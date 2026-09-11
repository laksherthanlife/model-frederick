from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.parse
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path("data/public_inputs.json")
CHUNK_SIZE = 1024 * 1024


class NetworkFailure(RuntimeError):
    pass


class HTTPSOnlyRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise NetworkFailure("refusing a public-input redirect away from HTTPS")
        return super().redirect_request(request, response, code, message, headers, newurl)


def relative_path(value):
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"invalid relative input path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value or value == ".":
        raise ValueError(f"input path is not canonical and relative: {value!r}")
    return Path(*path.parts)


def digest_text(value, length=64):
    if not isinstance(value, str) or len(value) != length or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"invalid {length}-character digest: {value!r}")
    return value


def nonsymlink_path(path):
    path = Path(path).absolute()
    for part in (path, *path.parents):
        if part.is_symlink():
            raise ValueError(f"public input path contains a symlink: {part}")
    return path


def external_directory(path, root):
    path = Path(path).expanduser().resolve()
    root = Path(root).resolve()
    if path.is_relative_to(root) or root.is_relative_to(path):
        raise ValueError("cache and report directories must be outside and not contain the checkout")
    return path


def verify_file(path, entry):
    path = nonsymlink_path(path)
    if not path.is_file():
        raise FileNotFoundError(f"mandatory public input is absent: {path}")
    size = path.stat().st_size
    mismatches = []
    if entry.get("size_bytes") is not None and size != entry["size_bytes"]:
        mismatches.append(f"byte count mismatch: expected {entry['size_bytes']}, got {size}")
    hashes = {"sha256": hashlib.sha256()}
    if "md5" in entry:
        hashes["md5"] = hashlib.md5(usedforsecurity=False)
    if "git_blob" in entry:
        hashes["git_blob"] = hashlib.sha1(f"blob {size}\0".encode(), usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(CHUNK_SIZE):
            for digest in hashes.values():
                digest.update(block)
    for name, digest in hashes.items():
        if digest.hexdigest() != entry[name]:
            mismatches.append(f"{name} mismatch: expected {entry[name]}, got {digest.hexdigest()}")
    if mismatches:
        raise ValueError(f"identity mismatch for {path}: " + "; ".join(mismatches))
    if "uncompressed_sha256" in entry:
        with gzip.open(path, "rb") as handle:
            actual = hashlib.file_digest(handle, "sha256").hexdigest()
        if actual != entry["uncompressed_sha256"]:
            raise ValueError(f"uncompressed SBML SHA-256 mismatch for {path}")
    return size


def source_url(source, entry):
    kind = source["kind"]
    if kind == "checkout":
        return None
    if kind == "checksum":
        url = source["urls"][entry["path"]]
        parsed = urllib.parse.urlsplit(url)
        if (parsed.scheme != "https" or not parsed.hostname or parsed.username is not None
                or parsed.password is not None or parsed.fragment):
            raise ValueError("checksum-pinned public inputs require an explicit credential-free HTTPS URL")
        return url
    if kind not in {"git", "zenodo", "kaggle"}:
        raise ValueError(f"unsupported public input source kind: {kind!r}")
    remote = relative_path(entry["remote_path"]).as_posix()
    if kind == "git":
        repository = relative_path(source["repository"]).as_posix()
        if len(PurePosixPath(repository).parts) != 2:
            raise ValueError("a GitHub source requires an owner/repository identity")
        commit = digest_text(source["commit"], 40)
        base = f"https://raw.githubusercontent.com/{repository}/{commit}/"
        suffix = urllib.parse.quote(remote, safe="/")
    elif kind == "zenodo":
        record = source["record"]
        if not isinstance(record, str) or not record.isdecimal() or int(record) <= 0:
            raise ValueError("a Zenodo source requires a version-specific numeric record")
        base = f"https://zenodo.org/api/records/{record}/files/"
        suffix = urllib.parse.quote(remote, safe="") + "/content"
    elif kind == "kaggle":
        dataset = relative_path(source["dataset"]).as_posix()
        version = source["version"]
        if len(PurePosixPath(dataset).parts) != 2 or type(version) is not int or version <= 0:
            raise ValueError("a Kaggle source requires a dataset identity and a positive version")
        base = f"https://www.kaggle.com/api/v1/datasets/download/{dataset}/"
        suffix = urllib.parse.quote(remote, safe="") + f"?datasetVersionNumber={version}"
    else:
        raise ValueError(f"unsupported public input source kind: {kind!r}")
    if source["base_url"] != base:
        raise ValueError("source URL does not match its immutable source identity")
    return base + suffix


def load_manifest(root=ROOT, manifest_path=None):
    root = Path(root).resolve()
    path = Path(manifest_path) if manifest_path is not None else root / MANIFEST
    manifest = json.loads(nonsymlink_path(path).read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported public-input manifest schema")
    if not isinstance(manifest.get("files"), list) or not manifest["files"]:
        raise ValueError("public-input manifest must enumerate mandatory files")
    if manifest.get("sources", {}).get("checkout", {}).get("kind") != "checkout":
        raise ValueError("the checkout source cannot be replaced by a downloader")
    for entry in manifest["files"]:
        validate_entry(entry, manifest["sources"])
    for collection in manifest.get("collections", []):
        relative_path(collection["path"])
        digest_text(collection["sha256"])
        if collection["format"] not in {"hog_sources", "portable_native"}:
            raise ValueError("unsupported public source collection")
    for name, value in manifest.get("environment", {}).items():
        if name not in {"YSTWIN_YEAST_GEM", "YSTWIN_EC_YEAST_GEM", "YSTWIN_THERMO"}:
            raise ValueError(f"unsupported public input resolver: {name}")
        if relative_path(value).parts[0] != "data":
            raise ValueError("public input resolver must point inside data/")
    return manifest


def validate_entry(entry, sources):
    path = relative_path(entry["path"])
    storage = entry.get("storage", "checkout")
    if storage not in {"checkout", "equilibrator_cache"}:
        raise ValueError(f"unknown input storage: {storage}")
    if storage == "checkout" and path.parts[0] != "data":
        raise ValueError("checkout public inputs must live under data/")
    if storage == "equilibrator_cache" and len(path.parts) != 1:
        raise ValueError("eQuilibrator cache input must be a single filename")
    digest_text(entry["sha256"])
    for name, length in (("md5", 32), ("git_blob", 40), ("uncompressed_sha256", 64)):
        if name in entry:
            digest_text(entry[name], length)
    source = sources[entry["source"]]
    size = entry.get("size_bytes")
    if size is not None and (type(size) is not int or size <= 0):
        raise ValueError("input byte counts must be positive integers")
    url = source_url(source, entry)
    if url is not None and size is None:
        raise ValueError("downloadable inputs require a byte count as well as SHA-256")


def inventory(root, manifest):
    entries = list(manifest["files"])
    errors = []
    for collection in manifest.get("collections", []):
        entry = {"path": collection["path"], "sha256": collection["sha256"], "source": "checkout"}
        entries.append(entry)
        try:
            path = root / relative_path(entry["path"])
            verify_file(path, entry)
            document = json.loads(path.read_text(encoding="utf-8"))
            if collection["format"] == "hog_sources":
                members = ((item["filename"], item["sha256"], item["bytes"])
                           for item in document["assets"].values())
            else:
                members = ((item["public_path"], item["public_sha256"], item["public_size_bytes"])
                           for item in document["records"])
            for name, digest, size in members:
                member = {"path": (relative_path(entry["path"]).parent / relative_path(name)).as_posix(),
                          "sha256": digest, "size_bytes": size, "source": "checkout"}
                validate_entry(member, manifest["sources"])
                entries.append(member)
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(f"{entry['path']}: {error}")
    identities = [(entry.get("storage", "checkout"), entry["path"]) for entry in entries]
    if len(identities) != len(set(identities)):
        raise ValueError("duplicate public input destinations")
    return entries, errors


def input_destination(root, entry):
    if entry.get("storage") == "equilibrator_cache":
        import pooch

        directory = Path(pooch.os_cache("equilibrator")).resolve()
        return nonsymlink_path(directory / relative_path(entry["path"]))
    return nonsymlink_path(root / relative_path(entry["path"]))


def tracked_paths(root):
    result = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, check=True)
    return {os.fsdecode(name) for name in result.stdout.split(b"\0") if name}


def authorize_destination(root, destination, tracked):
    if not destination.is_relative_to(root):
        return
    relative = destination.relative_to(root).as_posix()
    if relative in tracked:
        raise ValueError(f"refusing to provision over tracked input: {relative}")
    result = subprocess.run(["git", "check-ignore", "--quiet", "--", relative], cwd=root, check=False)
    if result.returncode != 0:
        raise ValueError(f"download destination must be explicitly gitignored: {relative}")


def copy_verified(source, destination, entry):
    nonsymlink_path(destination)
    if destination.exists():
        verify_file(destination, entry)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="public-input-", dir=destination.parent) as directory:
        staged = Path(directory) / "payload"
        shutil.copyfile(source, staged)
        verify_file(staged, entry)
        os.link(staged, destination)


def download_verified(url, destination, entry):
    nonsymlink_path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"Accept-Encoding": "identity", "User-Agent": "ystwin-public-inputs/1"})
    opener = urllib.request.build_opener(HTTPSOnlyRedirect())
    with tempfile.TemporaryDirectory(prefix="public-input-", dir=destination.parent) as directory:
        staged = Path(directory) / "payload"
        try:
            with opener.open(request, timeout=120) as response, staged.open("xb") as handle:
                total = 0
                while block := response.read(CHUNK_SIZE):
                    total += len(block)
                    if total > entry["size_bytes"]:
                        raise ValueError(f"download exceeds pinned byte count: {url}")
                    handle.write(block)
        except (urllib.error.URLError, OSError) as error:
            if isinstance(error, urllib.error.HTTPError):
                error.close()
            raise NetworkFailure(f"download stopped for {url}: {error}") from error
        verify_file(staged, entry)
        os.link(staged, destination)


def provision(root, manifest, cache_dir=None, *, verify_only=False):
    root = Path(root).resolve(strict=True)
    if not verify_only and cache_dir is None:
        raise ValueError("provisioning requires an explicit external cache directory")
    cache_dir = external_directory(cache_dir, root) if cache_dir is not None else None
    entries, errors = inventory(root, manifest)
    tracked = tracked_paths(root)
    for name in sorted(tracked):
        if name.startswith("data/"):
            try:
                path = nonsymlink_path(root / name)
                if not path.is_file():
                    raise FileNotFoundError(f"mandatory tracked public input is absent: {name}")
            except (OSError, ValueError) as error:
                errors.append(str(error))
    results = []
    network_failure = None
    for entry in entries:
        row = {"path": entry["path"], "storage": entry.get("storage", "checkout"),
               "sha256": entry["sha256"], "source": entry["source"]}
        source = manifest["sources"][entry["source"]]
        try:
            row["url"] = source_url(source, entry)
            destination = input_destination(root, entry)
            if destination.exists():
                verify_file(destination, entry)
                row["status"] = "verified"
                if not verify_only and source["kind"] != "checkout" and row["storage"] == "checkout":
                    copy_verified(destination, cache_dir / "sha256" / entry["sha256"], entry)
            elif verify_only or source["kind"] == "checkout":
                raise FileNotFoundError(f"mandatory public input is absent: {destination}")
            else:
                authorize_destination(root, destination, tracked)
                cached = (destination if row["storage"] == "equilibrator_cache"
                          else cache_dir / "sha256" / entry["sha256"])
                nonsymlink_path(cached)
                if cached.exists():
                    verify_file(cached, entry)
                else:
                    if network_failure is not None:
                        raise NetworkFailure(f"network operations already stopped: {network_failure}")
                    download_verified(row["url"], cached, entry)
                if cached != destination:
                    copy_verified(cached, destination, entry)
                row["status"] = "provisioned"
            row["size_bytes"] = destination.stat().st_size
        except (OSError, ValueError, ImportError, NetworkFailure) as error:
            row.update(status="failed", error=str(error))
            if isinstance(error, NetworkFailure) and network_failure is None:
                network_failure = str(error)
        results.append(row)
    for name, default in manifest.get("environment", {}).items():
        override = os.environ.get(name)
        if not override:
            continue
        expected = relative_path(default)
        for entry in entries:
            if entry.get("storage", "checkout") != "checkout":
                continue
            relative = relative_path(entry["path"])
            if relative != expected and relative.parent != expected:
                continue
            candidate = Path(override).expanduser()
            if relative.parent == expected:
                candidate = candidate / relative.name
            try:
                verify_file(candidate, entry)
            except (OSError, ValueError) as error:
                errors.append(f"{name}: {error}")
    return {"schema_version": 1, "mode": "verify_only" if verify_only else "provision",
            "passed": not errors and all(row["status"] != "failed" for row in results),
            "errors": errors, "inputs": results,
            "source_manifest_discrepancies": [
                {"path": entry["path"], "documented_uncompressed_sha256": entry["documented_uncompressed_sha256"],
                 "baseline_uncompressed_sha256": entry["uncompressed_sha256"], "evidence": entry["evidence"]}
                for entry in entries if "documented_uncompressed_sha256" in entry
                and entry["documented_uncompressed_sha256"] != entry["uncompressed_sha256"]],
            "optional_private_inputs": manifest.get("optional_private_inputs", [])}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Provision mandatory public inputs with pinned checksums; never accept refreshed bytes or hashes")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    output = None
    try:
        root = args.root.resolve(strict=True)
        output = external_directory(args.output_dir, root) if args.output_dir is not None else None
        manifest = load_manifest(root, args.manifest)
        report = provision(root, manifest, args.cache_dir, verify_only=args.verify_only)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        report = {"schema_version": 1, "passed": False, "errors": [str(error)], "inputs": []}
    if output is not None:
        output.mkdir(parents=True, exist_ok=True)
        with (output / "public-inputs.json").open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    for row in report["inputs"]:
        print(f"{row['status']}: {row['path']}" + (f": {row['error']}" if "error" in row else ""))
    for error in report["errors"]:
        print(f"failed: {error}")
    print(f"Public inputs: {len(report['inputs'])} checked; {'PASS' if report['passed'] else 'FAIL'}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
