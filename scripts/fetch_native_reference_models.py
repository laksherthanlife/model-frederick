from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import subprocess
import urllib.parse
import urllib.request
import zipfile

from ystwin.analysis.parameter_evidence import load_parameter_evidence


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIRECTORY = PurePosixPath("data/native_reference_models")


class AcquisitionError(RuntimeError):
    pass


def _digest(value, length=64):
    if not isinstance(value, str) or len(value) != length or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"expected a lowercase {length}-character digest")
    return value


def _relative(value):
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("expected a canonical relative asset path")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value or value == ".":
        raise ValueError("expected a canonical relative asset path")
    return path


def _url(value):
    parsed = urllib.parse.urlsplit(value)
    if (parsed.scheme != "https" or not parsed.hostname or parsed.username is not None
            or parsed.password is not None or parsed.fragment):
        raise ValueError("source acquisition requires a credential-free HTTPS URL")
    return value


def destination(root, relative):
    relative = _relative(relative)
    if not relative.is_relative_to(ASSET_DIRECTORY) or relative == ASSET_DIRECTORY:
        raise ValueError("writes are restricted to explicit data/native_reference_models/ assets")
    path = Path(root).absolute().joinpath(*relative.parts)
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError("source asset destinations must not contain symlinks")
    return path


def asset_entries(inventory, source_ids=None):
    data = inventory.to_dict()
    selected = set(data["sources"] if source_ids is None else source_ids)
    if not selected <= data["sources"].keys():
        raise ValueError("unknown native source id")
    entries = []
    for source_id, source in data["sources"].items():
        if source_id not in selected:
            continue
        for asset in source["local_artifacts"]:
            if "acquisition" not in asset:
                continue
            entry = {**asset, "source_id": source_id}
            path = _relative(entry["path"])
            if not path.is_relative_to(ASSET_DIRECTORY) or path == ASSET_DIRECTORY:
                raise ValueError("acquisition cannot target an existing non-native-reference asset")
            _digest(entry["sha256"])
            if type(entry["bytes"]) is not int or entry["bytes"] <= 0:
                raise ValueError("source assets require a positive pinned byte count")
            spec = entry["acquisition"]
            if spec["kind"] == "github":
                identity = source["identity"]
                repository = urllib.parse.urlsplit(_url(identity["repository"]))
                parts = _relative(repository.path.lstrip("/")).parts
                if repository.netloc != "github.com" or repository.query or len(parts) != 2:
                    raise ValueError("GitHub assets require an explicit owner/repository identity")
                entry["repository"] = "/".join(parts)
                entry["commit"] = _digest(identity["commit"], 40)
                _relative(spec["path"])
                _digest(spec["git_blob"], 40)
            elif spec["kind"] == "https":
                _url(spec["url"])
            elif spec["kind"] == "zip_member":
                _relative(spec["archive_path"])
                _relative(spec["member"])
            else:
                raise ValueError("unsupported native source acquisition kind")
            entries.append(entry)
    if not entries:
        raise ValueError("the selection has no pinned native-reference assets to verify or acquire")
    paths = [entry["path"] for entry in entries]
    if len(paths) != len(set(paths)):
        raise ValueError("duplicate native source asset destination")
    by_path = {entry["path"]: entry for entry in entries}
    for entry in entries:
        spec = entry["acquisition"]
        if spec["kind"] == "zip_member":
            archive = by_path.get(spec["archive_path"])
            if archive is None or archive["acquisition"]["kind"] != "https":
                raise ValueError("ZIP members require an explicitly registered HTTPS archive")
            if archive["source_id"] != entry["source_id"]:
                raise ValueError("ZIP members must belong to the same source as their archive")
    return entries


def verify_payload(payload, entry):
    if len(payload) != entry["bytes"]:
        raise ValueError(f"byte count mismatch for {entry['path']}")
    if hashlib.sha256(payload).hexdigest() != entry["sha256"]:
        raise ValueError(f"SHA-256 mismatch for {entry['path']}; never refresh the pinned hash")
    spec = entry["acquisition"]
    if spec["kind"] == "github":
        actual = hashlib.sha1(f"blob {len(payload)}\0".encode() + payload, usedforsecurity=False).hexdigest()
        if actual != spec["git_blob"]:
            raise ValueError(f"Git blob mismatch for {entry['path']}")


class _HTTPSRedirects(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, newurl):
        _url(newurl)
        return super().redirect_request(request, response, code, message, headers, newurl)


def download_https(url, expected_bytes):
    request = urllib.request.Request(_url(url), headers={
        "Accept-Encoding": "identity", "User-Agent": "ystwin-native-reference-acquisition/1",
    })
    # No authentication, alternative mirrors, retries, challenge solving, or hash refresh.
    with urllib.request.build_opener(_HTTPSRedirects()).open(request, timeout=60) as response:
        payload = response.read(expected_bytes + 1)
    if len(payload) != expected_bytes:
        raise AcquisitionError("HTTPS response differs from the pinned byte count (possibly an HTML challenge)")
    return payload


def download_github(entry):
    spec = entry["acquisition"]
    remote = urllib.parse.quote(spec["path"], safe="/")
    endpoint = f"repos/{entry['repository']}/contents/{remote}?ref={entry['commit']}"
    result = subprocess.run(["gh", "api", "--method", "GET", endpoint],
                            capture_output=True, check=True, timeout=60)
    response = json.loads(result.stdout)
    if (response.get("type") != "file" or response.get("path") != spec["path"]
            or response.get("sha") != spec["git_blob"] or response.get("size") != entry["bytes"]
            or response.get("encoding") != "base64"):
        raise AcquisitionError("GitHub commit/path/blob identity differs from the pinned source")
    return base64.b64decode("".join(response["content"].split()), validate=True)


def acquire(inventory, *, root=None, source_ids=None, fetch=False):
    root = ROOT if root is None else Path(root)
    entries = asset_entries(inventory, source_ids)
    by_path = {entry["path"]: entry for entry in entries}
    results = []
    # Archives precede their explicitly enumerated members; never extractall().
    ordered = sorted(entries, key=lambda entry: entry["acquisition"]["kind"] == "zip_member")
    for entry in ordered:
        path = destination(root, entry["path"])
        row = {"source_id": entry["source_id"], "path": entry["path"],
               "sha256": entry["sha256"], "bytes": entry["bytes"]}
        try:
            if path.exists():
                verify_payload(path.read_bytes(), entry)
                row["status"] = "verified_existing"
            elif not fetch:
                raise FileNotFoundError(f"missing pinned source asset: {entry['path']}")
            else:
                spec = entry["acquisition"]
                if spec["kind"] == "github":
                    payload = download_github(entry)
                elif spec["kind"] == "https":
                    payload = download_https(spec["url"], entry["bytes"])
                else:
                    archive_entry = by_path[spec["archive_path"]]
                    archive = destination(root, archive_entry["path"]).read_bytes()
                    verify_payload(archive, archive_entry)
                    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
                        if bundle.namelist().count(spec["member"]) != 1:
                            raise AcquisitionError("missing or duplicate ZIP member")
                        info = bundle.getinfo(spec["member"])
                        if info.is_dir() or info.file_size != entry["bytes"] or info.flag_bits & 1:
                            raise AcquisitionError("unexpected size, directory, or encrypted ZIP member")
                        payload = bundle.read(info)
                verify_payload(payload, entry)
                path.parent.mkdir(parents=True, exist_ok=True)
                destination(root, entry["path"])
                # Exclusive create: never overwrite, delete, restamp, or normalize source bytes.
                with path.open("xb") as handle:
                    handle.write(payload)
                row["status"] = "acquired_verified"
        except (OSError, ValueError, KeyError, AcquisitionError, subprocess.SubprocessError, zipfile.BadZipFile) as error:
            row.update(status="failed", error=f"{type(error).__name__}: {error}")
            results.append(row)
            # Stop acquisition on the first failed source instead of trying another route.
            if fetch:
                break
            continue
        results.append(row)
    return {"mode": "fetch" if fetch else "verify_only", "assets": results,
            "passed": len(results) == len(entries) and all(row["status"] != "failed" for row in results),
            "independent_real_validation_established": False}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Acquire only pinned native source assets; default is read-only verification")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--source", action="append", dest="source_ids")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--fetch", action="store_true", help="explicitly allow exclusive creation of missing pinned assets")
    mode.add_argument("--verify-only", action="store_false", dest="fetch", help="read-only verification (the default)")
    parser.set_defaults(fetch=False)
    args = parser.parse_args(argv)
    try:
        inventory = load_parameter_evidence(args.root / "data/parameter_evidence.json")
        report = acquire(inventory, root=args.root, source_ids=args.source_ids, fetch=args.fetch)
    except (OSError, ValueError, KeyError, TypeError) as error:
        report = {"passed": False, "error": f"{type(error).__name__}: {error}"}
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
