from __future__ import annotations

import ast
import csv
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import runpy
import shlex
import shutil
import site
import stat
import subprocess
import sys
import tempfile
import tomllib
import traceback
import uuid


REGISTRY_PATH = "data/artifact_registry.json"
LIFECYCLES = {"current", "historical_frozen", "local_exploratory", "superseded", "investigation_pending"}
DIMENSIONS = ("values", "missingness", "units", "row_keys")
PRIVATE_POLICY_PATH = "data/private_reproduction_policy.json"
_PRIVATE_PERMISSION = "explicit_env_file_read_only"
_PRIVATE_CI_SCOPE = "skip_unless_explicitly_reverified"
# Capture the launcher's OS temporary roots before the child redirects its own TMPDIR.
_PRIVATE_TEMPORARY_ROOTS = {Path(tempfile.gettempdir()).resolve(), Path("/tmp").resolve()}
# Source types are an allowlist, not a heuristic for excusing missing inputs.
# The policy additionally permits exact files/corpora for each registered recipe.
_PRIVATE_SOURCES = {
    "data/gen5": {"environment": "YSTWIN_GEN5_XPT", "suffix": ".xpt", "resolver": "gen5_xpt_dir"},
    "data/qpcr_raw": {"environment": "YSTWIN_QPCR_RAW", "suffix": ".xlsx", "resolver": "qpcr_raw_dir"},
    "data/qpcr_repaired": {"environment": "YSTWIN_QPCR", "suffix": ".xlsx", "resolver": "qpcr_dir",
                           "repository_root": "data/qpcr_repaired"},
}


class ArtifactError(ValueError):
    pass


def artifact_id(value: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or ":" in value:
        raise ArtifactError(f"not a repository-relative artifact identity: {value!r}")
    path = PurePosixPath(value)
    if (path.is_absolute() or len(path.parts) < 2 or path.as_posix() != value
            or any(part in {".", ".."} for part in value.split("/"))
            or any(character in value for character in "*?[]\x00")):
        raise ArtifactError(f"not a canonical full artifact identity: {value!r}")
    return value


def artifact_path(root: Path, identity: str) -> Path:
    root = Path(root).resolve()
    relative = artifact_id(identity)
    path = root
    for part in PurePosixPath(relative).parts:
        path = path / part
        if path.is_symlink():
            raise ArtifactError(f"symlink is not an artifact identity: {relative}")
    if not path.resolve().is_relative_to(root):
        raise ArtifactError(f"artifact escapes repository: {relative}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ArtifactError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _nonfinite(value):
    raise ArtifactError(f"non-finite JSON number: {value}")


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_object,
                      parse_constant=_nonfinite)


def canonical_digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def producer_digest(root: Path, paths: list[Path]) -> str:
    records = []
    for path in paths:
        relative = path.relative_to(root).as_posix()
        records.append({"path": artifact_id(relative), "sha256": sha256(artifact_path(root, relative))})
    return canonical_digest(sorted(records, key=lambda record: record["path"]))


def code_dependencies(root: Path, producer: dict) -> list[str]:
    pending = [producer["path"], *producer.get("libraries", [])]
    seen = set()
    while pending:
        relative = artifact_id(pending.pop())
        if relative in seen:
            continue
        seen.add(relative)
        path = artifact_path(root, relative)
        if not path.is_file():
            raise ArtifactError(f"missing producer code: {relative}")
        tree = ast.parse(path.read_text(encoding="utf-8"))
        dotted = ".".join(PurePosixPath(relative).with_suffix("").parts[1:]) if relative.startswith("src/") else ""
        package = dotted.removesuffix(".__init__") if path.name == "__init__.py" else dotted.rpartition(".")[0]
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = package.split(".")[:len(package.split(".")) - node.level + 1] if node.level else []
                prefix = ".".join([*base, *([node.module] if node.module else [])])
                names.add(prefix)
                names.update(f"{prefix}.{alias.name}" for alias in node.names)
            elif isinstance(node, ast.Call) and node.args:
                name = getattr(node.func, "attr", getattr(node.func, "id", ""))
                if name == "import_module" and isinstance(node.args[0], ast.Constant):
                    if isinstance(node.args[0].value, str):
                        names.add(node.args[0].value)
        for name in names:
            parts = name.split(".")
            if parts[0] != "ystwin":
                base = PurePosixPath("scripts", *parts)
                for candidate in (base.with_suffix(".py"), base / "__init__.py"):
                    if (root / candidate).is_file() and candidate.as_posix() not in seen:
                        pending.append(candidate.as_posix())
                continue
            for length in range(1, len(parts) + 1):
                base = PurePosixPath("src", *parts[:length])
                candidates = [base / "__init__.py"]
                if length == len(parts):
                    candidates.append(base.with_suffix(".py"))
                pending.extend(candidate.as_posix() for candidate in candidates
                               if (root / candidate).is_file() and candidate.as_posix() not in seen)
    return sorted(seen)


def _check_hash(value, label):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ArtifactError(f"{label} needs a full SHA-256 identity")


def git_tracking(root: Path) -> dict:
    process = subprocess.run(["git", "ls-files", "-z", "--", "outputs"], cwd=root, capture_output=True, text=True)
    available = process.returncode == 0
    tracked = set(process.stdout.split("\0")) - {""} if available else set()
    present = {path.relative_to(root).as_posix() for path in (root / "outputs").rglob("*") if path.is_file()}
    candidates = sorted(present - tracked)
    ignored = set()
    if available and candidates:
        process = subprocess.run(["git", "check-ignore", "-z", "--stdin"], cwd=root,
                                 input="\0".join(candidates) + "\0", capture_output=True, text=True)
        if process.returncode in (0, 1):
            ignored = set(process.stdout.split("\0")) - {""}
    return {"available": available, "tracked": tracked, "present": present, "ignored": ignored}


def json_pointer(payload, pointer: str):
    if pointer == "":
        return payload
    if not pointer.startswith("/"):
        raise ArtifactError(f"invalid JSON pointer: {pointer}")
    for token in pointer[1:].split("/"):
        key = token.replace("~1", "/").replace("~0", "~")
        payload = payload[int(key)] if isinstance(payload, list) else payload[key]
    return payload


def _expand_contract(contract: dict, schemas: dict, seen=()) -> dict:
    if not isinstance(contract, dict):
        raise ArtifactError("semantic contract must be an object")
    parent = contract.get("extends")
    base = {}
    if parent is not None:
        if parent in seen or parent not in schemas:
            raise ArtifactError(f"invalid semantic schema inheritance: {parent}")
        base = _expand_contract(schemas[parent], schemas, (*seen, parent))
    merged = {**base, **{key: value for key, value in contract.items() if key != "extends"}}
    units = dict(base.get("units", {}))
    units.update(contract.get("units", {}))
    for unit, columns in contract.get("unit_groups", {}).items():
        for column in columns:
            if column in units and units[column] != unit:
                raise ArtifactError(f"conflicting units for {column}")
            units[column] = unit
    merged["units"] = units
    merged.pop("unit_groups", None)
    return merged


def _contract(document: dict, identity: str, item: dict):
    direct = item.get("comparison")
    named = document.get("contracts", {}).get(identity)
    if direct is not None and named is not None:
        raise ArtifactError(f"two semantic contracts claim {identity}")
    contract = direct if direct is not None else named
    if isinstance(contract, str):
        if contract not in document.get("schemas", {}):
            raise ArtifactError(f"unknown semantic schema {contract}: {identity}")
        contract = document["schemas"][contract]
    if contract is None:
        return None
    contract = _expand_contract(contract, document.get("schemas", {}))
    if "column_variants" in contract:
        contract["column_variants"] = [document["column_sets"][value] if isinstance(value, str) else value
                                       for value in contract["column_variants"]]
    if isinstance(contract.get("columns"), str):
        contract["columns"] = document["column_sets"][contract["columns"]]
    return contract


def load_registry(root: Path, registry_path: str = REGISTRY_PATH) -> dict:
    document = read_json(artifact_path(root, registry_path))
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ArtifactError("artifact registry requires schema_version 1")
    if not isinstance(document.get("runs"), list) or not isinstance(document.get("frozen_bundles", []), list):
        raise ArtifactError("artifact registry requires explicit runs and frozen_bundles lists")
    runs, artifacts = {}, {}

    def register(identity, record):
        artifact_path(root, identity)
        if identity in artifacts:
            raise ArtifactError(f"artifact identity registered twice: {identity}")
        artifacts[identity] = record

    for run in document["runs"]:
        if not isinstance(run, dict) or not isinstance(run.get("id"), str) or not run["id"]:
            raise ArtifactError("every run needs an id")
        run = {**document.get("defaults", {}), **run}
        dependencies = list(run.get("dependencies", []))
        for name in run.get("dependency_sets", []):
            if name not in document.get("dependency_sets", {}):
                raise ArtifactError(f"unknown explicit dependency set: {name}")
            dependencies.extend(document["dependency_sets"][name])
        run["dependencies"] = dependencies
        if (run.get("model_identity") or {}).get("kind") == "declared_code_and_models":
            run["model_identity"] = {**run["model_identity"], "model_paths": [entry["path"] for entry in dependencies if entry.get("role") == "model"]}
        if run["id"] in runs:
            raise ArtifactError(f"duplicate run id: {run['id']}")
        if run.get("lifecycle") not in LIFECYCLES:
            raise ArtifactError(f"unknown lifecycle: {run['id']}")
        if run.get("status") not in {"investigation_pending", "verified"}:
            raise ArtifactError(f"explicit verification status required: {run['id']}")
        if not isinstance(run.get("artifacts"), list) or not run["artifacts"]:
            raise ArtifactError(f"run requires explicit artifact members: {run['id']}")
        run["artifacts"] = list(run["artifacts"])
        for index in run.get("artifact_indexes", []):
            path = artifact_path(root, index["path"])
            _check_hash(index.get("membership_sha256"), index["path"])
            if not path.is_file():
                run.setdefault("registry_issues", []).append(f"artifact membership index unavailable: {index['path']}")
                continue
            if sha256(path) != index["membership_sha256"]:
                raise ArtifactError(f"artifact membership index changed: {index['path']}")
            entries = read_json(path)
            if not isinstance(entries, list):
                raise ArtifactError(f"artifact membership index must be a list: {index['path']}")
            base = PurePosixPath(artifact_id(index["base_directory"]))
            for entry in entries:
                member = entry[index["path_field"]]
                if (not isinstance(member, str) or PurePosixPath(member).is_absolute()
                        or member != PurePosixPath(member).as_posix()):
                    raise ArtifactError(f"invalid indexed artifact: {member!r}")
                identity = artifact_id((base / member).as_posix())
                run["artifacts"].append({"path": identity, "comparison": index.get("comparison")})
        runs[run["id"]] = run
        expanded = []
        for item in run["artifacts"]:
            if not isinstance(item, dict):
                raise ArtifactError(f"artifact record must be an object: {run['id']}")
            identity = artifact_id(item.get("path"))
            item = {**item, "comparison": _contract(document, identity, item)}
            expanded.append(item)
            register(identity, {**item, "run": run["id"], "lifecycle": run["lifecycle"]})
        run["artifacts"] = expanded
    for bundle in document.get("frozen_bundles", []):
        manifest_id = artifact_id(bundle["manifest"])
        _check_hash(bundle.get("sha256"), manifest_id)
        manifest_path = artifact_path(root, manifest_id)
        if sha256(manifest_path) != bundle["sha256"]:
            raise ArtifactError(f"frozen manifest identity changed: {manifest_id}")
        manifest = read_json(manifest_path)
        if bundle.get("validator") == "native_v1_portable_evidence":
            from .analysis.portable_evidence import load_portable_evidence

            load_portable_evidence(root, manifest_id, bundle["sha256"])
        elif bundle.get("validator") is not None:
            raise ArtifactError(f"unknown frozen evidence validator: {bundle['validator']}")
        if not isinstance(manifest.get("records"), list):
            raise ArtifactError(f"frozen manifest requires exact records: {manifest_id}")
        frozen = {"lifecycle": "historical_frozen", "run": bundle["id"], "bundle": bundle}
        register(manifest_id, {**frozen, "sha256": bundle["sha256"]})
        origins = set()
        for item in manifest["records"]:
            origin = artifact_id(item["origin_path"])
            if origin in origins:
                raise ArtifactError(f"duplicate original identity in {manifest_id}: {origin}")
            origins.add(origin)
            public = (PurePosixPath(manifest_id).parent / item["public_path"]).as_posix()
            _check_hash(item.get("public_sha256"), public)
            _check_hash(item.get("original_sha256"), origin)
            register(public, {**frozen, "sha256": item["public_sha256"], "origin": origin})
            if origin.startswith("outputs/"):
                register(origin, {**frozen, "sha256": item["original_sha256"], "optional": True,
                                  "public_counterpart": public})
        if bundle.get("include_external_outputs") is True:
            for reference in manifest.get("external_references", []):
                origin = reference["origin_path"]
                if not origin.startswith("outputs/"):
                    continue
                _check_hash(reference.get("original_sha256"), origin)
                if origin in artifacts:
                    if artifacts[origin].get("sha256") != reference["original_sha256"]:
                        raise ArtifactError(f"conflicting frozen original identity: {origin}")
                    continue
                register(origin, {**frozen, "sha256": reference["original_sha256"], "optional": True,
                                  "public_counterpart": manifest_id, "integrity_scope": "original bytes if available; otherwise lineage only"})
    for local in document.get("local_runs", []):
        if local.get("lifecycle") != "local_exploratory" or not local.get("reason"):
            raise ArtifactError("local runs require an explicit exploratory disposition and reason")
        for identity in local.get("members", []):
            register(artifact_id(identity), {"lifecycle": "local_exploratory", "local_run": local,
                                             "run": local["id"], "retention": "local_only"})
    unknown_contracts = set(document.get("contracts", {})) - set(artifacts)
    if unknown_contracts:
        raise ArtifactError(f"contracts name unregistered artifacts: {sorted(unknown_contracts)}")
    return {"runs": runs, "artifacts": artifacts, "frozen_bundles": document.get("frozen_bundles", []),
            "receipts": document.get("receipts", {}), "reviews": document.get("reviews", {}), "document": document}


def inspect_contracts(root: Path, run_id: str | None = None, *, baseline: bool = False) -> list[dict]:
    registry = load_registry(root)
    result = []
    for run in registry["runs"].values():
        if run_id is not None and run["id"] != run_id:
            continue
        if run_id is None and run["lifecycle"] != "current":
            continue
        for item in run["artifacts"]:
            path = artifact_path(root, item["path"])
            record = {"run": run["id"], "path": item["path"], "producer": run.get("producer"),
                      "contract": item.get("comparison"), "present": path.is_file()}
            if path.is_file() and path.suffix == ".csv":
                with path.open(newline="", encoding="utf-8") as stream:
                    reader = csv.DictReader(stream)
                    record["columns"] = reader.fieldnames
                    samples = []
                    count = 0
                    for row in reader:
                        count += 1
                        if len(samples) < 2:
                            samples.append(row)
                    record.update(rows=count, examples=samples)
            if baseline:
                content, source = _baseline_bytes(root, registry, item, run)
                record["baseline_source"] = source
                record["baseline_sha256"] = hashlib.sha256(content).hexdigest() if content is not None else None
                if content is not None and path.suffix == ".csv":
                    reader = csv.DictReader(io.StringIO(content.decode("utf-8")))
                    record["baseline_columns"] = reader.fieldnames
                    record["baseline_examples"] = [row for _, row in zip(range(2), reader)]
            result.append(record)
    return result


def describe_local_run(root: Path, directory: str, metadata_path: str, producer: str,
                       assertions: dict, reason: str) -> dict:
    tracking = git_tracking(root)
    if not tracking["available"]:
        raise ArtifactError("local-only decisions require a Git index")
    directory_path = artifact_path(root, directory)
    members = sorted(path.relative_to(root).as_posix() for path in directory_path.rglob("*") if path.is_file())
    if not members or set(members) & tracking["tracked"]:
        raise ArtifactError("a local-only run must have explicit untracked members and no tracked artifacts")
    path = artifact_path(root, metadata_path)
    metadata = read_json(path)
    for pointer, expected in assertions.items():
        if json_pointer(metadata, pointer) != expected:
            raise ArtifactError(f"local run metadata does not establish {pointer}")
    artifact_id(producer)
    code_hashes = metadata.get("code_hashes", {})
    if producer not in code_hashes:
        raise ArtifactError("local run metadata does not identify the supplied producer")
    _check_hash(code_hashes[producer], producer)
    assertions = {**assertions, "/code_hashes/" + producer.replace("~", "~0").replace("/", "~1"): code_hashes[producer]}
    return {"id": directory, "lifecycle": "local_exploratory", "reason": reason,
            "metadata": {"path": metadata_path, "sha256": sha256(path)}, "assertions": assertions,
            "producer": {"path": producer, "identity_basis": "pinned run metadata, not current source"},
            "members": members}


def registered_writers(root: Path) -> dict[str, list[str]]:
    registry = load_registry(root)
    writers = {}
    for run in registry["runs"].values():
        producer = run.get("producer") or {}
        path = producer.get("path", "")
        if run["lifecycle"] != "current" or not path.startswith("scripts/"):
            continue
        destinations = writers.setdefault(path.removeprefix("scripts/"), [])
        destinations.extend(item["path"] for item in run["artifacts"] if item["path"].startswith("outputs/"))
    return {script: sorted(set(paths)) for script, paths in writers.items()}


def reproduction_command(run: dict) -> str:
    recipe = run.get("reproduction") or {}
    command = recipe.get("argv") or []
    if not command:
        return "unresolved: no responsible executable producer has been established"
    return shlex.join("python3" if value == "{python}" else value for value in command)


def run_metadata_issues(run: dict) -> list[str]:
    issues = list(run.get("registry_issues", []))
    producer = run.get("producer")
    if not isinstance(producer, dict) or not producer.get("path") or not producer.get("callable"):
        issues.append("responsible producer not established")
    else:
        artifact_id(producer["path"])
        for path in producer.get("libraries", []):
            artifact_id(path)
    dependencies = run.get("dependencies")
    if not isinstance(dependencies, list):
        issues.append("explicit input/output dependencies absent")
    else:
        seen = set()
        for dependency in dependencies:
            identity = artifact_id(dependency.get("path"))
            if identity in seen:
                issues.append(f"duplicate dependency: {identity}")
            seen.add(identity)
            if dependency.get("role") not in {"input", "model", "protocol", "upstream_output", "evidence"}:
                issues.append(f"dependency role absent: {identity}")
            if dependency.get("sha256") is not None:
                _check_hash(dependency["sha256"], identity)
        for item in run["artifacts"]:
            if item["path"] in seen:
                issues.append(f"input is also an output: {item['path']}")
    for key in ("parameters", "model_identity", "runtime"):
        if not isinstance(run.get(key), dict) or not run[key]:
            issues.append(f"{key} absent")
    recipe = run.get("reproduction")
    if not isinstance(recipe, dict) or not isinstance(recipe.get("argv"), list) or not recipe["argv"]:
        issues.append("reproduction command absent")
    elif not all(isinstance(value, str) and value for value in recipe["argv"]):
        issues.append("reproduction argv must contain nonempty strings")
    if not run.get("reason"):
        issues.append("lifecycle justification absent")
    if run["lifecycle"] == "superseded" and not run.get("superseded_by"):
        issues.append("superseded artifact needs its replacement identity")
    if run["lifecycle"] == "investigation_pending" and run["status"] != "investigation_pending":
        issues.append("an unresolved lifecycle cannot be verified")
    return issues


def _pinned_references(root: Path, references, label: str) -> list[str]:
    if not isinstance(references, list) or not references:
        return [f"{label} needs original full-path SHA-256 references"]
    issues = []
    for reference in references:
        if not isinstance(reference, dict) or "path" not in reference:
            issues.append(f"{label} reference is not a pinned identity")
            continue
        _check_hash(reference.get("sha256"), reference["path"])
        path = artifact_path(root, reference["path"])
        if not path.is_file() or sha256(path) != reference["sha256"]:
            issues.append(f"{label} identity not verified: {reference['path']}")
    return issues


def runtime_identity(root: Path) -> dict:
    result = {"python": platform.python_version(), "distributions": {}}
    project = root / "pyproject.toml"
    if project.is_file():
        payload = tomllib.loads(project.read_text(encoding="utf-8"))
        requirements = list(payload.get("project", {}).get("dependencies", []))
        for group in payload.get("project", {}).get("optional-dependencies", {}).values():
            requirements.extend(group)
        for requirement in requirements:
            name = requirement.split("==")[0].split("[")[0]
            try:
                result["distributions"][name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                result["distributions"][name] = None
        if result["distributions"].get("cobra") is not None:
            for name in ("optlang", "swiglpk"):
                try:
                    result["distributions"][name] = importlib.metadata.version(name)
                except importlib.metadata.PackageNotFoundError:
                    result["distributions"][name] = None
            cobra = importlib.import_module("cobra")
            result["cobra_solver"] = cobra.Configuration().solver.__name__
    return result


def _private_spec(root: Path, run: dict) -> tuple[dict, dict | None]:
    private = {entry["path"]: entry for entry in run["dependencies"] if any(
        PurePosixPath(entry["path"]).is_relative_to(base) for base in _PRIVATE_SOURCES)
        or (entry.get("source") or {}).get("environment") in {
            source["environment"] for source in _PRIVATE_SOURCES.values()}}
    if not private:
        return {}, None
    path = artifact_path(root, PRIVATE_POLICY_PATH)
    if not path.is_file():
        raise ArtifactError("private input policy unavailable; outside-repository input access is disabled")
    policy = read_json(path)
    if (not isinstance(policy, dict) or policy.get("schema_version") != 1
            or policy.get("permission") != _PRIVATE_PERMISSION
            or policy.get("ci_private_bytes") != _PRIVATE_CI_SCOPE
            or set(policy) != {"schema_version", "permission", "ci_private_bytes", "sources", "runs"}
            or not isinstance(policy["sources"], dict) or not isinstance(policy["runs"], dict)):
        raise ArtifactError("invalid explicit private input policy")
    for base, source in policy["sources"].items():
        if base not in _PRIVATE_SOURCES or source != _PRIVATE_SOURCES[base]:
            raise ArtifactError("private source type is outside the approved policy allowlist")
    permissions = policy["runs"].get(run["id"], {})
    if not isinstance(permissions, dict):
        raise ArtifactError("private recipe permissions must name exact logical inputs")
    if set(private) != set(permissions):
        raise ArtifactError("unapproved private input or private recipe membership changed: "
                            + ", ".join(sorted(set(private) ^ set(permissions))))
    rules = {}
    for name, permission in permissions.items():
        artifact_id(name)
        owners = [base for base in policy["sources"] if PurePosixPath(name).is_relative_to(base)]
        if len(owners) != 1 or not isinstance(permission, dict):
            raise ArtifactError(f"unapproved private input: {name}")
        base = owners[0]
        source = policy["sources"][base]
        dependency = private[name]
        kind = permission.get("kind")
        if kind == "file":
            if permission != {"kind": "file"} or PurePosixPath(name).suffix != source["suffix"]:
                raise ArtifactError(f"private input extension or permission is not approved: {name}")
        elif kind == "directory":
            if name != base or permission != {"kind": "directory", "glob": "*" + source["suffix"]}:
                raise ArtifactError(f"private corpus must select only top-level approved extensions: {name}")
        else:
            raise ArtifactError(f"private input kind is not approved: {name}")
        if dependency.get("kind", "file") != kind or dependency.get("role") != "input":
            raise ArtifactError(f"private input role/kind conflicts with policy: {name}")
        declared_source = {"environment": source["environment"], "resolver": source["resolver"]}
        if kind == "file":
            declared_source["relative_path"] = PurePosixPath(name).relative_to(base).as_posix()
        supplied_source = dependency.get("source", {})
        if not isinstance(supplied_source, dict) or any(
                key not in declared_source or value != declared_source[key] for key, value in supplied_source.items()):
            raise ArtifactError(f"private source conflicts with approved policy: {name}")
        rules[name] = {**permission, "logical_root": base, "source": declared_source}
    return rules, {"policy": {"path": PRIVATE_POLICY_PATH, "sha256": sha256(path)},
                   "permission": _PRIVATE_PERMISSION, "ci_private_bytes": _PRIVATE_CI_SCOPE,
                   "roots": sorted(rules)}


def private_input_declarations(root: Path, run: dict) -> dict:
    """Logical registry source declarations; never includes resolved private locations."""
    rules, _ = _private_spec(root, run)
    return {name: rule["source"] for name, rule in rules.items()}


def _no_symlinks(path: Path, label: str) -> Path:
    path = Path(path)
    if not path.is_absolute() or ".." in path.parts:
        raise ArtifactError(f"canonical absolute location required for {label}")
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        if current.is_symlink():
            raise ArtifactError(f"symlink access is disabled for {label}")
    return path


@contextmanager
def _private_open(path: Path, label: str, *, directory: bool = False):
    """Open read-only through no-follow directory descriptors, including every ancestor."""
    descriptors = []
    try:
        path = _no_symlinks(path, label)
        descriptor = os.open(path.anchor, os.O_RDONLY | os.O_DIRECTORY)
        descriptors.append(descriptor)
        for part in path.parts[1:-1]:
            descriptor = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=descriptor)
            descriptors.append(descriptor)
        descriptor = os.open(path.name, os.O_RDONLY | os.O_NOFOLLOW | (os.O_DIRECTORY if directory else os.O_NONBLOCK),
                             dir_fd=descriptor)
        descriptors.append(descriptor)
        info = os.fstat(descriptor)
        if not (stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode)):
            raise ArtifactError(f"private input is not a regular {'directory' if directory else 'file'}: {label}")
        if directory:
            yield descriptor
        else:
            with os.fdopen(os.dup(descriptor), "rb") as stream:
                yield stream
    except OSError:
        # No local root, OS exception filename, or raw workbook metadata may enter a receipt.
        raise ArtifactError(f"private input unavailable or unsafe: {label}") from None
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _private_sha256(path: Path, name: str) -> str:
    with _private_open(path, name) as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _private_locations(root: Path, rules: dict, env_file: Path | None) -> dict[str, Path]:
    if env_file is None:
        raise ArtifactError("private input access requires explicit --private-inputs-from-env <env-file>; ambient environment is not permission")
    env_file = Path(env_file)
    if not env_file.is_absolute():
        env_file = root / env_file
    names = {rule["source"]["environment"] for rule in rules.values()}
    settings = {}
    with _private_open(env_file, "private env file") as stream:
        for raw in stream:
            try:
                line = raw.decode("utf-8").strip()
            except UnicodeDecodeError:
                raise ArtifactError("private env file must be UTF-8 data, never shell code") from None
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[len("export "):].lstrip()
            name, separator, value = line.partition("=")
            name = name.strip()
            # Deliberately do not interpret, retain, or use unrelated settings.
            if name not in names:
                continue
            if not separator or name in settings:
                raise ArtifactError(f"private env file needs one assignment for {name}")
            value = value.strip()
            if value[:1] in {"'", '"'}:
                if len(value) < 2 or value[-1] != value[0] or value[0] in value[1:-1]:
                    raise ArtifactError(f"invalid quoted private source: {name}")
                value = value[1:-1]
            if (not value or any(character in value for character in "$`\\\x00\r\n\"'")
                    or not Path(value).is_absolute() or Path(value).as_posix() != value or ".." in Path(value).parts):
                raise ArtifactError(f"canonical literal absolute private source required: {name}")
            settings[name] = _no_symlinks(Path(value), name)
    locations = {}
    for rule in rules.values():
        base, environment = rule["logical_root"], rule["source"]["environment"]
        source = _PRIVATE_SOURCES[base]
        if environment in settings:
            locations[base] = settings[environment]
        elif source.get("repository_root"):
            locations[base] = artifact_path(root, source["repository_root"])
        else:
            raise ArtifactError(f"approved private source not configured: {environment}")
        with _private_open(locations[base], environment, directory=True):
            pass
    return locations


def _private_members(path: Path, name: str, rule: dict) -> dict[str, Path]:
    if rule["kind"] == "file":
        with _private_open(path, name):
            pass
        return {name: path}
    members = {}
    with _private_open(path, name, directory=True) as descriptor:
        # No recursive search, no inspection of unrelated files/subdirectories.
        names = sorted(value for value in os.listdir(descriptor) if value.endswith(_PRIVATE_SOURCES[rule["logical_root"]]["suffix"]))
        for member in names:
            logical = artifact_id((PurePosixPath(name) / member).as_posix())
            info = os.stat(member, dir_fd=descriptor, follow_symlinks=False)
            if stat.S_ISLNK(info.st_mode):
                raise ArtifactError(f"private corpus contains a symlink: {logical}")
            if not stat.S_ISREG(info.st_mode):
                raise ArtifactError(f"private corpus member is not a regular file: {logical}")
            members[logical] = path / member
    if not members:
        raise ArtifactError(f"required private input corpus is empty: {name}")
    return members


def _public_dependency_sources(root: Path, dependency: dict, *, staged: bool) -> dict[str, Path]:
    identity = artifact_id(dependency["path"])
    path = artifact_path(root, identity)
    source = dependency.get("source", {})
    supplied = os.environ.get(source.get("environment", ""))
    if not staged and source and (supplied or not path.exists()):
        if not supplied:
            raise ArtifactError(f"required in-repository input {identity} unavailable; supply an approved in-repository snapshot via {source.get('environment', 'registered source')}")
        supplied_path = Path(supplied).expanduser().absolute()
        if not supplied_path.resolve().is_relative_to(root.resolve()):
            raise ArtifactError(f"outside-repository input access is disabled: {identity}; supply an approved in-repository snapshot")
        path = _no_symlinks(supplied_path, identity)
        if source.get("relative_path"):
            relative = PurePosixPath(source["relative_path"])
            if relative.is_absolute() or ".." in relative.parts or relative.as_posix() != source["relative_path"]:
                raise ArtifactError(f"external input escapes its declared source: {identity}")
            path = _no_symlinks(path / relative, identity)
    if dependency.get("kind", "file") == "directory":
        if not path.is_dir():
            raise ArtifactError(f"required input directory unavailable: {identity}")
        members = {}
        for member in sorted(path.rglob("*")):
            if member.is_symlink():
                raise ArtifactError(f"input directory contains a symlink: {identity}")
            if member.is_file():
                logical = artifact_id((PurePosixPath(identity) / member.relative_to(path).as_posix()).as_posix())
                members[logical] = member
        if not members and not dependency.get("allow_empty", False):
            raise ArtifactError(f"required input directory is empty: {identity}")
        return members
    if not path.is_file():
        raise ArtifactError(f"required input unavailable: {identity}")
    return {identity: path}


def dependency_sources(root: Path, run: dict, *, staged: bool = False,
                       private_inputs_from_env: Path | None = None) -> dict[str, Path]:
    rules, _ = _private_spec(root, run)
    if staged and rules:
        _check_private_stage(root.parent)
        if root.name != "workspace" or not (root.parent / "execution_plan.json").is_file():
            raise ArtifactError("private staged identities require an isolated execution workspace")
    locations = _private_locations(root, rules, private_inputs_from_env) if rules and not staged else {}
    if private_inputs_from_env is not None and not rules:
        raise ArtifactError("recipe has no approved private inputs; env-file opt-in cannot authorize other sources")
    found = {}
    for dependency in run["dependencies"]:
        name = dependency["path"]
        rule = rules.get(name)
        if rule:
            path = artifact_path(root, name) if staged else locations[rule["logical_root"]] / rule["source"].get("relative_path", "")
            members = _private_members(path, name, rule)
        else:
            members = _public_dependency_sources(root, dependency, staged=staged)
        if set(members) & set(found):
            raise ArtifactError(f"input declared twice: {name}")
        if dependency.get("sha256") is not None:
            records = [{"path": key, "sha256": _private_sha256(value, key) if rule else sha256(value)}
                       for key, value in sorted(members.items())]
            actual = canonical_digest(records) if dependency.get("kind") == "directory" else records[0]["sha256"]
            if actual != dependency["sha256"]:
                raise ArtifactError(f"declared input/model identity changed: {name}")
        found.update(members)
    return found


def _dependency_owner(run: dict, name: str) -> dict:
    owners = [entry for entry in run["dependencies"] if name == entry["path"] or (
        entry.get("kind") == "directory" and PurePosixPath(name).is_relative_to(entry["path"]))]
    if len(owners) != 1:
        raise ArtifactError(f"ambiguous dependency role: {name}")
    return owners[0]


def _identity(root: Path, run: dict, dependencies: list, rules: dict, private: dict | None) -> dict:
    code = code_dependencies(root, run["producer"])
    if (root / "pyproject.toml").is_file():
        code.append("pyproject.toml")
    records = []
    for name in sorted(set(code)):
        path = root / name if name == "pyproject.toml" else artifact_path(root, name)
        records.append({"path": name, "sha256": sha256(path)})
    bound_roots = [{**entry, "source": rules[entry["path"]]["source"]} if entry["path"] in rules else entry
                   for entry in run["dependencies"]]
    result = {"producer": run["producer"], "code": records, "dependencies": sorted(dependencies, key=lambda entry: entry["path"]),
              "dependency_roots": bound_roots, "parameters": run["parameters"],
              "model_identity": run["model_identity"], "runtime": runtime_identity(root),
              "reproduction": run["reproduction"], "artifacts": run["artifacts"],
              "verifier": {"path": "src/ystwin/artifacts.py", "sha256": sha256(Path(__file__))}}
    if private is not None:
        result["private_inputs"] = private
    return result


def run_identity(root: Path, run: dict, *, staged: bool = False,
                 private_inputs_from_env: Path | None = None) -> dict:
    issues = run_metadata_issues(run)
    if issues:
        raise ArtifactError("; ".join(issues))
    rules, private = _private_spec(root, run)
    dependencies = []
    for name, source in sorted(dependency_sources(root, run, staged=staged, private_inputs_from_env=private_inputs_from_env).items()):
        owner = _dependency_owner(run, name)
        is_private = owner["path"] in rules
        entry = {"path": name, "role": owner["role"], "sha256": _private_sha256(source, name) if is_private else sha256(source)}
        if is_private:
            entry["private"] = True
        dependencies.append(entry)
    return _identity(root, run, dependencies, rules, private)


def _receipt_only_identity(root: Path, run: dict, recorded: dict) -> tuple[dict, list[dict]]:
    """Rehash public inputs/code, retaining ONLY policy-approved private pins from a receipt.

    This is not run_identity and must never be used to claim execution/reverification.
    Callers must first validate the durable receipt and its successful child binding.
    """
    rules, private = _private_spec(root, run)
    if not rules or recorded.get("private_inputs") != private:
        raise ArtifactError("private receipt policy identity changed or permission is absent")
    dependencies, skipped = [], []
    recorded_entries = recorded.get("dependencies", [])
    if not isinstance(recorded_entries, list):
        raise ArtifactError("private receipt needs an explicit input roster")
    for dependency in run["dependencies"]:
        name = dependency["path"]
        if name not in rules:
            for logical, source in _public_dependency_sources(root, dependency, staged=False).items():
                dependencies.append({"path": logical, "role": dependency["role"], "sha256": sha256(source)})
            continue
        rule = rules[name]
        members = [entry for entry in recorded_entries if isinstance(entry, dict) and isinstance(entry.get("path"), str)
                   and (entry["path"] == name or (rule["kind"] == "directory" and PurePosixPath(entry["path"]).is_relative_to(name)))]
        if not members or (rule["kind"] == "file" and len(members) != 1):
            raise ArtifactError(f"private receipt input membership is unavailable: {name}")
        seen = set()
        for entry in members:
            logical = artifact_id(entry["path"])
            if (logical in seen or entry.get("private") is not True or entry.get("role") != "input"
                    or set(entry) != {"path", "role", "sha256", "private"}
                    or (rule["kind"] == "directory" and (PurePosixPath(logical).parent.as_posix() != name
                        or PurePosixPath(logical).suffix != _PRIVATE_SOURCES[rule["logical_root"]]["suffix"]))):
                raise ArtifactError(f"unapproved private input in receipt: {logical}")
            _check_hash(entry.get("sha256"), logical)
            seen.add(logical)
        pin = dependency.get("sha256")
        if pin is not None:
            actual = members[0]["sha256"] if rule["kind"] == "file" else canonical_digest(
                sorted(({"path": entry["path"], "sha256": entry["sha256"]} for entry in members), key=lambda entry: entry["path"]))
            if pin != actual:
                raise ArtifactError(f"declared private input identity changed: {name}")
        dependencies.extend(members)
        skipped.extend(members)
    current = _identity(root, run, dependencies, rules, private)
    if len(dependencies) != len(recorded_entries):
        raise ArtifactError("receipt contains unknown private/public input paths")
    return current, sorted(skipped, key=lambda entry: entry["path"])


def _number_equal(left, right, rtol, atol):
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, int) and isinstance(right, int):
        return left == right
    try:
        a, b = Decimal(str(left)), Decimal(str(right))
    except (InvalidOperation, ValueError, TypeError):
        return left == right
    if not a.is_finite() or not b.is_finite():
        return not a.is_nan() and not b.is_nan() and a == b
    return abs(a - b) <= max(Decimal(str(atol)), Decimal(str(rtol)) * max(abs(a), abs(b)))


def _csv_table(path: Path, contract: dict):
    keys = contract["row_keys"]
    if not isinstance(keys, list) or not keys or not all(isinstance(key, str) for key in keys):
        raise ArtifactError("CSV row_keys must be a nonempty list of column names")
    missing = contract.get("missing_values", [""])
    nullable = set(contract.get("nullable_row_keys", []))
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream, strict=True)
        header = next(reader)
        for variant in contract.get("row_key_variants", []):
            if set(variant["when_columns"]) <= set(header):
                keys = variant["keys"]
                break
        if len(set(header)) != len(header) or not set(keys) <= set(header):
            raise ArtifactError("duplicate columns or missing row key columns")
        variants = contract.get("column_variants") or ([contract["columns"]] if "columns" in contract else [])
        if variants and set(header) not in [set(columns) for columns in variants]:
            raise ArtifactError(f"unregistered column schema: {header}")
        rows, groups = {}, {}
        for values in reader:
            if len(values) != len(header):
                raise ArtifactError("CSV row has wrong number of columns")
            row = dict(zip(header, values))
            numeric_keys = set(keys) & set(contract.get("numeric_columns", contract["units"]))
            key = tuple(str(Decimal(row[name]).normalize()) if name in numeric_keys and row[name] not in missing
                        else row[name] for name in keys)
            if any(row[name] in missing and name not in nullable for name in keys) or key in rows:
                raise ArtifactError(f"duplicate or missing row key: {key}")
            for name in contract.get("boolean_columns", []):
                if name in row and row[name] not in [*missing, "True", "False", "true", "false"]:
                    raise ArtifactError(f"invalid boolean: {name}")
            for name in contract.get("numeric_columns", contract["units"]):
                if name not in row or row[name] in missing:
                    continue
                try:
                    value = Decimal(row[name])
                except InvalidOperation as error:
                    raise ArtifactError(f"invalid numeric value in {name}: {row[name]}") from error
                if value.is_nan() or (not value.is_finite() and name not in contract.get("infinite_columns", [])):
                    raise ArtifactError(f"non-finite numeric value in {name}")
            for name, declaration in contract["units"].items():
                if name not in row or row[name] in missing or not isinstance(declaration, dict):
                    continue
                if "by" in declaration:
                    selector = row.get(declaration["by"])
                    if declaration.get("compound_separator") and isinstance(selector, str):
                        selector = selector.split(declaration["compound_separator"], 1)[0]
                    if selector == "" and any(row.get(column) != value for column, value in declaration.get("empty_selector_requires", {}).items()):
                        raise ArtifactError(f"empty unit selector is outside its declared context: {name}")
                    if selector not in declaration.get("cases", {}):
                        raise ArtifactError(f"undeclared conditional unit for {name}: {selector}")
                elif "column" in declaration and row.get(declaration["column"], "") in missing:
                    raise ArtifactError(f"missing unit for {name}")
            if contract.get("duplicate_keys") == "multiset":
                groups.setdefault(key, []).append(row)
            else:
                rows[key] = row
        for key, members in groups.items():
            for index, row in enumerate(sorted(members, key=lambda value: json.dumps(value, sort_keys=True))):
                rows[(*key, "occurrence", str(index))] = row
        if not rows and not contract.get("allow_empty", False):
            raise ArtifactError("empty result is not permitted by the semantic contract")
        return header, rows


def _negative_csv_changes(left, right, contract):
    changes = []
    for key, row in left.items():
        for column, rule in contract.get("negative_cases", {}).items():
            values = [str(value) for value in rule.get("values", [])]
            before = row.get(column)
            negative = before in values or (rule.get("nonempty") is True and before not in (None, ""))
            if rule.get("missing") is True and column in row and before in (None, ""):
                negative = True
            if rule.get("positive") is True and before not in (None, ""):
                negative = negative or Decimal(before) > 0
            after = right.get(key, {}).get(column)
            if negative and (key not in right or after != before):
                changes.append({"row_key": list(key), "column": column, "before": before,
                                "after": after, "row_retained": key in right})
    return changes


def _normalize_storage(payload, contract):
    payload = json.loads(json.dumps(payload, allow_nan=False))
    for pointer, rule in contract.get("storage_locations", {}).items():
        parent_pointer, _, last = pointer.rpartition("/")
        parent = json_pointer(payload, parent_pointer)
        key = int(last) if isinstance(parent, list) else last.replace("~1", "/").replace("~0", "~")
        value = parent[key]
        if not isinstance(value, str):
            raise ArtifactError(f"storage location is not a string: {pointer}")
        prefix = rule.get("prefix", "")
        if not value.startswith(prefix):
            raise ArtifactError(f"storage field has an unexpected prefix: {pointer}")
        location = value[len(prefix):]
        tail = ""
        if rule.get("suffix_separator"):
            location, separator, ending = location.rpartition(rule["suffix_separator"])
            if not separator:
                raise ArtifactError(f"storage field has no declared separator: {pointer}")
            tail = separator + ending
        suffix = rule["repo_relative_suffix"]
        if suffix not in {"outputs", "data", "figures"}:
            artifact_id(suffix)
        normalized = location.replace("\\", "/").rstrip("/")
        if normalized != suffix and not normalized.endswith("/" + suffix):
            raise ArtifactError(f"storage location does not identify {suffix}: {pointer}")
        parent[key] = prefix + suffix + tail
    return payload


def _pointer_matches(pointer: str, pattern: str) -> bool:
    actual, expected = pointer.split("/"), pattern.split("/")
    return len(actual) == len(expected) and all(a == b or b == "*" for a, b in zip(actual, expected))


def _json_row_keys(contract, pointer):
    rules = contract.get("row_keys", {})
    if not isinstance(rules, dict):
        return None
    matches = [keys for pattern, keys in rules.items() if _pointer_matches(pointer, pattern)]
    if len(matches) > 1:
        raise ArtifactError(f"ambiguous JSON row key rule: {pointer}")
    return matches[0] if matches else None


def _json_index(values, keys, contract, pointer):
    result = {}
    nullable = set(contract.get("nullable_json_row_keys", {}).get(pointer, []))
    for row in values:
        if not isinstance(row, dict):
            raise ArtifactError(f"JSON row is not an object: {pointer}")
        fields = []
        for key in keys:
            value = json_pointer(row, key) if key.startswith("/") else row[key]
            if value is None and key not in nullable:
                raise ArtifactError(f"missing JSON row key: {pointer}/{key}")
            fields.append(value)
        identity = json.dumps(fields, sort_keys=True, allow_nan=False)
        if identity in result:
            raise ArtifactError(f"duplicate JSON row key: {pointer}: {identity}")
        result[identity] = row
    return result


def _json_leaves(value, contract, pointer=""):
    if isinstance(value, list):
        keys = _json_row_keys(contract, pointer)
        value = _json_index(value, keys, contract, pointer) if keys else dict(enumerate(value))
    if isinstance(value, dict):
        leaves = {}
        for key, child in value.items():
            path = pointer + "/" + str(key).replace("~", "~0").replace("/", "~1")
            leaves.update(_json_leaves(child, contract, path))
        return leaves
    return {pointer: value}


def _negative_json_changes(left, right, contract):
    if not contract.get("negative_json_fields"):
        return []
    before, after = _json_leaves(left, contract), _json_leaves(right, contract)
    changes = []
    for pointer, value in before.items():
        for pattern, rule in contract.get("negative_json_fields", {}).items():
            negative = value in rule.get("values", []) or (rule.get("positive") is True and isinstance(value, (float, int)) and value > 0)
            if _pointer_matches(pointer, pattern) and negative and (pointer not in after or after[pointer] != value):
                changes.append({"json_pointer": pointer, "before": value, "after": after.get(pointer), "field_retained": pointer in after})
    return changes


def _metadata_values(payload, contract):
    values = {}
    for pointer, kind in contract.get("identity_fields", {}).items():
        value = json_pointer(payload, pointer)
        if kind in {"file_fingerprint", "file_fingerprints"}:
            entries = value if kind == "file_fingerprints" else [value]
            if not isinstance(entries, list):
                raise ArtifactError(f"identity field must be a fingerprint list: {pointer}")
            for entry in entries:
                if not isinstance(entry, dict) or set(entry) != {"path", "sha256", "bytes"}:
                    raise ArtifactError(f"identity field contains non-identity payload: {pointer}")
                if entry["path"] != "pyproject.toml":
                    artifact_id(entry["path"])
                _check_hash(entry["sha256"], pointer)
                if type(entry["bytes"]) is not int or entry["bytes"] < 0:
                    raise ArtifactError(f"invalid byte count in identity field: {pointer}")
        elif kind == "runtime_versions":
            if not isinstance(value, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in value.items()):
                raise ArtifactError(f"invalid runtime version map: {pointer}")
        elif kind == "python_version":
            if not isinstance(value, str):
                raise ArtifactError(f"invalid interpreter identity: {pointer}")
        else:
            raise ArtifactError(f"unknown identity-field schema: {kind}")
        values[pointer] = value
    return values


def _mask_metadata(payload, contract):
    values = _metadata_values(payload, contract)
    for pointer in values:
        parent_pointer, _, last = pointer.rpartition("/")
        parent = json_pointer(payload, parent_pointer)
        key = int(last) if isinstance(parent, list) else last.replace("~1", "/").replace("~0", "~")
        parent[key] = {"identity_field": contract["identity_fields"][pointer]}
    return payload, values


def candidate_metadata_issues(root: Path, path: Path, contract: dict, identity: dict) -> list[str]:
    if contract.get("format") != "json" or not contract.get("identity_fields"):
        return []
    values = _metadata_values(read_json(path), contract)
    pins = {entry["path"]: entry["sha256"] for entry in identity["code"] + identity["dependencies"]}
    issues = []
    for pointer, value in values.items():
        kind = contract["identity_fields"][pointer]
        if kind in {"file_fingerprint", "file_fingerprints"}:
            for entry in value if isinstance(value, list) else [value]:
                source = root / entry["path"]
                if pins.get(entry["path"]) != entry["sha256"] or not source.is_file() or source.stat().st_size != entry["bytes"]:
                    issues.append(f"candidate lineage does not match executed input/code: {entry['path']}")
        elif kind == "runtime_versions":
            if any(identity["runtime"]["distributions"].get(name) != version for name, version in value.items()):
                issues.append(f"candidate runtime versions do not match child: {pointer}")
        elif value.split()[0] != identity["runtime"]["python"]:
            issues.append("candidate interpreter identity does not match child")
    return issues


def semantic_comparison(original: Path, reproduced: Path, contract: dict | None) -> dict:
    result = {"matches": False, "dimensions": dict.fromkeys(DIMENSIONS, True), "differences": [],
              "difference_count": 0, "negative_case_changes": []}

    def fail(dimension, detail):
        result["dimensions"][dimension] = False
        result["difference_count"] += 1
        if len(result["differences"]) < 25:
            result["differences"].append(f"{dimension}: {detail}")

    if not isinstance(contract, dict) or "format" not in contract or "units" not in contract or "row_keys" not in contract:
        fail("units", "semantic contract pending; format, units and row_keys must be explicit")
        return result
    rtol, atol = contract.get("rtol", 0.0), contract.get("atol", 0.0)
    if not all(isinstance(value, (int, float)) and math.isfinite(value) and value >= 0 for value in (rtol, atol)):
        fail("values", "invalid comparison tolerances")
        return result
    if not isinstance(contract["units"], dict):
        fail("units", "unit declarations must be a mapping, including dimensionless quantities")
        return result
    if not original.is_file() or not reproduced.is_file():
        for dimension in DIMENSIONS:
            fail(dimension, "no complete baseline/candidate pair; equality is unverified")
        return result

    def compare_json(left, right, pointer=""):
        if left is None or right is None:
            if left is not None or right is not None:
                fail("missingness", pointer or "/")
        elif isinstance(left, dict) and isinstance(right, dict):
            if set(left) != set(right):
                fail("row_keys", f"object keys at {pointer or '/'}")
                fail("missingness", f"absent JSON members at {pointer or '/'}")
            for key in sorted(set(left) & set(right)):
                compare_json(left[key], right[key], pointer + "/" + key.replace("~", "~0").replace("/", "~1"))
        elif isinstance(left, list) and isinstance(right, list):
            keys = _json_row_keys(contract, pointer)
            if keys:
                compare_json(_json_index(left, keys, contract, pointer), _json_index(right, keys, contract, pointer), pointer)
            else:
                if len(left) != len(right):
                    fail("row_keys", f"array length at {pointer or '/'}")
                for index, (a, b) in enumerate(zip(left, right)):
                    compare_json(a, b, f"{pointer}/{index}")
        elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
            if not _number_equal(left, right, rtol, atol):
                fail("values", pointer or "/")
        elif type(left) is not type(right) or left != right:
            unit_field = pointer in contract["units"] or any(part in contract.get("unit_field_names", []) for part in pointer.split("/"))
            fail("units" if unit_field else "values", pointer or "/")

    try:
        if contract["format"] == "csv":
            missing = contract.get("missing_values", [""])
            left_header, left = _csv_table(original, contract)
            right_header, right = _csv_table(reproduced, contract)
            result["negative_case_changes"] = _negative_csv_changes(left, right, contract)
            result["row_counts"] = {"baseline": len(left), "candidate": len(right),
                                    "removed": len(set(left) - set(right)), "added": len(set(right) - set(left))}
            result["column_changes"] = {"removed": sorted(set(left_header) - set(right_header)),
                                        "added": sorted(set(right_header) - set(left_header))}
            result["changed_cell_samples"] = []
            if set(left_header) != set(right_header):
                fail("row_keys", "column identities differ")
                fail("missingness", "column presence differs")
                fail("values", "unpaired columns prevent value equality")
                if (set(left_header) ^ set(right_header)) & set(contract["units"]):
                    fail("units", "quantity/unit schema changed")
            optional = set(contract.get("optional_columns", []))
            variants = contract.get("column_variants", [])
            if variants:
                common = set(variants[0]).intersection(*(set(columns) for columns in variants[1:]))
                optional.update(set(contract["units"]) - common)
            if not (set(contract["units"]) - optional) <= set(left_header) & set(right_header):
                fail("units", "declared quantity column is absent")
            if set(left) != set(right):
                fail("row_keys", "row identities differ")
                fail("missingness", "row presence differs")
                fail("values", "unpaired rows prevent value equality")
            for key in sorted(set(left) & set(right)):
                for column in sorted(set(left_header) & set(right_header)):
                    a, b = left[key][column], right[key][column]
                    if (a in missing) != (b in missing):
                        fail("missingness", f"{key}: {column}")
                        if len(result["changed_cell_samples"]) < 25:
                            result["changed_cell_samples"].append({"row_key": list(key), "column": column, "before": a, "after": b, "dimension": "missingness"})
                    elif a not in missing:
                        numeric = column in contract.get("numeric_columns", contract["units"])
                        equal = _number_equal(a, b, rtol, atol) if numeric else a == b
                        if column in contract.get("boolean_columns", []):
                            equal = a.lower() == b.lower()
                        if column in contract.get("json_columns", []):
                            equal = json.loads(a) == json.loads(b)
                        if not equal:
                            dimension = "units" if column in contract.get("unit_columns", []) else "values"
                            fail(dimension, f"{key}: {column}")
                            if len(result["changed_cell_samples"]) < 25:
                                result["changed_cell_samples"].append({"row_key": list(key), "column": column, "before": a, "after": b, "dimension": dimension})
        elif contract["format"] == "json":
            left = _normalize_storage(read_json(original), contract)
            right = _normalize_storage(read_json(reproduced), contract)
            for pointer, expected in contract.get("invariants", {}).items():
                for payload in (left, right):
                    value = json_pointer(payload, pointer)
                    if type(value) is not type(expected) or value != expected:
                        raise ArtifactError(f"scientific scope invariant violated: {pointer}")
            result["negative_case_changes"] = _negative_json_changes(left, right, contract)
            left, left_metadata = _mask_metadata(left, contract)
            right, right_metadata = _mask_metadata(right, contract)
            result["identity_metadata_changes"] = [pointer for pointer in left_metadata
                                                    if left_metadata[pointer] != right_metadata[pointer]]
            compare_json(left, right)
            for pointer, expected in contract["units"].items():
                for payload in (left, right):
                    value = payload
                    for token in pointer.lstrip("/").split("/"):
                        key = token.replace("~1", "/").replace("~0", "~")
                        value = value[int(key)] if isinstance(value, list) else value[key]
                    if value != expected:
                        fail("units", f"{pointer} must equal {expected!r}")
        elif contract["format"] == "npz":
            import numpy as np

            with np.load(original, allow_pickle=False) as left, np.load(reproduced, allow_pickle=False) as right:
                for archive in (left, right):
                    if contract.get("members") and set(archive.files) != set(contract["members"]):
                        raise ArtifactError("unregistered NPZ member schema")
                    for axis in contract.get("row_keys", []):
                        if axis not in archive or archive[axis].ndim != 1 or len(set(archive[axis].tolist())) != len(archive[axis]):
                            raise ArtifactError(f"duplicate or absent NPZ axis identity: {axis}")
                if set(left.files) != set(right.files):
                    fail("row_keys", "NPZ member identities differ")
                for key in sorted(set(left.files) & set(right.files)):
                    a, b = left[key], right[key]
                    if a.shape != b.shape:
                        fail("row_keys", f"NPZ shape: {key}")
                    elif a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
                        if not np.array_equal(np.isnan(a), np.isnan(b)):
                            fail("missingness", key)
                        if not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
                            fail("values", key)
                    elif not np.array_equal(a, b):
                        fail("values", key)
        elif contract["format"] == "bytes":
            if sha256(original) != sha256(reproduced):
                for dimension in DIMENSIONS:
                    fail(dimension, "byte equality required; no format-specific equivalence asserted")
        else:
            raise ArtifactError(f"unsupported semantic format: {contract['format']}")
    except (OSError, ValueError, KeyError, TypeError, IndexError, InvalidOperation, StopIteration, csv.Error) as error:
        fail("row_keys", str(error))
    result["matches"] = all(result["dimensions"].values())
    return result


def _receipt_record(root: Path, registry: dict, run: dict, item: dict) -> tuple[dict, dict]:
    name = item["path"]
    proof = run.get("verification") or {}
    reference = proof.get("artifacts", {}).get(name)
    if reference is not None:
        receipt = registry["receipts"].get(reference.get("receipt_id"))
        if receipt is None or canonical_digest(receipt) != reference["receipt_id"]:
            raise ArtifactError("durable reproduction receipt identity mismatch")
    elif "path" in proof:
        path = artifact_path(root, proof["path"])
        if not path.is_file() or sha256(path) != proof.get("sha256"):
            raise ArtifactError("reproduction receipt identity mismatch")
        receipt = read_json(path)
        reference = {"mode": "equivalent", "accepted_sha256": receipt.get("working_sha256", receipt.get("retained_sha256", {})).get(name)}
    else:
        raise ArtifactError("successful reproduction receipt absent")
    return receipt, reference


def _receipt_issues(root: Path, registry: dict, run: dict, item: dict, current: dict) -> list[str]:
    name = item["path"]
    try:
        receipt, reference = _receipt_record(root, registry, run, item)
    except ArtifactError as error:
        return [str(error)]
    issues = []
    if not _execution_valid(receipt) or receipt.get("execution_complete") is not True:
        issues.append("producer execution or verified child binding failed")
    if receipt.get("run_id") != run["id"] or receipt.get("identity") != current:
        issues.append("input/model/code/parameters/runtime identity changed; numerical effect is untested")
    if _file_digest(root, name) != reference.get("accepted_sha256"):
        issues.append("retained artifact changed after semantic comparison or reviewed candidate not installed")
    if not receipt.get("candidate_validation", {}).get(name, {}).get("matches"):
        issues.append("candidate semantic contract was not verified")
    if reference.get("mode") == "equivalent":
        accepted = reference.get("accepted_sha256")
        witnessed = (accepted == receipt.get("candidate_sha256", {}).get(name) or (
            accepted == receipt.get("working_sha256", {}).get(name)
            and receipt.get("working_comparisons", {}).get(name, {}).get("matches") is True))
        if not accepted or not witnessed:
            issues.append("accepted output hash was not witnessed by this reproduction")
        comparison = receipt.get("comparisons", {}).get(name, {})
        if comparison.get("matches") is not True or comparison.get("dimensions") != dict.fromkeys(DIMENSIONS, True):
            issues.append("values/missingness/units/row_keys were not all verified")
    elif reference.get("mode") == "reviewed_change":
        if reference.get("accepted_sha256") != receipt.get("candidate_sha256", {}).get(name):
            issues.append("accepted output is not the exact reviewed candidate")
        review = registry["reviews"].get(reference.get("review_id"))
        if review is None or canonical_digest(review) != reference["review_id"] or not _review_valid(review, receipt, name):
            issues.append("scientific change has no exact-hash reviewed acceptance")
    else:
        issues.append("unknown reproduction acceptance mode")
    return issues


def _historical_runtime_issues(root: Path, run: dict) -> list[str]:
    runtime = run.get("frozen_runtime") or {}
    if runtime.get("registry"):
        reference = runtime["registry"]
        issues = _pinned_references(root, [reference], "original runtime registry")
        if issues:
            return issues
        recorded = read_json(artifact_path(root, reference["path"]))
        if recorded.get("code_ref") != runtime.get("git_ref") or not recorded.get("environment"):
            return ["original runtime/code reference does not match its registry"]
        for entry in runtime.get("code", []):
            process = subprocess.run(["git", "show", f"{runtime['git_ref']}:{artifact_id(entry['path'])}"], cwd=root, capture_output=True)
            if process.returncode or hashlib.sha256(process.stdout).hexdigest() != entry["sha256"]:
                issues.append(f"original Git producer bytes not verified: {entry['path']}")
        if not runtime.get("code"):
            issues.append("original producer code roster missing")
        return issues
    issues = []
    if not runtime.get("python") or not isinstance(runtime.get("packages"), dict):
        issues.append("original runtime identity is not recorded")
    issues.extend(_pinned_references(root, runtime.get("code"), "original code snapshot"))
    return issues


def _failed_export_issues(root: Path, run: dict, item: dict) -> list[str]:
    issues = _historical_runtime_issues(root, run)
    if _file_digest(root, item["path"]) != item.get("origin_sha256"):
        issues.append("original failed-export bytes changed")
    evidence = item.get("error_evidence") or {}
    issues.extend(_pinned_references(root, [evidence], "failed-export status evidence"))
    if not issues and json_pointer(read_json(artifact_path(root, evidence["path"])), evidence["pointer"]) != "failed_after_prediction_seal":
        issues.append("origin does not document failed_after_prediction_seal")
    return issues


def validation_scope(root: Path) -> tuple[set[str] | None, dict[str, str]]:
    if not (root / REGISTRY_PATH).is_file():
        return None, {}
    registry = load_registry(root)
    tracking = git_tracking(root)
    selected = set(tracking["tracked"])
    selected.update(name for name, item in registry["artifacts"].items() if item["lifecycle"] == "current")
    if not tracking["available"]:
        selected.update(tracking["present"])
    failed = {name: item["scientific_status"] for name, item in registry["artifacts"].items()
              if item.get("scientific_status") == "failed_after_prediction_seal"}
    return selected, failed


def audit_registry(root: Path) -> list[dict]:
    rows = []

    def row(identity, status, actual, detail, expected="explicit lineage and verified reproduction"):
        rows.append({"check": f"artifact provenance: {identity}", "expected": expected,
                     "actual": actual, "status": status, "detail": detail})

    try:
        registry = load_registry(root)
    except (OSError, ValueError, KeyError, TypeError) as error:
        row(REGISTRY_PATH, "FAIL", "investigation_pending", str(error))
        return rows
    tracking = git_tracking(root)
    present = tracking["present"]
    for identity in sorted((present | tracking["tracked"]) - set(registry["artifacts"])):
        local = tracking["available"] and identity not in tracking["tracked"]
        row(identity, "SKIP" if local else "FAIL", "local_decision_pending" if local else "investigation_pending",
            "Not in the authoritative artifact registry; no lifecycle or producer is inferred from the filename. "
            + ("Git does not retain this file; an explicit local decision remains visible, not a current-result failure." if local else "Tracked or unclassified without a Git index: retained ownership remains blocking."))
    identities, private_skips = {}, {}
    for identity, item in sorted(registry["artifacts"].items()):
        try:
            path = artifact_path(root, identity)
            if item.get("local_run"):
                local = item["local_run"]
                issues = _pinned_references(root, [local["metadata"]], "local run metadata")
                if not issues:
                    metadata = read_json(artifact_path(root, local["metadata"]["path"]))
                    issues = [f"local metadata assertion changed: {pointer}" for pointer, value in local.get("assertions", {}).items()
                              if json_pointer(metadata, pointer) != value]
                retained = identity in tracking["tracked"] or not tracking["available"]
                row(identity, "FAIL" if retained else "SKIP", "retained_disposition_required" if retained else "local_exploratory" if not issues else "local_decision_pending",
                    local["reason"] + "; " + "; ".join(issues) + "; producer=" + str(local.get("producer")),
                    "local-only run record, never an exemption for a tracked artifact")
                continue
            if item.get("bundle"):
                if not path.is_file() and item.get("optional"):
                    row(identity, "SKIP", "historical_frozen: original unavailable", f"Validate public counterpart {item['public_counterpart']}; original bytes not checked.")
                    continue
                matches = path.is_file() and sha256(path) == item["sha256"]
                bundle = item["bundle"]
                origin = item.get("origin", identity)
                producers = [binding for binding in bundle.get("original_producers", [])
                             if PurePosixPath(origin).is_relative_to(binding["run_directory"])]
                attribution = "; ".join(binding["entrypoint"] for binding in producers)
                row(identity, "PASS" if matches else "FAIL", "historical_frozen: bytes match" if matches else "frozen identity mismatch",
                    "Compared with the original/public manifest pin, never with current producer source. This is byte integrity, not a new fit or a runtime replay. "
                    + (f"Original producer: {attribution}. " if attribution else "") + bundle.get("reproduction_command", ""),
                    "pinned historical original/public integrity; no runtime execution asserted")
                continue
            run = registry["runs"][item["run"]]
            if item.get("scientific_status") == "failed_after_prediction_seal":
                issues = _failed_export_issues(root, run, item)
                row(identity, "FAIL" if issues else "SKIP", "failed_after_prediction_seal",
                    "; ".join(issues) if issues else "Original failed-export bytes and original code/runtime identity verified; this is not a valid current prediction and is not repaired or restamped.",
                    "preserve the registered historical scientific error and its origin SHA-256")
                continue
            if (run["lifecycle"] != "current" and run.get("retention") != "retained"
                    and tracking["available"] and identity not in tracking["tracked"]):
                row(identity, "SKIP", "local_decision_pending", run["reason"] + "; exact producer=" + reproduction_command(run),
                    "local run decision, not a retained-current provenance obligation")
                continue
            issues = run_metadata_issues(run)
            skipped_private = []
            if not path.is_file():
                issues.append("retained artifact missing")
            if run["lifecycle"] == "current" and not isinstance(item.get("comparison"), dict):
                issues.append("semantic contract pending: explicit values, missingness, units and row keys required")
            proof = run.get("verification") or {}
            has_proof = identity in proof.get("artifacts", {}) or "path" in proof
            if run["status"] == "investigation_pending" and not has_proof:
                issues.extend(run.get("issues", []) or ["successful reproduction and semantic comparison have not been recorded"])
            if run["lifecycle"] != "current":
                if not run.get("evidence"):
                    issues.append("noncurrent disposition lacks authoritative evidence")
                if not item.get("sha256"):
                    issues.append("noncurrent original byte identity unrecorded")
                elif path.is_file() and sha256(path) != item["sha256"]:
                    issues.append("noncurrent original byte identity changed")
                if run["lifecycle"] == "historical_frozen" and not run.get("frozen_runtime"):
                    issues.append("original frozen runtime/code identity absent")
                if run["status"] == "verified":
                    issues.extend(_pinned_references(root, run.get("evidence"), "lifecycle evidence"))
                    for replacement in run.get("superseded_by", []):
                        if not artifact_path(root, replacement).is_file():
                            issues.append(f"superseding artifact unavailable: {replacement}")
                    if run["lifecycle"] == "historical_frozen":
                        issues.extend(_historical_runtime_issues(root, run))
            elif not issues:
                rules, _ = _private_spec(root, run)
                if rules:
                    receipt, _ = _receipt_record(root, registry, run, item)
                    if not _execution_valid(receipt) or receipt.get("execution_complete") is not True:
                        raise ArtifactError("private receipt lacks successful, bound prior child execution")
                    key = (run["id"], canonical_digest(receipt["identity"]))
                    if key not in identities:
                        identities[key] = _receipt_only_identity(root, run, receipt["identity"])
                    current, skipped_private = identities[key]
                else:
                    if run["id"] not in identities:
                        identities[run["id"]] = run_identity(root, run)
                    current = identities[run["id"]]
                issues.extend(_receipt_issues(root, registry, run, item, current))
            detail = "; ".join(issues)
            if skipped_private and not issues:
                detail += "; durable prior private execution, output/code/public-input/runtime hashes verified; private bytes not reverified; not a public replay"
                for entry in skipped_private:
                    private_skips[(run["id"], entry["path"], entry["sha256"])] = entry
            producer = run.get("producer") or {}
            detail += f"; producer={producer.get('path', 'unresolved')}::{producer.get('callable', 'unresolved')}; reproduce={reproduction_command(run)}"
            if run["lifecycle"] == "current" and run.get("reproduction"):
                detail += f"; isolated check: python3 -B scripts/regenerate_current_artifacts.py --run {shlex.quote(run['id'])} --output-dir <new-temporary-directory>"
            actual = "investigation_pending" if issues else "verified_prior_private_execution" if skipped_private else run["lifecycle"]
            row(identity, "FAIL" if issues else "PASS", actual, detail.lstrip("; "))
        except (OSError, ValueError, KeyError, TypeError, SyntaxError) as error:
            row(identity, "FAIL", "investigation_pending", str(error))
    for (run_id, name, digest), _ in sorted(private_skips.items()):
        rows.append({"check": f"private input verification: {run_id}: {name}", "status": "SKIP",
                     "actual": "private_bytes_not_reverified", "expected": "explicit opt-in to reverify private input bytes",
                     "detail": f"Prior execution pinned {name} sha256={digest}; private bytes are unavailable to this non-opted-in audit and were NOT reverified. "
                               "This is receipt integrity, not a public replay. Re-execute with --private-inputs-from-env <env-file> to check the approved sources."})
    return rows


def plan_current_runs(root: Path) -> list[dict]:
    registry = load_registry(root)
    current = {name: run for name, run in registry["runs"].items() if run["lifecycle"] == "current"}
    graph = {name: {registry["artifacts"][entry["path"]]["run"] for entry in run["dependencies"]
                    if entry["path"] in registry["artifacts"] and registry["artifacts"][entry["path"]].get("run") in current}
             for name, run in current.items()}
    ordered, finished = [], set()
    while len(finished) < len(current):
        ready = sorted(name for name in current if name not in finished and graph[name] <= finished)
        if not ready:
            raise ArtifactError("registered current artifact dependency cycle")
        for name in ready:
            run = current[name]
            issues = run_metadata_issues(run)
            try:
                count = len(dependency_sources(root, run))
            except (OSError, ValueError) as error:
                count = None
                issues.append(str(error))
            ordered.append({"run_id": name, "upstream_runs": sorted(graph[name]), "producer": run["producer"],
                            "artifacts": [item["path"] for item in run["artifacts"]], "parameters": run["parameters"],
                            "model_identity": run["model_identity"], "dependency_count": count,
                            "input_issues": issues, "producer_recipe": run["reproduction"],
                            "staging_command": f"{shlex.quote(sys.executable)} -B scripts/regenerate_current_artifacts.py --run {shlex.quote(name)} --output-dir <new-temporary-directory>/{name}"})
            finished.add(name)
    return ordered


def _execution_argv(workspace: Path, run: dict) -> list[str]:
    replacements = {"python": sys.executable, "workspace": str(workspace), "output_dir": str(workspace / "outputs")}
    argv = [argument.format(**replacements) for argument in run["reproduction"]["argv"]]
    if argv[0] != sys.executable:
        raise ArtifactError("reproduction must use {python} with an explicit argv, not a shell")
    arguments = argv[2:] if len(argv) > 1 and argv[1] == "-B" else argv[1:]
    if not arguments or (workspace / arguments[0]).resolve() != artifact_path(workspace, run["producer"]["path"]):
        raise ArtifactError("reproduction argv does not execute the registered producer")
    for argument in argv[1:]:
        value = argument.partition("=")[2] if argument.startswith("--") and "=" in argument else argument
        path = Path(value)
        if (path.is_absolute() or ".." in path.parts) and not (workspace / path).resolve().is_relative_to(workspace):
            raise ArtifactError("reproduction argument escapes the isolated workspace")
    for value in run["reproduction"].get("environment", {}).values():
        path = Path(value.format(**replacements))
        if (path.is_absolute() or ".." in path.parts) and not (workspace / path).resolve().is_relative_to(workspace):
            raise ArtifactError("reproduction environment escapes the isolated workspace")
    return arguments


def _check_private_stage(directory: Path) -> None:
    directory = _no_symlinks(directory, "private staging directory")
    if not any(directory != base and directory.is_relative_to(base) for base in _PRIVATE_TEMPORARY_ROOTS):
        raise ArtifactError("private inputs may be staged only in an external temporary directory")
    if directory.exists() and (not directory.is_dir() or stat.S_IMODE(directory.stat().st_mode) != 0o700):
        raise ArtifactError("private staging directory must have restrictive 0700 permissions")


def _private_mkdir(directory: Path) -> None:
    missing, parent = [], directory
    while not parent.exists():
        missing.append(parent)
        parent = parent.parent
    _no_symlinks(parent, "private staging directory")
    for path in reversed(missing):
        path.mkdir(mode=0o700)


def _restrict_private_tree(directory: Path) -> None:
    for path in [directory, *directory.rglob("*")]:
        if path.is_symlink():
            raise ArtifactError("private staging may not contain symlinks")
        if not (path.is_file() or path.is_dir()):
            raise ArtifactError("private staging may contain only regular files and directories")
        path.chmod(0o700 if path.is_dir() else 0o600)


def _copy_private_input(source: Path, destination: Path, name: str) -> None:
    _private_mkdir(destination.parent)
    with _private_open(source, name) as reader:
        descriptor = os.open(destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(descriptor, "wb") as writer:
            shutil.copyfileobj(reader, writer)
    # No copystat/copy2: private instrument metadata and original permissions are not propagated.


def _comparison_record(result: dict, private: bool) -> dict:
    if not private:
        return result
    # Private producers may put raw workbook properties in warnings or malformed cells.
    # Persist structural verdicts/counts and exact detail hashes, never those free-text values.
    record = {key: result[key] for key in ("matches", "dimensions", "difference_count", "row_counts") if key in result}
    record["differences"] = [f"{key}: private comparison differs; inspect restricted staged outputs"
                             for key, matches in result["dimensions"].items() if not matches]
    record["details_sha256"] = canonical_digest(result)
    record["negative_case_changes"] = [{"change_sha256": canonical_digest(change), "row_retained": change.get("row_retained")}
                                        for change in result.get("negative_case_changes", [])]
    if "column_changes" in result:
        record["column_change_counts"] = {key: len(value) for key, value in result["column_changes"].items()}
    if "metadata_issues" in result:
        record["metadata_issues"] = ["private candidate lineage metadata invalid"]
    return record


def execute_plan(plan_path: Path) -> int:
    plan = read_json(plan_path)
    workspace = (plan_path.parent / "workspace").resolve()
    run, expected = plan["run"], plan["identity"]
    private = "private_inputs" in expected
    if private:
        os.umask(0o077)
        _check_private_stage(plan_path.parent)
        os.environ["TMPDIR"] = str(workspace / ".runtime/tmp")
        tempfile.tempdir = os.environ["TMPDIR"]
    expected_digest = canonical_digest(expected)
    before = run_identity(workspace, run, staged=True)
    reads, writes, outside_reads = set(), set(), set()
    runtime_roots = [Path(sys.prefix).resolve(), Path(sys.base_prefix).resolve(),
                     *[Path(value).resolve() for value in site.getsitepackages()],
                     Path(site.getusersitepackages()).resolve(),
                     *[Path(value).resolve() for value in ("/System/Library", "/Library/Fonts", "/usr/share", "/etc")]]
    active = [False]
    declared = {entry["path"] for entry in expected["code"] + expected["dependencies"]}
    immutable = declared | {PRIVATE_POLICY_PATH, "src/ystwin/artifacts.py"}
    if private:
        runtime_roots = [base for base in runtime_roots if base != Path("/etc")]

    def writable(value):
        if not isinstance(value, (str, bytes, os.PathLike)):
            raise ArtifactError("private producer may not mutate an unbound file descriptor")
        path = Path(os.fsdecode(value)).absolute()
        if not path.resolve().is_relative_to(workspace):
            raise ArtifactError("producer attempted to write outside its reproduction workspace")
        relative = path.relative_to(workspace).as_posix()
        if any(PurePosixPath(name).is_relative_to(relative) for name in immutable):
            raise ArtifactError("private producer attempted to modify staged input/code")

    def observe(event, args):
        if not active[0]:
            return
        if private:
            if event in {"socket.connect", "socket.connect_ex", "socket.getaddrinfo", "socket.sendto",
                         "subprocess.Popen", "os.system", "os.exec", "os.posix_spawn", "os.fork", "os.forkpty",
                         "os.symlink", "os.link"}:
                raise ArtifactError("private reproduction forbids network access, subprocesses and links")
            mutations = {"os.remove": (0,), "os.rmdir": (0,), "os.mkdir": (0,), "os.rename": (0, 1),
                         "os.chmod": (0,), "os.chown": (0,), "os.truncate": (0,), "os.utime": (0,)}
            if event in mutations:
                directory_arguments = {"os.remove": (1,), "os.rmdir": (1,), "os.mkdir": (2,), "os.rename": (2, 3),
                                       "os.chmod": (2,), "os.chown": (3,), "os.utime": (3,)}
                if any(args[index] not in (None, -1) for index in directory_arguments.get(event, ())):
                    raise ArtifactError("private producer may not mutate through directory descriptors")
                for index in mutations[event]:
                    writable(args[index])
        if event != "open" or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        path = Path(os.fsdecode(args[0])).absolute()
        mode, flags = args[1], args[2]
        writing = (isinstance(mode, str) and any(letter in mode for letter in "wax+")) or (
            isinstance(flags, int) and bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT)))
        if writing and not path.resolve().is_relative_to(workspace):
            raise ArtifactError("producer attempted to write outside its reproduction workspace")
        if private and writing:
            writable(path)
        if path.is_relative_to(workspace):
            relative = path.relative_to(workspace).as_posix()
            if writing:
                writes.add(relative)
            elif path.is_file() and relative not in writes and not relative.startswith(".runtime/"):
                reads.add(relative)
        elif not writing and (private or path.is_file()) and path != Path(os.devnull):
            resolved = path.resolve()
            if not any(resolved.is_relative_to(base) for base in runtime_roots):
                outside_reads.add(hashlib.sha256(str(resolved).encode()).hexdigest())
                raise ArtifactError("producer attempted to read an input outside the declared workspace; register and stage that dependency")

    sys.addaudithook(observe)
    code = 70
    failure = ""
    if before == expected:
        sys.argv = _execution_argv(workspace, run)
        sys.path.insert(0, str(artifact_path(workspace, run["producer"]["path"]).parent))
        active[0] = True
        try:
            runpy.run_path(str(artifact_path(workspace, run["producer"]["path"])), run_name="__main__")
            code = 0
        except SystemExit as error:
            code = 0 if error.code is None else int(error.code) if isinstance(error.code, int) else 1
            if not isinstance(error.code, (int, type(None))):
                print(error.code, file=sys.stderr)
        except BaseException as error:
            code = 130 if isinstance(error, KeyboardInterrupt) else 1
            failure = type(error).__name__
            traceback.print_exc()
        finally:
            active[0] = False
    else:
        failure = "staged input/code/runtime binding differs before execution"
    try:
        after_digest = canonical_digest(run_identity(workspace, run, staged=True))
    except (OSError, ValueError, KeyError, SyntaxError) as error:
        after_digest = None
        failure = f"post-execution binding unavailable: {type(error).__name__}"
    unknown_reads = sorted(reads - declared - {"src/ystwin/artifacts.py"})
    known_paths = immutable | {item["path"] for item in run["artifacts"]}

    def trace_paths(values):
        return sorted(value if not private or value in known_paths else "unregistered/" + hashlib.sha256(value.encode()).hexdigest()
                      for value in values)

    ledger = {"schema_version": 2, "run_id": run["id"], "exit_code": code,
              "before_binding_sha256": canonical_digest(before), "after_binding_sha256": after_digest,
              "expected_binding_sha256": expected_digest, "runtime": before["runtime"],
              "read_paths": trace_paths(reads), "written_paths": trace_paths(writes), "undeclared_read_paths": trace_paths(unknown_reads),
              "undeclared_external_reads": sorted(outside_reads), "failure": failure}
    (plan_path.parent / "child_execution.json").write_text(json.dumps(ledger, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return code


def _file_digest(root: Path, identity: str):
    path = artifact_path(root, identity)
    return sha256(path) if path.is_file() else None


def _baseline_bytes(root: Path, registry: dict, item: dict, run: dict):
    name = item["path"]
    accepted = (run.get("verification") or {}).get("artifacts", {}).get(name)
    if accepted:
        previous = registry["receipts"].get(accepted.get("receipt_id"))
        if previous is None or _receipt_issues(root, registry, run, item, previous["identity"]):
            raise ArtifactError(f"previous acceptance lacks valid output-bound evidence: {name}")
        path = artifact_path(root, name)
        if not path.is_file() or sha256(path) != accepted["accepted_sha256"]:
            raise ArtifactError(f"previously accepted output changed without review: {name}")
        return path.read_bytes(), {"kind": "accepted_receipt", "receipt_id": accepted["receipt_id"]}
    if item.get("origin_sha256"):
        path = artifact_path(root, name)
        if not path.is_file() or sha256(path) != item["origin_sha256"]:
            raise ArtifactError(f"registered original bytes unavailable: {name}")
        return path.read_bytes(), {"kind": "registered_original", "sha256": item["origin_sha256"]}
    reference = registry["document"].get("baseline_ref")
    if reference:
        if len(reference) != 40 or any(c not in "0123456789abcdef" for c in reference):
            raise ArtifactError("baseline_ref must be an immutable full Git commit")
        process = subprocess.run(["git", "show", f"{reference}:{name}"], cwd=root, capture_output=True)
        if process.returncode == 0:
            return process.stdout, {"kind": "git", "commit": reference, "path": name}
        if item.get("planned"):
            return None, {"kind": "new_artifact"}
        raise ArtifactError(f"no registered original in {reference}: {name}; an explicit origin identity is required")
    path = artifact_path(root, name)
    return (path.read_bytes(), {"kind": "working_tree"}) if path.is_file() else (None, {"kind": "new_artifact"})


def _execution_valid(receipt: dict) -> bool:
    execution = receipt.get("execution", {})
    digest = canonical_digest(receipt.get("identity"))
    return (receipt.get("schema_version") == 2 and type(receipt.get("exit_code")) is int
            and receipt["exit_code"] == 0 and receipt.get("inputs_unchanged") is True
            and execution.get("schema_version") == 2 and type(execution.get("exit_code")) is int
            and execution["exit_code"] == 0 and execution.get("run_id") == receipt.get("run_id")
            and execution.get("before_binding_sha256") == digest
            and execution.get("after_binding_sha256") == digest
            and execution.get("expected_binding_sha256") == digest
            and execution.get("runtime") == receipt["identity"]["runtime"]
            and receipt.get("execution_sha256") == canonical_digest(execution)
            and execution.get("undeclared_read_paths") == [] and execution.get("undeclared_external_reads") == []
            and not execution.get("failure"))


def reproduce_run(root: Path, run_id: str, output_dir: Path, *, timeout: float = 600.0,
                  private_inputs_from_env: Path | None = None) -> dict:
    root = root.resolve()
    if private_inputs_from_env is not None:
        _check_private_stage(Path(output_dir).absolute())
    output_dir = output_dir.resolve()
    if output_dir.is_relative_to(root) or root.is_relative_to(output_dir):
        raise ArtifactError("reproduction requires a new temporary directory outside the repository")
    registry = load_registry(root)
    if run_id not in registry["runs"]:
        raise ArtifactError(f"unknown run: {run_id}")
    run = registry["runs"][run_id]
    if run["lifecycle"] != "current":
        raise ArtifactError("only current artifacts may be regenerated; frozen/local/superseded evidence is never restamped")
    for item in run["artifacts"]:
        if not isinstance(item.get("comparison"), dict):
            raise ArtifactError(f"semantic contract must be completed before reproduction: {item['path']}")
    identity = run_identity(root, run, private_inputs_from_env=private_inputs_from_env)
    private = "private_inputs" in identity
    workspace = output_dir / "workspace"
    _execution_argv(workspace, run)
    baselines = {item["path"]: _baseline_bytes(root, registry, item, run) for item in run["artifacts"]}
    working = {item["path"]: _file_digest(root, item["path"]) for item in run["artifacts"]}
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ArtifactError("reproduction output directory must be new or empty")
    if private:
        _private_mkdir(output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
    workspace.mkdir(mode=0o700 if private else 0o777)
    for folder in ("src", "scripts"):
        if (root / folder).is_dir():
            shutil.copytree(root / folder, workspace / folder, symlinks=private,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    verifier = workspace / "src/ystwin/artifacts.py"
    verifier.parent.mkdir(parents=True, exist_ok=True)
    verifier.write_bytes(Path(__file__).read_bytes())
    if (root / "pyproject.toml").is_file():
        shutil.copy2(root / "pyproject.toml", workspace / "pyproject.toml")
    if private:
        _restrict_private_tree(output_dir)
        policy_destination = artifact_path(workspace, PRIVATE_POLICY_PATH)
        _copy_private_input(artifact_path(root, PRIVATE_POLICY_PATH), policy_destination, PRIVATE_POLICY_PATH)
    private_names = {entry["path"] for entry in identity["dependencies"] if entry.get("private") is True}
    sources = dependency_sources(root, run, private_inputs_from_env=private_inputs_from_env)
    if set(sources) != {entry["path"] for entry in identity["dependencies"]}:
        raise ArtifactError("input membership changed before staging; no changed corpus is authorized by the earlier identity")
    for name, source in sources.items():
        destination = artifact_path(workspace, name)
        if name in private_names:
            _copy_private_input(source, destination, name)
        else:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    (workspace / "outputs").mkdir(exist_ok=True)
    for name, (content, _) in baselines.items():
        artifact_path(workspace, name).parent.mkdir(parents=True, exist_ok=True)
        if content is not None:
            baseline = artifact_path(output_dir / "baseline", name)
            baseline.parent.mkdir(parents=True, exist_ok=True)
            baseline.write_bytes(content)
    for directory in ("tmp", "matplotlib", "cache", "home"):
        (workspace / ".runtime" / directory).mkdir(parents=True, exist_ok=True)
    plan_path = output_dir / "execution_plan.json"
    plan_path.write_text(json.dumps({"run": run, "identity": identity}, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    if private:
        _restrict_private_tree(output_dir)
        environment = {key: os.environ[key] for key in ("PATH", "LANG", "LC_ALL", "TZ", "SYSTEMROOT") if key in os.environ}
        environment.update({"HOME": str(workspace / ".runtime/home"), "PYTHONNOUSERSITE": "1"})
    else:
        environment = {key: value for key, value in os.environ.items() if not key.startswith("YSTWIN_") and key not in {"PYTHONPATH", "PYTHONHOME"}}
    environment.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": str(workspace / "src"),
                        "YSTWIN_OUTPUTS": str(workspace / "outputs"), "TMPDIR": str(workspace / ".runtime/tmp"),
                        "MPLCONFIGDIR": str(workspace / ".runtime/matplotlib"), "XDG_CACHE_HOME": str(workspace / ".runtime/cache")})
    replacements = {"python": sys.executable, "workspace": str(workspace), "output_dir": str(workspace / "outputs")}
    environment.update({key: value.format(**replacements) for key, value in run["reproduction"].get("environment", {}).items()})
    if private:
        environment["TMPDIR"] = tempfile.gettempdir()  # child validates launch scope before redirecting this
    bootstrap = "import importlib.util,sys,pathlib; s=importlib.util.spec_from_file_location('_artifact_child',sys.argv[1]); m=importlib.util.module_from_spec(s); sys.modules[s.name]=m; s.loader.exec_module(m); raise SystemExit(m.execute_plan(pathlib.Path(sys.argv[2])))"
    try:
        completed = subprocess.run([sys.executable, "-B", "-c", bootstrap, str(verifier), str(plan_path)],
                                   cwd=workspace, env=environment, capture_output=True, text=True, timeout=timeout)
        code, stdout, stderr = completed.returncode, completed.stdout, completed.stderr
    except subprocess.TimeoutExpired as error:
        code, stdout, stderr = 124, error.stdout or "", error.stderr or ""
        stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
    execution_path = output_dir / "child_execution.json"
    execution = read_json(execution_path) if execution_path.is_file() else {}
    comparisons, working_comparisons, validation, candidate_hashes = {}, {}, {}, {}
    for item in run["artifacts"]:
        name, contract = item["path"], item["comparison"]
        candidate = artifact_path(workspace, name)
        comparisons[name] = semantic_comparison(artifact_path(output_dir / "baseline", name), candidate, contract)
        working_comparisons[name] = semantic_comparison(artifact_path(root, name), candidate, contract)
        validation[name] = semantic_comparison(candidate, candidate, contract)
        if candidate.is_file() and validation[name]["matches"]:
            metadata_issues = candidate_metadata_issues(workspace, candidate, contract, identity)
            if metadata_issues:
                validation[name]["matches"] = False
                validation[name]["metadata_issues"] = metadata_issues
        candidate_hashes[name] = sha256(candidate) if candidate.is_file() else None
        comparisons[name] = _comparison_record(comparisons[name], private)
        working_comparisons[name] = _comparison_record(working_comparisons[name], private)
        validation[name] = _comparison_record(validation[name], private)
    try:
        unchanged = identity == run_identity(root, run, private_inputs_from_env=private_inputs_from_env) and working == {
            item["path"]: _file_digest(root, item["path"]) for item in run["artifacts"]}
    except (OSError, ValueError, KeyError, SyntaxError):
        unchanged = False
    receipt = {"schema_version": 2, "run_id": run_id, "identity": identity, "exit_code": code,
               "baseline_policy": registry["document"].get("baseline_ref"),
               "execution": execution, "execution_sha256": canonical_digest(execution),
               "working_sha256": working, "candidate_sha256": candidate_hashes,
               "retained_sha256": {name: hashlib.sha256(content).hexdigest() if content is not None else None for name, (content, _) in baselines.items()},
               "baseline_sources": {name: source for name, (_, source) in baselines.items()},
               "comparisons": comparisons, "working_comparisons": working_comparisons,
               "candidate_validation": validation, "inputs_unchanged": unchanged,
               "outputs_not_selected": sorted(name for name in execution.get("written_paths", [])
                                              if name.startswith("outputs/") and name not in candidate_hashes),
               "scope": "Staged execution and comparison only; no retained outputs or provenance stamps written"}
    receipt["execution_complete"] = _execution_valid(receipt) and all(result["matches"] for result in validation.values())
    receipt["success"] = receipt["execution_complete"] and all(result["matches"] for result in comparisons.values())
    if private:
        receipt["scope"] += "; explicit private-input execution, not a public replay; private stdout/stderr and raw comparison details suppressed"
        stdout = "Private producer stdout suppressed; consult the structured execution receipt.\n"
        stderr = "Private producer stderr suppressed; exception types and exit status are in the execution receipt.\n"
    (output_dir / "stdout.txt").write_text(stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(stderr, encoding="utf-8")
    (output_dir / "reproduction.json").write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    if private:
        _restrict_private_tree(output_dir)
    return receipt


def validate_reproduction(root: Path, receipt_path: Path, *,
                          private_inputs_from_env: Path | None = None) -> tuple[dict, dict, dict]:
    receipt_path = receipt_path.resolve(strict=True)
    receipt = read_json(receipt_path)
    registry = load_registry(root)
    run = registry["runs"].get(receipt.get("run_id"))
    if run is None or run["lifecycle"] != "current":
        raise ArtifactError("receipt does not identify a registered current run")
    if not _execution_valid(receipt) or receipt.get("execution_complete") is not True:
        raise ArtifactError("receipt lacks successful, bound child execution; reproduce with the current verifier")
    if receipt["identity"] != run_identity(root, run, private_inputs_from_env=private_inputs_from_env):
        raise ArtifactError("stale reproduction receipt: current code/input/model/parameters/runtime/contract binding changed")
    private = "private_inputs" in receipt["identity"]
    if private:
        _check_private_stage(receipt_path.parent)
    if receipt.get("baseline_policy") != registry["document"].get("baseline_ref"):
        raise ArtifactError("stale reproduction receipt: authoritative baseline policy changed")
    execution = read_json(receipt_path.parent / "child_execution.json")
    if execution != receipt["execution"] or canonical_digest(execution) != receipt["execution_sha256"]:
        raise ArtifactError("child execution record changed after reproduction")
    workspace = receipt_path.parent / "workspace"
    if workspace.resolve().is_relative_to(root.resolve()):
        raise ArtifactError("receipt workspace must be outside the repository")
    if run_identity(workspace, run, staged=True) != receipt["identity"]:
        raise ArtifactError("staged producer or input identities changed after reproduction")
    for item in run["artifacts"]:
        name, contract = item["path"], item["comparison"]
        if _file_digest(root, name) != receipt["working_sha256"][name]:
            raise ArtifactError(f"stale reproduction receipt: retained output changed after child execution: {name}")
        baseline = artifact_path(receipt_path.parent / "baseline", name)
        candidate = artifact_path(workspace, name)
        if (sha256(baseline) if baseline.is_file() else None) != receipt["retained_sha256"][name]:
            raise ArtifactError(f"recorded baseline changed: {name}")
        if not candidate.is_file() or sha256(candidate) != receipt["candidate_sha256"][name]:
            raise ArtifactError(f"staged candidate changed: {name}")
        if _comparison_record(semantic_comparison(baseline, candidate, contract), private) != receipt["comparisons"][name]:
            raise ArtifactError(f"recorded semantic comparison changed: {name}")
        if _comparison_record(semantic_comparison(artifact_path(root, name), candidate, contract), private) != receipt.get("working_comparisons", {}).get(name):
            raise ArtifactError(f"working-output comparison is not verified: {name}")
        if not semantic_comparison(candidate, candidate, contract)["matches"] or candidate_metadata_issues(workspace, candidate, contract, receipt["identity"]):
            raise ArtifactError(f"candidate contract or lineage is invalid: {name}")
    return registry, run, receipt


def review_template(root: Path, receipt_path: Path, *, private_inputs_from_env: Path | None = None) -> dict:
    _, run, receipt = validate_reproduction(root, receipt_path, private_inputs_from_env=private_inputs_from_env)
    entries = {}
    for item in run["artifacts"]:
        name = item["path"]
        comparison = receipt["comparisons"][name]
        if comparison["matches"]:
            continue
        entries[name] = {"baseline_sha256": receipt["retained_sha256"][name],
                         "candidate_sha256": receipt["candidate_sha256"][name],
                         "comparison_sha256": canonical_digest(comparison),
                         "observed_changed_dimensions": [key for key, equal in comparison["dimensions"].items() if not equal],
                         "accepted_dimensions": [], "reason": "",
                         "required_negative_case_changes_sha256": canonical_digest(comparison.get("negative_case_changes", [])),
                         "approved_negative_case_changes_sha256": None}
    return {"schema_version": 1, "approved": False, "reviewer": "", "reason": "",
            "reproduction_sha256": canonical_digest(receipt), "binding_sha256": canonical_digest(receipt["identity"]),
            "artifacts": entries}


def _review_valid(review: dict, receipt: dict, name: str) -> bool:
    if (review.get("schema_version") != 1 or review.get("approved") is not True
            or not review.get("reviewer") or not review.get("reason")
            or review.get("reproduction_sha256") != canonical_digest(receipt)
            or review.get("binding_sha256") != canonical_digest(receipt["identity"])):
        return False
    approval = review.get("artifacts", {}).get(name, {})
    comparison = receipt["comparisons"][name]
    changed = {key for key, equal in comparison["dimensions"].items() if not equal}
    return (approval.get("baseline_sha256") == receipt["retained_sha256"][name]
            and approval.get("candidate_sha256") == receipt["candidate_sha256"][name]
            and approval.get("comparison_sha256") == canonical_digest(comparison)
            and bool(approval.get("reason")) and changed <= set(approval.get("accepted_dimensions", []))
            and (not comparison.get("negative_case_changes") or approval.get("approved_negative_case_changes_sha256") == canonical_digest(comparison["negative_case_changes"])))


def adopt_reproduction(root: Path, receipt_path: Path, *, artifact_ids: list[str] | None = None,
                       accept_scientific_changes: bool = False, review_path: Path | None = None,
                       install: bool = False, private_inputs_from_env: Path | None = None) -> dict:
    root = root.resolve()
    registry_path = artifact_path(root, REGISTRY_PATH)
    original_registry = registry_path.read_bytes()
    registry, run, receipt = validate_reproduction(root, receipt_path, private_inputs_from_env=private_inputs_from_env)
    selected = set(artifact_ids or [item["path"] for item in run["artifacts"]])
    if not selected or not selected <= set(receipt["comparisons"]):
        raise ArtifactError("adoption must name exact artifacts from this reproduction")
    review = read_json(review_path) if review_path is not None else None
    workspace = receipt_path.resolve().parent / "workspace"
    decisions, contents = {}, {}
    receipt_id = canonical_digest(receipt)
    for name in sorted(selected):
        comparison = receipt["comparisons"][name]
        changed = not comparison["matches"]
        if changed and (not accept_scientific_changes or review is None or not _review_valid(review, receipt, name)):
            raise ArtifactError(f"scientific changes require explicit exact-hash reviewed acceptance, including negative cases: {name}")
        candidate = artifact_path(workspace, name)
        current = artifact_path(root, name)
        item = next(item for item in run["artifacts"] if item["path"] == name)
        if not changed and not install and not semantic_comparison(current, candidate, item["comparison"])["matches"]:
            raise ArtifactError(f"working output is not equivalent to the accepted baseline; explicit installation required: {name}")
        accepted = receipt["candidate_sha256"][name] if install or changed else receipt["working_sha256"][name]
        decisions[name] = {"receipt_id": receipt_id, "mode": "reviewed_change" if changed else "equivalent",
                           "accepted_sha256": accepted, "installed": install or accepted == _file_digest(root, name)}
        if changed:
            decisions[name]["review_id"] = canonical_digest(review)
        if install:
            contents[name] = candidate.read_bytes()
    document = registry["document"]
    document.setdefault("receipts", {})[receipt_id] = receipt
    if review is not None:
        document.setdefault("reviews", {})[canonical_digest(review)] = review
    original_run = next(entry for entry in document["runs"] if entry["id"] == run["id"])
    verification = original_run.setdefault("verification", {"artifacts": {}})
    verification.setdefault("artifacts", {}).update(decisions)
    original_run["status"] = "verified" if all(item["path"] in verification["artifacts"] for item in run["artifacts"]) else "investigation_pending"
    validate_reproduction(root, receipt_path, private_inputs_from_env=private_inputs_from_env)
    if registry_path.read_bytes() != original_registry:
        raise ArtifactError("artifact registry changed during adoption; retry against the new registry")
    for name, content in contents.items():
        target = artifact_path(root, name)
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(target.name + ".adopting-" + uuid.uuid4().hex)
        temporary.write_bytes(content)
        os.replace(temporary, target)
    temporary = registry_path.with_name(registry_path.name + ".adopting-" + uuid.uuid4().hex)
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, registry_path)
    return {"run_id": run["id"], "receipt_id": receipt_id, "decisions": decisions,
            "outputs_replaced": sorted(contents), "scope": "verified receipt adoption, not blind restamping"}
