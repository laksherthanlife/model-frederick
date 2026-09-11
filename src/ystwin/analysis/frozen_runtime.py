from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tempfile

from .portable_evidence import _relative
from .portable_replay import (
    CHECKPOINT_PATH,
    CHECKPOINT_SHA256,
    DEFAULT_MANIFEST_PATH,
    DEFAULT_MANIFEST_SHA256,
    FROZEN_CODE_PATHS,
    MODEL_SHA256,
    REPLAY_KIND,
    _VARIANTS,
    _load_bundle,
    _native_checkpoint_manifest,
    _relative_file,
    _same,
    _verify_native_checkpoint,
)
from .training_freeze import _parse_json, _sha256


FROZEN_CODE_REF = "84561f51502b5b4bf66c052dc88c853038917934"
RUNTIME_REGISTRY_PATH = "data/frozen_runtime.json"
RUNTIME_REGISTRY_SHA256 = "8e340b7d947d0a89b52d448fba9c20c6899eaf2d568e7a04dfc87426d1c90c0b"
HISTORICAL_AUDIT_KIND = "historical_source_bound_biology_learning_audit"


@dataclass(frozen=True)
class FrozenRuntimeSnapshot:
    root: Path
    registry: dict
    identities: dict
    supplied_data: list
    operation: str = "replay"


class HistoricalReplayError(RuntimeError):
    def __init__(self, reason, *, returncode=None, stdout="", stderr=""):
        super().__init__(reason)
        self.returncode = returncode
        self.report = {
            "schema_version": 1, "kind": "historical_portable_replay_failure",
            "verified": False, "code_ref": FROZEN_CODE_REF,
            "child_exit_code": returncode, "reason": reason,
            "stdout": stdout, "stderr": stderr,
            "new_fitting": False, "biological_validation": False,
        }


def load_frozen_runtime(root):
    root = Path(root).resolve(strict=True)
    data = _relative_file(root, RUNTIME_REGISTRY_PATH).read_bytes()
    if hashlib.sha256(data).hexdigest() != RUNTIME_REGISTRY_SHA256:
        raise ValueError("frozen runtime registry SHA-256 mismatch")
    registry = _parse_json(data)
    _same(registry["code_ref"], FROZEN_CODE_REF, "registered frozen Git commit")
    _same(registry["manifest"], {"path": DEFAULT_MANIFEST_PATH, "sha256": DEFAULT_MANIFEST_SHA256},
          "registered frozen manifest identity")
    _same(registry["checkpoint"], {
        "path": CHECKPOINT_PATH, "sha256": CHECKPOINT_SHA256, "model_sha256": MODEL_SHA256,
    }, "registered frozen checkpoint identity")
    return registry


def _identity(entry):
    path = _relative(entry["path"])
    if any(ord(character) < 32 for character in path):
        raise ValueError("frozen dependencies require canonical repository-relative paths")
    if entry.get("resolved_path", path) != path:
        raise ValueError("frozen dependency aliases are not allowed")
    size = entry["size_bytes"]
    if type(size) is not int or size <= 0:
        raise ValueError("invalid frozen dependency byte count")
    return {"path": path, "sha256": _sha256(entry["sha256"]), "size_bytes": size}


def _verified_content(data, identity, basis):
    if (hashlib.sha256(data).hexdigest() != identity["sha256"]
            or len(data) != identity["size_bytes"]):
        raise ValueError(f"{basis} SHA-256 mismatch for frozen dependency: {identity['path']}")
    return data


def _git(repository, *arguments, input=None):
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    environment.update(GIT_NO_REPLACE_OBJECTS="1", GIT_NO_LAZY_FETCH="1", GIT_TERMINAL_PROMPT="0")
    process = subprocess.run(
        ["git", "--no-replace-objects", "-C", str(repository), *arguments], input=input,
        env=environment, capture_output=True, timeout=60,
    )
    if process.returncode:
        raise ValueError(
            "frozen Git snapshot unavailable; provision the registered commit and its public source blobs "
            "in --source-repository, without substituting current code: "
            + process.stderr.decode("utf-8", errors="replace").strip()
        )
    return process.stdout


def _git_blobs(repository, code_ref, paths):
    if code_ref != FROZEN_CODE_REF:
        raise ValueError("historical replay requires the registered full frozen Git commit")
    top = Path(os.fsdecode(_git(repository, "rev-parse", "--show-toplevel").strip())).resolve(strict=True)
    if top != repository:
        raise ValueError("source repository must name its Git worktree root")
    if _git(repository, "cat-file", "-t", code_ref).strip() != b"commit":
        raise ValueError("the frozen code reference must identify a Git commit")
    paths = sorted(_relative(path) for path in paths)
    tree = _git(repository, "ls-tree", "-r", "-z", "--full-tree", code_ref, "--", *paths)
    objects = {}
    for record in tree.split(b"\0"):
        if not record:
            continue
        header, encoded_path = record.split(b"\t", 1)
        mode, kind, object_id = header.split()
        path = encoded_path.decode("utf-8")
        if path not in paths or path in objects or kind != b"blob" or mode not in {b"100644", b"100755"}:
            raise ValueError("frozen Git dependencies must be distinct regular files, never symlinks or submodules")
        objects[path] = object_id
    if set(objects) != set(paths):
        raise ValueError("frozen Git snapshot omits required files: " + ", ".join(sorted(set(paths) - objects.keys())))
    content = _git(repository, "cat-file", "--batch", input=b"\n".join(objects[path] for path in paths) + b"\n")
    result, offset = {}, 0
    for path in paths:
        end = content.index(b"\n", offset)
        object_id, kind, length = content[offset:end].split()
        size = int(length)
        start, offset = end + 1, end + 1 + size
        if object_id != objects[path] or kind != b"blob" or content[offset:offset + 1] != b"\n":
            raise ValueError("invalid frozen Git blob response")
        result[path] = content[start:offset]
        offset += 1
    if offset != len(content):
        raise ValueError("unexpected frozen Git blob response")
    return result


def _dependencies(bundle, registry, operation="replay"):
    frozen = bundle.fetch(CHECKPOINT_PATH, original_sha256=CHECKPOINT_SHA256).payload
    manifest = _native_checkpoint_manifest(frozen)
    _same(frozen["model"]["config"]["environment"], registry["environment"], "frozen runtime environment")
    entries = [*manifest["inputs"], *manifest["code"]]
    for variant in _VARIANTS:
        task = bundle.fetch(f"outputs/heldout_transfer_{variant}/task_manifest.json").payload
        score = bundle.fetch(f"outputs/heldout_score_{variant}/score.json").payload
        entries.extend(task["inputs"])
        entries.extend(score["label_manifest"]["inputs"])
    identities, data_inputs = {}, {}
    for entry in entries:
        identity = _identity(entry)
        path = identity["path"]
        if path in identities:
            _same(identity, identities[path], "shared frozen dependency identity")
        identities[path] = identity
        if entry["role"] not in {"training_code", "evaluation_code"}:
            data_inputs[path] = identity
    for entry in registry["replay_support"]:
        identity = _identity(entry)
        if identity["path"] in identities:
            raise ValueError("replay support must not replace a frozen scientific dependency")
        identities[identity["path"]] = identity
    if operation == "audit":
        for entry in registry["audit_support"]:
            identity = _identity(entry)
            path = identity["path"]
            if path in identities:
                raise ValueError("audit support must not replace a frozen scientific dependency")
            identities[path] = identity
            if entry["role"] in {"audit_provenance", "recorded_diagnostics"}:
                data_inputs[path] = identity
    return frozen, identities, data_inputs


def _supplied_data(root, inputs):
    checks, physical_files = [], set()
    for path, identity in sorted(inputs.items()):
        candidate = root / path
        if candidate.resolve() != candidate:
            raise ValueError(f"supplied replay data must use a nonsymlink canonical path: {path}")
        present = candidate.exists()
        if present:
            file = _relative_file(root, path)
            _verified_content(file.read_bytes(), identity, "supplied source data")
            stat = file.stat()
            physical = stat.st_dev, stat.st_ino
            if physical in physical_files:
                raise ValueError("supplied replay inputs must identify distinct physical files")
            physical_files.add(physical)
        checks.append({**identity, "status": "verified" if present else "absent_provisioned_from_verified_git"})
    return checks


def _verify_snapshot(snapshot):
    for path, identity in snapshot.identities.items():
        _verified_content(_relative_file(snapshot.root, path).read_bytes(), identity, "materialized snapshot")


@contextmanager
def materialize_frozen_runtime(root, *, source_repository=None, frozen_code_ref=None,
                               manifest_path=None, expected_sha256=None, operation="replay"):
    if operation not in {"replay", "audit"}:
        raise ValueError("frozen operation must be explicitly 'replay' or 'audit'")
    root = Path(root).resolve(strict=True)
    repository = root if source_repository is None else Path(source_repository).resolve(strict=True)
    registry = load_frozen_runtime(repository)
    if frozen_code_ref is not None and frozen_code_ref != registry["code_ref"]:
        raise ValueError("historical replay requires the registered full frozen Git commit")
    root, manifest_path, bundle = _load_bundle(root, manifest_path, expected_sha256)
    _same({"path": manifest_path, "sha256": bundle.manifest_sha256}, registry["manifest"],
          "registered frozen manifest identity")
    frozen, identities, data_inputs = _dependencies(bundle, registry, operation)
    supplied = _supplied_data(root, data_inputs)
    public_files = {manifest_path: {
        "path": manifest_path, "sha256": bundle.manifest_sha256,
        "size_bytes": _relative_file(root, manifest_path).stat().st_size,
    }}
    for record in bundle.manifest["records"]:
        path = (PurePosixPath(manifest_path).parent / record["public_path"]).as_posix()
        public_files[path] = _identity({
            "path": path, "sha256": record["public_sha256"], "size_bytes": record["public_size_bytes"],
        })
    if identities.keys() & public_files.keys():
        raise ValueError("public evidence storage must not alias executable dependencies")
    identities.update(public_files)
    blobs = _git_blobs(repository, registry["code_ref"], identities)
    for path, identity in identities.items():
        _verified_content(blobs[path], identity, "Git blob")
    for path, identity in public_files.items():
        _verified_content(_relative_file(root, path).read_bytes(), identity, "supplied public evidence")
    with tempfile.TemporaryDirectory(prefix="ystwin-frozen-replay-", dir=Path(tempfile.gettempdir()).resolve(strict=True)) as directory:
        snapshot = FrozenRuntimeSnapshot(Path(directory).resolve(strict=True), registry, identities, supplied, operation)
        for path, data in blobs.items():
            destination = snapshot.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("xb") as handle:
                handle.write(data)
        _verify_native_checkpoint(snapshot.root, frozen)
        yield snapshot
        _verify_snapshot(snapshot)
        _load_bundle(root, manifest_path, bundle.manifest_sha256)
        _same(_supplied_data(root, data_inputs), supplied, "post-replay supplied source identities")
        _same(load_frozen_runtime(repository), registry, "post-replay runtime registry")


_CHILD_REPLAY = r'''
import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import runpy
import sys

options = json.load(sys.stdin)
root = Path(options["root"]).resolve(strict=True)
source = root / "src"
if not sys.flags.isolated or not sys.dont_write_bytecode or os.environ.get("PYTHONPATH") != str(source):
    raise ValueError("historical replay requires isolated Python and its exact snapshot PYTHONPATH")
sys.path.insert(0, str(source))
identities = options["identities"]

def verify_sources():
    for path, identity in identities.items():
        file = root / path
        if file.resolve() != file or not file.is_file():
            raise ValueError("historical source alias: " + path)
        data = file.read_bytes()
        if hashlib.sha256(data).hexdigest() != identity["sha256"] or len(data) != identity["size_bytes"]:
            raise ValueError("historical source SHA-256 mismatch: " + path)

def verify_modules():
    loaded = []
    for name, module in sorted(sys.modules.items()):
        if name != "ystwin" and not name.startswith("ystwin."):
            continue
        filename = getattr(module, "__file__", None)
        if filename is None:
            raise ValueError("historical module lacks source identity: " + name)
        file = Path(filename)
        if file.resolve() != file or not file.is_relative_to(source):
            raise ValueError("historical module escaped the snapshot PYTHONPATH: " + name)
        relative = file.relative_to(root).as_posix()
        identity = identities.get(relative)
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        if identity is None or digest != identity["sha256"]:
            raise ValueError("historical module differs from its verified Git source: " + name)
        loaded.append({"module": name, "path": relative, "sha256": digest})
    return loaded

verify_sources()
verify_modules()
expected = options["environment"]
actual = {
    "python": platform.python_version(),
    "packages": {name: importlib.metadata.version(name) for name in expected["packages"]},
}
if actual != expected:
    raise ValueError("frozen runtime versions differ; expected " + json.dumps(expected, sort_keys=True)
                     + "; actual " + json.dumps(actual, sort_keys=True))
if options["operation"] == "replay":
    with contextlib.redirect_stdout(io.StringIO()) as output:
        namespace = runpy.run_path(str(root / options["entrypoint"]))
        exit_code = namespace["main"](options["arguments"])
    if exit_code != 0:
        raise SystemExit(exit_code)
    report = json.loads(output.getvalue())
elif options["operation"] == "audit":
    from ystwin.analysis import biology_learning_audit as audit

    request = options["audit_request"]
    sources = {"artifact_source": "portable", "manifest_path": options["manifest"]["path"],
               "expected_sha256": options["manifest"]["sha256"]}
    ledger = audit.collect_learning_evidence(root, **sources)
    supplied = request["ledger"]
    if supplied is not None:
        audit.require_claim(supplied, "portable_scientific_content", root=root, **sources)
    authorized, rejected = [], []
    for claim in request["required_claims"]:
        try:
            authorized.append(audit.require_claim(ledger, claim, root=root, **sources))
        except ValueError as exc:
            rejected.append({"claim": claim, "reason": str(exc)})
    diagnostics = {}
    origin = {"mode": "not_run_historical_source_audit", "numerical_diagnostics_rerun": False,
              "biological_validation_credit": False}
    if request["reuse_native_diagnostics"]:
        payload, origin = audit._reused_native_diagnostics(root)
        diagnostics = audit._json(payload)
    report = {
        "schema_version": 1, "kind": options["audit_kind"], "verified": True,
        "ledger": ledger, "integrity": ledger["integrity_scope"],
        "requested_claims": request["required_claims"], "authorized_claims": authorized,
        "rejected_required_claims": rejected, "required_claims_satisfied": not rejected,
        "diagnostics_origin": origin, "findings": audit._findings(ledger, diagnostics),
        "historical_diagnostic_code_reexecuted": False, "all_frozen_files_unchanged": True,
        "new_fitting": False, "new_selection": False, "independent_test": False,
        "biological_validation": False, "original_serialization_timestamp_verified": False,
    }
else:
    raise ValueError("unregistered historical operation")
verify_sources()
report["historical_execution"] = {
    "isolated_interpreter": True, "bytecode_disabled": True, "pythonpath": ["src"],
    "runtime_versions": actual, "loaded_modules": verify_modules(),
}
print(json.dumps(report, sort_keys=True, allow_nan=False))
'''


def _run_child(snapshot, python_executable, private_original_root, timeout, *, audit_request=None):
    arguments = ["--root", str(snapshot.root), "--manifest", snapshot.registry["manifest"]["path"],
                 "--manifest-sha256", snapshot.registry["manifest"]["sha256"]]
    if private_original_root is not None:
        arguments.extend(["--private-original-root", str(Path(private_original_root).resolve())])
    options = {
        "root": str(snapshot.root), "identities": snapshot.identities,
        "environment": snapshot.registry["environment"], "arguments": arguments,
        "entrypoint": snapshot.registry["replay_entrypoint"], "operation": snapshot.operation,
        "audit_request": audit_request, "manifest": snapshot.registry["manifest"],
        "audit_kind": HISTORICAL_AUDIT_KIND,
    }
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("PYTHON", "YSTWIN_", "GIT_"))}
    environment.update(PYTHONPATH=str(snapshot.root / "src"), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run(
        [os.fspath(python_executable), "-I", "-B", "-c", _CHILD_REPLAY],
        input=json.dumps(options, allow_nan=False), cwd=snapshot.root, env=environment,
        capture_output=True, text=True, timeout=timeout,
    )


def _child_report(process, snapshot):
    if process.returncode != 0:
        raise HistoricalReplayError("historical replay child failed", returncode=process.returncode,
                                    stdout=process.stdout, stderr=process.stderr)
    try:
        report = _parse_json(process.stdout)
        kind = HISTORICAL_AUDIT_KIND if snapshot.operation == "audit" else REPLAY_KIND
        if type(report["schema_version"]) is not int or report["schema_version"] != 1 or report["kind"] != kind or report["verified"] is not True:
            raise ValueError("child did not verify the recorded evidence")
        for field in ("new_fitting", "new_selection", "independent_test", "biological_validation",
                      "original_serialization_timestamp_verified"):
            if report[field] is not False:
                raise ValueError("child exceeded the historical reproduction scope")
        integrity = report["integrity"]
        _same(integrity["manifest_reference"], snapshot.registry["manifest"], "child manifest identity")
        _same(integrity["checkpoint_sha256"], CHECKPOINT_SHA256, "child checkpoint identity")
        _same(integrity["model_sha256"], MODEL_SHA256, "child model identity")
        for field in ("recorded_versions", "actual_versions"):
            _same(integrity["runtime"][field], snapshot.registry["environment"], "child runtime versions")
        execution = report["historical_execution"]
        if execution["isolated_interpreter"] is not True or execution["bytecode_disabled"] is not True:
            raise ValueError("child did not establish interpreter isolation")
        _same(execution["pythonpath"], ["src"], "child snapshot PYTHONPATH")
        _same(execution["runtime_versions"], snapshot.registry["environment"], "child isolated runtime versions")
        paths = set(FROZEN_CODE_PATHS) | {entry["path"] for entry in snapshot.registry["replay_support"]}
        if snapshot.operation == "audit":
            paths.update(entry["path"] for entry in snapshot.registry["audit_support"] if entry["role"] == "audit_code")
        expected = {path: snapshot.identities[path]["sha256"] for path in paths
                    if path.startswith("src/") and path.endswith(".py")}
        loaded = {row["path"]: row["sha256"] for row in execution["loaded_modules"]}
        _same(loaded, expected, "child loaded source identities")
        if len(loaded) != len(execution["loaded_modules"]):
            raise ValueError("child loaded aliased source modules")
    except (KeyError, TypeError, ValueError) as exc:
        raise HistoricalReplayError("invalid historical child report: " + str(exc), returncode=process.returncode,
                                    stdout=process.stdout, stderr=process.stderr) from exc
    return report


def _checked_child(snapshot, python_executable, private_original_root, timeout, *, audit_request=None):
    executable = sys.executable if python_executable is None else python_executable
    try:
        if audit_request is None:
            process = _run_child(snapshot, executable, private_original_root, timeout)
        else:
            process = _run_child(snapshot, executable, private_original_root, timeout, audit_request=audit_request)
    except subprocess.TimeoutExpired as exc:
        def text(value):
            return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""

        raise HistoricalReplayError("historical replay child timed out", stdout=text(exc.stdout),
                                    stderr=text(exc.stderr)) from exc
    except OSError as exc:
        raise HistoricalReplayError("historical replay child could not start: " + str(exc)) from exc
    return process, _child_report(process, snapshot)


def _execution_receipt(report, snapshot, process):
    report["historical_execution"].update({
        "operation": snapshot.operation,
        "code_ref": snapshot.registry["code_ref"], "child_exit_code": process.returncode,
        "registry_reference": {"path": RUNTIME_REGISTRY_PATH, "sha256": RUNTIME_REGISTRY_SHA256},
        "snapshot_files": [snapshot.identities[path] for path in sorted(snapshot.identities)],
        "supplied_data": snapshot.supplied_data, "git_bytes_verified_before_execution": True,
        "current_worktree_code_used": False, "snapshot_and_supplied_inputs_unchanged": True,
        "scope": snapshot.registry["audit_scope" if snapshot.operation == "audit" else "scope"],
    })
    return report


def replay_historical_evidence(root, *, source_repository=None, frozen_code_ref=None, python_executable=None,
                               manifest_path=None, expected_sha256=None, private_original_root=None, timeout=300):
    with materialize_frozen_runtime(
        root, source_repository=source_repository, frozen_code_ref=frozen_code_ref,
        manifest_path=manifest_path, expected_sha256=expected_sha256,
    ) as snapshot:
        process, report = _checked_child(snapshot, python_executable, private_original_root, timeout)
    return _execution_receipt(report, snapshot, process)


def _audit_report(process, snapshot, request):
    report = _child_report(process, snapshot)
    try:
        _same(report["requested_claims"], request["required_claims"], "historical requested claims")
        _same(report["ledger"]["integrity_scope"], report["integrity"], "historical audit integrity")
        _same(report["ledger"]["artifact_source"], "portable", "historical audit source")
        if (report["historical_diagnostic_code_reexecuted"] is not False
                or report["all_frozen_files_unchanged"] is not True):
            raise ValueError("historical audit exceeded its read-only source scope")
        authorized = report["authorized_claims"]
        rejected = report["rejected_required_claims"]
        assessed = [row["claim"] for row in (*authorized, *rejected)]
        _same(sorted(assessed), sorted(request["required_claims"]), "complete historical claim assessment")
        _same(report["required_claims_satisfied"], not rejected, "historical claim gate outcome")
        for row in authorized:
            if row["claim"] != "portable_scientific_content" or row["authorized"] is not True:
                raise ValueError("historical public content cannot authorize an original-byte or biological grade")
            _same(row["integrity_scope"], report["integrity"], "historical claim integrity")
        origin = report["diagnostics_origin"]
        if origin["numerical_diagnostics_rerun"] is not False or origin["biological_validation_credit"] is not False:
            raise ValueError("historical diagnostics cannot be refitted or upgraded")
        if request["reuse_native_diagnostics"]:
            _same(origin["mode"], "reused_unchanged_retrospective_evidence", "historical diagnostics mode")
            for field, path in (("source", "outputs/biology_learning_audit_01/native_diagnostics.json"),
                                ("original_plan", "outputs/biology_learning_audit_01/diagnostic_plan.json")):
                reference = origin[field]
                _same(reference["path"], path, "reused historical diagnostic path")
                _same(reference["sha256"], snapshot.identities[path]["sha256"],
                      "reused historical diagnostic identity")
        else:
            _same(origin["mode"], "not_run_historical_source_audit", "historical diagnostics mode")
    except (KeyError, TypeError, ValueError) as exc:
        raise HistoricalReplayError("invalid historical audit report: " + str(exc), returncode=process.returncode,
                                    stdout=process.stdout, stderr=process.stderr) from exc
    return report


def audit_historical_evidence(root, *, source_repository=None, frozen_code_ref=None, python_executable=None,
                              manifest_path=None, expected_sha256=None, required_claims=(), ledger=None,
                              reuse_native_diagnostics=False, output_dir=None, timeout=300):
    if (not isinstance(required_claims, (list, tuple))
            or any(not isinstance(claim, str) or not claim for claim in required_claims)
            or len(set(required_claims)) != len(required_claims)):
        raise ValueError("historical audit requires distinct explicit claim names")
    if type(reuse_native_diagnostics) is not bool:
        raise ValueError("reuse_native_diagnostics must be an explicit boolean")
    if ledger is not None and not isinstance(ledger, dict):
        raise ValueError("historical audit ledger must be an explicit JSON object")
    root = Path(root).resolve(strict=True)
    output = None
    if output_dir is not None:
        output = Path(output_dir)
        output = output if output.is_absolute() else root / output
        if output.resolve() != output or not output.parent.is_dir():
            raise ValueError("historical audit output requires an existing nonsymlink parent and a canonical path")
        if output.exists():
            raise FileExistsError("historical audit output is write-once; existing records cannot be overwritten")
    request = {"required_claims": list(required_claims), "ledger": ledger,
               "reuse_native_diagnostics": reuse_native_diagnostics}
    with materialize_frozen_runtime(
        root, source_repository=source_repository, frozen_code_ref=frozen_code_ref,
        manifest_path=manifest_path, expected_sha256=expected_sha256, operation="audit",
    ) as snapshot:
        process, _ = _checked_child(snapshot, python_executable, None, timeout, audit_request=request)
        report = _audit_report(process, snapshot, request)
        diagnostics = None
        if reuse_native_diagnostics:
            reference = report["diagnostics_origin"]["source"]
            diagnostics = _verified_content(_relative_file(snapshot.root, reference["path"]).read_bytes(),
                                            snapshot.identities[reference["path"]], "reused diagnostics")
    report = _execution_receipt(report, snapshot, process)
    if output is not None:
        if output.resolve() != output or not output.parent.is_dir():
            raise ValueError("historical audit output parent changed during execution")
        output.mkdir(exist_ok=False)
        if diagnostics is not None:
            with (output / "native_diagnostics.json").open("xb") as handle:
                handle.write(diagnostics)
        for name, value in (("ledger.json", report["ledger"]), ("audit_report.json", report)):
            with (output / name).open("x", encoding="utf-8") as handle:
                handle.write(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    return report
