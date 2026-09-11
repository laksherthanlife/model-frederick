"""Independent tests of the portable evidence boundary, not biological validation.

Run with --noconftest -p no:cacheprovider and PYTHONPATH=src:. so this review
never loads machine-local data configuration. All destructive probes must use
pytest temporary copies; retained originals and existing fixtures are read-only.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest

from _frozen_runtime_helpers import run_frozen_python


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_RELEASE = "data/frozen_evidence/native_v1"
NATIVE_CHECKPOINT = "outputs/native_training_run_02/checkpoint.json"
NATIVE_MANIFEST = "outputs/native_training_run_02/training_manifest.json"
POPULATION_FREEZE = (
    "outputs/native_population_development/native_training_01/development_prediction_freeze.json"
)
ORIGINAL_PINS = {
    NATIVE_CHECKPOINT: "384175e19b2451f542ab08e347ee9d9021e22772fbc205e0c404bbc28416ec94",
    POPULATION_FREEZE: "e03f6b9407b29e7e76da7a0be566ce9036c5cf27ac65bfaa9a4d3d7e009b5753",
}


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _pointer(document, pointer):
    for part in pointer.split("/")[1:]:
        key = part.replace("~1", "/").replace("~0", "~")
        document = document[int(key)] if isinstance(document, list) else document[key]
    return document


def _set_pointer(document, pointer, value):
    parent, _, key = pointer.rpartition("/")
    node = _pointer(document, parent)
    key = key.replace("~1", "/").replace("~0", "~")
    node[int(key) if isinstance(node, list) else key] = value


def _records(case):
    return {record["origin_path"]: record for record in case.manifest["records"]}


def _load(case, *, pin=None):
    from ystwin.analysis import portable_evidence

    return portable_evidence.load_portable_evidence(
        case.root, case.manifest_path, case.pin if pin is None else pin
    )


def _save_manifest(case):
    raw = _json_bytes(case.manifest)
    (case.root / case.manifest_path).write_bytes(raw)
    return _sha(raw)


def _envelope(case, origin):
    record = _records(case)[origin]
    return json.loads((case.directory / record["public_path"]).read_bytes())


def _rewrite_payload(case, origin, pointer, value):
    """Forge only temporary bytes plus all caller-computable content digests."""
    record = _records(case)[origin]
    envelope = _envelope(case, origin)
    _set_pointer(envelope["payload"], pointer, value)
    raw = _json_bytes(envelope)
    (case.directory / record["public_path"]).write_bytes(raw)
    record.update(public_sha256=_sha(raw), public_size_bytes=len(raw))
    science = _sha(_canonical(envelope["payload"]))
    record["scientific_payload"].update(original_sha256=science, public_sha256=science)
    for link in case.manifest["links"]:
        if link["target_origin_path"] == origin:
            link["target_public_sha256"] = record["public_sha256"]
    case.manifest["links"].sort(key=_canonical)
    return envelope


@pytest.fixture(scope="module")
def public_release():
    """Use only the distributable release; never rebuild from hidden originals."""
    directory = ROOT / PUBLIC_RELEASE
    manifest_raw = (directory / "manifest.json").read_bytes()
    manifest = json.loads(manifest_raw)
    files = {}
    for record in manifest["records"]:
        relative = PurePosixPath(record["public_path"])
        assert not relative.is_absolute() and ".." not in relative.parts
        assert relative.as_posix() == record["public_path"]
        path = directory / relative
        assert path.resolve() == path and path.is_file()
        raw = path.read_bytes()
        assert _sha(raw) == record["public_sha256"]
        files[record["public_path"]] = raw
    return manifest_raw, files


@pytest.fixture
def portable_copy(tmp_path, public_release):
    manifest_raw, files = public_release
    directory = tmp_path / "release"
    for relative, raw in {"manifest.json": manifest_raw, **files}.items():
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    return SimpleNamespace(
        root=tmp_path, directory=directory, manifest_path="release/manifest.json",
        manifest=json.loads(manifest_raw), pin=_sha(manifest_raw),
    )


@pytest.mark.parametrize("module", ["portable_evidence", "portable_replay"])
def test_portable_python_interfaces_can_be_imported(module):
    """Missing implementation is a failure, never a private-fixture skip."""
    assert importlib.import_module(f"ystwin.analysis.{module}") is not None


def test_public_replay_cli_exists():
    assert (ROOT / "scripts/replay_portable_evidence.py").is_file()


@pytest.mark.parametrize("mutation", [
    "coefficient", "unit", "condition_role", "censored_observation", "source_coordinate",
])
def test_exchange_replay_preserves_source_bound_values_units_roles_and_censoring(mutation):
    """Exercise the actual replay helper using only published native measurements."""
    from ystwin.analysis import native_physiology, portable_replay

    source = native_physiology.load_native_chemostat_data(
        ROOT / "data/physiology/chemostatData_VanHoek1998.tsv"
    )
    payload = json.loads(
        (ROOT / "outputs/native_training_run_02/native_exchange_fit.json").read_bytes()
    )
    assert portable_replay.restore_native_exchange_model(payload, source).to_dict() == payload
    changed = deepcopy(payload)
    if mutation == "coefficient":
        changed["knots"]["Acetate"][-1]["value"] += 0.125
    elif mutation == "unit":
        changed["input"]["unit"] = "1/min"
    elif mutation == "condition_role":
        changed["input"]["role"] = "scored_prediction"
    elif mutation == "censored_observation":
        assert changed["knots"]["Acetate"][0]["value"] is None
        changed["knots"]["Acetate"][0]["value"] = 0.0
    else:
        changed["knots"]["Acetate"][-1]["source_line"] += 1
    with pytest.raises(ValueError, match="differs"):
        portable_replay.restore_native_exchange_model(changed, source)


def test_interval_replay_keeps_failure_rows_denominators_and_explicit_units():
    from ystwin.analysis import portable_replay

    sealed = {
        "prediction_kind": "feasible_envelope",
        "prediction_scope": "synthetic_example_not_biological_validation",
        "task_manifest": {"quantity": "response", "unit": "a.u."},
        "predictions": [
            {"experiment_id": "a", "status": "ok", "lower": 1.0,
             "upper": 3.0, "reason": None},
            {"experiment_id": "b", "status": "ok", "lower": 2.0,
             "upper": 4.0, "reason": None},
            {"experiment_id": "c", "status": "infeasible", "lower": None,
             "upper": None, "reason": "declared constraint"},
            {"experiment_id": "d", "status": "unsupported", "lower": None,
             "upper": None, "reason": "outside declared support"},
        ],
    }
    observed = {
        identifier: {"value": value, "source_path": "data/example_observations.json"}
        for identifier, value in zip("abcd", [2.0, 6.0, 3.0, 4.0], strict=True)
    }
    result = portable_replay.recompute_interval_metrics(sealed, observed)
    assert result["unit"] == "a.u."
    assert result["n_total"] == len(result["rows"]) == 4
    assert result["n_scorable"] == 2
    assert result["n_contained"] == result["n_violated"] == 1
    assert result["complete"] is False
    assert result["containment_on_scorable"] == 0.5
    assert result["containment_fraction_total"] == 0.25
    assert result["mean_width_on_scorable"] == 2.0
    assert result["mean_violation_on_scorable"] == 1.0
    assert result["max_violation_on_scorable"] == 2.0
    assert result["status_counts"] == {
        "ok": 2, "failed": 0, "infeasible": 1, "unsupported": 1,
    }
    assert [row["violation"] for row in result["rows"]] == [0.0, 2.0, None, None]


def test_public_bundle_loads_without_originals_and_does_not_verify_their_old_sha(
    portable_copy, monkeypatch,
):
    original_read = Path.read_bytes
    opened = []

    def public_only(path):
        assert path.resolve().is_relative_to(portable_copy.root), (
            "A portable load attempted to consult files outside its relocated public root"
        )
        opened.append(path)
        return original_read(path)

    monkeypatch.setattr(Path, "read_bytes", public_only)
    bundle = _load(portable_copy)
    assert opened
    for origin, expected in ORIGINAL_PINS.items():
        assert not (portable_copy.root / origin).exists()
        artifact = bundle.fetch(origin, original_sha256=expected)
        scope = artifact.integrity_scope
        assert scope["original_sha256"] == expected
        assert scope["public_sha256"] != expected
        assert scope["public_content_verified"] is True
        assert scope["original_bytes_verified_on_load"] is False
        assert scope["new_scientific_claim_authorized"] is False
        assert scope["new_export_was_sealed_before_data"] is False
        assert bundle.fetch_public(scope["public_sha256"]).payload == artifact.payload


def test_plain_json_payloads_and_exact_byte_companions_have_distinct_honest_lineage(portable_copy):
    bundle = _load(portable_copy)
    for origin, record in _records(portable_copy).items():
        content = (portable_copy.directory / record["public_path"]).read_bytes()
        fetched = bundle.fetch(origin, original_sha256=record["original_sha256"])
        assert _sha(content) == record["public_sha256"]
        if record["format"] == "json_envelope":
            document = json.loads(content)
            assert set(document) == {"schema_version", "kind", "origin", "payload"}
            assert isinstance(document["payload"], dict)
            assert fetched.payload == document["payload"]
            assert document["origin"] == {"path": origin, "sha256": record["original_sha256"]}
            assert record["public_sha256"] != record["original_sha256"]
            assert fetched.integrity_scope["original_bytes_verified_on_load"] is False
            assert _sha(_canonical(document["payload"])) == record["scientific_payload"]["public_sha256"]
        else:
            assert record["format"] == "original_bytes"
            assert fetched.payload == content
            assert _sha(content) == record["original_sha256"]
            assert fetched.integrity_scope["original_bytes_verified_on_load"] is True


def test_returned_manifest_payloads_and_scopes_are_not_mutable_trust_state(portable_copy):
    bundle = _load(portable_copy)
    expected = bundle.fetch(NATIVE_CHECKPOINT)
    first = bundle.fetch(NATIVE_CHECKPOINT)
    first.payload["model"]["parameters"].clear()
    first.integrity_scope["public_content_verified"] = False
    first.integrity_scope["scientific_payload"]["public_sha256"] = "0" * 64
    manifest = bundle.manifest
    manifest["records"].clear()
    assert bundle.fetch(NATIVE_CHECKPOINT) == expected
    assert bundle.manifest == portable_copy.manifest


@pytest.mark.parametrize("forged_internal_digests", [False, True])
@pytest.mark.parametrize("mutation", ["coefficient", "prediction", "unit", "source_hash", "refpath"])
def test_mutated_public_science_cannot_be_laundered_by_internal_rehashing(
    portable_copy, mutation, forged_internal_digests,
):
    _load(portable_copy)
    if mutation == "coefficient":
        model = _envelope(portable_copy, NATIVE_CHECKPOINT)["payload"]["model"]
        key = sorted(model["parameters"])[0]
        pointer = "/model/parameters/" + key.replace("~", "~0").replace("/", "~1")
        _rewrite_payload(portable_copy, NATIVE_CHECKPOINT, pointer, model["parameters"][key] + 0.125)
    elif mutation == "prediction":
        value = _envelope(portable_copy, POPULATION_FREEZE)["payload"]["predictions"][0]["mean"]
        _rewrite_payload(portable_copy, POPULATION_FREEZE, "/predictions/0/mean", value + 0.25)
    elif mutation == "unit":
        _rewrite_payload(portable_copy, NATIVE_CHECKPOINT, "/model/config/unit", "changed-unit")
    elif mutation == "source_hash":
        _rewrite_payload(portable_copy, NATIVE_CHECKPOINT, "/training_manifest/inputs/0/sha256", "0" * 64)
    else:
        _rewrite_payload(portable_copy, NATIVE_CHECKPOINT, "/training_manifest/code/0/path", "../unapproved.py")
    if forged_internal_digests:
        assert _save_manifest(portable_copy) != portable_copy.pin
    with pytest.raises(ValueError, match="SHA-256|digest|identity"):
        _load(portable_copy)


@pytest.mark.parametrize("mutation", ["unknown_location", "wrong_operation", "wrong_replacement"])
def test_even_a_new_manifest_pin_cannot_authorize_unregistered_scientific_redaction(
    portable_copy, mutation,
):
    """A newly trusted fixture pin does not waive the transform schema contract."""
    _load(portable_copy)
    transform = _records(portable_copy)[NATIVE_CHECKPOINT]["transforms"][0]
    if mutation == "unknown_location":
        transform["json_pointer"] = "/model/config/unit"
    elif mutation == "wrong_operation":
        transform["operation"] = "replace_scientific_payload"
    else:
        transform["replacement"] = str(portable_copy.root / "private-context")
    new_pin = _save_manifest(portable_copy)
    with pytest.raises(ValueError, match="transform|redaction"):
        _load(portable_copy, pin=new_pin)


@pytest.mark.parametrize("mutation", ["public_path", "envelope_origin", "historical_hash_as_public"])
def test_copies_cannot_substitute_a_different_logical_origin_or_claim_old_byte_hashes(
    portable_copy, mutation,
):
    _load(portable_copy)
    record = _records(portable_copy)[NATIVE_CHECKPOINT]
    if mutation == "public_path":
        record["public_path"] = _records(portable_copy)[POPULATION_FREEZE]["public_path"]
    elif mutation == "envelope_origin":
        document = _envelope(portable_copy, NATIVE_CHECKPOINT)
        document["origin"]["path"] = POPULATION_FREEZE
        raw = _json_bytes(document)
        (portable_copy.directory / record["public_path"]).write_bytes(raw)
        record.update(public_sha256=_sha(raw), public_size_bytes=len(raw))
    else:
        record["public_sha256"] = record["original_sha256"]
    new_pin = _save_manifest(portable_copy)
    with pytest.raises(ValueError, match="storage|identity|SHA-256"):
        _load(portable_copy, pin=new_pin)


def test_symlinked_public_copy_is_refused_before_target_bytes_are_read(portable_copy, monkeypatch):
    record = _records(portable_copy)[NATIVE_CHECKPOINT]
    path = portable_copy.directory / record["public_path"]
    outside = portable_copy.root / "unexported-copy.json"
    outside.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(outside)
    original_read = Path.read_bytes

    def guarded(candidate):
        assert candidate.resolve() != outside.resolve(), "A symlink target was opened"
        return original_read(candidate)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    with pytest.raises(ValueError, match="symlink"):
        _load(portable_copy)


def test_lineage_only_nodes_are_not_loadable_payloads(portable_copy, monkeypatch):
    bundle = _load(portable_copy)

    def no_more_reads(path):
        raise AssertionError("A historical boundary reference attempted another file read")

    monkeypatch.setattr(Path, "read_bytes", no_more_reads)
    assert bundle.manifest["external_references"]
    for reference in bundle.manifest["external_references"]:
        with pytest.raises(ValueError):
            bundle.fetch(reference["origin_path"], original_sha256=reference["original_sha256"])


def test_exporter_rejects_rewritten_original_bytes_before_writing_any_export(tmp_path):
    from ystwin.analysis import portable_evidence

    path = tmp_path / NATIVE_CHECKPOINT
    path.parent.mkdir(parents=True)
    original = _json_bytes({
        "schema_version": 1, "kind": "native_checkpoint",
        "training_manifest": {"root": str(tmp_path)},
        "model": {"parameters": {"coefficient": 2.0}, "config": {"unit": "a.u."}},
    })
    path.write_bytes(original)
    destination = tmp_path / "public-export"
    with pytest.raises(ValueError, match="original checkpoint SHA-256 mismatch"):
        portable_evidence.export_portable_evidence(tmp_path, destination)
    assert path.read_bytes() == original
    assert not destination.exists()


def _forge_active_code_reference(case, new_path):
    """Keep every internal digest/link consistent while changing a typed active input."""
    manifest = _envelope(case, NATIVE_MANIFEST)["payload"]
    old = deepcopy(manifest["code"][0])
    assert old["role"] == "training_code"
    assert old["source"]["scope"] == "native_training_code"
    manifest["code"][0].update(path=new_path, resolved_path=new_path)
    manifest["sha256"] = _sha(_canonical({key: value for key, value in manifest.items() if key != "sha256"}))
    _rewrite_payload(case, NATIVE_MANIFEST, "/code", manifest["code"])
    _rewrite_payload(case, NATIVE_MANIFEST, "/sha256", manifest["sha256"])
    _rewrite_payload(case, NATIVE_CHECKPOINT, "/training_manifest", manifest)
    expected_pointers = {
        (NATIVE_CHECKPOINT, "/training_manifest/code/0/path"),
        (NATIVE_MANIFEST, "/code/0/path"),
    }
    affected = 0
    for link in case.manifest["links"]:
        if (link["source_origin_path"], link["json_pointer"]) in expected_pointers:
            assert link["target_origin_path"] == old["path"]
            link.update(target_origin_path=new_path, target_public_sha256=None)
            affected += 1
    assert affected == 2
    external = {
        (row["origin_path"], row["original_sha256"]): row
        for row in case.manifest["external_references"]
    }
    external[new_path, old["sha256"]] = {
        "origin_path": new_path, "original_sha256": old["sha256"],
        "integrity_scope": "lineage_only_not_exported_or_independently_verified",
    }
    needed = {
        (link["target_origin_path"], link["target_original_sha256"])
        for link in case.manifest["links"] if link["target_public_sha256"] is None
    }
    case.manifest["external_references"] = [external[key] for key in sorted(needed)]
    case.manifest["links"].sort(key=_canonical)
    return _save_manifest(case)


@pytest.mark.parametrize("path_kind", ["root_escape", "remote", "private_namespace", "unbound_local"])
def test_pinned_graph_must_not_reclassify_active_training_code_as_lineage_only(
    portable_copy, path_kind,
):
    """Semantic graph validation is still required even with a supplied fixture pin.

    This is not a claim that an attacker can replace the externally retained pin.
    The existing structured training_code role, rather than an unfamiliar filename,
    distinguishes these active dependencies from citation-only source metadata.
    """
    _load(portable_copy)
    new_path = {
        "root_escape": "../unapproved.py",
        "remote": "https://example.invalid/unapproved.py",
        "private_namespace": "data/native_law_v2/granados/sealed/raw_mixed/unapproved.py",
        "unbound_local": "unexported/active_input.py",
    }[path_kind]
    new_pin = _forge_active_code_reference(portable_copy, new_path)
    with pytest.raises(ValueError):
        _load(portable_copy, pin=new_pin)


def _materialize_public_runtime(case):
    """Copy public source and schema-classified dependencies, never frozen originals."""
    from ystwin.analysis import frozen_runtime, portable_replay

    bundle = _load(case)
    checkpoint = bundle.fetch(NATIVE_CHECKPOINT).payload
    registration = frozen_runtime.load_frozen_runtime(ROOT)
    # -I suppresses user-site discovery. Make only the recorded installed package
    # roots visible, without executing site customization or consulting local data config.
    dependency_roots = set()
    for name, version in checkpoint["model"]["config"]["environment"]["packages"].items():
        distribution = importlib.metadata.distribution(name)
        assert distribution.version == version
        dependency_roots.add(str(Path(distribution.locate_file("")).resolve()))
    case.dependency_roots = sorted(dependency_roots)
    dependencies = {}
    allowed_roles = {
        "fit", "prior", "provenance", "selection", "training_code",
        "evaluation_code", "task_structure", "task_conditions", "outcome_labels",
    }

    def add(entries):
        for entry in entries:
            assert entry["role"] in allowed_roles
            if entry["path"] in dependencies:
                assert dependencies[entry["path"]] == entry["sha256"]
            dependencies[entry["path"]] = entry["sha256"]

    add(checkpoint["training_manifest"]["inputs"])
    add(checkpoint["training_manifest"]["code"])
    variants = ("nad", "oxygen_peroxide", "oxygen_water")
    for variant in variants:
        task = bundle.fetch(f"outputs/heldout_transfer_{variant}/task_manifest.json").payload
        score = bundle.fetch(f"outputs/heldout_score_{variant}/score.json").payload
        add(task["inputs"])
        add(score["label_manifest"]["inputs"])
    label_manifest = "data/holdout_transfer/labels_manifest.json"
    dependencies[label_manifest] = _records(case)[label_manifest]["original_sha256"]
    allowed_paths = (
        set(portable_replay.NATIVE_INPUT_ROLES) | set(portable_replay.FROZEN_CODE_PATHS)
        | {label_manifest, "data/holdout_transfer/labels.json"}
        | {f"data/holdout_transfer/hypotheses/{kind}_{variant}.json"
           for variant in variants for kind in ("conditions", "chemistry")}
    )
    assert set(dependencies) <= allowed_paths
    for entry in (*registration["replay_support"], *registration["audit_support"]):
        assert entry["path"] not in dependencies
        dependencies[entry["path"]] = entry["sha256"]
    frozen = frozen_runtime._git_blobs(ROOT, registration["code_ref"], dependencies)
    for relative, digest in dependencies.items():
        raw = frozen[relative]
        assert _sha(raw) == digest
        destination = case.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(raw)
    for relative, digest in dependencies.items():
        assert _sha((case.root / relative).read_bytes()) == digest
    # Existing collector-owned provenance, already part of the public checkout.
    provenance = "data/gem/MANIFEST.md"
    assert _sha((case.root / provenance).read_bytes()) == dependencies[provenance]
    for origin in ORIGINAL_PINS:
        assert not (case.root / origin).exists()
    assert not (case.root / ".ystwin.env").exists()
    return bundle


_ISOLATED_CLI_RUNNER = """
import json
import os
from pathlib import Path
import runpy
import sys

original = Path(sys.argv[1]).resolve()
checkout = Path(sys.argv[2]).resolve()
dependency_roots = json.loads(sys.argv[3])
script = checkout / sys.argv[4]
sys.argv = [str(script), *sys.argv[5:]]
sys.path[:0] = [str(checkout / 'src'), *dependency_roots]

def forbid_original_checkout(event, arguments):
    if event == 'open' and isinstance(arguments[0], (str, bytes)):
        candidate = Path(os.fsdecode(arguments[0])).resolve()
        if candidate.is_relative_to(original):
            raise PermissionError('The original checkout is unavailable in this relocation test')

sys.addaudithook(forbid_original_checkout)
runpy.run_path(str(script), run_name='__main__')
"""


def _run_relocated_cli(case, arguments):
    return subprocess.run(
        [sys.executable, "-I", "-B", "-c", _ISOLATED_CLI_RUNNER,
         str(ROOT), str(case.root), json.dumps(case.dependency_roots),
         "scripts/replay_portable_evidence.py", *arguments],
        cwd=case.root, text=True, capture_output=True, timeout=120,
        env={"HOME": str(case.root / "isolated-user"), "LANG": "C.UTF-8",
             "PYTHONPATH": str(case.root / "src"), "PYTHONDONTWRITEBYTECODE": "1",
             "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "CI": "1",
             "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )


@pytest.mark.parametrize("mode", ["explicit_pin", "default_pin", "missing_private_originals"])
def test_relocated_cli_replays_real_metrics_without_private_originals_or_biology_upgrade(
    portable_copy, mode,
):
    bundle = _materialize_public_runtime(portable_copy)
    arguments = ["--root", str(portable_copy.root), "--output", "replayed.json"]
    if mode == "default_pin":
        destination = portable_copy.root / PUBLIC_RELEASE
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(portable_copy.directory, destination)
    else:
        arguments += ["--manifest", portable_copy.manifest_path,
                      "--manifest-sha256", portable_copy.pin]
    if mode == "missing_private_originals":
        empty = portable_copy.root / "unavailable-originals"
        empty.mkdir()
        arguments += ["--private-original-root", str(empty)]
    completed = _run_relocated_cli(portable_copy, arguments)
    assert completed.returncode == 0, completed.stderr
    report = json.loads(completed.stdout)
    assert json.loads((portable_copy.root / "replayed.json").read_bytes()) == report
    assert report["kind"] == "portable_replay_of_previously_recorded_content"
    assert report["verified"] is True
    for field in ("new_fitting", "new_selection", "independent_test", "biological_validation",
                  "original_serialization_timestamp_verified"):
        assert report[field] is False
    original = report["integrity"]["original_checkpoint_bytes"]
    assert original["expected_sha256"] == ORIGINAL_PINS[NATIVE_CHECKPOINT]
    assert original["status"] == "not_checked"
    assert original["actual_sha256"] is None and original["matches"] is None
    assert report["integrity"]["independent_biological_validation"] is False
    if mode == "missing_private_originals":
        assert len(report["private_original_checks"]) == len(bundle.manifest["records"])
        assert all(row["status"] == "not_checked" and row["actual_sha256"] is None
                   and row["matches"] is None for row in report["private_original_checks"])
    else:
        assert report["private_original_checks"] is None
    recorded_fit = bundle.fetch("outputs/native_training_run_02/native_hog_fit.json").payload
    assert report["native_hog"]["training_nrmse"] == pytest.approx(
        recorded_fit["diagnostics"]["training_nrmse"], rel=1e-10, abs=1e-14
    )
    assert len(report["native_hog"]["rows"]) == 82
    for name, metrics in report["native_hog"]["per_observable"].items():
        recorded = recorded_fit["diagnostics"]["per_observable_errors"][name]
        assert metrics["unit"] == recorded["unit"] == "mol/L"
        assert metrics["rows"] == recorded["rows"]
        for key in ("rmse", "nrmse"):
            assert metrics[key] == pytest.approx(recorded[key], rel=1e-10, abs=1e-14)
    assert report["native_exchange"]["computed_report"] == bundle.fetch(
        "outputs/native_training_run_02/native_exchange_validation.json"
    ).payload
    assert [row["variant"] for row in report["transfer_scores"]] == [
        "nad", "oxygen_peroxide", "oxygen_water",
    ]
    for score in report["transfer_scores"]:
        recorded = bundle.fetch(f"outputs/heldout_score_{score['variant']}/score.json").payload
        assert score["score_verified"] is True
        assert score["prediction_kernel_reexecuted"] is False
        for key, value in score["computed_metrics"].items():
            assert value == recorded[key], key
    population = report["population_metrics"]
    assert population["verified"] is True
    assert population["prediction_count"] == population["evaluated_frame_rows"] == 288
    assert population["unit"] == "dimensionless_ratio"
    assert population["raw_source_verification"] == "not_checked"
    for field in ("prediction_kernel_reexecuted", "selection_reexecuted",
                  "strong_claim_authorized", "final_test_defined"):
        assert population[field] is False
    population_summary = bundle.fetch(
        "outputs/native_population_development/native_training_01_phase_c_01_scoring/scoring_summary.json"
    ).payload
    assert set(population["models"]) == set(population_summary["models"])
    for model, metrics in population["models"].items():
        for field in ("mean_mse", "rmse", "mae", "scaled_rmse", "scaled_mae"):
            assert metrics[field] == pytest.approx(
                population_summary["models"][model][field], rel=1e-10, abs=1e-14
            )


def test_portable_collector_requires_a_distinct_grade_and_cannot_forge_original_integrity(
    portable_copy,
):
    _materialize_public_runtime(portable_copy)
    completed = run_frozen_python(portable_copy.root, """
        from copy import deepcopy
        import hashlib
        import pytest
        from ystwin.analysis import biology_learning_audit as audit

        options = json.load(sys.stdin)
        ledger = audit.collect_learning_evidence(root, **options)
        assert ledger["evidence_predicates"]["checkpoint_integrity"]["value"] is None
        assert ledger["integrity_scope"]["original_checkpoint_bytes"]["matches"] is None
        scoped = audit.require_claim(ledger, "portable_scientific_content", root=root, **options)
        assert scoped["authorized"] is True
        assert scoped["claim"] == "portable_scientific_content"
        for claim in ("integrity_established", "independent_mechanistic_evidence",
                      "independent_native_validation"):
            with pytest.raises(ValueError, match="unsupported claim"):
                audit.require_claim(ledger, claim, root=root, **options)
        forged = deepcopy(ledger)
        forged["evidence_predicates"]["checkpoint_integrity"]["value"] = True
        canonical = json.dumps(forged, sort_keys=True, separators=(",", ":"), allow_nan=False)
        forged["sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
        with pytest.raises(ValueError, match="authoritative collector"):
            audit.require_claim(forged, "integrity_established", root=root, **options)
    """, payload={
        "artifact_source": "portable", "manifest_path": portable_copy.manifest_path,
        "expected_sha256": portable_copy.pin,
    })
    assert completed.returncode == 0, completed.stderr


def test_public_population_score_arithmetic_and_all_frozen_predictions_are_unchanged(portable_copy):
    """Replay published aggregate rows, without claiming a new raw-source reanalysis."""
    bundle = _load(portable_copy)
    seal = bundle.fetch(POPULATION_FREEZE, original_sha256=ORIGINAL_PINS[POPULATION_FREEZE]).payload
    run = "outputs/native_population_development/native_training_01_phase_c_01_scoring"
    report = bundle.fetch(f"{run}/development_report.json").payload
    summary = bundle.fetch(f"{run}/scoring_summary.json").payload
    expected_mse = {
        "prefix_persistence": 0.0588983830872367,
        "training_constant": 0.01404106645401573,
        "prefix_offset": 0.013467129500997396,
        "prefix_affine_time": 0.013838253929103959,
        "first_order_response": 0.010543455061737508,
        "transient_response": 0.006773738644263618,
    }
    predictions = {
        (row["model_id"], row["group_id"], row["frame_1based"]): row
        for row in seal["predictions"]
    }
    assert len(predictions) == len(seal["predictions"]) == 288
    assert len(report["group_evaluations"]) == 12
    seen, group_losses = set(), {model: [] for model in expected_mse}
    for group in report["group_evaluations"]:
        assert group["complete"] is True and len(group["frames"]) == 24
        errors = []
        for frame in group["frames"]:
            key = frame["model_id"], frame["group_id"], frame["frame_1based"]
            assert key not in seen
            seen.add(key)
            prediction = predictions[key]
            assert prediction["status"] == "ok"
            assert frame["predicted_mean"] == prediction["mean"]
            error = prediction["mean"] - frame["observed_mean"]
            assert frame["squared_error"] == pytest.approx(error * error, abs=1e-14)
            errors.append(error * error)
        group_losses[group["model_id"]].append(math.fsum(errors) / len(errors))
    assert seen == set(predictions)
    for model, expected in expected_mse.items():
        assert len(group_losses[model]) == 2
        replayed = math.fsum(group_losses[model]) / 2
        assert replayed == pytest.approx(expected, rel=1e-10, abs=1e-14)
        assert summary["models"][model]["mean_mse"] == pytest.approx(replayed, abs=1e-14)
    assert summary["selected_development_model"] == report["selected_development_model"] == "transient_response"
    assert summary["parameter_bound_warnings"]["transient_response"]["A"] == "lower"
    assert summary["strong_claim_authorized"] is report["strong_claim_authorized"] is False
    assert summary["final_test_defined"] is report["final_test_defined"] is False


def _approved_storage_locations():
    native = "outputs/native_training_run_02"
    source = "data/holdout_transfer/source_index.json"
    prior = "outputs/ecmodel_prior_port_01/ecmodel_isoprenoid_prior.json"
    fix = (
        "outputs/native_population_development/native_training_01_phase_c_01_scoring_fix_v2/"
        "reusable_interface_fix_verification.json"
    )
    locations = {
        (NATIVE_CHECKPOINT, "/training_manifest/root"): ".",
        (NATIVE_MANIFEST, "/root"): ".",
        (f"{native}/contract.json", "/root"): ".",
        (source, "/original_assets/0/local_absolute_path"): "data/carotenoid/elizondo2025_steady_states.tsv",
        (source, "/original_assets/1/local_absolute_path"): "data/carotenoid/SOURCE.md",
        ("outputs/biology_readiness_01/readiness.json", "/claims/0/assessment/verification_root"): ".",
        (prior, "/source_code/0/path"): "scripts/ecmodel_isoprenoid_prior.py",
        (prior, "/source_code/1/path"): "src/ystwin/pathway/proteome.py",
        (prior, "/source_code/2/path"): "src/ystwin/pathway/enzyme_capacity.py",
        (prior, "/source_model/path"): "data/gem/ecYeastGEM_batch.xml.gz",
        (fix, "/authority_pins/authority_root"): ".",
    }
    for variant in ("nad", "oxygen_peroxide", "oxygen_water"):
        locations[f"outputs/heldout_transfer_{variant}/predictions.json", "/task_manifest/root"] = "."
        locations[f"outputs/heldout_transfer_{variant}/task_manifest.json", "/root"] = "."
        locations[f"outputs/heldout_score_{variant}/score.json", "/label_manifest/root"] = "."
    return locations


def test_storage_transforms_are_the_exact_schema_allowlist_not_generic_scientific_redaction(
    portable_copy,
):
    observed = {}
    for record in portable_copy.manifest["records"]:
        for transform in record["transforms"]:
            identity = record["origin_path"], transform["json_pointer"]
            assert identity not in observed
            assert transform["operation"] == "replace_storage_location"
            observed[identity] = transform["replacement"]
    assert observed == _approved_storage_locations()


@pytest.mark.parametrize("relative", [
    "src/ystwin/mech/hog.py", "data/hog2013/sources.json",
])
def test_relocated_cli_checks_active_code_and_input_bytes_before_replay(portable_copy, relative):
    _materialize_public_runtime(portable_copy)
    target = portable_copy.root / relative
    target.write_bytes(target.read_bytes() + b"\n")
    completed = _run_relocated_cli(portable_copy, [
        "--root", str(portable_copy.root), "--manifest", portable_copy.manifest_path,
        "--manifest-sha256", portable_copy.pin, "--output", "must-not-exist.json",
    ])
    assert completed.returncode != 0
    assert "SHA-256 mismatch" in completed.stderr
    assert not (portable_copy.root / "must-not-exist.json").exists()


def test_portable_export_test_suite_does_not_hide_private_fixture_dependencies_as_skips(
    portable_copy,
):
    """Run the actual exporter suite in a public-only checkout; inspect structured results."""
    _materialize_public_runtime(portable_copy)
    destination = portable_copy.root / PUBLIC_RELEASE
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(portable_copy.directory, destination)
    # The current exporter suite and CLI are the subjects here, not replay dependencies.
    for relative in ("tests/test_portable_evidence.py", "scripts/export_portable_evidence.py"):
        destination = portable_copy.root / relative
        destination.parent.mkdir(exist_ok=True)
        destination.write_bytes((ROOT / relative).read_bytes())
    test_path = portable_copy.root / "tests/test_portable_evidence.py"
    result_path = portable_copy.root / "public-test-results.xml"
    completed = subprocess.run(
        [sys.executable, "-B", "-m", "pytest", "--noconftest", "-p", "no:cacheprovider",
         "-q", "--basetemp", str(portable_copy.root / "pytest-temporary"),
         "--junitxml", str(result_path), str(test_path)],
        cwd=portable_copy.root, text=True, capture_output=True, timeout=120,
        env={"PYTHONPATH": str(portable_copy.root / "src"),
             "HOME": str(portable_copy.root / "isolated-user"), "LANG": "C.UTF-8",
             "PYTHONDONTWRITEBYTECODE": "1", "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "CI": "1",
             "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    suites = ET.fromstring(result_path.read_bytes()).findall("testsuite")
    assert suites
    skipped = sum(int(suite.attrib["skipped"]) for suite in suites)
    assert skipped == 0, (
        f"The portable exporter test suite skipped {skipped} cases when original private "
        "fixtures were unavailable. A public-only green run must not suppress those checks."
    )
