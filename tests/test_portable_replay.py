from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from ystwin.analysis import biology_learning_audit as audit
from ystwin.analysis import frozen_runtime
from ystwin.analysis import portable_replay as replay
from ystwin.analysis.portable_evidence import load_portable_evidence

from _frozen_runtime_helpers import run_frozen_python as _frozen_process


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def public_bundle():
    return load_portable_evidence(ROOT, replay.DEFAULT_MANIFEST_PATH, replay.DEFAULT_MANIFEST_SHA256)


@pytest.fixture(scope="module")
def relocated_root(tmp_path_factory, public_bundle):
    root = tmp_path_factory.mktemp("public_replay_relocated").resolve(strict=True)
    with frozen_runtime.materialize_frozen_runtime(ROOT, operation="audit") as snapshot:
        shutil.copytree(snapshot.root, root, dirs_exist_ok=True)
    support = {
        "scripts/run_native_reconciliation.py": "36a7869cabb76822d430c5a71131f09a238601c00d368569d0fd06d43ba39449",
        "src/ystwin/fba/native_reconciliation.py": "5a817180d0d1ecb09e898f884b7e187d7efdf29fad4f50423b56358296a24f53",
    }
    for relative, data in frozen_runtime._git_blobs(ROOT, frozen_runtime.FROZEN_CODE_REF, support).items():
        assert hashlib.sha256(data).hexdigest() == support[relative]
        destination = root / relative
        with destination.open("xb") as handle:
            handle.write(data)
    for record in public_bundle.manifest["records"]:
        if record["transforms"]:
            assert not (root / record["origin_path"]).exists()
    assert not (root / ".ystwin.env").exists()
    return root


@pytest.fixture(scope="module")
def replayed(relocated_root):
    return frozen_runtime.replay_historical_evidence(relocated_root, source_repository=ROOT)


def _runner():
    spec = importlib.util.spec_from_file_location("portable_native_reconciliation", ROOT / "scripts/run_native_reconciliation.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def _mutable_root(tmp_path, relocated_root):
    root = tmp_path / "mutated_public_copy"
    shutil.copytree(relocated_root, root)
    return root


def _save_json(path, payload):
    data = (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def _change_public_payload(root, origin, change):
    manifest_path = root / replay.DEFAULT_MANIFEST_PATH
    manifest = json.loads(manifest_path.read_bytes())
    record = next(item for item in manifest["records"] if item["origin_path"] == origin)
    path = manifest_path.parent / record["public_path"]
    envelope = json.loads(path.read_bytes())
    change(envelope["payload"])
    record["public_sha256"] = _save_json(path, envelope)
    record["public_size_bytes"] = path.stat().st_size
    digest = hashlib.sha256(replay._json_bytes(envelope["payload"])).hexdigest()
    record["scientific_payload"].update(original_sha256=digest, public_sha256=digest)
    for link in manifest["links"]:
        if link["target_origin_path"] == origin:
            link["target_public_sha256"] = record["public_sha256"]
    return _save_json(manifest_path, manifest)


def test_public_replay_closes_native_and_score_fixtures_without_originals(relocated_root, replayed):
    assert not (relocated_root / replay.CHECKPOINT_PATH).exists()
    assert replayed["kind"] == "portable_replay_of_previously_recorded_content"
    assert replayed["verified"]
    assert not replayed["new_fitting"] and not replayed["new_selection"]
    assert not replayed["independent_test"] and not replayed["biological_validation"]
    assert not replayed["original_serialization_timestamp_verified"]
    integrity = replayed["integrity"]
    assert integrity["input_count"] == 7 and integrity["code_count"] == 20
    assert len(integrity["dependencies"]) == 27
    assert all(row["matches"] for row in integrity["dependencies"])
    assert integrity["runtime"]["versions_match"]
    assert len(integrity["runtime"]["loaded_frozen_modules"]) == 17
    assert integrity["original_checkpoint_bytes"]["status"] == "not_checked"
    assert integrity["original_checkpoint_bytes"]["matches"] is None
    assert integrity["checkpoint_reference"]["sha256"] != replay.CHECKPOINT_SHA256
    assert integrity["model_sha256"] == replay.MODEL_SHA256
    hog = replayed["native_hog"]
    assert len(hog["rows"]) == 82
    assert hog["training_nrmse"] == pytest.approx(0.11016896451452807, rel=1e-10)
    assert not hog["parameters_fitted"] and not hog["observation_gains_fitted"]
    assert {row["unit"] for row in hog["rows"]} == {"mol/L"}
    exchange = replayed["native_exchange"]["computed_report"]
    assert exchange["n_readouts"] == 27 and exchange["n_scored"] == 18
    assert sum(row["prediction"] is None for row in exchange["readouts"]) == 9
    assert {row["growth_rate_per_h"] for row in exchange["readouts"]} == {0.15, 0.28, 0.35}
    assert {row["variant"] for row in replayed["transfer_scores"]} == set(replay._VARIANTS)
    for row in replayed["transfer_scores"]:
        score = row["computed_metrics"]
        assert row["score_verified"] and not row["prediction_kernel_reexecuted"]
        assert score["n_total"] == 6 and score["n_scorable"] == score["n_contained"] == 3
        assert score["status_counts"]["failed"] == 3
        assert score["containment_fraction_total"] == 0.5 and score["unit"] == "mg/gDW"
        assert sum(item["width"] is None for item in score["rows"]) == 3
    population = replayed["population_metrics"]
    assert population["verified"] and population["prediction_count"] == population["evaluated_frame_rows"] == 288
    assert population["unit"] == "dimensionless_ratio" and len(population["models"]) == 6
    assert all(len(row["groups"]) == 2 for row in population["models"].values())
    assert population["models"]["transient_response"]["mean_mse"] == pytest.approx(0.006773738644263623, rel=1e-10)
    assert population["raw_source_verification"] == "not_checked"
    assert not population["selection_reexecuted"] and not population["prediction_kernel_reexecuted"]
    assert not population["strong_claim_authorized"] and not population["final_test_defined"]
    data = json.dumps(replayed, allow_nan=False)
    assert str(ROOT) not in data and str(relocated_root) not in data


def test_portable_mode_never_calls_fitters_or_the_transfer_lp(relocated_root):
    process = _frozen_process(relocated_root, """
        from contextlib import ExitStack
        from unittest.mock import patch
        from ystwin.analysis import blind_transfer, hog_learning, native_physiology
        from ystwin.analysis import portable_replay as replay, biology_learning_audit as audit
        from scipy import optimize

        def forbidden(*args, **kwargs):
            raise AssertionError("portable replay cannot fit, select or silently rerun the transfer LP")

        with ExitStack() as stack:
            for owner, name in ((hog_learning, "fit_native_hog_parameters"),
                                (native_physiology, "fit_native_exchange_model"),
                                (native_physiology.NativeExchangeModel, "from_dict"),
                                (blind_transfer, "predict_balanced_envelope"),
                                (optimize, "least_squares")):
                stack.enter_context(patch.object(owner, name, forbidden))
            assert replay.replay_portable_evidence(root)["verified"]
            ledger = audit.collect_learning_evidence(root, artifact_source="portable")
            assert audit.require_claim(ledger, "portable_scientific_content", root=root,
                                       artifact_source="portable")["authorized"]
    """)
    assert process.returncode == 0, process.stderr


def test_audit_and_reconciliation_use_real_portable_sources_at_relocated_root(relocated_root):
    process = _frozen_process(relocated_root, """
        import runpy
        from ystwin.analysis import biology_learning_audit as audit, portable_replay as replay

        ledger = audit.collect_learning_evidence(root, artifact_source="portable")
        result = audit.require_claim(ledger, "portable_scientific_content", root=root, artifact_source="portable")
        assert result["authorized"] and result["verification_root"] == "."
        assert ledger["evidence_predicates"]["checkpoint_integrity"]["value"] is None
        for claim in ("integrity_established", "conditional_parameter_estimation", "empirical_interpolation",
                      "supplied_kinetic_law", "independent_native_validation", "independent_mechanistic_evidence"):
            try:
                audit.require_claim(ledger, claim, root=root, artifact_source="portable")
            except ValueError as exc:
                assert "unsupported claim" in str(exc)
            else:
                raise AssertionError("portable evidence cannot authorize " + claim)
        runner = runpy.run_path(str(root / "scripts/run_native_reconciliation.py"))
        scenarios, receipt = runner["_frozen_native_context"](root, artifact_source="portable")
        assert {row["nacl_molar"] for row in scenarios.values()} == {0.0, 0.4}
        assert receipt["original_checkpoint_bytes"]["status"] == "not_checked"
        assert receipt["original_checkpoint_bytes_verified"] is False
        assert receipt["all_frozen_inputs_and_code_verified"]
        assert receipt["checkpoint_sha256"] == replay.CHECKPOINT_SHA256
        assert receipt["sha256"] != replay.CHECKPOINT_SHA256
    """)
    assert process.returncode == 0, process.stderr


def test_explicit_original_mode_does_not_fall_back_to_public_bytes(relocated_root):
    with pytest.raises(ValueError, match="nonsymlink file"):
        audit.collect_learning_evidence(relocated_root, artifact_source="original")
    with pytest.raises(FileNotFoundError):
        _runner()._frozen_native_context(relocated_root, artifact_source="original")
    for operation in (audit.collect_learning_evidence, _runner()._frozen_native_context):
        with pytest.raises(ValueError, match="artifact_source"):
            operation(relocated_root, artifact_source="original", manifest_path=replay.DEFAULT_MANIFEST_PATH)
        with pytest.raises(ValueError, match="artifact_source"):
            operation(relocated_root, artifact_source="automatic")


def test_original_verification_only_claims_actual_matching_available_bytes(public_bundle, relocated_root, tmp_path):
    checked = replay._check_private_originals(public_bundle, relocated_root)
    checkpoint = next(row for row in checked if row["path"] == replay.CHECKPOINT_PATH)
    assert checkpoint["status"] == "not_checked" and checkpoint["actual_sha256"] is None
    native_model = next(row for row in checked if row["path"] == "data/hog2013/model_wt.xml")
    assert native_model["status"] == "verified" and native_model["matches"] is True
    assert native_model["actual_sha256"] == native_model["expected_sha256"]
    missing = replay._check_private_originals(public_bundle, tmp_path / "unavailable_originals")
    assert len(missing) == len(public_bundle.manifest["records"])
    assert all(row["status"] == "not_checked" and row["matches"] is None for row in missing)
    private = tmp_path / "not_an_original"
    path = private / replay.CHECKPOINT_PATH
    path.parent.mkdir(parents=True)
    _save_json(path, public_bundle.fetch(replay.CHECKPOINT_PATH).payload)
    with pytest.raises(ValueError, match="original byte digest mismatch"):
        replay._check_private_originals(public_bundle, private)
    with pytest.raises(ValueError, match="retained frozen checkpoint"):
        audit.collect_learning_evidence(private, artifact_source="original")


def test_manifest_and_export_bytes_are_externally_pinned(relocated_root, replayed, tmp_path):
    with pytest.raises(ValueError, match="externally retained"):
        replay.load_portable_native_checkpoint(relocated_root, manifest_path=replay.DEFAULT_MANIFEST_PATH)
    with pytest.raises(ValueError, match="manifest SHA-256"):
        replay.load_portable_native_checkpoint(relocated_root, expected_sha256="0" * 64)
    root = _mutable_root(tmp_path, relocated_root)
    path = root / replayed["integrity"]["checkpoint_reference"]["path"]
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="public artifact SHA-256"):
        replay.load_portable_native_checkpoint(root)


@pytest.mark.parametrize("field,value", [("operation", "replace_scientific_value"),
                                         ("json_pointer", "/model/parameters/hog.kv6_1"),
                                         ("replacement", "unapproved_root")])
def test_unknown_or_scientific_transforms_are_rejected_even_with_a_recomputed_manifest_pin(relocated_root, tmp_path, field, value):
    root = _mutable_root(tmp_path, relocated_root)
    path = root / replay.DEFAULT_MANIFEST_PATH
    manifest = json.loads(path.read_bytes())
    record = next(item for item in manifest["records"] if item["origin_path"] == replay.CHECKPOINT_PATH)
    record["transforms"][0][field] = value
    digest = _save_json(path, manifest)
    with pytest.raises(ValueError, match="transform|redaction"):
        replay.load_portable_native_checkpoint(root, manifest_path=replay.DEFAULT_MANIFEST_PATH, expected_sha256=digest)


@pytest.mark.parametrize("path", ["pyproject.toml", "src/ystwin/mech/hog.py", "data/hog2013/native_training_protocol.json",
                                  "data/physiology/chemostatData_VanHoek1998.tsv"])
def test_changed_relative_kernel_or_source_bytes_fail_closed(relocated_root, tmp_path, path):
    root = _mutable_root(tmp_path, relocated_root)
    file = root / path
    file.write_bytes(file.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        replay.load_portable_native_checkpoint(root)


@pytest.mark.parametrize("replacement", ["../unapproved.py", "https://example.invalid/unapproved.py",
                                         "data/native_law_v2/granados/sealed/raw_mixed/unapproved.py",
                                         "unexported/active_input.py"])
def test_active_code_cannot_be_reclassified_as_a_lineage_only_reference(public_bundle, relocated_root, replacement, monkeypatch):
    checkpoint = public_bundle.fetch(replay.CHECKPOINT_PATH).payload
    checkpoint["training_manifest"]["code"][0]["path"] = replacement
    reads = []
    original = replay._relative_file

    def record(root, path):
        reads.append(path)
        return original(root, path)

    monkeypatch.setattr(replay, "_relative_file", record)
    with pytest.raises(ValueError, match="executable dependency roster"):
        replay._verify_native_checkpoint(relocated_root, checkpoint)
    assert reads == []


def test_imported_runtime_kernel_is_checked_not_only_the_consumer_copy(relocated_root, tmp_path):
    changed = _mutable_root(tmp_path, relocated_root)
    path = changed / "src/ystwin/mech/hog.py"
    path.write_bytes(path.read_bytes() + b"\n")
    process = _frozen_process(changed, """
        from ystwin.mech import hog
        from ystwin.analysis import portable_replay as replay

        assert Path(hog.__file__) == root / "src/ystwin/mech/hog.py"
        replay.load_portable_native_checkpoint(Path(sys.argv[2]))
    """, str(relocated_root))
    assert process.returncode != 0
    assert "running frozen module differs from pinned code: ystwin.mech.hog" in process.stderr


def test_frozen_probe_isolates_real_imports_and_guard_state(relocated_root, public_bundle, monkeypatch):
    frozen = public_bundle.fetch(replay.CHECKPOINT_PATH).payload
    dependency = next(row for row in frozen["training_manifest"]["code"]
                      if row["path"] == "src/ystwin/fba/dynamic_rates.py")
    monkeypatch.setenv("YSTWIN_HOG_REFERENCE", "/not/a/frozen/source")
    monkeypatch.setenv("GIT_WORK_TREE", str(ROOT))
    process = _frozen_process(relocated_root, """
        import hashlib
        import importlib.metadata
        import os
        import platform
        import pytest

        # These names belong to the probe, not the bootstrap read/import guards.
        original = lambda *args: None
        source = "probe-local state"
        from ystwin.fba import dynamic_rates

        request = json.load(sys.stdin)
        assert os.getpid() != request["parent_pid"]
        assert sys.flags.isolated and sys.dont_write_bytecode
        assert Path(dynamic_rates.__file__) == root / request["dependency"]["path"]
        assert hashlib.sha256(Path(dynamic_rates.__file__).read_bytes()).hexdigest() == request["dependency"]["sha256"]
        expected = request["environment"]
        assert platform.python_version() == expected["python"]
        assert {name: importlib.metadata.version(name) for name in expected["packages"]} == expected["packages"]
        assert not any(key.startswith(("YSTWIN_", "GIT_")) for key in os.environ)
        with pytest.raises(PermissionError, match="current checkout is unavailable"):
            Path(sys.argv[2]).read_bytes()
    """, str(ROOT / dependency["path"]), payload={
        "parent_pid": os.getpid(), "dependency": dependency,
        "environment": frozen["model"]["config"]["environment"],
    })
    assert process.returncode == 0, process.stderr


@pytest.mark.parametrize("operation", ["checkpoint", "replay", "audit", "reconciliation"])
def test_current_process_cannot_authorize_a_materialized_frozen_tree(relocated_root, public_bundle, operation):
    from ystwin.fba import dynamic_rates

    frozen = public_bundle.fetch(replay.CHECKPOINT_PATH).payload
    dependency = next(row for row in frozen["training_manifest"]["code"]
                      if row["path"] == "src/ystwin/fba/dynamic_rates.py")
    assert hashlib.sha256((relocated_root / dependency["path"]).read_bytes()).hexdigest() == dependency["sha256"]
    current = Path(dynamic_rates.__file__)
    assert current == ROOT / dependency["path"]
    assert hashlib.sha256(current.read_bytes()).hexdigest() != dependency["sha256"]
    calls = {
        "checkpoint": lambda: replay.load_portable_native_checkpoint(relocated_root),
        "replay": lambda: replay.replay_portable_evidence(relocated_root),
        "audit": lambda: audit.collect_learning_evidence(relocated_root, artifact_source="portable"),
        "reconciliation": lambda: _runner()._frozen_native_context(relocated_root, artifact_source="portable"),
    }
    with pytest.raises(ValueError, match="running frozen module differs from pinned code"):
        calls[operation]()
    assert Path(dynamic_rates.__file__) == current


def test_runtime_versions_cannot_be_asserted_by_metadata(relocated_root, monkeypatch):
    monkeypatch.setattr(replay.platform, "python_version", lambda: "0.0.0")
    with pytest.raises(ValueError, match="runtime versions"):
        replay.load_portable_native_checkpoint(relocated_root)


@pytest.mark.parametrize("mutation", ["parameter", "unit", "role"])
def test_native_scientific_model_digest_is_independently_retained(public_bundle, relocated_root, mutation):
    checkpoint = public_bundle.fetch(replay.CHECKPOINT_PATH).payload
    changed = deepcopy(checkpoint)
    if mutation == "parameter":
        changed["model"]["parameters"]["hog.kv6_1"] *= 2
    elif mutation == "unit":
        changed["model"]["config"]["unit"] = "mg/L"
    else:
        changed["model"]["config"]["host_roles"]["oxygen_exchange_ids"] = ["r_1714"]
    with pytest.raises(ValueError, match="scientific model digest"):
        replay._verify_native_checkpoint(relocated_root, changed)


def test_recomputed_metrics_cannot_be_replaced_by_changed_pass_metadata(relocated_root, tmp_path):
    root = _mutable_root(tmp_path, relocated_root)
    digest = _change_public_payload(root, "outputs/heldout_score_nad/score.json",
                                    lambda score: score.update(mean_width_on_scorable=0.0))
    process = _frozen_process(root, """
        from ystwin.analysis import portable_replay as replay
        replay.replay_portable_evidence(root, manifest_path=replay.DEFAULT_MANIFEST_PATH, expected_sha256=sys.argv[2])
    """, digest)
    assert process.returncode != 0 and "recomputed interval score" in process.stderr


@pytest.mark.parametrize("mutation", ["metric", "frame_inventory"])
def test_population_aggregate_replay_rejects_changed_arithmetic_or_missing_frames(relocated_root, tmp_path, mutation):
    root = _mutable_root(tmp_path, relocated_root)
    run = "outputs/native_population_development/native_training_01_phase_c_01_scoring"
    if mutation == "metric":
        origin = f"{run}/scoring_summary.json"

        def change(payload):
            payload["models"]["transient_response"]["mean_mse"] = 0.0
    else:
        origin = f"{run}/development_report.json"

        def change(payload):
            payload["group_evaluations"][0]["frames"].pop()
    digest = _change_public_payload(root, origin, change)
    _, _, bundle = replay._load_bundle(root, replay.DEFAULT_MANIFEST_PATH, digest)
    with pytest.raises(ValueError, match="population.*differs"):
        replay._replay_population_metrics(bundle)


def test_outcome_roles_cannot_alias_native_input_provenance(relocated_root, tmp_path, public_bundle):
    root = _mutable_root(tmp_path, relocated_root)
    native_source = public_bundle.fetch(replay.CHECKPOINT_PATH).payload["training_manifest"]["inputs"][0]["source"]

    def overlap(score):
        source = score["label_manifest"]["inputs"][0]["source"]
        source.update(authority=native_source["authority"], source_id=native_source["source_id"])

    digest = _change_public_payload(root, "outputs/heldout_score_nad/score.json", overlap)
    with pytest.raises(ValueError, match="outcome labels overlap"):
        replay.replay_portable_evidence(root, manifest_path=replay.DEFAULT_MANIFEST_PATH, expected_sha256=digest)


@pytest.mark.parametrize("relative", ["src/ystwin/mech/hog.py", "data/hog2013/model_wt.xml"])
def test_current_strict_api_rejects_symlink_aliases_with_unchanged_bytes(relocated_root, tmp_path, relative):
    root = _mutable_root(tmp_path, relocated_root)
    path = root / relative
    target = path.with_name(path.name + ".original")
    path.rename(target)
    path.symlink_to(target.name)
    with pytest.raises(ValueError, match="nonsymlink file"):
        replay.load_portable_native_checkpoint(root)


def test_current_strict_api_rejects_distinct_roles_on_one_physical_code_file(relocated_root, tmp_path):
    root = _mutable_root(tmp_path, relocated_root)
    path = root / "src/ystwin/fba/__init__.py"
    path.rename(path.with_name("original_init.py"))
    os.link(root / "src/ystwin/analysis/__init__.py", path)
    with pytest.raises(ValueError, match="duplicate physical file in allowlist"):
        replay.load_portable_native_checkpoint(root)


def test_cli_replays_from_relocated_public_root_and_refuses_old_record_overwrites(relocated_root, tmp_path):
    environment = {key: value for key, value in os.environ.items() if not key.startswith("YSTWIN_")}
    environment.update(PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
    output = tmp_path.resolve(strict=True) / "new_replay.json"
    command = [sys.executable, str(ROOT / "scripts/replay_portable_evidence.py"),
               "--root", str(relocated_root), "--historical", "--source-repository", str(ROOT),
               "--output", str(output)]
    process = subprocess.run(command, cwd=relocated_root, env=environment, capture_output=True, text=True, timeout=90)
    assert process.returncode == 0, process.stderr
    report = json.loads(process.stdout)
    assert report == json.loads(output.read_bytes())
    assert report["verified"] and report["integrity"]["original_checkpoint_bytes"]["status"] == "not_checked"
    before = output.read_bytes()
    again = subprocess.run(command, cwd=relocated_root, env=environment, capture_output=True, text=True, timeout=30)
    assert again.returncode != 0 and "write-once" in again.stderr
    assert output.read_bytes() == before
