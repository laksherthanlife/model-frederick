from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ystwin.analysis import native_population_development as development

from test_population_historical_sources import historical_population_root as historical_population_root


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "data/native_law_v2/development_protocol.json"


@pytest.fixture(scope="module")
def protocol():
    return development.load_protocol(PROTOCOL_PATH)


def group_spec(replicate=1, rows=3):
    return {
        "group_id": f"1b38366b-7473-4ddf-9cf9-26e625533686:/rep{replicate}",
        "source_pointer": f"/rep{replicate}",
        "source_rows": rows,
        "split": "train" if replicate <= 3 else "development",
    }


def source_cell(group, row, frame, ratio, median=1.0, time=None):
    prefix = 28 <= frame <= 47
    return {
        "group_id": group["group_id"],
        "split": group["split"],
        "window": "prefix" if prefix else "response",
        "cell_index": row,
        "frame_1based": frame,
        "relative_time_min": float((frame - 49) * 2.5) if time is None else time,
        "max5": None if ratio is None else float(ratio * median),
        "median": float(median),
        "time_pointer": f"{group['source_pointer']}/general/times/{row}/{frame - 1}",
        "max5_pointer": f"{group['source_pointer']}/GFP/max5/{row}/{frame - 1}",
        "median_pointer": f"{group['source_pointer']}/GFP/median/{row}/{frame - 1}",
    }


def cells_for(group, response=True):
    frames = [*range(28, 48), *(range(50, 74) if response else ())]
    return [source_cell(group, row, frame, 2.0 + row + (frame - 28) * 0.1, median=1.0 + row) for row in range(group["source_rows"]) for frame in frames]


def synthetic_prepared(protocol):
    prefixes, frames, coordinates = [], [], []
    for index, spec in enumerate(development.group_specs(protocol)):
        mu = 1.2 + 0.2 * index
        prefixes.append({
            "group_id": spec["group_id"], "origin_1based": 49, "eligible_cell_indices": [0, 1],
            "mu": mu, "s": 0.0, "input_projection_sha256": "a" * 64,
            "operator_id": development.PREFIX_OPERATOR,
        })
        for frame in range(50, 74):
            time = (frame - 49) * 2.5
            row = {"group_id": spec["group_id"], "frame_1based": frame, "time_min": time}
            if index < 3:
                row.update(observed_mean=mu - 0.35 * (1.0 - np.exp(-time / 12.0)), prefix_cells=2, valid_cells=2, source_missing_cells=0, invalid_cells=0, cell_sample_sd=0.0)
                frames.append(row)
            else:
                coordinates.append(row)
    return development.finalize_prepared(prefixes, frames, coordinates, protocol)


def test_new_protocol_schema_and_mean_estimand_are_not_old_per_cell_mae(protocol):
    assert protocol["estimand"]["quantity"] == development.QUANTITY
    assert protocol["estimand"]["unit"] == "dimensionless_ratio"
    assert protocol["scoring_and_selection"]["selection_rule"]
    assert len(protocol["export_contract"]["schemas"]["$defs"]) >= 15
    assert [row["model_id"] for row in protocol["model_recipe"]["canonical_families"]] == list(development.MODEL_IDS)
    assert development.SCIENTIFIC_SCOPE == "retrospective_source_cohort_development"


def test_prefix_is_mean_of_individual_ratios_and_mean_frame_sample_sd():
    group = group_spec(rows=2)
    cells = [source_cell(group, row, frame, 2.0 + 2.0 * row + 10.0 * (frame - 28), median=1.0 + 2.0 * row) for row in range(2) for frame in range(28, 48)]
    summary, audit = development.derive_prefix(group, cells, "a" * 64)
    assert summary["mu"] == pytest.approx(98.0)
    assert summary["s"] == pytest.approx(np.sqrt(2.0))
    assert summary["eligible_cell_indices"] == [0, 1]
    assert audit["excluded_cells"] == []
    assert summary["mu"] != pytest.approx(np.mean([row["max5"] for row in cells]) / np.mean([row["median"] for row in cells]))


def test_prefix_cohort_is_fixed_and_response_poison_cannot_change_it():
    group = group_spec()
    cells = cells_for(group)
    missing = next(row for row in cells if row["cell_index"] == 2 and row["frame_1based"] == 30)
    missing["median"] = None
    before, audit = development.derive_prefix(group, cells, "a" * 64)
    poisoned = copy.deepcopy(cells)
    for row in poisoned:
        if row["window"] == "response":
            row["max5"] = 1e30
            row["median"] = 0.0
    after, _ = development.derive_prefix(group, poisoned, "a" * 64)
    assert before == after
    assert before["eligible_cell_indices"] == [0, 1]
    assert audit["excluded_cells"][0]["cell_index"] == 2
    assert audit["excluded_cells"][0]["frames"][0]["frame_1based"] == 30


def test_zero_prefix_dispersion_is_allowed_and_no_prefix_imputation():
    group = group_spec(rows=2)
    cells = [source_cell(group, row, frame, 3.0) for row in range(2) for frame in range(28, 48)]
    summary, _ = development.derive_prefix(group, cells, "a" * 64)
    assert summary["s"] == 0.0
    cells[0]["max5"] = None
    with pytest.raises(development.DevelopmentError, match="prefix"):
        development.derive_prefix(group, cells, "a" * 64)


def test_response_times_do_not_depend_on_intensity_validity_or_future_cell_count():
    group = group_spec()
    cells = cells_for(group)
    summary, _ = development.derive_prefix(group, cells, "a" * 64)
    rows = [row for row in cells if row["window"] == "response"]
    first = [row for row in rows if row["frame_1based"] == 50]
    first[0]["relative_time_min"] = 1.0
    first[1]["relative_time_min"] = 4.0
    first[2]["relative_time_min"] = 7.0
    first[0]["max5"] = None
    first[1]["median"] = 0.0
    frames = development.derive_frames(group, rows, rows, summary)
    row = frames[0]
    assert row["time_min"] == 4.0
    assert row["valid_cells"] == 1
    assert row["source_missing_cells"] == 1
    assert row["invalid_cells"] == 1
    assert row["cell_sample_sd"] is None
    assert row["observed_mean"] == pytest.approx(first[2]["max5"] / first[2]["median"])
    assert row["prefix_cells"] == 3


def test_scale_gauge_leaves_ratios_cohorts_targets_and_dispersion_unchanged():
    group = group_spec()
    original = cells_for(group)
    scaled = copy.deepcopy(original)
    for row in scaled:
        row["max5"] *= 2.0
        row["median"] *= 2.0
    a, _ = development.derive_prefix(group, original, "a" * 64)
    b, _ = development.derive_prefix(group, scaled, "a" * 64)
    assert a == b
    ra = [r for r in original if r["window"] == "response"]
    rb = [r for r in scaled if r["window"] == "response"]
    assert development.derive_frames(group, ra, ra, a) == development.derive_frames(group, rb, rb, b)


def test_source_identity_rows_windows_and_pointers_are_not_inferred_from_names():
    group = group_spec(rows=2)
    cells = cells_for(group)
    development.validate_group_records(group, cells, ("prefix", "response"))
    for change in ({"cell_index": 5}, {"frame_1based": 49}, {"max5_pointer": "/rep6/GFP/max5/0/27"}, {"group_id": group_spec(6)["group_id"]}):
        altered = copy.deepcopy(cells)
        altered[0].update(change)
        with pytest.raises(development.DevelopmentError):
            development.validate_group_records(group, altered, ("prefix", "response"))
    with pytest.raises(development.DevelopmentError):
        development.validate_group_records(group, cells[:-1] + [cells[0]], ("prefix", "response"))


def test_strict_json_rejects_duplicate_keys_and_nonfinite_tokens():
    for content in ('{"x":1,"x":2}', '{"x":NaN}', '{"x":1e999}'):
        with pytest.raises(development.DevelopmentError):
            development.strict_json(content)


def test_declared_output_map_is_used_before_loss_and_reports_binding():
    raw, mean, binds = development.evaluate_mean("prefix_offset", {"a": -3.0}, np.array([1.0, 4.0]), np.array([2.5, 60.0]))
    np.testing.assert_array_equal(raw, [-2.0, 1.0])
    np.testing.assert_array_equal(mean, [0.0, 1.0])
    np.testing.assert_array_equal(binds, [True, False])
    assert development.mean_mse(mean, np.array([2.0, 2.0])) == 2.5
    assert development.evaluate_mean("training_constant", {"c": 5.0}, np.array([1.0]), np.array([2.5]))[1][0] == 5.0


@pytest.mark.parametrize("model_id,parameters", [
    ("prefix_persistence", {}), ("training_constant", {"c": 2.0}),
    ("prefix_offset", {"a": -0.2}), ("prefix_affine_time", {"a": -0.1, "b": 0.3}),
    ("first_order_response", {"A": -0.4, "tau": 17.0}),
    ("transient_response", {"A": 0.5, "tau_off": 31.0, "tau_on": 8.0}),
])
def test_time_units_are_equivariant_for_every_canonical_family(model_id, parameters):
    time = np.array([2.5, 12.0, 60.0])
    mu = np.array([1.0, 2.0, 3.0])
    scaled_parameters = {key: value * 60.0 if key in ("tau", "tau_on", "tau_off") else value for key, value in parameters.items()}
    expected = development.evaluate_mean(model_id, parameters, mu, time)[1]
    actual = development.evaluate_mean(model_id, scaled_parameters, mu, time * 60.0, reference_time=3600.0)[1]
    np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-8)


def test_training_only_bound_scale_and_iqr_use_all_72_residual_means(protocol):
    prepared = synthetic_prepared(protocol)
    means = {row["group_id"]: row["mu"] for row in prepared["prefix_summaries"]}
    residual = [row["observed_mean"] - means[row["group_id"]] for row in prepared["training_frames"]]
    assert prepared["S_train"] == pytest.approx(np.quantile(residual, 0.75) - np.quantile(residual, 0.25))
    original_h, original_s = prepared["H_train"], prepared["S_train"]
    modified = copy.deepcopy(prepared)
    for prefix in modified["prefix_summaries"][-2:]:
        prefix["mu"] = 1e6
    recomputed = development.finalize_prepared(modified["prefix_summaries"], modified["training_frames"], modified["development_coordinates"], protocol)
    assert recomputed["H_train"] == original_h
    assert recomputed["S_train"] == original_s


def test_all_six_families_and_registered_starts_are_logged(protocol):
    fitted = development.fit_training(synthetic_prepared(protocol), protocol)
    assert list(fitted) == list(development.MODEL_IDS)
    assert len(fitted["prefix_persistence"]["attempts"]) == 1
    assert sum(len(row["attempts"]) for row in fitted.values()) <= 41
    for model_id, result in fitted.items():
        assert result["status"] == "ok"
        assert result["selected_start_index"] is not None
        if model_id != "prefix_persistence":
            assert len(result["attempts"]) == 8
        for index, attempt in enumerate(result["attempts"]):
            assert attempt["seed"] == 201801 + index
            assert 0 <= attempt["loss_evaluations"] <= 2000
            assert set(attempt["initial_parameters"]) == set(development.PARAMETERS[model_id])
            if attempt["status"] == "ok":
                assert attempt["training_mean_mse"] >= 0.0
    assert fitted["first_order_response"]["training_mean_mse"] < 1e-15
    assert fitted["first_order_response"]["parameters"]["tau"] == pytest.approx(12.0, rel=1e-4)


def test_json_object_order_is_not_a_model_dependency(protocol):
    prepared = synthetic_prepared(protocol)
    fitted = development.fit_training(prepared, protocol)
    restored = development.strict_json(development.canonical_bytes(fitted))
    assert development.predict_development(prepared, restored, protocol) == development.predict_development(prepared, fitted, protocol)


def test_failed_family_is_retained_with_all_48_prediction_slots(protocol, monkeypatch):
    prepared = synthetic_prepared(protocol)
    fitted = development.fit_training(prepared, protocol)
    fitted["transient_response"].update(status="failed", parameters=None, selected_start_index=None, reason="synthetic injected fitting failure")
    predictions, audit = development.predict_development(prepared, fitted, protocol)
    assert len(predictions) == 288
    failed = [row for row in predictions if row["model_id"] == "transient_response"]
    assert len(failed) == 48 and all(row["status"] == "failed" and row["mean"] is None for row in failed)
    assert len(audit) == 288
    for row in predictions:
        development.validate_document(protocol, row, "FramePrediction")


def test_incomplete_training_group_never_becomes_a_smaller_training_fit(protocol):
    prepared = synthetic_prepared(protocol)
    prepared["training_frames"][0]["observed_mean"] = None
    with pytest.raises(development.DevelopmentError, match="complete"):
        development.finalize_prepared(prepared["prefix_summaries"], prepared["training_frames"], prepared["development_coordinates"], protocol)


def test_prediction_schema_forbids_wrong_estimand_units_and_extra_hidden_inputs(protocol):
    prepared = synthetic_prepared(protocol)
    fitted = development.fit_training(prepared, protocol)
    row = development.predict_development(prepared, fitted, protocol)[0][0]
    for change in ({"quantity": "individual_cell_trajectory"}, {"unit": "phosphorylation_fraction"}, {"future_valid_cell_count": 50}):
        altered = {**row, **change}
        with pytest.raises(development.DevelopmentError):
            development.validate_document(protocol, altered, "FramePrediction")


def test_conditional_mean_mse_not_mean_cell_error_or_mae_selects_model():
    assert development.mean_mse(np.array([0.0, 0.0]), np.array([0.0, 2.0])) == 2.0
    assert development.mean_mse(np.array([1.0, 1.0]), np.array([0.0, 2.0])) == 1.0
    assert development.select_model({model: {"complete": True, "mean_mse": 1.0} for model in development.MODEL_IDS}) == "prefix_persistence"
    rows = {model: {"complete": True, "mean_mse": 1.0} for model in development.MODEL_IDS}
    rows["transient_response"] = {"complete": False, "mean_mse": None}
    assert development.select_model(rows) is None
    del rows["transient_response"]
    with pytest.raises(development.DevelopmentError):
        development.select_model(rows)


def sealed_synthetic(tmp_path, protocol, failed_family=False):
    prepared = synthetic_prepared(protocol)
    fitted = development.fit_training(prepared, protocol)
    if failed_family:
        fit = fitted["transient_response"]
        fit.update(status="failed", parameters=None, selected_start_index=None, reason="synthetic failed family")
        for attempt in fit["attempts"]:
            attempt.update(status="failed", training_mean_mse=None, reason="synthetic failed start")
    code = tmp_path / "synthetic_executor.py"
    code.write_text("pass\n")
    store = development.ArtifactStore(tmp_path, "outputs/native_population_development/synthetic")
    provenance = {"approval_sha256": "b" * 64, "training_projection_sha256": "a" * 64, "development_inputs_sha256": "a" * 64, "code_artifacts": [development.artifact_reference(tmp_path, code, "execution_code")], "operational_event_reference": "synthetic engineering fixture only"}
    summary = development.seal_training_freeze(protocol, prepared, fitted, provenance, store)
    freeze = development.strict_json((tmp_path / summary["prediction_freeze"]["path"]).read_bytes())
    return prepared, fitted, summary, freeze


def test_freeze_is_exact_schema_complete_non_authorizing_and_exclusive(tmp_path, protocol):
    prepared, fitted, summary, freeze = sealed_synthetic(tmp_path, protocol)
    assert development.validate_document(protocol, freeze, "DevelopmentPredictionFreeze")
    assert development.validate_freeze(protocol, freeze)
    assert len(freeze["attempts"]) == 41
    assert len(freeze["predictions"]) == 288
    assert set(summary["selected_checkpoints"]) == set(development.MODEL_IDS)
    assert freeze["strong_claim_authorized"] is False
    assert freeze["environment"]["scientific_scope"] == development.SCIENTIFIC_SCOPE
    audit = development.strict_json((tmp_path / summary["nonnegative_map_audit"]["path"]).read_bytes())
    assert len(audit["development"]) == 288
    assert len(audit["training"]) == 6 * 72
    with pytest.raises(FileExistsError):
        store = development.ArtifactStore(tmp_path, "outputs/native_population_development/synthetic")
        development.seal_training_freeze(protocol, prepared, fitted, freeze["environment"]["provenance"], store)


@pytest.mark.parametrize("defect", ["prediction", "duplicate_frame", "omitted_family", "source_prefix", "bound"])
def test_freeze_rejects_semantic_tampering_not_just_a_recomputed_hash(tmp_path, protocol, defect):
    _, _, _, freeze = sealed_synthetic(tmp_path, protocol)
    if defect == "prediction":
        freeze["predictions"][0]["mean"] += 0.5
    elif defect == "duplicate_frame":
        freeze["predictions"][1] = freeze["predictions"][0]
    elif defect == "omitted_family":
        freeze["attempts"] = [row for row in freeze["attempts"] if row["model_id"] != "transient_response"]
    elif defect == "source_prefix":
        freeze["prefix_summaries"][3]["input_projection_sha256"] = "c" * 64
    else:
        row = next(row for row in freeze["attempts"] if row["model_id"] == "first_order_response")
        row["parameters"]["tau"] = 100000.0
    with pytest.raises(development.DevelopmentError):
        development.validate_freeze(protocol, freeze)


def observed_development(prepared):
    prefix = {row["group_id"]: row for row in prepared["prefix_summaries"]}
    return [{**row, "observed_mean": prefix[row["group_id"]]["mu"], "prefix_cells": 2, "valid_cells": 2, "source_missing_cells": 0, "invalid_cells": 0, "cell_sample_sd": 0.0} for row in prepared["development_coordinates"]]


def test_scoring_uses_frozen_predictions_with_no_refit_and_keeps_all_families(tmp_path, protocol, monkeypatch):
    prepared, _, summary, freeze = sealed_synthetic(tmp_path, protocol)

    def forbidden(*args, **kwargs):
        raise AssertionError("no development refit")

    monkeypatch.setattr(development, "fit_training", forbidden)
    report = development.build_development_report(protocol, freeze, observed_development(prepared), summary["prediction_freeze"]["sha256"], "c" * 64)
    assert development.validate_document(protocol, report, "DevelopmentReport")
    assert report["selected_development_model"] == "prefix_persistence"
    assert len(report["group_evaluations"]) == 12
    assert all(len(row["frames"]) == 24 for row in report["group_evaluations"])
    assert report["strong_claim_authorized"] is False
    assert report["final_test_defined"] is False


def test_failed_comparator_prevents_a_winner_even_when_another_model_is_exact(tmp_path, protocol):
    prepared, _, summary, freeze = sealed_synthetic(tmp_path, protocol, failed_family=True)
    report = development.build_development_report(protocol, freeze, observed_development(prepared), summary["prediction_freeze"]["sha256"], "c" * 64)
    assert report["status"] == "incomplete_model_execution"
    assert report["selected_development_model"] is None
    assert len(report["group_evaluations"]) == 12
    assert all(row["mean_mse"] is None for row in report["group_evaluations"] if row["model_id"] == "transient_response")


def test_one_missing_response_frame_invalidates_complete_ranking_without_window_shrink(tmp_path, protocol):
    prepared, _, summary, freeze = sealed_synthetic(tmp_path, protocol)
    observed = observed_development(prepared)
    observed[0].update(observed_mean=None, valid_cells=0, source_missing_cells=2, cell_sample_sd=None)
    report = development.build_development_report(protocol, freeze, observed, summary["prediction_freeze"]["sha256"], "c" * 64)
    assert report["status"] == "incomplete_observation_coverage"
    assert report["selected_development_model"] is None
    assert sum(len(row["frames"]) for row in report["group_evaluations"]) == 288


def test_nonfinite_optimizer_failure_is_retained_as_finite_json(protocol, monkeypatch):
    def nonfinite(*args, **kwargs):
        return SimpleNamespace(x=np.array([np.nan, np.nan, np.nan]), fun=np.full(72, np.nan), success=False, status=-1, message="synthetic nonfinite optimizer")

    monkeypatch.setattr(development, "least_squares", nonfinite)
    fit = development.fit_family("transient_response", np.ones(72), np.tile(np.arange(1.0, 25.0), 3), np.ones(72), 1.0)
    assert fit["status"] == "failed" and len(fit["attempts"]) == 8
    assert development.canonical_bytes(fit)
    assert all(attempt["status"] == "failed" and attempt["training_mean_mse"] is None for attempt in fit["attempts"])


def test_blocked_custody_gate_does_not_read_any_projection_or_private_path(tmp_path, protocol, monkeypatch):
    protocol_path = tmp_path / "data/native_law_v2/development_protocol.json"
    protocol_path.parent.mkdir(parents=True)
    protocol_path.write_bytes(development.canonical_bytes(protocol))
    blocker_path = tmp_path / "data/native_law_v2/granados/development/operational_blocker.json"
    blocker_path.parent.mkdir(parents=True)
    blocker = {"document_type": "development_blocker_no_biological_verdict", "protocol_id": development.PROTOCOL_ID, "phase": "operational_approval", "reason_code": "custody_not_approved", "reason": "synthetic independent custody refusal", "affected_group_ids": [], "available_artifacts": [], "operational_event_reference": "synthetic fixture", "strong_claim_authorized": False, "further_exports_permitted": False}
    blocker_path.write_bytes(development.canonical_bytes(blocker))
    original = Path.read_bytes
    opened = []

    def metadata_only(path):
        opened.append(path)
        assert path in {protocol_path, blocker_path}
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", metadata_only)
    status = development.check_custody_gate(tmp_path, protocol_path)
    assert status["ready"] is False
    assert status["measurement_arrays_opened"] is False
    assert "synthetic independent custody refusal" in status["reason"]
    for path in ("data/native_law_v2/granados/sealed/raw_mixed/source.json", "data/native_law_v2/granados/development/private/source_coverage.json", "../outside.json"):
        with pytest.raises(development.DevelopmentError):
            development.safe_path(tmp_path, path)
    assert set(opened) == {protocol_path, blocker_path}


def test_actual_phase_b_manifest_is_data_ready_but_changed_code_needs_refresh(monkeypatch, historical_population_root):
    root = historical_population_root
    manifest = root / "data/native_law_v2/granados/development/release_manifest.json"
    forbidden = root / development.admission_contract()["projection_root"]
    original = Path.read_bytes

    def metadata_only(path):
        assert not path.resolve().is_relative_to(forbidden.resolve())
        assert not path.resolve().is_relative_to((root / "data/native_law_v2/granados/sealed").resolve())
        assert not path.resolve().is_relative_to((root / "data/native_law_v2/granados/development/private").resolve())
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", metadata_only)
    result = development.check_custody_gate(root, root / "data/native_law_v2/development_protocol.json", manifest)
    assert result["data_ready"] is True
    assert result["executor_ready"] is False and result["ready"] is False
    assert result["reason_code"] == "executor_admission_refresh_required"
    assert result["measurement_arrays_opened"] is False
    assert result["provenance"]["phase_b_manifest"]["sha256"] == "3b8c67c79191bf4dec0a1e7e0514b1a145359c30061868bbfdd136f0ee1cd3dd"


def synthetic_source_package(tmp_path, protocol, source_root):
    protocol_path = tmp_path / "data/native_law_v2/development_protocol.json"
    protocol_path.parent.mkdir(parents=True)
    protocol_path.write_bytes(development.canonical_bytes(protocol))
    for name in ("src/ystwin/analysis/native_population_development.py", "scripts/run_native_population_development.py", "tests/test_native_population_development.py"):
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((source_root / name).read_bytes())
    exporter = tmp_path / "scripts/synthetic_exporter.py"
    exporter.write_text("pass\n")
    request_summary = development.write_release_request(tmp_path, protocol_path, "outputs/native_population_development/request", [exporter], "synthetic_requester_not_a_real_release")
    public = tmp_path / "data/native_law_v2/granados/development"
    public.mkdir(parents=True)

    def write_metadata(name, payload, role):
        path = public / name
        path.write_bytes(development.canonical_bytes(payload))
        return development.artifact_reference(tmp_path, path, role)

    approval = {"document_type": "development_custody_decision_not_external_attestation", "protocol_id": development.PROTOCOL_ID, "request_sha256": request_summary["request"]["sha256"], "decision": "approve_development_exports", "custodian_identity": "synthetic_custodian_not_a_real_source_attestation", "operational_event_reference": "synthetic fixture only", "reviewed_source_asset_id": development.ASSET_ID, "reviewed_source_sha256": development.SOURCE_SHA256, "permitted_group_ids": [row["group_id"] for row in development.group_specs(protocol)], "conditions": ["synthetic fixture only"], "limitations": ["no actual source measurement or native approval"], "strong_claim_authorized": False}
    approval_ref = write_metadata("approval.json", approval, "operational_approval")
    audit_ref = write_metadata("synthetic_source_audit.json", {"fixture": "synthetic engineering only"}, "source_binding_audit")
    projections = {}
    for phase in ("training_values", "development_inputs", "development_responses"):
        groups = development.group_specs(protocol)[:3] if phase == "training_values" else development.group_specs(protocol)[3:]
        cells, times = [], []
        for spec in groups:
            for cell in range(spec["source_rows"]):
                baseline = 1.4 + cell * 0.001
                if phase != "development_responses":
                    cells.extend(source_cell(spec, cell, frame, baseline) for frame in range(28, 48))
                for frame in range(50, 74):
                    time = (frame - 49) * 2.5
                    record = source_cell(spec, cell, frame, baseline - 0.35 * (1.0 - np.exp(-time / 12.0)))
                    if phase == "development_inputs":
                        times.append({key: record[key] for key in ("group_id", "cell_index", "frame_1based", "relative_time_min", "time_pointer")})
                    else:
                        cells.append(record)
        projections[phase] = {"document_type": "source_bound_development_projection", "protocol_id": development.PROTOCOL_ID, "phase": phase, "approval_sha256": approval_ref["sha256"], "source_asset_id": development.ASSET_ID, "source_sha256": development.SOURCE_SHA256, "group_ids": [row["group_id"] for row in groups], "source_cells": cells, "source_times": times, "source_binding_audit": audit_ref}
    export_root = tmp_path / development.admission_contract()["projection_root"]
    export_root.mkdir(parents=True)
    references = {}
    for phase, role in (("training_values", "training_projection"), ("development_inputs", "development_inputs"), ("development_responses", "development_responses")):
        path = export_root / f"synthetic_{phase}.json"
        path.write_bytes(development.canonical_bytes(projections[phase]))
        references[role] = development.artifact_reference(tmp_path, path, role)

    def allowlist(phase, freeze_sha=None):
        roles = ("training_projection", "development_inputs") if phase == "training_and_conditioning_export" else ("development_responses",)
        artifacts = [references[role] for role in roles]
        checks = [{"artifact": reference, "source_scalar_count": development.admission_contract()["source_replay_scalar_counts"][reference["role"]], "all_source_scalars_equal": True, "nulls_preserved": True, "source_rows_and_pointers_verified": True, "group_window_inventory_verified": True} for reference in artifacts]
        receipt = {"document_type": "population_development_projection_source_replay", "protocol_id": development.PROTOCOL_ID, "protocol_sha256": request_summary["protocol"]["sha256"], "source_asset_id": development.ASSET_ID, "source_sha256": development.SOURCE_SHA256, "approval_sha256": approval_ref["sha256"], "custodian_identity": approval["custodian_identity"], "phase": phase, "operational_event_reference": "synthetic fixture only", "resolved_blocker_sha256": None, "prediction_freeze_sha256": freeze_sha, "artifact_checks": checks}
        replay_ref = write_metadata(f"{phase}_replay.json", receipt, "source_binding_audit")
        document = {"document_type": "population_development_projection_allowlist", "protocol_id": development.PROTOCOL_ID, "protocol_sha256": request_summary["protocol"]["sha256"], "phase": phase, "request": request_summary["request"], "approval": approval_ref, "configuration": request_summary["configuration"], "source_replay": replay_ref, "artifacts": artifacts, "custodian_identity": approval["custodian_identity"], "operational_event_reference": "synthetic fixture only", "resolved_blocker_sha256": None, "prediction_freeze_sha256": freeze_sha}
        path = public / f"{phase}_allowlist.json"
        path.write_bytes(development.canonical_bytes(document))
        return path

    return protocol_path, references, allowlist


def test_synthetic_phased_export_fit_freeze_release_and_fixed_scoring(tmp_path, protocol, monkeypatch, historical_population_root):
    protocol_path, references, allowlist = synthetic_source_package(tmp_path, protocol, historical_population_root)
    initial = allowlist("training_and_conditioning_export")
    response_path = tmp_path / references["development_responses"]["path"]
    raw_read = Path.read_bytes
    response_released = False
    accesses = []

    def guarded(path):
        accesses.append(path)
        if path == response_path and not response_released:
            raise AssertionError("development response read before prediction freeze and new release")
        return raw_read(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    assert development.check_custody_gate(tmp_path, protocol_path, initial)["ready"] is True
    assert not any(path == tmp_path / references["training_projection"]["path"] for path in accesses)
    summary = development.execute_training(tmp_path, protocol_path, initial, "outputs/native_population_development/run")
    assert summary["development_prediction_rows"] == 288
    assert summary["selected_development_model"] is None
    assert response_path not in accesses
    with pytest.raises(FileExistsError):
        development.execute_training(tmp_path, protocol_path, initial, "outputs/native_population_development/run")
    release = allowlist("development_response_export", summary["prediction_freeze"]["sha256"])
    response_released = True

    def no_refit(*args, **kwargs):
        raise AssertionError("development outcomes must never trigger refitting")

    monkeypatch.setattr(development, "fit_training", no_refit)
    result = development.execute_scoring(tmp_path, protocol_path, release, summary["prediction_freeze"])
    assert result["status"] == "complete_development_diagnostics"
    assert result["selected_development_model"] == "first_order_response"
    assert set(result["models"]) == set(development.MODEL_IDS)
    assert result["strong_claim_authorized"] is False
    assert result["scientific_scope"] == development.SCIENTIFIC_SCOPE
    assert response_path in accesses
    with pytest.raises(FileExistsError):
        development.execute_scoring(tmp_path, protocol_path, release, summary["prediction_freeze"])
