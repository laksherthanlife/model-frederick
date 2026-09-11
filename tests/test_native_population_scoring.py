from __future__ import annotations

from copy import deepcopy
from dataclasses import FrozenInstanceError, asdict, replace
from pathlib import Path

import pytest

from ystwin.analysis import native_population_development as frozen
from ystwin.analysis import native_population_scoring as scoring

from test_population_historical_sources import historical_population_root as historical_population_root


MANIFEST = Path("data/native_law_v2/granados/development/phase_c_01/release_manifest.json")
MANIFEST_SHA = "0e2d18b1d239f33fd21165219c512b52546635b733ec47d67c85f726a919ec82"


@pytest.fixture(scope="module")
def admitted_metadata(historical_population_root):
    return scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)


def test_valid_separate_training_executor_and_response_approvals_are_accepted_without_values(monkeypatch, historical_population_root):
    original = Path.read_bytes
    forbidden = historical_population_root / "data/native_law_v2/granados/development_exports"

    def no_measurements(path):
        assert not path.resolve().is_relative_to(forbidden.resolve())
        assert path.name != "coverage_audit.json"
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", no_measurements)
    bundle = scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)
    assert bundle.training.data_approval["sha256"] == "87536dc42cccde2733daf4cb52f54024661a7d58386acd180623adf266126533"
    assert bundle.training.executor_approval["sha256"] == "04ac928add5528beabd2a9f77a2055643dc120d3b76bff7315ffdc3f062f786e"
    assert bundle.release.manifest["custody_decision"]["sha256"] == "ccdf78fa1d986758a9f760747d755b02fd2d52c11098e14c3043ee7d8a47806e"
    assert len(bundle.prediction.data["predictions"]) == 288
    assert scoring.validate_lineage(bundle) is True


@pytest.mark.parametrize("defect", ["wrong_phase", "reused_training_approval", "reused_executor_approval", "wrong_request", "wrong_frozen_data", "wrong_code_lineage", "wrong_response_lineage", "changed_cohort", "changed_prediction", "missing_family"])
def test_typed_lineage_requires_actual_matching_links_not_reference_shapes(admitted_metadata, defect):
    bundle = deepcopy(admitted_metadata)
    manifest = bundle.release.manifest
    if defect == "wrong_phase":
        manifest["phase"] = "training_and_conditioning_export"
    elif defect == "reused_training_approval":
        manifest["custody_decision"] = bundle.training.data_approval
    elif defect == "reused_executor_approval":
        manifest["custody_decision"] = bundle.training.executor_approval
    elif defect == "wrong_request":
        bundle.release.decision["request_sha256"] = bundle.training.data_approval["sha256"]
    elif defect == "wrong_frozen_data":
        bundle.prediction.data["development_inputs_sha256"] = "0" * 64
    elif defect == "wrong_code_lineage":
        manifest["frozen_executor_artifacts"][0]["sha256"] = "0" * 64
    elif defect == "wrong_response_lineage":
        bundle.release.source_binding["response_release_approval_sha256"] = bundle.training.data_approval["sha256"]
    elif defect == "changed_cohort":
        key = next(iter(bundle.release.source_replay["frozen_prefix_cohorts_used"]))
        bundle.release.source_replay["frozen_prefix_cohorts_used"][key] = [0, 1]
    elif defect == "changed_prediction":
        bundle.prediction.data["predictions"][0]["mean"] += 0.5
    else:
        bundle.prediction.data["predictions"] = [row for row in bundle.prediction.data["predictions"] if row["model_id"] != "transient_response"]
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError)):
        scoring.validate_lineage(bundle)


def test_tampered_manifest_bytes_are_rejected_before_any_referenced_file(monkeypatch, historical_population_root):
    original = Path.read_bytes
    opened = []

    def tampered(path):
        opened.append(path)
        assert path.resolve() == (historical_population_root / MANIFEST).resolve()
        return original(path).replace(b'"development_response_export"', b'"training_and_conditioning_export"')

    monkeypatch.setattr(Path, "read_bytes", tampered)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="digest|byte"):
        scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)
    assert len(opened) == 1


def test_tampered_response_bytes_cannot_be_scored(admitted_metadata, monkeypatch, historical_population_root):
    target = historical_population_root / admitted_metadata.release.manifest["response_projection"]["path"]
    original = Path.read_bytes

    def changed(path):
        if path.resolve() == target.resolve():
            return b"{}"
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", changed)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="digest|byte"):
        scoring.load_released_observations(historical_population_root, admitted_metadata)


def test_tampered_bound_code_fails_before_response_access(admitted_metadata, monkeypatch, historical_population_root):
    target = historical_population_root / admitted_metadata.release.manifest["frozen_executor_artifacts"][0]["path"]
    original = Path.read_bytes
    projection_root = historical_population_root / "data/native_law_v2/granados/development_exports"

    def changed(path):
        assert not path.resolve().is_relative_to(projection_root.resolve())
        if path.resolve() == target.resolve():
            return b"pass\n"
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", changed)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="digest|byte"):
        scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)


def synthetic_means(bundle):
    prefixes = {row["group_id"]: row for row in bundle.prediction.data["prefix_summaries"]}
    return [{**coordinate, "observed_mean": prefixes[coordinate["group_id"]]["mu"], "prefix_cells": len(prefixes[coordinate["group_id"]]["eligible_cell_indices"]), "valid_cells": len(prefixes[coordinate["group_id"]]["eligible_cell_indices"]), "source_missing_cells": 0, "invalid_cells": 0, "cell_sample_sd": 0.0} for coordinate in bundle.prediction.conditioning["development_coordinates"]]


def test_scoring_uses_frozen_metric_and_prediction_helpers_without_refit(admitted_metadata, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("scoring cannot call any fitter or legacy same-approval executor")

    monkeypatch.setattr(frozen, "fit_training", forbidden)
    monkeypatch.setattr(frozen, "execute_scoring", forbidden)
    report, summary = scoring.score_frozen_means(admitted_metadata, synthetic_means(admitted_metadata))
    assert report["status"] == "complete_development_diagnostics"
    assert report["selected_development_model"] == "prefix_persistence"
    assert len(report["group_evaluations"]) == 12
    assert sum(len(row["frames"]) for row in report["group_evaluations"]) == 288
    assert set(summary["models"]) == set(frozen.MODEL_IDS)
    assert summary["strong_claim_authorized"] is False
    assert summary["scientific_scope"] == "retrospective_source_cohort_development"
    assert summary["parameter_bound_warnings"]["transient_response"]["A"] == "lower"


def test_missing_mean_cannot_be_removed_to_create_a_winner(admitted_metadata):
    means = synthetic_means(admitted_metadata)
    means[0].update(observed_mean=None, valid_cells=0, source_missing_cells=means[0]["prefix_cells"], cell_sample_sd=None)
    report, summary = scoring.score_frozen_means(admitted_metadata, means)
    assert report["status"] == "incomplete_observation_coverage"
    assert summary["selected_development_model"] is None
    assert len(report["group_evaluations"]) == 12
    assert all(len(group["frames"]) == 24 for group in report["group_evaluations"])


def test_no_new_output_can_overwrite_the_frozen_run_or_an_existing_score(tmp_path, admitted_metadata, historical_population_root):
    with pytest.raises(scoring.ScoringError):
        scoring.new_score_store(historical_population_root, admitted_metadata, "outputs/native_population_development/native_training_01")
    directory = "outputs/native_population_development/phase_c_scoring_engineering"
    store = scoring.new_score_store(tmp_path, admitted_metadata, directory)
    store.write("scoring_access_intent.json", {"synthetic_engineering": True}, "operational_event_log")
    with pytest.raises((scoring.ScoringError, FileExistsError)):
        scoring.new_score_store(tmp_path, admitted_metadata, directory)


def test_authority_pins_are_immutable_and_detached_from_cached_references(admitted_metadata):
    changed = deepcopy(admitted_metadata)
    pins = changed.admission_pins
    manifest_sha = pins.manifest.sha256
    with pytest.raises(FrozenInstanceError):
        pins.manifest.sha256 = "0" * 64
    with pytest.raises(FrozenInstanceError):
        pins.artifacts = ()
    changed.release.manifest_reference["sha256"] = "0" * 64
    changed.release.manifest["response_projection"]["sha256"] = "1" * 64
    assert pins.manifest.sha256 == manifest_sha == MANIFEST_SHA
    assert scoring.reload_admission(pins).admission_pins == admitted_metadata.admission_pins


def test_cached_projection_fields_are_ignored_in_favour_of_pinned_authority(admitted_metadata, monkeypatch, historical_population_root):
    expected = scoring.load_released_observations(historical_population_root, admitted_metadata)
    altered = deepcopy(admitted_metadata)
    forbidden = historical_population_root / "data/native_law_v2/granados/development/private/source_coverage.json"
    rebound = {**altered.release.manifest["response_projection"], "path": str(forbidden.relative_to(historical_population_root)), "sha256": "0" * 64}
    altered.release.manifest["response_projection"] = rebound
    altered.release.manifest["response_data_allowlist"] = [rebound]
    altered.release.source_replay["response_projection"] = rebound
    altered.release.source_binding["ordered_payload_sha256"] = "1" * 64
    original = Path.read_bytes

    def guarded(path):
        assert path.resolve() != forbidden.resolve()
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    assert scoring.load_released_observations(historical_population_root, altered) == expected


def test_cached_response_approval_cannot_replace_pinned_phase_c_decision(admitted_metadata):
    altered = deepcopy(admitted_metadata)
    altered.release.manifest["custody_decision"] = deepcopy(altered.training.data_approval)
    fresh = scoring.revalidate_admission(altered)
    assert fresh.release.manifest["custody_decision"] == admitted_metadata.release.manifest["custody_decision"]
    assert fresh.release.manifest["custody_decision"] != fresh.training.data_approval
    assert scoring.validate_lineage(fresh)


def test_changed_pinned_event_bytes_are_rejected_even_when_typed_links_still_match(admitted_metadata, monkeypatch, historical_population_root):
    event_path = historical_population_root / admitted_metadata.release.manifest["operational_event_reference"].split("#")[0]
    original = Path.read_bytes
    changed = []

    def substituted(path):
        raw = original(path)
        if path.resolve() == event_path.resolve():
            modified = raw.replace(b"2026-09-07T", b"2025-09-07T", 1)
            assert modified != raw and len(modified) == len(raw)
            changed.append(path)
            return modified
        return raw

    monkeypatch.setattr(Path, "read_bytes", substituted)
    with pytest.raises(scoring.ScoringError, match="pinned artifact.*digest"):
        scoring.revalidate_admission(admitted_metadata)
    assert changed


def test_canonical_json_cache_roundtrip_reloads_the_same_admission(admitted_metadata, historical_population_root):
    payload = frozen.strict_json(frozen.canonical_bytes(asdict(admitted_metadata.release)))
    restored = replace(admitted_metadata, release=scoring.PhaseCResponseRelease(**payload))
    fresh = scoring.revalidate_admission(restored)
    assert fresh == admitted_metadata
    assert scoring.load_released_observations(historical_population_root, fresh) == scoring.load_released_observations(historical_population_root, admitted_metadata)


def test_bit_identical_authorized_projection_replica_is_not_rejected_for_equal_values(tmp_path, admitted_metadata, historical_population_root):
    for key in ("response_projection", "previously_released_development_inputs", "coverage_audit"):
        reference = admitted_metadata.release.manifest[key]
        content = frozen.read_artifact(historical_population_root, reference, admitted_metadata.protocol, reference["role"], parse=False)
        target = tmp_path / reference["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    replica = scoring.load_released_observations(tmp_path, admitted_metadata)
    original = scoring.load_released_observations(historical_population_root, admitted_metadata)
    assert replica == original


def test_published_result_verification_does_not_fit_predict_or_select(tmp_path, admitted_metadata, monkeypatch, historical_population_root):
    def forbidden(*args, **kwargs):
        raise AssertionError("verification cannot run a fresh fit, prediction batch or selection")

    monkeypatch.setattr(frozen, "fit_training", forbidden)
    monkeypatch.setattr(frozen, "fit_family", forbidden)
    monkeypatch.setattr(frozen, "predict_development", forbidden)
    monkeypatch.setattr(frozen, "select_model", forbidden)
    monkeypatch.setattr(scoring, "score_once", forbidden)
    monkeypatch.setattr(scoring, "new_score_store", lambda root, bundle, directory: frozen.ArtifactStore(tmp_path, directory))
    result = scoring.verify_published_results(historical_population_root, "outputs/native_population_development/verification_fixture")
    assert result["published_metrics_exactly_unchanged"] is True
    assert result["new_fit_performed"] is False
    assert result["new_predictions_generated"] is False
    assert result["new_selection_performed"] is False
    record = frozen.read_artifact(tmp_path, result["verification_artifact"], admitted_metadata.protocol, "operational_event_log")
    assert record["issue"] == "SCORING001"
    assert record["original_adapter_was_bug_free"] is False
    assert record["historical_reports_overwritten"] is False
    assert len(record["original_adapter_snapshots"]) == 3
    assert len(record["fixed_adapter_snapshots"]) == 3
