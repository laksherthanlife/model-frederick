from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis import biology_learning_audit as audit
from ystwin.analysis.frozen_runtime import audit_historical_evidence, materialize_frozen_runtime
from ystwin.analysis.hog_data import load_native_hog_data
from ystwin.analysis.native_physiology import (
    evaluate_native_exchange_model,
    fit_native_exchange_model,
    load_native_chemostat_data,
)

from _frozen_runtime_helpers import run_frozen_python


ROOT = Path(__file__).resolve().parents[1]


def _require(ledger, claim, *, root):
    """Invoke the original gate, transporting only JSON and its real ValueError."""
    process = run_frozen_python(root, """
        from ystwin.analysis import biology_learning_audit as audit

        request = json.load(sys.stdin)
        try:
            result = audit.require_claim(request["ledger"], request["claim"], root=root,
                                         artifact_source="portable")
        except ValueError as exc:
            print(json.dumps({"error": str(exc)}))
        else:
            print(json.dumps({"result": result}, allow_nan=False))
    """, payload={"ledger": ledger, "claim": claim})
    assert process.returncode == 0, process.stderr
    response = json.loads(process.stdout)
    if "error" in response:
        raise ValueError(response["error"])
    return response["result"]


@pytest.fixture(scope="module")
def frozen_root():
    with materialize_frozen_runtime(ROOT, operation="audit") as snapshot:
        yield snapshot.root


@pytest.fixture(scope="module")
def native_source(frozen_root):
    return load_native_hog_data(frozen_root / "data/hog2013")


@pytest.fixture(scope="module")
def ledger(frozen_root):
    report = audit_historical_evidence(frozen_root, source_repository=ROOT)
    assert report["historical_execution"]["isolated_interpreter"]
    return report["ledger"]


def test_verified_native_ledger_distinguishes_estimation_from_discovery(ledger, frozen_root):
    assert ledger["checkpoint_sha256"] == audit.CHECKPOINT_SHA256
    assert all(item["matches"] for item in ledger["integrity_before"])
    assert ledger["artifact_source"] == "portable"
    assert "product outcomes" not in ledger["excluded"]
    assert "previously recorded outcomes and scores" in ledger["source_access_scope"]
    assert ledger["evidence_predicates"]["checkpoint_integrity"]["value"] is None
    assert ledger["integrity_scope"]["original_checkpoint_bytes"]["status"] == "not_checked"
    assert _require(ledger, "portable_scientific_content", root=frozen_root)["supported"]
    for claim in ("conditional_parameter_estimation", "empirical_interpolation", "supplied_kinetic_law", "integrity_established"):
        with pytest.raises(ValueError, match="checkpoint_integrity"):
            _require(ledger, claim, root=frozen_root)
    components = {item["id"]: item for item in ledger["components"]}
    hog = components["native_hog"]
    assert hog["native_rows"] == 82
    assert hog["unshocked_rows_used_for_fit"] == 28
    assert hog["kinetic_degrees_of_freedom"] == 3
    assert len(hog["absolute_parameter_overrides"]) == 4
    assert len(hog["estimated_observation_gains"]) == 2
    assert hog["remaining_fixed_parameter_count_including_protocol_inputs"] == hog["source_global_parameter_count"] - 4
    assert components["source_hog_structure"]["functional_forms_learned"] is False
    assert components["source_observation_and_growth_assumptions"]["OD_is_time_only_assignment"]
    assert components["native_gem_bridge"]["classification"] == "untested_coupling"
    assert components["native_gem_bridge"]["structure_verified"]
    assert components["independent_mechanistic_evidence"]["available"] is False
    json.dumps(ledger, allow_nan=False)


@pytest.mark.parametrize("claim", [
    "independent_native_validation", "independent_mechanistic_evidence",
    "teacher_recovery_is_biological_law_learning", "absolute_phosphorylation_fraction_learned",
    "individual_fast_rates_identified", "native_gem_coupling_validated",
    "genotype_copy_to_active_enzyme_validated", "realized_allocation_learned",
])
def test_native_artifacts_cannot_earn_unsupported_claims(ledger, claim, frozen_root):
    with pytest.raises(ValueError, match="unsupported claim"):
        _require(ledger, claim, root=frozen_root)


def test_hashes_and_guessed_kind_are_not_learning_evidence(ledger):
    evidence = {"source_kind": "measured", "independent": True, "training_error": 0.0,
                "checkpoint_integrity": ledger["evidence_predicates"]["checkpoint_integrity"]}
    for claim in ("conditional_parameter_estimation", "independent_mechanistic_evidence"):
        result = audit.assess_claim(claim, evidence)
        assert not result["supported"]
        assert any(item["state"] == "unresolved" for item in result["predicates"])
    advisory = audit.assess_claim("integrity_established", evidence)
    assert not advisory["metadata_consistent"]
    assert not advisory["supported"] and not advisory["authorized"]


def test_synthetic_supervision_is_rejected_even_with_perfect_recovery(ledger):
    evidence = deepcopy(ledger["evidence_predicates"])
    evidence.update(training_rmse=0.0, heldout_rmse=0.0, source_kind="independent_measurement")
    result = audit.assess_claim("teacher_recovery_is_biological_law_learning", evidence)
    assert result["status"] == "rejected"
    assert result["predicates"][0]["evidence"]["verification"] == "source_structure"


def test_missing_or_unbound_predicates_fail_closed(ledger):
    for mutation in ({"value": True}, True, {"value": "true", "verification": "source_kind"}):
        result = audit.assess_claim("integrity_established", {"checkpoint_integrity": mutation})
        assert not result["supported"]
    evidence = deepcopy(ledger["evidence_predicates"])
    evidence["checkpoint_integrity"]["evidence"][0]["sha256"] = "z" * 64
    assert not audit.assess_claim("integrity_established", evidence)["supported"]
    with pytest.raises(ValueError, match="unknown learning claim"):
        audit.assess_claim("biological_because_named_measured.csv", evidence)


def test_historical_exposure_and_unknown_split_timing_remain_visible(ledger):
    evidence = ledger["evidence_predicates"]
    assert evidence["historical_lineage_disjoint"]["value"] is False
    assert evidence["whole_intervention_groups_reserved"]["value"] is False
    assert evidence["split_before_all_selection"]["value"] is None
    assert evidence["independent_experimental_units"]["value"] is None
    inventory = ledger["native_observation_inventory"]
    assert inventory["imported_rows"] == 317
    assert inventory["unselected_rows"] == 235
    assert inventory["original_fit_status_counts"] == {"published_fitting_dataset": 317}
    host = next(item for item in ledger["components"] if item["id"] == "host_prior")
    assert host["historical_fit_lineage"] == "incomplete"
    assert host["upstream_fit_dataset_overlap"] == "unresolved, not assumed absent"


def test_native_binding_uses_cells_values_units_and_provenance(native_source):
    original = native_source.observations
    rows = original.loc[original.genotype.eq("wild_type")].copy()
    rows["split"] = "fit"
    assert audit.bind_native_observations(rows, original)["verified"]
    assert audit.bind_native_observations(rows.iloc[::-1], original)["verified"]
    changed = rows.copy()
    changed.iloc[0, changed.columns.get_loc("value")] *= 2
    result = audit.bind_native_observations(changed, original)
    assert not result["verified"]
    assert any(item.get("field") == "value" for item in result["issues"])
    assert not audit.bind_native_observations(pd.concat([rows, rows.iloc[:1]]), original)["verified"]


@pytest.mark.parametrize("field,replacement", [
    ("time_model_s", 12345.0), ("unit", "mmol/L"),
    ("measurement_provenance", "measured"), ("original_fit_status", "independent"),
    ("nacl_molar", 0.8), ("source_cell", "Z999"),
])
def test_stale_provenance_cannot_launder_changed_native_rows(native_source, field, replacement):
    original = native_source.observations
    rows = original.iloc[:13].copy()
    rows.iloc[0, rows.columns.get_loc(field)] = replacement
    result = audit.bind_native_observations(rows, original)
    assert not result["verified"]
    assert result["issues"]


def test_software_western_scale_identity_is_not_absolute_calibration(native_source):
    rows = native_source.observations.loc[native_source.observations.genotype.eq("wild_type")].copy()
    raw = np.linspace(0.01, 0.5, len(rows))
    initial, _, residual = audit._loss(rows, raw)
    mask = rows.observable_id.isin({"Hog1PP_measured", "Gpd1_measured"})
    rows.loc[mask, "value"] *= 2.0
    scaled, _, second_residual = audit._loss(rows, raw)
    np.testing.assert_allclose(second_residual, residual, rtol=0, atol=1e-12)
    for key in initial["gains"]:
        assert scaled["gains"][key] == pytest.approx(2 * initial["gains"][key])
    assert not audit.bind_native_observations(rows, native_source.observations)["verified"]


def test_software_time_baseline_cannot_fit_on_future_values(native_source):
    rows = native_source.observations.loc[native_source.observations.genotype.eq("wild_type")].copy()
    prefix = rows.loc[rows.time_model_s <= 7200].reset_index(drop=True)
    suffix = rows.loc[rows.time_model_s > 7200].reset_index(drop=True)
    first = audit._baselines(prefix, suffix)
    changed = suffix.copy()
    changed["value"] *= 10.0
    second = audit._baselines(prefix, changed)
    for name in first:
        assert first[name]["prediction"] == second[name]["prediction"]
        assert first[name]["coefficients"] == second[name]["coefficients"]
        assert first[name]["metrics"]["training_scales"] == second[name]["metrics"]["training_scales"]


def test_interpolator_recovers_shuffled_knots_without_mechanistic_credit(ledger, frozen_root):
    rows = load_native_chemostat_data(frozen_root / "data/physiology/chemostatData_VanHoek1998.tsv").split()["train"]
    changed = rows.copy()
    selected = rows.observable_id.eq("GlucoseUptake")
    changed.loc[selected, "value"] = rows.loc[selected, "value"].to_numpy()[::-1]
    changed.loc[selected, "reported_value"] = rows.loc[selected, "reported_value"].to_numpy()[::-1]
    model = fit_native_exchange_model(changed)
    result = evaluate_native_exchange_model(model, changed)
    assert result["metrics"]["GlucoseUptake"]["rmse"] == 0.0
    assert model.to_dict()["method"] == "piecewise_linear"
    assert result["n_scored"] < result["n_readouts"]
    assert not audit.assess_claim("independent_mechanistic_evidence", ledger["evidence_predicates"])["supported"]


def test_audit_refuses_existing_or_nonowned_destinations_without_writes():
    for path in (ROOT / "outputs", ROOT / "data" / "biology_learning_audit_test",
                 ROOT / "outputs" / "native_training_run_02", ROOT / "outputs" / "nested" / "biology_learning_audit_test"):
        with pytest.raises(ValueError, match="fresh direct"):
            audit.run_biology_learning_audit(ROOT, path)


def test_audit_plan_has_fixed_controls_and_no_goal_seeking_threshold():
    assert audit._PLAN["cutoff_model_time_s"] == 7200.0
    assert audit._PLAN["nulls"] == ["within_trace_shuffle", "within_trace_reverse_time"]
    assert audit._PLAN["common_fast_rate_factors"] == [0.2, 1.0, 5.0]
    assert audit._PLAN["performance_thresholds"] is None
    assert audit._PLAN["selection_from_diagnostic_results"] is False
    assert audit._PLAN["significance_or_confidence_intervals"] is None


def _forged_strong_evidence():
    return {name: {
        "value": True, "verification": "source_bound_values",
        "evidence": [{"path": "nonexistent/intervention.json", "sha256": "0" * 64,
                      "location": {"json_pointer": "/independent"}}],
    } for name in audit._CLAIMS["independent_mechanistic_evidence"]}


def test_forged_reference_shapes_cannot_authorize_mechanistic_evidence():
    evidence = _forged_strong_evidence()
    advisory = audit.assess_claim("independent_mechanistic_evidence", evidence)
    assert not advisory["supported"]
    with pytest.raises(ValueError, match="unsupported claim"):
        audit.require_claim({"evidence_predicates": evidence}, "independent_mechanistic_evidence",
                            root=ROOT, artifact_source="portable")


def test_tampered_real_refs_and_recomputed_internal_digest_do_not_authorize(ledger, frozen_root):
    import hashlib

    forged = deepcopy(ledger)
    for name in audit._CLAIMS["independent_mechanistic_evidence"]:
        forged["evidence_predicates"][name]["value"] = True
    forged["sha256"] = hashlib.sha256(json.dumps(forged, sort_keys=True).encode()).hexdigest()
    with pytest.raises(ValueError, match="unsupported claim"):
        _require(forged, "independent_mechanistic_evidence", root=frozen_root)


def test_advisory_metadata_never_authorizes_even_supported_native_grade(ledger):
    result = audit.assess_claim("portable_scientific_content", ledger["evidence_predicates"])
    assert not result["supported"]
    assert result["assessment_kind"] == "advisory_metadata"
    assert result["authorized"] is False
    assert result["metadata_consistent"] is True


def test_release_gate_rejects_tampered_reference_for_supported_grade(ledger, frozen_root):
    forged = deepcopy(ledger)
    forged["evidence_predicates"]["native_cells_bound"]["evidence"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="authoritative collector"):
        _require(forged, "portable_scientific_content", root=frozen_root)


def test_release_gate_rejects_tampered_parameter_payload(ledger, frozen_root):
    forged = deepcopy(ledger)
    hog = next(item for item in forged["components"] if item["id"] == "native_hog")
    hog["absolute_parameter_overrides"]["kv16f_1"] *= 2.0
    with pytest.raises(ValueError, match="authoritative collector"):
        _require(forged, "portable_scientific_content", root=frozen_root)


def test_release_gate_rechecks_current_frozen_bytes_without_modifying_files(ledger, frozen_root):
    process = run_frozen_python(frozen_root, r"""
        from unittest.mock import patch
        import pytest
        from ystwin.analysis import biology_learning_audit as audit

        ledger = json.load(sys.stdin)
        original = audit._bytes

        def changed_in_memory(root, path):
            payload = original(root, path)
            return payload + b"\n" if path == "src/ystwin/analysis/hog_learning.py" else payload

        with patch.object(audit, "_bytes", changed_in_memory):
            with pytest.raises(ValueError, match="frozen dependencies changed"):
                audit.require_claim(ledger, "portable_scientific_content", root=root,
                                    artifact_source="portable")
    """, payload=ledger)
    assert process.returncode == 0, process.stderr


def test_json_roundtrip_can_use_root_verified_gate_but_snapshot_is_not_authorization(ledger, frozen_root):
    report = json.loads(json.dumps(ledger))
    result = _require(report, "portable_scientific_content", root=frozen_root)
    assert result["assessment_kind"] == "verified_release_gate"
    assert result["authorized"] is True
    assert all(item["authorized"] is False for item in report["claims"])


def test_release_gate_never_opens_caller_supplied_reference_paths(ledger, frozen_root):
    process = run_frozen_python(frozen_root, """
        from copy import deepcopy
        from unittest.mock import patch
        import pytest
        from ystwin.analysis import biology_learning_audit as audit

        ledger = json.load(sys.stdin)
        original = audit._bytes
        reads = []

        def record(root, path):
            reads.append(path)
            return original(root, path)

        forged = deepcopy(ledger)
        forbidden = "nonexistent/never_open_outcomes.json"
        forged["evidence_predicates"]["native_cells_bound"]["evidence"][0]["path"] = forbidden
        with patch.object(audit, "_bytes", record):
            with pytest.raises(ValueError, match="authoritative collector"):
                audit.require_claim(forged, "portable_scientific_content", root=root,
                                    artifact_source="portable")
        assert forbidden not in reads
        assert audit.CHECKPOINT_PATH not in reads
        assert ledger["integrity_scope"]["checkpoint_reference"]["path"] in reads
    """, payload=ledger)
    assert process.returncode == 0, process.stderr


def test_portable_content_cannot_launder_the_original_checkpoint_byte_flag(ledger, frozen_root):
    forged = deepcopy(ledger)
    forged["evidence_predicates"]["checkpoint_integrity"]["value"] = True
    forged["integrity_scope"]["original_checkpoint_bytes"].update(
        status="verified", matches=True, actual_sha256=audit.CHECKPOINT_SHA256)
    with pytest.raises(ValueError, match="authoritative collector"):
        _require(forged, "integrity_established", root=frozen_root)


def test_release_gate_uses_exact_json_types_not_truthy_equivalence(ledger, frozen_root):
    forged = deepcopy(ledger)
    forged["evidence_predicates"]["native_cells_bound"]["value"] = 1
    with pytest.raises(ValueError, match="authoritative collector"):
        _require(forged, "portable_scientific_content", root=frozen_root)


@pytest.mark.parametrize("claim", [
    "independent_native_validation", "independent_mechanistic_evidence",
    "teacher_recovery_is_biological_law_learning", "absolute_phosphorylation_fraction_learned",
    "individual_fast_rates_identified", "native_gem_coupling_validated",
    "genotype_copy_to_active_enzyme_validated", "realized_allocation_learned",
])
def test_unimplemented_stronger_assessors_fail_closed_even_if_all_flags_true(ledger, claim, frozen_root):
    forged = deepcopy(ledger)
    for fact in forged["evidence_predicates"].values():
        fact["value"] = True
    forged["claims"] = [{"claim": claim, "supported": True, "authorized": True, "verified": True}]
    assert not audit.assess_claim(claim, forged["evidence_predicates"])["supported"]
    with pytest.raises(ValueError, match="unsupported claim.*implemented"):
        _require(forged, claim, root=frozen_root)


def test_parent_gauge_and_catalogue_remain_qualifications_not_validation(ledger, frozen_root):
    qualifications = ledger["additional_qualifications"]
    gauge = qualifications["gpd1_source_gauge"]
    assert gauge["source_model_identity_matches"]
    assert gauge["within_restricted_three_group_fit"] is False
    assert gauge["numerical_probe_rerun_by_collector"] is False
    assert gauge["reported_probe"]["maximum_absolute_rhs_covariance_error"] == 0.0
    assert gauge["biological_validation_credit"] is False
    assert qualifications["pool_counterflow"]["kinetic_parameter_equivalence_established"] is False
    candidate = qualifications["granados2018_candidate"]
    assert candidate["catalogue_status"] == "metadata_verified_measurement_arrays_not_opened"
    assert candidate["measurement_arrays_opened_by_audit"] is False
    assert candidate["validation_completed"] is False
    assert candidate["biological_validation_credit"] is False
    assert candidate["leakage_and_interpretation_gates"]
    with pytest.raises(ValueError, match="unsupported claim"):
        _require(ledger, "independent_mechanistic_evidence", root=frozen_root)


def test_reused_diagnostics_have_byte_binding_and_no_new_evidence_credit(frozen_root):
    import hashlib

    payload, receipt = audit._reused_native_diagnostics(frozen_root)
    assert receipt["source"]["sha256"] == hashlib.sha256(payload).hexdigest()
    assert receipt["mode"] == "reused_unchanged_retrospective_evidence"
    assert receipt["numerical_diagnostics_rerun"] is False
    assert receipt["biological_validation_credit"] is False
    assert receipt["checkpoint_sha256"] == audit.CHECKPOINT_SHA256
    assert set(json.loads(payload)) == {"hog", "native_exchange"}


def test_reused_diagnostics_reject_a_mismatched_prior_checkpoint(frozen_root, monkeypatch):
    original = audit._bytes

    def changed_plan(root, path):
        payload = original(root, path)
        if path == "outputs/biology_learning_audit_01/diagnostic_plan.json":
            plan = json.loads(payload)
            plan["checkpoint_sha256"] = "0" * 64
            return json.dumps(plan).encode()
        return payload

    monkeypatch.setattr(audit, "_bytes", changed_plan)
    with pytest.raises(ValueError, match="prior diagnostic plan"):
        audit._reused_native_diagnostics(frozen_root)
