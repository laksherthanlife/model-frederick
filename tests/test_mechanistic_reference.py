from __future__ import annotations

import ast
import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ystwin.analysis import mechanistic_reference as reference_module
from ystwin.analysis.mechanistic_reference import (
    AssayScale,
    FrozenReference,
    ReferenceDesign,
    ReferenceExperiment,
    SaltChange,
    default_reference_design,
    evaluate_native_reference,
    freeze_reference,
    generate_reference,
    run_synthetic_benchmark,
    score_parameter_recovery,
    score_synthetic_reference,
    sensitivity_identifiability,
    training_mean_baseline,
    verify_reference_freeze,
)
from ystwin.analysis.parameter_evidence import EvidenceGap, load_parameter_evidence, load_reference_model

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def full_reference():
    return generate_reference(freeze_reference(default_reference_design(seed=17), product_outcomes_seen=False))


@pytest.fixture
def small_design():
    times = (0, 20, 25, 40, 60)
    return ReferenceDesign(
        (
            ReferenceExperiment("train", "hog2013_wt", times, (SaltChange(20, 0.4),), "train"),
            ReferenceExperiment("test", "hog2013_wt", times, (SaltChange(20, 0.8),), "test", holdout_axis="environment"),
        ),
        (AssayScale("Hog1PP_measured", gain=1, background=0, noise_sd=0),),
        seed=31,
    )


def test_reference_contains_actual_source_models_and_all_four_challenge_axes(full_reference):
    truth = full_reference.evaluator_truth()
    assert len(truth["experiments"]) == 10
    assert {row["holdout_axis"] for row in truth["targets"]} == {
        "structure", "environment", "genotype", "exposure_history",
    }
    assert len({e["provenance"]["model_id"] for e in truth["experiments"]}) == 6
    for experiment in truth["experiments"]:
        assert len(experiment["species_ids"]) == 29
        assert len(experiment["reaction_ids"]) in (54, 58)
        assert experiment["provenance"]["imported_learner_coefficients"] is False
        assert experiment["provenance"]["identifiability"]["unresolved_directions"]
        assert np.isfinite(experiment["states"]).all()
    assert truth["full_biology_claim"] is False
    assert truth["independent_real_validation"] is False


def test_published_structure_replays_original_sbml_without_learned_values(full_reference):
    truth = full_reference.evaluator_truth()
    experiment = next(e for e in truth["experiments"] if e["design"]["experiment_id"] == "train_step")
    model = load_reference_model("hog2013_wt")
    original = model.simulate(experiment["times_s"], breakpoints=[3600, 3605, 4800, model.parameter_values["kv22_Hog1D_t"]])
    np.testing.assert_allclose(experiment["states"], original.states, rtol=2e-7, atol=1e-10)
    assert experiment["provenance"]["source_parameter_values"] == model.parameter_values
    assert experiment["provenance"]["reference_parameter_values"] == model.parameter_values
    assert experiment["provenance"]["structural_interventions"] == []


def test_structural_intervention_changes_only_the_declared_edge(full_reference):
    experiments = {e["design"]["experiment_id"]: e for e in full_reference.evaluator_truth()["experiments"]}
    original, altered = experiments["train_step"], experiments["test_structure"]
    original_parameters = original["provenance"]["reference_parameter_values"]
    altered_parameters = altered["provenance"]["reference_parameter_values"]
    assert {key for key in original_parameters if original_parameters[key] != altered_parameters[key]} == {"kv17f_1"}
    assert altered_parameters["kv17f_1"] == 0
    assert altered["provenance"]["structural_interventions"][0]["status"] == "declared_counterfactual_not_measured_mutant"
    index = altered["species_ids"].index("Gpd1")
    assert not np.allclose(np.asarray(altered["states"])[:, index], np.asarray(original["states"])[:, index], rtol=1e-5)


def test_exposure_history_preserves_memory_instead_of_resetting_state(full_reference):
    truth = full_reference.evaluator_truth()
    experiments = {e["design"]["experiment_id"]: e for e in truth["experiments"]}
    step, history = experiments["train_step"], experiments["test_history"]
    np.testing.assert_allclose(step["states"][:8], history["states"][:8], rtol=2e-6, atol=1e-9)
    assert not np.allclose(step["states"][-1], history["states"][-1], rtol=1e-5, atol=1e-8)
    assert "never reset state" in history["provenance"]["clock_policy"]
    query = next(row for row in full_reference.learner_inputs()["prediction_queries"] if row["experiment_id"] == "test_history" and row["time_s"] == 7200)
    assert query["nacl_molar"] == pytest.approx(0.4)


def test_learner_view_excludes_truth_coefficients_assay_gains_and_test_values(full_reference):
    view = full_reference.learner_inputs()
    assert set(view) == {"schema_version", "reference_id", "domain", "training_observations", "prediction_queries", "input_scope"}
    forbidden = {"states", "reaction_rates", "provenance", "source_parameter_values", "reference_parameter_values", "assay_scales", "gain", "structure", "holdout_axis", "expectation"}
    assert all(not (set(row) & forbidden) for row in view["prediction_queries"])
    assert all("value" not in row for row in view["prediction_queries"])
    assert all("known_od_input" in row for row in view["prediction_queries"])
    assert {r["experiment_id"] for r in view["training_observations"]}.isdisjoint({r["experiment_id"] for r in view["prediction_queries"]})
    view["training_observations"][0]["value"] = -999
    truth = full_reference.evaluator_truth()
    truth["experiments"][0]["states"][0][0] = -999
    assert full_reference.learner_inputs()["training_observations"][0]["value"] != -999
    assert full_reference.evaluator_truth()["experiments"][0]["states"][0][0] != -999


def test_synthetic_oracle_recovery_is_not_claimed_as_real_validation(full_reference):
    targets = full_reference.evaluator_truth()["targets"]
    predictions = [{"row_id": row["row_id"], "prediction": row["value"], "unit": row["unit"]} for row in targets]
    result = score_synthetic_reference(full_reference, predictions)
    assert result["coverage"] == 1
    assert result["evaluation_kind"] == "synthetic_recovery"
    assert not result["independent_real_validation"]
    assert all(metric["rmse"] == 0 for metric in result["metrics"])
    partial = score_synthetic_reference(full_reference, predictions[:2])
    assert partial["n_scored"] == 2
    assert partial["coverage"] == 2 / len(targets)
    assert all(metric["rmse"] is None for metric in partial["metrics"] if metric["n_scored"] == 0)


def test_parameter_coordinate_recovery_does_not_collapse_unidentified_directions(full_reference):
    experiment = next(e for e in full_reference.evaluator_truth()["experiments"] if e["design"]["experiment_id"] == "test_structure")
    provenance = experiment["provenance"]
    estimate = {
        "experiment_id": "test_structure", "parameter_id": "kv17f_1",
        "estimate": provenance["reference_parameter_values"]["kv17f_1"],
        "unit": provenance["units"]["parameters"]["kv17f_1"],
    }
    result = score_parameter_recovery(full_reference, [estimate])
    assert result["n_scored"] == 1 and result["n_total"] == 8 * 98
    assert not result["unique_parameter_recovery_established"]
    assert not result["independent_real_validation"]
    row = next(row for row in result["rows"] if row["scored"])
    assert row["point_error"] == 0 and row["reference_value"] == 0
    assert row["resolved_uncertainty"] is None
    assert row["identifiability"]["unresolved_directions"]
    with pytest.raises(ValueError, match="units must match"):
        score_parameter_recovery(full_reference, [{**estimate, "unit": "s^-1"}])
    with pytest.raises(ValueError, match="duplicate"):
        score_parameter_recovery(full_reference, [estimate, estimate])


def test_benchmark_callback_receives_only_the_learner_projection(full_reference):
    seen = []

    def learner(inputs):
        seen.append(inputs)
        assert "targets" not in inputs and "experiments" not in inputs
        return training_mean_baseline(inputs)

    result = run_synthetic_benchmark(full_reference, learner)
    assert len(seen) == 1 and result["coverage"] == 1
    assert any(metric["rmse"] > 0 for metric in result["metrics"])
    assert "generating parameter vector" in result["identifiability"]


@pytest.mark.parametrize("malformation", ["extra", "duplicate", "unknown", "unit", "nan"])
def test_scoring_rejects_misalignment_truth_injection_and_nonfinite_predictions(full_reference, malformation):
    row = full_reference.learner_inputs()["prediction_queries"][0]
    prediction = {"row_id": row["row_id"], "prediction": 0.0, "unit": row["unit"]}
    predictions = [prediction]
    if malformation == "extra":
        prediction["truth"] = 1
    elif malformation == "duplicate":
        predictions.append(dict(prediction))
    elif malformation == "unknown":
        prediction["row_id"] = "training-row-or-unknown"
    elif malformation == "unit":
        prediction["unit"] = "invented_molar_calibration"
    else:
        prediction["prediction"] = float("nan")
    with pytest.raises(ValueError):
        score_synthetic_reference(full_reference, predictions)


def test_assay_gain_background_and_noise_never_change_biological_truth(small_design):
    first = generate_reference(freeze_reference(small_design, product_outcomes_seen=False))
    scaled = replace(small_design, assays=(AssayScale("Hog1PP_measured", gain=3, background=0.7, noise_sd=0),))
    second = generate_reference(freeze_reference(scaled, product_outcomes_seen=False))
    for a, b in zip(first.evaluator_truth()["experiments"], second.evaluator_truth()["experiments"], strict=True):
        np.testing.assert_array_equal(a["states"], b["states"])
        assert a["provenance"]["reference_parameter_values"] == b["provenance"]["reference_parameter_values"]
    a = [row["value"] for row in first.learner_inputs()["training_observations"]]
    b = [row["value"] for row in second.learner_inputs()["training_observations"]]
    np.testing.assert_allclose(b, 3 * np.asarray(a) + 0.7)
    assert first.frozen.sha256 != second.frozen.sha256


def test_replay_and_noise_are_deterministic_and_order_independent(small_design):
    design = replace(small_design, assays=(AssayScale("Hog1PP_measured", 1, 0, 0.1),))
    frozen = freeze_reference(design, product_outcomes_seen=False)
    first, second = generate_reference(frozen), generate_reference(frozen)
    assert first.learner_inputs() == second.learner_inputs()
    assert first.evaluator_truth() == second.evaluator_truth()
    reordered = replace(design, experiments=tuple(reversed(design.experiments)))
    third = generate_reference(freeze_reference(reordered, product_outcomes_seen=False))
    assert first.learner_inputs()["training_observations"] == third.learner_inputs()["training_observations"]
    assert first.evaluator_truth()["targets"] == third.evaluator_truth()["targets"]


def test_freeze_roundtrip_rejects_tampering_and_runtime_changes(small_design, monkeypatch):
    frozen = freeze_reference(small_design, product_outcomes_seen=False)
    restored = FrozenReference.from_dict(frozen.to_dict())
    assert restored.sha256 == frozen.sha256
    verify_reference_freeze(restored, load_parameter_evidence())
    payload = frozen.to_dict()
    payload["design"]["seed"] += 1
    with pytest.raises(ValueError, match="digest mismatch"):
        FrozenReference.from_dict(payload)
    monkeypatch.setattr(reference_module, "_implementation_identity", lambda: {"files": {"changed": "changed"}})
    with pytest.raises(ValueError, match="changed after freeze"):
        generate_reference(frozen)
    with pytest.raises(ValueError, match="no-product-outcomes"):
        freeze_reference(small_design, product_outcomes_seen=True)


def test_design_rejects_learner_coefficients_and_random_timepoint_holdout(small_design):
    payload = small_design.to_dict()
    payload["learner_coefficients"] = {"k": 1}
    with pytest.raises(ValueError, match="not inputs"):
        ReferenceDesign.from_dict(payload)
    repeated = replace(small_design.experiments[1], salt_history=small_design.experiments[0].salt_history)
    with pytest.raises(ValueError, match="cross train/test"):
        replace(small_design, experiments=(small_design.experiments[0], repeated))
    unavailable = replace(small_design.experiments[1], model_key="jalihal2021")
    with pytest.raises(EvidenceGap, match="only the pinned"):
        freeze_reference(replace(small_design, experiments=(small_design.experiments[0], unavailable)), product_outcomes_seen=False)


@pytest.mark.parametrize("changes", [
    {"times_s": (1, 20, 30)}, {"times_s": (0, 20, 20)}, {"times_s": (0, float("nan"), 30)},
    {"salt_history": (SaltChange(20, 0.4), SaltChange(22, 0.0))},
    {"salt_history": (SaltChange(59, 0.4),)}, {"structure": "made_up_full_cell"},
])
def test_invalid_schedules_and_structures_fail_closed(small_design, changes):
    with pytest.raises(ValueError):
        replace(small_design.experiments[0], **changes)


def test_identifiability_preserves_null_directions_and_units_instead_of_inverting_them():
    matrix = np.array([[1, 1, 0], [2, 2, 1]], dtype=float)
    result = sensitivity_identifiability(matrix, ["biological_gain", "assay_gain", "turnover"], ["1", "RFU/molecule", "h^-1"])
    assert result["rank"] == 2 and result["unresolved_dimension"] == 1
    assert result["covariance"] is None
    vector = np.array([result["unresolved_directions"][0][key] for key in result["parameter_ids"]])
    np.testing.assert_allclose(matrix @ vector, 0, atol=1e-14)
    assert abs(vector[0]) > 0.5 and abs(vector[1]) > 0.5
    zero = sensitivity_identifiability(np.zeros((2, 3)), ["a", "b", "c"], ["1"] * 3)
    assert zero["rank"] == 0 and zero["unresolved_dimension"] == 3


def test_real_replay_is_separate_and_refuses_uncalibrated_protein_scales(full_reference):
    result = evaluate_native_reference(full_reference.frozen, evaluation_kind="published_fit_replay")
    assert result["evaluation_kind"] == "real_native_published_fit_replay"
    assert result["synthetic_recovery"] is False and result["independent_real_validation"] is False
    assert result["no_refit"] is True
    assert result["n_total"] == 54 and result["n_scored"] == 28
    assert len(result["refused_observations"]) == 26
    assert {row["observable_id"] for row in result["rows"] if row["scored"]} == {"glycerol_measured", "glycerol_e"}
    assert all(row["original_fit_status"] == "published_fitting_dataset" for row in result["rows"])
    with pytest.raises(EvidenceGap, match="independent real validation"):
        evaluate_native_reference(full_reference.frozen, evaluation_kind="independent_validation")


def test_reference_imports_no_fitted_learner_engine_or_other_mechanistic_coefficients():
    tree = ast.parse(Path(reference_module.__file__).read_text(encoding="utf-8"))
    imports = {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("ystwin")}
    assert imports == {
        "ystwin.analysis.parameter_evidence", "ystwin.mech.kinetic_sbml", "ystwin.analysis.hog_data",
    }
    assert not any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {"eval", "exec", "__import__"} for node in ast.walk(tree))


def _runner():
    spec = importlib.util.spec_from_file_location("benchmark_mechanistic_reference_runner", ROOT / "scripts/benchmark_mechanistic_reference.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cli_inventory_reports_blockers_and_full_coverage_gate(capsys):
    runner = _runner()
    exit_code = runner.main(["inventory"])
    report = json.loads(capsys.readouterr().out)
    assert exit_code == (0 if all(c["verified"] for c in report["source_checks"]) else 2)
    assert not report["full_biological_coverage"]
    assert runner.main(["inventory", "--require-complete"]) == 2
    capsys.readouterr()


def test_cli_exports_separate_roles_never_overwrites_and_scores_frozen_predictions(small_design, tmp_path, capsys):
    runner = _runner()
    design_path = tmp_path / "design.json"
    design_path.write_text(json.dumps(small_design.to_dict()), encoding="utf-8")
    output = tmp_path / "benchmark"
    args = ["synthetic", "--design", str(design_path), "--declare-no-product-outcomes", "--output-dir", str(output)]
    assert runner.main(args) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["learner"] == "training_mean_baseline_not_a_mechanistic_learner"
    inputs_path = output / "learner/inputs.json"
    original = inputs_path.read_bytes()
    inputs = json.loads(original)
    assert "targets" not in inputs and "experiments" not in inputs
    assert runner.main(args) == 2
    assert "restamp" in capsys.readouterr().err
    assert inputs_path.read_bytes() == original
    targets = json.loads((output / "evaluator/truth.json").read_text())["targets"]
    predictions_path = tmp_path / "predictions.json"
    predictions_path.write_text(json.dumps([{"row_id": row["row_id"], "prediction": row["value"], "unit": row["unit"]} for row in targets]))
    assert runner.main(["score-synthetic", "--freeze", str(output / "evaluator/freeze.json"), "--predictions", str(predictions_path)]) == 0
    scored = json.loads(capsys.readouterr().out)
    assert scored["coverage"] == 1 and all(m["rmse"] == 0 for m in scored["metrics"])
    estimates_path = tmp_path / "estimates.json"
    estimates_path.write_text("[]", encoding="utf-8")
    assert runner.main(["score-parameters", "--freeze", str(output / "evaluator/freeze.json"), "--estimates", str(estimates_path)]) == 0
    coordinates = json.loads(capsys.readouterr().out)
    assert coordinates["coverage"] == 0 and coordinates["n_total"] == 98
    assert not coordinates["unique_parameter_recovery_established"]
