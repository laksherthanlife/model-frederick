from __future__ import annotations

import importlib.util
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pytest

from ystwin.analysis import mechanistic_reference as reference_module
from ystwin.analysis.mechanistic_reference import (
    AssayScale,
    FrozenReference,
    MultiSourceReferenceDesign,
    NativeAssayScale,
    NativeReferenceExperiment,
    NutrientChange,
    NutrientInputs,
    ReferenceExperiment,
    SaltChange,
    SourceCounterfactual,
    default_multisource_reference_design,
    evaluate_native_reference,
    freeze_reference,
    generate_reference,
    reference_design_from_dict,
    run_synthetic_benchmark,
    score_parameter_recovery,
    score_synthetic_reference,
    training_mean_baseline,
    verify_reference_freeze,
)
from ystwin.analysis.parameter_evidence import (
    EvidenceGap,
    JALIHAL_NORMALIZED_INPUTS,
    content_sha256,
    load_parameter_evidence,
    load_reference_model,
    simulate_native_reference,
    validate_inventory,
)
from ystwin.mech.kinetic_sbml import KineticModel

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def multi_design():
    hog_times = (0, 5, 10)
    native_times = (0, 0.1, 1, 2, 3, 4, 5, 6)
    rich = NutrientInputs(0.5, 0.5, 1, 0, 0)
    shifted = NutrientInputs(0, 0.5, 1, 0, 0)
    native = [
        NativeReferenceExperiment("native_train", "jalihal2021", native_times, (NutrientChange(0, rich),), "train"),
        NativeReferenceExperiment("native_combination", "jalihal2021", native_times,
                                  (NutrientChange(0, NutrientInputs(0.25, 1, 0.2, 0.5, 0.8)),), "test", holdout_axis="input_combination"),
        NativeReferenceExperiment("native_history", "jalihal2021", native_times,
                                  (NutrientChange(0, rich), NutrientChange(2, shifted), NutrientChange(4, rich)),
                                  "test", holdout_axis="exposure_history"),
        NativeReferenceExperiment("native_parameter", "jalihal2021", native_times, (NutrientChange(0, rich),), "test",
                                  holdout_axis="parameter", counterfactuals=(SourceCounterfactual("gammasnf", 0.5),)),
    ]
    for index, name in enumerate(("w_pka_camp", "w_mig_snf", "w_sch9_torc")):
        native.append(NativeReferenceExperiment(
            f"native_structure_{index}", "jalihal2021", native_times, (NutrientChange(0, rich),), "test",
            holdout_axis="structure", counterfactuals=(SourceCounterfactual(name, 0, "structural_zero"),),
        ))
    return MultiSourceReferenceDesign(
        (
            ReferenceExperiment("hog_train", "hog2013_wt", hog_times, (SaltChange(0, 0.4),), "train"),
            ReferenceExperiment("hog_test", "hog2013_wt", hog_times, (SaltChange(0, 0.8),), "test", holdout_axis="environment"),
            *native,
        ),
        (
            AssayScale("Hog1PP_measured", 1, 0, 0),
            *(NativeAssayScale("jalihal2021", programme, observable) for programme, observables in (
                ("carbon_pka_snf1", ("PKA", "Snf1", "Mig1")), ("nitrogen_tor", ("TORC1", "Sch9")),
            ) for observable in observables),
        ),
        seed=712,
    )


@pytest.fixture(scope="module")
def multi_reference(multi_design):
    return generate_reference(freeze_reference(multi_design, product_outcomes_seen=False))


def _experiment(reference, experiment_id):
    return next(e for e in reference.evaluator_truth()["experiments"] if e["design"]["experiment_id"] == experiment_id)


def _keys(value):
    if isinstance(value, dict):
        return set(value) | {key for child in value.values() for key in _keys(child)}
    if isinstance(value, list):
        return {key for child in value for key in _keys(child)}
    return set()


def _small(design):
    return replace(design, experiments=tuple(e for e in design.experiments
                                             if e.experiment_id in ("hog_train", "hog_test", "native_train", "native_combination")))


def test_default_is_declared_multisource_with_independent_all_five_input_axes():
    design = default_multisource_reference_design(seed=91)
    assert isinstance(reference_design_from_dict(design.to_dict()), MultiSourceReferenceDesign)
    assert MultiSourceReferenceDesign.from_dict(design.to_dict()) == design
    assert len(design.experiments) == 27
    native = [e for e in design.experiments if isinstance(e, NativeReferenceExperiment)]
    assert len(native) == 17
    assert {e.holdout_axis for e in native if e.split == "test"} == {
        "input_combination", "exposure_history", "parameter", "structure",
    }
    baseline = next(e for e in native if e.experiment_id == "nutrient_train_control").input_history[0].inputs
    for axis in JALIHAL_NORMALIZED_INPUTS:
        train = next(e for e in native if e.experiment_id == f"nutrient_train_{axis}")
        history = next(e for e in native if e.experiment_id == f"nutrient_test_history_{axis}")
        assert {key for key in asdict(baseline) if getattr(train.input_history[0].inputs, key) != getattr(baseline, key)} == {axis}
        assert history.input_history[0].inputs == history.input_history[-1].inputs == baseline
        assert history.input_history[1].inputs == train.input_history[0].inputs
    training_inputs = {e.input_history for e in native if e.split == "train"}
    assert all(e.input_history not in training_inputs for e in native if e.holdout_axis == "input_combination")
    assert all("time_s" not in _keys(asdict(e)) and "times_s" not in asdict(e) for e in native)


def test_freeze_binds_design_sources_split_contracts_and_counterfactuals_before_simulation(multi_design, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("simulation is not allowed before the design/source/input split freeze")

    monkeypatch.setattr(reference_module, "simulate_native_reference", forbidden)
    monkeypatch.setattr(KineticModel, "simulate", forbidden)
    frozen = freeze_reference(multi_design, product_outcomes_seen=False)
    payload = frozen.to_dict()
    assert payload["schema_version"] == 2
    assert payload["evidence"]["source_ids"] == ["hog2013", "jalihal2021"]
    assert payload["input_split"]["frozen_before_generation"] is True
    assert set(payload["input_split"]["train"]).isdisjoint(payload["input_split"]["test"])
    assert payload["input_split"]["experiments_sha256"] == content_sha256(multi_design.to_dict()["experiments"])
    contract = payload["source_contracts"]["jalihal2021"]
    assert contract["times_field"] == "times_native" and contract["physical_seconds_per_native_unit"] is None
    assert contract["unit_audit"]["time"]["status"] == "conflicting_declarations"
    assert set(contract["input_axes"]) == set(JALIHAL_NORMALIZED_INPUTS)
    assert contract["physical_input_calibration"] is None
    assert payload["source_counterfactuals"]["native_structure_0"][0]["reference_value"] == 0
    assert verify_reference_freeze(FrozenReference.from_dict(payload), load_parameter_evidence()) == multi_design
    assert not payload["scope"]["full_environment_coverage"]
    assert set(payload["scope"]["unsupported_programmes"]) == {"oxidative_yap1", "heat_hsf1_hsp70", "upr_hac1", "ph"}
    assert "temperature" in payload["scope"]["unsupported_axes"]
    with pytest.raises(ValueError, match="no-product-outcomes"):
        freeze_reference(multi_design, product_outcomes_seen=True)


def test_actual_source_replay_matches_all_25_states_and_rates_without_rescaling(multi_reference):
    experiment = _experiment(multi_reference, "native_train")
    source = load_reference_model("jalihal2021")
    expected = simulate_native_reference("jalihal2021", experiment["times_native"],
                                         normalized_inputs=experiment["design"]["input_history"][0]["inputs"])
    assert experiment["species_ids"] == list(source.species_ids)
    assert experiment["reaction_ids"] == list(source.reaction_ids)
    assert np.shape(experiment["states"]) == (8, 25)
    np.testing.assert_array_equal(experiment["states"], expected.states)
    np.testing.assert_array_equal(experiment["reaction_rates"], expected.reaction_rates)
    provenance = experiment["provenance"]
    assert len(provenance["source_parameter_values"]) == 128
    assert len(provenance["reference_parameter_values"]) == 123
    assert provenance["source_parameter_values"] == source.parameter_values
    assert provenance["native_reference"]["source_bytes_modified"] is False
    assert provenance["execution_interface"] == "simulate_native_reference"
    assert provenance["imported_learner_coefficients"] is False
    for e in multi_reference.evaluator_truth()["experiments"]:
        assert np.isfinite(e["states"]).all() and np.isfinite(e["reaction_rates"]).all()


@pytest.mark.parametrize("experiment_id,parameter,observable,factor", [
    ("native_parameter", "gammasnf", "Snf1", 0.5),
    ("native_structure_0", "w_pka_camp", "PKA", 0),
    ("native_structure_1", "w_mig_snf", "Mig1", 0),
    ("native_structure_2", "w_sch9_torc", "Sch9", 0),
])
def test_actual_counterfactual_response_changes_only_the_declared_source_coordinate(multi_reference, experiment_id, parameter, observable, factor):
    base = _experiment(multi_reference, "native_train")
    altered = _experiment(multi_reference, experiment_id)
    source_values = base["provenance"]["reference_parameter_values"]
    reference_values = altered["provenance"]["reference_parameter_values"]
    assert {key for key in source_values if source_values[key] != reference_values[key]} == {parameter}
    assert reference_values[parameter] == source_values[parameter] * factor
    change = altered["provenance"]["source_counterfactuals"][0]
    assert change["status"] == "declared_source_counterfactual_not_biological_gene_intervention"
    assert change["unit"] == "unspecified" and change["factor_unit"] == "dimensionless"
    index = base["species_ids"].index(observable)
    assert np.max(np.abs(np.asarray(base["states"])[:, index] - np.asarray(altered["states"])[:, index])) > 1e-5
    model = load_reference_model("jalihal2021")
    parameters = {**base["design"]["input_history"][0]["inputs"], parameter: reference_values[parameter]}
    expected = model.simulate(altered["times_native"], parameters=parameters, method="Radau", rtol=1e-10, atol=1e-12)
    np.testing.assert_array_equal(altered["states"], expected.states)
    assert altered["provenance"]["identifiability"]["covariance"] is None


def test_history_carries_state_and_applies_right_continuous_inputs_without_coupling_atp(multi_reference):
    baseline = _experiment(multi_reference, "native_train")
    history = _experiment(multi_reference, "native_history")
    segments = history["provenance"]["segments"]
    for before, after in zip(segments, segments[1:]):
        np.testing.assert_array_equal(before["terminal_state"], after["initial_state"])
    np.testing.assert_allclose(history["states"][:4], baseline["states"][:4], rtol=1e-8, atol=1e-10)
    assert not np.allclose(history["states"][-1], baseline["states"][-1], rtol=1e-5, atol=1e-8)
    queries = [row for row in multi_reference.learner_inputs()["prediction_queries"] if row["experiment_id"] == "native_history"]
    at_change = next(row for row in queries if row["time_native"] == 2)
    assert at_change["normalized_inputs"]["Carbon"] == 0 and at_change["normalized_inputs"]["ATP"] == 0.5
    model = load_reference_model("jalihal2021")
    index = history["times_native"].index(2)
    expected_rates = model.reaction_rates(0, history["states"][index], at_change["normalized_inputs"])
    np.testing.assert_allclose(history["reaction_rates"][index], expected_rates, rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("axis,value", [("Carbon", 0), ("ATP", 0), ("Glutamine_ext", 0), ("NH4", 1), ("Proline", 1)])
def test_each_declared_nutrient_axis_executes_an_independent_heldout_history(multi_design, multi_reference, axis, value):
    train = next(e for e in multi_design.experiments if e.experiment_id == "native_train")
    baseline = train.input_history[0].inputs
    varied = replace(baseline, **{axis: value})
    test = replace(train, experiment_id=f"heldout_{axis}", split="test", holdout_axis="exposure_history",
                   input_history=(NutrientChange(0, baseline), NutrientChange(2, varied), NutrientChange(4, baseline)))
    design = replace(multi_design, experiments=(*multi_design.experiments[:2], train, test))
    reference = generate_reference(freeze_reference(design, product_outcomes_seen=False))
    experiment = _experiment(reference, test.experiment_id)
    original = _experiment(multi_reference, "native_train")
    assert np.max(np.abs(np.asarray(experiment["states"]) - np.asarray(original["states"]))) > 1e-7
    assert experiment["provenance"]["reference_parameter_values"] == original["provenance"]["reference_parameter_values"]
    rows = [r for r in reference.learner_inputs()["prediction_queries"] if r["experiment_id"] == test.experiment_id]
    for row in rows:
        expected = asdict(varied if 2 <= row["time_native"] < 4 else baseline)
        assert row["normalized_inputs"] == expected
        assert {key for key in expected if expected[key] != getattr(baseline, key)} <= {axis}


def test_off_grid_knots_and_segments_without_samples_still_preserve_state(multi_design):
    experiment = next(e for e in multi_design.experiments if e.experiment_id == "native_history")
    experiment = replace(experiment, times_native=(0, 0.1, 1, 6), input_history=(
        experiment.input_history[0], replace(experiment.input_history[1], time_native=2.1),
        replace(experiment.input_history[2], time_native=2.2),
    ))
    trajectory, provenance = reference_module._simulate_native_history(experiment, load_parameter_evidence())
    state = None
    for segment in provenance["segments"]:
        expected = simulate_native_reference("jalihal2021", [0, segment["end_native"] - segment["start_native"]],
                                             normalized_inputs=segment["normalized_inputs"], initial_state=state)
        np.testing.assert_allclose(segment["initial_state"], expected.states[0], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(segment["terminal_state"], expected.states[-1], rtol=1e-12, atol=1e-12)
        state = expected.states[-1]
    np.testing.assert_allclose(trajectory.states[-1], state, rtol=1e-12, atol=1e-12)


def test_native_output_clock_truth_firewall_and_frozen_observations_are_separate(multi_reference):
    inputs = multi_reference.learner_inputs()
    truth = multi_reference.evaluator_truth()
    forbidden = {"states", "reaction_rates", "provenance", "source_parameter_values", "reference_parameter_values", "factor", "gain",
                 "assay_scales", "counterfactuals", "source_counterfactuals", "expectation", "holdout_axis", "structure", "source_contracts"}
    assert not (_keys(inputs) & forbidden)
    assert all("value" not in row for row in inputs["prediction_queries"])
    assert set(row["experiment_id"] for row in inputs["training_observations"]).isdisjoint(row["experiment_id"] for row in inputs["prediction_queries"])
    assert truth["training_observations_sha256"] == content_sha256(inputs["training_observations"])
    assert truth["prediction_queries_sha256"] == content_sha256(inputs["prediction_queries"])
    for experiment in truth["experiments"]:
        if experiment["source_id"] == "jalihal2021":
            assert not (_keys(experiment) & {"times_s", "time_s", "time_h", "times_h", "hours"})
            assert experiment["time_unit"] == "source_native_unresolved"
            assert experiment["provenance"]["units"]["physical_seconds_per_native_unit"] is None
    native_rows = [row for row in inputs["prediction_queries"] if row["source_id"] == "jalihal2021"]
    assert all(row["time_unit"] == "source_native_unresolved" and "genotype" not in row for row in native_rows)
    assert all(set(row["normalized_inputs"]) == set(JALIHAL_NORMALIZED_INPUTS) for row in native_rows)
    inputs["training_observations"][0]["value"] = -999
    truth["experiments"][0]["states"][0][0] = -999
    assert multi_reference.learner_inputs()["training_observations"][0]["value"] != -999
    assert multi_reference.evaluator_truth()["experiments"][0]["states"][0][0] != -999


def test_multisource_observation_scales_and_noise_are_separate_reproducible_and_order_independent(multi_design):
    design = _small(multi_design)
    first = generate_reference(freeze_reference(design, product_outcomes_seen=False))
    noisy = replace(design, assays=tuple(replace(a, gain=3, background=0.7, noise_sd=0.01) if isinstance(a, NativeAssayScale) else a
                                         for a in design.assays))
    second = generate_reference(freeze_reference(noisy, product_outcomes_seen=False))
    third = generate_reference(freeze_reference(replace(noisy, experiments=tuple(reversed(noisy.experiments))), product_outcomes_seen=False))
    for a, b in zip(first.evaluator_truth()["experiments"], second.evaluator_truth()["experiments"], strict=True):
        np.testing.assert_array_equal(a["states"], b["states"])
    a = next(row for row in first.evaluator_truth()["targets"] if row["source_id"] == "jalihal2021")
    b = next(row for row in second.evaluator_truth()["targets"] if row["row_id"] == a["row_id"])
    assert b["expectation"] == pytest.approx(3 * a["expectation"] + 0.7)
    assert first.frozen.sha256 != second.frozen.sha256
    assert second.learner_inputs()["training_observations"] == third.learner_inputs()["training_observations"]
    assert second.evaluator_truth()["targets"] == third.evaluator_truth()["targets"]
    assert second.evaluator_truth() == generate_reference(second.frozen).evaluator_truth()


def test_per_programme_scoring_callback_and_missing_predictions_do_not_pool_units(multi_reference):
    targets = multi_reference.evaluator_truth()["targets"]
    predictions = [{"row_id": row["row_id"], "prediction": row["value"], "unit": row["unit"]} for row in targets]
    result = score_synthetic_reference(multi_reference, predictions)
    assert set(result["programme_scores"]) == {"osmotic_hog", "carbon_pka_snf1", "nitrogen_tor"}
    assert all(metric["rmse"] == 0 for metric in result["metrics"])
    assert all(score["coverage"] == 1 and "rmse" not in score for score in result["programme_scores"].values())
    assert all(score["evaluation_kind"] == "synthetic_recovery" for score in result["programme_scores"].values())
    assert result["full_biology_claim"] is False and result["independent_real_validation"] is False
    ids = {row["row_id"] for row in targets if row["programme"] == "nitrogen_tor"}
    partial = score_synthetic_reference(multi_reference, [p for p in predictions if p["row_id"] in ids])
    assert partial["programme_scores"]["nitrogen_tor"]["coverage"] == 1
    assert partial["programme_scores"]["carbon_pka_snf1"]["coverage"] == 0
    seen = []

    def learner(inputs):
        assert not (_keys(inputs) & {"targets", "experiments", "frozen", "counterfactuals"})
        seen.append(inputs)
        return training_mean_baseline(inputs)

    baseline = run_synthetic_benchmark(multi_reference, learner)
    assert len(seen) == 1 and baseline["coverage"] == 1
    assert any(metric["rmse"] > 0 for metric in baseline["metrics"])
    with pytest.raises(ValueError, match="units differ"):
        score_synthetic_reference(multi_reference, [{**predictions[-1], "unit": "mol/L"}])


def test_mean_baseline_never_borrows_observations_across_sources_or_programmes():
    rows = [
        {"source_id": "a", "programme": "p", "observable_id": "x", "unit": "u", "value": 2},
        {"source_id": "b", "programme": "p", "observable_id": "x", "unit": "u", "value": 40},
        {"source_id": "a", "programme": "q", "observable_id": "x", "unit": "u", "value": 8},
    ]
    queries = [{**{k: v for k, v in row.items() if k != "value"}, "row_id": str(index)} for index, row in enumerate(rows)]
    queries.append({**queries[0], "source_id": "unobserved", "row_id": "missing"})
    result = training_mean_baseline({"domain": "synthetic_native_reference", "training_observations": rows, "prediction_queries": queries})
    assert [r["prediction"] for r in result] == [2, 40, 8, None]


def test_native_parameter_diagnostic_excludes_imposed_inputs_and_preserves_unresolved_directions(multi_reference):
    altered = _experiment(multi_reference, "native_parameter")
    estimate = {"experiment_id": "native_parameter", "parameter_id": "gammasnf", "unit": "unspecified",
                "estimate": altered["provenance"]["reference_parameter_values"]["gammasnf"]}
    result = score_parameter_recovery(multi_reference, [estimate])
    assert result["n_total"] == 98 + 6 * 123 and result["n_scored"] == 1
    selected = next(row for row in result["rows"] if row["scored"])
    assert selected["point_error"] == 0 and selected["resolved_uncertainty"] is None
    assert selected["identifiability"]["unresolved_directions"]
    assert not result["unique_parameter_recovery_established"]
    with pytest.raises(ValueError, match="unknown"):
        score_parameter_recovery(multi_reference, [{**estimate, "parameter_id": "Carbon"}])


def test_real_replay_remains_hog_only_with_explicit_native_programme_refusals(multi_reference):
    real = evaluate_native_reference(multi_reference.frozen, evaluation_kind="published_fit_replay")
    assert real["evaluation_kind"] == "real_native_published_fit_replay"
    assert real["n_total"] == 54 and real["n_scored"] == 28
    assert real["synthetic_recovery"] is False and real["independent_real_validation"] is False
    assert set(real["programme_scores"]) == {"osmotic_hog"}
    assert set(real["refused_programmes"]) == {"carbon_pka_snf1", "nitrogen_tor"}
    assert all(row["original_fit_status"] == "published_fitting_dataset" for row in real["rows"])
    with pytest.raises(EvidenceGap, match="independent real validation"):
        evaluate_native_reference(multi_reference.frozen, evaluation_kind="independent_validation")


@pytest.mark.parametrize("field", ["design", "input_split", "source_contracts", "source_counterfactuals"])
def test_freeze_tampering_is_rejected(multi_reference, field):
    payload = multi_reference.frozen.to_dict()
    payload[field] = {}
    with pytest.raises(ValueError, match="digest mismatch"):
        FrozenReference.from_dict(payload)


def test_inventory_and_runtime_drift_refused_before_generation(multi_reference, monkeypatch):
    inventory = load_parameter_evidence().to_dict()
    inventory["sources"]["jalihal2021"]["verification"]["finding"] += " changed"
    with pytest.raises(ValueError, match="freeze no longer matches"):
        generate_reference(multi_reference.frozen, inventory=validate_inventory(inventory))
    monkeypatch.setattr(reference_module, "_implementation_identity", lambda: {"changed": True})
    with pytest.raises(ValueError, match="changed after freeze"):
        generate_reference(multi_reference.frozen)


@pytest.mark.parametrize("axis,value", [("Carbon", 20), ("ATP", True), ("NH4", -0.1), ("Glutamine_ext", float("nan")), ("Proline", "1")])
def test_invalid_normalized_inputs_fail_closed(axis, value):
    with pytest.raises(ValueError, match="normalized"):
        NutrientInputs(**{**asdict(NutrientInputs(0.5, 0.5, 0, 0, 0)), axis: value})


@pytest.mark.parametrize("field", ["learner_coefficients", "product_outcomes", "physical_environment", "genotype"])
def test_design_forbids_unseen_product_gene_or_learner_fields(multi_design, field):
    payload = multi_design.to_dict()
    payload[field] = {}
    with pytest.raises(ValueError, match="not inputs"):
        MultiSourceReferenceDesign.from_dict(payload)


@pytest.mark.parametrize("time_unit", ["s", "min", "h"])
def test_unresolved_source_clock_cannot_be_declared_physical(multi_design, time_unit):
    with pytest.raises(EvidenceGap, match="times_native"):
        replace(multi_design.experiments[2], time_unit=time_unit)


@pytest.mark.parametrize("source", ["zheng2016", "goulev2017", "unknown_source"])
def test_unsupported_source_dispatch_is_explicit_not_invented(multi_design, source):
    design = replace(multi_design,
                     experiments=tuple(replace(e, model_key=source) if isinstance(e, NativeReferenceExperiment) else e for e in multi_design.experiments),
                     assays=tuple(replace(a, model_key=source) if isinstance(a, NativeAssayScale) else a for a in multi_design.assays))
    with pytest.raises(EvidenceGap, match="no audited native nutrient dispatch"):
        freeze_reference(design, product_outcomes_seen=False)


@pytest.mark.parametrize("parameter,kind", [("invented_rate", "parameter_scale"), ("ATP", "parameter_scale"), ("gammasnf", "structural_zero")])
def test_unsupported_parameter_or_structural_coordinates_fail_before_simulation(multi_design, parameter, kind):
    changed = replace(multi_design.experiments[-1], counterfactuals=(SourceCounterfactual(parameter, 0, kind),),
                      holdout_axis="structure" if kind == "structural_zero" else "parameter")
    with pytest.raises(EvidenceGap):
        freeze_reference(replace(multi_design, experiments=(*multi_design.experiments[:-1], changed)), product_outcomes_seen=False)


def test_unsupported_assay_programme_and_physical_calibration_are_refused(multi_design):
    with pytest.raises(EvidenceGap, match="physical assay calibration"):
        NativeAssayScale("jalihal2021", "carbon_pka_snf1", "PKA", unit="mol/L")
    changed = replace(multi_design.assays[-1], programme="heat_hsf1_hsp70", observable_id="Hsf1")
    with pytest.raises(EvidenceGap, match="unsupported native programme"):
        freeze_reference(replace(multi_design, assays=(*multi_design.assays[:-1], changed)), product_outcomes_seen=False)


def test_physical_or_partial_input_dicts_and_hidden_timepoint_splits_are_rejected(multi_design):
    for inputs in ({"Carbon": 0.5}, {**asdict(NutrientInputs(0.5, 0.5, 0, 0, 0)), "glucose_g_L": 20}):
        payload = multi_design.to_dict()
        payload["experiments"][2]["input_history"][0]["inputs"] = inputs
        with pytest.raises(ValueError, match="all five normalized"):
            MultiSourceReferenceDesign.from_dict(payload)
    same = replace(multi_design.experiments[2], experiment_id="relabelled", times_native=(0, 0.2, 2, 5),
                   split="test", holdout_axis="input_combination")
    with pytest.raises(ValueError, match="cross train/test"):
        replace(multi_design, experiments=(*multi_design.experiments, same))
    rich = multi_design.experiments[2].input_history[0].inputs
    with pytest.raises(ValueError, match="redundant"):
        replace(multi_design.experiments[2], input_history=(NutrientChange(0, rich), NutrientChange(2, rich)))
    with pytest.raises(EvidenceGap, match="biological genotypes"):
        replace(same, holdout_axis="genotype")


def test_counterfactual_labels_or_order_cannot_hide_duplicate_source_experiments(multi_design):
    original = multi_design.experiments[-1]
    train = replace(original, split="train", holdout_axis=None, experiment_id="train_zero",
                    counterfactuals=(SourceCounterfactual("w_sch9_torc", 0),))
    with pytest.raises(ValueError, match="cross train/test"):
        replace(multi_design, experiments=(*multi_design.experiments, train))
    with pytest.raises(ValueError, match="must change"):
        SourceCounterfactual("w_pka_camp", 1)


def _runner():
    spec = importlib.util.spec_from_file_location("multisource_benchmark_runner", ROOT / "scripts/benchmark_mechanistic_reference.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cli_multisource_export_score_and_no_restamping(multi_design, tmp_path, capsys):
    runner = _runner()
    design = _small(multi_design)
    design_path = tmp_path / "design.json"
    design_path.write_text(json.dumps(design.to_dict()))
    output = tmp_path / "reference"
    args = ["synthetic", "--design", str(design_path), "--declare-no-product-outcomes", "--output-dir", str(output)]
    assert runner.main(args) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["source_ids"] == ["hog2013", "jalihal2021"]
    assert summary["learner"] == "training_mean_baseline_not_a_mechanistic_learner"
    assert summary["source_contracts"]["jalihal2021"]["times_field"] == "times_native"
    learner_path = output / "learner/inputs.json"
    original = learner_path.read_bytes()
    inputs = json.loads(original)
    assert "targets" not in _keys(inputs) and "counterfactuals" not in _keys(inputs)
    assert runner.main(args) == 2
    assert "restamp" in capsys.readouterr().err
    assert learner_path.read_bytes() == original
    targets = json.loads((output / "evaluator/truth.json").read_text())["targets"]
    predictions = tmp_path / "predictions.json"
    predictions.write_text(json.dumps([{"row_id": r["row_id"], "unit": r["unit"], "prediction": r["value"]} for r in targets]))
    assert runner.main(["score-synthetic", "--freeze", str(output / "evaluator/freeze.json"), "--predictions", str(predictions)]) == 0
    score = json.loads(capsys.readouterr().out)
    assert score["coverage"] == 1 and all(m["rmse"] == 0 for m in score["metrics"])
    assert set(score["programme_scores"]) == {"osmotic_hog", "carbon_pka_snf1", "nitrogen_tor"}


def test_cli_multisource_flag_selects_new_design_without_changing_legacy_default(multi_design, monkeypatch, capsys):
    runner = _runner()
    seen = []

    def selected(*, seed):
        seen.append(seed)
        return replace(_small(multi_design), seed=seed)

    monkeypatch.setattr(runner, "default_multisource_reference_design", selected)
    assert runner.main(["synthetic", "--multi-source", "--seed", "83", "--declare-no-product-outcomes"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert seen == [83] and result["source_ids"] == ["hog2013", "jalihal2021"]
    assert result["coverage"] == 1
