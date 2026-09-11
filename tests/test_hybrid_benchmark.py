from __future__ import annotations

import importlib.util
import json
from dataclasses import asdict, replace
from types import SimpleNamespace

import cobra
import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.analysis.hybrid_stress import fit_sensor_model
from ystwin.fba.product_panel import ProductTask
from ystwin.hybrid import Environment, ScenarioUnsupported, simulate_hybrid_case


@pytest.fixture(scope="module")
def fitted_sensor():
    frame = pd.read_csv(paths.REPO_ROOT / "outputs" / "sensor_characterisation.csv")
    plates = tuple(sorted(frame.plate.unique()))
    return fit_sensor_model(frame, training_plates=plates[:-1], seed=0)


@pytest.fixture
def metabolic_model():
    model = cobra.Model("hybrid_test")
    sugar_e = cobra.Metabolite("sugar_e", formula="C6H12O6", compartment="e")
    sugar_c = cobra.Metabolite("sugar_c", formula="C6H12O6", compartment="c")
    product = cobra.Metabolite("product_c", formula="C12H22O11", compartment="c")
    oxygen = cobra.Metabolite("oxygen_e", formula="O2", compartment="e")
    energy = cobra.Metabolite("energy_c", compartment="c")
    definitions = (
        ("r_1714", {sugar_e: -1}, (-10.0, 1000.0)),
        ("transport", {sugar_e: -1, sugar_c: 1}, (0.0, 1000.0)),
        ("r_2111", {sugar_c: -1}, (0.0, 1000.0)),
        ("synthesis", {sugar_c: -2, product: 1}, (0.0, 1000.0)),
        ("DM_product", {product: -1}, (0.0, 1000.0)),
        ("r_1992", {oxygen: -1}, (-20.0, 1000.0)),
        ("energy", {sugar_c: -1, energy: 1}, (0.0, 1000.0)),
        ("r_4046", {energy: -1}, (0.7, 1000.0)),
    )
    for name, metabolites, bounds in definitions:
        reaction = cobra.Reaction(name, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites(metabolites)
        model.add_reactions([reaction])
    model.objective = "r_2111"
    return model


@pytest.fixture
def product():
    return ProductTask("test_product", "DM_product", "product_c", 342.29648, 12.0,
                       "extra intracellular net accumulation", ("test network",))


def test_learned_state_drives_actual_fluxes_and_titre(fitted_sensor, metabolic_model, product):
    environment = Environment("DTT", stressor="DTT", dose_mM=1.0, hours=4.0)
    full = simulate_hybrid_case(fitted_sensor, metabolic_model, product, environment,
                               grid_points=5, steps=30)
    ablated = simulate_hybrid_case(fitted_sensor, metabolic_model, product, environment,
                                  couple_stress=False, grid_points=5, steps=30)

    assert full.stress.growth_retention < 0.95
    assert full.rate_info.growth_retention == full.stress.growth_retention
    assert ablated.rate_info.growth_retention == 1.0
    assert full.trajectory.titre_mg_per_l != pytest.approx(ablated.trajectory.titre_mg_per_l)
    assert "SIMULATED" in full.stress.source
    assert full.summary()["product_targets_used_for_fitting"] is False


def test_all_three_process_metrics_are_reported(fitted_sensor, metabolic_model, product):
    result = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                  Environment("baseline", hours=4.0), grid_points=5, steps=30)
    row = result.summary()

    assert row["titre_mg_per_l"] > 0
    assert row["mean_specific_rate_mmol_per_gdcw_h"] > 0
    assert row["yield_g_per_g"] > 0
    assert row["productivity_mg_per_l_h"] == pytest.approx(row["titre_mg_per_l"] / 4.0)
    assert 0 <= row["carbon_yield"] <= 1.0 + 1e-6
    assert row["prediction_basis"] == product.output_basis
    assert row["allocation_fraction"] == 0.25


def test_outcome_tables_are_not_opened_during_prediction(
        fitted_sensor, metabolic_model, product, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("prediction tried to read a table after training")

    monkeypatch.setattr(pd, "read_csv", forbidden)
    simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                         Environment("baseline", hours=2.0), grid_points=5, steps=20)


def test_out_of_calibration_stress_is_not_silently_used(fitted_sensor, metabolic_model, product):
    with pytest.raises(ScenarioUnsupported, match="calibrat"):
        simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                             Environment("lethal", stressor="DTT", dose_mM=10.0))


def test_inputs_and_predictions_are_repeatable(fitted_sensor, metabolic_model, product):
    bounds = {reaction.id: reaction.bounds for reaction in metabolic_model.reactions}
    environment = Environment("baseline", hours=2.0)
    first = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                 environment, grid_points=5, steps=20)
    second = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                  environment, grid_points=5, steps=20)

    np.testing.assert_array_equal(first.trajectory.product_mmol_per_l,
                                  second.trajectory.product_mmol_per_l)
    assert {reaction.id: reaction.bounds for reaction in metabolic_model.reactions} == bounds


def test_temperature_is_not_an_unused_input(fitted_sensor, metabolic_model, product):
    environment = Environment("baseline", hours=2.0)
    first = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                 environment, grid_points=5, steps=20)
    changed = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                   replace(environment, name="warm", temperature_c=37.0),
                                   grid_points=5, steps=20)

    assert first.rate_info.maximum_growth_per_h != changed.rate_info.maximum_growth_per_h
    assert first.trajectory.titre_mg_per_l != pytest.approx(changed.trajectory.titre_mg_per_l)


@pytest.mark.parametrize("fields", [
    {"hours": 0.0}, {"hours": float("nan")}, {"substrate_g_per_l": -1.0},
    {"initial_biomass_g_per_l": 0.0}, {"oxygen_uptake": float("inf")},
    {"carbon_uptake": -1.0}, {"carbon_source": "unknown"},
    {"stressor": None, "dose_mM": 1.0},
])
def test_invalid_environments_are_refused(fields):
    with pytest.raises(ValueError):
        Environment("invalid", **fields)


@pytest.fixture(scope="module")
def runner():
    spec = importlib.util.spec_from_file_location(
        "hybrid_benchmark_script", paths.REPO_ROOT / "scripts/run_hybrid_benchmark.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def external_targets():
    return pd.read_csv(paths.REPO_ROOT / "data/hybrid_benchmark/external_targets.csv")


def test_external_product_values_cannot_change_a_prediction(runner, external_targets, monkeypatch):
    targets = external_targets[(external_targets["product"] == "trehalose")
                               & (external_targets.case_id == "ph5_co2_0p04")].copy()
    task = ProductTask("trehalose", "demand", "trehalose", 342.29648, 12,
                       "extra accumulation", ())
    monkeypatch.setattr(runner, "install_product", lambda model, name: (model, task))
    monkeypatch.setattr(runner, "biomass_pool_quota", lambda model, name: 0.1)

    def rates(model, *args, **kwargs):
        assert "observed" not in kwargs
        return SimpleNamespace(at_uptake=lambda uptake: (kwargs["fixed_growth_per_h"], uptake, 0.002)), None

    monkeypatch.setattr(runner, "build_network_rates", rates)
    options = SimpleNamespace(allocation=0.25, growth_fraction=0.9, grid_points=5)
    first = runner._predict_external(None, targets, options)
    targets["observed"] *= 1000
    targets["uncertainty"] *= 1000
    second = runner._predict_external(None, targets, options)

    assert first.predicted.tolist() == second.predicted.tolist()
    expected = (0.1 + 0.002 * 0.97 / 0.026) * 342.29648
    assert first.predicted.iloc[0] == pytest.approx(expected)
    assert bool(first.comparable.iloc[0])


def test_missing_experimental_inputs_do_not_get_invented(runner, external_targets, monkeypatch):
    selected = external_targets[(external_targets.study_id == "hull2014")
                                | ((external_targets.study_id == "hakkaart2020")
                                   & external_targets.oxygen_uptake.isna())]

    def forbidden(*args, **kwargs):
        raise AssertionError("incomplete external case entered the metabolic solve")

    monkeypatch.setattr(runner, "build_network_rates", forbidden)
    options = SimpleNamespace(allocation=0.25, growth_fraction=0.9, grid_points=5)
    predictions = runner._predict_external(None, selected, options)

    assert len(predictions) == len(selected)
    assert predictions.prediction_status.eq("unsupported").all()
    assert predictions.predicted.isna().all()
    assert predictions.reason.str.len().gt(0).all()


def test_baselines_are_product_target_independent(runner, external_targets, monkeypatch):
    task = ProductTask("trehalose", "demand", "trehalose", 342.29648, 12,
                       "extra accumulation", ())
    monkeypatch.setattr(runner, "install_product", lambda model, name: (model, task))
    monkeypatch.setattr(runner, "biomass_pool_quota", lambda model, name: 0.1)
    baseline = runner._external_baselines(None, external_targets)

    assert baseline.loc[baseline["product"] == "glycerol", "predicted"].eq(0.0).all()
    assert baseline.loc[baseline["product"] == "trehalose", "predicted"].eq(34.229648).all()
    assert len(baseline) == len(external_targets)


def test_runner_does_not_overwrite_an_existing_run(runner, tmp_path):
    path = tmp_path / "sentinel.txt"
    path.write_text("keep this run")
    with pytest.raises(SystemExit, match="overwrite"):
        runner.main(["--output-dir", str(tmp_path)])
    assert path.read_text() == "keep this run"


def test_declared_environments_are_distinct_and_include_real_stress_calls(runner):
    environments = runner.environments(24.0)

    assert len(environments) == len({env.name for env in environments}) == 6
    assert {env.carbon_source for env in environments} == {"glucose", "ethanol"}
    assert {env.stressor for env in environments} == {None, "DTT", "H2O2"}


def test_supplied_observations_do_not_get_replaced_by_simulation(
        fitted_sensor, metabolic_model, product, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("observed features were replaced by simulated readings")

    monkeypatch.setattr(type(fitted_sensor), "predict_condition", forbidden)
    features = dict.fromkeys(fitted_sensor.feature_names, 1.0)
    result = simulate_hybrid_case(fitted_sensor, metabolic_model, product,
                                  Environment("observed", hours=1.0),
                                  observed_features=features, grid_points=5, steps=10)
    assert result.stress.source == "observed"


@pytest.fixture
def numerical_reference(runner, tmp_path):
    protocol = {
        "input_hashes": {"sensor": "data"}, "training_plates": ["one"],
        "heldout_plate": "two", "training_seed": 0, "allocation_fraction": 0.25,
        "growth_floor_fraction": 0.9, "environments": [{"name": "baseline"}],
        "code_hashes": {"src/ystwin/hybrid.py": "model", "scripts/run_hybrid_benchmark.py": "driver"},
    }
    frame = pd.DataFrame([{
        "product": "test", "environment": "baseline", "status": "predicted",
        "titre_mg_per_l": 10.0, "mean_specific_rate_mmol_per_gdcw_h": 0.5,
        "yield_g_per_g": 0.1,
    }])
    runner._write_json(tmp_path / "protocol.json", protocol)
    runner._write_json(tmp_path / "sensor_model.json", {"seed": 0})
    frame.to_csv(tmp_path / "scenarios.csv", index=False)
    return protocol, frame


def test_numerical_comparison_reports_changes_without_refitting(runner, numerical_reference, tmp_path):
    protocol, frame = numerical_reference
    frame["titre_mg_per_l"] *= 1.01
    destination = tmp_path / "comparison.csv"
    runner._compare_numerics(frame, protocol, tmp_path, destination)
    comparison = pd.read_csv(destination)

    assert comparison.titre_mg_per_l_change_pct.iloc[0] == pytest.approx(1.0)
    assert comparison.yield_g_per_g_change_pct.iloc[0] == 0.0


@pytest.mark.parametrize("field", ["allocation_fraction", "training_seed", "input_hashes"])
def test_numerical_comparison_rejects_changed_assumptions(
        runner, numerical_reference, tmp_path, field):
    protocol, frame = numerical_reference
    protocol[field] = {"sensor": "different"} if field == "input_hashes" else 0.5
    with pytest.raises(ValueError, match="comparison"):
        runner._compare_numerics(frame, protocol, tmp_path, tmp_path / "comparison.csv")


@pytest.fixture
def assessment_run(runner, tmp_path):
    directory = tmp_path / "reference"
    directory.mkdir()
    products = ["trehalose", "glycerol", "glutathione"]
    environment = asdict(runner.environments(24.0)[0])
    protocol = {
        "run_complete": True, "products": products, "environments": [environment],
        "scenario_rows": 3, "successful_scenarios": 3,
        "input_hashes": {"sensor_table": "sensor", "ec_model": "model", "proteome": "protein"},
        "code_hashes": {"src/ystwin/hybrid.py": "model-code",
                        "scripts/run_hybrid_benchmark.py": "old-driver"},
        "training_plates": ["training"], "heldout_plate": "heldout", "training_seed": 0,
        "allocation_fraction": 0.25, "growth_floor_fraction": 0.9,
        "product_targets_used_for_fitting": False, "grid_points": 5, "steps": 120,
    }
    runner._write_json(directory / "protocol.json", protocol)
    runner._write_json(directory / "sensor_model.json", {"seed": 0})
    pd.DataFrame([{
        "product": product, "environment": environment["name"], "status": "predicted",
        "titre_mg_per_l": 10.0, "mean_specific_rate_mmol_per_gdcw_h": 0.5,
        "yield_g_per_g": 0.1, "allocation_fraction": 0.25, "growth_fraction": 0.9,
        "couple_stress": True, "grid_points_requested": 5,
    } for product in products]).to_csv(directory / "scenarios.csv", index=False)
    return directory


def _assessment_external(runner, directory, model_values=(80.0, 160.0),
                         baseline_values=(30.0, 30.0), observations=(10.0, 20.0),
                         kinds=None, statuses=None):
    keys = ["study_id", "case_id", "product", "metric"]
    targets = pd.DataFrame([{
        "study_id": "external", "case_id": f"case_{index}", "product": "trehalose",
        "metric": "intracellular_content", "units": "mg/g_total_dry_biomass",
        "observed": observed, "detection_limit": np.nan,
        "observation_kind": "quantified" if kinds is None else kinds[index],
    } for index, observed in enumerate(observations)])
    frames = []
    for prefix, values in (("external", model_values), ("external_baseline", baseline_values)):
        predictions = targets[keys].copy()
        predictions["predicted"] = values
        predictions["prediction_units"] = targets.units
        predictions["comparable"] = True
        predictions["prediction_status"] = [
            "unsupported" if value is None else "predicted" for value in values]
        if prefix == "external" and statuses is not None:
            predictions["prediction_status"] = statuses
        predictions["prediction_basis"] = "frozen test assumption"
        predictions["reason"] = ""
        scores = runner.score_external_targets(targets, predictions)
        scores.to_csv(directory / f"{prefix}_scores.csv", index=False)
        runner.summarize_scores(scores).to_csv(directory / f"{prefix}_summary.csv", index=False)
        frames.append(scores)
    return frames


@pytest.fixture
def assessment_numerics(runner, assessment_run, tmp_path):
    directory = tmp_path / "refinement"
    directory.mkdir()
    protocol = json.loads((assessment_run / "protocol.json").read_text())
    protocol.update(products=["trehalose"], scenario_rows=1, successful_scenarios=1,
                    numerical_reference=str(assessment_run), grid_points=13, steps=480)
    protocol["code_hashes"]["scripts/run_hybrid_benchmark.py"] = "new-driver"
    runner._write_json(directory / "protocol.json", protocol)
    runner._write_json(directory / "sensor_model.json", {"seed": 0})
    frame = pd.read_csv(assessment_run / "scenarios.csv").iloc[:1].copy()
    for column in ("titre_mg_per_l", "mean_specific_rate_mmol_per_gdcw_h", "yield_g_per_g"):
        frame[column] *= 1.02
    frame["grid_points_requested"] = 13
    frame.to_csv(directory / "scenarios.csv", index=False)
    runner._compare_numerics(frame, protocol, assessment_run, directory / "numerical_refinement.csv")
    return directory


def test_assessment_completed_execution_is_not_scientific_validation(
        runner, assessment_run, monkeypatch):
    _assessment_external(runner, assessment_run)
    before = {path.name: path.read_bytes() for path in assessment_run.iterdir()}

    def forbidden(*args, **kwargs):
        raise AssertionError("assessment must not train, predict, or rewrite artifacts")

    for name in ("fit_sensor_model", "load_model", "_job_results", "_write_json"):
        monkeypatch.setattr(runner, name, forbidden)
    report = runner.assess_run(assessment_run)

    assert report["execution"]["run_complete"] is True
    assert report["execution"]["planned_scenarios"] == report["execution"]["actual_scenarios"] == 3
    assert report["validated_model"] is False
    assert report["scientific_result"] == "does_not_beat_baseline"
    product = report["external_validation"]["per_product"]["trehalose"]
    comparison = product["baseline_comparisons"][0]
    assert comparison["matched_scored_rows"] == 2
    assert comparison["model_mae"] == 105.0
    assert comparison["baseline_mae"] == 15.0
    assert comparison["model_median_abs_percent_error"] == 700.0
    assert comparison["result"] == "does_not_beat_baseline"
    assert report["external_validation"]["per_product"]["glutathione"]["status"] == "no_targets"
    assert report["external_validation"]["coverage_status"] == "incomplete"
    limitations = " ".join(report["evidence_limitations"])
    for phrase in ("allocation-conditioned", "empirical growth-cap", "gene", "ATP",
                   "paired observed sensor/product", "not confidence intervals"):
        assert phrase in limitations
    assert {path.name: path.read_bytes() for path in assessment_run.iterdir()} == before
    assert json.loads(json.dumps(report, allow_nan=False)) == report


def test_assessment_no_quantified_targets_is_not_a_zero_error_pass(runner, assessment_run):
    _assessment_external(runner, assessment_run, model_values=(0.0, 1.0),
                         baseline_values=(0.0, 0.0), observations=(0.0, 0.0),
                         kinds=["below_detection_limit"] * 2)
    report = runner.assess_run(assessment_run)
    product = report["external_validation"]["per_product"]["trehalose"]

    assert report["scientific_result"] == "insufficient_coverage"
    assert report["validated_model"] is False
    assert product["status"] == "no_quantified_targets"
    assert product["baseline_comparisons"][0]["matched_scored_rows"] == 0
    assert product["baseline_comparisons"][0]["model_mae"] is None
    assert product["baseline_comparisons"][0]["baseline_mae"] is None


def test_assessment_absent_run_and_missing_scenarios_are_explicit(runner, assessment_run, tmp_path):
    missing = runner.assess_run(tmp_path / "absent")
    assert missing["execution"]["status"] == "missing"
    assert missing["execution"]["run_complete"] is False
    assert missing["execution"]["planned_scenarios"] is None
    assert missing["external_validation"]["status"] == "missing"
    assert missing["scientific_result"] == "insufficient_coverage"
    (assessment_run / "scenarios.csv").unlink()
    incomplete = runner.assess_run(assessment_run)
    assert incomplete["execution"]["status"] == "incomplete"
    assert incomplete["execution"]["protocol_run_complete"] is True
    assert incomplete["execution"]["run_complete"] is False
    assert len(incomplete["execution"]["missing_scenarios"]) == 3


@pytest.mark.parametrize("partial_rows", [1, 3])
def test_assessment_incomplete_execution_is_not_inferred_complete(
        runner, assessment_run, partial_rows):
    protocol_path = assessment_run / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    protocol.update(run_complete=False, scenario_rows=partial_rows, successful_scenarios=partial_rows)
    runner._write_json(protocol_path, protocol)
    frame = pd.read_csv(assessment_run / "scenarios.csv").iloc[:partial_rows]
    frame.to_csv(assessment_run / "scenarios.csv", index=False)
    report = runner.assess_run(assessment_run)

    assert report["execution"]["status"] == "incomplete"
    assert report["execution"]["actual_scenarios"] == partial_rows
    assert report["execution"]["planned_scenarios"] == 3
    assert report["execution"]["run_complete"] is False


def test_assessment_joins_baselines_on_same_scored_target_keys(runner, assessment_run):
    _, baseline = _assessment_external(
        runner, assessment_run, model_values=(11.0, None, 1000.0),
        baseline_values=(15.0, 20.0, None), observations=(10.0, 10.0, 10.0))
    baseline.iloc[::-1].to_csv(assessment_run / "external_baseline_scores.csv", index=False)
    pd.DataFrame([{"product": "trehalose", "model_mae": 0, "baseline_mae": 999}]).to_csv(
        assessment_run / "external_baseline_comparison.csv", index=False)
    report = runner.assess_run(assessment_run)
    comparison = report["external_validation"]["per_product"]["trehalose"]["baseline_comparisons"][0]

    assert comparison["matched_scored_rows"] == 1
    assert comparison["model_only_scored_rows"] == comparison["baseline_only_scored_rows"] == 1
    assert comparison["model_mae"] == 1.0
    assert comparison["baseline_mae"] == 5.0
    assert comparison["result"] == "beats_baseline_on_scored_rows"
    assert comparison["matched_keys"][0]["case_id"] == "case_0"
    assert report["scientific_result"] == "insufficient_coverage"


def test_assessment_unknown_lod_and_infeasibility_are_not_double_counted(runner, assessment_run):
    _assessment_external(runner, assessment_run, model_values=(0.1, None),
                         baseline_values=(0.0, 0.0), observations=(0.0, 0.0),
                         kinds=["below_detection_limit"] * 2, statuses=["predicted", "infeasible"])
    report = runner.assess_run(assessment_run)
    external = report["external_validation"]

    assert external["target_rows"] == 2
    assert external["score_status_counts"]["censored"] == 1
    assert external["score_status_counts"]["infeasible"] == 1
    assert sum(external["score_status_counts"].values()) == 2
    assert external["observation_counts"]["below_detection_limit"] == 2
    assert external["observation_counts"]["unknown_detection_limit"] == 2
    assert external["summary_consistency"] == "verified"
    assert report["scientific_result"] == "insufficient_coverage"


def test_assessment_summary_only_cannot_establish_matched_baseline_success(runner, assessment_run):
    _assessment_external(runner, assessment_run, model_values=(10.0, 20.0))
    (assessment_run / "external_scores.csv").unlink()
    report = runner.assess_run(assessment_run)

    assert report["external_validation"]["status"] == "summary_only"
    assert report["external_validation"]["target_rows"] == 2
    assert report["external_validation"]["score_status_counts"]["scored"] == 2
    assert report["scientific_result"] == "insufficient_coverage"
    assert report["validated_model"] is False


@pytest.mark.parametrize("field", ["numerical_reference", "input_hashes", "allocation_fraction",
                                  "environments", "code_hashes", "training_seed"])
def test_assessment_rejects_unrelated_numerical_runs(
        runner, assessment_run, assessment_numerics, field):
    path = assessment_numerics / "protocol.json"
    protocol = json.loads(path.read_text())
    replacements = {
        "numerical_reference": str(assessment_run.parent / "unrelated"),
        "input_hashes": {"sensor_table": "unrelated"}, "allocation_fraction": 0.5,
        "environments": [{**protocol["environments"][0], "hours": 48.0}],
        "code_hashes": {"src/ystwin/hybrid.py": "different-model"}, "training_seed": 1,
    }
    protocol[field] = replacements[field]
    runner._write_json(path, protocol)

    with pytest.raises(runner.AssessmentError, match="numerical|seed"):
        runner.assess_run(assessment_run, numerical_run=assessment_numerics)


def test_assessment_numerical_maxima_do_not_override_bad_external_evidence(
        runner, assessment_run, assessment_numerics):
    _assessment_external(runner, assessment_run)
    report = runner.assess_run(assessment_run, numerical_run=assessment_numerics)
    numerics = report["numerical_refinement"]

    assert numerics["status"] == "assessed"
    assert numerics["compatibility"] == "verified"
    assert numerics["compared_scenarios"] == 1
    assert numerics["reference_scenarios_not_refined"] == 2
    assert numerics["max_abs_finite_change_pct"] == pytest.approx(2.0)
    assert numerics["any_failed"] is False
    assert numerics["any_zero_reference_changed"] is False
    assert numerics["resolution"]["reference"]["grid_points"] == [5]
    assert numerics["resolution"]["refined"]["steps"] == 480
    assert report["scientific_result"] == "does_not_beat_baseline"
    assert report["validated_model"] is False


def test_assessment_zero_reference_change_is_undefined_not_zero(runner, assessment_run, assessment_numerics):
    frame = pd.read_csv(assessment_run / "scenarios.csv")
    frame.loc[0, "titre_mg_per_l"] = 0.0
    frame.to_csv(assessment_run / "scenarios.csv", index=False)
    (assessment_numerics / "numerical_refinement.csv").unlink()
    report = runner.assess_run(assessment_run, numerical_run=assessment_numerics)
    numerics = report["numerical_refinement"]

    assert numerics["any_zero_reference_changed"] is True
    assert numerics["per_metric"]["titre_mg_per_l"]["max_abs_finite_change_pct"] is None
    assert numerics["per_metric"]["titre_mg_per_l"]["zero_reference_changes"] == 1
    assert json.loads(json.dumps(report, allow_nan=False)) == report


def test_assessment_failed_numerics_are_not_silently_dropped(runner, assessment_run, assessment_numerics):
    frame = pd.read_csv(assessment_numerics / "scenarios.csv")
    frame["status"] = "failed"
    frame.to_csv(assessment_numerics / "scenarios.csv", index=False)
    path = assessment_numerics / "protocol.json"
    protocol = json.loads(path.read_text())
    protocol["successful_scenarios"] = 0
    runner._write_json(path, protocol)
    (assessment_numerics / "numerical_refinement.csv").unlink()
    numerics = runner.assess_run(assessment_run, numerical_run=assessment_numerics)["numerical_refinement"]

    assert numerics["status"] == "incomplete"
    assert numerics["any_failed"] is True
    assert numerics["compared_scenarios"] == 0
    assert numerics["max_abs_finite_change_pct"] is None


@pytest.mark.parametrize("artifact", ["protocol", "scenario_keys", "summary", "baseline_targets", "refinement"])
def test_assessment_contradictory_artifacts_raise_named_errors(
        runner, assessment_run, assessment_numerics, artifact):
    scores, baseline = _assessment_external(runner, assessment_run)
    if artifact == "protocol":
        path = assessment_run / "protocol.json"
        protocol = json.loads(path.read_text())
        protocol["scenario_rows"] = 99
        runner._write_json(path, protocol)
    elif artifact == "scenario_keys":
        frame = pd.read_csv(assessment_run / "scenarios.csv")
        pd.concat([frame, frame.iloc[:1]]).to_csv(assessment_run / "scenarios.csv", index=False)
    elif artifact == "summary":
        path = assessment_run / "external_summary.csv"
        summary = pd.read_csv(path)
        summary["n_scored"] = 0
        summary.to_csv(path, index=False)
    elif artifact == "baseline_targets":
        baseline.loc[0, "observed"] = 1.0
        baseline.to_csv(assessment_run / "external_baseline_scores.csv", index=False)
    else:
        path = assessment_numerics / "numerical_refinement.csv"
        comparison = pd.read_csv(path)
        comparison["titre_mg_per_l_reference"] = 999.0
        comparison.to_csv(path, index=False)

    with pytest.raises(runner.AssessmentError):
        runner.assess_run(assessment_run, numerical_run=assessment_numerics)


def test_assessment_nonfinite_score_evidence_cannot_become_zero(runner, assessment_run):
    scores, _ = _assessment_external(runner, assessment_run)
    scores["abs_error"] = np.inf
    scores.to_csv(assessment_run / "external_scores.csv", index=False)
    (assessment_run / "external_summary.csv").unlink()
    report = runner.assess_run(assessment_run)
    comparison = report["external_validation"]["per_product"]["trehalose"]["baseline_comparisons"][0]

    assert comparison["matched_scored_rows"] == 0
    assert comparison["model_mae"] is None
    assert report["external_validation"]["nonfinite_scored_error_rows"] == 2
    assert report["scientific_result"] == "insufficient_coverage"
    assert json.loads(json.dumps(report, allow_nan=False)) == report


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_assessment_rejects_nonfinite_protocol_json(runner, assessment_run, value):
    path = assessment_run / "protocol.json"
    protocol = json.loads(path.read_text())
    protocol["allocation_fraction"] = value
    path.write_text(json.dumps(protocol))
    with pytest.raises(runner.AssessmentError, match="non.finite"):
        runner.assess_run(assessment_run)


def test_assessment_cli_writes_new_report_without_running_or_claiming_science_pass(
        runner, assessment_run, monkeypatch):
    _assessment_external(runner, assessment_run)

    def forbidden(*args, **kwargs):
        raise AssertionError("assessment CLI entered normal model execution")

    monkeypatch.setattr(runner, "fit_sensor_model", forbidden)
    monkeypatch.setattr(runner, "load_model", forbidden)
    assert runner.main(["--assess-run", str(assessment_run)]) == 0
    report_path = assessment_run / "assessment.json"
    contents = report_path.read_text()
    assert json.loads(contents)["scientific_result"] == "does_not_beat_baseline"
    with pytest.raises(SystemExit, match="overwrite"):
        runner.main(["--assess-run", str(assessment_run)])
    assert report_path.read_text() == contents


@pytest.mark.parametrize("destination", ["protocol.json", "scenarios.csv", "external_scores.csv"])
def test_assessment_output_cannot_overwrite_any_existing_artifact(runner, assessment_run, destination):
    _assessment_external(runner, assessment_run)
    path = assessment_run / destination
    before = path.read_bytes()
    with pytest.raises(SystemExit, match="overwrite"):
        runner.main(["--assess-run", str(assessment_run), "--assessment-output", str(path)])
    assert path.read_bytes() == before


def test_assessment_cli_accepts_numerics_and_an_explicit_new_output(
        runner, assessment_run, assessment_numerics, tmp_path):
    destination = tmp_path / "new-assessment.json"
    assert runner.main(["--assess-run", str(assessment_run),
                        "--numerical-run", str(assessment_numerics),
                        "--assessment-output", str(destination)]) == 0
    assert json.loads(destination.read_text())["numerical_refinement"]["compatibility"] == "verified"
    assert not (assessment_run / "assessment.json").exists()


def test_assessment_cli_does_not_relax_normal_run_output_requirement(runner):
    with pytest.raises(SystemExit):
        runner.main([])
    with pytest.raises(SystemExit):
        runner.main(["--numerical-run", "some-run"])


def test_assessment_missing_numerical_provenance_does_not_combine_runs(
        runner, assessment_run, assessment_numerics):
    path = assessment_numerics / "protocol.json"
    protocol = json.loads(path.read_text())
    del protocol["numerical_reference"]
    runner._write_json(path, protocol)
    numerics = runner.assess_run(assessment_run, numerical_run=assessment_numerics)["numerical_refinement"]

    assert numerics["status"] == numerics["compatibility"] == "unverified"
    assert numerics["max_abs_finite_change_pct"] is None
    assert numerics["any_failed"] is None
    assert numerics["compared_scenarios"] == 0


def test_assessment_historical_seed_fallback_does_not_invent_unrecorded_steps(
        runner, assessment_run, assessment_numerics):
    for directory in (assessment_run, assessment_numerics):
        path = directory / "protocol.json"
        protocol = json.loads(path.read_text())
        del protocol["steps"]
        del protocol["grid_points"]
        del protocol["training_seed"]
        runner._write_json(path, protocol)
    numerics = runner.assess_run(assessment_run, numerical_run=assessment_numerics)["numerical_refinement"]

    assert numerics["compatibility"] == "verified"
    assert numerics["resolution"]["status"] == "partially_recorded"
    assert numerics["resolution"]["reference"]["steps"] is None
    assert numerics["resolution"]["refined"]["steps"] is None
    assert numerics["resolution"]["reference"]["grid_points"] == [5]
    assert numerics["resolution"]["refined"]["grid_points"] == [13]


def test_assessment_nonfinite_numerical_values_are_flagged_and_json_safe(
        runner, assessment_run, assessment_numerics):
    path = assessment_numerics / "scenarios.csv"
    frame = pd.read_csv(path)
    frame["titre_mg_per_l"] = np.inf
    frame.to_csv(path, index=False)
    (assessment_numerics / "numerical_refinement.csv").unlink()
    report = runner.assess_run(assessment_run, numerical_run=assessment_numerics)
    numerics = report["numerical_refinement"]

    assert numerics["any_nonfinite_values"] is True
    assert numerics["per_metric"]["titre_mg_per_l"]["max_abs_finite_change_pct"] is None
    assert numerics["per_metric"]["titre_mg_per_l"]["nonfinite_value_rows"] == 1
    assert numerics["status"] == "assessed_with_undefined_changes"
    assert json.loads(json.dumps(report, allow_nan=False)) == report


def test_assessment_unquantified_and_invalid_rows_preserve_valid_numeric_rows(runner, assessment_run):
    scores, _ = _assessment_external(runner, assessment_run, model_values=(80.0, None),
                                    statuses=["predicted", "invalid"])
    scores["predicted"] = scores.predicted.astype(object)
    scores.loc[1, "predicted"] = "invalid numerical input"
    scores.to_csv(assessment_run / "external_scores.csv", index=False)
    report = runner.assess_run(assessment_run)

    assert report["external_validation"]["score_status_counts"]["invalid"] == 1
    assert report["external_validation"]["score_status_counts"]["scored"] == 1
    comparison = report["external_validation"]["per_product"]["trehalose"]["baseline_comparisons"][0]
    assert comparison["matched_scored_rows"] == 1
    assert comparison["model_mae"] == 70.0


def test_assessment_extreme_json_count_raises_a_named_error(runner, assessment_run):
    path = assessment_run / "protocol.json"
    protocol = json.loads(path.read_text())
    protocol["scenario_rows"] = 10 ** 400
    runner._write_json(path, protocol)

    with pytest.raises(runner.AssessmentError):
        runner.assess_run(assessment_run)


def test_assessment_missing_job_artifacts_prevent_execution_completion(runner, assessment_run):
    path = assessment_run / "protocol.json"
    protocol = json.loads(path.read_text())
    protocol["planned_jobs"] = 4
    runner._write_json(path, protocol)
    report = runner.assess_run(assessment_run)

    assert report["execution"]["actual_scenarios"] == report["execution"]["planned_scenarios"] == 3
    assert report["execution"]["actual_jobs"] == 3
    assert report["execution"]["job_census_status"] == "incomplete"
    assert report["execution"]["run_complete"] is False


def test_assessment_rejects_nonboolean_zero_reference_flags(runner, assessment_run, assessment_numerics):
    path = assessment_numerics / "numerical_refinement.csv"
    frame = pd.read_csv(path)
    frame["titre_mg_per_l_reference_zero_changed"] = 0
    frame.to_csv(path, index=False)

    with pytest.raises(runner.AssessmentError):
        runner.assess_run(assessment_run, numerical_run=assessment_numerics)


def test_assessment_output_refuses_dangling_symlinks(runner, assessment_run, tmp_path):
    destination = tmp_path / "report.json"
    target = tmp_path / "not-created.json"
    destination.symlink_to(target)
    with pytest.raises(SystemExit, match="overwrite"):
        runner.main(["--assess-run", str(assessment_run), "--assessment-output", str(destination)])
    assert destination.is_symlink()
    assert not target.exists()
