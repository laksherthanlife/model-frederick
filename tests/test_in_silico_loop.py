from __future__ import annotations

import importlib.util
import json
from dataclasses import fields, replace
from types import SimpleNamespace

import cobra
import numpy as np
import pandas as pd
import pytest

from ystwin import paths

from ystwin.fba.product_panel import ProductTask
from ystwin.in_silico import (
    InSilicoEnvironment,
    MetabolicControl,
    apply_synthetic_control,
    simulate_controlled_product,
)


@pytest.fixture
def network():
    model = cobra.Model("synthetic_control_test")
    carbon_e = cobra.Metabolite("carbon_e", formula="C6H12O6", compartment="e")
    carbon_c = cobra.Metabolite("carbon_c", formula="C6H12O6", compartment="c")
    product = cobra.Metabolite("product_c", formula="C12H22O11", compartment="c")
    oxygen = cobra.Metabolite("oxygen_e", formula="O2", compartment="e")
    energy = cobra.Metabolite("energy_c", compartment="c")
    protein = cobra.Metabolite("prot_pool[c]", compartment="c")
    definitions = (
        ("r_1714", {carbon_e: -1}, (-10.0, 1000.0)),
        ("transport", {carbon_e: -1, carbon_c: 1}, (0.0, 1000.0)),
        ("r_2111", {carbon_c: -1, protein: -0.05}, (0.0, 1000.0)),
        ("synthesis", {carbon_c: -2, product: 1, protein: -0.1}, (0.0, 1000.0)),
        ("DM_product", {product: -1}, (0.0, 1000.0)),
        ("r_1992", {oxygen: -1}, (-20.0, 1000.0)),
        ("energy", {carbon_c: -1, energy: 1}, (0.0, 1000.0)),
        ("r_4046", {energy: -1}, (0.7, 0.7)),
        ("prot_pool_exchange", {protein: 1}, (0.0, 0.2)),
    )
    for name, metabolites, bounds in definitions:
        reaction = cobra.Reaction(name, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites(metabolites)
        model.add_reactions([reaction])
    model.objective = "r_2111"
    return model


@pytest.fixture
def product():
    return ProductTask("test", "DM_product", "product_c", 342.29648, 12.0,
                       "synthetic variable intracellular pool", ())


def test_explicit_synthetic_control_changes_constraints_not_the_source(network):
    before = {reaction.id: reaction.bounds for reaction in network.reactions}
    control = MetabolicControl(0.8, 0.6, 1.4, 0.2)
    controlled, record = apply_synthetic_control(network, control, structural_mass_fraction=0.8)

    assert controlled.reactions.r_4046.bounds == pytest.approx((1.75, 1.75))
    assert controlled.reactions.prot_pool_exchange.upper_bound == pytest.approx(0.12)
    assert {reaction.id: reaction.bounds for reaction in network.reactions} == before
    assert record["control_source"] == "declared or learned in-silico controller"


@pytest.mark.parametrize("coefficient", [-2.0, -0.5, 0.5, 2.0])
def test_protein_budget_control_is_invariant_to_exchange_units(network, coefficient):
    pool = network.reactions.prot_pool_exchange
    metabolite = next(iter(pool.metabolites))
    pool.add_metabolites({metabolite: coefficient}, combine=False)
    pool.bounds = (0.0, 0.2 / coefficient) if coefficient > 0 else (0.2 / coefficient, 0.0)
    original = pool.bounds
    controlled, record = apply_synthetic_control(network, MetabolicControl(0.8, 0.6, 1.4, 0.2))
    changed = controlled.reactions.prot_pool_exchange
    endpoint = changed.upper_bound if coefficient > 0 else changed.lower_bound
    assert coefficient * endpoint == pytest.approx(0.12)
    assert record["protein_supply_applied"] == pytest.approx(0.12)
    assert pool.bounds == original


def test_identical_controls_reproduce_identical_teacher_and_student_rollouts(network, product):
    controls = np.array([[1.0, 1.0, 0.7, 0.2], [0.8, 0.8, 1.0, 0.3]])
    environment = InSilicoEnvironment("glucose", hours=4.0)
    first = simulate_controlled_product(network, product, [0.0, 2.0, 4.0], controls,
                                        environment, grid_points=5, output_dt=0.1)
    second = simulate_controlled_product(network, product, [0.0, 2.0, 4.0], controls,
                                         environment, grid_points=5, output_dt=0.1)

    np.testing.assert_array_equal(first.product_mmol_per_l, second.product_mmol_per_l)
    assert first.summary()["titre_mg_per_l"] > 0
    assert first.summary()["yield_g_per_g"] > 0
    assert first.summary()["mean_specific_rate_mmol_per_g_struct_h"] > 0
    assert first.summary()["content_mg_per_g_total_dw"] < 1000.0


def test_future_control_changes_do_not_change_the_past(network, product):
    controls = np.array([[1.0, 1.0, 0.7, 0.2], [0.8, 0.8, 1.0, 0.3]])
    changed = controls.copy()
    changed[1] = [0.5, 0.5, 1.4, 0.1]
    environment = InSilicoEnvironment("glucose", hours=4.0)
    first = simulate_controlled_product(network, product, [0, 2, 4], controls, environment,
                                        grid_points=5, output_dt=0.1)
    second = simulate_controlled_product(network, product, [0, 2, 4], changed, environment,
                                         grid_points=5, output_dt=0.1)
    prefix = first.times_h <= 2.0

    np.testing.assert_array_equal(first.product_mmol_per_l[prefix], second.product_mmol_per_l[prefix])
    assert first.product_mmol_per_l[-1] != pytest.approx(second.product_mmol_per_l[-1])
    from ystwin.generator.in_silico import observe_reporters

    latent = np.zeros((len(first.times_h), 3))
    observations = [observe_reporters(
        result.times_h, latent, result.growth_rate_per_h, result.structural_biomass_g_per_l,
        noise_cv=0.0,
    ) for result in (first, second)]
    # Reporter growth drivers are linearly interpolated. A new interval's rate
    # must not be interpolated backwards into observations used before its knot.
    np.testing.assert_allclose(observations[0][prefix], observations[1][prefix],
                               rtol=2e-7, atol=2e-8)


def test_control_knots_preserve_inventories_and_the_sampled_reporter_growth_convention(network, product):
    environment = InSilicoEnvironment("glucose", hours=2.0)
    controls = [[1.0, 1.0, 0.7, 0.2], [0.2, 0.5, 1.4, 0.1]]
    joined = simulate_controlled_product(network, product, [0.0, 1.0, 2.0], controls,
                                        environment, grid_points=5, output_dt=0.1)
    knot, = np.flatnonzero(joined.times_h == 1.0)
    continued = simulate_controlled_product(
        network, product, [0.0, 1.0], controls[1:],
        replace(environment, hours=1.0,
                initial_structural_biomass_g_per_l=joined.structural_biomass_g_per_l[knot],
                substrate_g_per_l=joined.substrate_mmol_per_l[knot]
                * joined.substrate_molar_mass_g_per_mol / 1000.0),
        initial_pool_mmol_per_l=joined.product_mmol_per_l[knot],
        grid_points=5, output_dt=0.1,
    )

    for name in ("structural_biomass_g_per_l", "substrate_mmol_per_l", "product_mmol_per_l"):
        np.testing.assert_allclose(getattr(joined, name)[knot:], getattr(continued, name),
                                   rtol=1e-12, atol=1e-12)
    previous = simulate_controlled_product(network, product, [0.0, 1.0], controls[:1],
                                           replace(environment, hours=1.0), grid_points=5, output_dt=0.1)
    np.testing.assert_array_equal(joined.growth_rate_per_h[:knot + 1], previous.growth_rate_per_h)
    np.testing.assert_allclose(joined.growth_rate_per_h[knot + 1:], continued.growth_rate_per_h[1:],
                               rtol=1e-12, atol=1e-12)


def test_intracellular_product_is_counted_once_in_total_dry_mass(network, product):
    result = simulate_controlled_product(network, product, [0, 2], [[1, 1, 0.7, 0.2]],
                                         InSilicoEnvironment("glucose", hours=2),
                                         grid_points=5, output_dt=0.1)
    expected = result.structural_biomass_g_per_l + result.product_mmol_per_l * product.molar_mass_g_per_mol / 1000

    np.testing.assert_allclose(result.total_dry_biomass_g_per_l, expected)
    assert result.product_mmol_per_l[0] == 0.0


def test_product_in_medium_is_not_added_to_cell_dry_mass(network, product):
    network.metabolites.product_c.compartment = "e"
    external = replace(product, output_basis="synthetic secreted product")
    result = simulate_controlled_product(network, external, [0, 2], [[1, 1, 0.7, 0.2]],
                                         InSilicoEnvironment("glucose", hours=2),
                                         grid_points=5, output_dt=0.1)
    np.testing.assert_array_equal(result.total_dry_biomass_g_per_l, result.structural_biomass_g_per_l)


def test_control_updates_preserve_the_carbon_inventory(network, product):
    result = simulate_controlled_product(network, product, [0, 2, 4],
                                         [[1, 1, 0.7, 0.2], [0.6, 0.8, 1.0, 0.4]],
                                         InSilicoEnvironment("glucose", hours=4),
                                         grid_points=5, output_dt=0.1)
    consumed = result.initial_substrate_mmol_per_l - result.substrate_mmol_per_l[-1]
    biomass_carbon = (result.structural_biomass_g_per_l[-1] - result.structural_biomass_g_per_l[0]) * 6
    product_carbon = result.product_mmol_per_l[-1] * 12

    assert biomass_carbon + product_carbon <= consumed * 6 + 1e-7
    assert result.total_substrate_consumed_mmol_per_l == pytest.approx(consumed, abs=1e-7)


@pytest.mark.parametrize("values", [
    [-0.1, 1, 0.7, 0.2], [1.1, 1, 0.7, 0.2], [1, 0, 0.7, 0.2],
    [1, 1, -1, 0.2], [1, 1, 0.7, 1.1], [1, 1, float("nan"), 0.2],
])
def test_invalid_control_domains_are_refused(values):
    with pytest.raises(ValueError):
        MetabolicControl(*values)


def test_malformed_schedule_is_refused(network, product):
    with pytest.raises(ValueError, match="schedule"):
        simulate_controlled_product(network, product, [0, 1, 1], [[1, 1, 0.7, 0.2]] * 2,
                                     InSilicoEnvironment("glucose", hours=1))


@pytest.fixture
def loop_runner():
    spec = importlib.util.spec_from_file_location(
        "in_silico_runner", paths.REPO_ROOT / "scripts/run_in_silico_loop.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_learning_only_command_reads_no_calibration_or_product_tables(loop_runner, tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("in-silico learning accessed an external table or metabolic model")

    monkeypatch.setattr(pd, "read_csv", forbidden)
    monkeypatch.setattr(loop_runner, "load_model", forbidden)
    destination = tmp_path / "run"
    result = loop_runner.main([
        "--output-dir", str(destination), "--train-episodes", "16", "--test-episodes", "2",
        "--skip-products", "--no-plots", "--skip-nuisance-audit",
    ])
    assert result == 0
    protocol = json.loads((destination / "protocol.json").read_text())
    assert protocol["run_complete"]
    assert not protocol["wet_lab_calibration_used"]
    assert not protocol["biological_validation"]
    assert (destination / "student.json").is_file()
    assert (destination / "split_manifest.csv").is_file()


def test_product_forecast_receives_only_the_observation_prefix(
        loop_runner, network, product, monkeypatch):
    seen = []

    class Student:
        def infer(self, times, observations):
            return np.zeros((len(times), 3))

        def initialize(self, estimates):
            return SimpleNamespace(state=estimates[-1], accepted=True, diagnostics=lambda: {})

        def forecast(self, times, inputs, prefix_times, prefix_observations):
            assert prefix_times[-1] <= 1.0 + 1e-12
            assert len(prefix_observations) == len(prefix_times)
            assert times[-1] == 2.0
            seen.append(len(prefix_times))
            return np.zeros((len(times), 3))

        def predict_controls(self, latent):
            return np.tile([1.0, 1.0, 0.7, 0.2], (len(latent), 1))

    monkeypatch.setattr(loop_runner, "_MODEL", network)
    monkeypatch.setattr(loop_runner, "_STUDENT", Student())
    monkeypatch.setattr(loop_runner, "_OPTIONS", SimpleNamespace(
        hours=2.0, dt=0.1, seed=0, noise_cv=0.001, prefix_hours=1.0,
        control_intervals=2, grid_points=5, production_seed=9000, production_protocol="pulse",
        modes=["student_prefix_forecast", "oracle_state_learned_control", "zero_state_baseline"]))
    monkeypatch.setattr(loop_runner, "install_product", lambda model, name: (model, product))
    monkeypatch.setattr(loop_runner, "separate_storage_biomass", lambda model: (
        model, SimpleNamespace(structural_mass_fraction=1.0, metadata={})))
    result = loop_runner._production_job(("test", "glucose_aerobic"))

    assert seen
    assert len(result["rows"]) == 4
    assert all(row["status"] == "predicted" for row in result["rows"])
    assert result["provenance"]["product_training_labels"] is False


@pytest.mark.parametrize("rejected", [False, True])
def test_assessment_counts_rejected_forecasts_instead_of_dropping_nan(loop_runner, rejected):
    states = pd.DataFrame([
        {"split": "test_protocol", "forecast_rmse": 0.001, "persistence_rmse": 0.1,
         "forecast_rejected": False},
        {"split": "test_protocol", "forecast_rmse": np.nan if rejected else 0.001,
         "persistence_rmse": 0.1, "forecast_rejected": rejected},
        {"split": "test_family_shift", "forecast_rmse": 0.1, "persistence_rmse": 0.2,
         "forecast_rejected": False},
    ])
    row = {"product": "test", "environment": "test", "status": "predicted"}
    errors = {f"{name}_absolute_percent_error": 1.0 for name in (
        "titre_mg_per_l", "mean_specific_rate_mmol_per_g_total_dw_h", "yield_g_per_g")}
    products = pd.DataFrame([
        {**row, "mode": "teacher"}, {**row, **errors, "mode": "student_prefix_forecast"},
    ])
    nuisance = pd.DataFrame([
        {"noise_cv": 0.0, "growth_profile": "slow", "state_rmse": 0.01, "status": "predicted"},
    ])
    assessment = loop_runner.assess_teacher_recovery(states, products, nuisance, 1)
    assert assessment["teacher_recovery_passed"] is (not rejected)
    assert assessment["cross_equation_generalization_passed"] is False
    assert assessment["biological_validation"] is False
    assert "in_silico_passed" not in assessment
    assert assessment["reference_forecasts_rejected"] == int(rejected)


def test_control_knots_use_bounded_estimates_and_no_future_observation(loop_runner):
    class Student:
        def infer(self, times, observations):
            return observations[:, 1:4] - 1.0

        def forecast(self, times, inputs, prefix_times, observations):
            return np.column_stack([times / 2, np.zeros(len(times)),
                                    np.where(times <= prefix_times[-1], -0.01, 0.0)])

        def initialize(self, estimates):
            return SimpleNamespace(state=np.maximum(estimates[-1], 0), accepted=True,
                                   diagnostics=lambda: {"prefix_projection_l2": 0.01})

        def predict_controls(self, latent):
            assert np.all(latent >= 0)
            return np.column_stack([np.ones(len(latent)), np.ones(len(latent)),
                                    0.7 + latent[:, 0] + latent[:, 2], np.full(len(latent), 0.2)])

    times = np.array([0, 0.5, 1.0, 1.5, 2.0])
    prefix_times = times[:3]
    observations = np.column_stack([np.ones(3), 1 + prefix_times / 2, np.ones(3),
                                    np.full(3, 0.99), np.ones(3)])
    _, controls, initializations = loop_runner._forecast_control_schedule(
        Student(), times, np.zeros((5, 2)), prefix_times, observations, np.array([0, 0.75, 1.0, 1.05, 1.5]))

    np.testing.assert_allclose(controls[:, 2], [0.7, 0.95, 1.2, 1.225, 1.45])
    assert initializations[1]["last_observation_h"] == 0.5
    assert len(initializations) == 3


def test_reassess_frozen_run_separates_teacher_recovery_without_refitting(loop_runner, tmp_path, monkeypatch):
    run = tmp_path / "frozen"
    run.mkdir()
    protocol = {"mode": "in_silico_teacher_student", "run_complete": True,
                "products": ["test"], "environments": ["test"], "biological_validation": False}
    original = json.dumps(protocol)
    (run / "protocol.json").write_text(original)
    pd.DataFrame([
        {"split": "test_protocol", "forecast_rmse": 0.001, "persistence_rmse": 0.1, "forecast_rejected": False},
        {"split": "test_family_shift", "forecast_rmse": 0.1, "persistence_rmse": 0.2, "forecast_rejected": False},
    ]).to_csv(run / "state_control_scores.csv", index=False)
    errors = {f"{name}_absolute_percent_error": 1.0 for name in (
        "titre_mg_per_l", "mean_specific_rate_mmol_per_g_total_dw_h", "yield_g_per_g")}
    pd.DataFrame([
        {"mode": "teacher", "product": "test", "environment": "test", "status": "predicted"},
        {"mode": "student_prefix_forecast", "product": "test", "environment": "test",
         "status": "predicted", **errors},
    ]).to_csv(run / "product_scores.csv", index=False)
    pd.DataFrame([{"noise_cv": 0.0, "growth_profile": "test", "state_rmse": 0.01,
                   "status": "predicted"}]).to_csv(run / "observation_invariance.csv", index=False)

    def forbidden(*args, **kwargs):
        raise AssertionError("review attempted to regenerate or refit the reference")

    monkeypatch.setattr(loop_runner, "fit_student", forbidden)
    monkeypatch.setattr(loop_runner, "generate_episodes", forbidden)
    output = tmp_path / "review.json"
    assert loop_runner.main(["--assess-run", str(run), "--assessment-output", str(output)]) == 0
    assessed = json.loads(output.read_text())
    assert assessed["teacher_recovery_passed"] is True
    assert assessed["cross_equation_generalization_passed"] is False
    assert assessed["biological_validation"] is False
    assert (run / "protocol.json").read_text() == original
    with pytest.raises(ValueError, match="overwrite"):
        loop_runner.reassess_run(run, output)


@pytest.mark.parametrize("allocation", [0.0, 0.2])
def test_reporting_separates_initial_inventory_from_new_production_and_names_its_basis(network, product, allocation):
    result = simulate_controlled_product(
        network, product, [0.0, 1.0], [[1.0, 1.0, 0.7, allocation]],
        InSilicoEnvironment("inventory", hours=1.0), initial_pool_mmol_per_l=0.01,
        grid_points=5, output_dt=0.1,
    )
    row = result.summary()
    net = result.product_mmol_per_l[-1] - result.product_mmol_per_l[0]
    assert result.product_mmol_per_l[0] == 0.01
    if allocation == 0.0:
        np.testing.assert_array_equal(result.product_mmol_per_l, 0.01)
        assert row["productivity_mg_per_l_h"] == row["yield_g_per_g"] == 0.0
    product_mass = product.molar_mass_g_per_mol  # g/mol is numerically mg/mmol.
    consumed_mg_l = result.total_substrate_consumed_mmol_per_l * result.substrate_molar_mass_g_per_mol
    exposure = result.structural_biomass_exposure_g_h_per_l
    total_exposure = exposure + np.trapezoid(result.product_mmol_per_l, result.times_h) * product_mass / 1000
    assert row["titre_mg_per_l"] == pytest.approx((0.01 + net) * product_mass)
    assert row["productivity_mg_per_l_h"] == pytest.approx(net * product_mass / result.times_h[-1])
    assert row["yield_g_per_g"] == pytest.approx(net * product_mass / consumed_mg_l)
    assert row["mean_specific_rate_mmol_per_g_struct_h"] == pytest.approx(net / exposure)
    assert row["mean_specific_rate_mmol_per_g_total_dw_h"] == pytest.approx(net / total_exposure)
    assert row["content_mg_per_g_total_dw"] == pytest.approx(
        row["titre_mg_per_l"] / result.total_dry_biomass_g_per_l[-1])
    assert row["simulation_route"] == "in_silico_gsmm_control_comparison"
    assert row["biological_validation"] is False
    assert row["prediction_basis"] == product.output_basis
    assert row["product_molar_mass_g_per_mol"] == product_mass
    assert row["substrate_molar_mass_g_per_mol"] == result.substrate_molar_mass_g_per_mol
    assert "initial product inventory excluded" in row["yield_basis"]
    assert "trapezoidal" in row["total_dry_mass_exposure_basis"]
    assert "preceding interval at control knots" in row["growth_rate_sampling"]
    assert result.stage_provenance[0]["integrator_route"] == "fixed_volume_comparison"
    assert result.stage_provenance[0]["specific_rate_units"] == "mmol per structural g per hour; growth per hour"


def test_actual_generator_student_control_and_network_interfaces_form_a_synthetic_comparison(
        loop_runner, network, product, monkeypatch):
    from ystwin.analysis import in_silico as learning
    from ystwin.fba import allocation
    from ystwin.generator import in_silico as generator
    from ystwin.mech import engine

    training = generator.generate_episodes(16, seed=417, hours=2.0, dt=0.1)
    student = learning.fit_student(training, seed=19)
    teacher, = generator.generate_episodes(1, seed=941, protocol="pulse", hours=2.0, dt=0.1, noise_cv=0.0)
    assert teacher.episode_id not in student.training_episode_ids
    assert tuple(field.name for field in fields(MetabolicControl)) == generator.CONTROL_NAMES
    assert student.metadata["targets"]["controller"] == list(generator.CONTROL_NAMES)
    assert student.metadata["biological_validation"] is False
    assert student.metadata["uses_product_targets"] is False
    metadata = json.dumps(student.metadata, sort_keys=True)
    environment = InSilicoEnvironment("synthetic interface test", hours=2.0, carbon_uptake=6.0)
    knots = np.array([0.0, 1.0, 2.0])
    teacher_controls = loop_runner._interpolate(knots[:-1], teacher.times_h, teacher.controls)
    reference = simulate_controlled_product(network, product, knots, teacher_controls,
                                            environment, grid_points=5, output_dt=0.1)
    observations = generator.observe_reporters(
        reference.times_h, loop_runner._interpolate(reference.times_h, teacher.times_h, teacher.latent),
        reference.growth_rate_per_h, reference.structural_biomass_g_per_l, noise_cv=0.0,
    )
    prefix = reference.times_h <= 1.0
    # The evaluator creates observations; only their prefix and declared imposed
    # inputs reach the fitted student. No latent, growth or control truth is passed.
    observations[~prefix] = np.nan

    def forbidden(*args, **kwargs):
        raise AssertionError("the fitted comparison must not refit, read targets, call a teacher, or enter the physical/allocation route")

    monkeypatch.setattr(learning, "fit_student", forbidden)
    monkeypatch.setattr(pd, "read_csv", forbidden)
    monkeypatch.setattr(engine, "simulate", forbidden)
    monkeypatch.setattr(allocation, "solve_allocation", forbidden)
    for name in ("generate_episodes", "simulate_teacher", "observe_reporters"):
        monkeypatch.setattr(loop_runner, name, forbidden)
        monkeypatch.setattr(generator, name, forbidden)
    forecast, controls, initializations = loop_runner._forecast_control_schedule(
        student, teacher.times_h, teacher.inputs, reference.times_h[prefix], observations[prefix], knots[:-1],
    )
    predicted = simulate_controlled_product(network, product, knots, controls,
                                            environment, grid_points=5, output_dt=0.1)
    assert np.isfinite(forecast).all()
    assert not np.array_equal(controls, teacher_controls)  # A fitted estimate, not copied labels.
    assert all(row["last_observation_h"] <= row["control_time_h"] for row in initializations)
    for control, stage in zip(controls, predicted.stage_provenance, strict=True):
        assert list(stage["control"]) == list(generator.CONTROL_NAMES)
        np.testing.assert_array_equal(list(stage["control"].values()), control)
        assert stage["rate_model"]["growth_retention"] == control[0]
        assert stage["ngam_applied"] == control[2]
        assert stage["rate_model"]["allocation_fraction"] == control[3]
    assert predicted.summary()["titre_mg_per_l"] > 0.0
    assert predicted.summary()["yield_g_per_g"] > 0.0
    assert predicted.summary()["biological_validation"] is False
    assert json.dumps(student.metadata, sort_keys=True) == metadata


def test_hybrid_comparison_reports_the_actual_coupling_and_keeps_its_dry_mass_basis(network, product):
    from ystwin import hybrid
    from ystwin.analysis.hybrid_stress import fit_sensor_model

    readings = pd.read_csv(paths.REPO_ROOT / "outputs/sensor_characterisation.csv")
    plates = tuple(sorted(readings.plate.unique()))
    sensor = fit_sensor_model(readings, training_plates=plates[:-1], seed=0)
    dose = sensor.metadata["calibrated_ranges"]["DTT"]["dose_mM"][1] / 2
    environment = hybrid.Environment("sensor comparison", stressor="DTT", dose_mM=dose, hours=1.0)
    for couple in (True, False):
        result = hybrid.simulate_hybrid_comparison(sensor, network, product, environment,
                                                   couple_stress=couple, grid_points=5, steps=10)
        row = result.summary()
        assert row["simulation_route"] == "gsmm_allocation_comparison"
        assert result.trajectory.simulation_route == "fixed_volume_comparison"
        assert row["biological_validation"] is False
        assert row["biomass_basis"] == "supplied network dry-mass basis; extra product is not added to biomass"
        assert row["applied_growth_retention"] == (result.stress.growth_retention if couple else 1.0)
        assert ("empirical growth-cap factor" in row["assumptions"]) is couple
        assert ("growth retention is not applied" in row["assumptions"]) is (not couple)
        assert row["titre_mg_per_l"] == result.trajectory.titre_mg_per_l
        assert row["mean_specific_rate_mmol_per_gdcw_h"] == result.trajectory.mean_specific_rate_mmol_per_gdcw_h
        assert row["yield_g_per_g"] == result.trajectory.yield_g_per_g
