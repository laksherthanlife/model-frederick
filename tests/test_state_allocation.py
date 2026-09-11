from __future__ import annotations

from dataclasses import replace
import importlib

import cobra
import numpy as np
import pytest

from ystwin.analysis.partial_orders import Source


SOURCE = Source("synthetic", "analytic allocation fixture", "constructed, not measured biology")


@pytest.fixture
def allocation():
    return importlib.import_module("ystwin.fba.allocation")


@pytest.fixture
def network():
    model = cobra.Model("allocation_toy")
    metabolites = {name: cobra.Metabolite(name, compartment=compartment)
                   for name, compartment in (("g_e", "e"), ("g_c", "c"),
                                              ("o_e", "e"), ("o_c", "c"), ("p_c", "c"))}
    definitions = [
        ("EX_G", {"g_e": -2.0}, (-1000.0, 1000.0)),
        ("EX_O", {"o_e": 1.0}, (0.0, 1000.0)),
        ("TG", {"g_e": -1.0, "g_c": 1.0}, (0.0, 1000.0)),
        ("TO", {"o_e": -1.0, "o_c": 1.0}, (0.0, 1000.0)),
        ("BIO", {"g_c": -1.0, "o_c": -1.0}, (0.0, 1000.0)),
        ("MAKE_P", {"g_c": -1.0, "p_c": 1.0}, (0.0, 1000.0)),
        ("PRODUCT", {"p_c": -1.0}, (0.0, 1000.0)),
    ]
    for name, stoichiometry, bounds in definitions:
        reaction = cobra.Reaction(name, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({metabolites[key]: value for key, value in stoichiometry.items()})
        model.add_reactions([reaction])
    model.objective = "BIO"
    return model


def identity(api, model):
    return api.ModelIdentity(model.id, "synthetic organism", "toy strain", "toy-v1", SOURCE)


def environment(api, oxygen=10.0, growth=2.0, condition="toy-condition"):
    return api.Environment(condition, "synthetic organism", "toy strain",
                           {"oxygen": oxygen, "growth": growth},
                           {"oxygen": "mmol/L", "growth": "1/h"}, SOURCE)


def state(api, name="growth", product_cap=8.0, mass=None):
    return api.StatePolicy(
        name=name, source=SOURCE, biomass_reaction="BIO",
        uptakes=(api.UptakePolicy("EX_G", 10.0, SOURCE),
                 api.UptakePolicy("EX_O", 10.0, SOURCE)),
        enzyme_capacities=(api.FluxCapacity("BIO", 2.0, 2.0, SOURCE),
                           api.FluxCapacity("MAKE_P", 0.0, product_cap, SOURCE)),
        tasks=(api.ObjectiveTask("PRODUCT", 1.0, SOURCE),),
        dry_mass_g_per_cell=mass, dry_mass_source=SOURCE if mass is not None else None,
    )


def test_state_solves_apply_uptake_enzyme_and_task_policy_without_mutating_model(allocation, network):
    before = {r.id: r.bounds for r in network.reactions}
    s = state(allocation, product_cap=3.0)
    result = allocation.solve_allocation(network, identity(allocation, network),
                                         environment(allocation), allocation.MechanisticPolicy(s))
    assert result.status == "optimal"
    assert result.fluxes["BIO"] == pytest.approx(2.0)
    assert result.fluxes["PRODUCT"] == pytest.approx(3.0)
    assert result.fluxes["EX_G"] == pytest.approx(-2.5)
    assert result.mass_balance_residual < 1e-8
    assert {r.id: r.bounds for r in network.reactions} == before
    assert network.objective.direction == "max"
    assert result.policy == "mechanistic"


def test_fractional_tasks_keep_final_maximization_distinct_from_source_minimization(allocation, network):
    policy = replace(state(allocation), tasks=(
        allocation.ObjectiveTask("BIO", 1.0, SOURCE),
        allocation.ObjectiveTask("PRODUCT", 0.5, SOURCE),
    ))
    result = allocation.solve_allocation(network, identity(allocation, network),
                                         environment(allocation), allocation.MechanisticPolicy(policy))
    assert result.status == "optimal"
    assert result.fluxes["BIO"] == pytest.approx(2.0)
    assert result.fluxes["PRODUCT"] == pytest.approx(8.0)
    with network:
        network.reactions.EX_G.lower_bound = -5.0
        network.reactions.EX_O.upper_bound = 10.0
        network.reactions.BIO.bounds = (2.0, 2.0)
        network.reactions.MAKE_P.upper_bound = 8.0
        for task in policy.tasks:
            network.objective = task.reaction
            maximum = network.optimize(objective_sense="maximize")
            assert maximum.status == "optimal"
            reaction = network.reactions.get_by_id(task.reaction)
            reaction.lower_bound = max(reaction.lower_bound, task.fraction * maximum.objective_value)
        minimum = network.optimize(objective_sense="minimize")
        assert minimum.status == "optimal"
        assert minimum.fluxes["BIO"] == pytest.approx(2.0)
        assert minimum.fluxes["PRODUCT"] == pytest.approx(4.0)


def test_state_specific_kinetics_respond_to_environment_at_fixed_growth(allocation, network):
    feature = allocation.Feature("oxygen", "mmol/L")
    s = replace(state(allocation), uptakes=(
        allocation.UptakePolicy("EX_G", 20.0, SOURCE, feature=feature, half_saturation=10.0),
        allocation.UptakePolicy("EX_O", 10.0, SOURCE),
    ))
    policy = allocation.MechanisticPolicy(s)
    low = allocation.solve_allocation(network, identity(allocation, network),
                                      environment(allocation, oxygen=2.0), policy)
    high = allocation.solve_allocation(network, identity(allocation, network),
                                       environment(allocation, oxygen=10.0), policy)
    assert low.growth_rate == high.growth_rate == pytest.approx(2.0)
    assert low.fluxes["PRODUCT"] == pytest.approx(20 / 6 - 2)
    assert high.fluxes["PRODUCT"] == pytest.approx(8.0)


def test_cell_fractions_are_converted_to_mass_weights_before_full_flux_blending(allocation, network):
    states = (state(allocation, "off", product_cap=0.0, mass=1e-12),
              state(allocation, "on", product_cap=8.0, mass=3e-12))
    distribution = allocation.FixedDistribution({"off": 0.5, "on": 0.5}, "cell", SOURCE)
    policy = allocation.MixturePolicy(states, distribution)
    result = allocation.solve_allocation(network, identity(allocation, network),
                                         environment(allocation), policy)
    assert result.status == "optimal"
    assert result.biomass_fractions == pytest.approx({"off": 0.25, "on": 0.75})
    assert result.cell_fractions == pytest.approx({"off": 0.5, "on": 0.5})
    assert result.fluxes["PRODUCT"] == pytest.approx(6.0)
    assert result.fluxes["EX_G"] == pytest.approx(-4.0)
    assert result.growth_rate == pytest.approx(2.0)
    assert set(result.fluxes) == {r.id for r in network.reactions}
    assert result.mass_balance_residual < 1e-8
    assert result.policy == "mixture"


def test_cell_fraction_blending_refuses_missing_cell_mass(allocation, network):
    policy = allocation.MixturePolicy((state(allocation),),
                                      allocation.FixedDistribution({"growth": 1.0}, "cell", SOURCE))
    with pytest.raises(ValueError, match="dry.mass"):
        allocation.solve_allocation(network, identity(allocation, network), environment(allocation), policy)


def test_infeasible_positive_weight_state_is_not_dropped_or_renormalized(allocation, network):
    infeasible = replace(state(allocation, "blocked"), uptakes=(
        allocation.UptakePolicy("EX_G", 0.0, SOURCE),
        allocation.UptakePolicy("EX_O", 10.0, SOURCE),
    ))
    policy = allocation.MixturePolicy((state(allocation, "ok"), infeasible),
                                      allocation.FixedDistribution({"ok": 0.9, "blocked": 0.1},
                                                                   "biomass", SOURCE))
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation), policy)
    assert result.status == "infeasible"
    assert result.growth_rate is None
    assert result.fluxes is None
    assert result.biomass_fractions["blocked"] == pytest.approx(0.1)


def test_zero_growth_is_a_feasible_state_not_an_infeasible_sentinel(allocation, network):
    s = replace(state(allocation), enzyme_capacities=(
        allocation.FluxCapacity("BIO", 0.0, 0.0, SOURCE),))
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                         allocation.MechanisticPolicy(s))
    assert result.status == "optimal"
    assert result.growth_rate == 0.0
    assert result.fluxes["PRODUCT"] == pytest.approx(10.0)


def training_rows(api):
    return [api.StateObservation(
        f"train-{i}", f"culture-{i}", "train", environment(api, oxygen=x),
        {"off": 1 - p, "on": p}, SOURCE,
    ) for i, (x, p) in enumerate(((0.0, 0.05), (1.0, 0.1), (9.0, 0.9), (10.0, 0.95)))]


def fit(api, rows):
    return api.fit_state_distribution(
        rows, features=(api.Feature("oxygen", "mmol/L"), api.Feature("growth", "1/h", "growth")),
        states=("off", "on"), projection_dimension=1, ridge=0.01, basis="biomass",
    )


def test_scaling_projection_and_distribution_are_training_only(allocation):
    rows = training_rows(allocation)
    test_row = allocation.StateObservation("held-out", "held-out-culture", "test",
                                          environment(allocation, oxygen=1e8),
                                          {"off": 1.0, "on": 0.0}, SOURCE)
    first = fit(allocation, rows)
    second = fit(allocation, [*rows, test_row])
    assert first.report() == second.report()
    np.testing.assert_allclose(first.center, [5.0, 2.0])
    low = first.predict(environment(allocation, oxygen=1.0))
    high = first.predict(environment(allocation, oxygen=9.0))
    assert high["on"] > low["on"] + 0.5
    assert sum(high.values()) == pytest.approx(1.0)
    assert set(first.training_ids) == {row.sample_id for row in rows}
    assert first.report()["biological_validation"] is False


def test_training_requires_environment_features_and_no_group_leakage(allocation):
    rows = training_rows(allocation)
    with pytest.raises(ValueError, match="environment"):
        allocation.fit_state_distribution(rows, features=(allocation.Feature("growth", "1/h", "growth"),),
                                          states=("off", "on"), projection_dimension=None, ridge=0.01)
    leaked = replace(rows[0], sample_id="test", split="test")
    with pytest.raises(ValueError, match="group"):
        fit(allocation, [*rows, leaked])


def test_prediction_refuses_wrong_environment_units_and_does_not_refit(allocation):
    fitted = fit(allocation, training_rows(allocation))
    snapshot = fitted.report()
    wrong = replace(environment(allocation), units={"oxygen": "percent", "growth": "1/h"})
    with pytest.raises(ValueError, match="unit"):
        fitted.predict(wrong)
    assert fitted.report() == snapshot


def test_host_and_strain_compatibility_are_explicit(allocation, network):
    policy = allocation.MechanisticPolicy(state(allocation))
    env = replace(environment(allocation), organism="CHO")
    with pytest.raises(ValueError, match="host"):
        allocation.solve_allocation(network, identity(allocation, network), env, policy)
    env = replace(environment(allocation), strain="another strain")
    with pytest.raises(ValueError, match="strain"):
        allocation.solve_allocation(network, identity(allocation, network), env, policy)
    result = allocation.solve_allocation(network, identity(allocation, network), env, policy,
                                         strain_transfer="synthetic cross-strain test, not a validation")
    assert result.status == "optimal"
    assert result.strain_transfer


def test_matched_comparison_refuses_different_base_networks(allocation, network):
    first = identity(allocation, network)
    second = replace(first, network="another-network")
    with pytest.raises(ValueError, match="network"):
        allocation.require_matched_models(first, second)


def test_fixed_objective_and_fixed_transition_are_separate_baselines(allocation, network):
    s = state(allocation)
    fixed = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                        allocation.FixedObjectivePolicy(s))
    assert fixed.policy == "fixed_objective"
    assert fixed.fluxes["PRODUCT"] == pytest.approx(0.0)
    transition = allocation.FixedTransition(("off", "on"), allocation.Feature("oxygen", "mmol/L"),
                                            midpoint=5.0, width=1.0, source=SOURCE)
    assert transition.predict(environment(allocation, oxygen=5.0)) == pytest.approx({"off": 0.5, "on": 0.5})


def test_protein_budget_is_never_automatically_calibrated_and_explicit_units_are_converted(allocation, network):
    pool = cobra.Metabolite("pool", compartment="c")
    enzyme = cobra.Metabolite("enzyme", compartment="c")
    supply = cobra.Reaction("protein_supply", lower_bound=-100.0, upper_bound=0.0)
    supply.add_metabolites({pool: -2.0})
    draw = cobra.Reaction("draw_enzyme", lower_bound=0.0, upper_bound=1000.0)
    draw.add_metabolites({pool: -50.0, enzyme: 1.0})
    network.add_reactions([supply, draw])
    network.reactions.MAKE_P.add_metabolites({enzyme: -1.0})
    original = supply.bounds
    native = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                          allocation.MechanisticPolicy(state(allocation)))
    explicit = replace(state(allocation), protein_budget=allocation.ProteinBudget(
        "protein_supply", 0.05, "g/gDW", "mg/gDW", SOURCE))
    constrained = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                               allocation.MechanisticPolicy(explicit))
    assert native.fluxes["PRODUCT"] == pytest.approx(4.0)
    assert constrained.fluxes["PRODUCT"] == pytest.approx(1.0)
    assert supply.bounds == original


def test_reactor_rhs_uses_dry_mass_and_physical_exchange_stoichiometry(allocation, network):
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                         allocation.MechanisticPolicy(state(allocation, product_cap=3.0)))
    rates = allocation.reactor_derivative(network, result, biomass_gdw_per_l=2.0,
                                          concentrations_mmol_per_l={"glucose": 5.0},
                                          exchanges={"glucose": "EX_G"}, dilution_per_h=0.1,
                                          feed_mmol_per_l={"glucose": 10.0})
    assert rates["glucose"] == pytest.approx(-9.5)
    assert rates["biomass_gDW_per_L_per_h"] == pytest.approx(3.8)


@pytest.mark.integration
def test_real_models_obey_measured_stoichiometric_uptake_without_pool_fitting(allocation, yeast_gem, ec_yeast_gem):
    from ystwin import paths
    from ystwin.bridge.regulation import exchange_coefficients
    from ystwin.fba.physiology import REFERENCE_AEROBIC_BATCH
    from ystwin.fba.solver import load_model

    reference = REFERENCE_AEROBIC_BATCH
    models = [yeast_gem, ec_yeast_gem]
    matched_path = paths.data_dir() / "gem" / "ecYeastGEM_yeast902.xml.gz"
    if matched_path.exists():
        models.append(load_model(matched_path)[0])
    for model in models:
        before = {reaction.id: reaction.bounds for reaction in model.reactions}
        declaration = allocation.ModelIdentity(model.id, "Saccharomyces cerevisiae", "model background",
                                                "declared separately, not a matched comparison", SOURCE)
        env = allocation.Environment("van-hoek-D0.4", declaration.organism, "DS28911", {}, {}, SOURCE)
        policy = allocation.MechanisticPolicy(allocation.StatePolicy(
            name="measured-medium", source=SOURCE,
            uptakes=(allocation.UptakePolicy("r_1714", reference.glucose_uptake, SOURCE, mode="fixed"),
                     allocation.UptakePolicy("r_1992", reference.oxygen_uptake, SOURCE, mode="fixed")),
        ))
        result = allocation.solve_allocation(model, declaration, env, policy,
                                             strain_transfer="engineering constraint check, not strain validation")
        if result.status != "optimal":
            assert model.id == "M_ecYeastGEM_batch_v8__46__3__46__4"
            assert result.status == "numerical_failure"
            assert result.fluxes is None and result.growth_rate is None
            assert result.states[0].mass_balance_residual > 10 * model.tolerance
        else:
            for exchange, expected in (("r_1714", reference.glucose_uptake), ("r_1992", reference.oxygen_uptake)):
                observed = sum(coefficient * result.fluxes[rid]
                               for rid, coefficient in exchange_coefficients(model, exchange).items())
                assert observed == pytest.approx(expected, abs=1e-6)
        assert {reaction.id: reaction.bounds for reaction in model.reactions} == before
        assert result.provenance["state_policies"][0]["protein_budget"] is None
        assert len(result.states[0].uptake_constraints) == 2


def test_allocation_cli_requires_explicit_outputs_and_labels_the_demo_synthetic(tmp_path):
    import json
    from pathlib import Path
    import subprocess
    import sys

    script = Path(__file__).resolve().parents[1] / "scripts" / "run_state_allocation.py"
    missing = subprocess.run([sys.executable, str(script), "--synthetic-demo"], capture_output=True, text=True)
    assert missing.returncode != 0
    output = tmp_path / "explicit-output"
    run = subprocess.run([sys.executable, str(script), "--synthetic-demo", "--output-dir", str(output)],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    data = json.loads((output / "state_allocation.json").read_text())
    assert data["biological_validation"] is False
    assert data["data_kind"] == "synthetic"
    assert {row["allocation"]["policy"] for row in data["results"]} == {
        "mixture", "mechanistic", "fixed_objective", "fixed_transition"}
    fitted = json.loads((output / "state_distribution.json").read_text())
    assert all(name.startswith("train-") for name in fitted["training_ids"])
    repeated = subprocess.run([sys.executable, str(script), "--synthetic-demo", "--output-dir", str(output)],
                              capture_output=True, text=True)
    assert repeated.returncode != 0


def test_config_driven_runner_loads_explicit_policies_and_training_rows(allocation, network, tmp_path):
    import json
    from pathlib import Path
    import runpy

    model_path = tmp_path / "model.xml"
    cobra.io.write_sbml_model(network, model_path)
    config = {
        "schema_version": 1, "model": "model.xml", "identity": allocation._report(identity(allocation, network)),
        "features": [{"name": "oxygen", "unit": "mmol/L"}, {"name": "growth", "unit": "1/h", "role": "growth"}],
        "states": [allocation._report(state(allocation, "off", product_cap=0.0)),
                   allocation._report(state(allocation, "on", product_cap=8.0))],
        "observations": [allocation._report(row) for row in training_rows(allocation)],
        "fit": {"projection_dimension": 1, "ridge": 0.01, "basis": "biomass"},
        "baselines": {"mechanistic_state": "on", "fixed_objective_state": "on"},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    script = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts" / "run_state_allocation.py"))
    output = tmp_path / "config-output"
    script["main"](["--config", str(path), "--output-dir", str(output)])
    result = json.loads((output / "state_allocation.json").read_text())
    assert result["config_sha256"]
    assert result["data_kind"] == "synthetic"
    assert all(row["allocation"]["status"] == "optimal" for row in result["results"])
    assert result["model_path"] == str(model_path)


def test_full_flux_biomass_growth_and_cell_growth_use_their_own_consistent_weights(allocation, network):
    off = state(allocation, "off", product_cap=0.0, mass=1e-12)
    on = replace(state(allocation, "on", product_cap=6.0, mass=3e-12), enzyme_capacities=(
        allocation.FluxCapacity("BIO", 4.0, 4.0, SOURCE),
        allocation.FluxCapacity("MAKE_P", 0.0, 6.0, SOURCE)))
    policy = allocation.MixturePolicy((off, on), allocation.FixedDistribution({"off": 0.5, "on": 0.5}, "cell", SOURCE))
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation), policy)
    assert result.growth_rate == pytest.approx(3.5)
    assert result.cell_population_growth_rate == pytest.approx(3.0)
    assert result.fluxes["PRODUCT"] == pytest.approx(4.5)
    assert result.fluxes["EX_G"] == pytest.approx(-4.0)
    assert result.mean_cell_dry_mass_g == pytest.approx(2e-12, abs=1e-20)
    assert result.mass_balance_residual < 1e-8


def test_allocation_report_is_strict_json_with_unbounded_native_exchange_capacity(allocation, network):
    import json

    network.reactions.EX_G.upper_bound = float("inf")
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation),
                                         allocation.MechanisticPolicy(state(allocation)))
    assert result.status == "optimal"
    json.dumps(result.report(), allow_nan=False)


def test_source_audit_records_verified_code_and_missing_tables():
    import json
    from pathlib import Path

    document = json.loads((Path(__file__).resolve().parents[1] / "data" / "cosmic_sources.json").read_text())
    implementation = document["public_implementation"]
    assert implementation["commit"] == "3d60d3a2ae67943d770d7f8fbbb8bcc16f391a86"
    assert implementation["sha256"] == "b82441d85b96c044e82a94f3bc27ab9e1c339ba9e5db0183133a8e7d21dd9b10"
    assert document["supplementary_and_training_data"]["missing"]
    assert document["local_implementation_boundary"]["biological_validation"] is False


def test_mechanistic_cap_is_not_applied_again_after_blending(allocation, network):
    states = (state(allocation, "a", product_cap=4.0), state(allocation, "b", product_cap=4.0))
    policy = allocation.MixturePolicy(states, allocation.FixedDistribution({"a": 0.5, "b": 0.5},
                                                                           "biomass", SOURCE))
    result = allocation.solve_allocation(network, identity(allocation, network), environment(allocation), policy)
    assert result.fluxes["PRODUCT"] == pytest.approx(4.0)
