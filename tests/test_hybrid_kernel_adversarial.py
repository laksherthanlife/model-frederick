from __future__ import annotations

import math

import cobra
import numpy as np
import pytest
from cobra.exceptions import OptimizationError
from swiglpk import GLP_DUAL, GLP_PRIMAL

from ystwin.fba import dynamic_rates
from ystwin.fba.dynamic import FeedProfile, simulate_batch
from ystwin.fba.dynamic_rates import (
    MAX_STRESS_NGAM,
    RESTING_NGAM,
    NetworkRateClosure,
    build_network_rates,
    build_rates,
)
from ystwin.fba.physiology import cap_uptake
from ystwin.fba.solver import PINNED_TOLERANCE


def _reaction(model, reaction_id, metabolites, bounds=(0.0, 1000.0)):
    reaction = cobra.Reaction(reaction_id, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites(metabolites)
    model.add_reactions([reaction])
    return reaction


def _exchange(model, reaction_id, metabolite, convention, coefficient=1.0):
    if convention == "split":
        _reaction(model, reaction_id, {metabolite: -coefficient}, (0.0, 0.0))
        return _reaction(model, f"{reaction_id}_REV", {metabolite: coefficient},
                         (0.0, 1000.0 / coefficient))
    sign = 1.0 if convention == "positive" else -1.0
    return _reaction(model, reaction_id, {metabolite: sign * coefficient},
                     (-1000.0 / coefficient, 1000.0 / coefficient))


def _network(convention="negative", coefficient=1.0, growth_atp=0.0, product_atp=0.0):
    model = cobra.Model("hybrid_adversarial")
    carbon_e = cobra.Metabolite("carbon_e", formula="C", compartment="e")
    carbon_c = cobra.Metabolite("carbon_c", formula="C", compartment="c")
    product = cobra.Metabolite("product_c", formula="C2", compartment="c")
    oxygen = cobra.Metabolite("oxygen_e", formula="O2", compartment="e")
    co2 = cobra.Metabolite("co2_e", formula="CO2", compartment="e")
    atp = cobra.Metabolite("atp_c", compartment="c")
    _exchange(model, "EX_carbon", carbon_e, convention, coefficient)
    _exchange(model, "EX_oxygen", oxygen, convention)
    _reaction(model, "transport", {carbon_e: -1.0, carbon_c: 1.0})
    _reaction(model, "growth", {carbon_c: -10.0, atp: -growth_atp})
    _reaction(model, "synthesis", {carbon_c: -2.0, product: 1.0, atp: -product_atp})
    _reaction(model, "DM_product", {product: -1.0})
    _reaction(model, "respiration", {carbon_c: -1.0, oxygen: -1.0, atp: 1.0, co2: 1.0})
    _reaction(model, "ATPM", {atp: -1.0})
    _reaction(model, "EX_co2", {co2: -1.0}, (-1000.0, 1000.0))
    model.objective = "growth"
    return model


def _build(model=None, **kwargs):
    options = {
        "product_reaction_id": "DM_product",
        "biomass_reaction_id": "growth",
        "growth_fraction": 0.8,
        "grid_points": 11,
        "kcat_per_s": 1.0,
        "enzyme_mmol_per_gdcw": 1.0,
    }
    options.update(kwargs)
    return build_rates(_network() if model is None else model,
                       "EX_carbon", "EX_oxygen", "ATPM", **options)


def _substrate_for_uptake(allowed, vmax=10.0, km=0.5):
    return km * allowed / (vmax - allowed)


def _uptake_reaction(model):
    if "EX_carbon_REV" in model.reactions:
        return model.reactions.get_by_id("EX_carbon_REV")
    return model.reactions.get_by_id("EX_carbon")


def _assert_flux_feasible(model, observed, ngam=RESTING_NGAM):
    scoped = model.copy()
    mu, uptake, product = observed
    boundary = _uptake_reaction(scoped)
    coefficient = next(iter(boundary.metabolites.values()))
    boundary.bounds = (uptake / coefficient, uptake / coefficient)
    scoped.reactions.growth.bounds = (mu, mu)
    scoped.reactions.DM_product.bounds = (product, product)
    scoped.reactions.ATPM.lower_bound = ngam
    scoped.reactions.EX_co2.lower_bound = 0.0
    scoped.objective = 0
    value = scoped.slim_optimize()
    assert scoped.solver.status == "optimal", (observed, ngam, scoped.solver.status)
    assert np.isfinite(value)
    fluxes = np.array([reaction.flux for reaction in scoped.reactions])
    stoichiometry = cobra.util.array.create_stoichiometric_matrix(scoped)
    np.testing.assert_allclose(stoichiometry @ fluxes, 0.0, atol=1e-7, rtol=0.0)
    for reaction, flux in zip(scoped.reactions, fluxes):
        assert reaction.lower_bound - 1e-7 <= flux <= reaction.upper_bound + 1e-7


def _snapshot(model):
    return (
        tuple((r.id, r.bounds, tuple((m.id, c) for m, c in r.metabolites.items()))
              for r in model.reactions),
        str(model.objective.expression),
        model.objective.direction,
        model.tolerance,
        type(model.solver),
        tuple((constraint.name, constraint.lb, constraint.ub, str(constraint.expression))
              for constraint in model.constraints),
    )


@pytest.mark.parametrize("substrate", [0.1, 1.0, 100.0])
def test_split_ec_exchanges_preserve_rates_and_oxygen_limitation(substrate):
    plain = _network(growth_atp=1.0, product_atp=1.0)
    split = _network("split", growth_atp=1.0, product_atp=1.0)
    expected, _ = _build(plain, oxygen_lower_bound=-0.9)
    actual, _ = _build(split, oxygen_lower_bound=-0.9)
    assert actual(substrate) == pytest.approx(expected(substrate), abs=1e-7)
    mu, uptake, product = actual(substrate)
    assert 0.0 < uptake <= 10.0 * substrate / (0.5 + substrate) + 1e-7
    assert RESTING_NGAM + mu + product <= 0.9 + 1e-7
    assert uptake == pytest.approx(11.0 * mu + 3.0 * product + RESTING_NGAM, abs=1e-7)


def test_reorienting_a_one_metabolite_exchange_does_not_open_uncapped_supply():
    expected, _ = _build(_network("negative"))
    actual, _ = _build(_network("positive"))
    substrate = 0.5
    assert actual(substrate)[1] <= 5.0 + 1e-7
    assert actual(substrate) == pytest.approx(expected(substrate), abs=1e-7)


@pytest.mark.parametrize("convention", ["negative", "split"])
@pytest.mark.parametrize("coefficient", [0.5, 2.0])
def test_growth_retention_is_invariant_to_exchange_stoichiometry_units(convention, coefficient):
    expected, _ = _build(_network(convention), growth_retention=0.5)
    actual, _ = _build(_network(convention, coefficient), growth_retention=0.5)
    assert actual(10.0) == pytest.approx(expected(10.0), abs=1e-7)


@pytest.mark.parametrize("convention", ["negative", "positive", "split"])
def test_only_the_declared_carbon_feed_can_supply_the_batch(convention):
    model = _network()
    alternate = cobra.Metabolite("alternate_e", formula="C2", compartment="e")
    _exchange(model, "EX_alternate", alternate, convention)
    _reaction(model, "alternate_transport", {
        alternate: -1.0, model.metabolites.carbon_c: 2.0,
    })
    expected, _ = _build()
    actual, _ = _build(model)
    assert actual(10.0) == pytest.approx(expected(10.0), abs=1e-7)


def test_an_open_intracellular_carbon_sink_cannot_bypass_the_feed_ledger():
    model = _network()
    sink = _reaction(model, "SK_carbon", {model.metabolites.carbon_c: -1.0},
                     (-1000.0, 1000.0))
    sink.annotation["sbo"] = "SBO:0000632"
    assert sink in model.boundary
    assert sink not in model.exchanges
    expected, _ = _build()
    actual, _ = _build(model)
    assert actual(10.0) == pytest.approx(expected(10.0), abs=1e-7)


@pytest.mark.parametrize("stress", [0.0, 0.5, 1.0])
def test_a_feasible_near_maintenance_state_is_not_replaced_with_dormancy(stress):
    rates, info = _build(stress=stress)
    allowed = info.ngam + 0.01
    expected = (0.8 * (allowed - info.ngam) / 10.0,
                allowed, 0.2 * (allowed - info.ngam) / 2.0)
    _assert_flux_feasible(_network(), expected, info.ngam)
    observed = rates(_substrate_for_uptake(allowed))
    assert observed == pytest.approx(expected, abs=1e-7)


def test_the_first_feasible_grid_node_is_itself_callable():
    rates, info = _build()
    index = int(np.flatnonzero(rates.feasible)[0])
    assert rates.uptake_grid[index] < rates.uptake_grid[-1]
    observed = rates(_substrate_for_uptake(float(rates.uptake_grid[index])))
    expected = tuple(float(value) for value in rates.grid_fluxes[index])
    assert observed == pytest.approx(expected, abs=1e-7)
    _assert_flux_feasible(_network(), observed, info.ngam)


@pytest.mark.parametrize("allowed", [1.05, 1.6, 2.0, 3.75, 5.4, 8.9])
def test_interpolated_rates_have_a_complete_feasible_flux_completion(allowed):
    model = _network(growth_atp=0.5, product_atp=0.25)
    model.reactions.synthesis.upper_bound = 0.1
    _reaction(model, "alternate_synthesis", {
        model.metabolites.carbon_c: -3.0,
        model.metabolites.product_c: 1.0,
        model.metabolites.co2_e: 1.0,
    }, (0.0, 0.7))
    pool = model.problem.Constraint(
        model.reactions.growth.flux_expression
        + model.reactions.alternate_synthesis.flux_expression,
        ub=0.55, name="shared_capacity",
    )
    model.add_cons_vars([pool])
    rates, info = _build(model, enzyme_mmol_per_gdcw=0.0002)
    observed = rates(_substrate_for_uptake(allowed))
    assert observed[0] > 0.0
    assert observed[2] > 0.0
    assert observed[1] <= allowed + 1e-7
    assert observed[2] <= info.q_product_max_mmol_per_gdcw_h * allowed / 10.0 + 1e-7
    _assert_flux_feasible(model, observed, info.ngam)


@pytest.mark.parametrize("primary_blocked", [False, True])
def test_alternate_product_routes_cannot_bypass_the_terminal_kinetic_cap(primary_blocked):
    model = _network()
    if primary_blocked:
        model.reactions.synthesis.bounds = (0.0, 0.0)
    _reaction(model, "alternate_synthesis", {
        model.metabolites.carbon_c: -3.0,
        model.metabolites.product_c: 1.0,
        model.metabolites.co2_e: 1.0,
    })
    rates, info = _build(model, enzyme_mmol_per_gdcw=0.0001, growth_fraction=0.5)
    substrate = 2.0
    observed = rates(substrate)
    cap = info.q_product_max_mmol_per_gdcw_h * substrate / (0.5 + substrate)
    assert observed[2] == pytest.approx(cap, abs=1e-7)
    assert observed[1] < 10.0 * substrate / (0.5 + substrate) - 0.1
    _assert_flux_feasible(model, observed, info.ngam)


def test_a_carbon_limited_product_flux_does_not_scale_with_a_nonbinding_kcat():
    normal, _ = _build()
    inflated, _ = _build(kcat_per_s=100.0, regulation_scale=10.0)
    assert normal(10.0)[2] > 0.0
    assert inflated(10.0) == pytest.approx(normal(10.0), abs=1e-7)


def test_an_infeasible_maintenance_network_is_not_a_successful_zero_titre_prediction():
    model = _network()
    model.reactions.respiration.bounds = (0.0, 0.0)
    with pytest.raises((ValueError, RuntimeError, OptimizationError)):
        _build(model)


def test_a_preexisting_maintenance_capacity_is_not_silently_relaxed():
    model = _network()
    model.reactions.ATPM.upper_bound = RESTING_NGAM / 2.0
    with pytest.raises((ValueError, RuntimeError, OptimizationError)):
        _build(model)


def test_an_unbounded_growth_objective_is_not_reinterpreted_as_zero_growth():
    model = _network()
    model.reactions.growth.add_metabolites({model.metabolites.carbon_c: 10.0})
    model.reactions.growth.upper_bound = math.inf
    assert not model.reactions.growth.metabolites
    with pytest.raises((ValueError, RuntimeError, OptimizationError)):
        _build(model)


def test_a_solver_nan_without_an_infeasibility_certificate_is_not_biological_dormancy(monkeypatch):
    monkeypatch.setattr(cobra.Model, "slim_optimize", lambda self, *args, **kwargs: math.nan)
    with pytest.raises((ValueError, RuntimeError, OptimizationError)):
        _build()


@pytest.mark.parametrize("name", [
    "kcat_per_s", "enzyme_mmol_per_gdcw", "stress", "regulation_scale",
    "uptake_vmax_mmol_per_gdcw_h", "uptake_km_mmol_per_l", "oxygen_lower_bound",
    "growth_fraction", "growth_retention", "maximum_growth_per_h",
])
@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_rate_inputs_are_rejected_before_solving(name, value, monkeypatch):
    def unexpected_solve(*args, **kwargs):
        pytest.fail("a nonfinite input reached the metabolic solver")

    monkeypatch.setattr(dynamic_rates, "growth_or_none", unexpected_solve)
    with pytest.raises(ValueError):
        _build(**{name: value})


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_a_nonfinite_grid_size_has_a_validation_error(value):
    with pytest.raises(ValueError, match="grid_points"):
        _build(grid_points=value)


def test_finite_factors_cannot_overflow_the_product_capacity_silently():
    with pytest.raises(ValueError, match="finite|capacity|kcat"):
        _build(kcat_per_s=1e308)


@pytest.mark.parametrize("substrate", [math.nan, math.inf, -math.inf, -0.1])
def test_nonfinite_or_negative_substrate_cannot_reach_a_cached_rate(substrate):
    rates, _ = _build()
    with pytest.raises(ValueError, match="substrate"):
        rates(substrate)


@pytest.mark.parametrize("observed", [
    (math.nan, 1.0, 0.1), (0.1, math.inf, 0.1),
    (0.1, 1.0, math.nan), (0.1, 1.0, -math.inf),
])
def test_nonfinite_rates_are_not_integrated(observed):
    with pytest.raises(ValueError, match="rates"):
        simulate_batch(lambda substrate: observed, 1.0, 10.0, 1.0, 24.0, steps=2)


@pytest.mark.parametrize("kwargs", [
    {"initial_biomass_g_per_l": math.inf},
    {"initial_substrate_mmol_per_l": math.inf},
    {"hours": math.nan},
    {"product_molar_mass_g_per_mol": math.nan},
    {"substrate_molar_mass_g_per_mol": math.inf},
])
def test_nonfinite_batch_inputs_fail_before_requesting_rates(kwargs):
    def unexpected_rates(substrate):
        pytest.fail("a nonfinite batch input reached the rate closure")

    options = dict(initial_biomass_g_per_l=1.0, initial_substrate_mmol_per_l=10.0,
                   hours=1.0, product_molar_mass_g_per_mol=24.0, steps=2)
    options.update(kwargs)
    with pytest.raises(ValueError):
        simulate_batch(unexpected_rates, **options)


@pytest.mark.parametrize("kwargs", [
    {"rate_mmol_per_l_h": math.inf},
    {"start_h": math.inf},
    {"end_h": math.nan},
    {"end_h": -math.inf},
])
def test_nonfinite_feed_parameters_do_not_create_unreported_carbon(kwargs):
    with pytest.raises(ValueError):
        FeedProfile(**kwargs)


def test_source_model_constraints_solver_settings_and_objective_are_isolated():
    model = _network("split")
    model.tolerance = 1e-6
    model.objective = model.problem.Objective(
        2.0 * model.reactions.synthesis.flux_expression, direction="min")
    model.add_cons_vars([model.problem.Constraint(
        model.reactions.transport.flux_expression, ub=8.0, name="source_capacity")])
    before = _snapshot(model)
    first, _ = _build(model, stress=0.0)
    first_value = first(10.0)
    _build(model, stress=0.4, growth_retention=0.3, regulation_scale=0.5)
    assert _snapshot(model) == before
    model.reactions.transport.upper_bound = 0.0
    assert first(10.0) == first_value


def test_source_model_is_unchanged_when_a_scenario_fails_after_copying():
    model = _network("split")
    before = _snapshot(model)
    with pytest.raises(KeyError):
        _build(model, product_reaction_id="missing_product")
    assert _snapshot(model) == before


def _changing_ratio_rates(substrate):
    if substrate <= 0.0:
        return 0.0, 0.0, 0.0
    return 0.0, RESTING_NGAM + 0.4 * substrate, 0.2 * substrate


def _changing_ratio_exact(hours):
    substrate = (1.0 + RESTING_NGAM / 0.4) * math.exp(-0.4 * hours) - RESTING_NGAM / 0.4
    product = (1.0 - substrate - RESTING_NGAM * hours) / 2.0
    return substrate, product


def test_substrate_rate_scaling_cannot_discount_mandatory_atp_maintenance():
    hours = 1.1
    exact_substrate, _ = _changing_ratio_exact(hours)
    assert exact_substrate > 0.0
    for substrate in (1.0, exact_substrate):
        _assert_flux_feasible(_network(), _changing_ratio_rates(substrate))
    result = simulate_batch(_changing_ratio_rates, 1.0, 1.0, hours, 24.0, steps=2)
    maintenance_carbon = result.substrate_consumed_mmol_per_l - 2.0 * result.product_mmol_per_l[-1]
    assert maintenance_carbon >= RESTING_NGAM * hours - 1e-9


def test_changing_flux_ratios_converge_to_the_analytic_metabolic_trajectory():
    hours = 1.1
    exact_substrate, exact_product = _changing_ratio_exact(hours)
    initial = _changing_ratio_rates(1.0)
    final = _changing_ratio_rates(exact_substrate)
    assert initial[2] / initial[1] > 10.0 * final[2] / final[1]
    errors = []
    for steps in (100, 400, 1600):
        result = simulate_batch(_changing_ratio_rates, 1.0, 1.0, hours, 24.0, steps=steps)
        errors.append(abs(result.product_mmol_per_l[-1] - exact_product))
    assert errors[1] < 0.3 * errors[0]
    assert errors[2] < 0.3 * errors[1]
    assert errors[-1] < 0.001 * exact_product


def test_stress_maintenance_is_present_in_every_positive_grid_flux():
    rates, info = _build(stress=0.8)
    assert info.ngam == pytest.approx(RESTING_NGAM + 0.8 * MAX_STRESS_NGAM)
    for observed in rates.grid_fluxes[rates.feasible]:
        _assert_flux_feasible(_network(), observed, info.ngam)
        mu, uptake, product = observed
        assert uptake == pytest.approx(10.0 * mu + 2.0 * product + info.ngam, abs=1e-7)


def _network_rates(model=None, **kwargs):
    if model is None:
        model = _network()
        model.reactions.ATPM.lower_bound = RESTING_NGAM
    options = dict(product_reaction_id="DM_product", biomass_reaction_id="growth",
                   growth_fraction=0.8, grid_points=11)
    options.update(kwargs)
    return build_network_rates(model, "EX_carbon", "EX_oxygen", "ATPM", **options)


def _protein_budget_model():
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    protein = cobra.Metabolite("protein_budget_c", compartment="c")
    _reaction(model, "protein_pool", {protein: 1.0}, (0.0, 0.5))
    model.reactions.growth.add_metabolites({protein: -1.0})
    model.reactions.synthesis.add_metabolites({protein: -1.0})
    return model


@pytest.mark.parametrize("allocation", [0.0, 0.25, 1.0])
def test_network_allocation_shares_the_supplied_protein_budget(allocation):
    model = _protein_budget_model()
    before = _snapshot(model)
    rates, info = _network_rates(model, allocation_fraction=allocation)
    observed = rates.at_uptake(10.0)
    expected_product = allocation * 0.1
    expected_growth = 0.5 - expected_product
    expected_uptake = 10.0 * expected_growth + 2.0 * expected_product + RESTING_NGAM
    assert observed == pytest.approx((expected_growth, expected_uptake, expected_product), abs=1e-7)
    _assert_flux_feasible(model, observed)
    assert _snapshot(model) == before
    assert info.attainable_growth_per_h == pytest.approx(0.5)
    assert info.allocation_fraction == allocation
    assert isinstance(info, NetworkRateClosure)
    assert "allocation-conditioned" in info.product_rate_basis
    assert "not measured" in info.summary()
    assert not hasattr(info, "kcat_per_s")
    assert not hasattr(info, "enzyme_mmol_per_gdcw")


def test_network_product_interval_minimum_includes_growth_coupled_secretion():
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    model.reactions.growth.add_metabolites({
        model.metabolites.carbon_c: -2.0,
        model.metabolites.product_c: 1.0,
    })
    rates, info = _network_rates(model, allocation_fraction=0.25)
    attainable = (10.0 - RESTING_NGAM) / 12.0
    floor = 0.8 * attainable
    q_min = floor
    q_max = floor + (10.0 - RESTING_NGAM - 12.0 * floor) / 2.0
    product = q_min + 0.25 * (q_max - q_min)
    growth = (10.0 - RESTING_NGAM - 2.0 * product) / 10.0
    assert q_min > 0.0
    assert info.attainable_growth_per_h == pytest.approx(attainable)
    assert rates.at_uptake(10.0) == pytest.approx((growth, 10.0, product), abs=1e-7)
    _assert_flux_feasible(model, rates.at_uptake(10.0))


def test_network_product_interval_retains_a_required_positive_demand():
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    model.reactions.DM_product.lower_bound = 0.2
    rates, _ = _network_rates(model, allocation_fraction=0.25)
    floor = 0.8 * (10.0 - RESTING_NGAM - 2.0 * 0.2) / 10.0
    q_max = (10.0 - RESTING_NGAM - 10.0 * floor) / 2.0
    expected = 0.2 + 0.25 * (q_max - 0.2)
    assert rates.at_uptake(10.0)[2] == pytest.approx(expected, abs=1e-7)
    assert rates.minimum_feasible_uptake_mmol_per_gdcw_h == pytest.approx(1.1)
    _assert_flux_feasible(model, rates.at_uptake(10.0))


def test_fixed_measured_growth_is_exact_and_reports_actual_uptake():
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    model.reactions.DM_product.lower_bound = 0.3
    rates, info = _network_rates(model, fixed_growth_per_h=0.2,
                                growth_fraction=0.1, allocation_fraction=0.25)
    q_max = (10.0 - RESTING_NGAM - 10.0 * 0.2) / 2.0
    product = 0.3 + 0.25 * (q_max - 0.3)
    uptake = 10.0 * 0.2 + 2.0 * product + RESTING_NGAM
    assert rates.at_uptake(10.0) == pytest.approx((0.2, uptake, product), abs=1e-7)
    assert uptake < 10.0
    assert info.physiology_conditioned
    assert info.fixed_growth_per_h == 0.2
    for observed in rates.grid_fluxes[rates.feasible]:
        assert observed[0] == pytest.approx(0.2, abs=1e-7)
        _assert_flux_feasible(model, observed)
    with pytest.raises(ValueError, match="fixed_growth_per_h.*infeasible"):
        rates.at_uptake(2.0)
    with pytest.raises(ValueError, match="fixed_growth_per_h.*infeasible"):
        rates(0.0)


def test_infeasible_measured_growth_is_not_silently_replaced_with_attainable_growth():
    with pytest.raises(ValueError, match="fixed_growth_per_h.*infeasible"):
        _network_rates(fixed_growth_per_h=1.0)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf, -0.1])
def test_invalid_fixed_growth_is_rejected_before_a_network_solve(value):
    with pytest.raises(ValueError, match="fixed_growth_per_h"):
        _network_rates(fixed_growth_per_h=value)


@pytest.mark.parametrize("value", [math.nan, math.inf, -0.1, 1.1])
def test_invalid_allocation_is_rejected_before_a_network_solve(value):
    with pytest.raises(ValueError, match="allocation_fraction"):
        _network_rates(allocation_fraction=value)


@pytest.mark.parametrize("convention", ["negative", "positive", "split"])
@pytest.mark.parametrize("coefficient", [0.5, 2.0])
def test_network_allocation_is_invariant_to_substrate_exchange_units(convention, coefficient):
    model = _network(convention, coefficient)
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    expected, _ = _network_rates(growth_retention=0.5)
    actual, _ = _network_rates(model, growth_retention=0.5)
    assert actual.at_uptake(10.0) == pytest.approx(expected.at_uptake(10.0), abs=1e-7)


def test_network_allocation_preserves_the_supplied_maintenance_requirement():
    model = _network()
    model.reactions.ATPM.bounds = (1.3, 2.0)
    before = _snapshot(model)
    rates, info = _network_rates(model)
    assert info.ngam == 1.3
    _assert_flux_feasible(model, rates.at_uptake(10.0), ngam=1.3)
    assert _snapshot(model) == before


def test_scalar_stress_cannot_lower_a_stronger_existing_maintenance_requirement():
    model = _network()
    model.reactions.ATPM.bounds = (1.3, 2.0)
    rates, info = _build(model, stress=0.0)
    assert info.ngam == 1.3
    _assert_flux_feasible(model, rates.at_uptake(10.0), ngam=1.3)


@pytest.mark.parametrize("allowed", [0.71, 1.1, 3.7, 5.6, 7.4, 9.9])
def test_network_interpolation_completes_to_a_full_protein_constrained_flux(allowed):
    model = _protein_budget_model()
    rates, _ = _network_rates(model)
    observed = rates.at_uptake(allowed)
    assert observed[1] <= allowed + PINNED_TOLERANCE
    assert observed[0] + observed[2] <= 0.5 + PINNED_TOLERANCE
    _assert_flux_feasible(model, observed)


def test_network_at_uptake_and_monod_queries_agree_at_every_feasible_node():
    rates, _ = _network_rates()
    for allowed, observed, feasible in zip(rates.uptake_grid, rates.grid_fluxes, rates.feasible):
        if not feasible:
            continue
        expected = tuple(float(value) for value in observed)
        assert rates.at_uptake(float(allowed)) == pytest.approx(expected, abs=1e-7)
        if allowed < rates.uptake_grid[-1]:
            assert rates(_substrate_for_uptake(allowed)) == pytest.approx(expected, abs=1e-7)


@pytest.mark.parametrize("allowed", [math.nan, math.inf, -math.inf, -0.1, 10.1])
def test_direct_molecular_uptake_queries_must_be_inside_the_configured_domain(allowed):
    rates, _ = _network_rates()
    with pytest.raises(ValueError, match="allowed_molecular_uptake"):
        rates.at_uptake(allowed)


def test_a_reverse_named_nonpartner_does_not_redirect_the_uptake_cap():
    model = _network()
    _reaction(model, "EX_carbon_REV", {model.metabolites.oxygen_e: 1.0}, (0.0, 1000.0))
    assert cap_uptake(model, "EX_carbon", 2.0) == "EX_carbon"
    assert model.reactions.EX_carbon.lower_bound == -2.0
    assert model.reactions.EX_carbon_REV.upper_bound == 1000.0


def test_capping_a_split_supply_also_closes_an_accidentally_reopened_export():
    model = _network("split")
    model.reactions.EX_carbon.lower_bound = -1000.0
    assert cap_uptake(model, "EX_carbon", 2.0) == "EX_carbon_REV"
    assert model.reactions.EX_carbon.lower_bound == 0.0
    model.objective = "growth"
    assert model.slim_optimize() == pytest.approx(0.2)


@pytest.mark.parametrize("magnitude", [math.nan, math.inf, -math.inf])
def test_nonfinite_exchange_caps_are_rejected_without_mutation(magnitude):
    model = _network("split")
    before = _snapshot(model)
    with pytest.raises(ValueError, match="finite"):
        cap_uptake(model, "EX_carbon", magnitude)
    assert _snapshot(model) == before


@pytest.mark.parametrize("violation", ["bounds", "balance", "shared_capacity"])
def test_an_optimal_status_cannot_hide_an_invalid_complete_flux(monkeypatch, violation):
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    model.add_cons_vars([model.problem.Constraint(
        model.reactions.growth.flux_expression, ub=0.5, name="shared_capacity")])
    solver_type = type(model.solver)
    original = solver_type.primal_values.fget
    growth_variable = model.reactions.growth.forward_variable.name
    supply_variable = model.reactions.EX_carbon.reverse_variable.name
    transport_variable = model.reactions.transport.forward_variable.name
    synthesis_variable = model.reactions.synthesis.forward_variable.name
    maintenance_variable = model.reactions.ATPM.forward_variable.name

    def invalid_primal(solver):
        values = dict(original(solver))
        if violation == "bounds":
            values[maintenance_variable] = -1.0
        elif violation == "balance":
            values[synthesis_variable] += 0.1
        else:
            values[growth_variable] += 0.01
            values[supply_variable] += 0.1
            values[transport_variable] += 0.1
        return values

    monkeypatch.setattr(solver_type, "primal_values", property(invalid_primal))
    match = {"bounds": "bound violation", "balance": r"S\*v=0",
             "shared_capacity": "constraint violated.*shared_capacity"}[violation]
    with pytest.raises(RuntimeError, match=match):
        _network_rates(model)
    assert model.tolerance == PINNED_TOLERANCE


def test_the_scalar_and_network_routes_share_the_same_nonbinding_capacity_limit():
    model = _network()
    model.reactions.ATPM.lower_bound = RESTING_NGAM
    scalar, scalar_info = _build(model)
    network, network_info = _network_rates(model, allocation_fraction=1.0)
    for allowed in (0.71, 1.5, 5.0, 10.0):
        assert scalar.at_uptake(allowed) == pytest.approx(network.at_uptake(allowed), abs=1e-7)
    for info in (scalar_info, network_info):
        assert "pinned GLPK" in info.solver_method
        assert "presolve off" in info.solver_method
        assert info.solver_working_feasibility_tolerance == 1e-10
        assert info.physical_validation_tolerance == PINNED_TOLERANCE


@pytest.mark.parametrize("solve_raises", [False, True])
def test_working_precision_does_not_change_the_pinned_configuration_or_lp(monkeypatch, solve_raises):
    model = _network()
    model.tolerance = 1e-6
    model.solver.configuration.presolve = True
    model.solver.configuration._smcp.meth = GLP_DUAL
    model.solver.configuration._smcp.it_lim = 23456
    before = _snapshot(model)
    original = dynamic_rates.growth_or_none

    def inspect_configuration(scoped):
        assert scoped.solver.configuration.tolerances.feasibility == 1e-10
        assert scoped.solver.configuration.presolve is False
        assert scoped.solver.configuration._smcp.meth == GLP_PRIMAL
        assert scoped.solver.configuration._smcp.it_lim == 10000
        assert scoped.tolerance == 1e-6
        if solve_raises:
            raise RuntimeError("failed inside solver")
        return original(scoped)

    monkeypatch.setattr(dynamic_rates, "growth_or_none", inspect_configuration)
    stoichiometry = cobra.util.array.create_stoichiometric_matrix(model)
    if solve_raises:
        with pytest.raises(RuntimeError, match="failed inside solver"):
            dynamic_rates._validated_optimize(model, stoichiometry, "configuration regression")
    else:
        result = dynamic_rates._validated_optimize(model, stoichiometry, "configuration regression")
        assert result is not None
    assert _snapshot(model) == before
    assert model.solver.configuration.tolerances.feasibility == 1e-6
    assert model.solver.configuration.presolve is True
    assert model.solver.configuration._smcp.meth == GLP_DUAL
    assert model.solver.configuration._smcp.it_lim == 23456


@pytest.mark.integration
def test_structural_glutathione_grid_avoids_primal_warm_start_stall(ec_yeast_gem, monkeypatch):
    from ystwin.fba.product_panel import install_product, prepare_panel_model
    from ystwin.fba.storage import separate_storage_biomass
    from ystwin.in_silico import MetabolicControl, apply_synthetic_control

    prepared, _ = prepare_panel_model(ec_yeast_gem)
    installed, task = install_product(prepared, "glutathione")
    structural, basis = separate_storage_biomass(installed)
    control = MetabolicControl(0.9744394316690723, 0.9794153513556764,
                               0.7176590019151449, 0.16296056539068168)
    model, _ = apply_synthetic_control(structural, control,
                                      structural_mass_fraction=basis.structural_mass_fraction)
    before = _snapshot(model)
    original = dynamic_rates._validated_optimize
    original_solve = dynamic_rates.growth_or_none
    stages, attempts = [], []

    def record_solve(scoped):
        method = scoped.solver.configuration._smcp.meth
        value = original_solve(scoped)
        attempts.append((method, scoped.solver.status))
        return value

    def bounded_solve(scoped, matrix, stage, **kwargs):
        parameters = scoped.solver.configuration._smcp
        limit = parameters.it_lim
        try:
            parameters.it_lim = 10000
            result = original(scoped, matrix, stage, **kwargs)
            stages.append(stage)
            return result
        finally:
            parameters.it_lim = limit

    monkeypatch.setattr(dynamic_rates, "_validated_optimize", bounded_solve)
    monkeypatch.setattr(dynamic_rates, "growth_or_none", record_solve)
    rates, info = build_network_rates(
        model, "r_1714", "r_1992", "r_4046", product_reaction_id=task.reaction_id,
        uptake_vmax_mmol_per_gdcw_h=20 / 6 / basis.structural_mass_fraction,
        oxygen_lower_bound=-1000 / basis.structural_mass_fraction,
        growth_retention=control.growth_retention, allocation_fraction=control.allocation_fraction,
        grid_points=5)

    assert "attainable growth at uptake" in stages
    assert "dual recovery" in info.solver_method
    assert (GLP_DUAL, "optimal") in attempts
    assert all(np.isfinite(rates.at_uptake(value)).all() for value in rates.uptake_grid)
    assert _snapshot(model) == before


@pytest.mark.integration
@pytest.mark.parametrize("product_name", [
    "squalene", "glycogen", "trehalose", "glutathione", "glycerol",
])
def test_frozen_ec_panel_grid_has_physically_valid_complete_solutions(
        ec_yeast_gem, monkeypatch, product_name):
    from ystwin.fba.product_panel import install_product, prepare_panel_model

    prepared, _ = prepare_panel_model(ec_yeast_gem)
    model, task = install_product(prepared, product_name)
    before = _snapshot(model)
    original = dynamic_rates._validated_optimize
    checked_stages = []
    physical_stoichiometry = None

    def check_complete_solution(scoped, stoichiometry, stage, **kwargs):
        nonlocal physical_stoichiometry
        if physical_stoichiometry is None:
            physical_stoichiometry = cobra.util.array.create_stoichiometric_matrix(
                scoped, array_type="lil").tocsr()
        result = original(scoped, stoichiometry, stage, **kwargs)
        if result is None:
            return None
        value, fluxes = result
        assert np.isfinite(value)
        assert scoped.solver.status == "optimal"
        assert scoped.tolerance == PINNED_TOLERANCE
        assert scoped.solver.configuration.tolerances.feasibility == PINNED_TOLERANCE
        np.testing.assert_allclose(
            physical_stoichiometry @ fluxes, 0.0, atol=PINNED_TOLERANCE, rtol=0.0,
            err_msg=stage)
        for reaction, flux in zip(scoped.reactions, fluxes):
            assert reaction.lower_bound - PINNED_TOLERANCE <= flux
            assert flux <= reaction.upper_bound + PINNED_TOLERANCE
        primal = scoped.solver.primal_values
        for variable in scoped.solver.variables:
            observed = primal[variable.name]
            assert np.isfinite(observed)
            assert variable.lb is None or observed >= variable.lb - PINNED_TOLERANCE
            assert variable.ub is None or observed <= variable.ub + PINNED_TOLERANCE
        for constraint in scoped.solver.constraints:
            observed = constraint.primal
            assert np.isfinite(observed)
            assert constraint.lb is None or observed >= constraint.lb - PINNED_TOLERANCE
            assert constraint.ub is None or observed <= constraint.ub + PINNED_TOLERANCE
        checked_stages.append(stage)
        return result

    monkeypatch.setattr(dynamic_rates, "_validated_optimize", check_complete_solution)
    rates, info = build_network_rates(
        model, "r_1714", "r_1992", "r_4046", product_reaction_id=task.reaction_id,
        uptake_vmax_mmol_per_gdcw_h=20.0 / 6.0, grid_points=5)
    assert _snapshot(model) == before
    assert info.attainable_growth_per_h > 0.0
    assert rates.at_uptake(20.0 / 6.0)[0] > 0.0
    assert rates.at_uptake(20.0 / 6.0)[2] > 0.0
    assert "reference attainable growth" in checked_stages
    assert checked_stages.count("minimum actual carbon uptake") == int(rates.feasible.sum())
    assert np.all(rates.grid_fluxes[rates.feasible, 1]
                  <= rates.uptake_grid[rates.feasible] + PINNED_TOLERANCE)


@pytest.mark.integration
def test_frozen_ec_genuinely_infeasible_maintenance_is_explicitly_refused(ec_yeast_gem):
    from ystwin.fba.product_panel import install_product, prepare_panel_model

    prepared, _ = prepare_panel_model(ec_yeast_gem)
    model, task = install_product(prepared, "squalene")
    maintenance = model.reactions.get_by_id("r_4046")
    assert maintenance.lower_bound > PINNED_TOLERANCE
    model.add_cons_vars([model.problem.Constraint(
        maintenance.flux_expression, ub=0.0, name="contradictory_maintenance_limit")])
    before = _snapshot(model)
    with pytest.raises(RuntimeError, match="metabolic solve failed.*infeasible"):
        build_network_rates(
            model, "r_1714", "r_1992", "r_4046", product_reaction_id=task.reaction_id,
            uptake_vmax_mmol_per_gdcw_h=20.0 / 6.0, grid_points=5)
    assert _snapshot(model) == before
