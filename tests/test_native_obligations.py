from __future__ import annotations

from copy import deepcopy
from dataclasses import fields, replace
import inspect
from pathlib import Path

import cobra
from cobra.core.solution import get_solution
from cobra.util.array import create_stoichiometric_matrix
import numpy as np
import pytest

from ystwin.fba.dynamic_rates import _validated_optimize
from ystwin.fba.native_obligations import (
    NativeCarbonObligation,
    NativeObligationRoles,
    impose_native_obligation,
)
from ystwin.fba.physiology import aerobic_batch_constraints
from ystwin.fba.solver import configure


def _reaction(model, rid, stoichiometry, bounds=(0., 1000.)):
    reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites(stoichiometry)
    model.add_reactions([reaction])
    return reaction


@pytest.fixture
def native_toy():
    model = cobra.Model("native_physiology")
    glucose = cobra.Metabolite("glucose_e", formula="C6H12O6", compartment="e", charge=0)
    triose = cobra.Metabolite("triose_c", formula="C3H6O3", compartment="c", charge=0)
    glycerol = cobra.Metabolite("glycerol_c", formula="C3H8O3", compartment="c", charge=0)
    glycerol_e = cobra.Metabolite("glycerol_e", formula="C3H8O3", compartment="e", charge=0)
    nadh = cobra.Metabolite("nadh", formula="C21H27N7O14P2", compartment="c", charge=-2)
    nad = cobra.Metabolite("nad", formula="C21H26N7O14P2", compartment="c", charge=-1)
    hydrogen = cobra.Metabolite("hydrogen", formula="H2", compartment="c", charge=0)
    proton = cobra.Metabolite("proton", formula="H", compartment="c", charge=1)
    protein = cobra.Metabolite("protein_pool", compartment="c")
    _reaction(model, "glucose_boundary", {glucose: -1}, (-10., 1000.))
    _reaction(model, "carbon_split", {glucose: -1, triose: 2})
    _reaction(model, "hydrogen_supply", {hydrogen: 1}, (0., 100.))
    _reaction(model, "reduction", {hydrogen: -1, nad: -1, nadh: 1, proton: 1})
    _reaction(model, "native_synthesis", {triose: -1, nadh: -1, proton: -1,
                                          glycerol: 1, nad: 1, protein: -2})
    _reaction(model, "native_transport", {glycerol: -1, glycerol_e: 1})
    _reaction(model, "native_export", {glycerol_e: -1})
    _reaction(model, "protein_supply", {protein: 1}, (0., 100.))
    _reaction(model, "growth", {triose: -1, protein: -1})
    model.objective = "growth"
    model.notes = {"native_context": {"carbon": "glucose"}}
    roles = NativeObligationRoles(
        native_metabolite_id="glycerol_c", glucose_metabolite_id="glucose_e",
        glucose_exchange_ids=("glucose_boundary",), native_export_ids=("native_export",),
        resource_metabolite_ids=("protein_pool",),
    )
    return model, roles


def _obligation(fraction=.25):
    return NativeCarbonObligation(
        carbon_fraction=fraction, source="native balance test coefficient",
        source_medium="glucose-only toy medium",
        metadata={"native_training_scope": "glucose_only", "native_parameters": {"k": 2.}},
    )


def _impose(model, roles, fraction=.25):
    return impose_native_obligation(model, roles=roles, obligation=_obligation(fraction),
                                    demand_reaction_id="native_commitment")


def _snapshot(model):
    return {
        "reactions": {r.id: (r.bounds, {m.id: c for m, c in r.metabolites.items()}) for r in model.reactions},
        "constraints": {c.name: (c.lb, c.ub, str(c.expression)) for c in model.constraints},
        "objective": str(model.objective.expression),
        "direction": model.objective.direction,
        "notes": deepcopy(model.notes),
        "solver": model.solver.interface.__name__,
        "feasibility": model.solver.configuration.tolerances.feasibility,
    }


def _fluxes(report, solution):
    uptake = sum(c * solution.fluxes[rid] for rid, c in report.glucose_flux_coefficients.items())
    native = sum(c * solution.fluxes[rid] for rid, c in report.terminal_flux_coefficients.items())
    return uptake, native


def test_native_obligation_is_a_real_carbon_and_cofactor_cost_in_the_same_solve(native_toy):
    model, roles = native_toy
    assert model.optimize().objective_value == pytest.approx(20.)
    before = _snapshot(model)
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(15.)
    uptake, native = _fluxes(report, solution)
    assert uptake == pytest.approx(10.)
    assert native == pytest.approx(5.)
    assert solution.fluxes["native_synthesis"] == pytest.approx(5.)
    assert solution.fluxes["reduction"] == pytest.approx(5.)
    assert solution.fluxes["protein_supply"] == pytest.approx(25.)
    np.testing.assert_allclose(create_stoichiometric_matrix(constrained) @ solution.fluxes.to_numpy(),
                               0., atol=1e-8)
    demand = constrained.reactions.get_by_id(report.demand_reaction_id)
    assert demand.metabolites == {constrained.metabolites.get_by_id("glycerol_c"): -1.}
    assert demand.lower_bound == 0
    assert report.native_moles_per_glucose_mole == pytest.approx(.5)
    assert report.metadata["native_training_scope"] == "glucose_only"
    assert report.metadata["fate_specific"] is False
    assert report.metadata["absolute_source_flux_identified"] is False
    assert _snapshot(model) == before
    constrained.notes["native_context"]["carbon"] = "changed"
    assert model.notes["native_context"]["carbon"] == "glucose"


def test_zero_coefficient_does_not_open_an_unforced_native_disposal_bypass(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("native_transport").bounds = (0., 0.)
    model.reactions.get_by_id("carbon_split").add_metabolites({
        model.metabolites.get_by_id("triose_c"): 1.,
        model.metabolites.get_by_id("glycerol_c"): 1.,
        model.metabolites.get_by_id("hydrogen"): -1.,
    }, combine=False)
    baseline = model.optimize().objective_value
    assert baseline == pytest.approx(0.)
    constrained, report = _impose(model, roles, fraction=0.)
    assert constrained.optimize().objective_value == pytest.approx(baseline)
    assert constrained.reactions.get_by_id(report.demand_reaction_id).bounds == (0., 0.)


def test_obligation_uses_utilized_glucose_not_an_unused_uptake_ceiling(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("growth").bounds = (8., 8.)
    constrained, report = _impose(model, roles)
    exchange = constrained.reactions.get_by_id("glucose_boundary")
    constrained.objective = -exchange.flux_expression
    constrained.objective.direction = "min"
    solution = constrained.optimize()
    uptake, native = _fluxes(report, solution)
    assert uptake == pytest.approx(8. / 1.5)
    assert native == pytest.approx(.5 * uptake)
    assert native < .5 * abs(exchange.lower_bound)
    assert solution.fluxes["growth"] == pytest.approx(8.)


def test_existing_native_exports_satisfy_total_commitment_without_double_charging(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("native_export").lower_bound = 5.
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(model.optimize().objective_value)
    assert solution.fluxes[report.demand_reaction_id] == pytest.approx(0.)
    assert _fluxes(report, solution) == pytest.approx((10., 5.))
    assert report.terminal_flux_coefficients == {"native_commitment": 1., "native_export": 1.}


def test_native_commitment_does_not_depend_on_an_export_route_being_open(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("native_transport").bounds = (0., 0.)
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(15.)
    assert solution.fluxes["native_export"] == pytest.approx(0.)
    assert solution.fluxes[report.demand_reaction_id] == pytest.approx(5.)
    assert "retention" not in report.metadata["demand_interpretation"].split()[0]


def test_shared_protein_and_other_hard_constraints_remain_binding(native_toy):
    model, roles = native_toy
    constraint = model.problem.Constraint(model.reactions.get_by_id("protein_supply").flux_expression,
                                          ub=10., name="hard_shared_protein")
    model.add_cons_vars([constraint])
    before = _snapshot(model)
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(6.)
    assert _fluxes(report, solution) == pytest.approx((4., 2.))
    assert solution.fluxes["protein_supply"] == pytest.approx(10.)
    assert constrained.constraints["hard_shared_protein"].ub == 10.
    assert _snapshot(model) == before


def test_native_redox_capacity_cannot_be_bypassed(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("hydrogen_supply").upper_bound = 1.
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(3.)
    assert _fluxes(report, solution) == pytest.approx((2., 1.))


def test_infeasible_native_obligation_is_refused_without_weakening_caller(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("growth").lower_bound = 15.
    model.reactions.get_by_id("native_synthesis").upper_bound = 1.
    before = _snapshot(model)
    with pytest.raises(ValueError, match="infeasible"):
        _impose(model, roles)
    assert _snapshot(model) == before
    assert model.optimize().status == "optimal"


@pytest.mark.parametrize("fraction", [-.01, 1.01, np.nan, np.inf, True])
def test_nonphysical_carbon_commitments_cannot_be_clipped_into_a_flux(fraction):
    with pytest.raises(ValueError):
        _obligation(fraction)


@pytest.mark.parametrize("rid", ["native_export", "carbon_reserve"])
def test_free_carbon_uptake_or_inventory_release_is_refused(native_toy, rid):
    model, roles = native_toy
    if rid == "native_export":
        model.reactions.get_by_id(rid).lower_bound = -1.
    else:
        _reaction(model, rid, {model.metabolites.get_by_id("triose_c"): 1.}, (0., 1.))
    before = _snapshot(model)
    with pytest.raises(ValueError, match="carbon|uptake|remobilization"):
        _impose(model, roles)
    assert _snapshot(model) == before


def test_other_carbon_media_are_not_silently_treated_as_glucose_only(native_toy):
    model, roles = native_toy
    acetate = cobra.Metabolite("acetate_e", formula="C2H3O2", compartment="e")
    _reaction(model, "acetate_boundary", {acetate: 1.}, (0., 3.))
    with pytest.raises(ValueError, match="carbon|medium"):
        _impose(model, roles)
    with pytest.raises(ValueError, match="carbon|medium"):
        _impose(model, replace(roles, resource_metabolite_ids=("protein_pool", "acetate_e")))


def test_unclassified_pseudo_resource_supply_is_not_assumed_carbon_free(native_toy):
    model, roles = native_toy
    with pytest.raises(ValueError, match="resource|formula"):
        _impose(model, replace(roles, resource_metabolite_ids=()))


def test_unknown_formula_is_not_evidence_of_carbon_free_input(native_toy):
    model, roles = native_toy
    ambiguous = cobra.Metabolite("unclassified_input", formula="R", compartment="e")
    _reaction(model, "unclassified_boundary", {ambiguous: 1.}, (0., 2.))
    with pytest.raises(ValueError, match="formula"):
        _impose(model, roles)


@pytest.mark.parametrize("coefficient", [-2., 2.])
def test_split_or_nonunit_glucose_boundaries_use_molecular_net_input(native_toy, coefficient):
    model, roles = native_toy
    original = model.reactions.get_by_id("glucose_boundary")
    original.bounds = (0., 0.)
    bounds = (-5., 0.) if coefficient < 0 else (0., 5.)
    _reaction(model, "glucose_input", {model.metabolites.get_by_id("glucose_e"): coefficient}, bounds)
    roles = replace(roles, glucose_exchange_ids=("glucose_boundary", "glucose_input"))
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(15.)
    assert _fluxes(report, solution) == pytest.approx((10., 5.))
    assert report.glucose_flux_coefficients["glucose_input"] == coefficient


def test_exchange_cycling_does_not_turn_gross_glucose_supply_into_utilized_carbon(native_toy):
    model, roles = native_toy
    model.reactions.get_by_id("glucose_boundary").bounds = (3., 3.)
    _reaction(model, "glucose_input", {model.metabolites.get_by_id("glucose_e"): 1.}, (8., 8.))
    roles = replace(roles, glucose_exchange_ids=("glucose_boundary", "glucose_input"))
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(7.5)
    assert _fluxes(report, solution) == pytest.approx((5., 2.5))


def test_nonunit_reversed_native_export_counts_molecules_not_reaction_direction(native_toy):
    model, roles = native_toy
    export = model.reactions.get_by_id("native_export")
    export.add_metabolites({model.metabolites.get_by_id("glycerol_e"): 2.}, combine=False)
    export.bounds = (-1000., 0.)
    constrained, report = _impose(model, roles)
    solution = constrained.optimize()
    assert solution.objective_value == pytest.approx(15.)
    assert _fluxes(report, solution) == pytest.approx((10., 5.))
    assert report.terminal_flux_coefficients["native_export"] == -2.


def test_missing_split_glucose_input_role_is_refused(native_toy):
    model, roles = native_toy
    _reaction(model, "additional_glucose", {model.metabolites.get_by_id("glucose_e"): 1.}, (0., 2.))
    with pytest.raises(ValueError, match="glucose|carbon"):
        _impose(model, roles)


def test_net_glucose_export_cannot_earn_a_negative_native_demand(native_toy):
    model, roles = native_toy
    constrained, report = _impose(model, roles)
    constrained.objective = constrained.reactions.get_by_id(report.demand_reaction_id)
    constrained.objective.direction = "min"
    constrained.reactions.get_by_id("glucose_boundary").lower_bound = 0.
    solution = constrained.optimize()
    assert solution.status == "optimal"
    assert _fluxes(report, solution) == pytest.approx((0., 0.))
    assert solution.fluxes["growth"] == pytest.approx(0.)
    assert constrained.constraints[report.uptake_constraint_id].lb == 0.


def test_roles_use_stoichiometry_and_chemical_identity_not_identifier_patterns(native_toy):
    model, roles = native_toy
    with pytest.raises(ValueError, match="boundary|one-metabolite"):
        _impose(model, replace(roles, native_export_ids=("native_transport",)))
    with pytest.raises(ValueError, match="glucose|formula"):
        _impose(model, replace(roles, glucose_metabolite_id="glycerol_e"))
    with pytest.raises(ValueError, match="native|glycerol|formula"):
        _impose(model, replace(roles, native_metabolite_id="triose_c"))
    with pytest.raises(ValueError, match="native|glycerol|formula"):
        _impose(model, replace(roles, native_export_ids=("glucose_boundary",)))


def test_irrelevant_labels_cannot_change_native_obligations(native_toy):
    model, roles = native_toy
    other = model.copy()
    other.id = "irrelevant_label"
    other.notes["unrelated_label"] = "another_label"
    for reaction in other.reactions:
        reaction.name = "unrelated display name"
    left, left_report = _impose(model, roles)
    right, right_report = _impose(other, roles)
    assert left.optimize().objective_value == pytest.approx(right.optimize().objective_value)
    assert left_report.glucose_flux_coefficients == right_report.glucose_flux_coefficients
    assert left_report.terminal_flux_coefficients == right_report.terminal_flux_coefficients
    assert "product" not in inspect.signature(impose_native_obligation).parameters
    assert all("product" not in field.name for field in fields(NativeObligationRoles))
    assert all("product" not in field.name for field in fields(NativeCarbonObligation))


def test_coefficient_provenance_medium_dependence_and_caller_metadata_are_preserved(native_toy):
    model, roles = native_toy
    obligation = _obligation()
    before = deepcopy(obligation.metadata)
    constrained, report = impose_native_obligation(model, roles=roles, obligation=obligation,
                                                    demand_reaction_id="native_commitment")
    assert report.metadata["source_medium"] == obligation.source_medium
    assert report.metadata["source"] == obligation.source
    assert "medium" in report.metadata["medium_dependence"]
    assert report.metadata["source_metadata"] == before
    report.metadata["source_metadata"]["native_parameters"]["k"] = -1
    assert obligation.metadata == before
    assert report.coupling_constraint_id in constrained.constraints
    with pytest.raises(ValueError, match="already|collision"):
        _impose(constrained, roles)


@pytest.mark.integration
def test_actual_native_ec_gem_pays_source_derived_glycerol_commitment():
    from ystwin.mech.hog import HogModel, HogProtocol, hog_glycerol_balance

    path = Path(__file__).resolve().parents[1] / "data" / "gem" / "ecYeastGEM_batch.xml.gz"
    if not path.exists():
        pytest.skip("native enzyme-constrained GEM artifact is absent")
    model = cobra.io.read_sbml_model(str(path))
    aerobic_batch_constraints(model)
    roles = NativeObligationRoles(
        native_metabolite_id="s_0765[c]", glucose_metabolite_id="s_0565[e]",
        glucose_exchange_ids=("r_1714", "r_1714_REV"), native_export_ids=("r_1808",),
        resource_metabolite_ids=("prot_pool[c]",),
    )
    native = HogModel.from_source()
    trajectory = native.simulate([0., 3600., 5400.], protocol=HogProtocol(.4))
    balance = hog_glycerol_balance(native, trajectory)
    fraction = float(balance.carbon_commitment()[-1])
    obligation = NativeCarbonObligation(fraction, "WT source balances at 5400 s", "W303 YPD",
                                        metadata=balance.metadata)
    before = _snapshot(model)
    constrained, report = impose_native_obligation(model, roles=roles, obligation=obligation,
                                                    demand_reaction_id="native_commitment")
    assert _snapshot(model) == before
    stoichiometry = create_stoichiometric_matrix(constrained, array_type="lil")
    _validated_optimize(constrained, stoichiometry, "native constrained growth")
    solution = get_solution(constrained)
    assert solution.status == "optimal"
    uptake, native_flux = _fluxes(report, solution)
    assert native_flux >= report.native_moles_per_glucose_mole * uptake - 1e-7
    assert native_flux > 0
    baseline = model.copy()
    configure(baseline)
    baseline_value, _ = _validated_optimize(
        baseline, create_stoichiometric_matrix(baseline, array_type="lil"), "native baseline growth")
    assert solution.objective_value < baseline_value - 1e-5
    assert constrained.reactions.get_by_id("prot_pool_exchange").bounds == (
        model.reactions.get_by_id("prot_pool_exchange").bounds)
    np.testing.assert_allclose(stoichiometry @ solution.fluxes.to_numpy(), 0., atol=1e-7)
