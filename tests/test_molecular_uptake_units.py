import math

import cobra
import pytest

from ystwin.fba.physiology import PhysiologyReference, cap_uptake, validate_physiology


_REPRESENTATIONS = (
    "negative", "positive", "split", "negative_split", "positive_split", "negative_source_split",
)


def _reaction(model, rid, metabolites, bounds):
    reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites(metabolites)
    model.add_reactions([reaction])
    return reaction


def _exchange(model, rid, metabolite, representation, coefficient=2.0):
    shapes = {
        "negative": ((-1, (-1.0, 1000.0)), None),
        "positive": ((1, (0.0, 1.0)), None),
        "split": ((-1, (0.0, 0.0)), (1, (0.0, 1.0))),
        "negative_split": ((1, (-1000.0, 0.0)), (-1, (-1.0, 0.0))),
        "positive_split": ((1, (0.0, 1.0)), (-1, (0.0, 1000.0))),
        "negative_source_split": ((-1, (-1.0, 0.0)), (1, (-1000.0, 0.0))),
    }
    forward, reverse = shapes[representation]
    _reaction(model, rid, {metabolite: forward[0] * coefficient}, forward[1])
    if reverse is not None:
        _reaction(model, f"{rid}_REV", {metabolite: reverse[0] * coefficient}, reverse[1])
    return f"{rid}_REV" if representation in {"split", "negative_split"} else rid


def _model(representation, coefficient=2.0, closed=False):
    model = cobra.Model(representation)
    extracellular = cobra.Metabolite("substrate_e", compartment="e")
    intracellular = cobra.Metabolite("substrate_c", compartment="c")
    supply_id = _exchange(model, "r_1714", extracellular, representation, coefficient)
    _reaction(model, "transport", {extracellular: -1.0, intracellular: 1.0}, (0.0, 1000.0))
    growth = _reaction(model, "r_2111", {intracellular: -1.0}, (0.0, 1000.0))
    _reaction(model, "maintenance", {intracellular: -1.0}, (0.0, 1000.0))
    model.objective = growth
    if closed:
        model.reactions.get_by_id(supply_id).bounds = (0.0, 0.0)
    return model, supply_id


def _bounds(model):
    return {reaction.id: reaction.bounds for reaction in model.reactions}


@pytest.mark.parametrize("representation", _REPRESENTATIONS)
@pytest.mark.parametrize("closed", [False, True])
@pytest.mark.parametrize("magnitude", [6.0, -6.0])
def test_molecular_cap_uses_the_actual_coefficient_and_supply_direction(
    representation, closed, magnitude,
):
    model, supply_id = _model(representation, closed=closed)

    assert cap_uptake(model, "r_1714", magnitude) == supply_id

    supply = model.reactions.get_by_id(supply_id)
    coefficient = next(iter(supply.metabolites.values()))
    bound = supply.upper_bound if coefficient > 0 else supply.lower_bound
    assert coefficient * bound == pytest.approx(abs(magnitude))
    solution = model.optimize(raise_error=True)
    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(abs(magnitude))
    assert coefficient * solution.fluxes[supply_id] == pytest.approx(abs(magnitude))


@pytest.mark.parametrize("coefficient", [0.5, 1.0, 2.0, 3.0])
def test_split_conversion_uses_the_supply_coefficient_not_the_export_coefficient(coefficient):
    model, supply_id = _model("split")
    supply = model.reactions.get_by_id(supply_id)
    supply.add_metabolites({model.metabolites.substrate_e: coefficient}, combine=False)

    assert cap_uptake(model, "r_1714", 6.0) == supply_id
    assert supply.upper_bound == pytest.approx(6.0 / coefficient)
    assert model.optimize(raise_error=True).objective_value == pytest.approx(6.0)


@pytest.mark.parametrize("representation", ["positive_split", "negative_source_split"])
def test_a_closed_export_partner_does_not_displace_an_open_supply(representation):
    model, supply_id = _model(representation)
    model.reactions.r_1714_REV.bounds = (0.0, 0.0)

    assert cap_uptake(model, "r_1714", 6.0) == supply_id
    assert model.reactions.r_1714_REV.bounds == (0.0, 0.0)
    assert model.optimize(raise_error=True).objective_value == pytest.approx(6.0)


def test_split_cap_closes_a_reopened_export_before_solving_in_molecular_units():
    model, supply_id = _model("split")
    model.reactions.r_1714.lower_bound = -1000.0

    assert cap_uptake(model, "r_1714", 6.0) == supply_id
    assert model.reactions.r_1714.bounds == (0.0, 0.0)
    assert model.optimize(raise_error=True).objective_value == pytest.approx(6.0)


@pytest.mark.parametrize("representation", _REPRESENTATIONS)
def test_zero_molecular_supply_can_be_feasible_with_zero_biomass(representation):
    model, _ = _model(representation)

    cap_uptake(model, "r_1714", 0.0)

    solution = model.optimize(raise_error=True)
    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(0.0)
    assert solution.fluxes.r_2111 == pytest.approx(0.0)


@pytest.mark.parametrize("representation", _REPRESENTATIONS)
def test_zero_supply_with_required_maintenance_is_infeasible_not_zero_biomass(representation):
    model, _ = _model(representation)
    model.reactions.maintenance.lower_bound = 1.0

    cap_uptake(model, "r_1714", 0.0)

    with pytest.raises(cobra.exceptions.Infeasible):
        model.slim_optimize(error_value=None)
    assert model.solver.status == "infeasible"


@pytest.mark.parametrize("magnitude", [
    math.nan, math.inf, -math.inf, True, False, None, "6", 6j,
    pytest.param(10**1000, id="unrepresentable-integer"),
])
def test_invalid_molecular_caps_fail_before_mutating_any_bounds(magnitude):
    model, _ = _model("split")
    model.reactions.r_1714.lower_bound = -1000.0
    before = _bounds(model)

    with pytest.raises(ValueError, match="finite"):
        cap_uptake(model, "r_1714", magnitude)

    assert _bounds(model) == before


def test_nonfinite_converted_reaction_flux_is_rejected_before_mutation():
    model, _ = _model("split", coefficient=0.5)
    model.reactions.r_1714.lower_bound = -1000.0
    before = _bounds(model)

    with pytest.raises(ValueError, match="finite"):
        cap_uptake(model, "r_1714", 1e308)

    assert _bounds(model) == before


@pytest.mark.parametrize("shape", ["empty", "transport"])
def test_uptake_requires_a_single_metabolite_boundary(shape):
    model, _ = _model("negative")
    reaction = model.reactions.r_1714
    if shape == "empty":
        reaction.subtract_metabolites(dict(reaction.metabolites))
    else:
        reaction.add_metabolites({model.metabolites.substrate_c: 1.0})
    before = _bounds(model)

    with pytest.raises(ValueError, match="one-metabolite boundary"):
        cap_uptake(model, "r_1714", 6.0)

    assert _bounds(model) == before


def test_a_reverse_named_transport_is_not_a_partner_even_when_it_shares_the_substrate():
    model, _ = _model("split")
    model.reactions.r_1714_REV.add_metabolites({model.metabolites.substrate_c: -2.0})
    before = model.reactions.r_1714_REV.bounds

    assert cap_uptake(model, "r_1714", 6.0) == "r_1714"
    assert model.reactions.r_1714_REV.bounds == before
    assert model.optimize(raise_error=True).objective_value == pytest.approx(6.0)


def test_a_reverse_transport_cannot_itself_be_capped_as_a_boundary():
    model, _ = _model("split")
    model.reactions.r_1714_REV.add_metabolites({model.metabolites.substrate_c: -2.0})
    before = _bounds(model)

    with pytest.raises(ValueError, match="one-metabolite boundary"):
        cap_uptake(model, "r_1714_REV", 6.0)

    assert _bounds(model) == before


@pytest.mark.parametrize("rid", ["r_1714", "r_1714_REV"])
@pytest.mark.parametrize("coefficient", [0.0, math.nan, math.inf, -math.inf])
def test_invalid_exchange_stoichiometry_cannot_open_an_alternative_supply(rid, coefficient):
    model, _ = _model("split")
    reaction = model.reactions.get_by_id(rid)
    reaction._metabolites[model.metabolites.substrate_e] = coefficient
    before = _bounds(model)

    with pytest.raises(ValueError, match="stoichiometry"):
        cap_uptake(model, "r_1714", 6.0)

    assert _bounds(model) == before


def test_a_required_flux_conflict_does_not_partially_close_the_split_pair():
    model, _ = _model("split")
    model.reactions.r_1714.lower_bound = -1000.0
    model.reactions.r_1714_REV.bounds = (4.0, 1000.0)
    before = _bounds(model)

    with pytest.raises(ValueError, match="required flux"):
        cap_uptake(model, "r_1714", 6.0)

    assert _bounds(model) == before


def _physiology_model(representation):
    model, supply_id = _model(representation)
    oxygen = cobra.Metabolite("oxygen_e", compartment="e")
    ethanol = cobra.Metabolite("ethanol_e", compartment="e")
    co2 = cobra.Metabolite("co2_e", compartment="e")
    oxygen_id = _exchange(model, "r_1992", oxygen, representation)
    for rid in (supply_id, oxygen_id):
        reaction = model.reactions.get_by_id(rid)
        coefficient = next(iter(reaction.metabolites.values()))
        if coefficient > 0:
            reaction.upper_bound = 2.0
        else:
            reaction.lower_bound = -2.0
    for rid, metabolite in (("r_1761", ethanol), ("r_1672", co2)):
        if representation == "positive":
            _reaction(model, rid, {metabolite: 2.0}, (-1000.0, 0.0))
        else:
            _reaction(model, rid, {metabolite: -2.0}, (0.0, 1000.0))
        if representation == "split":
            _reaction(model, f"{rid}_REV", {metabolite: 2.0}, (0.0, 0.0))
    model.reactions.r_2111.add_metabolites({oxygen: -1.0, ethanol: 2.0, co2: 3.0})
    reference = PhysiologyReference("toy", 4.0, 4.0, 4.0, 8.0, 12.0, "toy mass balance")
    return model, reference


@pytest.mark.parametrize("representation", ["negative", "positive", "split"])
@pytest.mark.parametrize("apply_constraints", [False, True])
def test_physiology_reports_molecular_rates_instead_of_raw_exchange_flux(
    representation, apply_constraints,
):
    model, reference = _physiology_model(representation)
    before = _bounds(model)

    report = validate_physiology(model, reference, apply_constraints=apply_constraints)

    assert report.observed == pytest.approx({
        "growth_rate": 4.0, "glucose": 4.0, "oxygen": 4.0, "ethanol": 8.0, "co2": 12.0,
    })
    assert report.passed
    assert _bounds(model) == before


def test_physiology_reports_the_net_molecular_rate_of_both_split_branches():
    model, reference = _physiology_model("split")
    model.reactions.r_1714.bounds = (1.0, 1.0)

    report = validate_physiology(model, reference, apply_constraints=False)

    assert report.observed == pytest.approx({
        "growth_rate": 2.0, "glucose": 2.0, "oxygen": 2.0, "ethanol": 4.0, "co2": 6.0,
    })


def test_physiology_does_not_report_an_infeasible_solution_as_zero_growth():
    model, reference = _physiology_model("split")
    model.reactions.maintenance.lower_bound = 10.0
    before = _bounds(model)

    with pytest.raises(cobra.exceptions.OptimizationError, match="infeasible"):
        validate_physiology(model, reference)

    assert _bounds(model) == before
