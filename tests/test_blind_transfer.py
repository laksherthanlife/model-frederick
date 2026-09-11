from __future__ import annotations

from dataclasses import replace

import cobra
import pytest

from ystwin.analysis.blind_transfer import BalancedCondition, HostFluxRoles, content_from_flux, predict_balanced_envelope
from ystwin.fba.chemical_task import ChemicalTask, MetaboliteDefinition, ReactionDefinition
from ystwin.fba.native_obligations import NativeCarbonObligation, NativeObligationRoles


class NativePhysiology:
    def predict(self, rates):
        return [{"observable_id": name, "prediction_status": "quantified", "prediction": value}
                for name, value in (("GlucoseUptake", 2.0), ("O2uptake", 1.0))]


def native_host():
    model = cobra.Model("neutral_native_host")
    model.solver = "glpk"
    metabolites = {name: cobra.Metabolite(name, formula=formula, charge=0, compartment="c")
                   for name, formula in (("G", "C6H12O6"), ("B", "C6H12O6"),
                                         ("H", "H2"), ("O", "O2"), ("Y", "C3H8O3"))}
    definitions = {
        "glucose_feed": {"G": 1}, "hydrogen_feed": {"H": 1}, "oxygen_feed": {"O": 1},
        "assimilate": {"G": -1, "B": 1}, "grow": {"B": -1},
        "native_synthesis": {"G": -1, "H": -2, "Y": 2},
    }
    for name, coefficients in definitions.items():
        reaction = cobra.Reaction(name, lower_bound=0, upper_bound=1000)
        reaction.add_metabolites({metabolites[key]: value for key, value in coefficients.items()})
        model.add_reactions([reaction])
    roles = HostFluxRoles(
        biomass_reaction_id="grow", native=NativeObligationRoles("Y", "G", ("glucose_feed",)),
        oxygen_exchange_ids=("oxygen_feed",), intracellular_compartment_ids=("c",),
        extracellular_compartment_ids=(),
    )
    return model, roles


def chemical_task(label="unseen_a"):
    return ChemicalTask((MetaboliteDefinition(label, "C6H12O6", 0, "c"),),
                        (ReactionDefinition(f"make_{label}", (("G", -1), (label, 1)), "independent chemistry"),), label)


def condition():
    return BalancedCondition("experiment_a", 0.1, "glucose", "balanced_growth_intracellular_content",
                             "explicit native-condition transfer assumption")


def test_native_obligations_reduce_a_name_blind_capacity_envelope():
    host, roles = native_host()
    records = []
    for fraction in (0.0, 0.1):
        result, details = predict_balanced_envelope(
            host, task=chemical_task(), condition=condition(), physiology=NativePhysiology(),
            obligations=(NativeCarbonObligation(fraction, "native evidence", "glucose"),), roles=roles)
        assert result["status"] == "ok", result
        assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.9 - 2 * fraction)
        records.append(result)
    assert records[1]["upper"] < records[0]["upper"]
    assert records[0]["lower"] == pytest.approx(0)
    assert "__chemical_task_output" not in host.reactions


def test_only_chemical_structure_not_identifier_controls_the_forecast():
    host, roles = native_host()
    kwargs = {"condition": condition(), "physiology": NativePhysiology(),
              "obligations": (NativeCarbonObligation(0.1, "native evidence", "glucose"),), "roles": roles}
    first, _ = predict_balanced_envelope(host, task=chemical_task("unseen_a"), **kwargs)
    second, _ = predict_balanced_envelope(host, task=chemical_task("different_label"), **kwargs)
    assert first == second


def test_product_mass_enters_the_total_dry_mass_denominator():
    assert content_from_flux(1.0, 0.1, 100.0) == pytest.approx(500.0)
    with pytest.raises(ValueError):
        content_from_flux(1.0, 0.0, 100.0)


def test_unsupported_inputs_and_bad_chemistry_are_retained_not_zero_filled():
    host, roles = native_host()
    kwargs = {"physiology": NativePhysiology(),
              "obligations": (NativeCarbonObligation(0.1, "native evidence", "glucose"),), "roles": roles}
    unsupported, _ = predict_balanced_envelope(host, task=chemical_task(),
                                               condition=replace(condition(), carbon_input="unsupported_feed"), **kwargs)
    assert unsupported["status"] == "unsupported" and unsupported["upper"] is None
    broken = replace(chemical_task(), metabolites=(MetaboliteDefinition("unseen_a", "C5H10O5", 0, "c"),))
    failed, _ = predict_balanced_envelope(host, task=broken, condition=condition(), **kwargs)
    assert failed["status"] == "failed" and failed["upper"] is None


def predict(host, roles, task=None):
    return predict_balanced_envelope(
        host, task=chemical_task() if task is None else task, condition=condition(),
        physiology=NativePhysiology(),
        obligations=(NativeCarbonObligation(0.1, "native evidence", "glucose"),), roles=roles,
    )


def set_boundary(reaction, coefficient, capacity):
    metabolite = next(iter(reaction.metabolites))
    reaction.add_metabolites({metabolite: coefficient}, combine=False)
    reaction.bounds = (0.0, capacity) if coefficient > 0 else (-capacity, 0.0)


def oxygen_task():
    return ChemicalTask(
        (MetaboliteDefinition("unseen_oxygenated", "C6H12O8", 0, "c"),),
        (ReactionDefinition("make_oxygenated", (("G", -1), ("O", -1), ("unseen_oxygenated", 1)),
                            "independent balanced chemistry"),),
        "unseen_oxygenated",
    )


def add_boundary(host, identifier, metabolite_id, coefficient, bounds):
    reaction = cobra.Reaction(identifier, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites({host.metabolites.get_by_id(metabolite_id): coefficient})
    host.add_reactions([reaction])


def test_learned_glucose_cap_preserves_a_tighter_hard_bound():
    host, roles = native_host()
    host.reactions.glucose_feed.upper_bound = 0.5
    result, details = predict(host, roles)
    assert result["status"] == "ok", result
    scenario = details["scenarios"][0]
    assert scenario["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(0.35)
    assert scenario["native_obligation"]["applied_boundary_input_bounds"]["glucose_feed"] == (0, 0.5)
    assert host.reactions.glucose_feed.bounds == (0, 0.5)


def test_closed_glucose_supply_is_not_reopened():
    host, roles = native_host()
    host.reactions.glucose_feed.bounds = (0, 0)
    result, _ = predict(host, roles)
    assert result["status"] in {"failed", "infeasible"}
    assert result["lower"] is None and result["upper"] is None
    assert host.reactions.glucose_feed.bounds == (0, 0)


@pytest.mark.parametrize("coefficient", [2.0, -2.0, 0.5])
def test_nonunit_glucose_boundary_has_the_same_molecular_budget(coefficient):
    host, roles = native_host()
    set_boundary(host.reactions.glucose_feed, coefficient, 1000)
    before = host.reactions.glucose_feed.bounds
    result, details = predict(host, roles)
    assert result["status"] == "ok", result
    scenario = details["scenarios"][0]
    assert scenario["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.7)
    assert scenario["native_obligation"]["applied_boundary_input_bounds"]["glucose_feed"] == before
    assert host.reactions.glucose_feed.bounds == before


@pytest.mark.parametrize("coefficient,capacity,expected", [
    (1.0, 0.2, 0.2), (2.0, 1000, 1.0), (-2.0, 1000, 1.0),
    (-2.0, 0.1, 0.2), (1.0, 0.0, 0.0),
])
def test_oxygen_caps_preserve_hard_bounds_and_molecular_units(coefficient, capacity, expected):
    host, roles = native_host()
    set_boundary(host.reactions.oxygen_feed, coefficient, capacity)
    before = host.reactions.oxygen_feed.bounds
    result, details = predict(host, roles, oxygen_task())
    assert result["status"] == "ok", result
    scenario = details["scenarios"][0]
    assert scenario["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(expected)
    assert host.reactions.oxygen_feed.bounds == before
    if capacity:
        assert scenario["native_obligation"]["applied_boundary_input_bounds"]["oxygen_feed"] == before


def test_split_glucose_boundaries_cap_net_input_not_gross_exchange():
    host, roles = native_host()
    host.reactions.glucose_feed.bounds = (3, 3)
    add_boundary(host, "glucose_export", "G", -1, (1, 1))
    roles = replace(roles, native=replace(
        roles.native, glucose_exchange_ids=("glucose_feed", "glucose_export")))
    result, details = predict(host, roles)
    assert result["status"] == "ok", result
    assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.7)
    assert host.reactions.glucose_feed.bounds == (3, 3)
    assert host.reactions.glucose_export.bounds == (1, 1)


def test_all_declared_oxygen_boundaries_share_one_molecular_cap():
    host, roles = native_host()
    set_boundary(host.reactions.oxygen_feed, 2, 1000)
    add_boundary(host, "additional_oxygen", "O", 0.5, (0, 1000))
    roles = replace(roles, oxygen_exchange_ids=("oxygen_feed", "additional_oxygen"))
    result, details = predict(host, roles, oxygen_task())
    assert result["status"] == "ok", result
    assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.0)


def test_split_oxygen_boundaries_cap_net_input_not_gross_exchange():
    host, roles = native_host()
    set_boundary(host.reactions.oxygen_feed, 2, 1)
    host.reactions.oxygen_feed.lower_bound = 1
    add_boundary(host, "oxygen_export", "O", -1, (0, 1000))
    roles = replace(roles, oxygen_exchange_ids=("oxygen_feed", "oxygen_export"))
    result, details = predict(host, roles, oxygen_task())
    assert result["status"] == "ok", result
    assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.0)
    assert host.reactions.oxygen_feed.bounds == (1, 1)
    assert host.reactions.oxygen_export.bounds == (0, 1000)


def test_undeclared_active_oxygen_boundary_cannot_bypass_the_cap():
    host, roles = native_host()
    add_boundary(host, "undeclared_oxygen", "O", 1, (0, 1000))
    result, _ = predict(host, roles, oxygen_task())
    assert result["status"] == "failed"
    assert result["upper"] is None
    assert "oxygen" in result["reason"] and "role" in result["reason"]


def test_other_carbon_supplies_require_explicit_medium_setup_not_silent_closure():
    host, roles = native_host()
    host.add_metabolites([cobra.Metabolite("other_carbon", formula="C2H4O2", charge=0, compartment="c")])
    add_boundary(host, "other_carbon_feed", "other_carbon", 1, (0, 1))
    result, _ = predict(host, roles)
    assert result["status"] == "failed"
    assert result["upper"] is None
    assert host.reactions.other_carbon_feed.bounds == (0, 1)


@pytest.mark.parametrize("compartment", ["e", "outside"])
def test_explicit_extracellular_compartment_aliases_are_equally_rejected(compartment):
    host, roles = native_host()
    host.add_metabolites([cobra.Metabolite("external_marker", formula="H2O", charge=0,
                                         compartment=compartment)])
    host.compartments = {"c": "cytosol", compartment: "extracellular space"}
    roles = replace(roles, extracellular_compartment_ids=(compartment,))
    task = chemical_task()
    task = replace(task, metabolites=(replace(task.metabolites[0], compartment=compartment),))
    result, _ = predict(host, roles, task)
    assert result["status"] == "failed"
    assert result["upper"] is None
    assert "extracellular" in result["reason"]


@pytest.mark.parametrize("compartment", ["inside", "e"])
def test_explicit_renamed_intracellular_compartment_is_accepted(compartment):
    host, roles = native_host()
    for metabolite in host.metabolites:
        metabolite.compartment = compartment
    host.compartments = {compartment: "cytosol"}
    roles = replace(roles, intracellular_compartment_ids=(compartment,))
    task = chemical_task()
    task = replace(task, metabolites=(replace(task.metabolites[0], compartment=compartment),))
    result, details = predict(host, roles, task)
    assert result["status"] == "ok", result
    assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(1.7)


def test_unclassified_output_compartment_is_not_assumed_intracellular():
    host, roles = native_host()
    host.add_metabolites([cobra.Metabolite("unknown_marker", formula="H2O", charge=0,
                                         compartment="unclassified")])
    task = chemical_task()
    task = replace(task, metabolites=(replace(task.metabolites[0], compartment="unclassified"),))
    result, _ = predict(host, roles, task)
    assert result["status"] == "failed"
    assert result["upper"] is None
    assert "compartment" in result["reason"] and "role" in result["reason"]


@pytest.mark.parametrize("field,value", [
    ("oxygen_exchange_ids", ()), ("oxygen_exchange_ids", "oxygen_feed"),
    ("oxygen_exchange_ids", ("oxygen_feed", "oxygen_feed")),
    ("intracellular_compartment_ids", ()), ("intracellular_compartment_ids", "c"),
    ("extracellular_compartment_ids", ("c",)),
])
def test_host_roles_require_explicit_unambiguous_boundary_and_compartment_sets(field, value):
    _, roles = native_host()
    with pytest.raises(ValueError):
        replace(roles, **{field: value})


@pytest.mark.parametrize("identifier", ["hydrogen_feed", "assimilate"])
def test_oxygen_boundary_roles_require_single_species_molecular_oxygen(identifier):
    host, roles = native_host()
    roles = replace(roles, oxygen_exchange_ids=(identifier,))
    result, _ = predict(host, roles)
    assert result["status"] == "failed"
    assert result["upper"] is None
    assert "oxygen" in result["reason"]


def test_final_solves_preserve_existing_bounds_and_full_stoichiometric_balance(monkeypatch):
    import numpy as np
    from cobra.core.solution import get_solution

    from ystwin.analysis import blind_transfer

    host, roles = native_host()
    set_boundary(host.reactions.glucose_feed, 2, 1000)
    set_boundary(host.reactions.oxygen_feed, -2, 0.1)
    original_bounds = {reaction.id: reaction.bounds for reaction in host.reactions}
    validated = blind_transfer._validated_optimize
    checked = []

    def verify(model, stoichiometry, context, **kwargs):
        for identifier, bounds in original_bounds.items():
            expected = (0.1, 0.1) if identifier == roles.biomass_reaction_id else bounds
            assert model.reactions.get_by_id(identifier).bounds == expected
        solution = validated(model, stoichiometry, context, **kwargs)
        assert solution is not None
        fluxes = get_solution(model).fluxes
        np.testing.assert_allclose(stoichiometry @ fluxes.to_numpy(), 0, rtol=0, atol=1e-8)
        glucose = 2 * fluxes["glucose_feed"]
        oxygen = -2 * fluxes["oxygen_feed"]
        assert -1e-8 <= glucose <= 2 + 1e-8
        assert -1e-8 <= oxygen <= 1 + 1e-8
        assert 2 * fluxes["native_synthesis"] >= 0.2 * glucose - 1e-8
        checked.append(model.objective.direction)
        return solution

    monkeypatch.setattr(blind_transfer, "_validated_optimize", verify)
    result, details = predict(host, roles, oxygen_task())
    assert result["status"] == "ok", result
    assert details["scenarios"][0]["flux_upper_mmol_per_base_gdw_h"] == pytest.approx(0.2)
    assert checked == ["min", "max"]
    assert {reaction.id: reaction.bounds for reaction in host.reactions} == original_bounds
