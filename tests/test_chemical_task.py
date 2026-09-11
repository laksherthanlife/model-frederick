from __future__ import annotations

from dataclasses import replace

import cobra
import pytest

from ystwin.fba.chemical_task import (
    ChemicalTask,
    EnzymeDefinition,
    MetaboliteDefinition,
    ReactionDefinition,
    install_chemical_task,
)


def host():
    model = cobra.Model("native_host")
    model.solver = "glpk"
    substrate = cobra.Metabolite("native_feed", formula="C6H12O6", charge=0, compartment="c")
    pool = cobra.Metabolite("native_protein_pool", compartment="c")
    feed = cobra.Reaction("feed", lower_bound=0, upper_bound=2)
    feed.add_metabolites({substrate: 1})
    protein = cobra.Reaction("protein_supply", lower_bound=0, upper_bound=0.01)
    protein.add_metabolites({pool: 1})
    model.add_reactions([feed, protein])
    return model


def task(label="opaque_a", *, enzyme=False):
    return ChemicalTask(
        metabolites=(MetaboliteDefinition(label, "C6H12O6", 0, "c"),),
        reactions=(ReactionDefinition(f"reaction_{label}", (("native_feed", -1), (label, 1)),
                                      "independent chemical definition",
                                      "catalyst_a" if enzyme else None, 1.0 if enzyme else None),),
        output_metabolite_id=label,
        enzymes=(EnzymeDefinition("catalyst_a", 10000, 0.0001, "independent enzyme assay"),) if enzyme else (),
        protein_pool_metabolite_id="native_protein_pool" if enzyme else None,
    )


def capacity(model, installed):
    model.objective = installed.reaction_id
    return model.slim_optimize(error_value=None)


def test_generic_chemical_task_consumes_real_substrate_without_mutating_host():
    original = host()
    before = tuple(original.reactions.list_attr("id"))
    model, installed = install_chemical_task(original, task())
    assert capacity(model, installed) == pytest.approx(2)
    assert tuple(original.reactions.list_attr("id")) == before
    assert installed.carbon_atoms == 6
    assert installed.metadata["prediction_scope"] == "feasible_capacity_not_realized_allocation"
    model.reactions.feed.upper_bound = 0
    assert capacity(model, installed) == pytest.approx(0)


def test_renaming_unseen_task_has_no_effect_on_predicted_capacity():
    first, a = install_chemical_task(host(), task("opaque_a"))
    second, b = install_chemical_task(host(), task("entirely_different_identifier"))
    assert capacity(first, a) == pytest.approx(capacity(second, b))


def test_independent_enzyme_inputs_limit_catalysis_inside_shared_pool():
    model, installed = install_chemical_task(host(), task(enzyme=True))
    assert capacity(model, installed) == pytest.approx(0.36)
    model.reactions.protein_supply.upper_bound = 0.0005
    assert capacity(model, installed) == pytest.approx(0.18)


def test_unbalanced_or_unidentified_chemical_participants_are_refused():
    spec = task()
    bad = replace(spec, metabolites=(replace(spec.metabolites[0], formula="C5H10O5"),))
    with pytest.raises(ValueError, match="balance"):
        install_chemical_task(host(), bad)
    reaction = replace(spec.reactions[0], stoichiometry=(("native_protein_pool", -1), ("opaque_a", 1)))
    with pytest.raises(ValueError, match="formula"):
        install_chemical_task(host(), replace(spec, reactions=(reaction,)))


def test_duplicate_reactions_and_missing_enzyme_evidence_are_refused():
    spec = task()
    with pytest.raises(ValueError, match="duplicate"):
        install_chemical_task(host(), replace(spec, reactions=spec.reactions * 2))
    constrained = task(enzyme=True)
    enzyme = replace(constrained.enzymes[0], source="")
    with pytest.raises(ValueError, match="source"):
        install_chemical_task(host(), replace(constrained, enzymes=(enzyme,)))


def test_input_schema_contains_no_product_fit_or_outcome_channel():
    spec = task()
    payload = spec.to_dict()
    assert ChemicalTask.from_dict(payload) == spec
    with pytest.raises(ValueError, match="unknown"):
        ChemicalTask.from_dict({**payload, "measured_titre": 1})
    with pytest.raises(ValueError, match="unknown"):
        ChemicalTask.from_dict({**payload, "allocation_fraction": 0.7})
