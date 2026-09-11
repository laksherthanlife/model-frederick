from __future__ import annotations

import importlib.util
import json
from dataclasses import fields

import cobra
import pytest

from ystwin import paths
from ystwin.fba.stress_pathways import PathwayItem, load_coverage, model_gene_coverage
from ystwin.generator.in_silico import INPUT_NAMES, LATENT_NAMES, TeacherParameters
from ystwin.generator.stress_panel import MODULES, REPORTERS


def item(key, orfs, layer="CASCADE"):
    return PathwayItem(layer, key, key, orfs, "NOT_METABOLIC", "-", "-", "-",
                       "Curated signaling identity, not a claim of dynamic coverage", "test fixture")


def test_gene_association_audit_uses_identifiers_not_display_names():
    model = cobra.Model("identity_fixture")
    metabolite = cobra.Metabolite("m", compartment="c")
    reaction = cobra.Reaction("metabolic_task")
    reaction.add_metabolites({metabolite: -1})
    reaction.gene_reaction_rule = "Y_REAL"
    model.add_reactions([reaction])
    model.genes.Y_REAL.name = "HOG1"

    result = model_gene_coverage(model, (item("hog", "YLR113W"),))[0]
    assert result.present_orfs == ()
    assert result.absent_orfs == ("YLR113W",)
    assert result.reaction_ids == ()
    assert result.evidence_scope == "gene_reaction_associations_only"


def test_associated_kinase_does_not_imply_a_signaling_model():
    model = cobra.Model("association_fixture")
    reaction = cobra.Reaction("generic_ATP_reaction")
    reaction.add_metabolites({cobra.Metabolite("ATP", compartment="c"): -1})
    reaction.gene_reaction_rule = "YLR113W"
    model.add_reactions([reaction])
    original = [(r.id, r.bounds, r.gene_reaction_rule, dict(r.metabolites)) for r in model.reactions]

    result = model_gene_coverage(model, (item("hog", "YLR113W|YJL128C"),))[0]
    assert result.model_id == model.id
    assert result.present_orfs == ("YLR113W",)
    assert result.absent_orfs == ("YJL128C",)
    assert result.reaction_ids == ("generic_ATP_reaction",)
    assert result.associations == (("YLR113W", ("generic_ATP_reaction",)),)
    assert result.evidence_scope == "gene_reaction_associations_only"
    assert [(r.id, r.bounds, r.gene_reaction_rule, dict(r.metabolites)) for r in model.reactions] == original


def test_coverage_is_queried_per_model_not_copied_from_reference_verdict():
    first, second = cobra.Model("first"), cobra.Model("second")
    reaction = cobra.Reaction("renamed_enzyme_reaction")
    reaction.gene_reaction_rule = "Y_GENE"
    first.add_reactions([reaction])
    declared = PathwayItem("ADAPTATION", "effector", "effector", "Y_GENE", "IN_GEM",
                           "reference_only_reaction", "-", "-", "reference annotation", "fixture")

    assert model_gene_coverage(first, (declared,))[0].reaction_ids == ("renamed_enzyme_reaction",)
    assert model_gene_coverage(second, (declared,))[0].reaction_ids == ()
    assert model_gene_coverage(second, (declared,))[0].absent_orfs == ("Y_GENE",)


def test_gene_free_items_and_duplicate_items_are_explicit():
    model = cobra.Model("empty")
    result = model_gene_coverage(model, (item("geometry", "-", "ADAPTATION"),))[0]
    assert result.declared_orfs == ()
    assert result.present_orfs == result.absent_orfs == result.reaction_ids == ()
    with pytest.raises(ValueError, match="duplicate"):
        model_gene_coverage(model, (item("same", "Y_A"), item("same", "Y_B")))


@pytest.fixture
def auditor():
    spec = importlib.util.spec_from_file_location("stress_map_auditor", paths.REPO_ROOT / "scripts/audit_stress_map.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_map_audit_does_not_infer_dynamics_or_optimize_metabolism(auditor, monkeypatch):
    mapping = json.loads(auditor.MAP_PATH.read_text())
    model = cobra.Model("empty_reference")

    def forbidden(*args, **kwargs):
        raise AssertionError("a gene-association audit must not perform metabolic optimization")

    monkeypatch.setattr(model, "slim_optimize", forbidden)
    summary, associations, branches = auditor.build_report(
        {"ec_predictor": model}, mapping, load_coverage())
    assert summary["current_student"]["imposed_inputs"] == list(INPUT_NAMES)
    assert summary["current_student"]["latent_coordinates"] == list(LATENT_NAMES)
    assert summary["current_student"]["independent_biological_validation"] is False
    assert summary["current_student"]["shared_metabolic_backend_in_reference_comparison"] is True
    assert summary["models"]["ec_predictor"]["items_with_gpr_associations"] == 0
    assert set(associations.evidence_scope) == {"gene_reaction_associations_only"}
    assert len(branches) == len(mapping["branches"])
    assert "signaling_covered" not in branches.columns


def test_every_teacher_constant_is_classified_as_an_input_assumption(auditor):
    inventory = auditor.teacher_assumptions()
    assert set(inventory.parameter) == {field.name for field in fields(TeacherParameters)}
    coefficients = inventory[inventory.role == "frozen_synthetic_control_law"]
    assert len(coefficients) > 0
    assert not coefficients.constructor_input.any()
    assert set(inventory.source) == {"declared reference-teacher assumption, not fitted biological evidence"}


def test_expanded_map_rejects_unresolved_evidence_references(auditor):
    mapping = json.loads(auditor.MAP_PATH.read_text())
    mapping["branches"][0]["coverage_items"].append("unverified.pathway")
    with pytest.raises(ValueError, match="unknown coverage item"):
        auditor.validate_map(mapping, load_coverage())


def test_atp_catalogue_matches_queen_measurand_and_yeast_source():
    assert REPORTERS["QUEEN-2m"].module == "atp"
    assert MODULES["atp"].transcription_factor == "cytosolic ATP"
    assert "30858198" in MODULES["atp"].source
