"""The heterologous pathway is built from the actual construct genes, balanced.

crtE / crtYB / crtI from Xanthophyllomyces dendrorhous plus truncated HMG1.
crtYB is bifunctional: phytoene synthase and lycopene cyclase are separate
reactions carrying the same gene.
"""

import pytest

pytestmark = pytest.mark.integration


def test_every_added_reaction_is_mass_and_charge_balanced(carotenoid_model):
    from ystwin.fba.carotenoid import PATHWAY_REACTION_IDS

    for rid in PATHWAY_REACTION_IDS:
        rxn = carotenoid_model.reactions.get_by_id(rid)
        if rxn.boundary:
            continue
        assert rxn.check_mass_balance() == {}, f"{rid} unbalanced: {rxn.check_mass_balance()}"


def test_the_pathway_draws_carbon_from_the_native_ggpp_pool(carotenoid_model):
    psy = carotenoid_model.reactions.get_by_id("CRTYB_PSY")

    assert carotenoid_model.metabolites.get_by_id("s_0189") in psy.reactants


def test_beta_carotene_can_carry_flux_when_it_is_the_objective(carotenoid_model):
    with carotenoid_model as m:
        m.reactions.get_by_id("r_1714").lower_bound = -10.0
        m.objective = "DM_betacarotene_c"
        assert m.slim_optimize() > 1e-6


def test_adding_a_drain_pathway_does_not_raise_maximum_growth(yeast_gem, carotenoid_model):
    base = yeast_gem.slim_optimize()

    assert carotenoid_model.slim_optimize() == pytest.approx(base, rel=1e-6)


def test_reactions_carry_the_real_construct_gene_names_not_generic_labels(carotenoid_model):
    rules = {
        rid: carotenoid_model.reactions.get_by_id(rid).gene_reaction_rule
        for rid in ("CRTE", "CRTYB_PSY", "CRTI", "CRTYB_LCY")
    }

    assert rules["CRTE"] == "crtE"
    assert rules["CRTI"] == "crtI"
    assert rules["CRTYB_PSY"] == "crtYB"
    assert rules["CRTYB_LCY"] == "crtYB"


def test_desaturation_is_oxidative_so_the_pathway_has_a_redox_cost(carotenoid_model):
    """Four desaturations take eight hydrogens off phytoene; something must accept them.

    This asserted FAD and FADH2 by name until 2026-08-30. The claim it was making -- that the
    step is oxidative and the pathway pays a redox cost -- is unchanged and still checked
    here; what changed is which species carries it. CrtI's flavin is a prosthetic group that
    is reoxidised, not a substrate consumed four times per turnover, and writing four free
    FADH2 as products made the reaction come out at +166.2 kJ/mol -- impossible, on a step
    Verwaal 2007 (PMID 17496128) measured running to completion. The acceptor is O2 and the
    two-electron product is hydrogen peroxide.
    """
    crti = carotenoid_model.reactions.get_by_id("CRTI")
    names = {m.name for m in crti.reactants} | {m.name for m in crti.products}

    assert "oxygen" in names and "hydrogen peroxide" in names
    assert not {"FAD", "FADH2"} & names


def test_the_redox_cost_is_paid_per_desaturation(carotenoid_model):
    """Eight hydrogens, two per O2, so four O2 in and four H2O2 out."""
    crti = carotenoid_model.reactions.get_by_id("CRTI")
    by_name = {m.name: c for m, c in crti.metabolites.items()}

    assert by_name["oxygen"] == -4.0
    assert by_name["hydrogen peroxide"] == 4.0


def test_intermediates_are_present_so_they_can_serve_as_independent_anchors(carotenoid_model):
    """Phytoene and lycopene come off the same HPLC run as beta-carotene."""
    for mid in ("phytoene_c", "lycopene_c", "betacarotene_c"):
        assert carotenoid_model.metabolites.has_id(mid)


def test_installing_the_pathway_twice_is_refused(carotenoid_model):
    from ystwin.fba.carotenoid import add_beta_carotene_pathway

    with pytest.raises(ValueError, match="already"):
        add_beta_carotene_pathway(carotenoid_model)


def test_the_source_model_is_left_untouched(yeast_gem, carotenoid_model):
    assert not yeast_gem.metabolites.has_id("betacarotene_c")
