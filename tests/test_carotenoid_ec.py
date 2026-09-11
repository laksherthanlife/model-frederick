"""The pathway must install on the enzyme-constrained model too.

GECKO exports carry the same metabolites under compartment-suffixed ids
(``s_0189[c]``) and drop charge annotation. The builder resolves ids rather than
assuming one convention, because the ec model is the one that reproduces the
measured phenotype and therefore the one any real bound has to come from.
"""

import pathlib

import pytest

from ystwin.fba.carotenoid import (
    _NEW_METABOLITES,
    PATHWAY_REACTION_IDS,
    PRODUCT_DEMAND_ID,
    add_beta_carotene_pathway,
)
from ystwin.fba.fva import product_flux_range

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def _ec_carotenoid_template(ec_yeast_gem_factory):
    return add_beta_carotene_pathway(ec_yeast_gem_factory())


@pytest.fixture
def ec_carotenoid(_ec_carotenoid_template, model_copy):
    return model_copy(_ec_carotenoid_template)


def test_installs_despite_the_compartment_suffixed_id_convention(ec_carotenoid):
    assert ec_carotenoid.metabolites.has_id("betacarotene_c")
    assert ec_carotenoid.reactions.has_id("CRTI")


def test_reactions_are_atom_balanced_on_the_ec_model(ec_carotenoid):
    for rid in PATHWAY_REACTION_IDS:
        rxn = ec_carotenoid.reactions.get_by_id(rid)
        if rxn.boundary:
            continue
        imbalance = {k: v for k, v in rxn.check_mass_balance().items() if k != "charge"}
        assert imbalance == {}, f"{rid} unbalanced: {imbalance}"


def test_it_draws_from_the_same_native_ggpp_pool(ec_carotenoid):
    psy = ec_carotenoid.reactions.get_by_id("CRTYB_PSY")

    assert any(m.id.startswith("s_0189") for m in psy.reactants)


def test_the_product_range_is_still_undetermined_under_enzyme_constraints(ec_carotenoid):
    """The ec model fixes the phenotype but not the missing pathway capacity."""
    rng = product_flux_range(
        ec_carotenoid, PRODUCT_DEMAND_ID, growth_fraction=0.90, glucose_uptake=None,
    )

    assert rng.minimum == pytest.approx(0.0, abs=1e-9)
    assert rng.relative_width > 0.99


def test_growth_is_not_inflated_by_installing_a_drain(ec_yeast_gem, ec_carotenoid):
    """A sink can only cost carbon, never create it. LP tolerance is ~1e-6 here."""
    base = ec_yeast_gem.slim_optimize()

    assert ec_carotenoid.slim_optimize() <= base * (1 + 1e-4)


def test_a_model_missing_a_required_precursor_is_refused_by_name(yeast_gem):
    with pytest.raises(KeyError, match="s_0189"):
        add_beta_carotene_pathway(yeast_gem, ggpp_id="s_0189_nonexistent")


class TestTheDeclaredIdentityIsTheSameInBothPlaces:
    """Which species `phytoene_c` IS, is written down twice, and twice is one too many.

    `_NEW_METABOLITES` here carries the KEGG and ModelSEED ids as cobra annotations, which
    is where every other consumer of a model looks -- `bridge/thermodynamic.py` and
    `bridge/equilibrator.py` both resolve energies through `kegg.compound` and neither
    knows this module exists. `data/thermo/heterologous_metabolites.tsv` carries the same
    pairs as data, which is what lets the gate resolve them for a model that does NOT have
    this pathway installed -- the case `predict.py` is in, gating from a declared
    stoichiometry with no model at all.

    Both are needed and a duplicated fact drifts, so this is the cross-check rather than a
    third copy. It is the same move `thermo_gate._check_pair_agrees` makes on the
    stoichiometry the spec declares beside the reaction id.
    """

    @pytest.fixture(scope="class")
    def declared(self):
        import csv

        path = (pathlib.Path(__file__).resolve().parents[1] / "data" / "thermo"
                / "heterologous_metabolites.tsv")
        with path.open(encoding="utf-8") as handle:
            rows = [line for line in handle if not line.startswith("#")]
        return {row["metabolite_id"]: row
                for row in csv.DictReader(rows, delimiter="\t")}

    def test_the_two_tables_name_the_same_three_species(self, declared):
        assert set(declared) == {mid for mid, *_ in _NEW_METABOLITES}

    def test_every_identifier_agrees_row_for_row(self, declared):
        from_module = {mid: (kegg, seed, formula)
                       for mid, _name, formula, _charge, kegg, seed in _NEW_METABOLITES}
        from_data = {mid: (row["kegg_id"], row["modelseed_id"], row["formula"])
                     for mid, row in declared.items()}

        assert from_module == from_data

    def test_every_declared_row_says_where_its_identifier_came_from(self, declared):
        """An identifier with no provenance is the thing this whole track was fixing."""
        for metabolite_id, row in declared.items():
            assert "cpd" in row["source"] and "KEGG" in row["source"], metabolite_id

    def test_the_installed_model_carries_the_annotation_the_energy_layer_reads(
            self, ec_carotenoid):
        """Not the declared table -- the model object. `bridge/equilibrator.py` resolves
        through `kegg.compound` and has no way to reach a TSV under data/thermo."""
        for metabolite_id, _name, _formula, _charge, kegg, seed in _NEW_METABOLITES:
            annotation = ec_carotenoid.metabolites.get_by_id(metabolite_id).annotation
            assert annotation["kegg.compound"] == kegg
            assert annotation["seed.compound"] == seed


@pytest.fixture(scope="module")
def _installed_template(yeast_gem_factory):
    """Plain Yeast9 with the pathway on it, for the chemistry rather than the EC bounds."""
    return add_beta_carotene_pathway(yeast_gem_factory())


@pytest.fixture
def installed(_installed_template, model_copy):
    return model_copy(_installed_template)

class TestCrtIHasTheRightElectronAcceptor:
    """It was written with four free FAD, which is impossible and cost the model its O2 need.

    The shipped form scored +166.2 +- 12.8 kJ/mol -- uphill by thirteen standard errors on a
    step Verwaal 2007 (PMID 17496128) measured running to completion. FAD in a bacterial-type
    phytoene desaturase is a prosthetic group that gets reoxidised, not a substrate consumed
    four times per turnover, and free FAD is a poor oxidant at E'0 = -0.219 V.
    """

    def test_it_consumes_oxygen(self, installed):
        """Carotenoid synthesis in yeast is strictly aerobic. The model had no O2 term."""
        reaction = installed.reactions.get_by_id("CRTI")
        consumed = {m.id for m, c in reaction.metabolites.items() if c < 0}

        assert "s_1275" in consumed

    def test_it_makes_hydrogen_peroxide_not_water(self, installed):
        """A flavin oxidase reduces O2 by two electrons at one site. The four-electron form
        is more favourable on paper (-676 against -293) and is not the chemistry."""
        reaction = installed.reactions.get_by_id("CRTI")
        produced = {m.id for m, c in reaction.metabolites.items() if c > 0}

        assert "s_0837" in produced
        assert "s_0803" not in produced

    def test_no_free_flavin_is_consumed(self, installed):
        """The defect, stated as the thing that must not come back."""
        reaction = installed.reactions.get_by_id("CRTI")
        touched = {m.id for m in reaction.metabolites}

        assert "s_0687" not in touched and "s_0689" not in touched

    def test_it_is_mass_and_charge_balanced(self, installed):
        """phytoene C40H64 + 4 O2 -> lycopene C40H56 + 4 H2O2. Eight H move, eight O in."""
        reaction = installed.reactions.get_by_id("CRTI")

        assert reaction.check_mass_balance() == {}

    def test_without_oxygen_no_beta_carotene_can_be_made(self, installed):
        """The consequence worth having. Before this the model would have predicted
        carotenoid in an anaerobic chemostat."""
        with installed as model:
            model.reactions.get_by_id("r_1992").lower_bound = 0.0
            model.objective = "DM_betacarotene_c"
            flux = model.slim_optimize()

        assert flux is None or abs(flux) < 1e-9

    def test_and_the_ceiling_rises_with_oxygen_until_it_saturates(self, installed):
        """A mechanism for the environment to reach the product, which the kinetic layer
        still does not use -- see docs/FINDINGS.md."""
        got = []
        for bound in (0.1, 1.0, 8.0):
            with installed as model:
                model.reactions.get_by_id("r_1992").lower_bound = -bound
                model.objective = "DM_betacarotene_c"
                got.append(model.slim_optimize() or 0.0)

        assert got[0] < got[1] < got[2] or got[1] == pytest.approx(got[2], rel=1e-6)
        assert got[0] < got[2]
