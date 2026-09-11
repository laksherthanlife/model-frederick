"""The oxidative cost of secretory folding, which yeast-GEM v9.0.2 does not carry.

Forming a disulfide bond in the ER consumes molecular oxygen and produces hydrogen peroxide.
The vendored model has no reaction that does it, no `ERO1` or `PDI1` gene, and no peroxide
species in the `er` compartment at all -- so the burden that makes a disulfide-rich secreted
protein expensive was absent rather than unparameterised.

That matters for this project specifically: `NativeYap1` and `AlteredYap1` were characterised
over an H2O2 ladder and `UPRE1` and `UPRE2` over a DTT ladder, and a secreted disulfide-rich
protein loads both. See `docs/research/SECRETED_PROTEIN.md`.
"""

from __future__ import annotations

import cobra
import pytest

from ystwin import paths
from ystwin.fba.secretion import (
    DISULFIDE_O2_PER_BOND,
    ER_PEROXIDE_ID,
    PEROXIDE_EXPORT_ID,
    add_oxidative_folding,
    folding_reaction,
)

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def base():
    if paths.yeast_gem() is None:
        pytest.skip("needs yeast-GEM; see docs/REPRODUCING.md")
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cobra.io.read_sbml_model(str(paths.yeast_gem()))


@pytest.fixture(scope="module")
def installed(base):
    return add_oxidative_folding(base)


def _protein(model, formula="C10H20N2O5S6"):
    metabolite = cobra.Metabolite("testprot_er", name="test protein", formula=formula,
                                  charge=0, compartment="er")
    model.add_metabolites([metabolite])
    return metabolite


class TestWhatTheVendoredModelWasMissing:
    """Pinned so that a future GEM version adding these is noticed rather than duplicated."""

    def test_the_stock_model_has_no_er_peroxide(self, base):
        assert not base.metabolites.has_id(ER_PEROXIDE_ID)
        peroxide = [m.compartment for m in base.metabolites if m.formula == "H2O2"]

        assert "er" not in peroxide
        assert set(peroxide) == {"c", "m", "n", "p"}

    def test_and_no_reaction_that_forms_a_protein_disulfide(self, base):
        """Four name matches and all four are cofactor recycling, not substrate oxidation."""
        named = [r.id for r in base.reactions
                 if r.name and "disulfide" in r.name.lower()]

        assert set(named) <= {"r_1806", "r_4186"}

    def test_and_no_ero1_or_pdi1_gene(self, base):
        names = {(g.name or "").upper() for g in base.genes}

        assert "ERO1" not in names and "PDI1" not in names


class TestInstallingIt:
    def test_it_adds_er_peroxide(self, installed):
        peroxide = installed.metabolites.get_by_id(ER_PEROXIDE_ID)

        assert peroxide.compartment == "er"
        assert peroxide.formula == "H2O2"

    def test_it_exports_to_the_cytosol_so_the_pool_is_not_a_dead_end(self, installed):
        export = installed.reactions.get_by_id(PEROXIDE_EXPORT_ID)
        moved = {m.id: c for m, c in export.metabolites.items()}

        assert moved == {ER_PEROXIDE_ID: -1.0, "s_0837": 1.0}

    def test_the_cytosol_can_actually_dispose_of_it(self, installed):
        """Catalase. Without a sink the folding reaction cannot carry flux at all."""
        catalase = installed.reactions.get_by_id("r_0255")

        assert installed.metabolites.get_by_id("s_0837") in catalase.metabolites

    def test_the_input_model_is_not_modified(self, base, installed):
        assert not base.metabolites.has_id(ER_PEROXIDE_ID)
        assert installed is not base

    def test_installing_twice_is_refused(self, installed):
        with pytest.raises(ValueError, match="already installed"):
            add_oxidative_folding(installed)

    def test_a_model_without_er_oxygen_is_refused_by_name(self, base):
        stripped = base.copy()
        stripped.remove_metabolites([stripped.metabolites.get_by_id("s_1276")])

        with pytest.raises(ValueError, match="s_1276"):
            add_oxidative_folding(stripped)


class TestTheCostItself:
    def test_one_oxygen_and_one_peroxide_per_disulfide(self, installed):
        model = installed.copy()
        _protein(model)
        reaction = folding_reaction(model, "testprot_er", 3)
        moved = {m.id: c for m, c in reaction.metabolites.items()}

        assert moved["s_1276"] == -3.0 * DISULFIDE_O2_PER_BOND
        assert moved[ER_PEROXIDE_ID] == 3.0

    def test_the_folded_protein_is_two_hydrogens_lighter_per_bond(self, installed):
        """`2 R-SH + O2 -> R-S-S-R + H2O2` moves two hydrogens off the protein per bond.
        Giving the folded species the same formula leaves the reaction unbalanced by
        exactly 2n hydrogen, which is what the first version of this did."""
        model = installed.copy()
        _protein(model, formula="C10H20N2O5S6")
        folding_reaction(model, "testprot_er", 3)

        assert model.metabolites.get_by_id("testprot_er_folded").formula == "C10H14N2O5S6"

    def test_and_the_reaction_balances(self, installed):
        model = installed.copy()
        _protein(model)
        reaction = folding_reaction(model, "testprot_er", 3)

        assert reaction.check_mass_balance() == {}

    def test_zero_disulfides_is_legal_and_free(self, installed):
        """Not every secreted protein has one, and refusing zero would make every caller
        special-case it."""
        model = installed.copy()
        _protein(model)
        reaction = folding_reaction(model, "testprot_er", 0)
        moved = {m.id: c for m, c in reaction.metabolites.items()}

        assert "s_1276" not in moved and ER_PEROXIDE_ID not in moved
        assert reaction.check_mass_balance() == {}

    def test_a_negative_count_is_refused(self, installed):
        model = installed.copy()
        _protein(model)

        with pytest.raises(ValueError, match="disulfides must be"):
            folding_reaction(model, "testprot_er", -1)

    def test_a_formula_with_too_few_hydrogens_is_refused(self, installed):
        """Either the formula or the count is wrong, and guessing which would be worse."""
        model = installed.copy()
        _protein(model, formula="C10H2N2O5S6")

        with pytest.raises(ValueError, match="fewer than"):
            folding_reaction(model, "testprot_er", 3)

    def test_folding_without_installing_first_is_refused(self, base):
        model = base.copy()
        _protein(model)

        with pytest.raises(ValueError, match="add_oxidative_folding"):
            folding_reaction(model, "testprot_er", 3)


class TestItCanActuallyCarryFlux:
    def test_the_folding_reaction_reaches_a_sink(self, installed):
        """The whole point of the export: a folded protein with nowhere to put its peroxide
        is a reaction that optimises to zero and looks like a capacity limit."""
        model = installed.copy()
        _protein(model)
        folding_reaction(model, "testprot_er", 3)
        model.add_boundary(model.metabolites.get_by_id("testprot_er"), type="sink")
        model.add_boundary(model.metabolites.get_by_id("testprot_er_folded"), type="demand")
        model.objective = "DM_testprot_er_folded"

        assert (model.slim_optimize() or 0.0) > 1e-9

    def test_and_stops_without_oxygen(self, installed):
        """Disulfide formation is aerobic. If this passes anaerobically the cost is not
        actually being paid."""
        model = installed.copy()
        _protein(model)
        folding_reaction(model, "testprot_er", 3)
        model.add_boundary(model.metabolites.get_by_id("testprot_er"), type="sink")
        model.add_boundary(model.metabolites.get_by_id("testprot_er_folded"), type="demand")
        model.objective = "DM_testprot_er_folded"
        model.reactions.get_by_id("r_1992").lower_bound = 0.0

        assert (model.slim_optimize() or 0.0) < 1e-9
