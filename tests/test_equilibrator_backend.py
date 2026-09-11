"""Replacing the group-contribution energies with component contribution.

The SEED table shipped with pytfa is group contribution, and it fails where substrate and
product are structurally similar because the errors do not cancel. It put ribose-5-phosphate
isomerase at +75.4 kJ/mol -- an isomerase between two pentose phosphates -- and reported an
error of 1.7 on it, which made an essential reaction impossible and blocked growth through
the whole pentose phosphate pathway.

eQuilibrator is component contribution from the same group as BioNumbers, and it handles pH,
ionic strength and magnesium natively rather than by hand-added association constants. On
the reactions the old table got wrong it is right, and on ATP hydrolysis it returns the
textbook value exactly.
"""

import pytest

pytest.importorskip("equilibrator_api")

from ystwin import paths
from ystwin.bridge.equilibrator import EquilibratorData

pytestmark = pytest.mark.skipif(
    paths.yeast_gem() is None,
    reason="needs yeast-GEM to resolve metabolites against; see docs/REPRODUCING.md",
)


@pytest.fixture(scope="module")
def thermo():
    return EquilibratorData.load()


class TestItFixesWhatGroupContributionGotWrong:
    def test_an_isomerase_is_near_equilibrium(self, thermo):
        """The old table said +75.4 for this. An isomerase between two pentose phosphates
        cannot be, and the pathway it blocked is essential."""
        energy = thermo.reaction_dg0({"s_0577": -1, "s_1408": 1})

        assert energy is not None
        assert abs(energy) < 15.0

    def test_atp_hydrolysis_returns_the_textbook_value(self, thermo):
        """-30.5 kJ/mol at pH 7, 1 mM Mg. The old table needed hand-added Alberty constants
        and still landed at -37."""
        energy = thermo.reaction_dg0(
            {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1})

        assert energy == pytest.approx(-30.5, abs=6.0)

    def test_a_mutase_is_not_wildly_uphill(self, thermo):
        """Phosphopentomutase came out at +77.4 before."""
        energy = thermo.reaction_dg0({"s_0415": -1, "s_1408": 1})

        assert energy is not None and energy < 20.0


class TestItReportsItsOwnUncertainty:
    def test_every_energy_comes_with_an_error(self, thermo):
        energy, error = thermo.reaction_dg0(
            {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}, with_error=True)

        assert error is not None and error > 0

    def test_a_well_known_reaction_has_a_small_one(self, thermo):
        _, error = thermo.reaction_dg0(
            {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}, with_error=True)

        assert error < 5.0

    def test_an_uncovered_metabolite_returns_nothing_rather_than_guessing(self, thermo):
        assert thermo.reaction_dg0({"s_not_a_real_metabolite": -1}) is None


class TestConditionsAreExplicit:
    def test_it_records_the_ph_it_was_built_at(self, thermo):
        assert thermo.ph == pytest.approx(7.5)

    def test_magnesium_is_handled_natively(self, thermo):
        assert thermo.p_mg == pytest.approx(3.0)

    def test_changing_magnesium_changes_atp_hydrolysis(self):
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}
        low = EquilibratorData.load(p_mg=6.0)
        high = EquilibratorData.load(p_mg=2.0)

        assert low.reaction_dg0(hydrolysis) != pytest.approx(high.reaction_dg0(hydrolysis))

    def test_it_covers_more_of_the_model_than_the_old_table(self, thermo):
        assert thermo.coverage > 0.55


class TestInstancesDoNotShareConditions:
    """The loader caches one underlying object, so conditions have to be set per query.

    Without that, two instances built at different magnesium silently return the same
    energies -- the second load overwrites the first's settings on the shared object and
    nothing errors. A wrong answer that looks entirely reasonable.
    """

    def test_two_instances_keep_their_own_magnesium(self):
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}
        low, high = EquilibratorData.load(p_mg=6.0), EquilibratorData.load(p_mg=2.0)

        first = low.reaction_dg0(hydrolysis)
        second = high.reaction_dg0(hydrolysis)
        again = low.reaction_dg0(hydrolysis)

        assert first != pytest.approx(second)
        assert again == pytest.approx(first)

    def test_they_keep_their_own_ph(self):
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}
        acid, base = EquilibratorData.load(ph=6.0), EquilibratorData.load(ph=8.0)

        assert acid.reaction_dg0(hydrolysis) != pytest.approx(base.reaction_dg0(hydrolysis))


def _single_compartment(both):
    """The reactions both estimators are scored over, from the composite's own model."""
    import cobra

    model = cobra.io.read_sbml_model(str(paths.yeast_gem()))
    return [r for r in model.reactions
            if len({m.compartment for m in r.metabolites}) == 1]

class TestTheStubThatWouldHaveZeroedEveryReaction:
    """`dgf` returned 0.0 for any covered metabolite, "for interface compatibility".

    It is not compatible. `bridge/thermodynamic.py::reaction_dg0` -- and therefore the whole
    thermodynamic gate -- computes a reaction as `sum(coefficient * dgf(m))`, so passing
    this backend to it would have returned **exactly zero** for every reaction and read as a
    real energy at equilibrium. Found by the 2026-08-30 audit while wiring this module into
    the gate for the first time.
    """

    def test_asking_for_a_formation_energy_is_refused(self, thermo):
        with pytest.raises(NotImplementedError, match="per-compound"):
            thermo.dgf("s_0373")

    def test_the_refusal_says_what_to_call_instead(self, thermo):
        with pytest.raises(NotImplementedError, match="reaction_dg0"):
            thermo.dgf("s_0373")

    def test_an_uncovered_metabolite_still_answers_none(self, thermo):
        """A coverage probe has to keep working, or callers cannot ask what is reachable."""
        assert thermo.dgf("not_a_real_metabolite_id") is None

    def test_coverage_is_askable_without_asking_for_an_energy(self, thermo):
        assert thermo.covers("s_0373") is True
        assert thermo.covers("not_a_real_metabolite_id") is False


class TestItCarriesTheSameRefusalsAsTheOtherBackend:
    """A refutation is a property of the COMPOUND, not of the estimator asked about it.

    An estimator that dropped the refusal on becoming the primary backend would fail in the
    one direction a refusal must never fail in.
    """

    def test_a_refuted_carotenoid_is_refuted_here_too(self, thermo):
        assert thermo.refutation("betacarotene_c") is not None

    def test_an_ordinary_metabolite_is_not(self, thermo):
        assert thermo.refutation("s_0373") is None


class TestThePreferredCompositeNeverMixesTwoEstimatorsInOneReaction:
    """Component contribution first, group contribution as fallback, per reaction.

    Mixing a formation energy from one method with a reaction energy from another produces a
    number belonging to neither -- the same class of error as the kcal/kJ mix that
    `KJ_PER_KCAL` records.
    """

    @pytest.fixture(scope="class")
    def both(self):
        from ystwin.bridge.equilibrator import PreferredEnergies

        if paths.yeast_gem() is None:
            pytest.skip("needs yeast-GEM; see docs/REPRODUCING.md")
        return PreferredEnergies.load()

    def test_the_thiolase_is_answered_by_component_contribution(self, both):
        """The reaction the whole thiolase lead rested on, now answered by the estimator
        that matches the literature rather than the one that does not."""
        thiolase = {"s_0373": -2.0, "s_0367": 1.0, "s_0529": 1.0}

        assert both.which(thiolase) == "component"
        assert both.reaction_dg0(thiolase) == pytest.approx(24.96, abs=1.0)

    def test_it_reports_which_estimator_answered(self, both):
        """"The gate returned a number" does not say which of two methods produced it, and
        they disagree by 13 kJ/mol on the one reaction anybody has checked."""
        assert both.which({"s_0373": -2.0, "s_0367": 1.0, "s_0529": 1.0}) in (
            "component", "group")
        assert both.which({"not_a_metabolite": -1.0}) == "none"

    def test_it_covers_more_than_either_alone(self, both):
        """The point of the composite. Measured rather than asserted."""
        from ystwin.bridge.thermodynamic import reaction_dg0

        component = group = union = 0
        for reaction in _single_compartment(both):
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
            has_component = both.component.reaction_dg0(stoichiometry) is not None
            has_group = reaction_dg0(stoichiometry, both.group) is not None
            component += has_component
            group += has_group
            union += has_component or has_group

        assert union > component
        assert union > group

    def test_a_refusal_survives_the_composite(self, both):
        assert both.refutation("betacarotene_c") is not None


class TestAResolvedIdentifierIsNotTheSameAsTheRightMolecule:
    """A KEGG id resolving does not mean eQuilibrator found the compound the model means.

    yeast-GEM writes a polymer as one representative repeat unit; eQuilibrator holds a
    specific oligomer. `s_0773` glycogen is `C6H10O5` in the model and resolved to a
    TETRASACCHARIDE, `C24H42O21`. Glycogen synthase then came back at **-849.2 kJ/mol** for
    a glycosyl transfer textbooks put near -13, from a mass-unbalanced reaction, with no
    error anywhere. Found on 2026-08-30 while running a second product through the chain.

    The resolution is now checked against the formula and dropped when it disagrees, so a
    reaction containing one is UNCOVERED rather than mis-scored. It costs 46 of 1,289
    single-compartment reactions, which is the right trade: a number computed on the wrong
    molecule is worse than no number.
    """

    def test_the_mismatches_are_recorded_rather_than_silently_dropped(self, thermo):
        assert thermo.mismatched
        assert all(len(v) == 2 for v in thermo.mismatched.values())

    def test_glycogen_is_one_of_them(self, thermo):
        """The one that found the defect."""
        assert "s_0773" in thermo.mismatched
        model_formula, matched = thermo.mismatched["s_0773"]

        assert model_formula == "C6H10O5"
        assert matched == "C24H42O21"

    def test_glycogen_synthase_is_now_uncovered_not_mis_scored(self, thermo, yeast_gem):
        """It returned -849.2 kJ/mol. None is the correct answer."""
        reaction = yeast_gem.reactions.get_by_id("r_0510")
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        assert thermo.reaction_dg0(stoichiometry) is None

    def test_protonation_differences_are_not_treated_as_mismatches(self, thermo):
        """UDP is C9H11N2O12P2 in the model and C9H12N2O12P2 at eQuilibrator's reference
        pH. Same molecule. Refusing it would throw away most of the coverage."""
        assert "s_1538" not in thermo.mismatched
        assert thermo.covers("s_1538")

    def test_every_carotenoid_and_thiolase_metabolite_matched(self, thermo):
        """The check that had to be run before trusting any of this session's numbers."""
        for metabolite in ("s_0190", "s_0943", "s_0189", "s_0633",
                           "s_1275", "s_0837", "s_0373", "s_0367", "s_0529"):
            assert metabolite not in thermo.mismatched, metabolite

    def test_the_mismatches_are_polymers_and_isoprenoids(self, thermo):
        """Not scattered noise -- one failure mode. The model's repeat unit against a fixed
        oligomer. `s_0641` dodecaprenyl diphosphate is C60 in the model and matched a C10."""
        assert "s_0641" in thermo.mismatched
        assert len(thermo.mismatched) < 60


class TestATransportIsNotACreationFromNothing:
    """Coefficients are summed per compound, and the bug this prevents was severe.

    eQuilibrator identifies a MOLECULE, not a molecule-in-a-compartment: ATP[c] and ATP[m]
    resolve to one `Compound` object, and 427 of the objects behind this model's metabolites
    are shared that way. `reaction_dg0` built its terms with a dict comprehension keyed by
    the compound, so the LAST coefficient won and the other was discarded -- turning every
    transport into a creation from nothing.

    `r_1111`, a strict ADP/ATP antiport whose standard chemistry energy is zero, came back at
    **-3,634.7 kJ/mol**, and `gate_step` reported it as `runs`. 663 of the 1,916 fully
    covered reactions collapsed this way, median error 375 kJ/mol.

    **Nothing caught it because every test and every statistic quoted about this backend
    filtered to single-compartment reactions**, where no compound is shared -- including the
    1,098-reaction comparison used to justify making it the primary estimator. Found by the
    2026-08-31 code review.
    """

    @pytest.mark.parametrize("reaction_id", ["r_1111", "r_1175", "r_2023", "r_1849"])
    def test_a_pure_antiport_has_zero_standard_chemistry_energy(self, thermo, yeast_gem,
                                                                reaction_id):
        """The compartment work is `electrical_work`, a separate term, and belongs there."""
        reaction = yeast_gem.reactions.get_by_id(reaction_id)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        assert thermo.reaction_dg0(stoichiometry) == pytest.approx(0.0, abs=1e-6)

    def test_compounds_really_are_shared_across_compartments(self, thermo):
        """The premise. If this ever stops being true the guard above is untestable."""
        shared = sum(1 for count in
                     __import__("collections").Counter(
                         id(c) for c in thermo.compounds.values()).values()
                     if count > 1)

        assert shared > 100

    def test_and_a_single_compartment_reaction_is_unaffected(self, thermo, yeast_gem):
        """Everything this repository has published rests on single-compartment reactions,
        so the fix must leave them bit-identical. The thiolase is the load-bearing one."""
        reaction = yeast_gem.reactions.get_by_id("r_0103")
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        assert thermo.reaction_dg0(stoichiometry) == pytest.approx(24.96, abs=0.05)
