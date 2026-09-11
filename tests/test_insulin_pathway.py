"""The insulin precursor is a PROTEIN, and almost every test here exists because of that.

`test_carotenoid_pathway.py` can check a four-reaction chain by naming its metabolites. It
cannot check a stoichiometry with 41 terms computed from an 86-residue string, and the
interesting failure modes are different: a mistranscribed residue, a formula that does not
correspond to any real molecule, an energy cost silently omitted. So the checks here are
mostly ANCHORS -- each one reaches an authority outside `fba/insulin.py` and asks whether
the module still agrees with it.

The strongest is `TestTheFormulaReproducesPublishedInsulin`. Run the module's own machinery
on the two mature insulin chains instead of proinsulin and it has to return the published
molecular formula of human insulin, a number this repository did not choose and cannot
adjust.
"""

from __future__ import annotations

import pathlib

import pytest

pytestmark = pytest.mark.integration

_GLUCOSE_EXCHANGE = "r_1714"
_GROWTH = "r_2111"

#: Cytosolic free amino acid names in Yeast9, for the independent formula route. Names, not
#: ids, because a name is what the SBML and a reader agree on; the ids are resolved from
#: them inside the test so a renumbering fails loudly instead of silently checking nothing.
_FREE_AMINO_ACIDS = {
    "A": "L-alanine", "R": "L-arginine", "N": "L-asparagine", "D": "L-aspartate",
    "C": "L-cysteine", "Q": "L-glutamine", "E": "L-glutamate", "G": "L-glycine",
    "H": "L-histidine", "I": "L-isoleucine", "L": "L-leucine", "K": "L-lysine",
    "M": "L-methionine", "F": "L-phenylalanine", "P": "L-proline", "S": "L-serine",
    "T": "L-threonine", "W": "L-tryptophan", "Y": "L-tyrosine", "V": "L-valine",
}


@pytest.fixture(scope="module")
def _insulin_model_template(yeast_gem_factory):
    """One install for the whole file: copying Yeast9 is seconds, and this is done twelve
    times otherwise."""
    from ystwin.fba.insulin import add_insulin_precursor_pathway

    return add_insulin_precursor_pathway(yeast_gem_factory())


@pytest.fixture
def insulin_model(_insulin_model_template, model_copy):
    return model_copy(_insulin_model_template)


class TestTheSequenceIsTheOneUniProtHas:
    """P01308 is the authority and the module quotes it in four pieces. These check the
    pieces are the ones the feature table delimits and that they reassemble."""

    def test_the_precursor_is_eighty_six_residues(self):
        from ystwin.fba.insulin import PROINSULIN_SEQUENCE

        assert len(PROINSULIN_SEQUENCE) == 86

    def test_the_chains_are_the_lengths_the_feature_table_gives(self):
        """SIGNAL 1..24, PEPTIDE 25..54 (B), PROPEP 57..87 (C), PEPTIDE 90..110 (A). The
        two-residue gaps at 55-56 and 88-89 are the dibasic Kex2/PC sites."""
        from ystwin.fba.insulin import PROINSULIN_SEQUENCE

        b_chain, rr, c_peptide, kr, a_chain = (
            PROINSULIN_SEQUENCE[:30], PROINSULIN_SEQUENCE[30:32],
            PROINSULIN_SEQUENCE[32:63], PROINSULIN_SEQUENCE[63:65],
            PROINSULIN_SEQUENCE[65:])

        assert (len(b_chain), len(c_peptide), len(a_chain)) == (30, 31, 21)
        assert (rr, kr) == ("RR", "KR")
        assert a_chain == "GIVEQCCTSICSLYQLENYCN"
        assert b_chain.startswith("FVNQHL") and b_chain.endswith("YTPKT")

    def test_all_six_cysteines_are_paired(self):
        """DISULFID 31..96, 43..109, 95..100 -- three bonds, six cysteines, none left over."""
        from ystwin.fba.insulin import (DISULFIDE_BONDS, PROINSULIN_SEQUENCE,
                                        residue_composition)

        assert residue_composition(PROINSULIN_SEQUENCE)["C"] == 2 * DISULFIDE_BONDS == 6

    def test_proinsulin_contains_no_methionine_and_no_tryptophan(self):
        """Not trivia. It is why two of `r_4047`'s twenty aminoacyl-tRNA nodes appear with
        coefficient zero below, and it is a cheap check that the sequence was not padded
        with an initiator Met that the signal peptide actually carries."""
        from ystwin.fba.insulin import PROINSULIN_SEQUENCE, residue_composition

        composition = residue_composition(PROINSULIN_SEQUENCE)

        assert "M" not in composition and "W" not in composition
        assert len(composition) == 18

    def test_a_residue_with_no_trna_in_the_model_is_refused(self):
        """Selenocysteine has no tRNA in `r_4047`, so there is nothing to write it against
        and a substitution would be an invented residue in a product formula."""
        from ystwin.fba.insulin import residue_composition

        with pytest.raises(ValueError, match="no tRNA"):
            residue_composition("GIVEQU")


class TestTheStoichiometryComesFromTheHostModel:
    """The alternative is a hardcoded table of twenty `s_` ids, right for one GEM release."""

    def test_it_draws_on_the_same_species_the_protein_pseudoreaction_does(self, insulin_model):
        translation = insulin_model.reactions.get_by_id("INSPRE_TRANSLATION")
        pseudo = insulin_model.reactions.get_by_id("r_4047")
        charged = {m.id for m, c in pseudo.metabolites.items() if c < 0}

        drawn = {m.id for m, c in translation.metabolites.items() if c < 0}

        assert drawn & charged, "translation does not touch the model's aminoacyl-tRNA pool"
        assert drawn - charged <= {"s_0785", "s_0803"}, "unexpected extra substrate"

    def test_the_two_absent_residues_have_no_term_at_all(self, insulin_model):
        """Met-tRNA (s_1148) and Trp-tRNA (s_1527) are in `r_4047` and must not be here."""
        translation = insulin_model.reactions.get_by_id("INSPRE_TRANSLATION")
        ids = {m.id for m in translation.metabolites}

        assert "s_1148" not in ids and "s_1527" not in ids

    def test_the_residue_coefficients_are_the_sequence_counts(self, insulin_model):
        from ystwin.fba.insulin import PROINSULIN_SEQUENCE, residue_composition

        translation = insulin_model.reactions.get_by_id("INSPRE_TRANSLATION")
        by_id = {m.id: c for m, c in translation.metabolites.items()}
        counts = residue_composition(PROINSULIN_SEQUENCE)

        assert by_id["s_1077"] == -counts["L"] == -12    # Leu, the commonest residue
        assert by_id["s_0432"] == -counts["D"] == -1     # Asp, the rarest
        assert by_id["s_0542"] == -counts["C"] == -6     # Cys, and there are six

    def test_a_host_without_a_protein_pseudoreaction_is_refused_by_name(self, yeast_gem):
        from ystwin.fba.insulin import add_insulin_precursor_pathway

        with pytest.raises(KeyError, match="r_9999"):
            add_insulin_precursor_pathway(yeast_gem, protein_pseudoreaction_id="r_9999")


class TestTheProductFormulaIsDerivedTwiceAndAgrees:
    """The chain formula is built by summing (aminoacyl-tRNA minus tRNA) over the sequence.
    Nothing forces that to be a real molecule, so it is rebuilt by a route that shares no
    code and no metabolites with it."""

    def test_it_equals_the_free_amino_acids_minus_the_peptide_waters(self, insulin_model):
        """The textbook polypeptide: sum of the free amino acids, less one water per bond.

        The two routes share only the model. If the tRNA-difference arithmetic, the
        termination water or the per-residue proton were wrong, this is where it shows."""
        from ystwin.fba.insulin import PROINSULIN_SEQUENCE, residue_composition

        counts = residue_composition(PROINSULIN_SEQUENCE)
        length = sum(counts.values())
        expected: dict[str, float] = {}
        for code, n in counts.items():
            matches = [m for m in insulin_model.metabolites
                       if m.name == _FREE_AMINO_ACIDS[code] and m.compartment == "c"]
            assert len(matches) == 1, f"{_FREE_AMINO_ACIDS[code]} is not unique in the cytosol"
            for symbol, count in matches[0].elements.items():
                expected[symbol] = expected.get(symbol, 0.0) + n * count
        expected["H"] -= 2 * (length - 1)
        expected["O"] -= (length - 1)

        reduced = insulin_model.metabolites.get_by_id("insulin_precursor_red_c")

        assert reduced.elements == pytest.approx(expected)

    def test_the_folded_form_is_six_hydrogens_lighter(self, insulin_model):
        """Three disulfides, two hydrogens each, and nothing else changes."""
        reduced = insulin_model.metabolites.get_by_id("insulin_precursor_red_c")
        folded = insulin_model.metabolites.get_by_id("insulin_precursor_er")

        assert folded.charge == reduced.charge
        assert {s: v for s, v in folded.elements.items() if s != "H"} == \
               {s: v for s, v in reduced.elements.items() if s != "H"}
        assert reduced.elements["H"] - folded.elements["H"] == 6


class TestTheFormulaReproducesPublishedInsulin:
    """THE ONE ANCHOR OUTSIDE THIS REPOSITORY, and the reason to believe the rest.

    Mature human insulin -- the A and B chains, three disulfides, C peptide excised -- has a
    published molecular formula, C257H383N65O77S6, and an average mass of 5807.57 Da. That
    number was not chosen here and cannot be adjusted here. Feed the module's own
    composition machinery the two chains instead of proinsulin and it must land on it.

    The one correction is a protonation state and it is the interesting part. Yeast9 writes
    charged species, so what comes back is the ANION at pH 7, charge -2: four glutamates
    against one arginine and one lysine, plus two N-termini and two C-termini. Add the two
    protons back and the formula is the published one exactly. Getting that off by even one
    hydrogen would mean the per-residue proton bookkeeping in `_chain_composition` is wrong,
    which would put a wrong mass on every gram-per-gram yield this module reports.
    """

    def test_the_two_mature_chains_give_the_published_formula(self, yeast_gem):
        import cobra

        from ystwin.fba.insulin import (_chain_composition, _hill_formula,
                                        _translation_species, residue_composition)

        a_chain, b_chain = "GIVEQCCTSICSLYQLENYCN", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
        charged, free = _translation_species(yeast_gem, "r_4047")
        elements: dict[str, float] = {}
        charge = 0.0
        for chain in (a_chain, b_chain):
            chain_elements, chain_charge = _chain_composition(
                residue_composition(chain), charged, free)
            for symbol, count in chain_elements.items():
                elements[symbol] = elements.get(symbol, 0.0) + count
            charge += chain_charge
        elements["H"] -= 6.0                       # three disulfides
        assert charge == -2.0
        elements["H"] -= charge                    # protonate to the neutral molecule

        assert _hill_formula(elements) == "C257H383N65O77S6"
        assert cobra.Metabolite("x", formula="C257H383N65O77S6").formula_weight == \
            pytest.approx(5807.57, abs=0.05)

    def test_proinsulin_lands_on_its_own_published_mass(self, insulin_model):
        """Same arithmetic, one chain, C peptide included: about 9389 Da neutral."""
        import cobra

        folded = insulin_model.metabolites.get_by_id("insulin_precursor_er")
        neutral = dict(folded.elements)
        neutral["H"] -= folded.charge

        assert cobra.Metabolite("x", formula="C410H638N114O127S6").elements == \
            pytest.approx(neutral)
        assert cobra.Metabolite("x", formula="C410H638N114O127S6").formula_weight == \
            pytest.approx(9388.5, abs=0.5)


class TestTheEnergyCost:
    """Where the ATP and the GTP come from, and what happens when the assumption is dropped."""

    def test_the_charging_atp_is_paid_by_the_host_not_asserted_here(self, insulin_model):
        """Writing insulin off the FREE amino acid pool would skip the synthetases and the
        2 ATP-equivalents per residue they cost. Drawing on the CHARGED pool means the host
        pays it and this module asserts nothing -- but only if every one of the twenty
        really is written ATP to AMP plus diphosphate, which is two high-energy phosphate
        bonds and not one. Checked on all twenty rather than on alanine, because "the model
        already charges it" is the load-bearing half of Martin 2025's split and one
        synthetase written ATP-to-ADP would quietly halve it.
        """
        from ystwin.fba.insulin import _translation_species

        charged, _ = _translation_species(insulin_model, "r_4047")
        for code, metabolite in charged.items():
            makers = [r for r in metabolite.reactions
                      if r.id != "r_4047" and metabolite in r.products]
            assert len(makers) == 1, f"{code}: {[r.id for r in makers]}"
            names = {m.name for m in makers[0].metabolites}
            assert {"ATP", "AMP", "diphosphate"} <= names, f"{code}: {makers[0].reaction}"
            assert makers[0].check_mass_balance() == {}

    def test_two_gtp_per_residue_are_charged(self, insulin_model):
        """Martin 2025 (PMID 40562331): 4 ATP per peptide bond, "2 from PPi formation at
        aminoacyl tRNA synthesis ... and 1 GTP each for the two elongation factors"."""
        from ystwin.fba.insulin import GTP_PER_RESIDUE, PROINSULIN_SEQUENCE

        translation = insulin_model.reactions.get_by_id("INSPRE_TRANSLATION")
        by_id = {m.id: c for m, c in translation.metabolites.items()}
        expected = GTP_PER_RESIDUE * len(PROINSULIN_SEQUENCE)

        assert by_id["s_0785"] == -expected == -172.0     # GTP
        assert by_id["s_0739"] == expected                # GDP
        assert by_id["s_1322"] == expected                # phosphate

    def test_the_assumption_can_be_switched_off_and_then_no_gtp_appears(self, yeast_gem):
        """yeast-GEM's own convention: elongation energy inside the biomass GAM and nowhere
        else. A heterologous protein never reaches it, so this setting charges nothing."""
        from ystwin.fba.insulin import add_insulin_precursor_pathway

        free_model = add_insulin_precursor_pathway(yeast_gem, gtp_per_residue=0.0)
        translation = free_model.reactions.get_by_id("INSPRE_TRANSLATION")

        assert "s_0785" not in {m.id for m in translation.metabolites}
        assert translation.check_mass_balance() == {}


class TestOxidativeFolding:
    """Three disulfides are not free, and what they cost is hydrogen peroxide."""

    def test_one_oxygen_in_and_one_peroxide_out_per_disulfide(self, insulin_model):
        """Gross 2006 (PMID 16407158) on recombinant yeast Ero1p: "reduction of molecular
        oxygen by Ero1p yielded stoichiometric hydrogen peroxide"."""
        from ystwin.fba.insulin import DISULFIDE_BONDS

        folding = insulin_model.reactions.get_by_id("INSPRE_FOLDING")
        by_name = {m.name: c for m, c in folding.metabolites.items()}

        assert by_name["oxygen"] == -float(DISULFIDE_BONDS) == -3.0
        assert by_name["hydrogen peroxide"] == float(DISULFIDE_BONDS) == 3.0

    def test_the_pathway_therefore_cannot_run_without_oxygen(self, insulin_model):
        """A property of the reaction as written, not a measurement of the cell -- and the
        reason it is worth asserting is that Gross 2006 exists because Ero1p's ANAEROBIC
        electron acceptor was unknown."""
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        with insulin_model as m:
            m.reactions.get_by_id("r_1992").lower_bound = 0.0     # no oxygen uptake
            m.objective = PRODUCT_DEMAND_ID

            assert m.slim_optimize() == pytest.approx(0.0, abs=1e-9)

    def test_more_disulfides_than_cysteines_is_refused(self, yeast_gem):
        from ystwin.fba.insulin import add_insulin_precursor_pathway

        with pytest.raises(ValueError, match="cysteines"):
            add_insulin_precursor_pathway(yeast_gem, disulfide_bonds=4)


class TestTheInstallContract:
    """The same contract `fba/carotenoid.py` keeps, checked the same way."""

    def test_every_added_reaction_is_mass_and_charge_balanced(self, insulin_model):
        from ystwin.fba.insulin import PATHWAY_REACTION_IDS

        for rid in PATHWAY_REACTION_IDS:
            reaction = insulin_model.reactions.get_by_id(rid)
            if reaction.boundary:
                continue
            assert reaction.check_mass_balance() == {}, \
                f"{rid} unbalanced: {reaction.check_mass_balance()}"

    def test_installing_the_pathway_twice_is_refused(self, insulin_model):
        from ystwin.fba.insulin import add_insulin_precursor_pathway

        with pytest.raises(ValueError, match="already"):
            add_insulin_precursor_pathway(insulin_model)

    def test_the_source_model_is_left_untouched(self, yeast_gem, insulin_model):
        assert not yeast_gem.metabolites.has_id("insulin_precursor_er")

    def test_adding_a_drain_pathway_does_not_raise_maximum_growth(self, yeast_gem,
                                                                  insulin_model):
        assert insulin_model.slim_optimize() == pytest.approx(yeast_gem.slim_optimize(),
                                                              rel=1e-6)

    def test_insulin_can_carry_flux_when_it_is_the_objective(self, insulin_model):
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        with insulin_model as m:
            m.reactions.get_by_id(_GLUCOSE_EXCHANGE).lower_bound = -10.0
            m.objective = PRODUCT_DEMAND_ID

            assert m.slim_optimize() > 1e-6


class TestTheCeiling:
    """The measured numbers, pinned so they cannot drift silently.

    Every one of these is a BOUND on a model, not a prediction about a cell. There is no
    insulin measurement anywhere in this repository to compare them with, which is stated in
    `TestWhatThisCannotDo` and is the honest end of this exercise.

    Conditions are yeast-GEM's own defaults -- glucose exchange at its shipped lower bound
    of -1 mmol/gDCW/h, aerobic, no oxygen cap -- so these are the same conditions
    `docs/FINDINGS.md` reports growth of 0.08584 under.
    """

    def test_maximum_insulin_flux_with_growth_switched_off(self, insulin_model):
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        with insulin_model as m:
            m.objective = PRODUCT_DEMAND_ID

            assert m.slim_optimize() == pytest.approx(0.009898, rel=1e-3)

    def test_the_mass_yield_on_glucose_is_about_half_a_gram_per_gram(self, insulin_model):
        """0.516 g insulin / g glucose, i.e. 68% of the glucose CARBON ends up in one
        86-residue polypeptide. No process does this, which is the point of a ceiling."""
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        molar_mass = insulin_model.metabolites.get_by_id("insulin_precursor_er").formula_weight
        with insulin_model as m:
            m.objective = PRODUCT_DEMAND_ID
            flux = m.slim_optimize()
        glucose = -insulin_model.reactions.get_by_id(_GLUCOSE_EXCHANGE).lower_bound

        assert flux * molar_mass / (glucose * 180.156) == pytest.approx(0.516, rel=1e-2)

    def test_growth_and_insulin_trade_off_and_the_exchange_rate_is_about_one_to_one(
            self, insulin_model):
        """Roughly a gram of insulin per gram of biomass forgone, which is what a carbon
        budget should give for two things made of the same elements. It is the sanity check
        on the whole construction more than a result."""
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        molar_mass = insulin_model.metabolites.get_by_id("insulin_precursor_er").formula_weight
        mu_max = insulin_model.slim_optimize()
        with insulin_model as m:
            m.objective = PRODUCT_DEMAND_ID
            at_zero_growth = m.slim_optimize()
            m.reactions.get_by_id(_GROWTH).lower_bound = 0.5 * mu_max
            at_half_growth = m.slim_optimize()

        assert at_half_growth < at_zero_growth
        # Fluxes are mmol/gDCW/h against a molar mass in g/mol, so the 1e-3 is the mmol.
        grams_per_gdcw = (
            (at_zero_growth - at_half_growth) * 1e-3 * molar_mass / (0.5 * mu_max))

        assert grams_per_gdcw == pytest.approx(1.07, rel=0.05)

    def test_charging_the_elongation_gtp_costs_about_a_tenth_of_the_ceiling(
            self, yeast_gem, insulin_model):
        """How load-bearing the one assumption in the module is, stated as a number rather
        than as a caveat. Dropping it raises the ceiling by 9%, so it changes the answer and
        does not change the conclusion."""
        from ystwin.fba.insulin import PRODUCT_DEMAND_ID, add_insulin_precursor_pathway

        free_model = add_insulin_precursor_pathway(yeast_gem, gtp_per_residue=0.0)
        with free_model as m:
            m.objective = PRODUCT_DEMAND_ID
            without = m.slim_optimize()
        with insulin_model as m:
            m.objective = PRODUCT_DEMAND_ID
            with_gtp = m.slim_optimize()

        assert without / with_gtp == pytest.approx(1.101, rel=1e-2)


class TestWhatThisCannotDo:
    """The findings that are refusals. They are tests because a refusal that is only prose
    stops being true the moment somebody adds a file, and nothing says so."""

    def test_the_ceiling_cannot_make_any_flux_happen(self, insulin_model):
        """`fba/carotenoid.py` says it in its own docstring: cap a reaction at the measured
        magnitude and FVA returns [0, cap] -- the ceiling moves and the floor stays at
        exactly zero, because an upper bound cannot make a flux mandatory. It is worth
        re-checking on a protein because the failure it guards against is reading a ceiling
        as a forecast, which is the whole risk of this module."""
        from cobra.flux_analysis import flux_variability_analysis

        from ystwin.fba.insulin import PRODUCT_DEMAND_ID

        with insulin_model as m:
            m.reactions.get_by_id(PRODUCT_DEMAND_ID).upper_bound = 0.005
            interval = flux_variability_analysis(
                m, reaction_list=[PRODUCT_DEMAND_ID], fraction_of_optimum=0.0)

        assert interval.loc[PRODUCT_DEMAND_ID, "minimum"] == pytest.approx(0.0, abs=1e-9)
        assert interval.loc[PRODUCT_DEMAND_ID, "maximum"] == pytest.approx(0.005, rel=1e-6)

    def test_there_is_no_insulin_measurement_anywhere_in_this_repository(self):
        """The reason nothing here is calibrated. Checked rather than asserted, because "no
        data exists" is the kind of claim that quietly stops being true.

        If this fails, somebody has added insulin data, and the right response is to find
        out what units it is in -- a titre in mg/L is NOT convertible to the mmol/gDCW/h
        this module reports without a paired biomass concentration.
        """
        repo = pathlib.Path(__file__).resolve().parents[1]
        directories = [repo / "data", repo / "outputs", repo / "figures"]
        hits = [p for d in directories if d.is_dir() for p in d.rglob("*")
                if p.is_file() and "insulin" in p.name.lower()]

        assert hits == []

    def test_there_is_deliberately_no_pathway_spec_for_insulin(self):
        """Three reasons, each independently sufficient, and the weakest one has expired.

        This docstring used to LEAD with the fourth: `pathway/spec.py` refused
        `fate = "secreted"` at load, and a recombinant insulin precursor is secreted by
        design -- expressed behind the MFalpha1 leader precisely so that it leaves the cell.
        On 2026-09-01 that refusal moved to solve time and a secreted terminal node became
        representable, so the reason that was stated first is no longer a reason at all.

        The absence stands on the other three, which the export outlet does not touch:

        1. There is no measurement of any kind behind insulin in this repository -- no
           titre, no content, no rate -- so there is nothing to fit and nothing to score.
        2. A `PathwaySpec` is a chain of metabolite nodes. A polypeptide is not one.
        3. Its entry "enzyme" would be the ribosome, and `entry_enzyme` names a gene whose
           expression drives the flux law. That map does not exist for translation.

        Kastberg 2025 does not lift any of the three: it measures the insulin precursor in a
        chemostat and reports only relative fold changes, with no absolute quantity anywhere.

        This asserts the absence so that the reasoning in `fba/insulin.py`'s docstring and
        the state of `data/pathways/` cannot drift apart.
        """
        from ystwin.pathway.spec import Fate, available_pathways

        assert "insulin" not in available_pathways()
        assert Fate.SECRETED in Fate.ALL, "the fate this product has is still declarable"
