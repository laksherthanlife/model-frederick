"""The solver's own concentrations, put to the layer that can say whether they work.

`pathway/solve.py` produces a content for every node and `bridge/thermodynamic.py` consumes
concentrations, and until now nothing joined them -- so a pathway could solve a pool for a
step whose driving force points the other way and report it as an ordinary float. These
tests pin the join.

Three things are being checked and they are different kinds of claim.

**The arithmetic**, on invented energies. Anything named ``made_up`` below is a Tier 0
number chosen to make a case sharp; it is not a measurement and nothing about the chemistry
of any real pathway follows from it. What follows is that ``dG = dG0 + RT ln Q``, that the
comparison against zero is the right way round, and that the three states stay three.

**The third state.** ``cannot_say`` must never be reported as ``runs``. The tables cover
44.5% of yeast-GEM's REACTIONS -- 1840 of 4131, so ``cannot_say`` is the verdict for 55.5%
of them and is the COMMON case, not the exceptional one. (54.7% is metabolite coverage,
1536/2806. Quoting that one for a gate that gates reactions overstates the reach by ten
points.) A gate that cleared what it could not compute would be worse than no gate at all --
it would hand back a clearance it never earned.

**The one case with a known answer**, on the vendored ModelSEED energies: yeast-GEM's
cytosolic thiolase `r_0103`, +15.6 kJ/mol, threshold 221 uM acetyl-CoA, glucose below at
about 10 uM and ethanol above at about 425 uM. `tests/test_thiolase_threshold.py` pins
those numbers against a hand-written formula; these pin that the gate reaches the same ones
through the API a solver would call, on a real solved `PathwaySolution`. Both files pin
ARITHMETIC. Neither pins that this is the mechanism: the two concentrations come from two
laboratories with two quench protocols, and `scripts/thiolase_threshold.py` names the
measurement that would settle it.

**And the unit gap.** mmol/gDCW becomes molar only by dividing by a cytosolic volume nobody
here has measured, so the tests below show what the verdict does across the whole plausible
range rather than at one value -- including one case that flips inside it, which is what a
conclusion resting on the convention looks like.
"""

from __future__ import annotations

import inspect
from copy import deepcopy

import pytest

from ystwin import paths
from ystwin.bridge.thermodynamic import (
    KJ_PER_KCAL,
    NULL_CUE_UNCERTAINTY_KCAL,
    RefutedEnergy,
    ThermodynamicData,
    dg_prime,
    reaction_dg0,
    reaction_dg0_uncertainty,
)
from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS
from ystwin.pathway.solve import NodeKinetics, solve_pathway
from ystwin.pathway.spec import Node, PathwaySpec, RateLaw, load_pathway
from ystwin.pathway.thermo_gate import (
    CYTOSOLIC_VOLUMES_ML_PER_GDCW,
    Feasibility,
    ThermodynamicallyBlocked,
    Unresolved,
    content_from_molar,
    gate_across_volumes,
    gate_pathway,
    gate_step,
    molar_from_content,
    require_feasible,
    substrate_threshold_m,
    volume_dependent_steps,
)

# Cofactor and precursor levels the beta-carotene solve does not pin. Order-of-magnitude
# placeholders and nothing rests on their values: the steps they belong to are refused on
# their energies, before any quotient is formed.
CAROTENOID_BACKGROUND = {
    "s_0189": 1e-5,   # GGPP, the precursor
    "s_0633": 1e-4,   # diphosphate
    "s_0687": 1e-4,   # FAD
    "s_0689": 1e-4,   # FADH2
    "s_0794": 1e-7,   # proton, pH 7
}

# yeast-GEM's cytosolic thiolase, 2 acetyl-CoA <=> acetoacetyl-CoA + CoA. Pinned against the
# SBML in tests/test_thiolase_threshold.py::test_the_reaction_is_the_one_claimed, and written
# out here so the cheap tests do not need an SBML parse to reach the shape of it.
THIOLASE = {"s_0373": -2.0, "s_0367": 1.0, "s_0529": 1.0}
ACETYL_COA, ACETOACETYL_COA, COA = "s_0373", "s_0367", "s_0529"

# The cofactor levels scripts/thiolase_threshold.py holds the other two participants at, and
# the two measured acetyl-CoA pools in umol/gDW -- glucose- and ethanol-limited CEN.PK113-7D
# at D = 0.05. Two laboratories and two quench protocols; that script carries the caveats.
COA_M = 100e-6
ACETOACETYL_COA_M = 1e-6
GLUCOSE_UMOL_PER_GDW = 0.0199
ETHANOL_UMOL_PER_GDW = 0.85


def _umol_per_gdw(amount: float) -> float:
    """A measurement in umol/gDW as the mmol/gDCW the solver would have produced."""
    return amount * 1e-3


class _Metabolite:
    """The three attributes `reaction_energy` reads off a cobra metabolite, and no more.

    A stub rather than the real model, so the one case that needs a MEMBRANE -- and there
    is no membrane anywhere in a cytosolic heterologous pathway -- does not drag an SBML
    parse into the cheap half of this file.
    """

    def __init__(self, identifier: str, charge: int, compartment: str):
        self.id = identifier
        self.charge = charge
        self.compartment = compartment


class _Reaction:
    def __init__(self, identifier: str, metabolites: dict):
        self.id = identifier
        self.metabolites = metabolites


@pytest.fixture
def made_up_energies():
    """Formation energies chosen to make a case, not measured.

    ``m_prod`` at +5 kJ/mol against zero for the rest makes the condensation below uphill by
    exactly 5, which puts its threshold inside the pools the synthetic chain solves. Real
    energies come from the vendored tables and appear only in the classes that say so.
    """
    return ThermodynamicData({"m_sub": 0.0, "m_prod": 5.0, "m_side": 0.0}, coverage=1.0)


# One uphill condensation: two substrate, one product, one spectator. The thiolase's shape,
# which is what makes the substrate threshold a square root.
CONDENSATION = {"m_sub": -2.0, "m_prod": 1.0, "m_side": 1.0}
SIDE_M = 100e-6


@pytest.fixture
def chain():
    """Two nodes, the first of which holds a pool, so the solve pins a real concentration.

    ``proportional`` rather than ``passthrough`` on purpose: a passthrough node solves to
    exactly zero and the whole point of this fixture is a step gated on the solver's own
    number.
    """
    return PathwaySpec(
        product="product_x", organism="Saccharomyces cerevisiae",
        nodes=(Node("substrate", RateLaw.PROPORTIONAL, enzyme="E1"),
               Node("product_x", RateLaw.PASSTHROUGH, molar_mass_g_per_mol=100.0)),
        entry_enzyme="E1", precursor_metabolite="m_x")


def _solved(chain, entry_flux, growth_rate=0.1):
    return solve_pathway(chain, entry_flux, growth_rate,
                         {"substrate": NodeKinetics(rate_constant=1.0)})


def _gate(solution, thermo, volume, **overrides):
    arguments = dict(steps={"product_x": CONDENSATION},
                     metabolites={"substrate": "m_sub", "product_x": "m_prod"},
                     thermo=thermo, background_m={"m_side": SIDE_M})
    arguments.update(overrides)
    return gate_pathway(solution, cytosolic_volume_ml_per_gdcw=volume, **arguments)


@pytest.fixture(scope="module")
def measured_energies(request):
    """The vendored ModelSEED energies, or a skip.

    Not a module-level ``importorskip``: pytfa is needed for these three classes and for
    nothing else in this file, and skipping the whole file on a machine without it would
    take the arithmetic tests -- which need no data at all -- down with it.
    """
    pytest.importorskip("pytfa")
    return request.getfixturevalue("thermo_data")


class TestTheUnitGapIsAConventionAndTheApiSaysSo:
    def test_a_content_becomes_a_concentration_by_one_division(self):
        assert molar_from_content(1.0, 2.0) == pytest.approx(0.5)

    def test_the_round_trip_returns_the_content(self):
        assert content_from_molar(molar_from_content(3.7e-4, 2.7), 2.7) == pytest.approx(
            3.7e-4, rel=1e-15)

    def test_a_cell_with_no_volume_is_refused(self):
        """Zero volume is an infinite concentration, which would report every uphill step
        as running. The refusal names the range because there is no right value to default
        to -- see CYTOSOLIC_VOLUMES_ML_PER_GDCW."""
        with pytest.raises(ValueError, match="must be positive"):
            molar_from_content(1.0, 0.0)

    def test_a_negative_content_is_refused(self):
        with pytest.raises(ValueError, match="non-negative"):
            molar_from_content(-1e-6, 2.0)

    def test_the_glucose_measurement_is_the_ten_micromolar_that_gets_quoted(self):
        """0.0199 umol/gDW is only "about 10 uM" once a volume is assumed, and this is the
        assumption. Two thirds of the distance to the threshold is in this line."""
        molar = molar_from_content(_umol_per_gdw(GLUCOSE_UMOL_PER_GDW), 2.0)

        assert molar * 1e6 == pytest.approx(9.95, rel=1e-3)

    def test_the_volumes_are_the_ones_the_thiolase_script_reports_across(self):
        """Both files convert the same measurements. A range that drifted apart would let
        the script and the gate disagree about the same state with nothing to point at."""
        assert CYTOSOLIC_VOLUMES_ML_PER_GDCW == (1.0, 2.0, 2.7)

    def test_the_volume_has_no_default_anywhere_a_caller_can_reach(self):
        """The one thing this API must not do is let the convention in silently. If a
        default ever appears here, every verdict downstream acquires an unstated assumption
        and no caller has to notice."""
        parameter = inspect.signature(gate_pathway).parameters[
            "cytosolic_volume_ml_per_gdcw"]

        assert parameter.default is inspect.Parameter.empty


class TestThreeStatesAndTheThirdIsNotTheFirst:
    def test_a_downhill_step_runs(self):
        downhill = ThermodynamicData({"a": 0.0, "b": -20.0}, coverage=1.0)

        step = gate_step({"a": -1.0, "b": 1.0}, {"a": 1e-4, "b": 1e-4}, downhill)

        assert step.feasibility == Feasibility.RUNS

    def test_an_uphill_step_at_a_thin_pool_cannot_run(self, made_up_energies):
        step = gate_step(CONDENSATION,
                         {"m_sub": 1e-6, "m_prod": 1e-4, "m_side": SIDE_M},
                         made_up_energies)

        assert step.feasibility == Feasibility.CANNOT_RUN

    def test_the_same_step_runs_once_the_substrate_pool_is_large_enough(self,
                                                                       made_up_energies):
        step = gate_step(CONDENSATION,
                         {"m_sub": 1e-2, "m_prod": 1e-4, "m_side": SIDE_M},
                         made_up_energies)

        assert step.feasibility == Feasibility.RUNS

    def test_a_step_sitting_exactly_at_equilibrium_cannot_run(self):
        """dG = 0 is no net flux in either direction, and the balance that produced the
        pool assumed one. Reporting equilibrium as forward flux would put the boundary in
        the wrong place by exactly the case that defines it."""
        flat = ThermodynamicData({"a": 0.0, "b": 0.0}, coverage=1.0)

        step = gate_step({"a": -1.0, "b": 1.0}, {"a": 1e-4, "b": 1e-4}, flat)

        assert step.feasibility == Feasibility.CANNOT_RUN

    def test_an_uncovered_participant_is_not_reported_as_running(self, made_up_energies):
        """The failure this module would be worthless without. `reaction_dg0` returns None
        rather than a partial sum when a metabolite is uncovered; if that None became a
        verdict of `runs`, 45% of yeast-GEM would be cleared by not being in the tables."""
        step = gate_step({"m_sub": -1.0, "not_in_the_tables": 1.0},
                         {"m_sub": 1e-4, "not_in_the_tables": 1e-4}, made_up_energies)

        assert step.feasibility == Feasibility.CANNOT_SAY

    def test_and_it_names_the_metabolite_the_tables_do_not_reach(self, made_up_energies):
        step = gate_step({"m_sub": -1.0, "not_in_the_tables": 1.0},
                         {"m_sub": 1e-4, "not_in_the_tables": 1e-4}, made_up_energies)

        assert "not_in_the_tables" in step.note

    def test_an_uncovered_step_carries_no_energy_at_all(self, made_up_energies):
        """Not a partial sum dressed as a number. A caller reading `dg0_kj_per_mol` off a
        `cannot_say` step must get None, not the covered metabolites' total."""
        step = gate_step({"m_sub": -1.0, "not_in_the_tables": 1.0},
                         {"m_sub": 1e-4, "not_in_the_tables": 1e-4}, made_up_energies)

        assert step.dg0_kj_per_mol is None

    def test_a_transport_step_written_without_its_counter_ion_is_not_reported_as_running(
            self):
        """The second way to have no energy, and it is a different fact from the first.
        `reaction_energy` refuses a bare uniport of a highly charged species because
        charging it the full zF*dPsi gives energies no carrier delivers -- +77 kJ/mol for
        PRPP transport. Every metabolite here is covered, so this would sail through a gate
        that only asked whether the tables reached them."""
        crossing = ThermodynamicData({"x_c": 0.0, "x_m": 0.0}, coverage=1.0)
        uniport = _Reaction("t_1", {_Metabolite("x_c", -4, "c"): -1.0,
                                    _Metabolite("x_m", -4, "m"): 1.0})

        step = gate_step({"x_c": -1.0, "x_m": 1.0}, {"x_c": 1e-4, "x_m": 1e-4},
                         crossing, reaction=uniport)

        assert step.feasibility == Feasibility.CANNOT_SAY

    def test_a_participant_with_no_concentration_is_not_reported_as_running(
            self, made_up_energies):
        step = gate_step(CONDENSATION, {"m_sub": 1e-2, "m_prod": 1e-4}, made_up_energies)

        assert step.feasibility == Feasibility.CANNOT_SAY

    def test_and_the_note_says_where_that_concentration_would_come_from(
            self, made_up_energies):
        """A refusal that does not say what would lift it is a dead end for the caller."""
        step = gate_step(CONDENSATION, {"m_sub": 1e-2, "m_prod": 1e-4}, made_up_energies)

        assert "m_side" in step.note and "background_m" in step.note

    def test_a_verdict_carries_the_concentrations_it_was_reached_at(self, made_up_energies):
        """A verdict without them is not reproducible: a threshold is a statement about all
        of the participants at once, so the same step is feasible or not depending on
        numbers that would otherwise never be recorded anywhere."""
        step = gate_step(CONDENSATION,
                         {"m_sub": 1e-2, "m_prod": 1e-4, "m_side": SIDE_M},
                         made_up_energies)

        assert step.concentrations_m == {"m_sub": 1e-2, "m_prod": 1e-4, "m_side": SIDE_M}

    def test_the_energy_is_the_one_the_formula_gives(self, made_up_energies):
        """dG0 + RT ln Q, computed independently through `dg_prime`. If the gate ever
        applied a factor of its own -- a sign, a temperature, a unit -- this is what would
        catch it."""
        concentrations = {"m_sub": 1e-3, "m_prod": 1e-4, "m_side": SIDE_M}
        step = gate_step(CONDENSATION, concentrations, made_up_energies)

        assert step.dg_kj_per_mol == pytest.approx(
            dg_prime(5.0, CONDENSATION, concentrations), rel=1e-12)


class TestTheGateReadsTheSolversOwnNumbers:
    def test_a_thin_solve_cannot_run_its_own_last_step(self, chain, made_up_energies):
        """The contradiction this module exists to find: the balance moved carbon into a
        pool through a step that cannot carry it at the concentration the balance itself
        returned."""
        report = _gate(_solved(chain, 1e-3), made_up_energies, 2.0)

        assert report.step("product_x").feasibility == Feasibility.CANNOT_RUN

    def test_the_same_step_runs_when_the_solve_puts_more_in_the_pool(self, chain,
                                                                    made_up_energies):
        report = _gate(_solved(chain, 1e-1), made_up_energies, 2.0)

        assert report.step("product_x").feasibility == Feasibility.RUNS

    def test_the_concentration_gated_is_the_solved_content_over_the_volume(
            self, chain, made_up_energies):
        """Nothing else may enter. If the gate ever reached for a nominal or a literature
        concentration for a node the solve pinned, it would stop testing the solve."""
        solution = _solved(chain, 1e-1)
        report = _gate(solution, made_up_energies, 2.0)

        assert report.step("product_x").concentrations_m["m_sub"] == pytest.approx(
            solution.node("substrate").content_mmol_per_gdcw / 2.0, rel=1e-15)

    def test_a_background_that_clashes_with_a_solved_content_is_refused(
            self, chain, made_up_energies):
        """Silently preferring either one would answer a different question than the one
        asked. The refusal names the metabolite and both ways out."""
        with pytest.raises(ValueError, match="m_sub"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0,
                  background_m={"m_side": SIDE_M, "m_sub": 1e-3})

    def test_nodes_with_no_declared_step_are_reported_rather_than_omitted(
            self, chain, made_up_energies):
        """A gate that said "nothing is blocked" while looking at one step of two has not
        cleared the pathway, and the caller's only way to know is for the report to say
        what it never looked at."""
        report = _gate(_solved(chain, 1e-1), made_up_energies, 2.0)

        assert report.ungated == ("substrate",)

    def test_a_metabolite_keyed_by_something_that_is_not_a_node_is_refused(
            self, chain, made_up_energies):
        """The quiet way the guarantee above breaks. A mistyped node name maps nothing, so
        the solved content never reaches Q and `background_m` fills the gap instead -- and
        the clash check cannot fire, because no solved value claimed the metabolite. The
        gate would then report a verdict about an assumed level while saying it tested the
        solver's own. Here that flips `product_x` from runs to cannot_run."""
        with pytest.raises(KeyError, match="substarte"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0,
                  metabolites={"substarte": "m_sub", "product_x": "m_prod"},
                  background_m={"m_side": SIDE_M, "m_sub": 1e-6})

    def test_two_nodes_declared_to_be_the_same_species_are_refused(self, chain,
                                                                   made_up_energies):
        """One content would silently replace the other, and which one is whichever the
        spec happens to list last. A node IS a species, so this is a mapping that cannot be
        true rather than one that is merely ambiguous."""
        with pytest.raises(ValueError, match="both declared to be"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0,
                  metabolites={"substrate": "m_sub", "product_x": "m_sub"})

    def test_a_step_keyed_by_something_that_is_not_a_node_is_refused(self, chain,
                                                                     made_up_energies):
        with pytest.raises(KeyError, match="typo"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0,
                  steps={"typo": CONDENSATION})

    def test_a_reaction_id_with_no_model_to_resolve_it_is_refused(self, chain,
                                                                 made_up_energies):
        """And the refusal names both ways out, because the second one -- give the
        stoichiometry directly -- is how a step gets gated without an SBML parse."""
        with pytest.raises(ValueError, match="model"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0,
                  steps={"product_x": "r_0103"})

    def test_a_step_that_is_neither_an_id_nor_a_stoichiometry_is_refused(self, chain,
                                                                        made_up_energies):
        with pytest.raises(TypeError, match="product_x"):
            _gate(_solved(chain, 1e-1), made_up_energies, 2.0, steps={"product_x": 17})


class TestAPassthroughNodeCannotBeGatedFromTheSolveAlone:
    """The shipped PHB spec declares both CoA-thioester intermediates `passthrough`, so the
    solver returns exactly zero for each and there is no concentration to gate on. The
    energies here are invented; what is being pinned is which state that produces.
    """

    @pytest.fixture
    def solved_phb(self):
        return solve_pathway(load_pathway("phb"), 9.93e-2, 0.1)

    @pytest.fixture
    def energies(self):
        return ThermodynamicData(
            {ACETYL_COA: 0.0, ACETOACETYL_COA: 5.0, COA: 0.0}, coverage=1.0)

    def _gate_phb(self, solution, energies, **overrides):
        arguments = dict(steps={"acetoacetyl_coa": THIOLASE},
                         metabolites={"acetoacetyl_coa": ACETOACETYL_COA},
                         thermo=energies, cytosolic_volume_ml_per_gdcw=2.0)
        arguments.update(overrides)
        return gate_pathway(solution, **arguments)

    def test_the_entry_step_cannot_be_gated(self, solved_phb, energies):
        report = self._gate_phb(solved_phb, energies)

        assert report.step("acetoacetyl_coa").feasibility == Feasibility.CANNOT_SAY

    def test_and_the_note_says_the_pool_is_declared_absent_rather_than_unknown(
            self, solved_phb, energies):
        """Two different facts that would otherwise read the same. "Nobody told me" is
        fixed by supplying a number; "the spec asserts there is no pool" is a claim the
        spec makes, and the fix is to change the rate law or to admit the assumption."""
        note = self._gate_phb(solved_phb, energies).step("acetoacetyl_coa").note

        assert "passthrough" in note and "exactly zero" in note

    def test_a_measured_level_for_the_declared_zero_lifts_it(self, solved_phb, energies):
        """Background is consulted exactly where the solve pins nothing, which is what
        makes a passthrough node gateable from a measurement without the gate ever
        overriding a number the solver did produce."""
        report = self._gate_phb(
            solved_phb, energies,
            background_m={ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M})

        assert report.step("acetoacetyl_coa").feasibility != Feasibility.CANNOT_SAY

    def test_two_of_the_three_nodes_are_not_gated_at_all(self, solved_phb, energies):
        assert self._gate_phb(solved_phb, energies).ungated == (
            "r3_hydroxybutyryl_coa", "phb")


class TestWhatChangesAcrossTheVolumeRange:
    """The conversion is a convention, so a verdict that depends on which value was chosen
    is an artefact of the choice and must be visible as one. Q here is volume-dependent
    because one participant is a fixed background level while two are solved contents --
    which is the ordinary case, since cofactors are measured in molar and pools are not.
    """

    def test_the_sweep_returns_one_report_per_volume(self, chain, made_up_energies):
        reports = gate_across_volumes(
            _solved(chain, 1e-1), steps={"product_x": CONDENSATION},
            metabolites={"substrate": "m_sub", "product_x": "m_prod"},
            thermo=made_up_energies, background_m={"m_side": SIDE_M})

        assert [r.cytosolic_volume_ml_per_gdcw for r in reports] == [1.0, 2.0, 2.7]

    def test_a_verdict_that_holds_across_the_whole_range_is_not_flagged(self, chain,
                                                                       made_up_energies):
        reports = tuple(_gate(_solved(chain, 1e-1), made_up_energies, v)
                        for v in CYTOSOLIC_VOLUMES_ML_PER_GDCW)

        assert volume_dependent_steps(reports) == ()

    def test_a_verdict_that_flips_inside_the_range_is_named(self, chain, made_up_energies):
        """The case the sweep exists for. At this flux the step runs at 1.0 mL/gDCW and
        cannot run at 2.7, so the "verdict" is a statement about the convention."""
        reports = tuple(_gate(_solved(chain, 1.32e-2), made_up_energies, v)
                        for v in CYTOSOLIC_VOLUMES_ML_PER_GDCW)

        assert volume_dependent_steps(reports) == ("product_x",)

    def test_and_the_flip_is_a_real_disagreement_between_the_two_ends(self, chain,
                                                                     made_up_energies):
        thin = _gate(_solved(chain, 1.32e-2), made_up_energies, 1.0)
        thick = _gate(_solved(chain, 1.32e-2), made_up_energies, 2.7)

        assert (thin.step("product_x").feasibility,
                thick.step("product_x").feasibility) == (Feasibility.RUNS,
                                                         Feasibility.CANNOT_RUN)

    def test_sweeping_no_volumes_at_all_is_refused(self, chain, made_up_energies):
        with pytest.raises(ValueError, match="at least one"):
            gate_across_volumes(
                _solved(chain, 1e-1), volumes=(), steps={"product_x": CONDENSATION},
                metabolites={"substrate": "m_sub", "product_x": "m_prod"},
                thermo=made_up_energies, background_m={"m_side": SIDE_M})


class TestTheThresholdIsAStatementAboutEveryParticipant:
    def test_at_the_threshold_the_energy_is_exactly_zero(self, made_up_energies):
        """The definition, checked by putting the answer back through the other formula. A
        threshold that did not zero dG would be a number with no meaning."""
        held = {"m_prod": 1e-4, "m_side": SIDE_M}
        threshold = substrate_threshold_m(CONDENSATION, "m_sub", held, 5.0)

        assert dg_prime(5.0, CONDENSATION, {**held, "m_sub": threshold}) == pytest.approx(
            0.0, abs=1e-9)

    def test_a_substrate_runs_forward_above_its_threshold(self, made_up_energies):
        held = {"m_prod": 1e-4, "m_side": SIDE_M}
        threshold = substrate_threshold_m(CONDENSATION, "m_sub", held, 5.0)

        assert gate_step(CONDENSATION, {**held, "m_sub": threshold * 1.01},
                         made_up_energies).feasibility == Feasibility.RUNS

    def test_and_cannot_run_below_it(self, made_up_energies):
        held = {"m_prod": 1e-4, "m_side": SIDE_M}
        threshold = substrate_threshold_m(CONDENSATION, "m_sub", held, 5.0)

        assert gate_step(CONDENSATION, {**held, "m_sub": threshold * 0.99},
                         made_up_energies).feasibility == Feasibility.CANNOT_RUN

    def test_a_product_threshold_points_the_other_way(self, made_up_energies):
        """The coefficient divides, so the inequality flips sign with it. A gate that
        assumed every threshold was a floor would have this exactly backwards for the half
        of participants that are products."""
        held = {"m_sub": 1e-2, "m_side": SIDE_M}
        threshold = substrate_threshold_m(CONDENSATION, "m_prod", held, 5.0)

        assert gate_step(CONDENSATION, {**held, "m_prod": threshold * 1.01},
                         made_up_energies).feasibility == Feasibility.CANNOT_RUN

    def test_a_metabolite_that_does_not_take_part_has_no_threshold(self):
        with pytest.raises(ValueError, match="does not take part"):
            substrate_threshold_m(CONDENSATION, "m_elsewhere",
                                  {"m_prod": 1e-4, "m_side": SIDE_M}, 5.0)

    def test_a_missing_other_participant_is_refused_rather_than_assumed(self):
        """Assuming a level for it would make the threshold a number about that assumption,
        and nothing in the returned float would say so."""
        with pytest.raises(KeyError, match="m_side"):
            substrate_threshold_m(CONDENSATION, "m_sub", {"m_prod": 1e-4}, 5.0)


class TestRefusingIsSeparateFromReporting:
    def test_a_clean_report_does_not_raise(self, chain, made_up_energies):
        require_feasible(_gate(_solved(chain, 1e-1), made_up_energies, 2.0))

    def test_a_blocked_step_raises_and_names_it(self, chain, made_up_energies):
        with pytest.raises(ThermodynamicallyBlocked, match="product_x"):
            require_feasible(_gate(_solved(chain, 1e-3), made_up_energies, 2.0))

    def test_the_refusal_says_what_would_lift_it(self, chain, made_up_energies):
        """Three different things lift it and they are three different claims -- a larger
        pool, a measured cofactor level, or a different volume convention. A refusal naming
        none of them leaves the caller with a raise and no move."""
        with pytest.raises(ThermodynamicallyBlocked) as raised:
            require_feasible(_gate(_solved(chain, 1e-3), made_up_energies, 2.0))

        assert "substrate_threshold_m" in str(raised.value)
        assert "mL/gDCW" in str(raised.value)

    def test_the_refusal_quotes_the_energies_it_refused_on(self, chain, made_up_energies):
        with pytest.raises(ThermodynamicallyBlocked, match=r"dG = \+"):
            require_feasible(_gate(_solved(chain, 1e-3), made_up_energies, 2.0))

    def test_a_step_nobody_could_compute_is_not_refused(self, chain, made_up_energies):
        """Deliberate, and the opposite of the `cannot_say` rule elsewhere in this file. A
        gate must not CLEAR what it could not compute; refusing it instead would be a
        refusal with no energy behind it, which is a fabricated result shaped like caution.
        The step is in `report.cannot_say` for a caller who wants it treated as a failure.
        """
        report = _gate(_solved(chain, 1e-3), made_up_energies, 2.0, background_m={})

        require_feasible(report)

        assert len(report.cannot_say) == 1

    def test_the_summary_states_the_volume_the_verdicts_are_conditional_on(
            self, chain, made_up_energies):
        summary = _gate(_solved(chain, 1e-3), made_up_energies, 2.0).summary()

        assert "1 cannot run" in summary and "2 mL/gDCW" in summary


@pytest.mark.integration
@pytest.mark.skipif(
    paths.yeast_gem() is None or paths.thermo_dir() is None,
    reason="needs yeast-GEM and the thermodynamic tables; see docs/REPRODUCING.md")
class TestTheOneCaseWithAKnownAnswer:
    """r_0103 on the vendored ModelSEED energies, through the API a solver would call.

    Everything above runs on invented numbers. This class runs on the tables, and it is the
    only place in this file where a verdict is a statement about real chemistry -- and even
    here it is a statement about the arithmetic of a lead, not about a mechanism. The two
    acetyl-CoA measurements come from two laboratories with two quench protocols.
    """

    def test_the_standard_energy_is_the_one_the_lead_rests_on(self, measured_energies):
        step = gate_step(THIOLASE,
                         {ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M},
                         measured_energies)

        assert step.dg0_kj_per_mol == pytest.approx(38.07, abs=0.5)

    def test_the_threshold_is_the_nineteen_millimolar_that_is_claimed(self,
                                                                      measured_energies):
        """Reproduced through this module rather than through the script's closed form.
        `substrate_threshold_m` solves the general case, and a coefficient of -2 is what
        turns it into the square root the script writes out."""
        step = gate_step(THIOLASE,
                         {ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M},
                         measured_energies)

        threshold = substrate_threshold_m(THIOLASE, ACETYL_COA, step.concentrations_m,
                                          step.dg0_kj_per_mol)

        assert threshold * 1e6 == pytest.approx(19042.0, rel=0.05)

    def test_the_threshold_in_the_units_the_solver_would_have_to_reach(self,
                                                                      measured_energies):
        """A threshold in molar is not actionable by a layer that works in mmol/gDCW, and
        the number that is depends on the volume convention."""
        step = gate_step(THIOLASE,
                         {ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M},
                         measured_energies)
        threshold = substrate_threshold_m(THIOLASE, ACETYL_COA, step.concentrations_m,
                                          step.dg0_kj_per_mol)

        assert content_from_molar(threshold, 2.0) == pytest.approx(3.81e-2, rel=0.05)

    @pytest.mark.parametrize("volume", CYTOSOLIC_VOLUMES_ML_PER_GDCW)
    def test_glucose_cannot_run_it_at_any_plausible_volume(self, measured_energies,
                                                           volume):
        step = gate_step(
            THIOLASE,
            {ACETYL_COA: molar_from_content(_umol_per_gdw(GLUCOSE_UMOL_PER_GDW), volume),
             ACETOACETYL_COA: ACETOACETYL_COA_M, COA: COA_M},
            measured_energies)

        assert step.feasibility == Feasibility.CANNOT_RUN

    @pytest.mark.parametrize("volume", CYTOSOLIC_VOLUMES_ML_PER_GDCW)
    def test_ethanol_cannot_run_it_either_at_any_plausible_volume(self, measured_energies,
                                                                  volume):
        """This asserted ``RUNS`` until 2026-08-30 and the inversion is the record.

        The lead was that the two feeds straddle the threshold. They do not: the +15.6
        kJ/mol it rested on was a unit error, and at the corrected +38.07 the threshold is
        19 mM against a largest plausible ethanol pool of 850 uM. It is not a near miss and
        no volume convention closes it."""
        step = gate_step(
            THIOLASE,
            {ACETYL_COA: molar_from_content(_umol_per_gdw(ETHANOL_UMOL_PER_GDW), volume),
             ACETOACETYL_COA: ACETOACETYL_COA_M, COA: COA_M},
            measured_energies)

        assert step.feasibility == Feasibility.CANNOT_RUN

    def test_the_gem_reaction_gives_the_same_energy_as_its_stoichiometry(
            self, measured_energies, yeast_gem):
        """The path a caller takes with a real model, against the one the cheap tests take.
        A disagreement would mean the reaction id resolves to something other than the
        reaction these tests are written about."""
        concentrations = {ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M}
        reaction = yeast_gem.reactions.get_by_id("r_0103")

        through_the_model = gate_step(
            {m.id: c for m, c in reaction.metabolites.items()}, concentrations,
            measured_energies, reaction=reaction)

        assert through_the_model.dg0_kj_per_mol == pytest.approx(
            gate_step(THIOLASE, concentrations, measured_energies).dg0_kj_per_mol,
            rel=1e-12)

    @pytest.mark.parametrize(
        "amount,expected",
        [(GLUCOSE_UMOL_PER_GDW, Feasibility.CANNOT_RUN),
         (ETHANOL_UMOL_PER_GDW, Feasibility.CANNOT_RUN)])
    def test_a_solved_phb_pathway_is_gated_end_to_end_through_the_model(
            self, measured_energies, yeast_gem, amount, expected):
        """Everything joined: the shipped PHB spec, solved, its entry step resolved from
        yeast-GEM by id, gated at the measured acetyl-CoA for each feed. This is the call a
        solver would make.

        The ethanol row read ``RUNS`` until 2026-08-30. Both feeds are now below the
        threshold, which is the refutation of the thiolase lead rather than a change in what
        this test exercises: the wiring it pins is unaffected."""
        solution = solve_pathway(load_pathway("phb"), 9.93e-2, 0.1)

        report = gate_pathway(
            solution, steps={"acetoacetyl_coa": "r_0103"},
            metabolites={"acetoacetyl_coa": ACETOACETYL_COA}, thermo=measured_energies,
            cytosolic_volume_ml_per_gdcw=2.0, model=yeast_gem,
            background_m={ACETYL_COA: molar_from_content(_umol_per_gdw(amount), 2.0),
                          ACETOACETYL_COA: ACETOACETYL_COA_M, COA: COA_M})

        assert report.step("acetoacetyl_coa").feasibility == expected

    def test_and_the_glucose_state_is_refused_with_the_reaction_named(
            self, measured_energies, yeast_gem):
        solution = solve_pathway(load_pathway("phb"), 9.93e-2, 0.1)
        report = gate_pathway(
            solution, steps={"acetoacetyl_coa": "r_0103"},
            metabolites={"acetoacetyl_coa": ACETOACETYL_COA}, thermo=measured_energies,
            cytosolic_volume_ml_per_gdcw=2.0, model=yeast_gem,
            background_m={
                ACETYL_COA: molar_from_content(_umol_per_gdw(GLUCOSE_UMOL_PER_GDW), 2.0),
                ACETOACETYL_COA: ACETOACETYL_COA_M, COA: COA_M})

        with pytest.raises(ThermodynamicallyBlocked, match="r_0103"):
            require_feasible(report)


class TestTheGuardsOnHowTheGateIsCalled:
    """Six ways a caller could get a confident verdict out of a question the gate never asked.

    Every one of these was reachable, and each produces a NUMBER rather than an error, which
    is what makes them worth pinning: a gate that refuses is doing its job, and a gate that
    answers the wrong question in the right format is indistinguishable from one that works.
    """

    def test_a_stoichiometry_and_reaction_that_disagree_are_refused(self, made_up_energies):
        """dG0 comes from the reaction and Q from the mapping, so a mismatched pair is neither.

        The sharpest form: hexokinase's energy against the thiolase's quotient. Nothing about
        the pair is individually wrong, and the report would name the thiolase throughout.
        """
        hexokinase = _Reaction("r_0534", {_Metabolite("m_sub", 0, "c"): -1.0,
                                          _Metabolite("m_prod", 0, "c"): 1.0})
        with pytest.raises(ValueError, match="different reactions"):
            gate_step(THIOLASE, {ACETYL_COA: 1e-3, ACETOACETYL_COA: 1e-6, COA: 1e-4},
                      made_up_energies, reaction=hexokinase)

    def test_the_same_pair_is_accepted_when_it_does_agree(self, made_up_energies):
        """The guard has to let the honest call through, or it is just a ban on `reaction=`."""
        matching = _Reaction("r_x", {_Metabolite("m_sub", 0, "c"): -2.0,
                                     _Metabolite("m_prod", 0, "c"): 1.0,
                                     _Metabolite("m_side", 0, "c"): 1.0})
        step = gate_step(CONDENSATION, {"m_sub": 1e-3, "m_prod": 1e-6, "m_side": SIDE_M},
                         made_up_energies, reaction=matching)
        assert step.feasibility in (Feasibility.RUNS, Feasibility.CANNOT_RUN)

    def test_a_step_written_in_the_consuming_direction_is_refused(self, made_up_energies):
        """dG0's sign is the direction's sign, and yeast-GEM writes r_0103 REVERSIBLY.

        Bounds (-1000, 1000) carry no claim about which way it runs, so the +15.606 kJ/mol
        the thiolase lead rests on is inherited from an SBML authoring convention. Written
        the other way it is -15.606 and the glucose verdict reverses.
        """
        with pytest.raises(ValueError, match="CONSUMES it"):
            gate_step(CONDENSATION, {"m_sub": 1e-3, "m_prod": 1e-6, "m_side": SIDE_M},
                      made_up_energies, node="substrate", produces="m_sub")

    def test_a_mistyped_product_metabolite_is_refused_rather_than_ignored(
            self, made_up_energies):
        """A typo leaves the solved content out of Q and `background_m` fills the gap."""
        with pytest.raises(ValueError, match="does not appear"):
            gate_step(CONDENSATION, {"m_sub": 1e-3, "m_prod": 1e-6, "m_side": SIDE_M},
                      made_up_energies, node="product_x", produces="m_prodd")

    def test_a_volume_in_litres_rather_than_millilitres_is_refused(self):
        """The slip is a factor of 1000 and ln Q never objects to the result."""
        with pytest.raises(ValueError, match="pure water"):
            molar_from_content(1e-3, 2.0e-9)

    def test_the_ceiling_is_above_every_pool_this_repository_solves(self):
        """The bound must not be reachable by an honest number, or it is a bug of its own."""
        assert molar_from_content(50.0, 1.0) == pytest.approx(50.0)

    def test_requiring_feasibility_of_a_report_that_gated_nothing_raises(self, chain,
                                                                        made_up_energies):
        """Zero blocked steps out of zero gated steps is not a pass."""
        report = _gate(_solved(chain, 1e-3), made_up_energies, 2.0, steps={})
        assert report.steps == () and len(report.ungated) == 2
        with pytest.raises(ThermodynamicallyBlocked, match="nothing was gated"):
            require_feasible(report)

    def test_a_single_volume_sweep_cannot_report_volume_independence(self, chain,
                                                                    made_up_energies):
        """`()` from one report reads as 'nothing rests on the convention'. It cannot know."""
        one = gate_across_volumes(_solved(chain, 1e-3), volumes=(2.0,),
                                  steps={"product_x": CONDENSATION},
                                  metabolites={"substrate": "m_sub", "product_x": "m_prod"},
                                  thermo=made_up_energies, background_m={"m_side": SIDE_M})
        with pytest.raises(ValueError, match="at least two volumes"):
            volume_dependent_steps(one)


class TestWhatTheJoinIsActuallyExercisedOn:
    """The gate's reason to exist is putting SOLVED concentrations to real energies.

    It does not do that anywhere yet, and this pins the gap rather than papering it. On the
    synthetic `chain` the solved pool does reach Q -- the plumbing works -- but the energies
    are invented. On `phb`, the only real-energy pathway gated here, both thioesters are
    `passthrough` and solve to exactly zero and acetyl-CoA is not a node at all, so every
    concentration in that verdict comes from `background_m`. The two halves never meet.
    """

    def test_on_invented_energies_the_solved_pool_does_reach_the_quotient(
            self, chain, made_up_energies):
        """Doubling the solved content moves dG by RT ln 2, which is the join working."""
        low = _gate(_solved(chain, 1e-3), made_up_energies, 2.0).step("product_x")
        high = _gate(_solved(chain, 2e-3), made_up_energies, 2.0).step("product_x")
        assert high.concentrations_m["m_sub"] == pytest.approx(
            2 * low.concentrations_m["m_sub"])
        assert high.dg_kj_per_mol < low.dg_kj_per_mol

    def test_but_on_the_real_pathway_every_concentration_is_an_assumption(self):
        """phb's own nodes contribute nothing to Q, so the real-energy verdict is background.

        Not a failure of the gate -- a fact about the spec. `phb` declares both thioesters
        `passthrough`, which asserts they hold no appreciable pool, and a node at exactly
        zero is `cannot_say` by construction. Until a spec keeps a pool on a step the tables
        cover, the join has a seam and no join.
        """
        solution = solve_pathway(load_pathway("phb"), 9.93e-2, 0.1)
        assert [s.content_mmol_per_gdcw for s in solution.nodes[:2]] == [0.0, 0.0]
        # The precursor the chain is pulled off, and not a node of it -- an exact
        # name, because "acetoacetyl_coa" contains "acetyl" and a substring test here
        # passes for the wrong reason.
        assert "acetyl_coa" not in {s.name for s in solution.nodes}


# --------------------------------------------------------------------------------------
# The two kinds of silence.
#
# `cannot_say` was one state carrying a sentence, and two of the things it meant were being
# read as the same thing. "The tables do not cover this species" is closed by supplying an
# identifier. "The tables cover it and the number must not be used" is closed by nothing
# except a measurement. Every C40 carotenoid reported the first and was in the second, and
# `data/pathways/beta_carotene.toml` carried the wrong one as prose for as long as the spec
# has existed. These classes pin that the states stay apart.
# --------------------------------------------------------------------------------------

# Two molecules that share almost all of their structure -- the shape a group-contribution
# error cancels in, and the reason a reaction's uncertainty is not the quadrature of its
# metabolites'. `big` is a scaffold both sides carry; `tail` is the one group that changes.
CUE_ERRORS_KCAL = {"big": 10.0, "tail": 0.5}
SHARED_SCAFFOLD_CUES = {"a": {"big": 1.0, "tail": 1.0}, "b": {"big": 1.0, "tail": 2.0}}
ISOMERISATION = {"a": -1.0, "b": 1.0}


def _energies(by_metabolite, *, cues=None, refuted=None, cue_errors=None):
    """A `ThermodynamicData` built by hand, for a case that needs no tables at all."""
    return ThermodynamicData(
        dict(by_metabolite), coverage=1.0,
        cues_by_metabolite=dict(cues or {}),
        cue_uncertainty_kcal=dict(cue_errors if cue_errors is not None else CUE_ERRORS_KCAL),
        refuted=dict(refuted or {}))


REFUSED = RefutedEnergy("cpd99999", "invented compound",
                        "two estimators disagree by more than they quote",
                        "tests/test_thermo_gate.py, not a measurement")


class TestAnEnergyThatIsRefusedIsNotAnEnergyThatIsAbsent:
    """The distinction the carotenoid pathway needed and did not have.

    Everything here runs on invented energies and an invented refutation; what is pinned is
    which STATE each case produces, not any chemistry. The real rows and their citations are
    in `data/thermo/refuted_energies.tsv` and the class below reads them.
    """

    def test_an_uncovered_metabolite_is_reported_as_no_energy(self):
        step = gate_step(ISOMERISATION, {"a": 1e-4, "b": 1e-4}, _energies({"a": 0.0}))

        assert step.feasibility == Feasibility.CANNOT_SAY
        assert step.reason == Unresolved.NO_ENERGY
        assert step.dg0_kj_per_mol is None

    def test_the_no_energy_note_says_to_check_for_a_naming_gap_first(self):
        """The error this repository actually made. The tables carried all three C40
        carotenoids and the gap was that nothing had declared which species the invented
        ids were, so `no_energy` was read as "no such measurement exists" for months."""
        step = gate_step(ISOMERISATION, {"a": 1e-4, "b": 1e-4}, _energies({"a": 0.0}))

        assert "NAMING gap" in step.note

    def test_a_refused_metabolite_is_a_different_state_from_an_absent_one(self):
        step = gate_step(ISOMERISATION, {"a": 1e-4, "b": 1e-4},
                         _energies({"a": 0.0, "b": -20.0}, refuted={"b": REFUSED}))

        assert step.feasibility == Feasibility.CANNOT_SAY
        assert step.reason == Unresolved.REFUTED_UNCERTAINTY

    def test_and_it_reports_the_energy_it_declined_to_use(self):
        """Withholding the number as well would make the refusal uncheckable. -20 kJ/mol
        is a large driving force and a reader is entitled to see that that is what was set
        aside, and how far from zero it was."""
        step = gate_step(ISOMERISATION, {"a": 1e-4, "b": 1e-4},
                         _energies({"a": 0.0, "b": -20.0}, refuted={"b": REFUSED}))

        assert step.dg0_kj_per_mol == pytest.approx(-20.0)
        assert step.dg_kj_per_mol is None

    def test_the_note_carries_the_row_s_own_reason_and_source(self):
        """A refusal with no citation is a result nobody can overturn."""
        note = gate_step(ISOMERISATION, {"a": 1e-4, "b": 1e-4},
                         _energies({"a": 0.0, "b": -20.0},
                                   refuted={"b": REFUSED})).note

        assert REFUSED.reason in note and REFUSED.source in note
        assert "cover it" in note

    def test_a_refusal_is_not_lifted_by_supplying_the_concentrations(self):
        """The order matters. If the concentration check ran first, a caller who supplied a
        background would get a verdict computed off the refused number, and the note would
        say the pathway had been gated."""
        with_everything = gate_step(ISOMERISATION, {"a": 1e-3, "b": 1e-9},
                                    _energies({"a": 0.0, "b": -20.0},
                                              refuted={"b": REFUSED}))

        assert with_everything.reason == Unresolved.REFUTED_UNCERTAINTY

    def test_a_refused_step_does_not_make_require_feasible_raise(self):
        """`cannot_say` never raises, and this is the case that makes that matter: the
        beta-carotene chain is one Elizondo 2025 measured running, and refusing it on an
        extrapolated energy would be a refutation of a strain rather than of a model."""
        solution = solve_pathway(load_pathway("phb"), 9.93e-2, 0.1)
        report = gate_pathway(
            solution, steps={"acetoacetyl_coa": THIOLASE},
            metabolites={"acetoacetyl_coa": ACETOACETYL_COA},
            thermo=_energies({ACETYL_COA: 0.0, ACETOACETYL_COA: 5.0, COA: 0.0},
                             refuted={ACETOACETYL_COA: REFUSED}),
            cytosolic_volume_ml_per_gdcw=2.0,
            background_m={ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M})

        require_feasible(report)  # does not raise
        assert report.because(Unresolved.REFUTED_UNCERTAINTY) == report.steps

    def test_asking_for_a_reason_that_does_not_exist_is_refused(self):
        """`()` from a typo reads as "no step failed that way", which is a positive claim
        from a query incapable of finding one."""
        report = gate_pathway(
            _solved(load_pathway("phb"), 9.93e-2), steps={}, metabolites={},
            thermo=_energies({}), cytosolic_volume_ml_per_gdcw=2.0)
        with pytest.raises(ValueError, match="not a reason"):
            report.because("no_enrgy")


class TestTheUncertaintyIsPropagatedOverCuesAndNotOverFormationEnergies:
    """A reaction's error is not the quadrature of its metabolites' errors.

    Group contribution builds a formation energy out of the groups a molecule is made of, so
    two molecules sharing a scaffold share most of their error and it cancels in the
    difference. Counting it twice is how a well-determined reaction acquires an uncertainty
    larger than itself -- see the real numbers in the class below. `pytfa` propagates over
    net cues for this reason and so does this.
    """

    def test_a_shared_scaffold_cancels_out_of_the_reaction(self):
        sigma = reaction_dg0_uncertainty(
            ISOMERISATION, _energies({"a": 0.0, "b": 1.0}, cues=SHARED_SCAFFOLD_CUES))

        # Only the one `tail` group survives the difference: 0.5 kcal/mol.
        assert sigma == pytest.approx(0.5 * 4.184)

    def test_the_quadrature_of_the_two_formation_errors_would_be_twenty_times_larger(self):
        """The mistake, written out arithmetically so its size is on the page.

        Each molecule's formation error is dominated by the `big` scaffold both of them
        carry: sqrt(10^2 + 0.5^2) and sqrt(10^2 + 2*0.5^2) kcal/mol. Adding those in
        quadrature counts the scaffold twice and gives 14.2 kcal/mol for a difference in
        which it is not there at all.
        """
        import math

        substrate = math.hypot(1 * 10.0, 1 * 0.5)
        product = math.hypot(1 * 10.0, 2 * 0.5)
        naive = math.hypot(substrate, product)
        over_cues = reaction_dg0_uncertainty(
            ISOMERISATION, _energies({"a": 0.0, "b": 1.0},
                                     cues=SHARED_SCAFFOLD_CUES)) / KJ_PER_KCAL

        assert naive == pytest.approx(14.2, abs=0.1)
        assert naive > 20 * over_cues

    def test_an_exactly_cancelling_cue_vector_reports_pytfa_s_floor_not_zero(self):
        """Zero would let any threshold pass. A decomposition that sees no difference
        between the two sides is reporting its own resolution, and 2 kcal/mol is the value
        `pytfa`'s own `prepare(null_error_override=2)` uses for it."""
        sigma = reaction_dg0_uncertainty(
            ISOMERISATION,
            _energies({"a": 0.0, "b": 1.0},
                      cues={"a": {"big": 1.0}, "b": {"big": 1.0}}))

        assert sigma == pytest.approx(NULL_CUE_UNCERTAINTY_KCAL * KJ_PER_KCAL)

    def test_a_metabolite_with_no_decomposition_gives_no_uncertainty_at_all(self):
        """The same refusal `reaction_dg0` makes on an uncovered metabolite and for the
        same reason: an error summed over some of a reaction is a different reaction's."""
        assert reaction_dg0_uncertainty(
            ISOMERISATION,
            _energies({"a": 0.0, "b": 1.0}, cues={"a": {"big": 1.0}})) is None

    def test_the_error_is_reported_beside_an_ordinary_verdict_too(self):
        """Not only beside the withdrawn ones. A `runs` inside its own error bar is the
        case a reader most needs to see, and no gate can see it for them."""
        step = gate_step(ISOMERISATION, {"a": 1e-3, "b": 1e-6},
                         _energies({"a": 0.0, "b": -20.0}, cues=SHARED_SCAFFOLD_CUES))

        assert step.feasibility == Feasibility.RUNS
        assert step.dg0_uncertainty_kj_per_mol == pytest.approx(0.5 * 4.184)


class TestTheUncertaintyThresholdIsOptInAndSaysWhy:
    """Withdrawing a verdict smaller than its own error is a claim about the result, not a
    default. That decision belongs to whoever owns the result, which is the whole reason
    the threshold is opt-in.

    This docstring used to add that `phb`'s entry step "runs at -3.3 kJ/mol on ethanol
    against a 4.5 kJ/mol propagated error, so a one-sigma rule applied at load time would
    silently delete this module's only real result". Withdrawn 2026-09-04: the kcal/kJ
    correction in commit 8809344 moved that step to +38.07 kJ/mol against sigma 4.55, so it
    returns `cannot_run` at 8.4 sigma uphill and a one-sigma rule would delete nothing. The
    tests below never depended on it -- they use a synthetic -1 kJ/mol step -- which is why
    the claim survived here unchallenged. The opt-in default stands on the principle alone.
    """

    def _small_driving_force(self, **kwargs):
        # dG0 = -1 kJ/mol with everything at 1 M, so dG = -1 against a 2.1 kJ/mol error.
        return gate_step(ISOMERISATION, {"a": 1.0, "b": 1.0},
                         _energies({"a": 0.0, "b": -1.0}, cues=SHARED_SCAFFOLD_CUES),
                         **kwargs)

    def test_by_default_the_verdict_stands_and_the_error_is_reported_beside_it(self):
        step = self._small_driving_force()

        assert step.feasibility == Feasibility.RUNS
        assert step.dg0_uncertainty_kj_per_mol > abs(step.dg_kj_per_mol)

    def test_asking_for_one_sigma_withdraws_it(self):
        step = self._small_driving_force(uncertainty_sigma=1.0)

        assert step.feasibility == Feasibility.CANNOT_SAY
        assert step.reason == Unresolved.UNCERTAINTY_EXCEEDS_EFFECT

    def test_and_the_withdrawn_step_still_carries_both_numbers(self):
        step = self._small_driving_force(uncertainty_sigma=1.0)

        assert step.dg_kj_per_mol == pytest.approx(-1.0)
        assert step.dg0_uncertainty_kj_per_mol == pytest.approx(0.5 * 4.184)

    def test_a_driving_force_well_outside_the_error_survives_the_threshold(self):
        step = gate_step(ISOMERISATION, {"a": 1.0, "b": 1.0},
                         _energies({"a": 0.0, "b": -50.0}, cues=SHARED_SCAFFOLD_CUES),
                         uncertainty_sigma=1.0)

        assert step.feasibility == Feasibility.RUNS

    def test_a_reaction_with_no_error_at_all_is_not_withdrawn_by_the_threshold(self):
        """No decomposition means no error, and no error is not a small error. Withdrawing
        here would make an uncomputable uncertainty behave like an infinite one."""
        step = gate_step(ISOMERISATION, {"a": 1.0, "b": 1.0},
                         _energies({"a": 0.0, "b": -1.0}), uncertainty_sigma=1.0)

        assert step.feasibility == Feasibility.RUNS


@pytest.mark.integration
@pytest.mark.skipif(
    paths.yeast_gem() is None or paths.thermo_dir() is None,
    reason="needs yeast-GEM and the thermodynamic tables; see docs/REPRODUCING.md")
class TestTheCarotenoidGapWasANameAndNotACoverageGap:
    """On the real tables. The claim `docs/research/CAROTENOID_ENERGIES.md` corrects.

    Before `data/thermo/heterologous_metabolites.tsv` existed, `dgf` returned None for all
    three C40 species and the gate said `no_energy` -- and `beta_carotene.toml` wrote that
    down as the tables carrying no energy for a carotenoid. They carry all three.
    """

    @pytest.mark.parametrize("metabolite,seed",
                             [("phytoene_c", "cpd03205"), ("lycopene_c", "cpd03217"),
                              ("betacarotene_c", "cpd01420")])
    def test_each_declared_species_now_resolves_to_its_modelseed_compound(
            self, measured_energies, metabolite, seed):
        assert measured_energies.dgf(metabolite) is not None
        assert measured_energies.seed_by_metabolite[metabolite] == seed

    def test_they_resolve_without_the_pathway_being_installed_on_any_model(
            self, measured_energies):
        """The gate is asked about a spec before a model has the chemistry in it -- see
        `predict.py`, which gates from the declared stoichiometry with `model=None`. An
        annotation on a metabolite object cannot help there, which is why the declaration
        is data and not only an annotation."""
        assert set(measured_energies.declared) == {
            "phytoene_c", "lycopene_c", "betacarotene_c"}

    def test_declaring_them_does_not_change_what_the_tables_cover(self, measured_energies):
        """Coverage is over the model's own metabolites. Counting three species the model
        does not have would report a reach the tables do not have."""
        assert measured_energies.coverage == pytest.approx(0.547, abs=0.001)

    def test_the_three_steps_are_refused_for_the_stated_reason_not_for_a_missing_name(self):
        """The whole track in one assertion. Same verdict as before -- `cannot_say` on all
        three -- and a different fact underneath it."""
        spec = load_pathway("beta_carotene")
        report = gate_pathway(
            solve_pathway(spec, 4.2e-4, 0.18, kinetics=BETA_CAROTENE_KINETICS),
            steps=spec.thermo_steps, metabolites=spec.thermo_metabolites,
            thermo=ThermodynamicData.load(), cytosolic_volume_ml_per_gdcw=2.0,
            background_m=CAROTENOID_BACKGROUND)

        assert set(report.verdicts.values()) == {Feasibility.CANNOT_SAY}
        assert set(report.reasons.values()) == {Unresolved.REFUTED_UNCERTAINTY}

    def test_and_every_one_of_them_now_has_an_energy_to_report(self):
        spec = load_pathway("beta_carotene")
        report = gate_pathway(
            solve_pathway(spec, 4.2e-4, 0.18, kinetics=BETA_CAROTENE_KINETICS),
            steps=spec.thermo_steps, metabolites=spec.thermo_metabolites,
            thermo=ThermodynamicData.load(), cytosolic_volume_ml_per_gdcw=2.0,
            background_m=CAROTENOID_BACKGROUND)

        assert all(s.dg0_kj_per_mol is not None for s in report.steps)
        assert all(s.dg0_uncertainty_kj_per_mol is not None for s in report.steps)

    def test_the_thiolase_error_is_the_one_that_makes_its_verdict_readable(
            self, measured_energies):
        """+38.07 kJ/mol against 4.5, so the standard energy is 8.5 sigma from zero. It was
        3.4 sigma at the pre-2026-08-30 +15.6, and the ethanol verdict then sat inside its
        own error bar at dG = -3.3. It no longer does: correcting the unit moved the energy
        away from zero, so the verdict is further from the noise, not nearer it."""
        step = gate_step(THIOLASE,
                         {ACETYL_COA: 425e-6, ACETOACETYL_COA: ACETOACETYL_COA_M,
                          COA: COA_M},
                         measured_energies)

        assert step.dg0_uncertainty_kj_per_mol == pytest.approx(4.5, abs=0.2)
        assert abs(step.dg0_kj_per_mol) > 8 * step.dg0_uncertainty_kj_per_mol

    def test_the_naive_quadrature_over_formation_errors_would_be_eight_times_that(
            self, measured_energies):
        """36.2 kJ/mol on a reaction whose energy is +38, produced entirely by counting one
        CoA moiety three times. The reason `reaction_dg0_uncertainty` propagates over cues.
        This number is a property of the tabulated formation errors, which are in kJ/mol in
        the source and were never on the mixed scale, so the unit correction did not move
        it -- which is why it is still 36.2 here."""
        import math

        errors = {"s_0373": 14.79, "s_0367": 14.98, "s_0529": 14.44}  # kJ/mol, tabulated
        naive = math.sqrt(sum((THIOLASE[m] * e) ** 2 for m, e in errors.items()))

        assert naive == pytest.approx(36.2, abs=0.5)
        assert naive > 7 * reaction_dg0_uncertainty(THIOLASE, measured_energies)


def _prepare_energy_template(energies):
    from ystwin.bridge.equilibrator import PreferredEnergies

    component = energies.component if isinstance(energies, PreferredEnergies) else energies
    for compound in set(component.compounds.values()):
        tuple(compound.microspecies)
    return _copy_energy_template(energies)


def _copy_energy_template(energies):
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker

    from ystwin.bridge.equilibrator import PreferredEnergies

    component = energies.component if isinstance(energies, PreferredEnergies) else energies
    cache = component.contribution.ccache
    engine = create_engine(cache.engine.url)
    session = sessionmaker(bind=engine)()
    return deepcopy(energies, {id(cache.engine): engine, id(cache.session): session})


@pytest.mark.integration
@pytest.mark.skipif(
    paths.yeast_gem() is None or paths.thermo_dir() is None,
    reason="needs yeast-GEM and the thermodynamic tables; see docs/REPRODUCING.md")
class TestTheTwoEstimatorsDisagreeFarBeyondTheirQuotedError:
    """The evidence `data/thermo/refuted_energies.tsv` cites, recomputed rather than quoted.

    A refusal that rests on a number in a hand-written row is a refusal nobody re-derives.
    These are the two estimators this repository already ships -- ModelSEED group
    contribution through `bridge/thermodynamic.py`, eQuilibrator component contribution
    through `bridge/equilibrator.py` -- on the three steps of `beta_carotene.toml`.
    """

    @pytest.fixture(scope="class")
    def _both_template(self, yeast_gem_factory):
        pytest.importorskip("equilibrator_api")
        from ystwin.bridge.equilibrator import EquilibratorData
        from ystwin.fba.carotenoid import add_beta_carotene_pathway

        installed = add_beta_carotene_pathway(yeast_gem_factory())
        return (ThermodynamicData.load(),
                _prepare_energy_template(EquilibratorData.load(installed)), installed)

    @pytest.fixture
    def both(self, _both_template, model_copy):
        seed, component, model = _both_template
        component = _copy_energy_template(component)
        try:
            yield deepcopy(seed), component, model_copy(model)
        finally:
            component.contribution.ccache.close()

    @pytest.mark.parametrize("reaction_id,quoted,ratio", [("CRTYB_PSY", 4.9, 4.5),
                                                          ("CRTI", 17.4, 5.6),
                                                          ("CRTYB_LCY", 15.5, 8.0)])
    def test_the_disagreement_is_far_beyond_what_equilibrator_quotes(
            self, both, reaction_id, quoted, ratio):
        """The ratios were 16-22x before the 2026-08-30 unit correction and are 4.5-17.6x
        after it, because part of the old disagreement was the unit error rather than the
        chemistry. The conclusion is unchanged and the reason is
        `test_the_two_estimators_do_not_even_agree_on_the_sign`: an estimator that is 22
        kJ/mol out on a step whose quoted error is 4.9, and 225 kJ/mol out with the wrong
        sign on the next, is not an estimator whose number may be used."""
        seed, equilibrator, model = both
        reaction = model.reactions.get_by_id(reaction_id)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
        component, error = equilibrator.reaction_dg0(stoichiometry, with_error=True)
        group = reaction_dg0(stoichiometry, seed)

        assert error == pytest.approx(quoted, abs=0.6)
        assert abs(group - component) == pytest.approx(ratio * error, rel=0.15)
        assert abs(group - component) > 4 * error

    def test_they_now_agree_on_the_sign_of_every_step(self, both):
        """**This asserted the opposite until 2026-08-30, and the inversion is the point.**

        It read `group < 0 < component` on CRTI -- -59 against +166 -- and called that the
        decisive reason to refuse both numbers. The straddle was not a fact about
        carotenoids. `fba/carotenoid.py` wrote CrtI with four free FAD as the terminal
        acceptor, which is thermodynamically impossible (+166 on a step Verwaal 2007 runs to
        completion) and is not what a flavin oxidase does. Written with O2 the same
        chemistry is -292.7, and the sign disagreement is gone.

        What is left is a magnitude disagreement, which is a weaker claim and the honest
        one. Both estimators now put every step of the pathway well below zero.
        """
        seed, equilibrator, model = both
        for reaction_id in ("CRTE", "CRTYB_PSY", "CRTI", "CRTYB_LCY"):
            reaction = model.reactions.get_by_id(reaction_id)
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

            group = reaction_dg0(stoichiometry, seed)
            component = equilibrator.reaction_dg0(stoichiometry)

            assert group is not None and component is not None, reaction_id
            assert group < 0 and component < 0, reaction_id

    def test_the_carotenoids_are_what_equilibrator_is_least_sure_of(self, both):
        """The uncertainty rises down the chain as the conjugation does, and CRTE -- which
        touches no C40 species at all -- is the one it is surest of.

        The exact ordering moved on 2026-08-30 when CrtI's terminal acceptor was corrected
        from four free FAD to O2: CRTI's error went 12.8 -> 17.4 and overtook CRTYB_LCY's
        15.5, because four O2 and four H2O2 carry their own component uncertainty. So the
        assertion is on CRTE being lowest and the C40 steps being several times it, which is
        the claim, rather than on a total order that a stoichiometry change reshuffles.
        """
        _seed, equilibrator, model = both
        errors = {}
        for reaction_id in ("CRTE", "CRTYB_PSY", "CRTI", "CRTYB_LCY"):
            stoichiometry = {
                m.id: c
                for m, c in model.reactions.get_by_id(reaction_id).metabolites.items()}
            errors[reaction_id] = equilibrator.reaction_dg0(
                stoichiometry, with_error=True)[1]

        assert errors["CRTE"] < errors["CRTYB_PSY"]
        assert errors["CRTE"] < min(errors["CRTI"], errors["CRTYB_LCY"])
        assert min(errors["CRTI"], errors["CRTYB_LCY"]) > 3 * errors["CRTE"]


class TestADeclarationIsCheckedRatherThanTrusted:
    """A hand-written table that attaches one molecule's energy to another molecule's name
    is the failure `_same_skeleton` already exists to catch on the alias route, and a
    declaration is the same route with a shorter path to a typo.
    """

    @pytest.fixture
    def tables(self, tmp_path):
        pytest.importorskip("pytfa")
        return tmp_path

    def _write(self, path, seed, formula):
        path.write_text(
            "metabolite_id\tmodelseed_id\tkegg_id\tname\tformula\tinstalled_by\tsource\n"
            f"x_c\t{seed}\tC00000\tx\t{formula}\ttest\tnot a source\n", encoding="utf-8")
        return path

    def test_a_formula_that_disagrees_with_the_database_is_refused(self, tables):
        bad = self._write(tables / "declared.tsv", "cpd03205", "C6H12O6")
        with pytest.raises(ValueError, match="disagrees"):
            ThermodynamicData.load(heterologous=bad)

    def test_an_id_the_database_does_not_carry_is_refused_by_name(self, tables):
        """A declaration that resolves to nothing is worse than none: it reads as a
        closed gap."""
        bad = self._write(tables / "declared.tsv", "cpd99999999", "C40H64")
        with pytest.raises(ValueError, match="cpd99999999"):
            ThermodynamicData.load(heterologous=bad)


class TestARefusedMagnitudeIsCheckedOnDGAndNotOnDG0:
    """The rule that lets a refused step still return a verdict, and why it has no threshold.

    **The first version of this was a bandaid and this class is the counterexample that
    killed it.** It compared the two estimators' dG0 and called the verdict robust when
    "the gap is smaller than the smaller estimate". A verdict is not about dG0. Component
    -60 and group -35 pass that rule, and at a concentration ratio worth +40 kJ/mol they
    give -20 and +5 -- opposite verdicts from a check that called them robust.

    What replaced it asks the only question that matters: put BOTH estimators at the SAME
    concentrations and see whether they give the same sign. No threshold, because the
    quantity being thresholded was the wrong one.
    """

    @pytest.fixture(scope="class")
    def _both_template(self, yeast_gem_factory):
        pytest.importorskip("equilibrator_api")
        from ystwin.bridge.equilibrator import PreferredEnergies
        from ystwin.fba.carotenoid import add_beta_carotene_pathway

        installed = add_beta_carotene_pathway(yeast_gem_factory())
        return _prepare_energy_template(PreferredEnergies.load(model=installed)), installed

    @pytest.fixture
    def both(self, _both_template, model_copy):
        energies, model = _both_template
        energies = _copy_energy_template(energies)
        try:
            yield energies, model_copy(model)
        finally:
            energies.component.contribution.ccache.close()

    def _cyclase(self, model):
        reaction = model.reactions.get_by_id("CRTYB_LCY")
        return {m.id: c for m, c in reaction.metabolites.items()}

    def test_the_dg0_rule_would_have_passed_a_case_it_should_not(self):
        """Two lines of arithmetic, kept as the reason the rule changed."""
        component, group = -60.0, -35.0
        gap = abs(component - group)

        assert gap < min(abs(component), abs(group))          # the old rule said robust
        assert (component + 40.0 < 0) != (group + 40.0 < 0)   # and the verdicts differ

    def test_at_physiological_concentrations_the_verdict_stands(self, both):
        energies, model = both
        step = gate_step(self._cyclase(model),
                         {"lycopene_c": 1e-6, "betacarotene_c": 1e-6}, energies)

        assert step.feasibility is Feasibility.RUNS
        assert "MAGNITUDE is refused and the VERDICT is not" in step.note

    def test_and_the_note_carries_both_numbers_so_a_reader_can_check(self, both):
        energies, model = both
        step = gate_step(self._cyclase(model),
                         {"lycopene_c": 1e-6, "betacarotene_c": 1e-6}, energies)

        assert "kJ/mol" in step.note and "both estimators give the same sign" in step.note

    def test_where_the_two_estimators_actually_differ_it_refuses(self, both):
        """Driven to a quotient extreme enough that the refused 124 kJ/mol decides the
        sign. Not a physiological state -- the point is that the check FIRES, which the
        dG0 rule would not have."""
        energies, model = both
        step = gate_step(self._cyclase(model),
                         {"lycopene_c": 1e-30, "betacarotene_c": 1.0}, energies)

        assert step.feasibility is Feasibility.CANNOT_SAY
        assert step.reason is Unresolved.REFUTED_UNCERTAINTY
        assert "opposite directions" in step.note

    def test_a_single_estimator_backend_still_refuses_outright(self, measured_energies):
        """The conservative default is unchanged for every caller that does not deliberately
        supply two estimators. One method agreeing with itself is not agreement."""

        assert not hasattr(measured_energies, "both_dg0")


class TestTheDisagreementCheckAppliesToEveryReaction:
    """It ran only on reactions containing a hand-listed compound, and that was the defect.

    `data/thermo/refuted_energies.tsv` was written after looking at the carotenoids and
    nothing else. The standard it applied was not the one the rest of the model was held to:

      CRTYB_PSY  refused by name        gap  22.2 kJ/mol   -- BELOW the model median
      r_4710     reported `runs`        gap 964.4 kJ/mol   -- estimators at -963.0 and +1.4

    A curated list cannot be a standard. Over the 1,098 single-compartment reactions both
    estimators cover, the gap runs 7.3 / 25.2 / 46.8 / 81.8 kJ/mol at the 25th / 50th / 75th
    / 90th percentile, so refusing a step for 22 would refuse most of yeast-GEM. The check is
    now a measured property of the reaction being gated, and the compound list is what
    records why those three were investigated -- not what triggers the check.
    """

    @pytest.fixture(scope="class")
    def _both_template(self, yeast_gem_factory):
        pytest.importorskip("equilibrator_api")
        from ystwin.bridge.equilibrator import PreferredEnergies

        yeast_gem = yeast_gem_factory()
        return _prepare_energy_template(PreferredEnergies.load(model=yeast_gem)), yeast_gem

    @pytest.fixture
    def both(self, _both_template, model_copy):
        energies, model = _both_template
        energies = _copy_energy_template(energies)
        try:
            yield energies, model_copy(model)
        finally:
            energies.component.contribution.ccache.close()

    def test_a_reaction_with_no_listed_compound_is_still_checked(self, both):
        """r_4710: the two estimators put it at -963.0 and +1.4. It used to return `runs`."""
        energies, model = both
        reaction = model.reactions.get_by_id("r_4710")
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        assert all(energies.refutation(m) is None for m in stoichiometry)

        step = gate_step(stoichiometry, {m: 1e-4 for m in stoichiometry}, energies)

        assert step.feasibility is Feasibility.CANNOT_SAY
        assert "opposite directions" in step.note

    def test_a_large_gap_that_agrees_on_direction_is_still_reported(self, both):
        """r_2141: -1281.5 against -410.9, an 871 kJ/mol gap and the same verdict. The check
        is on the direction, not the size, so this is a `runs` and should stay one."""
        energies, model = both
        reaction = model.reactions.get_by_id("r_2141")
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        step = gate_step(stoichiometry, {m: 1e-4 for m in stoichiometry}, energies)

        assert step.feasibility is Feasibility.RUNS

    def test_the_carotenoid_steps_are_not_special(self, both):
        """Two of the four sit at the model-wide median. The refutation implied otherwise."""
        import numpy as np

        energies, model = both
        gaps = []
        for reaction in model.reactions:
            if len({m.compartment for m in reaction.metabolites}) != 1:
                continue
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
            pair = energies.both_dg0(stoichiometry)
            if pair is not None:
                gaps.append(abs(pair[0] - pair[1]))

        assert len(gaps) > 900
        assert 20.0 < float(np.median(gaps)) < 30.0
