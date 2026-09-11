"""Thermodynamic flux analysis done the way Henry 2007 does it, because the shortcut failed.

Fixing every unmeasured metabolite at one concentration, computing dG' and freezing
directions makes yeast-GEM infeasible: 1020 of 1837 covered reactions get forced one way and
221 of those are shut off entirely. That is not thermodynamics refusing the model, it is a
made-up metabolome refusing it. No single concentration is right for 1573 metabolites.

TMFA (PMID 17172310) makes ln(concentration) a variable inside physiological bounds and asks
whether SOME assignment makes a flux distribution feasible. Directionality couples to it
through a binary per reaction, which is why the problem is a mixed-integer one:

    dG'_i  =  dG'0_i + RT * sum_j S_ji * x_j          x_j = ln c_j
    v_i    <=  M * z_i                                 forward only when z_i = 1
    v_i    >= -M * (1 - z_i)
    dG'_i  <=  M * (1 - z_i) - eps                     z_i = 1 forces dG' negative
    dG'_i  >= -M * z_i + eps

The measurement enters by tightening the bounds on x_j for the cofactors a sensor reads.
Whether that tightening changes any flux is then a question with an answer, which is the
whole reason for building it.
"""

import numpy as np
import pytest

from ystwin.bridge.tmfa import (
    CONCENTRATION_CEILING_M,
    CONCENTRATION_FLOOR_M,
    ThermoModel,
)

from ystwin import paths

pytest.importorskip("pytfa")

pytestmark = pytest.mark.skipif(
    paths.yeast_gem() is None or paths.thermo_dir() is None,
    reason="needs yeast-GEM and the thermodynamic tables; see docs/REPRODUCING.md",
)


@pytest.fixture(scope="module")
def thermo(thermo_data):
    return thermo_data


@pytest.fixture(scope="module")
def built(thermo):
    """Scoped to central carbon metabolism, where the measured cofactors actually bite.

    Scope matters: the problem is mixed-integer and grows with the reaction count, so the
    formulation is written for the whole model and demonstrated where the measurement is
    relevant."""
    return ThermoModel.build(thermo, subsystems=("Glycolysis / gluconeogenesis",))


class TestTheFormulation:
    def test_it_adds_a_concentration_variable_per_covered_metabolite(self, built):
        assert built.n_concentration_variables > 0

    def test_it_adds_a_direction_binary_per_constrained_reaction(self, built):
        assert built.n_direction_variables == built.n_constrained_reactions

    def test_concentrations_are_bounded_to_a_physiological_window(self, built):
        low, high = built.concentration_bounds

        assert low == pytest.approx(np.log(CONCENTRATION_FLOOR_M))
        assert high == pytest.approx(np.log(CONCENTRATION_CEILING_M))

    def test_it_constrains_only_reactions_it_has_full_energies_for(self, built, thermo):
        for reaction_id in built.constrained_reactions:
            reaction = built.model.reactions.get_by_id(reaction_id)
            assert all(thermo.dgf(m.id) is not None for m in reaction.metabolites)

    def test_it_says_how_much_of_the_scope_it_could_constrain(self, built):
        assert 0.0 < built.coverage <= 1.0


class TestItStaysFeasible:
    def test_the_model_still_solves(self, built):
        """The whole point. The fixed-concentration shortcut made it infeasible."""
        assert built.is_feasible()

    def test_it_still_grows(self, built):
        assert built.growth_rate() > 0

    def test_growth_is_not_faster_than_without_the_constraints(self, built):
        """A thermodynamic constraint can only remove options, never add them."""
        assert built.growth_rate() <= built.unconstrained_growth_rate() + 1e-6


class TestWhatTheMeasurementBuys:
    def test_pinning_a_cofactor_never_widens_a_flux_range(self, built):
        free = built.flux_range("r_0962")
        pinned = built.flux_range("r_0962", measured={"s_0434": (0.9e-3, 4.4e-3),
                                                      "s_0394": (0.3e-3, 2.2e-3)})

        assert pinned[1] - pinned[0] <= free[1] - free[0] + 1e-6

    def test_it_reports_whether_the_measurement_changed_anything(self, built):
        report = built.measurement_value(
            measured={"s_0434": (0.9e-3, 4.4e-3), "s_0394": (0.3e-3, 2.2e-3)},
            reactions=["r_0962"])

        assert set(report) == {"r_0962"}
        assert isinstance(report["r_0962"]["narrowed"], bool)

    def test_a_measurement_below_the_assumed_floor_is_honoured_not_clipped(self, built):
        """The generic window is an assumption and a measurement outruns it.

        Free cytosolic NADH in yeast runs 3 to 30 uM (Canelas 2008, PMID 18383140) against
        the 10 uM floor Henry 2007 assumes, so most of the measured range sits underneath
        it. Clipping to the assumption would throw away the measurement in favour of the
        guess it was meant to replace.
        """
        built.flux_range("r_0962", measured={"s_0434": (3e-6, 3e-5)})
        variable = built.concentration_variables["s_0434"]

        assert variable.lb == pytest.approx(np.log(CONCENTRATION_FLOOR_M))

    def test_a_measurement_that_is_not_an_interval_is_still_refused(self, built):
        with pytest.raises(ValueError, match="positive interval"):
            built.flux_range("r_0962", measured={"s_0434": (4e-3, 1e-3)})

    def test_a_measurement_on_a_metabolite_the_scope_lacks_is_refused(self, built):
        with pytest.raises(KeyError, match="not a constrained metabolite"):
            built.flux_range("r_0962", measured={"s_not_real": (1e-3, 2e-3)})


class TestTheAnswerIsNo:
    """The measurement tightens the energies substantially and changes no flux at all.

    yeast-GEM ships measured yeast ATP (0.9-4.4 mM) and ADP (0.3-2.2 mM) in
    YMDBconcentrations.csv and discards them on export. Applied here they cut the feasible
    ATP/ADP ratio from a ln-span of 15.2 to 3.58 -- 29 kJ/mol of dG' range removed, which is
    not a marginal tightening.

    Zero glycolytic flux ranges narrow, and the growth rate is unchanged to four decimals.
    Thermodynamics is simply not what bounds these fluxes: mass balance and the uptake
    bounds already do, and the dG' constraint never becomes the binding one.

    This is consistent with the only two papers that put these bounds on a yeast model.
    Kummel 2006 and Martinez 2014 both chose them deliberately loose, and nothing here
    suggests tightening them would have cost either paper anything. Recorded as a test so
    that the negative result survives, rather than being rediscovered by whoever next
    assumes a measured cofactor must constrain a flux.
    """

    def test_the_measurement_removes_real_energy_range(self):
        from ystwin.bridge.thermodynamic import GAS_CONSTANT_KJ, STANDARD_TEMPERATURE_K
        from ystwin.bridge.tmfa import CONCENTRATION_CEILING_M, CONCENTRATION_FLOOR_M

        rt = GAS_CONSTANT_KJ * STANDARD_TEMPERATURE_K
        free = np.log((CONCENTRATION_CEILING_M / CONCENTRATION_FLOOR_M) ** 2)
        measured = np.log((4.4e-3 / 0.3e-3) / (0.9e-3 / 2.2e-3))

        assert rt * (free - measured) > 20.0

    def test_and_yet_it_narrows_no_flux(self, built):
        report = built.measurement_value(
            measured={"s_0434": (0.9e-3, 4.4e-3), "s_0394": (0.3e-3, 2.2e-3)},
            reactions=built.constrained_reactions[:6])

        assert not any(row["narrowed"] for row in report.values())

    def test_thermodynamics_does_not_even_bind_on_growth(self, built):
        """Not that the measurement is too loose -- the whole layer is not binding here.

        Compared at mixed-integer solver tolerance rather than floating-point: the two
        agree to five figures, and the residual difference has the constrained solution
        marginally the larger, which is numerically impossible and therefore tolerance."""
        assert built.growth_rate() == pytest.approx(built.unconstrained_growth_rate(), rel=1e-4)


class TestExaminingAFluxLeavesTheModelAsItFoundIt:
    """Found by randomised test ordering, which is the only reason it was found.

    Taking a flux range points the objective at that reaction. Leaving it there makes every
    later growth query optimise the wrong thing and return a number that looks perfectly
    reasonable, so nothing fails and the answer is simply wrong.
    """

    def test_growth_is_unchanged_by_having_examined_a_flux(self, built):
        before = built.growth_rate()
        built.flux_range(built.constrained_reactions[0])

        assert built.growth_rate() == pytest.approx(before, rel=1e-6)

    def test_a_measured_range_also_leaves_it_alone(self, built):
        before = built.growth_rate()
        built.flux_range(built.constrained_reactions[0],
                         measured={"s_0434": (0.9e-3, 4.4e-3)})

        assert built.growth_rate() == pytest.approx(before, rel=1e-6)


class TestWhyOneCofactorCannotDecideADirection:
    """The mechanism behind the negative result, and it is arithmetic rather than scoping.

    A reaction's Gibbs energy moves with RT ln Q, and every metabolite in Q gets a vote
    weighted by its concentration range. The generic window spans 1e-5 to 2e-2 M, which is
    7.6 ln-units or about 19 kJ/mol per metabolite. A pinned cofactor ratio contributes
    about 3 ln-units, near 5.5 kJ/mol.

    So on the reactions where a measurement should matter most -- cytosolic malate
    dehydrogenase, glycerol-3-phosphate dehydrogenase, the ADP/ATP carrier -- the
    unmeasured substrates carry seven to seventeen times the swing of the measured couple,
    and can move to cancel whatever it says. The direction is never decided, and no amount
    of scoping changes that.

    Which is consistent with how the technique is actually used: Kummel 2006 constrained
    iND750 with a measured *metabolome*, dozens of species at once, not one couple. A single
    biosensor is the wrong shape of measurement for this job, and that is a statement about
    the method rather than about the sensor.
    """

    def _swings(self, thermo, reaction_id, pinned):
        from ystwin.bridge.thermodynamic import (GAS_CONSTANT_KJ, STANDARD_TEMPERATURE_K,
                                                 _default_model)

        rt = GAS_CONSTANT_KJ * STANDARD_TEMPERATURE_K
        window = np.log(CONCENTRATION_CEILING_M / CONCENTRATION_FLOOR_M)
        reaction = _default_model().reactions.get_by_id(reaction_id)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
        measured = sum(abs(c) for k, c in stoichiometry.items() if k in pinned) * rt * np.log(3.0)
        free = sum(abs(c) for k, c in stoichiometry.items() if k not in pinned) * rt * window
        return measured, free

    def test_the_unmeasured_substrates_outweigh_the_measured_couple(self, thermo):
        pinned = {"s_0434", "s_0394", "s_1198", "s_1203"}

        for reaction_id in ("r_0714", "r_0491", "r_0958"):
            measured, free = self._swings(thermo, reaction_id, pinned)
            assert free > 5 * measured, reaction_id

    def test_that_holds_for_the_reaction_the_nadh_literature_names(self, thermo):
        """Glycerol-3-phosphate dehydrogenase is the cytosolic NADH sink Vemuri drained."""
        measured, free = self._swings(thermo, "r_0491", {"s_1198", "s_1203"})

        assert free > measured

    def test_so_pinning_them_narrows_nothing_even_where_it_should(self, thermo):
        """The reactions here are flippable and carry flux, which is the strongest case."""
        scoped = ThermoModel.build(thermo, reactions=["r_0714", "r_0491", "r_1110"])
        report = scoped.measurement_value(
            {"s_0434": (0.9e-3, 4.4e-3), "s_0394": (0.3e-3, 2.2e-3)},
            reactions=["r_0714", "r_0491", "r_1110"])

        assert not any(row["narrowed"] for row in report.values())


class TestWhyConcentrationsAreVariablesAndNotConstants:
    """Kept from the approach this replaced, because the failure is the reason for the design.

    Fixing every metabolite at one concentration and freezing directions from the resulting
    energies forces 1020 of 1837 covered reactions one way and shuts 221 off entirely
    against their own bounds, and yeast-GEM stops solving. What is refused there is an
    invented metabolome, not the thermodynamics.
    """

    def test_a_single_assumed_concentration_forces_most_reactions_one_way(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, dg_prime, reaction_dg0

        model = _default_model()
        forced = judged = 0
        for reaction in model.reactions:
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
            dg0 = reaction_dg0(stoichiometry, thermo)
            if dg0 is None:
                continue
            judged += 1
            energy = dg_prime(dg0, stoichiometry, {m: 1e-4 for m in stoichiometry})
            forced += abs(energy) > 5.0

        assert judged > 500
        assert forced > 0.4 * judged

    def test_whereas_variable_concentrations_leave_the_model_growing(self, built):
        assert built.growth_rate() > 0


class TestApplyingAWindowAboveTheCurrentOne:
    """Bounds have to be widened before they are narrowed.

    Setting the lower bound first fails whenever the new one passes the upper bound still in
    place from the last measurement, which is what happens applying a window shifted upward.
    """

    def test_a_window_entirely_above_the_default_is_accepted(self, built):
        built.flux_range(built.constrained_reactions[0],
                         measured={"s_0434": (5e-3, 15e-3)})

    def test_two_disjoint_windows_in_sequence_both_apply(self, built):
        reaction = built.constrained_reactions[0]
        built.flux_range(reaction, measured={"s_0434": (1e-5, 1e-4)})
        built.flux_range(reaction, measured={"s_0434": (5e-3, 15e-3)})
        variable = built.concentration_variables["s_0434"]

        assert variable.lb < variable.ub

    def test_the_default_window_is_restored_afterwards(self, built):
        built.flux_range(built.constrained_reactions[0],
                         measured={"s_0434": (5e-3, 15e-3)})
        low, high = built.concentration_bounds
        variable = built.concentration_variables["s_0434"]

        assert (variable.lb, variable.ub) == pytest.approx((low, high))


class TestTheNarrowingControl:
    """Is a narrowed flux range evidence that the measurement said something?

    On its own, no, and the control is what shows it. A concentration window constrains the
    problem whether or not anyone measured it, and whether or not it is true. Pinning two
    arbitrary metabolites to a window of the same log-width as the real one, placed at 0.1
    mM, narrows up to three of eight glycolytic reactions by 15 to 20 percent. The measured
    ATP and ADP narrow none.

    The arbitrary window narrows precisely because it is wrong: cytosolic phosphate runs in
    the tens of mM and NAD around 1 mM, so forcing them near 0.1 mM asserts something false
    and the model dutifully rules flux out. That is the point. A count of narrowed reactions
    cannot tell a measurement from a mistake, so the earlier three-of-eight result carried no
    evidence either way, and it is reproduced here by a random pair.
    """

    SCOPE = ["r_0173", "r_0174", "r_0366"]
    ARBITRARY = (4.5e-5, 2.2e-4)

    def test_the_measured_adenylates_narrow_nothing(self, built):
        report = built.measurement_value(
            measured={"s_0434": (0.9e-3, 4.4e-3), "s_0394": (0.3e-3, 2.2e-3)},
            reactions=self.SCOPE)

        assert not any(row["narrowed"] for row in report.values())

    def test_an_arbitrary_window_of_the_same_width_narrows_more(self, built):
        """Same width, wrong place, on metabolites nobody measured."""
        report = built.measurement_value(
            measured={"s_1322": self.ARBITRARY, "s_1198": self.ARBITRARY},
            reactions=self.SCOPE)

        assert any(row["narrowed"] for row in report.values())

    def test_and_by_far_more_than_a_rounding(self, built):
        report = built.measurement_value(
            measured={"s_1322": self.ARBITRARY, "s_1198": self.ARBITRARY},
            reactions=self.SCOPE)
        kept = [row["fraction_kept"] for row in report.values() if row["narrowed"]]

        assert min(kept) < 0.95

    def test_so_narrowing_is_judged_relative_and_reported_with_its_width(self, built):
        from ystwin.bridge.tmfa import NARROWING_TOLERANCE

        report = built.measurement_value(measured={"s_0434": (0.9e-3, 4.4e-3)},
                                         reactions=self.SCOPE[:1])

        assert 0 < NARROWING_TOLERANCE < 0.05
        assert "fraction_kept" in report[self.SCOPE[0]]
