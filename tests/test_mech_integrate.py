"""The stiff driver, the tau/T audit, and the two things the ratio test cannot decide.

Four claims are pinned here and each one is a number this suite computes rather than repeats:

1. **The stiffness is real.** The catalogue spans 1.1e8 in time constant, and an explicit method
   on a two-state system carrying the fastest measured rate costs four orders of magnitude more
   right-hand-side evaluations than the pinned implicit one -- for an answer that agrees to seven
   figures. Cost, not accuracy, which is exactly what stiffness means.

2. **The window is a declared input and the two reduction sets are not nested.** Seven states the
   plate read must integrate are algebra in a fed-batch; the one state a fed-batch must integrate
   is frozen on a plate. Neither set contains the other.

3. **Window-awareness earns its place under criterion (a).** Applying the fed-batch reduction set
   to a plate read eliminates reporter dilution and biases the time-integral of the reporter by
   51.7-62.7%, against a MEASURED 14.6% floor -- 3.5 to 4.3 times the noise of the assay that
   would measure it.

4. **The ratio test is necessary and not sufficient.** Crz1 passes it by a factor of 18 and the
   audit still refuses, because its frequency is the signal.

Nothing here needs the wet-lab exports except the one class that re-derives the growth band from
the committed gate-1 tables, and that skips when they are absent.
"""

from __future__ import annotations

import functools
import math
import re

import numpy as np
import pytest

from ystwin import paths
from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE, OBSERVED_ACTIVITY_CV
from ystwin.mech.integrate import (
    ATOL,
    FALLBACK_METHOD,
    GROWTH_FLOOR_PER_H,
    PINNED_METHOD,
    REDUCTION_RATIO_CEILING,
    RETAINED_RATIO_CEILING,
    RTOL,
    Justification,
    ReductionRefused,
    Verdict,
    audit_timescales,
    freeze_bias,
    integrate_window,
    qss_bias,
    reduce_for_window,
    require_reducible,
)
from ystwin.mech.state import (
    FEDBATCH_5D,
    PLATE_GROWTH_BAND_PER_H,
    PLATE_READ_4H,
    TIMESCALE_CATALOGUE,
    Encoding,
    MechState,
    StateUnidentifiable,
    StateVar,
    TimescaleUnmeasured,
    Window,
    catalogue_span_hours,
)

CATALOGUE = TIMESCALE_CATALOGUE

# Petelenz-Kurdziel 2013, PLoS Comput Biol 9(6):e1003084, PMID 23762021: after 0.4 M NaCl
# glycerol is the dominant glycolytic sink by ~15 min and Fps1 has reopened by 30 min.
OSMOTIC_EXCURSION_H = 0.5


def _one_state(name: str) -> MechState:
    return MechState((CATALOGUE[name],))


class TestTheStateVectorRefusesRatherThanInventing:
    """Criterion (b), enforced at construction rather than in review."""

    def test_a_state_with_no_constraint_and_no_sweep_axis_is_refused(self):
        with pytest.raises(StateUnidentifiable) as excinfo:
            StateVar(name="ghost", units="a.u.", source="nowhere", tau_h=(1.0, 1.0))

        assert "sweep axis" in str(excinfo.value)

    def test_a_declared_sweep_axis_is_enough(self):
        assert StateVar(name="axis", units="a.u.", source="s", tau_h=(1.0, 1.0),
                        sweep_axis=True).sweep_axis

    def test_a_refused_time_constant_raises_on_float_and_names_the_nearest_source(self):
        """The idiom: an object whose ``float()`` is a traceback naming the missing
        measurement, not a plausible number in a table."""
        refused = CATALOGUE["snf1_p"].tau_h

        with pytest.raises(TimescaleUnmeasured) as excinfo:
            float(refused)

        assert "Snf1" in str(excinfo.value)
        assert "PMID 11782433" in str(excinfo.value)
        assert "WRONG kinase" in str(excinfo.value)

    def test_asking_a_refused_state_for_its_ratio_raises_rather_than_defaulting(self):
        with pytest.raises(TimescaleUnmeasured):
            CATALOGUE["ire1_cluster"].taus_in(PLATE_READ_4H)

    def test_exactly_one_of_a_fixed_interval_and_a_growth_multiple(self):
        with pytest.raises(ValueError):
            StateVar(name="both", units="a.u.", source="s", sweep_axis=True,
                     tau_h=(1.0, 1.0), tau_growth_multiple=1.0)
        with pytest.raises(ValueError):
            StateVar(name="neither", units="a.u.", source="s", sweep_axis=True)

    def test_an_unordered_interval_is_an_error(self):
        with pytest.raises(ValueError):
            StateVar(name="backwards", units="a.u.", source="s", sweep_axis=True,
                     tau_h=(2.0, 1.0))


class TestTheStateVectorIsNamedOrderedAndUnitTagged:
    def test_pack_and_unpack_round_trip(self):
        system = MechState((CATALOGUE["biomass"], CATALOGUE["reporter_mature"]))
        values = {"biomass": 0.5, "reporter_mature": 12.0}

        assert system.unpack(system.pack(values)) == values

    def test_pack_refuses_a_missing_name(self):
        with pytest.raises(KeyError):
            _one_state("biomass").pack({})

    def test_pack_refuses_a_name_that_is_not_a_state(self):
        with pytest.raises(KeyError):
            _one_state("biomass").pack({"biomass": 1.0, "smell": 2.0})

    def test_duplicate_names_are_refused(self):
        with pytest.raises(ValueError):
            MechState((CATALOGUE["biomass"], CATALOGUE["biomass"]))

    def test_every_catalogued_state_carries_units_and_a_source(self):
        for v in CATALOGUE:
            assert v.units, v.name
            assert v.source, v.name

    def test_every_time_constant_names_a_pmid_a_bnid_or_its_own_arithmetic(self):
        """House rule 3, as a gate. A number attributed to a paper that does not contain it is
        the worst outcome available here, so every row has to point somewhere checkable."""
        pattern = re.compile(r"PMID \d{6,8}|BNID \d{4,6}|REFUSED|DEFINITION|MEASURED per window|"
                             r"Genes Dev 9:1559|eLife 3:e05031")
        unsourced = [v.name for v in CATALOGUE if not pattern.search(v.source)]

        assert unsourced == []

    def test_every_state_is_constrained_by_an_assay_or_declared_as_an_axis(self):
        for v in CATALOGUE:
            assert bool(v.constrained_by) != v.sweep_axis or v.sweep_axis, v.name


class TestTheStiffnessIsMeasuredNotAsserted:
    """The argument for a stiff driver, computed rather than quoted from the design doc."""

    def test_the_catalogue_spans_eight_orders_of_magnitude(self):
        fastest, slowest = catalogue_span_hours()

        assert fastest == pytest.approx(1.0 / 160.0 / 3600.0)     # Ypd1 -> Ssk1, 160 /s
        assert slowest == pytest.approx(193.0, rel=1e-3)          # plasmid loss at mu = 0.101
        assert slowest / fastest == pytest.approx(1.11e8, rel=0.02)

    def test_an_explicit_method_costs_four_orders_more_than_the_pinned_one(self):
        stiff, explicit = _stiff_against_explicit()

        assert stiff.n_rhs_evals < 1_000
        assert explicit.n_rhs_evals > 1_000_000
        assert explicit.n_rhs_evals / stiff.n_rhs_evals > 1e4

    def test_the_two_methods_agree_so_the_cost_is_not_accuracy(self):
        """This is what makes it stiffness rather than a hard problem: the step size is set by
        chemistry the observable cannot see, and the answers are the same."""
        stiff, explicit = _stiff_against_explicit()

        assert stiff.time_integral("reporter_mature") == pytest.approx(
            explicit.time_integral("reporter_mature"), rel=1e-5)


class TestTheDriverIsPinned:
    def test_the_method_is_bdf_with_lsoda_behind_it(self):
        assert (PINNED_METHOD, FALLBACK_METHOD) == ("BDF", "LSODA")

    def test_the_default_run_reports_which_method_actually_ran(self):
        system, rhs, y0 = _phosphorelay_and_reporter()

        assert integrate_window(rhs, y0, PLATE_READ_4H, system).method == PINNED_METHOD

    def test_tightening_the_tolerance_a_hundredfold_does_not_move_the_observable(self):
        """RTOL and ATOL are ASSERTED numerical settings, not model parameters. This is the
        check that says so: they are pinned for reproducibility, not fitted for an answer."""
        system = _one_state("reporter_mature")
        rhs = _step_driven_reporter(mu=0.3)

        loose = integrate_window(rhs, {"reporter_mature": 0.0}, PLATE_READ_4H, system)
        tight = integrate_window(rhs, {"reporter_mature": 0.0}, PLATE_READ_4H, system,
                                 rtol=RTOL / 100, atol=ATOL / 100)

        moved = abs(tight.time_integral("reporter_mature")
                    - loose.time_integral("reporter_mature"))
        assert moved / loose.time_integral("reporter_mature") < OBSERVED_ACTIVITY_CV / 1000

    def test_a_per_state_tolerance_may_be_given_by_name(self):
        system, rhs, y0 = _phosphorelay_and_reporter()

        run = integrate_window(rhs, y0, PLATE_READ_4H, system,
                               atol={"ssk1_p": 1e-12, "reporter_mature": 1e-8})

        assert run.method == PINNED_METHOD

    def test_a_tolerance_missing_a_state_is_an_error(self):
        system, rhs, y0 = _phosphorelay_and_reporter()

        with pytest.raises(KeyError):
            integrate_window(rhs, y0, PLATE_READ_4H, system, atol={"ssk1_p": 1e-12})


class TestNonNegativityIsEnforcedInsideTheRightHandSide:
    """Clipping the state would make the trajectory discontinuous and defeat the implicit
    method's Jacobian. Clipping the argument keeps the flow smooth."""

    def test_the_right_hand_side_never_sees_a_negative_state(self):
        seen: list[float] = []
        system = _one_state("reporter_mature")

        def rhs(t, y):
            seen.append(float(y[0]))
            return [-5.0]

        integrate_window(rhs, {"reporter_mature": 1.0}, PLATE_READ_4H, system)

        assert min(seen) >= 0.0

    def test_it_can_be_turned_off_and_then_the_state_goes_negative(self):
        seen: list[float] = []
        system = _one_state("reporter_mature")

        def rhs(t, y):
            seen.append(float(y[0]))
            return [-5.0]

        integrate_window(rhs, {"reporter_mature": 1.0}, PLATE_READ_4H, system, nonnegative=False)

        assert min(seen) < 0.0


class TestTheRatioTestIsNecessaryAndNotSufficient:
    """Cai, Dalal & Elowitz 2008, PMID 18818649, MEASURED: calcium sets the burst frequency and
    not its duration, so the mean of the train does not determine the output."""

    def test_crz1_passes_the_ratio_by_a_wide_margin(self):
        row = audit_row("crz1_nuclear", PLATE_READ_4H)

        assert row.ratio_high < REDUCTION_RATIO_CEILING / 10

    def test_and_the_audit_refuses_it_anyway(self):
        assert audit_row("crz1_nuclear", PLATE_READ_4H).verdict is Verdict.FREQUENCY_ENCODED

    def test_the_refusal_says_why_and_names_the_measurement(self):
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["crz1_nuclear"], PLATE_READ_4H)

        message = str(excinfo.value)
        assert "frequency-encoded" in message
        assert "PMID 18818649" in message

    def test_msn2_is_refused_for_the_same_reason_on_a_different_measurement(self):
        """Hao & O'Shea 2012, PMID 22179789: same TF, frequency-modulated under glucose
        limitation and amplitude-modulated under oxidative stress. The stressor picks."""
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["msn2_nuclear"], FEDBATCH_5D)

        assert "PMID 22179789" in str(excinfo.value)

    def test_a_mean_field_reduction_is_the_only_route_and_must_carry_a_reason(self):
        crz1 = CATALOGUE["crz1_nuclear"]

        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(crz1, PLATE_READ_4H, Justification.MEAN_FIELD)

        assert "requires a reason" in str(excinfo.value)

    def test_a_declared_mean_field_reduction_is_granted_and_recorded_as_such(self):
        row = require_reducible(
            CATALOGUE["crz1_nuclear"], PLATE_READ_4H, Justification.MEAN_FIELD,
            reason="the reporter integrates over 4.14 h and the burst train is stationary here")

        assert row.verdict is Verdict.FREQUENCY_ENCODED

    def test_mean_field_is_not_a_justification_for_a_level_encoded_state(self):
        """It is a claim about a carrier. Offering it for a relaxing state is a category error
        and would let any state be reduced by asserting the wrong thing about it."""
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["reporter_mature"], PLATE_READ_4H,
                              Justification.MEAN_FIELD, reason="because")

        assert "level-encoded" in str(excinfo.value)

    def test_every_granted_reduction_records_which_justification_it_rests_on(self):
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            for row in audit_timescales(CATALOGUE, window):
                if row.reduces:
                    assert row.justification is Justification.RATIO, row.state.name
                else:
                    assert row.justification is None, row.state.name


class TestTheAuditRefusesABadReduction:
    def test_it_refuses_a_state_inside_the_irreducible_band(self):
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["reporter_mature"], PLATE_READ_4H)

        message = str(excinfo.value)
        assert "irreducible band" in message
        assert "14.6%" in message

    def test_it_refuses_a_straddling_interval_rather_than_picking_an_end(self):
        """mRNA half-life spans 9x across methods (Chan 2018's 3.6 min against a 32 min
        transcription shutoff). One end licenses the reduction and the other refuses it."""
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["mrna"], PLATE_READ_4H)

        message = str(excinfo.value)
        assert "PMID 30192227" in message and "BNID 114181" in message

    def test_the_reporter_maturation_straddle_names_both_yeast_papers(self):
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["reporter_immature"], PLATE_READ_4H)

        message = str(excinfo.value)
        assert "PMID 35180343" in message and "PMID 17237792" in message

    def test_it_refuses_a_state_whose_time_constant_nobody_has_measured(self):
        with pytest.raises(ReductionRefused) as excinfo:
            require_reducible(CATALOGUE["slt2_p"], PLATE_READ_4H)

        assert "no measured time constant" in str(excinfo.value)

    def test_it_grants_a_state_that_actually_passes(self):
        row = require_reducible(CATALOGUE["ssk1_p"], PLATE_READ_4H)

        assert (row.verdict, row.justification) == (Verdict.ELIMINATE, Justification.RATIO)

    def test_a_straddling_state_is_integrated_rather_than_dropped(self):
        """The safe direction when the criterion cannot decide is to keep the equation."""
        plate = reduce_for_window(CATALOGUE, PLATE_READ_4H)

        assert "mrna" in plate.integrated.names
        assert "reporter_immature" in plate.integrated.names


class TestTheIrreducibleSetsAreNotNested:
    """The structural finding, and it is arithmetic on measured time constants."""

    def test_neither_window_set_contains_the_other(self):
        plate = set(reduce_for_window(CATALOGUE, PLATE_READ_4H).integrated.names)
        fedbatch = set(reduce_for_window(CATALOGUE, FEDBATCH_5D).integrated.names)

        assert not plate <= fedbatch
        assert not fedbatch <= plate

    def test_the_plate_integrates_seven_states_the_fedbatch_reduces_to_algebra(self):
        plate = reduce_for_window(CATALOGUE, PLATE_READ_4H)
        fedbatch = reduce_for_window(CATALOGUE, FEDBATCH_5D)

        only_on_the_plate = set(plate.integrated.names) - set(fedbatch.integrated.names)

        assert only_on_the_plate == {
            "cross_protection", "glutathione", "mrna", "reporter_immature",
            "reporter_mature", "trehalose", "upr_output",
        }

    def test_plasmid_loss_is_frozen_on_a_plate_and_integrated_in_a_fedbatch(self):
        """The one process that flips from negligible to dominant between the two vessels."""
        assert audit_row("plasmid_bearing", PLATE_READ_4H).verdict is Verdict.FREEZE
        assert audit_row("plasmid_bearing", FEDBATCH_5D).verdict is Verdict.KEEP

    def test_cross_protection_the_marquee_emergence_result_is_algebra_in_a_fedbatch(self):
        """Recorded because it is a conclusion the architecture states the premise of and never
        draws: on a 5-day run the cross-protection state reduces away under its own rule."""
        assert audit_row("cross_protection", PLATE_READ_4H).verdict is Verdict.KEEP
        assert audit_row("cross_protection", FEDBATCH_5D).verdict is Verdict.ELIMINATE

    def test_the_fedbatch_keeps_only_the_population_block(self):
        fedbatch = reduce_for_window(CATALOGUE, FEDBATCH_5D)
        measured = {r.state.name for r in fedbatch.audit
                    if r.verdict is not Verdict.REFUSED
                    and r.state.encoding is not Encoding.FREQUENCY
                    and r.state.name in fedbatch.integrated.names}

        assert measured == {"biomass", "product", "generations", "plasmid_bearing"}

    def test_an_integrating_state_is_never_eliminated_because_it_has_no_steady_state(self):
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            for row in audit_timescales(CATALOGUE, window):
                if row.state.encoding is Encoding.INTEGRATING:
                    assert row.verdict is not Verdict.ELIMINATE, (row.state.name, window.name)


class TestWindowAwarenessEarnsItsPlace:
    """Criterion (a). Removing the piece must move a bioproduction-relevant observable by more
    than the noise of the assay that would measure it."""

    def test_using_the_fedbatch_set_on_a_plate_biases_the_reporter_past_its_floor(self):
        """ABLATION. The fed-batch set eliminates reporter dilution (tau/T = 0.083 at mu = 0.101);
        the plate set keeps it (tau/T = 0.67-0.99). Do it anyway and the time-integral of the
        reporter -- what a stable FP records over the read -- moves by 3.5 to 4.3 times the
        MEASURED 14.6% plate CV. Observable: reporter activity. Floor:
        panel_experiment.OBSERVED_ACTIVITY_CV."""
        system = _one_state("reporter_mature")
        biases = []

        for mu in PLATE_GROWTH_BAND_PER_H:
            kept = integrate_window(_step_driven_reporter(mu), {"reporter_mature": 0.0},
                                    PLATE_READ_4H, system)
            true_integral = kept.time_integral("reporter_mature")
            eliminated_integral = (1.0 / mu) * PLATE_READ_4H.duration_h
            biases.append((eliminated_integral - true_integral) / eliminated_integral)

        assert biases == pytest.approx([0.6270, 0.5166], abs=5e-4)
        assert min(biases) / OBSERVED_ACTIVITY_CV > 3.5

    def test_the_closed_form_agrees_with_the_integration(self):
        system = _one_state("reporter_mature")
        for mu in PLATE_GROWTH_BAND_PER_H:
            kept = integrate_window(_step_driven_reporter(mu), {"reporter_mature": 0.0},
                                    PLATE_READ_4H, system)
            eliminated_integral = PLATE_READ_4H.duration_h / mu
            integrated_bias = 1.0 - kept.time_integral("reporter_mature") / eliminated_integral
            assert qss_bias(1.0 / mu, PLATE_READ_4H.duration_h) == pytest.approx(
                integrated_bias, abs=5e-4)

    def test_using_the_plate_set_on_a_fedbatch_erases_a_quarter_of_the_plasmid_loss(self):
        """The reverse direction. The plate set freezes the plasmid-bearing fraction; over 120 h
        at the setpoint the culture runs 17.5 generations and 46% of it goes plasmid-free.

        Scored as an EFFECT SIZE and not against a floor, deliberately. The assay is replica
        plating on selective against non-selective medium, and this repository has no measured
        CV for it -- borrowing the plate reader's activity CV for a colony count is the exact
        substitution that scoring rule exists to prevent."""
        system = _one_state("plasmid_bearing")
        tau = CATALOGUE["plasmid_bearing"].taus_in(FEDBATCH_5D)[0]

        run = integrate_window(lambda t, y: [-y[0] / tau], {"plasmid_bearing": 1.0},
                               FEDBATCH_5D, system)

        assert run.time_average("plasmid_bearing") == pytest.approx(0.7447, abs=1e-3)
        assert run.of("plasmid_bearing")[-1] == pytest.approx(0.5370, abs=1e-3)
        assert freeze_bias(tau, FEDBATCH_5D.duration_h) == pytest.approx(0.2553, abs=1e-3)


class TestTheBiasFormulaIsExactAndTheRatioBoundsIt:
    def test_the_quasi_steady_bias_matches_a_numerical_integration(self):
        tau, window, mu = 2.0, 4.14, 0.5
        system = _one_state("reporter_mature")

        run = integrate_window(lambda t, y: [(1.0 - y[0]) / tau], {"reporter_mature": 0.0},
                               Window("probe", window, mu, mu, "probe"), system)
        numerical = (window - run.time_integral("reporter_mature")) / window

        assert numerical == pytest.approx(qss_bias(tau, window), rel=1e-5)

    def test_freezing_and_eliminating_are_the_two_halves_of_one_integral(self):
        assert qss_bias(3.0, 4.14) + freeze_bias(3.0, 4.14) == pytest.approx(1.0)

    def test_the_ratio_the_criterion_tests_is_an_upper_bound_on_the_bias(self):
        """So the criterion is conservative, and the audit's bias column says by how much."""
        for tau in (0.01, 0.1, 0.5, 1.0, 4.0):
            assert qss_bias(tau, 4.14) <= tau / 4.14 + 1e-12

    def test_a_non_positive_window_is_an_error(self):
        with pytest.raises(ValueError):
            qss_bias(1.0, 0.0)


class TestADecayingDriverChangesTheDenominator:
    """The tau/T bias formula is derived for a STEP input. When a source measures how long the
    driver actually persists and that is shorter than the run, the run length is the wrong
    denominator and the audit scores against the driver instead."""

    def test_hog1_is_eliminated_against_the_whole_read(self):
        assert audit_row("hog1_p", PLATE_READ_4H).verdict is Verdict.ELIMINATE

    def test_and_kept_against_the_measured_osmotic_excursion(self):
        """Petelenz-Kurdziel 2013, PMID 23762021: glycerol is the dominant sink by ~15 min and
        Fps1 has reopened by 30 min, so the driver is gone long before the read ends."""
        row = audit_row("hog1_p", PLATE_READ_4H,
                        driver_persists_h={"hog1_p": OSMOTIC_EXCURSION_H})

        assert row.effective_window_h == OSMOTIC_EXCURSION_H
        assert row.verdict is Verdict.KEEP

    def test_nothing_moves_when_no_driver_is_declared(self):
        plain = audit_timescales(CATALOGUE, PLATE_READ_4H)
        empty = audit_timescales(CATALOGUE, PLATE_READ_4H, driver_persists_h={})

        assert [r.verdict for r in plain] == [r.verdict for r in empty]

    def test_a_driver_longer_than_the_window_does_not_extend_it(self):
        row = audit_row("hog1_p", PLATE_READ_4H, driver_persists_h={"hog1_p": 999.0})

        assert row.effective_window_h == PLATE_READ_4H.duration_h

    def test_a_non_positive_driver_duration_is_an_error(self):
        with pytest.raises(ValueError):
            audit_row("hog1_p", PLATE_READ_4H, driver_persists_h={"hog1_p": 0.0})


class TestTheAuditTableSaysWhatItCannotDecide:
    def test_it_states_that_the_ratio_test_is_not_sufficient(self):
        table = reduce_for_window(CATALOGUE, PLATE_READ_4H).audit_table()

        assert "NECESSARY AND NOT SUFFICIENT" in table
        assert "PMID 18818649" in table

    def test_it_states_that_a_systematic_bias_is_not_a_random_cv(self):
        table = reduce_for_window(CATALOGUE, FEDBATCH_5D).audit_table()

        assert "SYSTEMATIC" in table and "RANDOM" in table

    def test_it_names_the_window_and_the_growth_band_it_was_scored_at(self):
        table = reduce_for_window(CATALOGUE, FEDBATCH_5D).audit_table()

        assert "FEDBATCH_5D" in table and "T = 120 h" in table and "0.101" in table

    def test_every_state_gets_a_row(self):
        plate = reduce_for_window(CATALOGUE, PLATE_READ_4H)

        assert len(plate.audit) == len(CATALOGUE)
        assert {r.state.name for r in plate.audit} == set(CATALOGUE.names)

    def test_the_three_dispositions_partition_the_system(self):
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            reduced = reduce_for_window(CATALOGUE, window)
            parts = (set(reduced.integrated.names) | set(reduced.eliminated.names)
                     | set(reduced.frozen.names))

            assert parts == set(CATALOGUE.names)
            assert (len(reduced.integrated) + len(reduced.eliminated)
                    + len(reduced.frozen)) == len(CATALOGUE)


class TestTheFloorsAreImportedNotRestated:
    """A criterion that no longer matches the instrument it was derived from is worse than none,
    so the digits live in one place."""

    def test_the_reduction_criterion_is_the_measured_plate_cv(self):
        assert REDUCTION_RATIO_CEILING is OBSERVED_ACTIVITY_CV
        assert REDUCTION_RATIO_CEILING == 0.146

    def test_the_other_half_of_the_inequality_carries_no_new_number(self):
        assert RETAINED_RATIO_CEILING == pytest.approx(1.0 / 0.146)

    def test_the_growth_floor_is_the_measured_growth_standard_error(self):
        assert GROWTH_FLOOR_PER_H is MEASURED_GROWTH_RATE_SE
        assert GROWTH_FLOOR_PER_H == 0.0117


class TestTheGrowthBandIsMeasuredFromTheRealPlates:
    """House rule: mu is measured off the biomass series, never asserted from a pump. This
    re-derives the band from the committed gate-1 tables and fails if they move."""

    def test_the_plate_band_is_the_iqr_of_the_committed_gate_one_tables(self):
        pandas = pytest.importorskip("pandas")
        tables = sorted((paths.REPO_ROOT / "outputs").glob("g1_2026*.csv"))
        wanted = [t for t in tables if "NewProtocol" in t.name or "Replicate3" in t.name]
        if len(wanted) != 3:
            pytest.skip("the three NewProtocol gate-1 tables are not in this checkout")

        rates = pandas.concat([
            np.log(pandas.read_csv(t).query("passed")["od_fold_change"].astype(float))
            / PLATE_READ_4H.duration_h
            for t in wanted])

        assert len(rates) == 183
        assert (rates.quantile(0.25), rates.quantile(0.75)) == pytest.approx(
            PLATE_GROWTH_BAND_PER_H, rel=1e-12)

    def test_the_fedbatch_setpoint_is_the_state_the_vessel_reproduces(self):
        """0.254 /h is the upper Elizondo state and also that strain's measured mu_max, so
        `fba/fedbatch.py` refuses it: at the ceiling the feed no longer sets mu."""
        from ystwin.fba.fedbatch import REFUSES_THE_UPPER_ELIZONDO_STATE

        assert FEDBATCH_5D.growth_rate_low_per_h == 0.101
        assert "0.254" in REFUSES_THE_UPPER_ELIZONDO_STATE

    def test_the_plate_window_is_the_length_of_the_real_read(self):
        from ystwin.generator.plate import PlateConditions

        assert PLATE_READ_4H.duration_h == PlateConditions().duration_h


class TestTheWindowRefusesAnImpossibleDeclaration:
    def test_a_zero_length_window_is_an_error(self):
        with pytest.raises(ValueError):
            Window("nothing", 0.0, 0.3, 0.3, "s")

    def test_an_unordered_growth_band_is_an_error(self):
        with pytest.raises(ValueError):
            Window("backwards", 4.0, 0.4, 0.3, "s")


def audit_row(name, window, **kwargs):
    """The audit row for one catalogued state, for the tests that want just one."""
    rows = audit_timescales(_one_state(name), window, **kwargs)
    return rows[0]


def _step_driven_reporter(mu):
    """dR/dt = 1 - mu*R: a promoter switched on at t = 0 feeding a stable FP diluted by growth."""
    return lambda t, y: [1.0 - mu * y[0]]


def _phosphorelay_and_reporter():
    """The two-state system the stiffness argument is made on.

    Ypd1 -> Ssk1 at a MEASURED 160 /s (Janiak-Spens 2005, PMID 15628880) driving a reporter that
    dilutes at mu. Nothing about the reporter can see 6.25 ms; an explicit method's step size
    can see nothing else.
    """
    system = MechState((CATALOGUE["ssk1_p"], CATALOGUE["reporter_mature"]))
    tau_fast = CATALOGUE["ssk1_p"].taus_in(PLATE_READ_4H)[0]
    mu = 0.3

    def rhs(t, y):
        return [(0.5 - y[0]) / tau_fast, y[0] - mu * y[1]]

    return system, rhs, {"ssk1_p": 0.0, "reporter_mature": 0.0}


@functools.lru_cache(maxsize=1)
def _stiff_against_explicit():
    """The pinned method and RK45 on the same system, computed once.

    Cached because the explicit run really does cost five million right-hand-side evaluations,
    which is the finding -- and paying for it twice would put 37 s on the suite to re-measure a
    number that does not change between two assertions about it.
    """
    system, rhs, y0 = _phosphorelay_and_reporter()
    return (integrate_window(rhs, y0, PLATE_READ_4H, system, method=PINNED_METHOD),
            integrate_window(rhs, y0, PLATE_READ_4H, system, method="RK45"))


def test_the_module_docstring_numbers_are_the_ones_the_code_produces():
    """The claim in the module docstring above -- 1.1e8 span, 51.7-62.7% ablation -- is the
    output of the two tests above. This asserts they have not drifted apart in prose."""
    fastest, slowest = catalogue_span_hours()

    assert math.log10(slowest / fastest) == pytest.approx(8.05, abs=0.05)
    low, high = PLATE_GROWTH_BAND_PER_H
    assert qss_bias(1.0 / high, PLATE_READ_4H.duration_h) == pytest.approx(0.517, abs=1e-3)
    assert qss_bias(1.0 / low, PLATE_READ_4H.duration_h) == pytest.approx(0.627, abs=1e-3)
