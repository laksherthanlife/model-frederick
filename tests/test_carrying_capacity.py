"""The generator's carrying capacity is a function of the dose, and what that did not fix.

`simulate_culture` integrated every well against `params.carrying_capacity` -- one asserted
number per strain that no dose could move -- while `culture.substrate_limited_capacity` sat
beside it deriving the same quantity from van Hoek 1998's measured biomass yield and was
imported by nothing but `tests/test_crabtree_transition.py`. The capacity is now derived by
default, from the glucose charge and the yield at the realised rate, so a dose moves it.

Half of this file pins that wiring. The other half pins the result of scoring it, because the
result is a negative and negatives are what rot. `REVISED_BUILD_LIST.md` ranks this change
first on the whole build list on the strength of a 203-well high-dose R_late/R_early violation
that "a logistic with a dose-dependent K reproduces 0.593 of", with a fitted
K_dosed/K_control of 0.245 at 2-5 mM. The measured yield curve cannot produce a capacity ratio
below 1 at any glucose charge: van Hoek's yield *rises* as the rate falls, so slowing a
culture makes its substrate budget bigger. Scored on the committed late/early table the fix
closes about 5% of the gap, not 59%, and the sign of the dose effect is wrong below 20 g/L.

That is not an argument against the change -- the change removes four asserted scalars, adds
one declared one, and clears the growth assay's own noise floor. It is an argument against the
reason the build list gave for it, and it belongs in the suite rather than in a report.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pytest

from ystwin.generator.culture import (
    STANDARD_GLUCOSE_G_PER_L,
    CultureParameters,
    measured_growth_rate_range,
    resolved_capacity,
    simulate_culture,
    substrate_limited_capacity,
)
from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE
from ystwin.generator.plate import DEFAULT_PANEL

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "carrying_capacity_check.py"

# The default plate's inoculum in g/L: 0.16 OD at 0.42 gDCW/OD, both from PlateConditions.
INOCULUM_G_PER_L = 0.0672

FAST = CultureParameters(
    mu_max=0.62, growth_ic50=0.8, growth_hill=1.5,
    promoter_basal=1.0e-3, promoter_peak=2.4e-3, promoter_ec50=0.9, promoter_hill=2.0,
    carrying_capacity=1.4,
)
"""A strain whose unstressed rate is outside the curve van Hoek measured.

0.62 /h is not exotic -- it is the rate `tests/test_generator_culture.py` has always used,
and the median `mu_max` fitted to the real wells is 0.567 /h. Every one of those is above the
0.40 /h ceiling of the measured yield, which is why the fallback is a live path and not a
corner case.
"""


@pytest.fixture(scope="module")
def check():
    spec = importlib.util.spec_from_file_location("carrying_capacity_check", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheCapacityIsDerivedFromTheMeasuredYield:
    def test_a_dose_moves_the_capacity_the_well_is_integrated_against(self):
        """The whole defect in one assertion: the incumbent returned 3.0 at every dose."""
        params = DEFAULT_PANEL["UPRE1"]
        times = np.linspace(0.0, 4.14, 25)

        control = simulate_culture(times, 0.0, params, INOCULUM_G_PER_L)
        dosed = simulate_culture(times, 5.0, params, INOCULUM_G_PER_L)

        assert control["carrying_capacity"] != dosed["carrying_capacity"]
        assert control["carrying_capacity"] != params.carrying_capacity

    def test_it_is_exactly_the_charge_times_the_measured_yield(self):
        """Derived, not fitted, and identical to the orphan it wires in.

        Written as an identity rather than a number so that a change to the vendored van Hoek
        table moves both sides at once instead of leaving a stale literal here.
        """
        params = DEFAULT_PANEL["UPRE1"]
        rate = params.growth_rate_at(2.0)

        capacity, source = resolved_capacity(params, rate, INOCULUM_G_PER_L)

        assert capacity == substrate_limited_capacity(
            INOCULUM_G_PER_L, STANDARD_GLUCOSE_G_PER_L, rate)
        assert "van Hoek" in source and "9797269" in source

    def test_the_numbers_the_default_panel_now_runs_at(self):
        """Pinned because every synthetic plate in the repository moved when they changed.

        UPRE1 was integrated against 3.0 g/L at every rung and is now integrated against
        4.55 at 0 mM and 9.77 at 5 mM. A reader comparing an old table to a new one needs
        the two numbers written down somewhere that fails when they drift.
        """
        params = DEFAULT_PANEL["UPRE1"]

        control, _ = resolved_capacity(params, params.growth_rate_at(0.0), INOCULUM_G_PER_L)
        top, _ = resolved_capacity(params, params.growth_rate_at(5.0), INOCULUM_G_PER_L)

        assert control == pytest.approx(4.5472, abs=1e-3)
        assert top == pytest.approx(9.7732, abs=1e-3)

    def test_the_result_carries_which_capacity_it_used(self):
        """Two routes that disagree by more than a factor of two need a label, not a value.

        `scripts/run_environment_sweep.py` already wrote `carrying_capacity_source` into every
        row it produced; this is the same bookkeeping moved to where the choice is made.
        """
        out = simulate_culture(np.linspace(0.0, 4.14, 25), 0.0, DEFAULT_PANEL["UPRE1"],
                               INOCULUM_G_PER_L)

        assert isinstance(out["capacity_source"], str)
        assert out["carrying_capacity"] > 0.0


class TestTheAssertedCapacityIsTheDeclaredOptOut:
    """A capacity-limited well is a real thing to want, so the old behaviour stays reachable."""

    def test_declining_the_substrate_route_restores_the_asserted_capacity(self):
        params = DEFAULT_PANEL["UPRE1"]
        times = np.linspace(0.0, 4.14, 25)

        out = simulate_culture(times, 5.0, params, INOCULUM_G_PER_L, glucose_g_per_L=None)

        assert out["carrying_capacity"] == params.carrying_capacity
        assert "declined" in out["capacity_source"]

    def test_the_opt_out_still_saturates_where_it_always_did(self):
        """The pre-fix behaviour, unchanged: a long read closes on `carrying_capacity`."""
        times = np.linspace(0.0, 40.0, 200)

        out = simulate_culture(times, 0.0, FAST, 0.05, glucose_g_per_L=None)

        assert out["biomass"][-1] == pytest.approx(FAST.carrying_capacity, rel=0.05)

    def test_a_substrate_limited_well_ends_up_somewhere_else_entirely(self):
        """Which is why the opt-out has to be explicit rather than a default nobody notices."""
        times = np.linspace(0.0, 40.0, 200)
        params = DEFAULT_PANEL["UPRE1"]

        limited = simulate_culture(times, 0.0, params, INOCULUM_G_PER_L)
        capped = simulate_culture(times, 0.0, params, INOCULUM_G_PER_L, glucose_g_per_L=None)

        assert limited["biomass"][-1] > capped["biomass"][-1] * 1.4


class TestARefusedYieldFallsBackAndSaysWhy:
    """`chemostat_physiology` refuses outside 0.025-0.40 /h, and a refusal is not a licence."""

    @pytest.mark.parametrize("rate", [0.62, 1.2, 0.41, 0.01])
    def test_a_rate_off_the_measured_curve_takes_the_asserted_capacity(self, rate):
        capacity, source = resolved_capacity(FAST, rate, 0.05)

        assert capacity == FAST.carrying_capacity
        assert "outside the measured yield curve" in source

    def test_the_fallback_names_the_range_it_left(self):
        """So the reader's next move is obvious: stay inside it, or measure more of it."""
        low, high = measured_growth_rate_range()

        _, source = resolved_capacity(FAST, 0.62, 0.05)

        assert f"{low:g}" in source and f"{high:g}" in source

    @pytest.mark.parametrize("charge, inoculum", [(0.0, 0.05), (-1.0, 0.05), (20.0, 0.0)])
    def test_a_bad_charge_or_inoculum_raises_instead_of_falling_back(self, charge, inoculum):
        """The fallback is for a refused *rate* only, and the distinction is load-bearing.

        A capacity that quietly became the asserted one because the charge was zero would be
        the asserted capacity wearing the measured route's label, which is exactly the kind of
        silent substitution the fallback's own message exists to prevent.
        """
        with pytest.raises(ValueError, match="positive charge and inoculum"):
            resolved_capacity(FAST, 0.30, inoculum, charge)

    def test_the_fallback_is_not_a_corner_case_on_the_real_plates(self, check):
        """91% of the committed conditions grow faster than van Hoek's fastest chemostat.

        Worth an assertion because it bounds what this fix can ever do to a *calibrated*
        panel: fitted to the real plates, most wells never reach the substrate route at all.
        """
        observed = check.measured()
        if observed is None:
            pytest.skip("outputs/sensor_characterisation.csv is absent")
        _, ceiling = measured_growth_rate_range()

        above = float((observed.mu_max > ceiling).mean())

        assert above > 0.85
        assert observed.mu_max.median() > ceiling


class TestTheSignTheYieldCurveGivesIsNotTheSignTheFitsWant:
    """The load-bearing refutation, and the reason the fix cannot close what it was aimed at."""

    def test_slowing_a_culture_raises_its_substrate_budget(self):
        """van Hoek's yield goes 0.20 g/g at 0.40 /h to 0.48 at 0.25, so K goes up with dose."""
        params = DEFAULT_PANEL["UPRE1"]

        control, _ = resolved_capacity(params, params.growth_rate_at(0.0), INOCULUM_G_PER_L)
        dosed, _ = resolved_capacity(params, params.growth_rate_at(5.0), INOCULUM_G_PER_L)

        assert dosed > control
        assert dosed / control == pytest.approx(2.149, rel=0.01)

    @pytest.mark.parametrize("charge", [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 100.0])
    def test_no_glucose_charge_makes_the_capacity_fall_with_dose(self, charge):
        """`REVISED_BUILD_LIST.md` quotes a fitted K_dosed/K_control of 0.245 at 2-5 mM.

        The charge is the only free quantity in the derived capacity, and it multiplies both
        ends of the ratio, so it cannot change the ratio's side of 1. The fitted drop is
        therefore a different claim from "the measured yield curve predicts it", and this
        parametrisation is the proof rather than a paragraph asserting it.
        """
        params = DEFAULT_PANEL["UPRE1"]

        control, _ = resolved_capacity(
            params, params.growth_rate_at(0.0), INOCULUM_G_PER_L, charge)
        dosed, _ = resolved_capacity(
            params, params.growth_rate_at(2.0), INOCULUM_G_PER_L, charge)

        assert dosed / control > 1.0


class TestTheAblationClearsTheGrowthAssaysOwnFloor:
    """Criterion (a), on the growth channel and its measured floor -- not the reporter CV.

    `MEASURED_GROWTH_RATE_SE` = 0.0117 /h is the standard error of a growth rate estimated
    from one of these traces, median over 210 real wells. `OBSERVED_ACTIVITY_CV` = 0.146 is a
    reporter floor and is deliberately not used here: nothing in this comparison reads a
    reporter.
    """

    @staticmethod
    def _worst_shift(check, duration_h: float, seeds) -> float:
        worst = []
        for seed in seeds:
            on = check.simulated(STANDARD_GLUCOSE_G_PER_L, duration_h, seed)
            off = check.simulated(None, duration_h, seed)
            worst.append(float((on.mu_late - off.mu_late).abs().max()))
        return float(np.mean(worst))

    def test_the_short_protocol_moves_the_late_growth_rate_past_the_floor(self, check):
        """Barely, and the margin is the point: 4.14 h from OD 0.16 is not substrate-limited.

        The orphan's own docstring scoped the effect to a plate "that runs to substrate
        exhaustion". This one does not, so a fix aimed at it should be expected to sit near
        the floor rather than above it, and it does.
        """
        shift = self._worst_shift(check, 4.14, (0, 1, 2))

        assert shift > MEASURED_GROWTH_RATE_SE
        assert shift / MEASURED_GROWTH_RATE_SE == pytest.approx(1.23, abs=0.15)

    def test_a_read_that_reaches_the_capacity_moves_it_far_past_the_floor(self, check):
        shift = self._worst_shift(check, 24.0, (0, 1))

        assert shift / MEASURED_GROWTH_RATE_SE > 4.0


class TestItDoesNotCloseTheHighDoseDeceleration:
    """The negative half, pinned so nobody re-derives it or re-claims the 0.593.

    Scored on `outputs/sensor_characterisation.csv` rather than on the raw exports, because
    `paths.biosensor_plates()` returns None on this machine and that table is the only
    committed carrier of `mu_late` beside `mu_max`.
    """

    def test_the_measured_drop_is_what_the_committed_table_says(self, check):
        observed = check.measured()
        if observed is None:
            pytest.skip("outputs/sensor_characterisation.csv is absent")

        control, high, drop = check._bands(observed)

        # Current sensor-characterisation replay, with one row per condition (not per well).
        assert len(observed) == 112
        assert observed["band"].value_counts().to_dict() == {
            "control": 16, "mid": 64, "high": 32}
        assert control == pytest.approx(0.4511440571683947, abs=0.002)
        assert high == pytest.approx(0.1299948836096634, abs=0.002)
        assert drop == pytest.approx(0.3211491735587313, abs=0.003)

    def test_the_generator_reproduces_a_small_fraction_of_it_either_way(self, check):
        """Both capacities under-predict the dosed well's deceleration by about 0.23-0.24.

        The claim under test is `REVISED_BUILD_LIST.md`'s "reproduces 0.593 of it". It does
        not: the derived capacity closes single-digit percent of what the asserted one leaves.
        """
        observed = check.measured()
        if observed is None:
            pytest.skip("outputs/sensor_characterisation.csv is absent")
        _, _, target = check._bands(observed)

        seeds = check.SEEDS[:3]
        fixed = np.mean([check._bands(check.simulated(None, 4.14, s))[2] for s in seeds])
        derived = np.mean(
            [check._bands(check.simulated(STANDARD_GLUCOSE_G_PER_L, 4.14, s))[2] for s in seeds])

        assert abs(target - derived) > 0.20
        assert (derived - fixed) / (target - fixed) < 0.15
