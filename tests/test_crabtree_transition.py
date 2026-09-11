"""The generator's growth rate has to carry the metabolism that goes with it.

A stressor in this generator does one thing to physiology: it lowers the growth rate. In a real
aerobic glucose culture that is not one thing, it is a switch. Below a critical rate yeast
respires completely -- no ethanol, biomass yield near 0.48 g/g, oxygen uptake rising with rate.
Above it the cells overflow into ethanol, the yield collapses towards 0.20, and oxygen uptake
*falls*. van Hoek 1998 (PMID 9797269, Table 1) measured every one of those numbers across ten
dilution rates and put the switch at 0.28 /h.

The generator's untreated glucose culture grows at 0.40 /h, which is above the switch, and a
mid-range stress dose takes it below. So the switch is not an exotic corner of the sweep; it is
in the middle of every dose ladder the project runs. Before this, nothing in the generator knew
it existed: biomass yield was a per-strain constant, so a slowed culture ended up sparser when
the measured yield says it ends up denser on the same glucose.

These tests pin the curve against the vendored table, and pin the refusal that keeps it from
being extrapolated into growth rates nobody measured.
"""

import numpy as np
import pytest

from ystwin.generator.culture import (
    CRITICAL_GROWTH_RATE_PER_H,
    CultureParameters,
    chemostat_physiology,
    measured_growth_rate_range,
    simulate_culture,
    substrate_limited_capacity,
)


class TestTheMeasuredCurveIsReproduced:
    def test_the_critical_growth_rate_is_the_one_van_hoek_measured(self):
        """0.28 /h is where Table 1 first detects ethanol, not a round number chosen for it.

        `tests/test_physiology_validation.py` already asserts the critical rate is near 0.28 on
        this paper's authority. Naming it in the generator too means the two cannot drift apart
        without one of them failing.
        """
        assert CRITICAL_GROWTH_RATE_PER_H == 0.28

    def test_no_ethanol_is_produced_below_the_critical_rate(self):
        """Every dilution rate up to 0.25 in Table 1 is below the assay's detection limit.

        A generator that let a trace of ethanol leak in below the switch would make the
        transition a slope rather than a threshold, and the threshold is the finding.
        """
        for rate in (0.025, 0.05, 0.10, 0.15, 0.20, 0.25):
            assert chemostat_physiology(rate).ethanol == 0.0
            assert not chemostat_physiology(rate).fermentative

    def test_ethanol_appears_at_the_critical_rate_and_climbs_steeply(self):
        assert chemostat_physiology(0.28).ethanol == pytest.approx(0.11)
        assert chemostat_physiology(0.30).ethanol == pytest.approx(2.30)
        assert chemostat_physiology(0.40).ethanol == pytest.approx(13.90)

    def test_the_biomass_yield_collapses_across_the_switch(self):
        """0.48 to 0.20 g/g, a factor of 2.4, over a growth-rate change of 0.15 /h.

        This is the fact the generator's fixed `carrying_capacity` contradicted. It is worth an
        assertion of its own because the *direction* is counterintuitive: the faster culture
        makes less biomass per gram of sugar.
        """
        respiratory = chemostat_physiology(0.25).biomass_yield
        fermentative = chemostat_physiology(0.40).biomass_yield

        assert respiratory == pytest.approx(0.48)
        assert fermentative == pytest.approx(0.20)
        assert respiratory / fermentative == pytest.approx(2.4, rel=0.01)

    def test_oxygen_uptake_is_not_monotonic_in_growth_rate(self):
        """It peaks at the switch and then falls, which is the trap in this dataset.

        Code that assumes oxygen demand rises with growth rate is right on six of van Hoek's ten
        rows and wrong on the four that matter. The peak is 7.4 at D = 0.28 and the value at the
        highest rate measured, 0.40, is *half* of it.
        """
        rates = np.linspace(0.025, 0.40, 60)
        oxygen = np.array([chemostat_physiology(r).oxygen_uptake for r in rates])

        assert oxygen.max() == pytest.approx(7.4, rel=0.01)
        assert oxygen[-1] < oxygen.max() / 1.9
        assert not np.all(np.diff(oxygen) >= 0)

    def test_the_respiratory_quotient_stays_near_one_until_fermentation_starts(self):
        """The quantity a wet lab watches to see the switch happen.

        van Hoek reports RQ within 8% of unity below the critical rate. If this ever drifts, the
        transcription of q_CO2 or q_O2 has drifted with it.
        """
        for rate in (0.05, 0.10, 0.15, 0.20, 0.25):
            assert chemostat_physiology(rate).respiratory_quotient == pytest.approx(1.0, abs=0.1)
        assert chemostat_physiology(0.40).respiratory_quotient > 4.0

    def test_a_rate_between_measured_rows_is_interpolated_and_says_so(self):
        between = chemostat_physiology(0.32)

        assert between.interpolated
        assert 0.30 < between.growth_rate < 0.35
        assert chemostat_physiology(0.30).ethanol < between.ethanol
        assert between.ethanol < chemostat_physiology(0.35).ethanol

    def test_a_rate_that_landed_on_a_measured_row_is_not_called_interpolated(self):
        assert not chemostat_physiology(0.30).interpolated


class TestItRefusesRatherThanExtrapolating:
    """The single most tempting bad fix here is to clamp to the nearest measured row.

    A stressed culture in this generator goes to zero growth, so the low edge is hit constantly
    and clamping is silent. It would report a biomass yield of 0.45 g/g for a poisoned well, on
    the authority of a chemostat held at 0.025 /h -- a completely different physiological state
    that happens to share a growth rate with nothing.
    """

    def test_the_measured_range_is_the_range_the_paper_ran(self):
        assert measured_growth_rate_range() == (0.025, 0.40)

    @pytest.mark.parametrize("rate", [0.0, 0.01, 0.024, 0.45, 1.0])
    def test_a_growth_rate_outside_the_measured_curve_is_refused(self, rate):
        with pytest.raises(ValueError, match="outside the measured chemostat curve"):
            chemostat_physiology(rate)

    def test_the_refusal_names_the_range_and_the_source(self):
        """A refusal a caller cannot act on is an obstruction.

        The two things the caller needs are where the boundary is and who measured it, so both
        are in the message rather than in a comment somewhere else.
        """
        with pytest.raises(ValueError) as raised:
            chemostat_physiology(0.9)

        message = str(raised.value)
        assert "0.025" in message and "0.4" in message
        assert "9797269" in message

    def test_the_endpoints_themselves_are_inside_the_range(self):
        assert chemostat_physiology(0.025).biomass_yield == pytest.approx(0.45)
        assert chemostat_physiology(0.40).biomass_yield == pytest.approx(0.20)


class TestCarryingCapacityFollowsTheMeasuredYield:
    """What turns the switch from a table lookup into something the OD channel can see."""

    def test_slowing_a_culture_across_the_switch_raises_the_biomass_it_can_reach(self):
        """The sign the old fixed capacity got wrong.

        On 20 g/L glucose a culture at 0.40 /h can make about 4 g/L of biomass and one at
        0.25 /h about 9.6 g/L, because the slower one is not throwing carbon away as ethanol.
        A generator whose capacity is a per-strain constant says a stressed well ends up at the
        same density as an unstressed one, and a generator that lowered the capacity with dose
        would have the sign backwards.
        """
        fast = substrate_limited_capacity(0.07, glucose_g_per_L=20.0, growth_rate=0.40)
        slow = substrate_limited_capacity(0.07, glucose_g_per_L=20.0, growth_rate=0.25)

        assert slow > fast
        assert slow / fast == pytest.approx(2.3, rel=0.05)

    def test_it_is_the_inoculum_plus_what_the_glucose_buys(self):
        capacity = substrate_limited_capacity(0.5, glucose_g_per_L=10.0, growth_rate=0.40)

        assert capacity == pytest.approx(0.5 + 10.0 * 0.20)

    def test_it_inherits_the_refusal_for_an_unmeasured_growth_rate(self):
        with pytest.raises(ValueError, match="outside the measured chemostat curve"):
            substrate_limited_capacity(0.07, glucose_g_per_L=20.0, growth_rate=0.005)

    @pytest.mark.parametrize("bad", [0.0, -1.0])
    def test_a_well_with_no_glucose_in_it_is_refused(self, bad):
        with pytest.raises(ValueError, match="glucose_g_per_L"):
            substrate_limited_capacity(0.07, glucose_g_per_L=bad, growth_rate=0.30)

    def test_a_well_with_no_cells_in_it_is_refused(self):
        with pytest.raises(ValueError, match="initial_biomass"):
            substrate_limited_capacity(0.0, glucose_g_per_L=20.0, growth_rate=0.30)


BASE = CultureParameters(
    mu_max=0.40, growth_ic50=0.8, growth_hill=1.5,
    promoter_basal=1.0e-3, promoter_peak=2.5e-3, promoter_ec50=0.7, promoter_hill=2.0,
    carrying_capacity=4.0,
)


class TestReporterMaturationIsNowReachable:
    """`reporter.py` has modelled chromophore maturation all along and `culture.py` never used it.

    Every simulated well in this repository has therefore assumed a reporter that is fluorescent
    the instant it is translated. That is an assumption, not a measurement -- Nagai 2002
    (PMID 11753368) identifies chromophore oxidation as the rate-limiting step of YFP maturation
    and shows a single substitution changing it, so the rate belongs to the exact variant. It
    stays `None` by default so no unmeasured number enters, and it is now reachable so the cost
    of the assumption can be measured.
    """

    def test_the_default_is_still_instantaneous_maturation(self):
        assert BASE.k_mat is None

    def test_a_maturing_reporter_lags_an_instantaneous_one(self):
        """The bias runs one way, and it is the way that matters on a short read.

        An unmatured reporter is invisible, so a 4.14 h read of a slow-maturing construct
        understates the promoter. Reading the same culture for longer closes the gap, which is
        why the effect is a property of the read window and not a calibration offset.
        """
        times = np.linspace(0.0, 4.14, 25)
        instant = simulate_culture(times, 1.0, BASE, 0.07)
        maturing = simulate_culture(
            times, 1.0, CultureParameters(**{**BASE.__dict__, "k_mat": 1.0}), 0.07)

        assert maturing["reporter"][-1] < instant["reporter"][-1]

    def test_a_faster_maturing_reporter_lags_less(self):
        times = np.linspace(0.0, 4.14, 25)
        slow = simulate_culture(
            times, 1.0, CultureParameters(**{**BASE.__dict__, "k_mat": 0.5}), 0.07)
        fast = simulate_culture(
            times, 1.0, CultureParameters(**{**BASE.__dict__, "k_mat": 4.0}), 0.07)

        assert slow["reporter"][-1] < fast["reporter"][-1]

    @pytest.mark.parametrize("bad", [0.0, -1.0])
    def test_a_non_positive_maturation_rate_is_refused_at_construction(self, bad):
        """Caught where the parameters are built, not deep inside an integrator.

        Zero is the dangerous value: it reads as "no maturation" and means "never matures", so
        it would silently produce a permanently dark reporter.
        """
        with pytest.raises(ValueError, match="k_mat"):
            CultureParameters(**{**BASE.__dict__, "k_mat": bad})
