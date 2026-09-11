"""The maintenance measurement, and the twenty-fold gap it opens under the latent branch.

`bridge/latent_bridge.py` scales its ATP maintenance constraint by
``_MAX_STRESS_MAINTENANCE = 6.5`` mmol ATP/gDW/h, which its own comment calls "an envelope,
not their measured value". These tests pin the measurement that envelope stood in for, and
pin the disagreement, because the disagreement is the reason the module exists.

The sign question is deliberately NOT resolved here. The measurement compares a strain that
can mount the general stress response against one that cannot, and finds the responder
spends LESS -- which is a claim about what the GSR buys, not about what stress costs. Any
test asserting a direction would be asserting more than the paper does.
"""

from __future__ import annotations

import pytest

from ystwin.bridge.latent_bridge import _MAX_STRESS_MAINTENANCE, _RESTING_MAINTENANCE
from ystwin.bridge.maintenance_calibration import (
    ATP_PER_GLUCOSE_AEROBIC,
    MEASURED_MAINTENANCE,
    stress_maintenance_increment,
)


class TestTheMeasurementIsCarriedAsPublished:
    def test_both_states_are_present(self):
        assert set(MEASURED_MAINTENANCE) == {"gsr_intact", "gsr_absent"}

    def test_the_values_are_the_published_ones(self):
        assert MEASURED_MAINTENANCE["gsr_intact"].m_s_mmol_glucose_per_gdw_h == 0.066
        assert MEASURED_MAINTENANCE["gsr_absent"].m_s_mmol_glucose_per_gdw_h == 0.082

    def test_removing_the_stress_response_RAISES_maintenance(self):
        """The paper's actual finding, and the reason the sign is left open: the strain that
        CANNOT mount the response spends MORE, so the response lowers the cost of homeostasis."""
        assert (MEASURED_MAINTENANCE["gsr_absent"].m_s_mmol_glucose_per_gdw_h
                > MEASURED_MAINTENANCE["gsr_intact"].m_s_mmol_glucose_per_gdw_h)

    def test_the_two_values_overlap_within_one_standard_deviation(self):
        """Which is why the increment is returned as an interval and never as a point."""
        intact = MEASURED_MAINTENANCE["gsr_intact"]
        absent = MEASURED_MAINTENANCE["gsr_absent"]
        gap = absent.m_s_mmol_glucose_per_gdw_h - intact.m_s_mmol_glucose_per_gdw_h

        assert gap < intact.standard_deviation
        assert gap < absent.standard_deviation

    def test_every_measurement_names_its_source(self):
        for measurement in MEASURED_MAINTENANCE.values():
            assert "40181231" in measurement.source


class TestTheConversionToATP:
    def test_it_scales_linearly_with_the_assumed_yield(self):
        intact = MEASURED_MAINTENANCE["gsr_intact"]

        assert intact.as_atp(32.0) == pytest.approx(2 * intact.as_atp(16.0))

    def test_the_interval_spans_the_published_uncertainty_in_the_yield(self):
        low, high = MEASURED_MAINTENANCE["gsr_intact"].atp_interval()

        assert low == pytest.approx(0.066 * ATP_PER_GLUCOSE_AEROBIC[0])
        assert high == pytest.approx(0.066 * ATP_PER_GLUCOSE_AEROBIC[1])
        assert low < high

    def test_a_non_positive_yield_is_refused(self):
        with pytest.raises(ValueError, match="must be positive"):
            MEASURED_MAINTENANCE["gsr_intact"].as_atp(0.0)

    def test_total_maintenance_lands_around_one_mmol_atp(self):
        """The scale sanity check. Anything that returns tens here has a unit error."""
        for measurement in MEASURED_MAINTENANCE.values():
            low, high = measurement.atp_interval()

            assert 0.5 < low < 3.0
            assert 0.5 < high < 3.0


class TestTheGapWithTheAssertedEnvelope:
    """The finding. These will fail if the envelope is ever changed, which is the point --
    a change to `_MAX_STRESS_MAINTENANCE` should have to come past this file."""

    def test_the_asserted_increment_is_at_least_ten_times_the_measured_one(self):
        low, high = stress_maintenance_increment()

        assert _MAX_STRESS_MAINTENANCE / high > 10.0

    def test_the_asserted_full_stress_total_exceeds_any_measured_total(self):
        asserted_total = _RESTING_MAINTENANCE + _MAX_STRESS_MAINTENANCE
        highest_measured = max(m.atp_interval()[1] for m in MEASURED_MAINTENANCE.values())

        assert asserted_total > 4 * highest_measured

    def test_the_increment_is_an_interval_not_a_point(self):
        low, high = stress_maintenance_increment()

        assert low < high

    def test_pinning_the_yield_collapses_the_interval(self):
        """A caller who is willing to assert a P/O gets a point, and has asserted it
        explicitly rather than inherited it."""
        low, high = stress_maintenance_increment(atp_per_glucose=18.0)

        assert low == high == pytest.approx((0.082 - 0.066) * 18.0)


class TestTheCheckCanActuallyFail:
    def test_a_measurement_equal_to_the_envelope_would_close_the_gap(self):
        """Guards against the gap test passing for the wrong reason -- it must be measuring
        the ratio, not just that two numbers exist."""
        from ystwin.bridge.maintenance_calibration import MaintenanceMeasurement

        hypothetical = MaintenanceMeasurement(
            strain="x", genotype="x",
            m_s_mmol_glucose_per_gdw_h=_MAX_STRESS_MAINTENANCE / 18.0,
            standard_deviation=0.0, source="hypothetical")

        assert hypothetical.as_atp(18.0) == pytest.approx(_MAX_STRESS_MAINTENANCE)
