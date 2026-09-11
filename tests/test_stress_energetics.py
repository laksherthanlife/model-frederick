"""The joined chain: reporter -> stress -> maintenance -> growth, and the gates on it.

`bridge/stress_energetics.py` composes four modules that were built separately. These tests
pin the composition, and more importantly they pin the two things that make it refuse:

  * the maintenance scale it inherits from `latent_bridge` is ~20x the measured effect, and
  * the sign of the effect is not established by the measurement that fixes the scale.

A test that asserted a direction here would be asserting more than PMID 40181231 does, so
none of them do. What they assert is that the module does not either.
"""

from __future__ import annotations

import pytest

from ystwin.bridge.latent_bridge import _ACTIVITY_FULL_SCALE, _MAX_STRESS_MAINTENANCE
from ystwin.bridge.maintenance_calibration import MEASURED_MAINTENANCE
from ystwin.bridge.stress_energetics import (
    SIGN_IS_UNRESOLVED,
    maintenance_from_reporter,
)

FULL = _ACTIVITY_FULL_SCALE


class TestTheChainComposes:
    def test_a_dark_reporter_gives_resting_maintenance(self):
        result = maintenance_from_reporter(0.0, FULL)

        assert result.stress_fraction == 0.0
        assert (result.ngam_low, result.ngam_high) == (
            MEASURED_MAINTENANCE["gsr_intact"].atp_interval())

    def test_a_full_scale_reporter_adds_the_whole_measured_increment(self):
        resting_high = MEASURED_MAINTENANCE["gsr_intact"].atp_interval()[1]
        result = maintenance_from_reporter(FULL, FULL)

        assert result.stress_fraction == 1.0
        assert result.ngam_high > resting_high

    def test_maintenance_rises_monotonically_with_the_reading(self):
        readings = [0.0, 0.25 * FULL, 0.5 * FULL, FULL]
        highs = [maintenance_from_reporter(r, FULL).ngam_high for r in readings]

        assert highs == sorted(highs)

    def test_a_reading_above_full_scale_saturates_rather_than_extrapolating(self):
        """Full scale is the top of the construct's range, not a slope to continue past."""
        at_full = maintenance_from_reporter(FULL, FULL)
        over = maintenance_from_reporter(10 * FULL, FULL)

        assert over.stress_fraction == 1.0
        assert over.ngam_high == at_full.ngam_high

    def test_pinning_the_po_ratio_collapses_the_interval_to_a_point(self):
        result = maintenance_from_reporter(0.5 * FULL, FULL, atp_per_glucose=18.0)

        assert result.ngam_low == pytest.approx(result.ngam_high)


class TestTheTwentyFoldGapSurvivesComposition:
    """The finding, restated at the point of use. If these ever pass trivially the module
    has stopped carrying the disagreement it was written to carry."""

    def test_the_asserted_scale_gives_a_far_higher_ngam_at_full_stress(self):
        measured = maintenance_from_reporter(FULL, FULL)
        asserted = maintenance_from_reporter(FULL, FULL, use_measured_scale=False)

        assert asserted.ngam_high > 4 * measured.ngam_high

    def test_the_two_scales_disagree_about_the_direction_at_zero_stress(self):
        """A quieter surprise, and worth pinning: the MEASURED resting maintenance is HIGHER
        than the asserted resting constant, so the envelope is not uniformly conservative --
        it is low at rest and ~20x high under stress."""
        measured = maintenance_from_reporter(0.0, FULL)
        asserted = maintenance_from_reporter(0.0, FULL, use_measured_scale=False)

        assert measured.ngam_low > asserted.ngam_low

    def test_the_measured_increment_is_a_small_fraction_of_the_asserted_one(self):
        measured = maintenance_from_reporter(FULL, FULL)
        resting_high = MEASURED_MAINTENANCE["gsr_intact"].atp_interval()[1]
        increment = measured.ngam_high - resting_high

        assert increment < 0.1 * _MAX_STRESS_MAINTENANCE

    def test_each_result_names_which_scale_produced_it(self):
        assert "measured" in maintenance_from_reporter(0.0, FULL).scale_source
        assert "asserted" in maintenance_from_reporter(
            0.0, FULL, use_measured_scale=False).scale_source


class TestBothGatesHold:
    def test_a_result_is_not_usable_by_default(self):
        assert not maintenance_from_reporter(0.5 * FULL, FULL).usable_for_a_number

    def test_validating_the_reporter_alone_is_not_enough(self):
        result = maintenance_from_reporter(0.5 * FULL, FULL, reporter_validated=True)

        assert not result.usable_for_a_number

    def test_resolving_the_sign_alone_is_not_enough(self):
        result = maintenance_from_reporter(0.5 * FULL, FULL, sign_resolved=True)

        assert not result.usable_for_a_number

    def test_both_together_lift_it(self):
        result = maintenance_from_reporter(0.5 * FULL, FULL,
                                           reporter_validated=True, sign_resolved=True)

        assert result.usable_for_a_number

    def test_the_summary_says_which_gate_is_down(self):
        text = maintenance_from_reporter(0.5 * FULL, FULL, sign_resolved=True).summary()

        assert "REFUSED" in text
        assert "reporter" in text
        assert "sign unresolved" not in text

    def test_a_clean_summary_carries_no_refusal(self):
        text = maintenance_from_reporter(0.5 * FULL, FULL,
                                         reporter_validated=True,
                                         sign_resolved=True).summary()

        assert "REFUSED" not in text

    def test_the_sign_question_names_what_would_settle_it(self):
        """A refusal that does not say what would lift it is an obstacle, not a finding."""
        assert "same cultures" in SIGN_IS_UNRESOLVED


class TestItRefusesBadInputs:
    @pytest.mark.parametrize("scale", [0.0, -1.0])
    def test_a_non_positive_full_scale_is_refused(self, scale):
        with pytest.raises(ValueError, match="must be positive"):
            maintenance_from_reporter(100.0, scale)

    def test_the_full_scale_refusal_says_why_there_is_no_default(self):
        with pytest.raises(ValueError, match="construct and the instrument"):
            maintenance_from_reporter(100.0, 0.0)

    def test_a_negative_activity_is_refused(self):
        with pytest.raises(ValueError, match="non-negative"):
            maintenance_from_reporter(-1.0, FULL)

    def test_full_scale_is_required_rather_than_defaulted_to_the_local_constructs(self):
        """1305.0 belongs to four constructs on one instrument. A signature that defaulted
        to it would silently apply them to a different reporter."""
        import inspect

        signature = inspect.signature(maintenance_from_reporter)

        assert signature.parameters["activity_full_scale"].default is inspect.Parameter.empty


@pytest.mark.integration
class TestTheClaimsAboutTheGEMAreTrue:
    """The docstrings of both modules make quantitative claims about Yeast9. They are
    load-bearing -- the whole argument for preferring the measured scale is what the
    difference does to growth -- so they are checked against the model, not asserted.

    An earlier version of `maintenance_calibration`'s docstring said the asserted envelope
    made the constraint "saturate" and quoted a no-growth threshold of 19.19 as though it
    were a constant. Both were wrong: the response is exactly linear, and the threshold moves
    tenfold with the carbon supply. These tests exist because that went unnoticed.
    """

    @staticmethod
    def _growth(model, ngam):
        """None when the constraint is infeasible.

        The warning filter is deliberate and narrow: the suite runs with
        ``filterwarnings = ["error"]``, and probing past the no-growth threshold is the
        POINT of the bisection below, so cobra's 'Solver status is infeasible' warning is an
        expected result here rather than a defect. It is caught around the solve alone.
        """
        import warnings

        with model:
            model.reactions.get_by_id("r_4046").bounds = (ngam, 1000.0)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "Solver status is", UserWarning)
                solution = model.optimize()
            return (solution.objective_value
                    if solution.status == "optimal" and solution.objective_value is not None
                    else None)

    @pytest.fixture
    def aerobic(self, yeast_gem):
        def configure(glucose):
            yeast_gem.reactions.get_by_id("r_1714").lower_bound = glucose
            yeast_gem.reactions.get_by_id("r_1992").lower_bound = -1000.0
            return yeast_gem
        return configure

    @pytest.mark.parametrize("glucose", [-1.0, -1.5, -10.0])
    def test_growth_is_linear_in_the_ngam_bound(self, aerobic, glucose):
        """Not merely monotonic. Linearity is why one ratio describes the error at every
        carbon supply, and it is what makes 'saturates' the wrong word."""
        model = aerobic(glucose)
        base = self._growth(model, 0.0)
        slopes = [(base - self._growth(model, v)) / v for v in (1.0, 2.0, 5.0, 7.2, 10.0)]

        assert max(slopes) - min(slopes) < 1e-9 * max(slopes)

    @pytest.mark.parametrize("glucose", [-1.0, -1.5, -10.0])
    def test_the_slope_is_the_same_at_every_carbon_supply(self, aerobic, glucose):
        model = aerobic(glucose)
        base = self._growth(model, 0.0)

        assert (base - self._growth(model, 10.0)) / 10.0 == pytest.approx(0.004642, abs=1e-6)

    @pytest.mark.parametrize("glucose", [-1.0, -1.5, -10.0])
    def test_the_growth_penalty_ratio_is_4_39_and_does_not_move_with_carbon(self, aerobic,
                                                                           glucose):
        """The headline number, and the one that follows the TOTAL rather than the increment.
        The increment ratio is ~20x; what reaches growth is 4.39x."""
        model = aerobic(glucose)
        base = self._growth(model, 0.0)
        asserted = base - self._growth(model, 7.2)
        measured = base - self._growth(model, 1.64)

        assert asserted / measured == pytest.approx(4.39, abs=0.01)

    def test_the_no_growth_threshold_is_not_a_constant(self, aerobic):
        """19.19 was quoted as if it were one. It is a property of the carbon supply, and it
        moves by an order of magnitude across the range this project actually solves at."""
        def threshold(glucose):
            model = aerobic(glucose)
            low, high = 0.0, 400.0
            while high - low > 0.01:
                mid = (low + high) / 2
                g = self._growth(model, mid)
                low, high = (low, mid) if g is None or g < 1e-9 else (mid, high)
            return low

        assert threshold(-1.0) == pytest.approx(19.19, abs=0.05)
        assert threshold(-10.0) > 10 * threshold(-1.0)
