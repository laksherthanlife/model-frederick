"""Lethality is a growth measurement, and reading it off the reporter was the error.

The first calibration fitted a toxic collapse to the reporter trace and got a lethal dose
of 7.2 mM for DTT. Growth on the same wells says otherwise: the culture is at 24% of the
control's rate by 2 mM. The reporter is an indirect, model-dependent readout of a quantity
optical density measures directly, and where a direct measurement exists it wins.

The same measurement settles a second thing the first pass got wrong. Both ladders were
fitted across their full range, but induction only remains interpretable while the culture
is healthy -- the "peak" that anchored the DTT half-maximal estimate sat at a dose where
growth was already down to 71% of control, and the peroxide peak at 23%. Inside the healthy
range neither curve saturates, so neither EC50 is identifiable, and saying so is the result.
"""

import numpy as np
import pytest

from ystwin.generator.panel_calibration import (
    healthy_doses,
    induction_is_identifiable,
    lethal_dose_from_growth,
)


class TestLethalityFromGrowth:
    def test_it_finds_where_growth_halves(self):
        doses = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
        growth = np.array([0.40, 0.40, 0.30, 0.20, 0.05])

        assert lethal_dose_from_growth(doses, growth) == pytest.approx(2.0, abs=0.6)

    def test_a_ladder_that_never_halves_growth_returns_nothing(self):
        doses = np.array([0.0, 0.5, 1.0, 2.0])
        growth = np.array([0.40, 0.39, 0.38, 0.37])

        assert lethal_dose_from_growth(doses, growth) is None

    def test_it_interpolates_between_bracketing_doses(self):
        doses = np.array([0.0, 1.0, 2.0])
        growth = np.array([0.40, 0.40, 0.20])

        assert 1.0 < lethal_dose_from_growth(doses, growth) <= 2.0

    def test_it_uses_the_undosed_well_as_the_reference(self):
        doses = np.array([0.0, 1.0, 2.0])
        slow = lethal_dose_from_growth(doses, np.array([0.10, 0.10, 0.05]))
        fast = lethal_dose_from_growth(doses, np.array([0.40, 0.40, 0.20]))

        assert slow == pytest.approx(fast)

    def test_it_refuses_a_ladder_with_no_control(self):
        with pytest.raises(ValueError, match="undosed"):
            lethal_dose_from_growth([0.5, 1.0, 2.0], [0.4, 0.3, 0.2])


class TestTheHealthyRange:
    def test_it_keeps_wells_growing_near_the_control(self):
        keep = healthy_doses([0.0, 0.5, 1.0], [0.40, 0.38, 0.36])

        assert list(keep) == [True, True, True]

    def test_it_drops_a_badly_slowed_culture(self):
        keep = healthy_doses([0.0, 1.0, 2.0], [0.40, 0.28, 0.10])

        assert list(keep) == [True, True, False]

    def test_the_threshold_is_relative_not_absolute(self):
        """A slow strain is not a stressed one; what matters is the fall from its own control."""
        keep = healthy_doses([0.0, 1.0], [0.08, 0.075])

        assert list(keep) == [True, True]


class TestIdentifiabilityInsideTheHealthyRange:
    def test_a_curve_still_rising_at_the_healthy_edge_is_not_identifiable(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0])
        activity = np.array([1000.0, 1100.0, 1200.0, 1350.0, 1500.0])

        assert not induction_is_identifiable(doses, activity)

    def test_a_curve_that_flattens_before_the_edge_is_identifiable(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0])
        activity = np.array([1000.0, 1300.0, 1450.0, 1495.0, 1500.0])

        assert induction_is_identifiable(doses, activity)

    def test_a_flat_ladder_is_not_identifiable(self):
        doses = np.array([0.0, 0.5, 1.0])
        activity = np.array([1000.0, 1002.0, 999.0])

        assert not induction_is_identifiable(doses, activity)
