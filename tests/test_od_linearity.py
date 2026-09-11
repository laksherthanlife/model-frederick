"""Finding where the reader stops responding linearly, from a dilution series.

`OpticalQualityGate.od_linear_max` is a placeholder, and the code says so: it is the single
most important number to measure for the platform. Everything downstream divides by OD, so
a growth rate taken above the linear range is wrong by however much the reader is
compressing, and no amount of modelling recovers it.

A dilution series settles it. A known dilution factor gives the OD that should have been
read; the reader gives what was read; the two agree until absorbance starts saturating.
What is returned is the highest OD still within tolerance -- so the answer is a measurement
with a stated tolerance rather than a number carried over from a cuvette.
"""

import numpy as np
import pytest

from ystwin.calib.od_linearity import LinearityFit, fit_linear_range
from ystwin.readings import CorrectedOD


def series(true_od, saturating_at=None):
    """Reader response: linear, then compressing above a threshold."""
    od = np.asarray(true_od, dtype=float)
    if saturating_at is None:
        return od
    return np.where(od <= saturating_at, od, saturating_at + (od - saturating_at) * 0.35)


class TestFitting:
    def test_a_perfectly_linear_reader_is_linear_to_the_top(self):
        true = np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6])
        fit = fit_linear_range(true, CorrectedOD(series(true)))

        assert fit.od_linear_max >= true.max()

    def test_it_finds_a_known_break_point(self):
        true = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6])
        fit = fit_linear_range(true, CorrectedOD(series(true, saturating_at=0.6)))

        assert fit.od_linear_max == pytest.approx(0.6, abs=0.25)

    def test_a_tighter_tolerance_never_reports_a_wider_range(self):
        true = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6])
        measured = series(true, saturating_at=0.6)

        tight = fit_linear_range(true, CorrectedOD(measured), tolerance=0.02).od_linear_max
        loose = fit_linear_range(true, CorrectedOD(measured), tolerance=0.15).od_linear_max

        assert tight <= loose

    def test_it_reports_the_slope_it_calibrated_against(self):
        true = np.array([0.05, 0.1, 0.2, 0.4, 0.8])
        fit = fit_linear_range(true, CorrectedOD(0.5 * series(true)))

        assert fit.slope == pytest.approx(0.5, rel=0.1)

    def test_it_survives_noise(self):
        rng = np.random.default_rng(0)
        true = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6])
        measured = series(true, 0.6) * rng.normal(1.0, 0.01, size=len(true))
        fit = fit_linear_range(true, CorrectedOD(measured))

        assert fit.od_linear_max == pytest.approx(0.6, abs=0.35)

    def test_it_returns_a_fit_that_names_its_tolerance(self):
        true = np.array([0.05, 0.1, 0.2, 0.4, 0.8])
        fit = fit_linear_range(true, CorrectedOD(series(true)), tolerance=0.07)

        assert isinstance(fit, LinearityFit)
        assert fit.tolerance == 0.07


class TestItRefusesWhatItCannotAnswer:
    def test_a_series_too_short_to_fit_is_refused(self):
        with pytest.raises(ValueError, match="at least"):
            fit_linear_range([0.1, 0.2], CorrectedOD([0.1, 0.2]))

    def test_mismatched_lengths_are_refused(self):
        with pytest.raises(ValueError, match="same length"):
            fit_linear_range([0.1, 0.2, 0.4, 0.8], CorrectedOD([0.1, 0.2, 0.4]))

    def test_a_series_that_never_leaves_the_linear_range_says_so(self):
        true = np.array([0.02, 0.04, 0.08, 0.12])
        fit = fit_linear_range(true, CorrectedOD(series(true)))

        assert fit.saturation_observed is False

    def test_a_series_that_does_leave_it_says_so(self):
        true = np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2])
        fit = fit_linear_range(true, CorrectedOD(series(true, saturating_at=0.4)))

        assert fit.saturation_observed is True

    def test_it_will_not_extrapolate_past_the_highest_dilution_run(self):
        true = np.array([0.02, 0.04, 0.08, 0.12])
        fit = fit_linear_range(true, CorrectedOD(series(true)))

        assert fit.od_linear_max <= true.max()
