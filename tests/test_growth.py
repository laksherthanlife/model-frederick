import numpy as np
import pytest

from ystwin.growth import specific_growth_rate
from ystwin.readings import CorrectedOD


def test_recovers_the_true_rate_for_pure_exponential_growth():
    t = np.linspace(0, 10, 120)
    od = 0.05 * np.exp(0.35 * t)

    mu = specific_growth_rate(t, CorrectedOD(od))

    assert mu == pytest.approx(0.35, abs=0.005)


def test_rate_declines_monotonically_as_a_logistic_culture_saturates():
    t = np.linspace(0, 24, 200)
    od = 2.0 / (1 + (2.0 / 0.05 - 1) * np.exp(-0.4 * t))

    mu = specific_growth_rate(t, CorrectedOD(od))

    interior = mu[5:-5]
    assert np.all(np.diff(interior) < 1e-9)
    assert interior[0] > 0.3
    assert interior[-1] < 0.02


def test_stays_close_to_truth_under_realistic_reader_noise():
    rng = np.random.default_rng(0)
    t = np.linspace(0, 10, 120)
    clean = 0.05 * np.exp(0.3 * t)
    od = clean + rng.normal(0, 0.005, clean.shape)

    mu = specific_growth_rate(t, CorrectedOD(od), window_h=2.0)

    assert np.median(mu) == pytest.approx(0.3, abs=0.02)


def test_nonpositive_optical_density_is_rejected_rather_than_silently_logged():
    t = np.linspace(0, 5, 30)
    od = np.linspace(-0.01, 0.5, 30)

    with pytest.raises(ValueError, match="positive"):
        specific_growth_rate(t, CorrectedOD(od))


def test_returns_one_rate_per_timepoint():
    t = np.linspace(0, 8, 64)
    od = 0.05 * np.exp(0.25 * t)

    assert specific_growth_rate(t, CorrectedOD(od)).shape == t.shape


def test_uneven_sampling_is_handled_by_resampling_onto_a_regular_grid():
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 10, 90))
    od = 0.05 * np.exp(0.3 * t)

    mu = specific_growth_rate(t, CorrectedOD(od))

    assert np.median(mu) == pytest.approx(0.3, abs=0.02)


class TestRobustMaximumGrowthRate:
    """The maximum of a noisy derivative is not a growth rate.

    These cultures are inoculated at OD ~0.11 against a blank of ~0.11, so the
    blank-corrected density starts near zero and the early log-derivative is
    dominated by reader noise. Taking a plain max there returns the largest noise
    spike, which is how a stressed well ends up "growing" faster than an untreated
    one.
    """

    def test_it_recovers_the_true_rate_for_clean_exponential_growth(self):
        from ystwin.growth import max_specific_growth_rate

        t = np.linspace(0, 10, 120)
        od = 0.05 * np.exp(0.35 * t)

        assert max_specific_growth_rate(t, CorrectedOD(od)) == pytest.approx(0.35, abs=0.02)

    def test_a_noisy_near_zero_start_does_not_inflate_the_estimate(self):
        from ystwin.growth import max_specific_growth_rate

        rng = np.random.default_rng(0)
        t = np.linspace(0, 12, 145)
        od = 0.004 * np.exp(0.25 * t) + rng.normal(0, 0.004, t.size)
        od = np.clip(od, 1e-4, None)

        plain_max = float(np.nanmax(specific_growth_rate(t, CorrectedOD(od))))
        robust = max_specific_growth_rate(t, CorrectedOD(od), min_od=0.02)

        assert plain_max > 0.6, "the naive statistic should be visibly inflated here"
        assert robust < 0.45

    def test_a_slower_culture_reports_a_lower_rate_than_a_faster_one(self):
        from ystwin.growth import max_specific_growth_rate

        t = np.linspace(0, 12, 145)
        fast = max_specific_growth_rate(t, CorrectedOD(0.05 * np.exp(0.40 * t)))
        slow = max_specific_growth_rate(t, CorrectedOD(0.05 * np.exp(0.18 * t)))

        assert fast > slow

    def test_a_culture_that_never_clears_the_detection_floor_returns_nan(self):
        from ystwin.growth import max_specific_growth_rate

        t = np.linspace(0, 12, 145)
        od = np.full_like(t, 0.005)

        assert np.isnan(max_specific_growth_rate(t, CorrectedOD(od), min_od=0.02))

    def test_the_window_excludes_points_below_the_detection_floor(self):
        from ystwin.growth import growth_window

        t = np.linspace(0, 12, 100)
        od = 0.002 * np.exp(0.4 * t)

        mask = growth_window(t, CorrectedOD(od), min_od=0.02)

        assert not mask[0]
        assert mask[-1]
        assert od[mask].min() >= 0.02

    def test_using_a_quantile_rather_than_the_maximum_resists_a_single_spike(self):
        from ystwin.growth import max_specific_growth_rate

        t = np.linspace(0, 12, 145)
        od = 0.05 * np.exp(0.30 * t)
        spiked = od.copy()
        spiked[70] *= 1.25

        clean = max_specific_growth_rate(t, CorrectedOD(od))
        with_spike = max_specific_growth_rate(t, CorrectedOD(spiked))

        assert abs(with_spike - clean) < 0.08
