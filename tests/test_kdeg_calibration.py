"""Move 2: measure the reporter loss rate instead of assuming it is zero.

D2 showed the data demand a loss rate of 0.049-0.074 /h that the stable-YFP
assumption cannot supply. Loss has two candidate sources with different signatures:

  * proteolysis / maturation-independent degradation -- proportional to *time*
  * photobleaching -- proportional to the *number of reads*

They are separable only by reading two otherwise identical plates at different
intervals. Without that design, all a single plate yields is the total loss rate,
which is what the chase below measures.
"""

import numpy as np
import pytest

from ystwin.calib.kdeg import (
    fit_decay_rate,
    fit_kdeg_from_chase,
    minimum_consistent_kdeg,
    partition_loss,
)
from ystwin.readings import RawOD, RawRFU


def _chase(k_total=0.30, hours=6.0, n=37, r0=8000.0, noise=0.0, seed=0):
    """After translation shutoff the per-cell signal decays at the total loss rate."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, hours, n)
    signal = r0 * np.exp(-k_total * t)
    return t, signal * (1 + rng.normal(0, noise, t.shape))


def test_recovers_a_known_decay_rate_from_a_chase():
    t, signal = _chase(k_total=0.30)

    fit = fit_decay_rate(t, signal)

    assert fit.rate == pytest.approx(0.30, rel=0.01)


def test_the_standard_error_shrinks_as_the_chase_is_sampled_more_densely():
    sparse = fit_decay_rate(*_chase(n=8, noise=0.03, seed=1))
    dense = fit_decay_rate(*_chase(n=80, noise=0.03, seed=1))

    assert dense.stderr < sparse.stderr


def test_a_signal_that_does_not_decay_returns_a_rate_near_zero():
    t = np.linspace(0, 6, 37)

    fit = fit_decay_rate(t, np.full_like(t, 5000.0))

    assert fit.rate == pytest.approx(0.0, abs=1e-6)


def test_a_rising_signal_reports_a_negative_rate_rather_than_being_clamped():
    """Translation was not actually blocked. Say so instead of hiding it."""
    t = np.linspace(0, 6, 37)

    fit = fit_decay_rate(t, 5000.0 * np.exp(0.1 * t))

    assert fit.rate < 0


def test_residual_growth_during_the_chase_is_subtracted():
    """Cycloheximide arrests growth, but not always completely."""
    t, signal = _chase(k_total=0.30)
    od = 0.4 * np.exp(0.05 * t)

    fit = fit_kdeg_from_chase(t, rfu=RawRFU(signal * od), optical_density=RawOD(od), od_blank=0.0)

    assert fit.k_deg == pytest.approx(0.30 - 0.05, rel=0.03)


def test_a_photobleaching_control_is_subtracted_from_the_chase():
    t, signal = _chase(k_total=0.30)

    fit = fit_kdeg_from_chase(
        t, rfu=RawRFU(signal), optical_density=RawOD(np.ones_like(t)), od_blank=0.0,
        bleaching_rate=0.08,
    )

    assert fit.k_deg == pytest.approx(0.22, rel=0.03)


def test_the_chase_fit_carries_its_own_quality_measure():
    t, signal = _chase(noise=0.02, seed=4)

    fit = fit_kdeg_from_chase(t, rfu=RawRFU(signal), optical_density=RawOD(np.ones_like(t)), od_blank=0.0)

    assert fit.r_squared > 0.95


class TestSeparatingBleachingFromDegradation:
    def test_two_read_intervals_resolve_both_components(self):
        k_deg, bleach_per_read = 0.05, 0.004
        fast_dt, slow_dt = 1 / 6, 1 / 2  # every 10 min vs every 30 min
        loss_fast = k_deg + bleach_per_read / fast_dt
        loss_slow = k_deg + bleach_per_read / slow_dt

        split = partition_loss(
            [(fast_dt, loss_fast), (slow_dt, loss_slow)]
        )

        assert split.k_deg == pytest.approx(0.05, rel=1e-6)
        assert split.bleach_per_read == pytest.approx(0.004, rel=1e-6)

    def test_a_single_read_interval_cannot_resolve_them(self):
        with pytest.raises(ValueError, match="two distinct read intervals"):
            partition_loss([(1 / 6, 0.074)])

    def test_identical_intervals_are_refused_however_many_plates(self):
        with pytest.raises(ValueError, match="two distinct read intervals"):
            partition_loss([(1 / 6, 0.07), (1 / 6, 0.08), (1 / 6, 0.075)])

    def test_the_partition_predicts_a_read_interval_it_was_not_fitted_on(self):
        """What the split is FOR. Two cadences resolve the two components; the point of
        resolving them is to answer what a third cadence would cost, which is the question
        an experimenter actually has -- "if I read every five minutes, how much more signal
        do I lose?" Fitting on 10 and 30 min and predicting 5 is the only assertion here
        that is not a restatement of the inputs."""
        k_deg, bleach = 0.05, 0.004
        fitted = partition_loss([(dt, k_deg + bleach / dt) for dt in (1 / 6, 1 / 2)])

        unseen = 1 / 12  # every five minutes, twice as often as the fastest fitted plate
        assert fitted.total_at_interval(unseen) == pytest.approx(k_deg + bleach / unseen,
                                                                 rel=1e-6)

    def test_read_rarely_enough_and_only_degradation_is_left(self):
        """The limit that says the two components are separated and not just co-fitted:
        bleaching is per read, so amortised over a long enough interval it vanishes and the
        loss rate falls to k_deg. A partition that had absorbed bleaching into k_deg would
        not have this limit."""
        fitted = partition_loss([(1 / 6, 0.05 + 0.004 * 6), (1 / 2, 0.05 + 0.004 * 2)])

        assert fitted.total_at_interval(1e6) == pytest.approx(fitted.k_deg, rel=1e-6)
        assert fitted.total_at_interval(1 / 6) > fitted.total_at_interval(1 / 2)

    def test_pure_degradation_shows_no_dependence_on_read_frequency(self):
        split = partition_loss([(1 / 6, 0.05), (1 / 2, 0.05)])

        assert split.bleach_per_read == pytest.approx(0.0, abs=1e-9)
        assert split.k_deg == pytest.approx(0.05)


def test_the_minimum_consistent_rate_matches_the_bound_d2_reports():
    """Same quantity, reachable without running a chase at all."""
    from ystwin.diagnostics.dilution_confound import dilution_confound_report
    from ystwin.reporter import ReporterKinetics

    t = np.linspace(0, 24, 900)
    mu = np.full_like(t, 0.25)
    od = 0.05 * np.exp(0.25 * t)
    from ystwin.reporter import simulate_reporter

    k_syn = 1.0 / (1 + np.exp((t - 8.0) / 0.35))
    per_cell = simulate_reporter(t, k_syn, mu, ReporterKinetics(k_deg=0.30), r0=1 / 0.55)
    rfu = per_cell * od

    direct = minimum_consistent_kdeg(t, RawOD(od), RawRFU(rfu))
    via_report = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert direct == pytest.approx(via_report.implied_min_k_deg, rel=1e-6)
    assert direct == pytest.approx(0.30, abs=0.05)
