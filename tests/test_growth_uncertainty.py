"""How accurately a growth rate can be estimated, derived rather than fitted.

The noise split was picked by grid search: two free parameters tuned until a simulated
plate reproduced one measured summary statistic. Two parameters against one number is a
ridge of equally good answers, and the point chosen off it carried no information. Trying
to identify them from the plates instead fails outright -- with two plates per condition
each CV is a one-degree-of-freedom estimate, and regressing CV squared on 1/mu squared
returns an R2 of 0.07.

The growth term does not need to be fitted at all. A specific growth rate is the slope of
log optical density against time, so its standard error is the standard error of that
slope: the residual spread about the fit, divided by the spread of the time points and the
root of their number. Every one of those is measurable on the trace itself.

That leaves one quantity to fit rather than two, and the remaining one is fixed by
subtraction from the measured total instead of by search.
"""

import numpy as np
import pytest

from ystwin.growth import growth_rate_uncertainty, specific_growth_rate
from ystwin.readings import CorrectedOD


def trace(times, mu, od0=0.05, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    clean = od0 * np.exp(mu * np.asarray(times, dtype=float))
    return clean * rng.lognormal(0.0, noise, size=len(times)) if noise else clean


class TestTheStandardErrorOfAGrowthRate:
    def test_a_noiseless_trace_has_no_uncertainty(self):
        times = np.linspace(0, 6, 13)

        assert growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.35))) == pytest.approx(0.0, abs=1e-9)

    def test_a_noisier_trace_is_less_certain(self):
        times = np.linspace(0, 6, 13)
        quiet = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.35, noise=0.01, seed=1)))
        loud = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.35, noise=0.10, seed=1)))

        assert loud > quiet

    def test_a_longer_window_is_more_certain(self):
        short = np.linspace(0, 2, 5)
        long = np.linspace(0, 8, 17)
        brief = growth_rate_uncertainty(short, CorrectedOD(trace(short, 0.35, noise=0.05, seed=2)))
        extended = growth_rate_uncertainty(long, CorrectedOD(trace(long, 0.35, noise=0.05, seed=2)))

        assert extended < brief

    def test_it_matches_the_spread_of_repeated_estimates(self):
        """The claim that makes it a standard error and not a formula: estimates drawn from
        independent noise must scatter by about what it predicts."""
        times = np.linspace(0, 6, 13)
        estimates = [float(np.max(specific_growth_rate(times, CorrectedOD(trace(times, 0.35, noise=0.05,
                                                                   seed=s)))))
                     for s in range(200)]
        predicted = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.35, noise=0.05, seed=0)))

        assert predicted == pytest.approx(float(np.std(estimates)), rel=1.5)

    def test_it_refuses_a_trace_too_short_to_fit(self):
        with pytest.raises(ValueError, match="at least"):
            growth_rate_uncertainty([0.0, 1.0], CorrectedOD([0.05, 0.07]))

    def test_it_refuses_mismatched_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            growth_rate_uncertainty([0.0, 1.0, 2.0], CorrectedOD([0.05, 0.07]))

    def test_a_dead_culture_gives_no_finite_estimate(self):
        times = np.linspace(0, 6, 13)

        assert not np.isfinite(growth_rate_uncertainty(times, CorrectedOD(np.zeros(13))))


class TestWhatItPropagatesInto:
    def test_the_relative_error_grows_as_growth_slows(self):
        """Why the term has dose structure: the same absolute error is a larger relative one
        in a culture that has been slowed."""
        times = np.linspace(0, 6, 13)
        sigma = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.35, noise=0.05, seed=3)))

        assert sigma / 0.10 > sigma / 0.35

    def test_it_is_an_absolute_error_not_a_relative_one(self):
        """A slow culture is not measured proportionally worse; it is measured equally badly
        in absolute terms, which is what makes the relative error blow up."""
        times = np.linspace(0, 6, 13)
        fast = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.40, noise=0.05, seed=4)))
        slow = growth_rate_uncertainty(times, CorrectedOD(trace(times, 0.10, noise=0.05, seed=4)))

        assert slow == pytest.approx(fast, rel=0.3)
