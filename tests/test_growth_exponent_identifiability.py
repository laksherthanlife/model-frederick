"""How well does the data pin the growth exponent in ``Vmax = capacity * mu``?

This file exists because the answer was asserted three times in this repository and the
first two were wrong in opposite directions.

**First claim: unidentifiable.** "Two dilution rates cannot establish a functional form."
Wrong, and wrong for an interesting reason -- see TestWhereTheIdentifiabilityComesFrom.

**Second claim: identifiable, 0.85 to 1.20.** Also wrong, and wrong in a way that is easier
to make and harder to notice. That range came from re-fitting the point estimate with one
state, then one strain, held out. **A leave-one-out spread measures how stable a fit is
under resampling. It is not a confidence interval and it is systematically far too narrow**
-- each refit still sees most of the data, so the estimates cluster whether or not the
parameter is well determined.

**What the data actually says.** Profile likelihood, F-test on the profiled sum of squares
with n=6 and k=3: the exponent is 1.04 with a 95% interval of **[0.53, 1.78]**. That is 3.6x
the resampling spread. So:

- the exponent IS bounded -- 0 and 3 are excluded, and the model is not indifferent to it;
- the exponent is NOT pinned to 1 -- mu^0.6 and mu^1.7 are both inside the interval;
- "every other exponent fits worse, and not narrowly" is true of the POINT ESTIMATES and
  false as a statement about what is ruled out.

Both facts are tested below, and the resampling spread is kept, clearly labelled as what it
is, because the contrast between the two numbers is the lesson.

`docs/research/KINETIC_FIT.md` said "two dilution rates cannot establish a functional form;
`vmax ~ mu` and any other monotone curve through the same pair are indistinguishable here."
That reasoning is correct for a fit whose only leverage is the growth-rate axis: two points,
two parameters, no residual freedom.

It is the wrong description of THIS fit. The calibration set is three strains at two
dilution rates, and the three strains carry three different pathway fluxes. The cyclase is
Michaelis-Menten, so the map from flux to product rate is *curved* -- and the same mu
observed at three flux levels on a curved response is not the same measurement three times.
The curvature is the leverage, and it is on the flux axis rather than the growth axis.

Fitted freely the exponent lands at 1.04 with three residual degrees of freedom, and it
survives dropping any single state (0.90-1.12) and any whole strain (1.01-1.09).

None of that rescues the MECHANISM, which is a separate question and still open: Elizondo's
own RT-qPCR has CrtYB mRNA falling by about a third from the low dilution rate to the high
one, which is the wrong direction for "capacity rises because there is more enzyme". What
these tests establish is that the exponent is measured rather than assumed. Why it is 1 is
in KINETIC_FIT.md section 6 and is not settled.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.kinetic.carotenoid import ELIZONDO2025

CALIBRATION = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"

# Where the point estimate lands under resampling. This is NOT the parameter's uncertainty
# and must not be read as one -- see TestTheIntervalIsWiderThanTheResamplingSpread below,
# which is the number to quote.
RESAMPLING_RANGE = (0.85, 1.20)

# The 95% profile-likelihood interval, from an F-test on the profiled sum of squares with
# n=6, k=3. This is the honest statement of what the data pins down.
PROFILE_INTERVAL = (0.53, 1.78)


@pytest.fixture(scope="module")
def measured():
    if not CALIBRATION.exists():
        raise AssertionError(f"{CALIBRATION} is tracked; a missing file is a bad checkout")
    return pd.read_csv(CALIBRATION, sep="\t")


def _predict(capacity, km, exponent, mu, flux):
    """Beta-carotene specific rate under ``Vmax = capacity * mu**exponent``.

    The same closed-form branch solve as ``solve_branch_from_flux``, restated here on
    arrays. It is restated rather than imported because the shipped function pins the
    exponent at 1 by construction, and a test of whether 1 is right cannot be written
    through a function that cannot express anything else.
    """
    vmax = capacity * mu ** exponent
    a = mu
    b = mu * km + vmax - flux
    c = -flux * km
    lycopene = (-b + np.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
    return vmax * lycopene / (km + lycopene)


def _fit(mu, flux, observed, exponent=None):
    """Least squares in log space. ``exponent=None`` fits it as a third free parameter."""
    from scipy.optimize import least_squares

    if exponent is None:
        result = least_squares(
            lambda p: np.log(np.maximum(_predict(p[0], p[1], p[2], mu, flux), 1e-12))
            - np.log(observed),
            x0=[2.2e-3, 5e-4, 1.0], bounds=([1e-8, 1e-8, -3.0], [1e2, 1e2, 5.0]))
        return result.x, float(np.sqrt(np.mean(result.fun ** 2)))
    result = least_squares(
        lambda p: np.log(np.maximum(_predict(p[0], p[1], exponent, mu, flux), 1e-12))
        - np.log(observed),
        x0=[2.2e-3, 5e-4], bounds=([1e-8, 1e-8], [1e2, 1e2]))
    return np.append(result.x, exponent), float(np.sqrt(np.mean(result.fun ** 2)))


def _axes(frame):
    return (frame.mu_per_h.to_numpy(),
            (frame.q_lycopene + frame.q_betacarotene).to_numpy(),
            frame.q_betacarotene.to_numpy())


class TestTheExponentIsMeasuredAndNotAssumed:
    def test_fitting_it_freely_lands_on_one(self, measured):
        mu, flux, observed = _axes(measured)

        (_, _, exponent), _ = _fit(mu, flux, observed)

        assert RESAMPLING_RANGE[0] < exponent < RESAMPLING_RANGE[1]

    def test_the_shipped_calibration_agrees_with_the_free_fit(self, measured):
        """`ELIZONDO2025` was fitted with the exponent pinned at 1. If the free fit had
        come back at 2, the shipped numbers would describe a different model."""
        mu, flux, observed = _axes(measured)

        (capacity, km, _), _ = _fit(mu, flux, observed, exponent=1.0)

        assert capacity == pytest.approx(ELIZONDO2025.capacity_mmol_per_gdcw, rel=0.15)
        assert km == pytest.approx(ELIZONDO2025.km_mmol_per_gdcw, rel=0.30)

    @pytest.mark.parametrize("exponent", [0.0, 0.5, 1.5, 2.0, 3.0])
    def test_every_other_exponent_fits_worse(self, measured, exponent):
        """The claim in one line. If these tied, the exponent would be unidentifiable and
        `Vmax ~ mu` would be an assumption; they do not tie, so it is a measurement."""
        mu, flux, observed = _axes(measured)

        _, rmse_one = _fit(mu, flux, observed, exponent=1.0)
        _, rmse_other = _fit(mu, flux, observed, exponent=exponent)

        assert rmse_other > rmse_one

    def test_the_margin_is_not_marginal(self, measured):
        """The nearest competitors are mu^0.5 and mu^1.5, and both are more than half as
        bad again. A 5% edge would not be worth calling identifiable."""
        mu, flux, observed = _axes(measured)

        _, one = _fit(mu, flux, observed, exponent=1.0)
        neighbours = [_fit(mu, flux, observed, exponent=e)[1] for e in (0.5, 1.5)]

        assert min(neighbours) > 1.5 * one


class TestNoSingleObservationCarriesIt:
    def test_dropping_any_one_state_leaves_the_exponent_near_one(self, measured):
        mu, flux, observed = _axes(measured)

        for held_out in range(len(mu)):
            keep = np.setdiff1d(np.arange(len(mu)), [held_out])
            (_, _, exponent), _ = _fit(mu[keep], flux[keep], observed[keep])

            assert RESAMPLING_RANGE[0] < exponent < RESAMPLING_RANGE[1], (
                f"dropping state {measured.condition[held_out]} moves the exponent to "
                f"{exponent:.3f}, so it rested on that one state")

    def test_dropping_a_whole_strain_leaves_it_near_one(self, measured):
        """The harder test. Leave-one-state-out still keeps both dilution rates for every
        strain; leave-one-strain-out removes a whole flux level, which is the axis the
        identifiability actually comes from."""
        mu, flux, observed = _axes(measured)

        for strain in measured.strain.unique():
            keep = (measured.strain != strain).to_numpy()
            (_, _, exponent), _ = _fit(mu[keep], flux[keep], observed[keep])

            assert RESAMPLING_RANGE[0] < exponent < RESAMPLING_RANGE[1], (
                f"without {strain} the exponent is {exponent:.3f}")


class TestWhereTheIdentifiabilityComesFrom:
    def test_it_is_the_flux_axis_and_not_the_growth_axis(self, measured):
        """The direct statement of why the earlier reasoning was wrong.

        Collapse the three strains to their mean flux and the fit has two points on the
        growth axis and nothing else -- which IS degenerate, exactly as
        `KINETIC_FIT.md` said. Keep the strains apart and it is not. The leverage was never
        on the growth axis."""
        mu, flux, observed = _axes(measured)

        collapsed = np.full_like(flux, flux.mean())
        _, degenerate_one = _fit(mu, collapsed, observed, exponent=1.0)
        _, degenerate_two = _fit(mu, collapsed, observed, exponent=2.0)

        assert degenerate_two == pytest.approx(degenerate_one, rel=0.05)

    def test_the_three_strains_really_do_sit_at_different_fluxes(self, measured):
        """The premise of the test above. If the strains had similar fluxes there would be
        no spread to get leverage from."""
        by_strain = (measured.q_lycopene + measured.q_betacarotene).groupby(
            measured.strain).mean()

        assert by_strain.max() / by_strain.min() > 3.0

    def test_and_the_response_to_flux_is_curved_rather_than_proportional(self, measured):
        """The other half of the premise: a straight line through the origin would rescale
        away and leave the degeneracy intact. Michaelis-Menten saturation is what does not
        rescale, so doubling the flux must give less than double the rate."""
        mu = np.array([0.18])
        capacity, km = ELIZONDO2025.capacity_mmol_per_gdcw, ELIZONDO2025.km_mmol_per_gdcw

        low = _predict(capacity, km, 1.0, mu, np.array([3.0e-4]))[0]
        high = _predict(capacity, km, 1.0, mu, np.array([6.0e-4]))[0]

        assert high < 2.0 * low


class TestTheRegisteredPredictionDoesNotTurnOnIt:
    def test_letting_the_exponent_float_barely_moves_it(self, measured):
        """The practical consequence. If the registered D = 0.18 /h prediction had swung on
        the exponent, it would be a prediction about a modelling choice rather than about
        the strains."""
        mu, flux, observed = _axes(measured)
        (cap_free, km_free, exponent), _ = _fit(mu, flux, observed)
        (cap_one, km_one, _), _ = _fit(mu, flux, observed, exponent=1.0)

        registered = np.array([0.18])
        for strain, group in measured.groupby("strain"):
            strain_flux = np.array([float((group.q_lycopene + group.q_betacarotene).mean())])
            pinned = _predict(cap_one, km_one, 1.0, registered, strain_flux)[0]
            floating = _predict(cap_free, km_free, exponent, registered, strain_flux)[0]

            assert floating == pytest.approx(pinned, rel=0.05), strain


class TestTheIntervalIsWiderThanTheResamplingSpread:
    """The correction, and the reason this class is last rather than first: everything above
    is about the point estimate, and a point estimate is not a claim until it has an
    interval on it."""

    @staticmethod
    def _profile(mu, flux, observed):
        """95% profile-likelihood interval on the exponent, by F-test on profiled SSE."""
        from scipy import stats

        def sse(exponent):
            _, rmse = _fit(mu, flux, observed, exponent=exponent)
            return rmse ** 2 * len(mu)

        grid = np.linspace(-1.0, 4.0, 1001)
        profiled = np.array([sse(e) for e in grid])
        best = profiled.min()
        n, k = len(mu), 3
        threshold = best * (1.0 + stats.f.ppf(0.95, 1, n - k) / (n - k))
        inside = grid[profiled <= threshold]
        return float(inside.min()), float(inside.max()), float(grid[profiled.argmin()])

    def test_the_ninety_five_percent_interval_is_the_one_to_quote(self, measured):
        mu, flux, observed = _axes(measured)

        low, high, point = self._profile(mu, flux, observed)

        assert point == pytest.approx(1.04, abs=0.05)
        assert low == pytest.approx(PROFILE_INTERVAL[0], abs=0.05)
        assert high == pytest.approx(PROFILE_INTERVAL[1], abs=0.05)

    def test_it_is_several_times_wider_than_the_resampling_spread(self, measured):
        """The whole point. If these were close, leave-one-out would be a usable proxy for
        uncertainty and this repository could go on quoting it."""
        mu, flux, observed = _axes(measured)

        low, high, _ = self._profile(mu, flux, observed)
        resampling_width = RESAMPLING_RANGE[1] - RESAMPLING_RANGE[0]

        assert (high - low) > 3.0 * resampling_width

    def test_one_half_is_inside_the_interval_even_though_it_fits_worse(self, measured):
        """The claim that had to be withdrawn. mu^0.5 fits worse than mu^1 by a wide margin
        on the point estimates AND sits inside the 95% interval. Both are true, and only the
        second bears on what has been ruled out."""
        mu, flux, observed = _axes(measured)

        low, high, _ = self._profile(mu, flux, observed)
        _, rmse_half = _fit(mu, flux, observed, exponent=0.5)
        _, rmse_one = _fit(mu, flux, observed, exponent=1.0)

        assert rmse_half > 1.5 * rmse_one
        assert low < 0.6

    def test_but_zero_and_three_are_still_excluded(self, measured):
        """Bounded is not nothing. A growth-rate-independent Vmax and a cubic one are both
        outside, so the data does constrain the law even if it does not pin it."""
        mu, flux, observed = _axes(measured)

        low, high, _ = self._profile(mu, flux, observed)

        assert low > 0.0
        assert high < 3.0
