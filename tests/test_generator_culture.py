"""S0: the mechanistic culture model the synthetic data comes from.

Two coupled effects, both established and both visible in the real plates:

  growth      mu(d)      = mu0 / (1 + (d/IC50)^m)          stressor slows growth
  promoter    k_synth(d) = basal + (peak-basal) * hill(d)  stressor induces the reporter
  reporter    dR/dt      = k_synth - (mu + k_deg) R        dilution sets the gain

The third equation is why synthetic data generated from a *constant* promoter still
shows an apparent dose response: growth inhibition alone produces one. A generator
that skipped it would produce data on which the dilution correction had nothing to
find, and would validate the analysis against a strawman.
"""

import numpy as np
import pytest

from ystwin.generator.culture import CultureParameters, simulate_culture

BASE = CultureParameters(
    mu_max=0.62, growth_ic50=0.8, growth_hill=1.5,
    promoter_basal=1.0e-3, promoter_peak=2.4e-3, promoter_ec50=0.9, promoter_hill=2.0,
    carrying_capacity=1.4, k_deg=0.0,
)


def test_an_unstressed_culture_grows_at_the_stated_maximum_rate():
    t = np.linspace(0, 4.14, 25)

    out = simulate_culture(t, dose=0.0, params=BASE, initial_biomass=0.05)

    mu = np.gradient(np.log(out["biomass"]), t)
    assert mu[2] == pytest.approx(0.62, rel=0.1)


def test_the_stressor_slows_growth_monotonically():
    t = np.linspace(0, 4.14, 25)
    rates = []
    for dose in (0.0, 0.2, 0.5, 1.0, 2.0):
        out = simulate_culture(t, dose=dose, params=BASE, initial_biomass=0.05)
        rates.append(out["biomass"][-1])

    assert all(a > b for a, b in zip(rates, rates[1:]))


def test_growth_is_halved_at_the_stated_ic50():
    t = np.linspace(0, 1.0, 40)
    free = simulate_culture(t, 0.0, BASE, 0.01)
    at_ic50 = simulate_culture(t, BASE.growth_ic50, BASE, 0.01)

    mu_free = np.gradient(np.log(free["biomass"]), t)[2]
    mu_half = np.gradient(np.log(at_ic50["biomass"]), t)[2]
    assert mu_half == pytest.approx(mu_free / 2, rel=0.12)


def test_the_promoter_is_induced_by_the_stressor():
    t = np.linspace(0, 4.14, 25)

    low = simulate_culture(t, 0.1, BASE, 0.05)
    high = simulate_culture(t, 1.0, BASE, 0.05)

    assert high["promoter_activity"][-1] > low["promoter_activity"][-1] * 1.5


def test_growth_saturates_at_the_carrying_capacity():
    t = np.linspace(0, 40, 200)

    out = simulate_culture(t, 0.0, BASE, 0.05)

    assert out["biomass"][-1] == pytest.approx(BASE.carrying_capacity, rel=0.05)


class TestTheDilutionEffectIsInTheGeneratedData:
    def test_a_constant_promoter_still_produces_an_apparent_dose_response(self):
        """The confound the analysis is meant to find must exist in the synthetic data."""
        flat = CultureParameters(
            mu_max=0.62, growth_ic50=0.8, growth_hill=1.5,
            promoter_basal=1.0e-3, promoter_peak=1.0e-3,   # no induction at all
            promoter_ec50=0.9, promoter_hill=2.0,
            carrying_capacity=1.4, k_deg=0.0,
        )
        t = np.linspace(0, 4.14, 25)

        control = simulate_culture(t, 0.0, flat, 0.05)
        stressed = simulate_culture(t, 2.0, flat, 0.05)

        specific_control = control["reporter"][-1]
        specific_stressed = stressed["reporter"][-1]
        assert specific_stressed > specific_control * 1.2

    def test_and_the_promoter_activity_it_reports_is_genuinely_flat(self):
        flat = CultureParameters(
            mu_max=0.62, growth_ic50=0.8, growth_hill=1.5,
            promoter_basal=1.0e-3, promoter_peak=1.0e-3,
            promoter_ec50=0.9, promoter_hill=2.0,
            carrying_capacity=1.4, k_deg=0.0,
        )
        t = np.linspace(0, 4.14, 25)

        control = simulate_culture(t, 0.0, flat, 0.05)
        stressed = simulate_culture(t, 2.0, flat, 0.05)

        assert control["promoter_activity"][-1] == pytest.approx(
            stressed["promoter_activity"][-1], rel=1e-9
        )


def test_a_degron_tagged_reporter_reaches_a_lower_steady_state():
    t = np.linspace(0, 20, 200)
    stable = simulate_culture(t, 1.0, BASE, 0.05)
    tagged = simulate_culture(
        t, 1.0, CultureParameters(**{**BASE.__dict__, "k_deg": 1.0}), 0.05
    )

    assert tagged["reporter"][-1] < stable["reporter"][-1]


def test_negative_parameters_are_refused():
    with pytest.raises(ValueError, match="mu_max"):
        CultureParameters(
            mu_max=-0.1, growth_ic50=0.8, growth_hill=1.5,
            promoter_basal=1e-3, promoter_peak=2e-3, promoter_ec50=0.9,
            promoter_hill=2.0, carrying_capacity=1.4,
        )


def test_every_returned_series_is_on_the_requested_grid():
    t = np.linspace(0, 4.14, 25)

    out = simulate_culture(t, 0.5, BASE, 0.05)

    for key in ("biomass", "reporter", "promoter_activity", "growth_rate"):
        assert out[key].shape == t.shape


class TestHillFunctionsAreNumericallySafe:
    """Both curves are bounded, so no parameter value should overflow.

    Written as `d**n / (EC50**n + d**n)` the intermediate blows up long before the
    result would: at d=5, n=1000 the numerator overflows a float even though the
    answer is exactly 1. Augmentation draws wide, so it hits this.
    """

    def _params(self, **over):
        base = dict(
            mu_max=0.4, growth_ic50=0.8, growth_hill=1.5,
            promoter_basal=1e-3, promoter_peak=2e-3, promoter_ec50=0.9,
            promoter_hill=2.0, carrying_capacity=3.0,
        )
        base.update(over)
        return CultureParameters(**base)

    @pytest.mark.parametrize("hill", [0.5, 2.0, 50.0, 500.0, 5000.0])
    def test_induction_stays_finite_at_any_exponent(self, hill):
        params = self._params(promoter_hill=hill)

        value = params.promoter_activity_at(5.0)

        assert np.isfinite(value)
        assert params.promoter_basal <= value <= params.promoter_peak

    @pytest.mark.parametrize("hill", [0.5, 2.0, 50.0, 500.0, 5000.0])
    def test_inhibition_stays_finite_at_any_exponent(self, hill):
        params = self._params(growth_hill=hill)

        value = params.growth_rate_at(5.0)

        assert np.isfinite(value)
        assert 0.0 <= value <= params.mu_max

    def test_a_steep_curve_becomes_a_step_rather_than_a_nan(self):
        params = self._params(promoter_hill=5000.0, promoter_ec50=1.0)

        assert params.promoter_activity_at(0.5) == pytest.approx(params.promoter_basal)
        assert params.promoter_activity_at(2.0) == pytest.approx(params.promoter_peak)

    def test_zero_dose_gives_the_basal_and_unrestricted_values(self):
        params = self._params()

        assert params.promoter_activity_at(0.0) == pytest.approx(params.promoter_basal)
        assert params.growth_rate_at(0.0) == pytest.approx(params.mu_max)

    def test_the_stable_form_agrees_with_the_naive_one_where_that_is_computable(self):
        params = self._params(promoter_hill=2.0, promoter_ec50=0.9)
        d = 1.7
        naive_sat = d**2.0 / (0.9**2.0 + d**2.0)
        expected = params.promoter_basal + (params.promoter_peak - params.promoter_basal) * naive_sat

        assert params.promoter_activity_at(d) == pytest.approx(expected, rel=1e-12)
