"""The reporter is a dynamic state, not an instantaneous readout of stress.

For a stable fluorescent protein, per-cell signal relaxes to k_synth/(mu + k_deg).
Growth rate therefore sets the gain. Any 'stress state' inferred from raw RFU/OD
without inverting that gain is partly, and sometimes mostly, a measurement of 1/mu.
"""

import numpy as np
import pytest

from ystwin.reporter import (
    ReporterKinetics,
    naive_specific_fluorescence,
    promoter_activity,
    promoter_activity_from_total,
    simulate_reporter,
)
from ystwin.readings import CorrectedOD, CorrectedRFU, SpecificFluorescence


def _sigmoid_mu(t, hi=0.40, lo=0.05, t_mid=15.0, width=2.0):
    return lo + (hi - lo) / (1 + np.exp((t - t_mid) / width))


def test_steady_state_specific_signal_equals_synthesis_over_dilution():
    t = np.linspace(0, 60, 1200)
    mu = np.full_like(t, 0.30)
    k_syn = np.full_like(t, 1.2)

    r = simulate_reporter(t, k_syn, mu, ReporterKinetics(k_deg=0.0), r0=1.2 / 0.30)

    assert r[-1] == pytest.approx(1.2 / 0.30, rel=1e-3)


def test_specific_fluorescence_climbs_when_growth_slows_at_constant_promoter_activity():
    """The confound, stated as a test: this is a fake 'stress induction'."""
    t = np.linspace(0, 40, 800)
    mu = _sigmoid_mu(t)
    k_syn = np.ones_like(t)

    r = simulate_reporter(t, k_syn, mu, ReporterKinetics(k_deg=0.0), r0=1.0 / 0.40)

    fold_change = r[-1] / r[0]
    assert fold_change > 3.0, "expected the dilution artefact to look like induction"


def test_promoter_activity_recovers_constant_input_across_the_growth_transition():
    """The fix: inverting the dilution gain removes the artefact."""
    t = np.linspace(0, 40, 800)
    mu = _sigmoid_mu(t)
    k_syn = np.ones_like(t)
    r = simulate_reporter(t, k_syn, mu, ReporterKinetics(k_deg=0.0), r0=1.0 / 0.40)

    recovered = promoter_activity(t, SpecificFluorescence(r), mu, ReporterKinetics(k_deg=0.0))

    interior = recovered[20:-20]
    assert np.max(np.abs(interior - 1.0)) < 0.05


def test_forward_then_inverse_round_trips_a_time_varying_promoter_input():
    t = np.linspace(0, 30, 900)
    mu = _sigmoid_mu(t, hi=0.35, lo=0.10, t_mid=12.0)
    k_syn = 1.0 + 0.8 * np.exp(-((t - 18.0) ** 2) / (2 * 2.5**2))
    kin = ReporterKinetics(k_deg=0.05)
    r = simulate_reporter(t, k_syn, mu, kin, r0=k_syn[0] / (mu[0] + kin.k_deg))

    recovered = promoter_activity(t, SpecificFluorescence(r), mu, kin)

    interior = slice(30, -30)
    assert np.allclose(recovered[interior], k_syn[interior], atol=0.03)


def test_maturation_delays_the_observed_signal_behind_the_promoter_pulse():
    t = np.linspace(0, 30, 900)
    mu = np.full_like(t, 0.25)
    k_syn = 1.0 + 2.0 * np.exp(-((t - 10.0) ** 2) / (2 * 1.0**2))
    slow = ReporterKinetics(k_deg=0.3, k_mat=1.0)
    fast = ReporterKinetics(k_deg=0.3, k_mat=60.0)

    r_slow = simulate_reporter(t, k_syn, mu, slow, r0=k_syn[0] / (mu[0] + slow.k_deg))
    r_fast = simulate_reporter(t, k_syn, mu, fast, r0=k_syn[0] / (mu[0] + fast.k_deg))

    assert t[np.argmax(r_slow)] > t[np.argmax(r_fast)] + 0.3


def test_inversion_accounting_for_maturation_recovers_a_promoter_pulse():
    t = np.linspace(0, 30, 1800)
    mu = _sigmoid_mu(t, hi=0.35, lo=0.08, t_mid=14.0)
    k_syn = 1.0 + 1.5 * np.exp(-((t - 18.0) ** 2) / (2 * 2.0**2))
    kin = ReporterKinetics(k_deg=0.2, k_mat=2.0)
    r = simulate_reporter(t, k_syn, mu, kin, r0=k_syn[0] / (mu[0] + kin.k_deg))

    recovered = promoter_activity(t, SpecificFluorescence(r), mu, kin)

    interior = slice(80, -80)
    assert np.allclose(recovered[interior], k_syn[interior], atol=0.08)


def test_naive_specific_fluorescence_is_just_the_ratio_it_claims_to_be():
    rfu = np.array([1000.0, 2000.0, 4000.0])
    od = np.array([0.1, 0.2, 0.5])

    specific = naive_specific_fluorescence(CorrectedRFU(rfu), CorrectedOD(od))

    # `.values` rather than comparing the reading itself: the return is a
    # SpecificFluorescence now, and unwrapping it deliberately takes a keystroke so
    # that every place the label is dropped is visible in the diff.
    assert specific.values == pytest.approx([10000, 10000, 8000])


def test_zero_maturation_rate_is_rejected_as_a_degenerate_reporter():
    with pytest.raises(ValueError, match="maturation"):
        ReporterKinetics(k_deg=0.1, k_mat=0.0)


class TestTotalSignalForm:
    """The mu-free route to the same quantity, and why it is more accurate.

    k_synth = (dF/dt)/X + k_deg (F/X) follows from R = F/X by substitution, and mu
    cancels exactly. That removes the growth-rate estimate -- itself a smoothed
    derivative of a noisy optical trace -- from the calculation entirely.
    """

    def _culture(self, activity=800.0, cv=0.0, seed=0, n=145):
        """A decelerating culture with a known constant promoter activity."""
        rng = np.random.default_rng(seed)
        t = np.linspace(0.0, 24.0, n)
        mu = 0.35 * np.exp(-0.06 * t)
        biomass = 0.05 * np.exp(np.cumsum(mu) * (t[1] - t[0]))
        per_cell = np.empty_like(t)
        per_cell[0] = activity / mu[0]
        for i in range(1, t.size):
            steady = activity / mu[i]
            per_cell[i] = steady + (per_cell[i - 1] - steady) * np.exp(-mu[i] * (t[i] - t[i - 1]))
        total = per_cell * biomass
        if cv:
            total = total * (1 + rng.normal(0, cv, t.size))
            biomass = biomass * (1 + rng.normal(0, cv, t.size))
        return t, total, biomass, mu, activity

    def test_recovers_a_known_activity_without_being_given_growth(self):
        t, total, biomass, _, truth = self._culture()
        recovered = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass))
        mid = slice(20, -20)
        assert np.allclose(recovered[mid], truth, rtol=0.05)

    def test_beats_the_per_cell_route_under_noise(self):
        """The point of the change: fewer error sources, not a different answer."""
        from ystwin.growth import specific_growth_rate

        t, total, biomass, _, truth = self._culture(cv=0.02, seed=0)
        estimated_mu = specific_growth_rate(t, CorrectedOD(biomass))
        per_cell_route = promoter_activity(t, SpecificFluorescence(total / biomass), estimated_mu)
        total_route = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass))

        mid = slice(20, -20)
        err = lambda a: float(np.sqrt(np.mean(((a[mid] - truth) / truth) ** 2)))  # noqa: E731
        assert err(total_route) < err(per_cell_route)

    def test_agrees_with_the_per_cell_route_when_growth_is_known_exactly(self):
        """Both invert the same equation, so with exact mu they must agree."""
        t, total, biomass, mu, _ = self._culture()
        per_cell_route = promoter_activity(t, SpecificFluorescence(total / biomass), mu)
        total_route = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass))
        mid = slice(20, -20)
        assert np.allclose(per_cell_route[mid], total_route[mid], rtol=0.02)

    def test_the_cancellation_holds_for_a_degron_reporter(self):
        t, _, biomass, mu, _ = self._culture()
        k_deg = 0.6
        kinetics = ReporterKinetics(k_deg=k_deg)
        truth = 800.0
        per_cell = simulate_reporter(t, np.full(t.size, truth), mu, kinetics)
        total = per_cell * biomass

        recovered = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass), kinetics)
        mid = slice(20, -20)
        assert np.allclose(recovered[mid], truth, rtol=0.10)

    def test_maturation_is_refused_rather_than_silently_wrong(self):
        t, total, biomass, _, _ = self._culture()
        with pytest.raises(ValueError, match="maturation is not recoverable"):
            promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass),
                                         ReporterKinetics(k_mat=2.0))

    def test_non_positive_biomass_is_refused(self):
        t, total, biomass, _, _ = self._culture()
        biomass[5] = 0.0
        with pytest.raises(ValueError, match="strictly positive"):
            promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass))


@pytest.mark.parametrize("origin", [0.0, 10000.0])
def test_total_signal_derivative_uses_each_irregular_timestamp(origin):
    elapsed = np.array([0.0, 0.07, 0.2, 0.51, 0.9, 1.4, 2.1, 3.0, 4.2, 6.0, 7.8])
    t = origin + elapsed
    elapsed = t - origin
    total = 2.0 + 3.0 * elapsed + 0.5 * elapsed**2
    biomass = 0.2 + 0.03 * elapsed
    kinetics = ReporterKinetics(k_deg=0.1)
    expected = (3.0 + elapsed + kinetics.k_deg * total) / biomass

    got = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass),
                                      kinetics, window_h=2.0)

    np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("kinetics", [ReporterKinetics(k_deg=0.1),
                                      ReporterKinetics(k_deg=0.1, k_mat=2.0,
                                                       k_deg_immature=0.2)])
def test_per_cell_and_maturation_derivatives_use_irregular_timestamps(kinetics):
    t = np.array([0.0, 0.07, 0.2, 0.51, 0.9, 1.4, 2.1, 3.0, 4.2, 6.0, 7.8])
    mature = 2.0 + 3.0 * t + 0.5 * t**2
    mu = np.full_like(t, 0.2)
    balance = 3.0 + t + (mu + kinetics.k_deg) * mature
    expected = balance
    if kinetics.has_maturation:
        immature = balance / kinetics.k_mat
        derivative = (1.0 + (mu + kinetics.k_deg) * (3.0 + t)) / kinetics.k_mat
        expected = derivative + (mu + kinetics.k_mat + kinetics.immature_loss) * immature

    got = promoter_activity(t, SpecificFluorescence(mature), mu, kinetics, window_h=2.0)

    np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-10)


def test_uniform_grid_keeps_the_existing_savgol_derivative():
    from scipy.signal import savgol_filter

    t = np.linspace(0.0, 4.0, 25)
    signal = 100.0 + np.sin(t) + np.random.default_rng(42).normal(0.0, 0.01, t.size)
    expected = savgol_filter(signal, window_length=7, polyorder=3, deriv=1,
                             delta=float(np.median(np.diff(t))))

    got = promoter_activity_from_total(t, CorrectedRFU(signal), CorrectedOD(np.ones_like(t)),
                                      window_h=1.0)

    np.testing.assert_array_equal(got, expected)
