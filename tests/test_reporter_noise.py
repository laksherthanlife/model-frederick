"""Can the dilution inversion survive real reader noise?

The inversion differentiates the signal once (no maturation) or twice (with
maturation). Differentiation amplifies noise, so an inversion that is exact on
clean simulated traces can still be useless on a plate reader. If this fails, the
promoter-activity estimate needs a state-space filter, not a smoothing derivative.
"""

import numpy as np
import pytest

from ystwin.reporter import ReporterKinetics, promoter_activity, simulate_reporter
from ystwin.readings import SpecificFluorescence

READER_CV = 0.01  # 1% coefficient of variation is typical for a Synergy H1


def _noisy_case(seed, kinetics, cv=READER_CV, n=145, hours=24.0):
    """Induction pulse on a saturating culture, observed through a noisy reader."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, hours, n)
    mu = 0.05 + 0.30 / (1 + np.exp((t - 10.0) / 2.0))
    k_syn = 1.0 + 2.0 / (1 + np.exp(-(t - 12.0) / 1.5))
    clean = simulate_reporter(t, k_syn, mu, kinetics, r0=k_syn[0] / (mu[0] + kinetics.k_deg))
    noisy = clean * (1 + rng.normal(0, cv, t.shape))
    return t, mu, k_syn, clean, noisy


def _relative_error(recovered, truth, interior=slice(20, -20)):
    return float(np.median(np.abs(recovered[interior] - truth[interior]) / truth[interior]))


def test_inversion_without_maturation_tolerates_one_percent_reader_noise():
    t, mu, k_syn, _, noisy = _noisy_case(0, ReporterKinetics(k_deg=0.0))

    recovered = promoter_activity(t, SpecificFluorescence(noisy), mu, ReporterKinetics(k_deg=0.0), window_h=4.0)

    assert _relative_error(recovered, k_syn) < 0.15


def test_inversion_still_tracks_the_shape_of_the_induction_under_noise():
    t, mu, k_syn, _, noisy = _noisy_case(1, ReporterKinetics(k_deg=0.0))

    recovered = promoter_activity(t, SpecificFluorescence(noisy), mu, ReporterKinetics(k_deg=0.0), window_h=4.0)
    interior = slice(20, -20)

    r = np.corrcoef(recovered[interior], k_syn[interior])[0, 1]
    assert r > 0.9


def test_the_smoothing_window_trades_noise_against_bias_as_expected():
    t, mu, k_syn, _, noisy = _noisy_case(2, ReporterKinetics(k_deg=0.0))
    kin = ReporterKinetics(k_deg=0.0)

    narrow = promoter_activity(t, SpecificFluorescence(noisy), mu, kin, window_h=0.6)
    wide = promoter_activity(t, SpecificFluorescence(noisy), mu, kin, window_h=4.0)

    assert np.std(np.diff(wide)) < np.std(np.diff(narrow))


def test_double_differentiation_for_maturation_is_far_noisier():
    """Documented cost: recovering through a maturation step amplifies noise."""
    kin = ReporterKinetics(k_deg=0.1, k_mat=2.0)
    t, mu, k_syn, _, noisy = _noisy_case(3, kin)

    with_mat = promoter_activity(t, SpecificFluorescence(noisy), mu, kin, window_h=4.0)
    without_mat = promoter_activity(
        t, SpecificFluorescence(noisy), mu, ReporterKinetics(k_deg=0.1), window_h=4.0
    )

    assert np.std(np.diff(with_mat)) > np.std(np.diff(without_mat))


def test_noise_free_recovery_is_much_better_than_noisy_recovery():
    """Sanity: the estimator is limited by the noise, not by a coding error."""
    kin = ReporterKinetics(k_deg=0.0)
    t, mu, k_syn, clean, noisy = _noisy_case(4, kin)

    exact = promoter_activity(t, SpecificFluorescence(clean), mu, kin, window_h=4.0)
    approx = promoter_activity(t, SpecificFluorescence(noisy), mu, kin, window_h=4.0)

    assert _relative_error(exact, k_syn) < _relative_error(approx, k_syn)


def test_the_default_smoothing_window_scales_with_the_run_not_a_fixed_constant():
    from ystwin.reporter import default_activity_window_h

    short = default_activity_window_h(np.linspace(0, 12, 73))
    long = default_activity_window_h(np.linspace(0, 48, 289))

    assert long > short
    assert short == pytest.approx(2.0, rel=0.01)
    assert long == pytest.approx(8.0, rel=0.01)


def test_the_default_window_never_falls_below_a_few_sampling_intervals():
    from ystwin.reporter import default_activity_window_h

    sparse = default_activity_window_h(np.linspace(0, 6, 7))  # hourly sampling

    assert sparse >= 4 * 1.0


def test_the_default_window_recovers_activity_to_a_few_percent_at_reader_noise():
    """The quality claim the default has to earn, measured not assumed."""
    kin = ReporterKinetics(k_deg=0.0)
    errors = []
    for seed in range(8):
        t, mu, k_syn, _, noisy = _noisy_case(seed, kin)
        recovered = promoter_activity(t, SpecificFluorescence(noisy), mu, kin)  # no explicit window
        errors.append(_relative_error(recovered, k_syn))

    assert np.median(errors) < 0.03


def test_an_explicit_window_still_overrides_the_default():
    kin = ReporterKinetics(k_deg=0.0)
    t, mu, _, _, noisy = _noisy_case(5, kin)

    tight = promoter_activity(t, SpecificFluorescence(noisy), mu, kin, window_h=0.6)
    default = promoter_activity(t, SpecificFluorescence(noisy), mu, kin)

    assert np.std(np.diff(tight)) > np.std(np.diff(default))
