"""A declared observation distribution, not a variance chosen after the outcome arrives."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from ystwin.estimator import Observation, ParticleFilter, TwinPriors
from ystwin.observation import ReporterOptics, observe_od, observe_rfu


PRIORS = TwinPriors(
    biomass=(1.0, 0.0), reporter=(2.0, 0.0), promoter_activity=(0.3, 0.0),
    growth_rate=(0.0, 0.0), od_rel_sigma=0.2, rfu_rel_sigma=0.3,
)
OPTICS = ReporterOptics(2.0, 0.5, 0.3, 0.4)


def _filter(**overrides):
    filt = ParticleFilter(replace(PRIORS, **overrides), OPTICS, 0.5, 0.2,
                          n_particles=3, seed=71, resample_threshold=0.0)
    filt._x = np.array([[0.5, 2.0, 0.3, 0.0], [1.5, 3.0, 0.5, 0.1], [4.0, 1.0, 0.7, -0.1]])
    filt._w = np.array([0.2, 0.5, 0.3])
    return filt


def _expected(filt, channel, carotenoid=0.7):
    if channel == "od":
        return np.array([observe_od(b, 0.0, filt.calibration, filt.gdcw_per_od) for b in filt._x[:, 0]])
    return np.array([observe_rfu(r, b, carotenoid, filt.optics) for b, r in filt._x[:, :2]])


@pytest.mark.parametrize("channel", ["od", "rfu"])
def test_changing_the_outcome_cannot_change_the_prior_predictive_variance(channel):
    filt = _filter()
    particles, weights = filt._x.copy(), filt._w.copy()
    rng = deepcopy(filt._rng.bit_generator.state)
    moments = [filt._predictive(filt._x, Observation(1.0, carotenoid=0.7, **{channel: y}))[0]
               for y in (-20.0, 0.0, 0.1, 3.0, 200.0)]
    for key in ("predicted", "state_variance", "measurement_variance"):
        assert len({moment[key] for moment in moments}) == 1, key
    np.testing.assert_array_equal(filt._x, particles)
    np.testing.assert_array_equal(filt._w, weights)
    assert filt._rng.bit_generator.state == rng


@pytest.mark.parametrize("channel", ["od", "rfu"])
def test_predictive_variance_is_the_weighted_gaussian_mixture_variance(channel):
    filt = _filter()
    expected = _expected(filt, channel)
    sigma = getattr(filt.priors, f"{channel}_rel_sigma") * np.maximum(np.abs(expected), 1e-6)
    mean = float(filt._w @ expected)
    state_variance = float(filt._w @ (expected - mean)**2)
    measurement_variance = float(filt._w @ sigma**2)
    _, innovations = filt.update_with_innovations(Observation(0.0, carotenoid=0.7, **{channel: -1.0}))
    innovation, = innovations
    assert innovation.predicted == pytest.approx(mean)
    assert innovation.state_variance == pytest.approx(state_variance)
    assert innovation.measurement_variance == pytest.approx(measurement_variance)
    assert innovation.total_variance == pytest.approx(state_variance + measurement_variance)
    assert measurement_variance != pytest.approx((getattr(filt.priors, f"{channel}_rel_sigma") * mean)**2)


@pytest.mark.parametrize("channel", ["od", "rfu"])
@pytest.mark.parametrize("observed", [-1.0, 0.0, 4.0])
def test_likelihood_is_a_normalized_density_for_each_particle(channel, observed):
    filt = _filter()
    expected = _expected(filt, channel)
    sigma = getattr(filt.priors, f"{channel}_rel_sigma") * np.maximum(np.abs(expected), 1e-6)
    observation = Observation(0.0, carotenoid=0.7, **{channel: observed})
    np.testing.assert_allclose(filt._log_likelihood(filt._x, observation),
                               norm.logpdf(observed, loc=expected, scale=sigma), rtol=1e-14)


def test_joint_likelihood_sums_the_two_conditionally_independent_normalized_channels():
    filt = _filter()
    observation = Observation(0.0, od=2.5, rfu=-0.4, carotenoid=0.7)
    expected_logpdf = np.zeros(filt.n_particles)
    for channel in observation.channels:
        expected = _expected(filt, channel)
        sigma = getattr(filt.priors, f"{channel}_rel_sigma") * expected
        expected_logpdf += norm.logpdf(getattr(observation, channel), loc=expected, scale=sigma)
    np.testing.assert_allclose(filt._log_likelihood(filt._x, observation), expected_logpdf)


@pytest.mark.parametrize("channel", ["od", "rfu"])
def test_each_particles_likelihood_integrates_to_one_over_possible_readings(channel):
    filt = _filter()
    expected = _expected(filt, channel)
    sigmas = getattr(filt.priors, f"{channel}_rel_sigma") * expected
    for index, (mean, sigma) in enumerate(zip(expected, sigmas)):
        def density_in_standard_units(z):
            observation = Observation(0.0, carotenoid=0.7, **{channel: float(mean + sigma * z)})
            return float(np.exp(filt._log_likelihood(filt._x, observation)[index]) * sigma)

        integral, _ = quad(density_in_standard_units, -12.0, 12.0, epsabs=1e-10)
        assert integral == pytest.approx(1.0, abs=1e-9)


def test_equal_standardized_errors_do_not_give_equal_weights_when_sigmas_differ():
    priors = replace(PRIORS, od_rel_sigma=0.5)
    filt = ParticleFilter(priors, OPTICS, 1.0, 0.0, n_particles=2, seed=0, resample_threshold=0.0)
    filt._x[:, 0] = [1.0, 2.0]
    observation = Observation(0.0, od=4.0 / 3.0)
    posterior = filt.update(observation)
    np.testing.assert_allclose(filt._w, [2.0 / 3.0, 1.0 / 3.0], rtol=1e-14)
    assert posterior.posterior_ess == pytest.approx(1.8)
    assert posterior.mean["biomass"] == pytest.approx(4.0 / 3.0)


@pytest.mark.parametrize("channel", ["od", "rfu"])
def test_the_explicit_numerical_scale_floor_handles_zero_expected_signal(channel):
    priors = replace(PRIORS, observation_scale_floor=1e-4)
    optics = ReporterOptics(1.0, 0.0, 0.0, 0.0)
    filt = ParticleFilter(priors, optics, 1.0, 0.0, n_particles=3, seed=0, resample_threshold=0.0)
    filt._x[:, 0], filt._x[:, 1] = [0.0, 1e-8, 1.0], 1.0
    filt._w = np.array([0.2, 0.5, 0.3])
    expected = np.array([0.0, 1e-8, 1.0])
    sigma = getattr(priors, f"{channel}_rel_sigma") * np.array([1e-4, 1e-4, 1.0])
    obs = Observation(0.0, **{channel: -1e-5})
    np.testing.assert_allclose(filt._log_likelihood(filt._x, obs),
                               norm.logpdf(-1e-5, loc=expected, scale=sigma))
    moment, = filt._predictive(filt._x, obs)
    assert moment["measurement_variance"] == pytest.approx(float(filt._w @ sigma**2))
    assert TwinPriors.__dataclass_fields__["observation_scale_floor"].default == 1e-6


@pytest.mark.parametrize("value", [0.0, -1e-6, np.nan, np.inf, -np.inf])
def test_a_nonpositive_or_nonfinite_numerical_floor_is_refused(value):
    with pytest.raises(ValueError, match="observation_scale_floor"):
        replace(PRIORS, observation_scale_floor=value)


def test_an_undeclared_noise_mode_is_refused():
    with pytest.raises(ValueError, match="observation_noise_mode"):
        replace(PRIORS, observation_noise_mode="automatic")


@pytest.mark.parametrize("channel", ["od", "rfu"])
@pytest.mark.parametrize("observed", [-2.0, 0.0, 8.0])
def test_observed_scale_legacy_explicitly_freezes_the_old_unnormalized_weighting(channel, observed):
    filt = _filter(observation_noise_mode="observed_scale_legacy")
    expected = _expected(filt, channel)
    sigma = getattr(filt.priors, f"{channel}_rel_sigma") * max(abs(observed), 1e-6)
    observation = Observation(0.0, carotenoid=0.7, **{channel: observed})
    np.testing.assert_array_equal(filt._log_likelihood(filt._x, observation),
                                  -0.5 * ((observed - expected) / sigma)**2)
    moment, = filt._predictive(filt._x, observation)
    assert moment["measurement_variance"] == sigma**2
    assert PRIORS.observation_noise_mode == "predicted_scale"


@pytest.mark.parametrize("mode", ["predicted_scale", "observed_scale_legacy"])
def test_diagnostic_collection_preserves_weights_particles_ancestry_and_rng(mode):
    plain, diagnosed = (_filter(observation_noise_mode=mode) for _ in range(2))
    for filt in (plain, diagnosed):
        filt.resample_threshold = 1.0
    for observation in (Observation(0.0, od=2.0, rfu=3.0), Observation(0.5),
                        Observation(1.0, od=-0.2, rfu=-0.1)):
        expected = plain.update(observation)
        actual, _ = diagnosed.update_with_innovations(observation)
        assert actual == expected
        np.testing.assert_array_equal(plain._x, diagnosed._x)
        np.testing.assert_array_equal(plain._w, diagnosed._w)
        np.testing.assert_array_equal(plain._ancestors, diagnosed._ancestors)
        assert plain._rng.bit_generator.state == diagnosed._rng.bit_generator.state
        assert plain.n_resamples == diagnosed.n_resamples
    assert plain.n_resamples > 0


@pytest.mark.parametrize("dt", [24.0, 1.0, 1.0 / 6.0])
@pytest.mark.parametrize("positive_noise_mode", ["mean_preserving", "median_preserving_legacy"])
def test_observation_noise_choice_cannot_change_24h_process_noise_or_random_draws(dt, positive_noise_mode):
    current = _filter(positive_noise_mode=positive_noise_mode)
    legacy = _filter(positive_noise_mode=positive_noise_mode, observation_noise_mode="observed_scale_legacy")
    a, b = current._x.copy(), legacy._x.copy()
    for _ in range(round(24.0 / dt)):
        a = current._propagate(a, dt)
        b = legacy._propagate(b, dt)
    np.testing.assert_array_equal(a, b)
    assert current._rng.bit_generator.state == legacy._rng.bit_generator.state
