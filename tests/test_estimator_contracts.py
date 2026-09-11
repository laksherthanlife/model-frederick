"""Estimator contracts: conditioning, genealogy, physical noise and the forward equation."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest

from ystwin.estimator import Observation, ParticleFilter, TwinPriors
from ystwin.observation import ReporterOptics, observe_rfu
from ystwin.reporter import ReporterKinetics, simulate_reporter


OPTICS = ReporterOptics(1.0, 0.0, 0.0, 0.0)
PRIORS = TwinPriors(
    biomass=(1.0, 0.0), reporter=(2.0, 0.0), promoter_activity=(3.0, 0.0),
    growth_rate=(0.0, 0.0), growth_walk=0.0, activity_walk=0.0, biomass_walk=0.0,
)


def _filter(priors=PRIORS, *, n=64, seed=3, optics=OPTICS, threshold=0.5):
    return ParticleFilter(priors, optics, gdcw_per_od=1.0, od_blank=0.0,
                          n_particles=n, seed=seed, resample_threshold=threshold)


def _collapse_filter(observation_noise_mode="predicted_scale"):
    filt = _filter(replace(PRIORS, od_rel_sigma=1e-6, observation_noise_mode=observation_noise_mode))
    filt._x[:, 0] = np.arange(1, filt.n_particles + 1)
    return filt


@pytest.mark.parametrize("observation_noise_mode", ["predicted_scale", "observed_scale_legacy"])
def test_deliberate_collapse_is_visible_after_uniform_weights_are_restored(observation_noise_mode):
    filt = _collapse_filter(observation_noise_mode)
    posterior, innovations = filt.update_with_innovations(Observation(0.0, od=1.0))
    innovation, = innovations

    assert posterior.effective_sample_size == pytest.approx(1.0)
    assert filt.effective_sample_size == pytest.approx(filt.n_particles)
    for record in (posterior, innovation):
        assert record.prior_ess == pytest.approx(filt.n_particles)
        assert record.posterior_ess == pytest.approx(1.0)
        assert record.post_resample_ess == pytest.approx(filt.n_particles)
        assert record.resampled is True
        assert record.prior_unique_particles == filt.n_particles
        assert record.unique_particles == 1
        assert record.unique_ancestors == 1
    assert innovation.state_variance > 100.0
    assert innovation.predicted == pytest.approx(32.5)


@pytest.mark.parametrize("observation_noise_mode", ["predicted_scale", "observed_scale_legacy"])
def test_a_no_channel_update_does_not_relabel_yesterdays_collapse_as_todays(observation_noise_mode):
    filt = _collapse_filter(observation_noise_mode)
    filt.update(Observation(0.0, od=1.0))
    posterior, innovations = filt.update_with_innovations(Observation(1.0))

    assert innovations == []
    assert posterior.resampled is False
    assert posterior.prior_ess == pytest.approx(filt.n_particles)
    assert posterior.posterior_ess == pytest.approx(filt.n_particles)
    assert posterior.post_resample_ess == pytest.approx(filt.n_particles)
    assert posterior.unique_particles == 1
    assert posterior.unique_ancestors == 1


@pytest.mark.parametrize("observation_noise_mode", ["predicted_scale", "observed_scale_legacy"])
def test_process_noise_can_restore_particle_spread_but_not_lost_ancestry(observation_noise_mode):
    filt = _collapse_filter(observation_noise_mode)
    filt.update(Observation(0.0, od=1.0))
    filt.priors = replace(filt.priors, biomass_walk=0.1)
    posterior = filt.update(Observation(1.0))

    assert posterior.unique_particles == filt.n_particles
    assert posterior.unique_ancestors == 1


def test_without_resampling_posterior_and_post_resampling_ess_agree():
    filt = _filter(replace(PRIORS, biomass=(1.0, 0.4)), threshold=0.0)
    posterior, innovations = filt.update_with_innovations(Observation(0.0, od=1.0))
    assert 1.0 < posterior.posterior_ess < posterior.prior_ess
    assert posterior.posterior_ess == posterior.post_resample_ess
    assert posterior.resampled is False
    assert innovations[0].posterior_ess == posterior.posterior_ess


@pytest.mark.parametrize("observation_noise_mode", ["predicted_scale", "observed_scale_legacy"])
def test_diagnostics_forecasts_and_genealogy_do_not_spend_the_filter_rng(observation_noise_mode):
    plain, diagnosed = _collapse_filter(observation_noise_mode), _collapse_filter(observation_noise_mode)
    for observation in (Observation(0.0, od=1.0), Observation(1.0), Observation(2.0, od=1.0)):
        expected = plain.update(observation)
        actual, _ = diagnosed.update_with_innovations(observation)
        before = deepcopy(diagnosed._rng.bit_generator.state)
        diagnosed.posterior()
        diagnosed.forecast(2.0)
        assert diagnosed._rng.bit_generator.state == before
        assert plain._rng.bit_generator.state == diagnosed._rng.bit_generator.state
        np.testing.assert_array_equal(plain._x, diagnosed._x)
        np.testing.assert_array_equal(plain._w, diagnosed._w)
        assert actual == expected


@pytest.mark.parametrize("horizon", [-1.0, -1e-14, np.nan, np.inf, -np.inf])
def test_invalid_forecast_horizons_are_refused_without_mutation(horizon):
    filt = _filter()
    filt.update(Observation(0.0))
    before = deepcopy(filt._rng.bit_generator.state)
    particles = filt._x.copy()
    with pytest.raises(ValueError, match="finite|non-negative"):
        filt.forecast(horizon)
    np.testing.assert_array_equal(filt._x, particles)
    assert filt._time == 0.0
    assert filt._rng.bit_generator.state == before


@pytest.mark.parametrize("dt", [-0.1, np.nan, np.inf, -np.inf])
def test_the_internal_propagator_refuses_invalid_elapsed_times(dt):
    filt = _filter()
    with pytest.raises(ValueError, match="finite|non-negative"):
        filt._propagate(filt._x, dt)


@pytest.mark.parametrize("mu,k_deg", [(-0.3, 0.1), (-0.1, 0.1), (0.0, 0.0), (0.3, 0.1)])
def test_particle_reporter_matches_the_canonical_equation_including_nonpositive_loss(mu, k_deg):
    filt = _filter(replace(PRIORS, growth_rate=(mu, 0.0), k_deg=k_deg))
    t = np.linspace(0.0, 2.0, 13)
    canonical = simulate_reporter(t, np.full_like(t, 3.0), np.full_like(t, mu),
                                  ReporterKinetics(k_deg=k_deg), r0=2.0)
    projected = filt._propagate(filt._x, float(t[-1]))
    np.testing.assert_allclose(projected[:, 1], canonical[-1], rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("dt", [24.0, 1.0, 1.0 / 6.0])
@pytest.mark.parametrize("mode", [None, "median_preserving_legacy"])
def test_24h_activity_and_biomass_walk_moments_match_the_analytic_lognormal(dt, mode):
    """Mean and variance, not a particular draw, at three observation cadences."""
    priors = replace(PRIORS, promoter_activity=(1.0, 0.0), activity_walk=0.15,
                     biomass_walk=0.02)
    if mode is not None:
        priors = replace(priors, positive_noise_mode=mode)
    n = 100_000
    filt = _filter(priors, n=n, seed=1701)
    state = filt._x.copy()
    for _ in range(round(24.0 / dt)):
        state = filt._propagate(state, dt)
    for column, sigma in ((0, priors.biomass_walk), (2, priors.activity_walk)):
        variance_log = sigma**2 * 24.0
        mean = np.exp(variance_log / 2.0) if mode else 1.0
        variance = mean**2 * np.expm1(variance_log)
        assert np.mean(state[:, column]) == pytest.approx(mean, abs=5 * np.sqrt(variance / n))
        assert np.var(state[:, column]) == pytest.approx(variance, rel=0.06)
        assert np.all(state[:, column] > 0.0)


def test_the_legacy_mode_is_explicit_and_uses_the_same_random_draws():
    mean = _filter(replace(PRIORS, activity_walk=0.15, biomass_walk=0.02))
    legacy = _filter(replace(mean.priors, positive_noise_mode="median_preserving_legacy"))
    a, b = mean._propagate(mean._x, 24.0), legacy._propagate(legacy._x, 24.0)
    for column, sigma in ((0, mean.priors.biomass_walk), (2, mean.priors.activity_walk)):
        np.testing.assert_allclose(b[:, column] / a[:, column], np.exp(sigma**2 * 24.0 / 2.0))
    assert mean._rng.bit_generator.state == legacy._rng.bit_generator.state
    with pytest.raises(ValueError, match="positive_noise_mode"):
        replace(PRIORS, positive_noise_mode="unspecified")


def test_custom_photophysics_reaches_both_the_predictive_moments_and_likelihood():
    optics = ReporterOptics(gain=7.0, background=4.0, autofluorescence=5.0,
                            inner_filter_coeff=0.8, detector_max=25.0)
    filt = _filter(optics=optics, threshold=0.0)
    filt._x[:, 0] = np.linspace(0.5, 5.0, filt.n_particles)
    observation = Observation(0.0, rfu=20.0, carotenoid=1.5)
    predicted = np.array([observe_rfu(r, b, observation.carotenoid, optics)
                          for b, r in filt._x[:, :2]])
    sigma = filt.priors.rfu_rel_sigma * np.maximum(np.abs(predicted), 1e-6)
    np.testing.assert_allclose(filt._log_likelihood(filt._x, observation),
                               -0.5 * ((observation.rfu - predicted) / sigma)**2
                               - np.log(sigma) - 0.5 * np.log(2.0 * np.pi))
    _, innovations = filt.update_with_innovations(observation)
    assert innovations[0].predicted == pytest.approx(predicted.mean())
    assert innovations[0].state_variance == pytest.approx(predicted.var())
    assert innovations[0].measurement_variance == pytest.approx(np.mean(sigma**2))


@pytest.mark.parametrize("name", ["activity_walk", "biomass_walk", "growth_walk", "k_deg"])
@pytest.mark.parametrize("value", [-0.1, np.nan, np.inf])
def test_invalid_process_parameters_are_rejected_at_construction(name, value):
    with pytest.raises(ValueError, match=name):
        replace(PRIORS, **{name: value})
