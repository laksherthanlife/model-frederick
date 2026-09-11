"""The optional growth-rate deceleration term, tested on a specified synthetic culture.

This is a mechanism check, not provenance for the former docstring's real-well sweep.
The term is off by default and consumes no additional random numbers. The explicit
process- and observation-noise legacy modes preserve historical conventions without
restoring the clipped-loss bug; current process noise has separate analytic 24 h tests.
"""

import numpy as np
import pytest

from ystwin.estimator import Observation, ParticleFilter, TwinPriors
from ystwin.observation import ReporterOptics

OPTICS = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                        inner_filter_coeff=0.0)
PRIORS = TwinPriors(
    biomass=(0.05, 0.01), reporter=(4e-3, 1e-3),
    promoter_activity=(1.0e-3, 5e-4), growth_rate=(0.30, 0.08), k_deg=0.05,
)


def _filter(deceleration=0.0, seed=11, n_particles=256, **overrides):
    priors = TwinPriors(**{**PRIORS.__dict__, "growth_deceleration": deceleration, **overrides})
    return ParticleFilter(priors=priors, optics=OPTICS, gdcw_per_od=0.42, od_blank=0.09,
                          n_particles=n_particles, seed=seed)


def _pre_drift_propagate(filt, x, dt):
    """The dynamics as they stood before the deceleration term, re-stated here.

    Written out rather than pinned as literals so the comparison keeps its meaning if the
    observation model or the priors move: what is being asserted is that the *dynamics* are
    unchanged, and this is the code they were. The draws are taken from the filter's own
    generator, in the order the filter takes them, so the two streams line up.
    """
    biomass, reporter, activity, mu = x.T
    scale = np.sqrt(dt)
    mu_next = mu + filt._rng.normal(0, filt.priors.growth_walk * scale, mu.shape)
    activity_next = activity * np.exp(
        filt._rng.normal(0, filt.priors.activity_walk * scale, activity.shape)
    )
    loss = np.clip(mu_next + filt.priors.k_deg, 1e-6, None)
    steady = activity_next / loss
    reporter_next = steady + (reporter - steady) * np.exp(-loss * dt)
    biomass_next = biomass * np.exp(mu_next * dt) * np.exp(
        filt._rng.normal(0, filt.priors.biomass_walk * scale, biomass.shape)
    )
    return np.column_stack([
        np.clip(biomass_next, 1e-9, None),
        np.clip(reporter_next, 0.0, None),
        np.clip(activity_next, 0.0, None),
        mu_next,
    ])


def _observations(hours=6.0, n=25, mu=0.35):
    t = np.linspace(0.1, hours, n)
    biomass = 0.05 * np.exp(mu * t)
    return [
        Observation(time_h=float(time), od=float(b / 0.42 + 0.09), rfu=float(250.0 + 3e5 * b))
        for time, b in zip(t, biomass)
    ]


def test_the_deceleration_term_is_off_unless_it_is_asked_for():
    assert PRIORS.growth_deceleration == 0.0


def test_zero_deceleration_preserves_the_legacy_walk_draws_at_positive_loss():
    """Walks remain bit-identical only in the named legacy mode.

    The reporter's cancellation-safe exact step is numerically, not bitwise, equivalent
    at positive loss. Zero and negative loss now follow the canonical equation instead
    of the old floor; test_estimator_contracts checks those cases independently.
    """
    shipped = _filter(0.0, positive_noise_mode="median_preserving_legacy",
                      observation_noise_mode="observed_scale_legacy")
    reference = _filter(0.0, positive_noise_mode="median_preserving_legacy",
                        observation_noise_mode="observed_scale_legacy")
    state = shipped._x.copy()

    for dt in (0.1667, 0.5, 1.0):
        moved = shipped._propagate(state, dt)
        expected = _pre_drift_propagate(reference, state, dt)
        unchanged_loss = moved[:, 3] + shipped.priors.k_deg >= 1e-6
        assert unchanged_loss.any()
        np.testing.assert_array_equal(moved[:, [0, 2, 3]], expected[:, [0, 2, 3]])
        np.testing.assert_allclose(moved[unchanged_loss, 1], expected[unchanged_loss, 1], rtol=2e-14)
        assert shipped._rng.bit_generator.state == reference._rng.bit_generator.state
        state = moved


def test_the_deceleration_term_draws_no_random_numbers():
    """Why the bit-for-bit claim survives a *sequence* of steps, not just one.

    A term that consumed a draw would leave the generator one number further along, and
    every later step in the run would differ even at zero deceleration.
    """
    off, on = _filter(0.0), _filter(0.4)
    off._propagate(off._x, 0.5)
    on._propagate(on._x, 0.5)

    assert off._rng.bit_generator.state == on._rng.bit_generator.state


def test_a_larger_deceleration_leaves_the_growth_rate_lower():
    """Monotone in the parameter, compared under common random numbers so the ordering is
    the term and not the noise."""
    rates = [0.0, 0.05, 0.2, 0.6, 1.5]
    means = []
    for rate in rates:
        filt = _filter(rate)
        state = filt._propagate(filt._x, 1.0)
        means.append(float(np.mean(state[:, 3])))

    assert means == sorted(means, reverse=True)
    assert means[-1] < means[0]


def test_the_deceleration_cannot_drive_a_positive_growth_rate_negative():
    """The term is multiplicative, so it can only shrink the rate toward zero.

    An additive drift would have been simpler and is wrong: over a long enough step it
    drives the rate through zero into a culture that is shrinking because the filter said
    so, and the biomass state then reports a die-off nobody measured.

    At an absurd rate over an absurd interval the factor underflows and the rate lands on
    exactly zero. That is the arrested culture, which is the right limit; the point is that
    nothing passes through it.
    """
    filt = _filter(5.0, growth_walk=0.0)
    filt._x[:, 3] = np.linspace(1e-4, 0.8, filt.n_particles)

    for dt in (0.1, 1.0, 24.0, 1000.0):
        moved = filt._propagate(filt._x.copy(), dt)
        assert np.all(moved[:, 3] >= 0.0)
        if dt <= 24.0:
            assert np.all(moved[:, 3] > 0.0)


def test_the_deceleration_does_not_turn_a_shrinking_culture_into_a_growing_one():
    """A lysing well carries a negative rate. Relaxation toward zero must approach it from
    below, never cross it: a sign flip would turn measured death into growth."""
    filt = _filter(3.0, growth_walk=0.0)
    filt._x[:, 3] = np.linspace(-0.4, -1e-4, filt.n_particles)

    moved = filt._propagate(filt._x.copy(), 2.0)

    assert np.all(moved[:, 3] < 0.0)
    assert np.all(moved[:, 3] > filt._x[:, 3])


def test_a_half_step_taken_twice_equals_a_full_step():
    """The dt scaling, which an additive-per-step drift would get wrong.

    exp(-r*dt/2) applied twice is exp(-r*dt) exactly, so the term is a property of elapsed
    time and not of how often the plate happened to be read. The walk is silenced here
    because two half steps draw twice and a full step draws once, and that difference is
    the noise model rather than the drift.
    """
    full = _filter(0.35, growth_walk=0.0)
    halved = _filter(0.35, growth_walk=0.0)
    state = full._x.copy()

    one = full._propagate(state, 0.8)
    two = halved._propagate(halved._propagate(state, 0.4), 0.4)

    assert one[:, 3] == pytest.approx(two[:, 3], rel=1e-12)


def test_the_same_seed_gives_the_same_answer_with_the_deceleration_on():
    """The repository's determinism requirement, re-checked on the new code path."""
    runs = []
    for _ in range(2):
        filt = _filter(0.3)
        posteriors = [filt.update(obs) for obs in _observations()]
        runs.append([p.mean["growth_rate"] for p in posteriors]
                    + [filt.forecast(2.0).mean["biomass"]])

    assert runs[0] == runs[1]


def test_a_different_seed_gives_a_different_answer_with_the_deceleration_on():
    """The half of determinism people leave out. A deterministic drift applied to every
    particle alike could mask a seed that stopped being threaded."""
    finals = []
    for seed in (11, 12):
        filt = _filter(0.3, seed=seed)
        for obs in _observations():
            filt.update(obs)
        finals.append(filt.posterior().mean["growth_rate"])

    assert finals[0] != finals[1]


def test_the_deceleration_reaches_the_forecast():
    """The term has to be in the forecast or it is decoration: anticipating the slowing is
    the only thing it buys that a reweighting could not."""
    straight, decelerating = _filter(0.0), _filter(0.5)
    for filt in (straight, decelerating):
        for obs in _observations():
            filt.update(obs)

    assert decelerating.forecast(4.0).mean["growth_rate"] < \
           straight.forecast(4.0).mean["growth_rate"]


def test_a_negative_deceleration_is_refused_rather_than_clamped():
    """Naming what is wrong beats silently taking the absolute value, which would turn a
    sign error in a caller into a plausible run."""
    with pytest.raises(ValueError, match="growth_deceleration"):
        TwinPriors(**{**PRIORS.__dict__, "growth_deceleration": -0.1})


def test_the_deceleration_whitens_the_innovations_of_a_decelerating_culture():
    """An OD-only implementation check with known state and independent reader noise.

    The old fixture supplied a reporter signal hundreds of times outside its prior,
    collapsing the joint filter. Here the growth state and endpoint dynamics are
    known and noiseless: assimilation cannot hide an omitted deterministic term in
    a growth-rate walk. Only reader noise remains. This is a mechanism check, not
    identification of the defect on real plates or validation of the endpoint solver.
    """
    t = np.linspace(0.0, 8.0, 40)
    rate = 0.45 * np.exp(-0.3 * t)
    biomass = 0.05 * np.exp(np.r_[0.0, np.cumsum(np.diff(t) * rate[1:])])
    rng = np.random.default_rng(91)
    od = (biomass / 0.42 + 0.09) * (1.0 + rng.normal(0.0, PRIORS.od_rel_sigma, t.size))
    observations = [Observation(time_h=float(time), od=float(reading))
                    for time, reading in zip(t, od)]

    def lag1_od_autocorrelation(deceleration):
        filt = _filter(deceleration, n_particles=64, growth_rate=(0.45, 0.0),
                       biomass=(0.05, 0.0), growth_walk=0.0, biomass_walk=0.0, activity_walk=0.0)
        residuals = []
        for obs in observations:
            _, innovations = filt.update_with_innovations(obs)
            residuals += [i.standardised for i in innovations if i.channel == "od"]
        x = np.asarray(residuals)
        x = x - x.mean()
        return float((x[:-1] @ x[1:]) / (x @ x))

    assert abs(lag1_od_autocorrelation(0.0)) > 0.5
    assert abs(lag1_od_autocorrelation(0.3)) < 2.0 / np.sqrt(t.size)
