"""Particle filter over the joint biomass / reporter / promoter-activity / growth state.

Promoter activity is a state with a posterior, not a number read off RFU/OD. Growth is inferred
jointly, since the dilution term depends on it. A missing channel contributes no likelihood term
and is never zero-filled, which would assert a dark cell rather than an absent measurement.

A particle filter rather than a Kalman variant because activity cannot be negative, and a
Gaussian approximation there puts mass on impossible states.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .calib.od import ODCalibration
from .observation import ReporterOptics, observe_od, observe_rfu

__all__ = ["Innovation", "Observation", "ParticleFilter", "Posterior", "TwinPriors"]


@dataclass(frozen=True)
class Innovation:
    """What the filter predicted, what arrived, and how sure it claimed to be.

    The residual between a prediction and the reading that follows it can be checked
    without knowing the latent truth. Correct predictive first and second moments give
    its standardised square expectation one. A Gaussian predictive distribution also
    gives a chi-squared(1) reference, but a particle mixture is not generally Gaussian
    and sequential residuals need not be independent. Departures diagnose mismatch:
    above one the filter claims more certainty than it earns, below one it is hedging.

    That test -- Normalized Innovation Squared, from target tracking -- is what makes
    a posterior falsifiable on real plates where no ground truth exists. It needs no
    held-out data and no producing strain.

    Args:
        time_h: When the reading was taken.
        channel: ``"od"`` or ``"rfu"``.
        observed: The instrument reading.
        predicted: Prior predictive mean over the propagated particles, before the
            reading was used. Using the posterior mean here would compare the reading
            with a quantity fitted to it and always look good.
        state_variance: Spread of the particles' predicted readings -- what the filter
            thinks it does not know about the state.
        measurement_variance: What the filter assumes about the instrument.
        effective_sample_size: Compatibility name for ``posterior_ess``: weight ESS
            AFTER conditioning on all measured channels, BEFORE resampling. The old
            prior-only diagnostic could not reveal a collapse followed by resampling.
        prior_ess: Weight ESS of the ensemble used to predict this reading.
        post_resample_ess: Weight ESS after the optional resample. Uniform weights do
            not imply distinct particles, so this is not a diversity certificate.
        resampled: Whether this update actually resampled, not whether ESS equals N.
        prior_unique_particles: Exact distinct state vectors used for the prediction.
        unique_particles: Exact distinct state vectors after the update/resample.
        unique_ancestors: Distinct initial particle labels still represented. Process
            noise can restore state spread without restoring lost ancestry.

    New diagnostic fields default to ``None`` for manually constructed legacy records;
    the filter always supplies them. Neither weight ESS nor genealogy alone establishes
    the validity of a Gaussian predictive interval or an independent chi-squared null.
    """

    time_h: float
    channel: str
    observed: float
    predicted: float
    state_variance: float
    measurement_variance: float
    effective_sample_size: float
    prior_ess: float | None = None
    post_resample_ess: float | None = None
    resampled: bool = False
    prior_unique_particles: int | None = None
    unique_particles: int | None = None
    unique_ancestors: int | None = None

    @property
    def posterior_ess(self) -> float:
        """Weight ESS after conditioning, before any resampling."""
        return self.effective_sample_size

    @property
    def total_variance(self) -> float:
        """The filter's own claim about how far off it might be."""
        return self.state_variance + self.measurement_variance

    @property
    def nis(self) -> float:
        """Squared residual / prior variance; expectation one under correct moments."""
        total = self.total_variance
        if not np.isfinite(total) or total <= 0:
            return float("nan")
        return float((self.observed - self.predicted) ** 2 / total)

    @property
    def standardised(self) -> float:
        """Residual in units of the filter's own claimed spread.

        Whiteness is tested on this, not on ``nis``: squaring discards the sign, and a
        filter whose residuals are all the same sign is missing a term in its dynamics
        even when their magnitudes are perfectly calibrated.
        """
        total = self.total_variance
        if not np.isfinite(total) or total <= 0:
            return float("nan")
        return float((self.observed - self.predicted) / np.sqrt(total))

_STATES = ("biomass", "reporter", "promoter_activity", "growth_rate")


@dataclass(frozen=True)
class Observation:
    """What was measured at one time. ``None`` means not measured.

    Finite signed OD/RFU values, including negative blank-corrected readings, are
    valid under additive Gaussian observation noise. They are not physical biomass
    or a denominator for a concentration ratio. Nonfinite readings must be handled
    explicitly by the caller; :meth:`ParticleFilter.update` refuses them.
    """

    time_h: float
    od: float | None = None
    rfu: float | None = None
    carotenoid: float = 0.0

    @property
    def channels(self) -> tuple[str, ...]:
        return tuple(n for n in ("od", "rfu") if getattr(self, n) is not None)


@dataclass(frozen=True)
class TwinPriors:
    """Initial beliefs and process noise for the state.

    Each entry is ``(mean, standard deviation)``. The random-walk scales say how
    fast a quantity is allowed to move between observations; they are the knob that
    trades responsiveness against smoothing, and they are stated rather than tuned
    against the outcome being tested.

    ``observation_noise_mode="predicted_scale"`` declares conditionally independent
    additive Gaussian OD/RFU readings given a particle. Each channel has expected
    reading ``m_i`` from the calibrated observation function and standard deviation
    ``sigma_i = rel_sigma * max(abs(m_i), observation_scale_floor)``. The likelihood
    includes ``-log(sigma_i) - log(2*pi)/2``; dropping that term would favour particles
    merely for predicting wider noise. Predictive measurement variance is the prior
    weighted mean of ``sigma_i**2``, not noise evaluated at the observed outcome or
    at the ensemble mean. The relative scales are assumptions, not fitted NIS knobs.

    ``observation_scale_floor`` is a numerical reference-scale floor, expressed in
    each channel's reading units, not a measured instrument noise floor. Its default
    1e-6 keeps zero expected signals nondegenerate without estimating reader noise.
    ``observed_scale_legacy`` instead uses the absolute observed value and omits the
    Gaussian normalization, preserving the old particle weights bit for bit when
    the other historical conventions are also selected. It is an outcome-dependent
    pseudo-likelihood, not a coherent pre-observation probability distribution.

    ``positive_noise_mode="mean_preserving"`` is the default. Activity and biomass
    receive ``exp(N(-sigma**2 * dt / 2, sigma * sqrt(dt)))`` multipliers. Their
    expectation is one and their variance is ``exp(sigma**2 * dt) - 1``. This removes
    noise-induced drift in the physical-state mean, not biological growth or the
    separate uncertainty induced by the growth-rate walk.

    ``positive_noise_mode="median_preserving_legacy"`` explicitly selects the old
    zero-log-mean multiplier. It preserves the median but inflates the physical mean
    by ``exp(sigma**2 * t / 2)`` and has variance
    ``exp(sigma**2 * t) * (exp(sigma**2 * t) - 1)``. At the default scales, the
    analytic 24 h mean inflation is about 31% for activity and 0.5% for biomass.
    Those are consequences of the distribution, not measured biological effects.
    ``tests/test_estimator_contracts.py`` checks both moments at multiple cadences.

    The change of default requires regeneration of estimator-derived results. The
    legacy mode names only the noise convention: it does not restore the former
    clipped reporter-loss equation or invalid ESS diagnostics.

    ``growth_deceleration`` is the one term here that is not noise. A batch culture
    slows systematically as substrate depletes, and a zero-mean walk on the growth
    rate models that as noise: the filter cannot anticipate the slowing, so it is
    dragged down by each reading in turn and its residuals stay time-correlated. This
    is the rate, in 1/h, at which the specific growth rate relaxes toward zero between
    observations. It defaults to zero, which is the pre-drift filter exactly.

    Three forms were considered and two rejected.

    * *Logistic deceleration*, the textbook batch model, makes the growth rate a
      deterministic function of biomass approaching a carrying capacity. That
      contradicts growth being a free inferred state -- the filter would assert the
      deceleration rather than track it -- and the carrying capacity is a per-well
      constant nothing in this project measures.
    * *Drift proportional to the inferred stress activity* is the project's own
      thesis and the most interesting candidate, but the activity state is a
      per-gDCW synthesis rate whose scale is set by the prior, so its coefficient
      has no portable units: the same coefficient means different things on two
      plates read at different gain. It also feeds the reporter channel into the
      biomass dynamics, so a mis-specified optics model would corrupt the OD
      prediction -- exactly the fault the innovation diagnostic exists to detect,
      made invisible.
    * *Relaxation toward a lower asymptote* generalises what is used here, and the
      asymptote is a second constant nothing measures. Zero is the asymptote a
      closed batch culture actually has, since with no feed the rate goes to zero
      when the substrate does, so relaxing to zero is the choice that refuses an
      invented number rather than the choice that is convenient.

    WHAT IS REPRODUCIBLE ABOUT THE SWEEP BEHIND THE DEFAULT, AND WHAT IS NOT. This
    docstring used to report a sweep of this term across two dozen real wells with four
    figures attached: an OD median NIS falling by a quarter, half the wells in band
    becoming three fifths, the reporter median NIS halving, and the reporter's lag-1
    autocorrelation not moving in the third decimal. Those four numbers are withdrawn:
    nothing in this repository identifies the script, exact wells and configuration
    that produced them. A prefix of a dose-ordered well list is not a representative
    sample. The NIS runner now accepts an explicit deceleration, selection and output
    directory and records a manifest, but a newly configured run is not provenance
    for those historical figures.

    On a simulated culture that really does decelerate, the term reduces OD innovation
    autocorrelation: ``tests/test_estimator_drift.py`` tests that mechanism. This is
    synthetic evidence, not identification of a real-plate defect. Autocorrelation can
    reveal misspecification without uniquely locating it in growth, reporter optics,
    or noise. The default stays zero absent a separately supported dynamic model.
    """

    biomass: tuple[float, float]
    reporter: tuple[float, float]
    promoter_activity: tuple[float, float]
    growth_rate: tuple[float, float]
    k_deg: float = 0.0
    activity_walk: float = 0.15
    growth_walk: float = 0.06
    biomass_walk: float = 0.02
    od_rel_sigma: float = 0.02
    rfu_rel_sigma: float = 0.03
    growth_deceleration: float = 0.0
    positive_noise_mode: Literal["mean_preserving", "median_preserving_legacy"] = "mean_preserving"
    observation_noise_mode: Literal["predicted_scale", "observed_scale_legacy"] = "predicted_scale"
    observation_scale_floor: float = 1e-6

    def __post_init__(self) -> None:
        for name in _STATES:
            mean, std = getattr(self, name)
            if not np.isfinite(mean) or not np.isfinite(std) or std < 0:
                raise ValueError(f"{name} prior must have a finite mean and non-negative finite std")
            if name != "growth_rate" and mean < 0:
                raise ValueError(f"{name} prior mean must be non-negative")
        for name in ("activity_walk", "biomass_walk", "growth_walk", "k_deg"):
            value = getattr(self, name)
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative")
        for name in ("od_rel_sigma", "rfu_rel_sigma", "observation_scale_floor"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.observation_noise_mode not in ("predicted_scale", "observed_scale_legacy"):
            raise ValueError("observation_noise_mode must be predicted_scale or observed_scale_legacy")
        if self.positive_noise_mode not in ("mean_preserving", "median_preserving_legacy"):
            raise ValueError("positive_noise_mode must be mean_preserving or median_preserving_legacy")
        if not np.isfinite(self.growth_deceleration) or self.growth_deceleration < 0.0:
            raise ValueError(
                f"growth_deceleration must be a finite rate at or above zero, "
                f"got {self.growth_deceleration}; a negative value would accelerate "
                f"growth without bound, which is not a culture"
            )


@dataclass(frozen=True)
class Posterior:
    """Filtered belief at one time.

    Args:
        time_h: When the belief is held.
        mean: Weighted physical-state mean, not a log-space mean or median.
        std: Weighted standard deviation per state.
        effective_sample_size: Compatibility name for ``posterior_ess``, the ESS
            before resampling. All ESS and diversity fields have the same meaning
            as in :class:`Innovation`. A forecast does not condition or resample:
            its three ESS values are the current weights' ESS, and its distinct
            particle count describes the projected ensemble.
        n_channels_used: How many channels the update conditioned on.
    """

    time_h: float
    mean: dict[str, float]
    std: dict[str, float]
    effective_sample_size: float
    n_channels_used: int
    prior_ess: float | None = None
    post_resample_ess: float | None = None
    resampled: bool = False
    prior_unique_particles: int | None = None
    unique_particles: int | None = None
    unique_ancestors: int | None = None

    @property
    def posterior_ess(self) -> float:
        """Weight ESS after conditioning, before any resampling."""
        return self.effective_sample_size

    def interval(self, name: str, z: float = 1.96) -> tuple[float, float]:
        """Approximate central credible interval for one state."""
        return (self.mean[name] - z * self.std[name], self.mean[name] + z * self.std[name])

    def summary(self) -> str:
        return (
            f"t={self.time_h:6.2f} h  "
            f"X={self.mean['biomass']:.4f}+/-{self.std['biomass']:.4f} g/L  "
            f"mu={self.mean['growth_rate']:.3f}+/-{self.std['growth_rate']:.3f} /h  "
            f"k_syn={self.mean['promoter_activity']:.3e}+/-{self.std['promoter_activity']:.1e}  "
            f"[{self.n_channels_used} channels, ESS {self.effective_sample_size:.0f}]"
        )


class ParticleFilter:
    """Sequential Monte Carlo over the twin's observable state.

    Args:
        priors: Initial beliefs and random-walk scales.
        optics: Calibrated reporter measurement parameters.
        gdcw_per_od: Dry weight per linearised OD unit.
        od_blank: Media-only OD.
        n_particles: Sample size.
        od_calibration: Reader response. Defaults to a linear reader, which is only
            appropriate below the measured linear limit.
        seed: Random seed.
        resample_threshold: Resample when the effective sample size falls below
            this fraction of ``n_particles``.

    Prior, posterior/pre-resampling and post-resampling ESS are reported separately.
    :attr:`effective_sample_size` reads current weights; :attr:`last_reweight_ess`
    preserves the conditioning ESS even when resampling resets the weights. Exact
    particle counts and persistent initial-ancestor labels disclose the collapse
    that uniform weights alone would hide. Counting diversity draws no random numbers.
    """

    def __init__(
        self,
        priors: TwinPriors,
        optics: ReporterOptics,
        gdcw_per_od: float,
        od_blank: float,
        n_particles: int = 1000,
        od_calibration: ODCalibration | None = None,
        seed: int | None = None,
        resample_threshold: float = 0.5,
    ) -> None:
        if not np.isfinite(n_particles) or n_particles < 1 or int(n_particles) != n_particles:
            raise ValueError("n_particles must be a positive integer")
        if not np.isfinite(resample_threshold) or not 0 <= resample_threshold <= 1:
            raise ValueError("resample_threshold must be a finite fraction in [0, 1]")
        if not np.isfinite(gdcw_per_od) or gdcw_per_od <= 0:
            raise ValueError("gdcw_per_od must be finite and positive")
        self.priors = priors
        self.optics = optics
        self.gdcw_per_od = float(gdcw_per_od)
        self.calibration = od_calibration or ODCalibration(
            saturation_k=0.0, top_true_od=1.0, blank=od_blank, gdcw_per_od=gdcw_per_od
        )
        self.n_particles = int(n_particles)
        self.resample_threshold = float(resample_threshold)
        self._rng = np.random.default_rng(seed)
        self.n_resamples = 0
        # Diversity left by the last conditioning, read before any resample. Uniform
        # until the first reading arrives, which is what n_particles states.
        self.last_reweight_ess = float(self.n_particles)
        self.last_prior_ess = float(self.n_particles)
        self.last_resampled = False
        self._ancestors = np.arange(self.n_particles)
        self._time: float | None = None
        self._pending_innovations: list[Innovation] | None = None

        draw = self._rng.normal
        self._x = np.column_stack([
            np.abs(draw(*priors.biomass, self.n_particles)),
            np.abs(draw(*priors.reporter, self.n_particles)),
            np.abs(draw(*priors.promoter_activity, self.n_particles)),
            draw(*priors.growth_rate, self.n_particles),
        ])
        self._w = np.full(self.n_particles, 1.0 / self.n_particles)
        self._prior_unique_particles = self.unique_particles

    # ---- state access -------------------------------------------------

    @property
    def effective_sample_size(self) -> float:
        """Diversity of the weights as they stand. This is what the resample rule reads.

        It is NOT a degeneracy diagnostic: :meth:`_resample` sets the weights back to
        uniform, so immediately after a resample this is ``n_particles`` however
        degenerate the ensemble was a line earlier. What survives that reset is
        :attr:`last_reweight_ess`, and that is what :class:`Posterior` and
        :class:`Innovation` carry.
        """
        return float(1.0 / np.sum(self._w**2))

    @property
    def unique_particles(self) -> int:
        """Number of exact distinct state vectors, irrespective of weight uniformity."""
        return int(np.unique(self._x, axis=0).shape[0])

    @property
    def unique_ancestors(self) -> int:
        """Initial particle labels retained through all resampling events."""
        return int(np.unique(self._ancestors).size)

    def _diagnostics(self) -> dict:
        return dict(
            prior_ess=self.last_prior_ess, post_resample_ess=self.effective_sample_size,
            resampled=self.last_resampled, prior_unique_particles=self._prior_unique_particles,
            unique_particles=self.unique_particles, unique_ancestors=self.unique_ancestors,
        )

    def posterior(self, time_h: float | None = None, n_channels: int = 0) -> Posterior:
        mean = {n: float(np.sum(self._w * self._x[:, i])) for i, n in enumerate(_STATES)}
        std = {
            n: float(np.sqrt(max(np.sum(self._w * (self._x[:, i] - mean[n]) ** 2), 0.0)))
            for i, n in enumerate(_STATES)
        }
        return Posterior(
            time_h=self._time if time_h is None else time_h,
            mean=mean, std=std,
            effective_sample_size=self.last_reweight_ess,
            n_channels_used=n_channels, **self._diagnostics(),
        )

    # ---- dynamics -----------------------------------------------------

    def _positive_log_drift(self, sigma: float, dt: float) -> float:
        return -0.5 * sigma**2 * dt if self.priors.positive_noise_mode == "mean_preserving" else 0.0

    def _propagate(self, x: np.ndarray, dt: float) -> np.ndarray:
        """Step the canonical ``dR/dt = activity - (mu + k_deg) R`` equation.

        Zero and negative net loss are valid, as in ``reporter.simulate_reporter``.
        ``expm1`` and the exact zero-loss limit avoid dividing by an artificial floor.
        The log-walk comments below name both modes; the default adds the
        mean-preserving drift in log space. The growth
        endpoint approximation is unchanged, not an exact stochastic integrator.
        """
        if not np.isfinite(dt) or dt < 0:
            raise ValueError("propagation horizon must be finite and non-negative")
        if dt == 0:
            return x
        biomass, reporter, activity, mu = x.T
        scale = np.sqrt(dt)
        # Systematic deceleration first, then the zero-mean walk about it. Multiplicative,
        # so it shrinks the rate toward zero without ever changing its sign: a positive
        # rate cannot be driven negative, and a lysing well's negative rate is not driven
        # positive. exp(-r*dt/2) applied twice is exp(-r*dt), so a half step taken twice
        # is a full step with no discretisation error.
        #
        # At the zero default the factor is exactly 1.0 and multiplying by 1.0 is the
        # identity in IEEE arithmetic, so this line is bit for bit the pre-drift line.
        # The term draws no random numbers, so the generator stream is unmoved too, and a
        # zero-drift run reproduces every existing result rather than restating it.
        decelerated = mu * np.exp(-self.priors.growth_deceleration * dt)
        mu_next = decelerated + self._rng.normal(0, self.priors.growth_walk * scale, mu.shape)
        # A walk in LOG space, so the default subtracts sigma**2*dt/2 to hold the MEAN.
        # Legacy mode omits it: median-preserving, mean +31% over 24 h at this default.
        activity_next = activity * np.exp(
            self._rng.normal(0, self.priors.activity_walk * scale, activity.shape)
            + self._positive_log_drift(self.priors.activity_walk, dt)
        )
        loss = mu_next + self.priors.k_deg
        # The biomass and reporter steps hold the end-of-interval rate over the whole
        # interval, unchanged from before. Under deceleration that is the lowest rate in
        # the interval, so the biomass increment is understated by a relative term of
        # order deceleration*dt/2 -- under a percent at the plate's ten-minute cadence.
        # Left alone deliberately: integrating the decaying rate exactly would be a second
        # discretisation change riding along with the term being tested.
        # Exact step for dR/dt = k_syn - loss * R with both held over the interval.
        accumulated = np.divide(
            -np.expm1(-loss * dt), loss, out=np.full_like(loss, dt), where=loss != 0,
        )
        reporter_next = reporter * np.exp(-loss * dt) + activity_next * accumulated
        # The same log-space walk at biomass_walk, where legacy mode's uncorrected mean
        # runs +0.5% over 24 h.
        biomass_next = biomass * np.exp(mu_next * dt) * np.exp(
            self._rng.normal(0, self.priors.biomass_walk * scale, biomass.shape)
            + self._positive_log_drift(self.priors.biomass_walk, dt)
        )
        return np.column_stack([
            biomass_next,
            np.clip(reporter_next, 0.0, None),
            np.clip(activity_next, 0.0, None),
            mu_next,
        ])

    # ---- measurement --------------------------------------------------

    def _expected_measurement(self, x: np.ndarray, channel: str, carotenoid: float) -> np.ndarray:
        """Calibrated particle means; no realised OD or RFU is an input."""
        if channel == "od":
            return np.array([
                observe_od(b, 0.0, self.calibration, self.gdcw_per_od) for b in x[:, 0]
            ])
        return np.array([
            observe_rfu(r, b, carotenoid, self.optics) for b, r in x[:, :2]
        ])

    def _measurement_sigma(self, predicted: np.ndarray, channel: str) -> np.ndarray:
        """Conditional Gaussian SD, shared by likelihood and prior mixture moments."""
        relative = getattr(self.priors, f"{channel}_rel_sigma")
        return relative * np.maximum(np.abs(predicted), self.priors.observation_scale_floor)

    def _legacy_measurement_sigma(self, observed: float, channel: str) -> float:
        """Historical outcome-dependent scale; not a prior observation distribution."""
        relative = getattr(self.priors, f"{channel}_rel_sigma")
        return relative * max(abs(observed), self.priors.observation_scale_floor)

    def _log_likelihood(self, x: np.ndarray, obs: Observation) -> np.ndarray:
        total = np.zeros(x.shape[0])
        legacy = self.priors.observation_noise_mode == "observed_scale_legacy"
        for channel in obs.channels:
            predicted = self._expected_measurement(x, channel, obs.carotenoid)
            observed = getattr(obs, channel)
            if legacy:
                sigma = self._legacy_measurement_sigma(observed, channel)
                total += -0.5 * ((observed - predicted) / sigma) ** 2
            else:
                sigma = self._measurement_sigma(predicted, channel)
                total += (-0.5 * ((observed - predicted) / sigma) ** 2
                          - np.log(sigma) - 0.5 * np.log(2.0 * np.pi))
        return total

    # ---- filtering ----------------------------------------------------

    def _predictive(self, x: np.ndarray, obs: Observation) -> list[dict]:
        """Prior predictive moments per measured channel, before reweighting.

        Computed from the propagated ensemble and the current weights, so it is what
        the filter believed *before* seeing the reading. Reuses the same observation
        model the likelihood uses, so a residual cannot disagree with the weight it
        produced.

        Returns the moments rather than finished :class:`Innovation` objects because
        one field of an ``Innovation`` is not knowable yet: the effective sample size
        it carries is the diversity *after* the reading is weighted in, which has not
        happened at this point in the step. :meth:`update` stamps it on. Constructing
        the object here with a placeholder would be the older defect restated -- the
        field would hold whatever was convenient at construction time.

        By the law of total variance, the Gaussian mixture has state variance
        ``Var_w(m_i)`` plus measurement variance ``E_w(sigma_i**2)``. Neither uses
        the realised reading. The explicit ``observed_scale_legacy`` exception
        reproduces the old outcome-dependent diagnostic alongside its old weighting;
        it must not be interpreted as a genuinely prior-predictive variance.
        """
        found: list[dict] = []
        legacy = self.priors.observation_noise_mode == "observed_scale_legacy"
        for channel in obs.channels:
            predicted = self._expected_measurement(x, channel, obs.carotenoid)
            mean = float(np.sum(self._w * predicted))
            if legacy:
                sigma = self._legacy_measurement_sigma(getattr(obs, channel), channel)
                measurement_variance = float(sigma**2)
            else:
                sigma = self._measurement_sigma(predicted, channel)
                measurement_variance = float(np.sum(self._w * sigma**2))
            found.append(dict(
                time_h=obs.time_h, channel=channel, observed=float(getattr(obs, channel)), predicted=mean,
                state_variance=float(np.sum(self._w * (predicted - mean) ** 2)),
                measurement_variance=measurement_variance,
            ))
        return found

    def update_with_innovations(
        self, observation: Observation
    ) -> tuple[Posterior, list[Innovation]]:
        """Update, and also return what the filter predicted before it saw the reading.

        The diagnostic sibling of :meth:`update`. Identical filtering behaviour --
        the innovations are read off the propagated ensemble on the way past and
        change nothing -- so a run instrumented this way and a run without it produce
        the same posterior. The one field that is not read off the ensemble on the way
        past is the effective sample size: that is the diversity *after* the reading is
        weighted in, so :meth:`update` stamps it on once the reweighting has happened
        and before any resample resets it.

        Args:
            observation: The reading, with ``None`` for any channel not measured.

        Returns:
            The posterior, and one :class:`Innovation` per channel that was measured.
            An observation carrying no channels yields an empty list rather than a
            fabricated residual.
        """
        self._pending_innovations = None
        posterior = self.update(observation, _collect_innovations=True)
        found = self._pending_innovations or []
        self._pending_innovations = None
        return posterior, found

    def update(self, observation: Observation, _collect_innovations: bool = False) -> Posterior:
        """Propagate to the observation time and condition on whatever was measured.

        Args:
            observation: The reading, with ``None`` for any channel not measured.
            _collect_innovations: Internal. Set by
                :meth:`update_with_innovations`; stashes the prior predictive moments
                on the way past. Filtering is unaffected either way.

        An update carrying no channels reweights nothing, so there is no concentration
        to report and :attr:`last_reweight_ess` is set to the diversity of the weights
        as they stand -- ``n_particles`` if the previous step resampled. It is not left
        holding the previous step's value, which would attribute to this step a
        collapse that happened at another one.

        Raises:
            ValueError: if the observation predates the last update. The filter is
                a forward pass; admitting an earlier reading would use information
                the forecast at that time did not have.
        """
        if not np.isfinite(observation.time_h):
            raise ValueError("observation time must be finite")
        for channel in observation.channels:
            if not np.isfinite(getattr(observation, channel)):
                raise ValueError(f"{channel} observation must be finite, or None when missing")
        if not np.isfinite(observation.carotenoid) or observation.carotenoid < 0:
            raise ValueError("carotenoid must be finite and non-negative")
        if self._time is not None and observation.time_h < self._time:
            raise ValueError(
                f"observation at {observation.time_h} h predates the filter's "
                f"position at {self._time} h; the filter does not run backwards"
            )
        dt = 0.0 if self._time is None else observation.time_h - self._time
        self._x = self._propagate(self._x, dt)
        self._time = observation.time_h
        self.last_prior_ess = self.effective_sample_size
        self._prior_unique_particles = self.unique_particles
        self.last_resampled = False

        channels = observation.channels
        # Read off the prior predictive before the weights move. After reweighting the
        # ensemble has already been pulled toward the reading, and a residual measured
        # against that is a comparison of the reading with something fitted to it.
        moments = self._predictive(self._x, observation) if _collect_innovations else []
        if channels:
            log_w = np.log(np.clip(self._w, 1e-300, None))
            log_w += self._log_likelihood(self._x, observation)
            log_w -= log_w.max()
            w = np.exp(log_w)
            total = w.sum()
            # If every particle is incompatible, keep the prior rather than dividing by zero.
            self._w = w / total if np.isfinite(total) and total > 0 else self._w
        # The one point in the step where degeneracy is visible: the reading is in and
        # `_resample` has not yet reset the weights to uniform.
        self.last_reweight_ess = self.effective_sample_size
        if channels and self.last_reweight_ess < self.resample_threshold * self.n_particles:
            self._resample()
        if _collect_innovations:
            diagnostics = self._diagnostics()
            self._pending_innovations = [
                Innovation(**moment, effective_sample_size=self.last_reweight_ess, **diagnostics)
                for moment in moments
            ]

        return self.posterior(observation.time_h, len(channels))

    def _resample(self) -> None:
        """Systematic resampling: lower variance than multinomial at the same cost."""
        positions = (self._rng.random() + np.arange(self.n_particles)) / self.n_particles
        idx = np.searchsorted(np.cumsum(self._w), positions, side="right")
        idx = np.clip(idx, 0, self.n_particles - 1)
        self._x = self._x[idx]
        self._ancestors = self._ancestors[idx]
        self._w = np.full(self.n_particles, 1.0 / self.n_particles)
        self.n_resamples += 1
        self.last_resampled = True

    def forecast(self, horizon_h: float) -> Posterior:
        """Propagate the current belief forward with no further observations.

        Uses only the state as of the last update, which is what the firewall
        requires of a forecast.

        **A forecast leaves the filter exactly as it found it, including the random
        stream.** ``_propagate`` draws from ``self._rng``, and the ``.copy()`` below only
        protects the particle array -- until 2026-08-31 the draws were still taken from
        the shared generator, so a forecast advanced it and every subsequent ``update``
        got different numbers than it would have. Two consequences, both silent:

        * calling this changed the filtered posterior of the readings that came after it.
          On 200 particles over nine observations, forecasting once in the middle moved
          the final biomass by 6.5% and the promoter activity by 12%.
        * calling it twice on one filter returned two different answers, for a method
          whose entire contract is that it reads the state and does not touch it.

        Neither was reachable from ``src/`` -- nothing here forecasts mid-run -- and
        ``audit_determinism.py`` calls it last, which is exactly where the fault hides.
        So this was a trap for the next caller rather than a wrong number in ``outputs/``.

        The state is saved and restored rather than the draws being taken from a forked
        generator, because restoring is what makes the claim testable in one line:
        ``bit_generator.state`` before and after must be equal, and
        ``tests/test_estimator.py`` asserts exactly that. A fork would also be correct and
        would not let a test say so as directly.
        """
        if not np.isfinite(horizon_h) or horizon_h < 0:
            raise ValueError("forecast horizon must be finite and non-negative")
        if self._time is None:
            raise ValueError("cannot forecast before the first observation")
        entry_state = self._rng.bit_generator.state
        try:
            projected = self._propagate(self._x.copy(), horizon_h)
        finally:
            self._rng.bit_generator.state = entry_state
        mean = {n: float(np.sum(self._w * projected[:, i])) for i, n in enumerate(_STATES)}
        std = {
            n: float(np.sqrt(max(np.sum(self._w * (projected[:, i] - mean[n]) ** 2), 0.0)))
            for i, n in enumerate(_STATES)
        }
        return Posterior(
            time_h=self._time + horizon_h, mean=mean, std=std,
            effective_sample_size=self.effective_sample_size, n_channels_used=0,
            prior_ess=self.effective_sample_size, post_resample_ess=self.effective_sample_size,
            resampled=False, prior_unique_particles=self.unique_particles,
            unique_particles=int(np.unique(projected, axis=0).shape[0]),
            unique_ancestors=self.unique_ancestors,
        )
