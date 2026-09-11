from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import numpy as np

from ..estimator import Innovation
from . import engine
from .contracts import (
    Genotype, ObservationModel, PhysicalState, Protocol, ScientificRefusal,
    SimulationResult, finite, frozen_mapping,
)
from .params import Param

__all__ = ["PhysicalEstimate", "PhysicalParticle", "PhysicalParticleFilter"]

_CHANNELS = ("od", "rfu")


@dataclass(frozen=True)
class PhysicalParticle:
    parameters: engine.EngineParameters
    ancestor: int
    history: tuple[PhysicalState, ...]
    mature_history_mmol: tuple[float, ...]
    simulation: SimulationResult | None = field(default=None, compare=False, repr=False)

    @property
    def state(self) -> PhysicalState:
        return self.history[-1]

    @property
    def mature_mmol(self) -> float:
        return self.mature_history_mmol[-1]

    @property
    def expected_readings(self) -> Mapping[str, float]:
        if self.simulation is None:
            raise ValueError("expected readings require a positive shared-engine protocol interval")
        return frozen_mapping({name: float(self.simulation.observations.values[name][-1])
                               for name in _CHANNELS})


@dataclass(frozen=True)
class PhysicalEstimate:
    particles: tuple[PhysicalParticle, ...]
    weights: tuple[float, ...]
    prior_ess: float
    posterior_ess: float
    post_resample_ess: float
    resampled: bool
    prior_unique_particles: int
    unique_particles: int
    unique_ancestors: int
    parent_indices: tuple[int, ...]
    n_channels_used: int
    innovations: tuple[Innovation, ...]

    @property
    def time_h(self) -> float:
        return self.particles[0].state.time_h

    @property
    def n_particles(self) -> int:
        return len(self.particles)

    @property
    def conditional_validity(self) -> str:
        return (
            "Finite-ensemble inference conditional on the supplied physical states, parameters, "
            "genotype, protocol, optical calibration and independent Gaussian measurement SDs. "
            "Missing channels supply no likelihood; informative missingness is not modeled. "
            "Reweighting does not identify unrepresented biology or independently validate the model."
        )


def _ess(weights: np.ndarray) -> float:
    return float(1.0 / np.sum(weights**2))


def _unique_states(particles: tuple[PhysicalParticle, ...]) -> int:
    distinct = []
    for particle in particles:
        if particle.state not in distinct:
            distinct.append(particle.state)
    return len(distinct)


class PhysicalParticleFilter:
    def __init__(
        self,
        ensemble: Sequence[PhysicalState],
        *,
        parameters: engine.EngineParameters | Sequence[engine.EngineParameters],
        genotype: Genotype,
        observation_model: ObservationModel,
        measurement_noise: Mapping[str, float],
        weights: Sequence[float],
        initial_mature_mmol: Sequence[float] | None = None,
        seed: int = 0,
        resample_threshold: float = 0.5,
        rtol: float = 1e-6,
        atol: float | Mapping[str, float] | None = None,
        max_step_h: float = 0.1,
        method: str = engine.PINNED_METHOD,
    ) -> None:
        states = tuple(ensemble)
        if not states or any(not isinstance(state, PhysicalState) for state in states):
            raise ValueError("ensemble must contain complete PhysicalState instances")
        if not isinstance(genotype, Genotype):
            raise ValueError("genotype must be an explicit Genotype")
        if not isinstance(observation_model, ObservationModel):
            raise ScientificRefusal("an explicit, fully specified ObservationModel is required; the measurement map is unidentified")
        genes = {gene.name: gene for gene in genotype.genes}
        reporter = genes.get(observation_model.reporter) if isinstance(observation_model.reporter, str) else None
        if reporter is None or reporter.role != "reporter":
            raise ScientificRefusal("the measurement map must name a reporter gene in the supplied genotype")
        if any(set(state.expression) != set(genes) for state in states):
            raise ValueError("complete expression states for exactly the supplied genotype are required")
        if any(state.time_h != states[0].time_h for state in states):
            raise ValueError("all physical particles must share one absolute clock")
        if any(state.event_receipts != states[0].event_receipts for state in states):
            raise ValueError("all physical particles must share the same event receipts")
        parameter_sets = ((parameters,) * len(states) if isinstance(parameters, engine.EngineParameters)
                          else tuple(parameters))
        if len(parameter_sets) != len(states) or any(
            not isinstance(p, engine.EngineParameters) for p in parameter_sets
        ):
            raise ValueError("parameters must be an EngineParameters or one explicit set per particle")
        if not isinstance(measurement_noise, Mapping) or set(measurement_noise) != set(_CHANNELS):
            raise ScientificRefusal("measurement_noise requires explicit Gaussian SDs for exactly od and rfu, in OD and RFU units")
        noise = {name: finite(measurement_noise[name], f"{name} measurement SD", positive=True)
                 for name in _CHANNELS}
        self._variance = {name: finite(sd * sd, f"{name} measurement variance", positive=True)
                          for name, sd in noise.items()}
        w = np.array([finite(value, "particle weight", minimum=0.0) for value in weights])
        if w.shape != (len(states),) or not np.any(w > 0):
            raise ValueError("weights require one nonnegative value per particle and positive total mass")
        w /= w.max()
        w /= w.sum()
        mature = ((observation_model.initial_mature_mmol,) * len(states)
                  if initial_mature_mmol is None else tuple(initial_mature_mmol))
        if len(mature) != len(states):
            raise ValueError("initial_mature_mmol must provide one amount per particle")
        mature = tuple(finite(value, "initial_mature_mmol", minimum=0.0) for value in mature)
        if any(value > state.expression[reporter.name].protein_mmol
               for value, state in zip(mature, states, strict=True)):
            raise ValueError("initial mature reporter cannot exceed total reporter")
        if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        self.resample_threshold = finite(resample_threshold, "resample_threshold", minimum=0.0, maximum=1.0)
        if method not in ("BDF", "Radau"):
            raise ValueError("the shared stiff clock supports BDF or Radau, without silent fallback")
        self._solver_options = dict(
            rtol=finite(rtol, "rtol", positive=True),
            atol=frozen_mapping(atol) if isinstance(atol, Mapping) else atol,
            max_step_h=finite(max_step_h, "max_step_h", positive=True), method=method,
        )
        self.genotype = genotype
        self.observation_model = observation_model
        self.measurement_noise = frozen_mapping(noise)
        expectation_values = dict(observation_model.values)
        for name in ("od_noise_sd", "rfu_noise_sd", "missing_probability"):
            expectation_values[name] = Param.asserted(
                name, 0.0, observation_model.values[name].units,
                "Conditional expectation evaluation: synthetic instrument draws disabled, not a fitted calibration",
            )
        self._expected_model = replace(observation_model, values=expectation_values)
        self._particles = tuple(PhysicalParticle(p, i, (state,), (value,))
                                for i, (p, state, value) in enumerate(
                                    zip(parameter_sets, states, mature, strict=True)))
        self._w = w
        self._rng = np.random.default_rng(int(seed))
        self.n_resamples = 0
        self.last_result: PhysicalEstimate | None = None

    @property
    def particles(self) -> tuple[PhysicalParticle, ...]:
        return self._particles

    @property
    def weights(self) -> tuple[float, ...]:
        return tuple(float(value) for value in self._w)

    @property
    def time_h(self) -> float:
        return self._particles[0].state.time_h

    @property
    def effective_sample_size(self) -> float:
        return _ess(self._w)

    @property
    def unique_particles(self) -> int:
        return _unique_states(self._particles)

    @property
    def unique_ancestors(self) -> int:
        return len({particle.ancestor for particle in self._particles})

    def _check_protocol(self, protocol: Protocol) -> None:
        if not isinstance(protocol, Protocol):
            raise ValueError("an explicit Protocol is required")
        if protocol.times_h[0] != self.time_h:
            raise ValueError("protocol must start at the particles' current absolute clock; no backward time or gaps")

    def _propagate(self, protocol: Protocol) -> tuple[PhysicalParticle, ...]:
        propagated = []
        for particle in self._particles:
            simulation = engine.simulate(
                protocol, self.genotype, particle.state, particle.parameters,
                observation=replace(self._expected_model, initial_mature_mmol=particle.mature_mmol),
                **self._solver_options,
            )
            kernel = engine._Kernel(particle.parameters, self.genotype)
            names = simulation.diagnostics["stoichiometric_state_names"]
            samples = tuple(
                kernel.unpack(
                    float(time), np.array([simulation.truth[name][j] for name in names]),
                    tuple(event for event in simulation.final_state.event_receipts if event.time_h <= time),
                )
                for j, time in enumerate(simulation.times_h[:-1])
            ) + (simulation.final_state,)
            mature = tuple(float(value) for value in simulation.observations.metadata["mature_mmol"])
            propagated.append(replace(
                particle, history=particle.history[:-1] + samples,
                mature_history_mmol=particle.mature_history_mmol[:-1] + mature,
                simulation=simulation,
            ))
        return tuple(propagated)

    def _result(self, propagated, posterior_weights, moments, n_channels, *, resample):
        prior_ess = self.effective_sample_size
        posterior_ess = _ess(posterior_weights)
        prior_unique = _unique_states(propagated)
        if resample:
            positions = (self._rng.random() + np.arange(len(propagated))) / len(propagated)
            positions = np.minimum(positions, np.nextafter(1.0, 0.0))
            cumulative = np.cumsum(posterior_weights)
            cumulative /= cumulative[-1]
            cumulative[-1] = 1.0
            indices = tuple(int(index) for index in np.searchsorted(cumulative, positions, side="right"))
            particles = tuple(propagated[index] for index in indices)
            weights = np.full(len(particles), 1.0 / len(particles))
        else:
            indices = tuple(range(len(propagated)))
            particles, weights = propagated, posterior_weights
        diagnostics = dict(
            prior_ess=prior_ess, post_resample_ess=_ess(weights), resampled=resample,
            prior_unique_particles=prior_unique, unique_particles=_unique_states(particles),
            unique_ancestors=len({particle.ancestor for particle in particles}),
        )
        innovations = tuple(Innovation(**moment, effective_sample_size=posterior_ess, **diagnostics)
                            for moment in moments)
        return PhysicalEstimate(
            particles, tuple(float(value) for value in weights),
            posterior_ess=posterior_ess, parent_indices=indices, n_channels_used=n_channels,
            innovations=innovations, **diagnostics,
        )

    def step(self, protocol: Protocol, measurements: Mapping[str, float | None], *,
             collect_innovations: bool = True) -> PhysicalEstimate:
        self._check_protocol(protocol)
        if not isinstance(measurements, Mapping):
            raise ValueError("measurements must be a channel mapping; use None or omit missing channels")
        unknown = tuple(name for name in measurements if name not in _CHANNELS)
        if unknown:
            raise ScientificRefusal(f"unidentified measurement maps for {unknown!r}; supported channels are {_CHANNELS!r}")
        observed = {name: finite(measurements[name], f"{name} observation (use None when missing)")
                    for name in _CHANNELS if name in measurements and measurements[name] is not None}
        if not isinstance(collect_innovations, bool):
            raise ValueError("collect_innovations must be a bool")
        propagated = self._propagate(protocol)
        moments = []
        log_weights = np.full(len(propagated), -np.inf)
        positive = self._w > 0
        log_weights[positive] = np.log(self._w[positive])
        for name, reading in observed.items():
            predicted = np.array([particle.expected_readings[name] for particle in propagated])
            if not np.isfinite(predicted).all():
                raise ScientificRefusal(f"the supplied {name} measurement map has nonfinite expectations")
            mean = float(self._w @ predicted)
            with np.errstate(over="ignore", invalid="ignore"):
                variance = float(self._w @ (predicted - mean)**2)
            if not np.isfinite(mean) or not np.isfinite(variance):
                raise ScientificRefusal("prior predictive moments exceed numerical range")
            if collect_innovations:
                moments.append(dict(
                    time_h=protocol.times_h[-1], channel=name, observed=reading, predicted=mean,
                    state_variance=variance, measurement_variance=self._variance[name],
                ))
            with np.errstate(over="ignore", invalid="ignore"):
                log_weights -= 0.5 * ((reading - predicted) / self.measurement_noise[name])**2
        weights = self._w.copy()
        if observed:
            maximum = log_weights.max()
            if not np.isfinite(maximum):
                raise ScientificRefusal("no particle has numerically resolvable observation likelihood; refusing a prior-weight fallback")
            weights = np.exp(log_weights - maximum)
            weights /= weights.sum()
        resample = bool(observed) and _ess(weights) < self.resample_threshold * len(propagated)
        result = self._result(propagated, weights, moments, len(observed), resample=resample)
        self._particles = result.particles
        self._w = np.array(result.weights)
        self.n_resamples += int(result.resampled)
        self.last_result = result
        return result

    def forecast(self, protocol: Protocol, *, horizon_h: float | None = None) -> PhysicalEstimate:
        horizon = None if horizon_h is None else finite(horizon_h, "forecast horizon_h", minimum=0.0)
        self._check_protocol(protocol)
        if horizon == 0.0:
            receipts = self._particles[0].state.event_receipts
            if any(event.time_h == self.time_h and event not in receipts for event in protocol.events):
                raise ValueError("a zero-horizon forecast cannot apply a pending event; supply a positive shared-engine interval")
            propagated = self._particles
        else:
            if horizon is not None:
                stop = finite(self.time_h + horizon, "forecast stop", positive=True)
                if stop <= self.time_h or stop > protocol.times_h[-1]:
                    raise ValueError("forecast horizon must advance the clock within the supplied protocol")
                if stop != protocol.times_h[-1]:
                    protocol = replace(
                        protocol, times_h=tuple(time for time in protocol.times_h if time < stop) + (stop,),
                        controls=tuple(control for control in protocol.controls if control.time_h <= stop),
                        events=tuple(event for event in protocol.events if event.time_h <= stop),
                    )
            propagated = self._propagate(protocol)
        return self._result(propagated, self._w.copy(), [], 0, resample=False)
