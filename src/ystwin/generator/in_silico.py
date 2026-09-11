from __future__ import annotations

import json
import operator
from dataclasses import asdict, dataclass, field, replace

import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.special import expit

from ..reporter import ReporterKinetics, simulate_reporter

INPUT_NAMES = ("imposed_upr_stress", "imposed_oxidative_stress")
LATENT_NAMES = ("upr", "oxidative", "burden")
CONTROL_NAMES = (
    "growth_retention", "enzyme_budget_scale", "ngam_mmol_per_gdcw_h", "allocation_fraction"
)
OBSERVATION_NAMES = ("cell_density", "reporter_1", "reporter_2", "reporter_3", "reporter_4")
CONTROL_FEATURE_NAMES = (
    "1", "upr", "oxidative", "burden", "upr^2", "oxidative^2", "burden^2",
    "upr*oxidative", "upr*burden", "oxidative*burden",
)

__all__ = [
    "INPUT_NAMES", "LATENT_NAMES", "CONTROL_NAMES", "OBSERVATION_NAMES", "CONTROL_FEATURE_NAMES",
    "TeacherParameters", "SyntheticEpisode", "teacher_rhs", "teacher_controls",
    "simulate_teacher", "observe_reporters", "generate_episodes",
]


def _finite_array(value, name, shape=None):
    try:
        if np.iscomplexobj(value):
            raise ValueError(f"{name} must be real")
        array = np.asarray(value, dtype=float, order="C")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain finite real numbers") from exc
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite real numbers")
    return array


def _scalar(value, name, *, positive=False):
    result = float(_finite_array(value, name, ()))
    if result < 0.0 or (positive and result == 0.0):
        raise ValueError(f"{name} must be {'positive' if positive else 'non-negative'}")
    return result


def _integer(value, name):
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a non-negative integer") from exc
    if isinstance(value, (bool, np.bool_)) or result < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return int(result)


def _grid(times_h):
    times = _finite_array(times_h, "times_h")
    if times.ndim != 1 or times.size < 5:
        raise ValueError("times_h must be one-dimensional with at least 5 timepoints")
    if not np.all(np.diff(times) > 0.0):
        raise ValueError("times_h must be strictly increasing")
    return times


def _unit_interval(array, name):
    if np.any((array < 0.0) | (array > 1.0)):
        raise ValueError(f"{name} must lie in [0, 1]")
    return array


def _family(value):
    if value not in ("bilinear", "saturating"):
        raise ValueError("family must be 'bilinear' or 'saturating'")
    return value


@dataclass(frozen=True)
class TeacherParameters:
    family: str = "bilinear"
    tau_u_h: float = 0.8
    tau_o_h: float = 0.5
    tau_b_h: float = 2.5
    cross_uo_per_h: float = 0.12
    burden_weights: tuple[float, float] = (0.45, 0.40)
    saturation_half: float = 0.35
    saturation_power: int = 2
    reporter_basal_synthesis_per_h: tuple[float, ...] = (0.30, 0.25, 0.28, 0.22)
    reporter_loadings: tuple[tuple[float, ...], ...] = (
        (1.00, 0.10, 0.05),
        (0.05, 1.00, 0.10),
        (0.10, 0.15, 1.00),
        (0.45, 0.40, 0.50),
    )
    reporter_kinetics: ReporterKinetics = ReporterKinetics(k_deg=0.05)
    assay_growth_rate_per_h: float = 0.35
    initial_cell_density: float = 0.1
    solver_rtol: float = 1e-9
    solver_atol: float = 1e-11
    solver_max_step_h: float = 0.1
    solver_method: str = field(default="DOP853", init=False)
    default_noise_cv: float = field(default=0.001, init=False)
    noise_model: str = field(default="independent_mean_preserving_lognormal", init=False)
    control_feature_names: tuple[str, ...] = field(default=CONTROL_FEATURE_NAMES, init=False)
    growth_cost_coefficients: tuple[float, ...] = field(
        default=(0.0, 0.18, 0.14, 0.25, 0.06, 0.05, 0.08, 0.04, 0.0, 0.0), init=False
    )
    enzyme_cost_coefficients: tuple[float, ...] = field(
        default=(0.0, 0.07, 0.06, 0.22, 0.02, 0.02, 0.04, 0.02, 0.0, 0.0), init=False
    )
    ngam_coefficients: tuple[float, ...] = field(
        default=(0.7, 0.10, 0.12, 0.18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), init=False
    )
    allocation_logit_coefficients: tuple[float, ...] = field(
        default=(-1.6, -0.25, -0.20, -0.35, -0.08, -0.06, -0.12, -0.08, -0.05, -0.04),
        init=False,
    )

    def __post_init__(self):
        _family(self.family)
        for name in (
            "tau_u_h", "tau_o_h", "tau_b_h", "saturation_half", "initial_cell_density",
            "solver_rtol", "solver_atol", "solver_max_step_h",
        ):
            object.__setattr__(self, name, _scalar(getattr(self, name), name, positive=True))
        for name in ("cross_uo_per_h", "assay_growth_rate_per_h"):
            object.__setattr__(self, name, _scalar(getattr(self, name), name))
        power = _integer(self.saturation_power, "saturation_power")
        if power == 0 or self.saturation_half > 1.0 or self.saturation_half ** power == 0.0:
            raise ValueError("saturation requires power >= 1 and a representable half-power in (0, 1]")
        object.__setattr__(self, "saturation_power", power)
        weights = _finite_array(self.burden_weights, "burden_weights", (2,))
        if np.any(weights < 0.0) or weights.sum() > 1.0:
            raise ValueError("burden_weights must be non-negative and sum to at most 1")
        object.__setattr__(self, "burden_weights", tuple(weights.tolist()))
        basal = _finite_array(self.reporter_basal_synthesis_per_h, "reporter_basal_synthesis_per_h", (4,))
        loadings = _finite_array(self.reporter_loadings, "reporter_loadings", (4, 3))
        if np.any(basal <= 0.0) or np.any(loadings < 0.0) or np.linalg.matrix_rank(loadings) != 3:
            raise ValueError("reporters require positive basal synthesis and non-negative rank-3 loadings")
        object.__setattr__(self, "reporter_basal_synthesis_per_h", tuple(basal.tolist()))
        object.__setattr__(self, "reporter_loadings", tuple(tuple(row) for row in loadings.tolist()))
        if not isinstance(self.reporter_kinetics, ReporterKinetics):
            raise ValueError("reporter_kinetics must be ReporterKinetics")
        kin = self.reporter_kinetics
        loss = _scalar(kin.k_deg, "reporter k_deg", positive=True)
        mat = None if kin.k_mat is None else _scalar(kin.k_mat, "reporter k_mat", positive=True)
        immature_loss = None if kin.k_deg_immature is None else _scalar(kin.k_deg_immature, "reporter k_deg_immature")
        object.__setattr__(self, "reporter_kinetics", ReporterKinetics(loss, mat, immature_loss))

    @property
    def equations(self):
        return {
            "integration": "solve each known piecewise-linear input interval separately with solve_ivp; check the invariant at every returned knot without projection or clipping",
            "bilinear": [
                "du/dt = (input1*(1-u)-u)/tau_u_h + cross_uo_per_h*o*(1-u)",
                "do/dt = (input2*(1-o)-o)/tau_o_h",
                "db/dt = (burden_weights[0]*u + burden_weights[1]*o - b)/tau_b_h",
            ],
            "saturating": [
                "H(x) = (1 + saturation_half**saturation_power)*x**saturation_power/(saturation_half**saturation_power + x**saturation_power)",
                "du/dt = (H(input1)*(1-u)-u)/tau_u_h + cross_uo_per_h*H(o)*(1-u)",
                "do/dt = (H(input2)*(1-o)-o)/tau_o_h",
                "db/dt = (burden_weights[0]*H(u) + burden_weights[1]*H(o) - b)/tau_b_h",
            ],
            "reporter": [
                "synthesis_j = reporter_basal_synthesis_per_h[j] + reporter_loadings[j] @ latent",
                "instantaneous maturation: dR_j/dt = synthesis_j - (growth_rate_per_h + k_deg)*R_j",
                "finite maturation: dI_j/dt = synthesis_j - (growth_rate_per_h + k_mat + immature_loss)*I_j; dR_j/dt = k_mat*I_j - (growth_rate_per_h + k_deg)*R_j",
                "initial I_j and R_j are quasi-steady at the first supplied growth and synthesis",
                "synthesis and growth are piecewise linear on times_h; integrate with ystwin.reporter.simulate_reporter",
                "reporter_loadings and basal synthesis use normalized synthetic instrument units per cell per hour",
            ],
            "density": [
                "default assay growth = assay_growth_rate_per_h * growth_retention(latent)",
                "dX/dt = growth_rate_per_h(t)*X; X(t0) = initial_cell_density",
                "growth is piecewise linear on times_h; X = X(t0)*exp(cumulative_trapezoid(growth, times_h))",
                "supplied density is used unchanged before observation noise; supplied growth controls reporter dilution",
            ],
            "controls": [
                "phi = [1,u,o,b,u*u,o*o,b*b,u*o,u*b,o*b]",
                "growth_retention = exp(-phi @ growth_cost_coefficients)",
                "enzyme_budget_scale = exp(-phi @ enzyme_cost_coefficients)",
                "ngam_mmol_per_gdcw_h = phi @ ngam_coefficients",
                "allocation_fraction = sigmoid(phi @ allocation_logit_coefficients)",
                "common frozen declared synthetic control laws; not fitted to any product amount",
            ],
            "noise": [
                "observation = clean_observation * exp(sigma*epsilon - sigma*sigma/2)",
                "sigma = sqrt(log(1 + noise_cv**2)); epsilon independently standard normal by time and channel",
                "noise_cv is the public function argument; default_noise_cv is its declared default",
            ],
        }

    @property
    def provenance(self):
        return {
            "parameter_source": "declared synthetic mathematical teacher and library priors",
            "biological_validation": False,
            "product_targets_used": False,
            "calibration_data_used": False,
            "truth_available_only_in_simulation": True,
            "model_scope": {
                "purpose": "specified_synthetic_teacher_recovery",
                "imposed_inputs": list(INPUT_NAMES),
                "latent_coordinates": list(LATENT_NAMES),
                "control_outputs": list(CONTROL_NAMES),
                "state_semantics": "designed aggregate coordinates, not identified molecular signaling states",
                "virtual_sensor_panel": True,
                "full_yeast_stress_network": False,
                "control_law_source": "fixed synthetic reference coefficients, not learned biological regulation",
            },
            "parameters": asdict(self),
            "equations": self.equations,
            "library_versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        }

    def to_json(self):
        return json.dumps(self.provenance, allow_nan=False, sort_keys=True)


@dataclass(frozen=True)
class SyntheticEpisode:
    episode_id: str
    times_h: np.ndarray
    inputs: np.ndarray
    observations: np.ndarray
    latent: np.ndarray
    latent_derivative: np.ndarray
    controls: np.ndarray
    growth_rate_per_h: np.ndarray
    metadata: dict

    def __post_init__(self):
        if not isinstance(self.episode_id, str) or not self.episode_id:
            raise ValueError("episode_id must be a nonempty string")
        times = _grid(self.times_h)
        n = times.size
        for name, shape in (
            ("times_h", (n,)), ("inputs", (n, 2)), ("observations", (n, 5)),
            ("latent", (n, 3)), ("latent_derivative", (n, 3)),
            ("controls", (n, 4)), ("growth_rate_per_h", (n,)),
        ):
            array = _finite_array(getattr(self, name), name, shape).copy()
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        _unit_interval(self.inputs, "inputs")
        _unit_interval(self.latent, "latent")
        if np.any(self.observations <= 0.0) or np.any(self.growth_rate_per_h < 0.0):
            raise ValueError("observations must be positive and growth must be non-negative")
        if not isinstance(self.metadata, dict):
            raise ValueError("metadata must be a dictionary")
        object.__setattr__(self, "metadata", json.loads(json.dumps(self.metadata, allow_nan=False)))


def _parameters(parameters):
    if parameters is None:
        return TeacherParameters()
    if not isinstance(parameters, TeacherParameters):
        raise ValueError("parameters must be TeacherParameters or None")
    return parameters


def _hill(value, p):
    power = np.asarray(value) ** p.saturation_power
    half = p.saturation_half ** p.saturation_power
    return (1.0 + half) * power / (half + power)


def _rhs(latent, inputs, p):
    u, o, b = latent
    input1, input2 = inputs
    signal_u, signal_o = u, o
    if p.family == "saturating":
        input1, input2 = _hill(inputs, p)
        signal_u, signal_o = _hill(np.array([u, o]), p)
    return np.array([
        (input1 * (1.0 - u) - u) / p.tau_u_h + p.cross_uo_per_h * signal_o * (1.0 - u),
        (input2 * (1.0 - o) - o) / p.tau_o_h,
        (p.burden_weights[0] * signal_u + p.burden_weights[1] * signal_o - b) / p.tau_b_h,
    ])


def teacher_rhs(latent, inputs, parameters=TeacherParameters()) -> np.ndarray:
    state = _unit_interval(_finite_array(latent, "latent", (3,)), "latent")
    imposed = _unit_interval(_finite_array(inputs, "inputs", (2,)), "inputs")
    return _rhs(state, imposed, _parameters(parameters))


_CONTROL_PARAMETERS = TeacherParameters()


def teacher_controls(latent) -> np.ndarray:
    state = _finite_array(latent, "latent")
    if state.ndim < 1 or state.shape[-1] != 3:
        raise ValueError("latent must have shape (..., 3)")
    _unit_interval(state, "latent")
    u, o, b = np.moveaxis(state, -1, 0)
    phi = np.stack((np.ones_like(u), u, o, b, u*u, o*o, b*b, u*o, u*b, o*b), axis=-1)
    p = _CONTROL_PARAMETERS
    return np.stack((
        np.exp(-phi @ p.growth_cost_coefficients),
        np.exp(-phi @ p.enzyme_cost_coefficients),
        phi @ p.ngam_coefficients,
        expit(phi @ p.allocation_logit_coefficients),
    ), axis=-1)


def _culture_series(growth_rate_per_h, cell_density, n):
    growth = _finite_array(growth_rate_per_h, "growth_rate_per_h", (n,))
    density = _finite_array(cell_density, "cell_density", (n,))
    if np.any(growth < 0.0) or np.any(density <= 0.0):
        raise ValueError("growth_rate_per_h must be non-negative and cell_density strictly positive")
    return growth, density


def _noise_log_variance(noise_cv):
    return float(np.logaddexp(0.0, 2.0 * np.log(noise_cv))) if noise_cv else 0.0


def observe_reporters(
    times_h, latent, growth_rate_per_h, cell_density, *, seed=0, noise_cv=0.001, parameters=None
) -> np.ndarray:
    times = _grid(times_h)
    state = _unit_interval(_finite_array(latent, "latent", (times.size, 3)), "latent")
    growth, density = _culture_series(growth_rate_per_h, cell_density, times.size)
    p = _parameters(parameters)
    cv = _scalar(noise_cv, "noise_cv")
    seed = _integer(seed, "seed")
    synthesis = np.asarray(p.reporter_basal_synthesis_per_h) + state @ np.asarray(p.reporter_loadings).T
    observed = np.column_stack([
        density,
        *[simulate_reporter(times, synthesis[:, j], growth, p.reporter_kinetics) for j in range(4)],
    ])
    if cv:
        variance = _noise_log_variance(cv)
        observed *= np.random.default_rng(seed).lognormal(-variance / 2.0, np.sqrt(variance), observed.shape)
    if not np.all(np.isfinite(observed)) or np.any(observed <= 0.0):
        raise RuntimeError("reporter integration or observation noise produced nonpositive/nonfinite signals")
    return observed


def simulate_teacher(
    times_h, inputs, *, parameters=None, episode_id="episode", seed=0, noise_cv=0.001,
    family=None, initial_latent=None, growth_rate_per_h=None, cell_density=None,
) -> SyntheticEpisode:
    times = _grid(times_h)
    imposed = _unit_interval(_finite_array(inputs, "inputs", (times.size, len(INPUT_NAMES))), "inputs")
    p = _parameters(parameters)
    if family is not None:
        p = replace(p, family=_family(family))
    seed = _integer(seed, "seed")
    cv = _scalar(noise_cv, "noise_cv")
    initial = np.zeros(3) if initial_latent is None else _finite_array(initial_latent, "initial_latent", (3,))
    _unit_interval(initial, "initial_latent")
    if not isinstance(episode_id, str) or not episode_id:
        raise ValueError("episode_id must be a nonempty string")
    if cell_density is not None and growth_rate_per_h is None:
        raise ValueError("supplied cell_density requires growth_rate_per_h from the same culture")
    if growth_rate_per_h is not None:
        growth = _finite_array(growth_rate_per_h, "growth_rate_per_h", (times.size,))
        if np.any(growth < 0.0):
            raise ValueError("growth_rate_per_h must be non-negative")
        if cell_density is not None:
            growth, cell_density = _culture_series(growth, cell_density, times.size)

    def rhs(time, state):
        stress = np.array([np.interp(time, times, imposed[:, j]) for j in range(2)])
        return _rhs(state, stress, p)

    latent = np.empty((times.size, 3))
    latent[0] = initial
    for index in range(1, times.size):
        solution = solve_ivp(
            rhs, (times[index - 1], times[index]), latent[index - 1], t_eval=[times[index]],
            method=p.solver_method, rtol=p.solver_rtol, atol=p.solver_atol,
            max_step=min(p.solver_max_step_h, float(times[index] - times[index - 1])),
        )
        if not solution.success:
            raise RuntimeError(f"teacher integration failed: {solution.message}")
        if solution.y.shape != (3, 1) or not np.all(np.isfinite(solution.y)) or np.any((solution.y < 0.0) | (solution.y > 1.0)):
            raise RuntimeError("teacher integration violated the [0, 1] state invariant; no clipping is applied")
        latent[index] = solution.y[:, 0]
    derivative = np.array([_rhs(z, stress, p) for z, stress in zip(latent, imposed)])
    controls = teacher_controls(latent)
    if growth_rate_per_h is None:
        growth = p.assay_growth_rate_per_h * controls[:, 0]
    if cell_density is None:
        with np.errstate(over="raise", invalid="raise"):
            density = p.initial_cell_density * np.exp(cumulative_trapezoid(growth, times, initial=0.0))
    else:
        density = cell_density
    observed = observe_reporters(times, latent, growth, density, seed=seed, noise_cv=cv, parameters=p)
    metadata = {
        **p.provenance,
        "family": p.family,
        "family_selection": "an explicit family argument overrides parameters.family; otherwise the parameter object's family is used",
        "episode_id": episode_id,
        "seed": seed,
        "initial_latent": initial.tolist(),
        "initial_cell_density": float(density[0]),
        "input_names": list(INPUT_NAMES),
        "input_units": "normalized imposed stress levels in [0, 1]; not mM",
        "input_interpolation": "piecewise_linear",
        "latent_names": LATENT_NAMES,
        "observation_names": OBSERVATION_NAMES,
        "observation_units": [
            "normalized synthetic cell-density proxy",
            *["normalized synthetic instrument units per cell"] * 4,
        ],
        "control_names": CONTROL_NAMES,
        "control_units": ["fraction", "fraction", "mmol per gDCW per hour", "fraction"],
        "control_status": "declared synthetic control laws; no product-amount fitting",
        "growth_source": "declared_synthetic_assay" if growth_rate_per_h is None else "supplied_same_culture",
        "density_source": "integrated_growth" if cell_density is None else "supplied_same_culture",
        "supplied_culture_contract": "caller supplies aligned growth and a proportional density proxy for the same culture; numerical grid and domains checked, biological identity not inferred",
        "noise": {"model": p.noise_model, "cv": cv, "log_variance": _noise_log_variance(cv), "seed": seed},
        "split_unit": "episode_id",
    }
    return SyntheticEpisode(episode_id, times, imposed, observed, latent, derivative, controls, growth, metadata)


def _protocol_inputs(times, protocol, rng):
    duration = float(times[-1])
    if protocol == "random":
        knots = np.linspace(0.0, duration, 9)
        levels = rng.uniform(0.0, 1.0, (knots.size, 2))
        values = np.column_stack([np.interp(times, knots, levels[:, j]) for j in range(2)])
        details = {"knot_times_h": knots.tolist(), "knot_values": levels.tolist()}
    elif protocol == "pulse":
        values = np.zeros((times.size, 2))
        details = {"pulses": []}
        for j in range(2):
            start = duration * rng.uniform(0.1, 0.3)
            width = duration * rng.uniform(0.15, 0.35)
            rise = width * 0.2
            amplitude = rng.uniform(0.35, 1.0)
            knots = [0.0, start, start + rise, start + width, start + width + rise, duration]
            values[:, j] = np.interp(times, knots, [0.0, 0.0, amplitude, amplitude, 0.0, 0.0])
            details["pulses"].append({"start_h": start, "width_h": width, "rise_h": rise, "amplitude": amplitude})
    elif protocol == "ramp":
        start = rng.uniform(0.0, 0.2, 2)
        end = rng.uniform(0.6, 1.0, 2)
        values = start + (end - start) * (times / duration)[:, None]
        details = {"start": start.tolist(), "end": end.tolist()}
    else:
        amplitude = rng.uniform(0.2, 0.45, 2)
        cycles = rng.uniform(1.5, 3.5, 2)
        phase = rng.uniform(0.0, 2.0 * np.pi, 2)
        values = 0.5 + amplitude * np.sin(2.0 * np.pi * (times / duration)[:, None] * cycles + phase)
        details = {"offset": 0.5, "amplitude": amplitude.tolist(), "cycles": cycles.tolist(), "phase_radians": phase.tolist()}
    return values, details


def generate_episodes(
    n, *, seed=0, family="bilinear", protocol="random", hours=8.0, dt=0.1, noise_cv=0.001,
) -> tuple[SyntheticEpisode, ...]:
    n = _integer(n, "n")
    seed = _integer(seed, "seed")
    family = _family(family)
    if protocol not in ("random", "pulse", "ramp", "periodic"):
        raise ValueError("protocol must be random, pulse, ramp, or periodic")
    hours = _scalar(hours, "hours", positive=True)
    dt = _scalar(dt, "dt", positive=True)
    cv = _scalar(noise_cv, "noise_cv")
    if not np.isfinite(hours / dt):
        raise ValueError("hours/dt must be finite")
    times = dt * np.arange(int(np.floor(hours / dt)) + 1)
    if np.isclose(times[-1], hours, rtol=1e-12, atol=0.0):
        times[-1] = hours
    else:
        times = np.append(times, hours)
    _grid(times)
    episodes = []
    for index, child in enumerate(np.random.SeedSequence(seed).spawn(n)):
        episode_seed = int(child.generate_state(1, dtype=np.uint64)[0])
        protocol_stream, observation_stream = np.random.SeedSequence(episode_seed).spawn(2)
        protocol_seed = int(protocol_stream.generate_state(1, dtype=np.uint64)[0])
        noise_seed = int(observation_stream.generate_state(1, dtype=np.uint64)[0])
        inputs, details = _protocol_inputs(times, protocol, np.random.default_rng(protocol_seed))
        episode = simulate_teacher(
            times, inputs, episode_id=f"{family}/{protocol}/{seed}/{index}",
            seed=noise_seed, family=family, noise_cv=cv,
        )
        metadata = {
            **episode.metadata,
            "protocol": protocol,
            "protocol_parameters": details,
            "protocol_sampling": "sample the protocol on times_h, then impose piecewise-linear inputs on that grid",
            "batch_seed": seed,
            "episode_index": index,
            "episode_seed": episode_seed,
            "protocol_seed": protocol_seed,
            "seed_scheme": "batch SeedSequence children yield uint64 episode seeds; each episode seed spawns separate uint64 protocol and noise seeds",
            "requested_dt_h": dt,
        }
        episodes.append(replace(episode, metadata=metadata))
    return tuple(episodes)
