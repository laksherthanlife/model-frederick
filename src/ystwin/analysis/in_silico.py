from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import lsq_linear
from scipy.special import expit, logit
from scipy.stats import chi2


__all__ = ["InSilicoStudent", "fit_student", "score_student"]

LATENT_NAMES = ("upr", "oxidative", "burden")
CONTROL_NAMES = (
    "growth_retention",
    "enzyme_budget_scale",
    "ngam_mmol_per_gdcw_h",
    "allocation_fraction",
)
_INPUT_NAMES = ("imposed_upr_stress", "imposed_oxidative_stress")
_REPORTER_NAMES = tuple(f"reporter_{index}" for index in range(1, 5))
_HISTORY_WINDOW_H = 0.8
_HISTORY_POLYNOMIAL_DEGREE = 3
_LINK_EPSILON = 1e-6
_INITIALIZATION_TAIL_PROBABILITY = 1e-6
_INITIALIZATION_RESIDUAL_LIMIT = float(np.sqrt(chi2.isf(_INITIALIZATION_TAIL_PROBABILITY, df=3)))


def _array(value, name, shape=None):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must contain only real values")
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numeric array") from exc
    if shape is not None and result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {result.shape}")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain only finite values")
    return result


def _times(value, name="times_h", minimum=1):
    result = _array(value, name)
    if result.ndim != 1 or result.size < minimum:
        raise ValueError(f"{name} must be a 1-D array with at least {minimum} timepoints")
    if np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _observations(value, times, name="observations"):
    result = _array(value, name, (len(times), 5))
    if np.any(result[:, 0] <= 0):
        raise ValueError(f"{name} require strictly positive cell_density")
    if np.any(result[:, 1:] < 0):
        raise ValueError(f"{name} require nonnegative per-cell reporter readings")
    return result


def _inputs(value, shape):
    result = _array(value, "inputs", shape)
    if np.any(result < 0):
        raise ValueError("known imposed inputs must be nonnegative")
    return result


def _latent(value, name="latent", bounded=False):
    result = _array(value, name)
    if result.ndim == 0 or result.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3 in the order {LATENT_NAMES}")
    if bounded and np.any((result < 0) | (result > 1)):
        raise ValueError(f"{name} must lie in the invariant state domain [0, 1]")
    return result


def _controls(value, length):
    result = _array(value, "controls", (length, 4))
    if np.any(result <= 0) or np.any(result[:, [0, 1, 3]] > 1):
        raise ValueError("controls require positive values and fractions/scales no greater than 1")
    return result


def _episode_ids(episodes):
    identifiers = []
    for episode in episodes:
        identifier = getattr(episode, "episode_id", None)
        if not isinstance(identifier, str) or not identifier.strip():
            raise ValueError("each episode must have a nonempty string episode_id")
        if identifier in identifiers:
            raise ValueError(f"duplicate episode_id {identifier!r}; episodes must be unique")
        identifiers.append(identifier)
    return tuple(identifiers)


def _kinetics(episode):
    metadata = getattr(episode, "metadata", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("episode metadata must be a mapping")
    for key, names in (("latent_names", LATENT_NAMES), ("control_names", CONTROL_NAMES)):
        if key in metadata and tuple(metadata[key]) != names:
            raise ValueError(f"{key} must match the declared student coordinate order {names}")
    parameters = metadata.get("parameters", {})
    if not isinstance(parameters, Mapping):
        raise ValueError("declared teacher parameters must be a mapping")
    kinetics = parameters.get("reporter_kinetics", {})
    if not isinstance(kinetics, Mapping):
        raise ValueError("declared reporter kinetics must be a mapping")
    if (metadata.get("reporter_maturation", "instantaneous") != "instantaneous"
            or kinetics.get("k_mat") is not None):
        raise ValueError("only declared instantaneous reporter maturation is supported")
    if metadata.get("input_interpolation", "piecewise_linear") != "piecewise_linear":
        raise ValueError("known input interpolation must be piecewise_linear")
    loss = _array(
        metadata.get("reporter_degradation_per_h", kinetics.get("k_deg", 0.05)),
        "reporter degradation",
    )
    if loss.ndim == 0:
        loss = np.full(4, float(loss))
    if loss.shape != (4,) or np.any(loss < 0):
        raise ValueError("reporter degradation must be a nonnegative scalar or length-four array")
    if "k_deg" in kinetics:
        declared = _array(kinetics["k_deg"], "declared reporter degradation", ())
        if not np.all(loss == declared):
            raise ValueError("conflicting declared reporter degradation kinetics")
    return loss.copy()


def _readonly(value):
    result = np.array(value, dtype=float, copy=True)
    result.flags.writeable = False
    return result


def _scale(values):
    spread = np.std(values, axis=0)
    return np.where(spread > 1e-12, spread, 1.0)


def _quadratic(values):
    matrix = np.asarray(values, dtype=float)
    return np.column_stack([
        matrix,
        *[matrix[:, first] * matrix[:, second]
          for first, second in combinations_with_replacement(range(matrix.shape[1]), 2)],
    ])


def _quadratic_names(names):
    return tuple(names) + tuple(
        f"{first}*{second}" for first, second in combinations_with_replacement(names, 2)
    )


def _causal_features(times, observations, degradation):
    series = np.column_stack([np.log(observations[:, 0]), observations[:, 1:]])
    features = np.empty((len(times), 4))
    features[0] = degradation * observations[0, 1:]
    for index in range(1, len(times)):
        tolerance = 16 * np.finfo(float).eps * max(1.0, abs(times[index]))
        start = int(np.searchsorted(times, times[index] - _HISTORY_WINDOW_H - tolerance))
        start = min(start, index - 1)
        local_times = times[start:index + 1]
        span = times[index] - local_times[0]
        order = min(_HISTORY_POLYNOMIAL_DEGREE, max(1, (len(local_times) - 1) // 2))
        design = np.vander((local_times - times[index]) / span, N=order + 1, increasing=True)
        coefficients = np.linalg.lstsq(design, series[start:index + 1], rcond=None)[0]
        growth = coefficients[1, 0] / span
        reporter = coefficients[0, 1:]
        derivative = coefficients[1, 1:] / span
        features[index] = derivative + (growth + degradation) * reporter
    if not np.isfinite(features).all():
        raise ValueError("causal reporter inversion produced nonfinite features")
    return features


@dataclass(frozen=True)
class _LinearMap:
    coefficients: np.ndarray = field(repr=False)
    intercept: np.ndarray = field(repr=False)
    feature_mean: np.ndarray = field(repr=False)
    feature_scale: np.ndarray = field(repr=False)
    feature_names: tuple[str, ...]

    def predict(self, features):
        result = ((features - self.feature_mean) / self.feature_scale) @ self.coefficients
        return result + self.intercept

    def description(self):
        return {
            "coefficients": self.coefficients.tolist(),
            "intercept": self.intercept.tolist(),
            "feature_mean": self.feature_mean.tolist(),
            "feature_scale": self.feature_scale.tolist(),
            "feature_names": list(self.feature_names),
            "coefficient_coordinates": "training-standardized features; unpenalized target mean",
        }


def _fit_linear(features, targets, names, ridge):
    mean = features.mean(axis=0)
    scale = _scale(features)
    intercept = targets.mean(axis=0)
    design = (features - mean) / scale
    augmented = np.vstack([design, np.sqrt(len(features) * ridge) * np.eye(design.shape[1])])
    response = np.vstack([targets - intercept, np.zeros((design.shape[1], targets.shape[1]))])
    weights = np.linalg.lstsq(augmented, response, rcond=None)[0]
    return _LinearMap(
        _readonly(weights), _readonly(intercept), _readonly(mean), _readonly(scale), tuple(names),
    )


@dataclass(frozen=True)
class _InvariantDynamics:
    coefficients: np.ndarray = field(repr=False)
    feature_scale: np.ndarray = field(repr=False)
    optimality: tuple[float, ...]

    def predict(self, latent, inputs):
        base = np.concatenate([
            np.ones((*latent.shape[:-1], 1)), latent, 1.0 - latent, inputs,
        ], axis=-1)
        activation = base @ self.coefficients[:9]
        relaxation = base @ self.coefficients[9:]
        return (1.0 - latent) * activation - latent * relaxation

    def description(self):
        base_names = ("1", *LATENT_NAMES, *(f"1-{name}" for name in LATENT_NAMES), *_INPUT_NAMES)
        return {
            "coefficients": self.coefficients.tolist(),
            "feature_scale": self.feature_scale.tolist(),
            "feature_names": [f"activation:{name}" for name in base_names]
                             + [f"relaxation:{name}" for name in base_names],
            "output_names": list(LATENT_NAMES),
            "polynomial_degree": 2,
            "equation": "dx_j/dt = (1-x_j) * [1,x,1-x,u] @ W_activation[:,j] - x_j * [1,x,1-x,u] @ W_relaxation[:,j]",
            "constraint": "all rate coefficients >= 0; x in [0,1]^3, known inputs >= 0",
            "coefficient_coordinates": "physical unscaled affine rates, per hour",
            "solver_optimality": list(self.optimality),
        }


def _fit_dynamics(latent, inputs, derivative, ridge):
    base = np.column_stack([np.ones(len(latent)), latent, 1.0 - latent, inputs])
    weights = np.empty((2 * base.shape[1], 3))
    scales = np.empty_like(weights)
    optimality = []
    for coordinate in range(3):
        state = latent[:, coordinate:coordinate + 1]
        design = np.column_stack([(1.0 - state) * base, -state * base])
        scale = np.sqrt(np.mean(design ** 2, axis=0))
        scale = np.where(scale > 1e-12, scale, 1.0)
        augmented = np.vstack([
            design / scale,
            np.sqrt(len(latent) * ridge) * np.eye(design.shape[1]),
        ])
        response = np.concatenate([derivative[:, coordinate], np.zeros(design.shape[1])])
        fitted = lsq_linear(augmented, response, bounds=(0.0, np.inf), method="bvls", tol=1e-11)
        if not fitted.success:
            raise RuntimeError(f"constrained dynamics fit failed: {fitted.message}")
        weights[:, coordinate] = fitted.x / scale
        scales[:, coordinate] = scale
        optimality.append(float(fitted.optimality))
    return _InvariantDynamics(_readonly(weights), _readonly(scales), tuple(optimality))


def _control_links(controls):
    linked = np.log(controls)
    linked[:, :2] *= -1.0
    linked[:, 3] = logit(np.clip(controls[:, 3], _LINK_EPSILON, 1 - _LINK_EPSILON))
    return linked


def _inverse_control_links(linked):
    result = np.empty_like(linked)
    result[..., :2] = np.exp(-np.maximum(linked[..., :2], 0.0))
    result[..., 3] = expit(linked[..., 3])
    try:
        with np.errstate(over="raise", invalid="raise"):
            result[..., 2] = np.exp(linked[..., 2])
    except FloatingPointError as exc:
        raise ValueError("latent query extrapolates beyond finite positive control predictions") from exc
    if not np.isfinite(result).all() or np.any(result <= 0):
        raise ValueError("latent query extrapolates beyond finite positive control predictions")
    return result


@dataclass(frozen=True)
class _ControlMap:
    log_cost_coefficients: np.ndarray = field(repr=False)
    log_cost_feature_scale: np.ndarray = field(repr=False)
    positive_and_allocation: _LinearMap = field(repr=False)

    def predict(self, features):
        basis = np.column_stack([np.ones(len(features)), features])
        return np.column_stack([
            basis @ self.log_cost_coefficients,
            self.positive_and_allocation.predict(features),
        ])

    def description(self):
        other = self.positive_and_allocation
        weights = other.coefficients / other.feature_scale[:, None]
        intercept = other.intercept - other.feature_mean @ weights
        return {
            "coefficients": np.column_stack([
                self.log_cost_coefficients, np.vstack([intercept, weights]),
            ]).tolist(),
            "feature_names": ["1", *other.feature_names],
            "output_names": list(CONTROL_NAMES),
            "coefficient_coordinates": "physical unscaled quadratic features, including constant row",
            "links": ["negative_rectified_log_cost", "negative_rectified_log_cost", "log", "logit"],
            "log_cost_feature_scale": self.log_cost_feature_scale.tolist(),
            "ngam_allocation_standardized_fit": other.description(),
            "constraint": "retention log-cost coefficients >= 0; exp(-cost) is in (0,1] for latent in [0,1]^3",
            "extrapolation_rule": "negative costs from out-of-domain encoder estimates are rectified to zero; scoring reports their count and magnitude",
        }


def _fit_controller(features, controls, ridge):
    linked = _control_links(controls)
    basis = np.column_stack([np.ones(len(features)), features])
    scale = np.sqrt(np.mean(basis ** 2, axis=0))
    scale = np.where(scale > 1e-12, scale, 1.0)
    penalty = np.eye(basis.shape[1])
    penalty[0, 0] = 0.0
    design = np.vstack([basis / scale, np.sqrt(len(features) * ridge) * penalty])
    costs = np.empty((basis.shape[1], 2))
    for coordinate in range(2):
        response = np.concatenate([linked[:, coordinate], np.zeros(basis.shape[1])])
        fitted = lsq_linear(design, response, bounds=(0.0, np.inf), method="bvls", tol=1e-11)
        if not fitted.success:
            raise RuntimeError(f"constrained retention controller fit failed: {fitted.message}")
        costs[:, coordinate] = fitted.x / scale
    other = _fit_linear(features, linked[:, 2:], _quadratic_names(LATENT_NAMES), ridge)
    return _ControlMap(_readonly(costs), _readonly(scale), other)


def _rmse(prediction, target):
    return np.sqrt(np.mean((prediction - target) ** 2, axis=0))


def _excursions(latent):
    excess = np.maximum(np.maximum(-latent, latent - 1.0), 0.0)
    return {
        "count": int(np.count_nonzero(excess)),
        "max_violation": float(np.max(excess, initial=0.0)),
    }


@dataclass(frozen=True)
class _InitialState:
    raw_state: np.ndarray = field(repr=False)
    state: np.ndarray = field(repr=False)
    standardized_residual: float
    solver_optimality: float

    @property
    def accepted(self):
        return self.standardized_residual <= _INITIALIZATION_RESIDUAL_LIMIT

    def diagnostics(self):
        projection = self.state - self.raw_state
        excursion = _excursions(self.raw_state)
        result = {
            "prefix_initialization_accepted": self.accepted,
            "prefix_raw_bound_excursion_count": excursion["count"],
            "prefix_raw_max_bound_violation": excursion["max_violation"],
            "prefix_projection_count": int(np.count_nonzero(projection)),
            "prefix_projection_l2": float(np.linalg.norm(projection)),
            "prefix_max_projection": float(np.max(np.abs(projection))),
            "prefix_standardized_residual": self.standardized_residual,
            "prefix_residual_limit": _INITIALIZATION_RESIDUAL_LIMIT,
            "prefix_initialization_optimality": self.solver_optimality,
        }
        for index, name in enumerate(LATENT_NAMES):
            result[f"prefix_raw_{name}"] = float(self.raw_state[index])
            result[f"prefix_initial_{name}"] = float(self.state[index])
            result[f"prefix_projection_{name}"] = float(projection[index])
        return result


@dataclass(frozen=True)
class InSilicoStudent:
    encoder: _LinearMap = field(repr=False)
    startup_encoder: _LinearMap = field(repr=False)
    dynamics: _InvariantDynamics = field(repr=False)
    controller: _ControlMap = field(repr=False)
    reporter_degradation_per_h: np.ndarray = field(repr=False)
    training_episode_ids: tuple[str, ...]
    control_target_scale: np.ndarray = field(repr=False)
    encoder_residual_covariance: np.ndarray = field(repr=False)
    startup_residual_covariance: np.ndarray = field(repr=False)
    _metadata: dict = field(repr=False)

    @property
    def metadata(self):
        return deepcopy(self._metadata)

    def infer(self, times_h, observations) -> np.ndarray:
        times = _times(times_h)
        readings = _observations(observations, times)
        activities = _causal_features(times, readings, self.reporter_degradation_per_h)
        latent = self.encoder.predict(activities)
        latent[0] = self.startup_encoder.predict(_quadratic(readings[:1, 1:]))[0]
        if not np.isfinite(latent).all():
            raise ValueError("observations extrapolate beyond finite latent estimates")
        return latent

    def initialize(self, inferred_prefix):
        estimates = _latent(inferred_prefix, "inferred_prefix")
        if estimates.ndim != 2 or len(estimates) == 0:
            raise ValueError("inferred_prefix must have shape (positive timepoints, 3)")
        raw = estimates[-1]
        state = raw.copy()
        optimality = 0.0
        residual = 0.0
        if _excursions(raw)["count"]:
            covariance = (self.startup_residual_covariance if len(estimates) == 1
                          else self.encoder_residual_covariance)
            whitening = np.linalg.solve(np.linalg.cholesky(covariance), np.eye(3))
            fitted = lsq_linear(
                whitening, whitening @ raw, bounds=(0.0, 1.0), method="bvls", tol=1e-12,
            )
            if not fitted.success or not np.isfinite(fitted.x).all():
                raise RuntimeError(f"bounded prefix initialization failed: {fitted.message}")
            state = fitted.x
            optimality = float(fitted.optimality)
            residual = float(np.linalg.norm(whitening @ (state - raw)))
            if not np.isfinite(residual):
                raise ValueError("prefix initialization produced a nonfinite residual")
        return _InitialState(_readonly(raw), _readonly(state), residual, optimality)

    def predict_controls(self, latent) -> np.ndarray:
        states = _latent(latent)
        flat = states.reshape((-1, 3))
        linked = self.controller.predict(_quadratic(flat))
        return _inverse_control_links(linked).reshape((*states.shape[:-1], 4))

    def predict_derivative(self, latent, inputs) -> np.ndarray:
        states = _latent(latent, bounded=True)
        imposed = _inputs(inputs, (*states.shape[:-1], 2))
        return self.dynamics.predict(states, imposed)

    def forecast(self, times_h, inputs, prefix_times_h, prefix_observations) -> np.ndarray:
        times = _times(times_h)
        imposed = _inputs(inputs, (len(times), 2))
        prefix_times = _times(prefix_times_h, "prefix_times_h")
        readings = _observations(prefix_observations, prefix_times, "prefix_observations")
        if times[0] < prefix_times[0]:
            raise ValueError("forecast queries cannot precede the first prefix observation")
        prefix_end = prefix_times[-1]
        future = times > prefix_end
        if future.any() and times[0] > prefix_end:
            raise ValueError("known input grid must cover the gap from prefix end to forecast queries")
        inferred = self.infer(prefix_times, readings)
        result = np.empty((len(times), 3))
        past_indices = np.searchsorted(prefix_times, times[~future], side="right") - 1
        result[~future] = inferred[past_indices]
        if not future.any():
            return result
        initialization = self.initialize(inferred)
        if not initialization.accepted:
            excursion = _excursions(initialization.raw_state)
            raise ValueError(
                "prefix terminal encoder estimate is outside the invariant state domain [0, 1]; "
                f"max violation {excursion['max_violation']:.6g}; "
                f"standardized initialization residual {initialization.standardized_residual:.6g} "
                f"exceeds {_INITIALIZATION_RESIDUAL_LIMIT:.6g}; no forecast performed"
            )
        initial = initialization.state

        def rhs(time, state):
            known = np.array([
                np.interp(time, times, imposed[:, coordinate]) for coordinate in range(2)
            ])
            return self.dynamics.predict(state, known)

        start = prefix_end
        for index in np.flatnonzero(future):
            stop = times[index]
            solution = solve_ivp(
                rhs, (start, stop), initial, t_eval=[stop],
                method="DOP853", rtol=1e-7, atol=1e-9, max_step=min(0.1, float(stop - start)),
            )
            if not solution.success or solution.y.shape != (3, 1) or not np.isfinite(solution.y).all():
                raise RuntimeError(f"learned dynamics integration failed: {solution.message}")
            initial = solution.y[:, -1]
            excursion = _excursions(initial)
            if excursion["max_violation"] > 1e-7:
                raise RuntimeError(f"numerical forecast violated the invariant state domain: {excursion}")
            result[index] = initial
            start = stop
        return result


def fit_student(training_episodes, *, seed=0, ridge=1e-6) -> InSilicoStudent:
    ridge = float(_array(ridge, "ridge", ()))
    if ridge <= 0:
        raise ValueError("ridge must be a finite positive scalar")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    episodes = tuple(training_episodes)
    if not episodes:
        raise ValueError("training episodes cannot be empty")
    identifiers = _episode_ids(episodes)
    degradation = _kinetics(episodes[0])
    activities, initial_readings, initial_states = [], [], []
    latent_arrays, input_arrays, derivatives, controls, encoded_targets = [], [], [], [], []
    for episode in episodes:
        times = _times(episode.times_h, minimum=2)
        readings = _observations(episode.observations, times)
        imposed = _inputs(episode.inputs, (len(times), 2))
        latent = _latent(episode.latent, bounded=True)
        if latent.shape != (len(times), 3):
            raise ValueError("training latent must have shape (timepoints, 3)")
        derivative = _array(episode.latent_derivative, "latent_derivative", latent.shape)
        control = _controls(episode.controls, len(times))
        if not np.array_equal(_kinetics(episode), degradation):
            raise ValueError("all training episodes must share declared reporter degradation kinetics")
        activities.append(_causal_features(times, readings, degradation)[1:])
        encoded_targets.append(latent[1:])
        initial_readings.append(readings[0, 1:])
        initial_states.append(latent[0])
        latent_arrays.append(latent)
        input_arrays.append(imposed)
        derivatives.append(derivative)
        controls.append(control)
    encoder_features = np.concatenate(activities)
    encoder_targets = np.concatenate(encoded_targets)
    initial_features = _quadratic(np.array(initial_readings))
    initial_targets = np.array(initial_states)
    latent = np.concatenate(latent_arrays)
    imposed = np.concatenate(input_arrays)
    derivative = np.concatenate(derivatives)
    control = np.concatenate(controls)
    encoder = _fit_linear(
        encoder_features, encoder_targets,
        tuple(f"causal_activity_{index}" for index in range(1, 5)),
        ridge,
    )
    startup = _fit_linear(initial_features, initial_targets, _quadratic_names(_REPORTER_NAMES), ridge)
    dynamics = _fit_dynamics(latent, imposed, derivative, ridge)
    control_features = _quadratic(latent)
    controller = _fit_controller(control_features, control, ridge)
    control_scale = _readonly(_scale(control))
    state_fit = encoder.predict(encoder_features)
    initial_fit = startup.predict(initial_features)
    control_fit = _inverse_control_links(controller.predict(control_features))
    residual = state_fit - encoder_targets
    initial_residual = initial_fit - initial_targets
    encoder_covariance = residual.T @ residual / len(residual) + ridge * np.eye(3)
    startup_covariance = initial_residual.T @ initial_residual / len(initial_residual) + ridge * np.eye(3)
    metadata = {
        "schema_version": 2,
        "training_episode_ids": list(identifiers),
        "n_training_episodes": len(episodes),
        "n_training_timepoints": len(latent),
        "uses_product_targets": False,
        "biological_validation": False,
        "product_target_claim": "No product titres, yields, fluxes, wet-lab tables, or GEM optima are read or used as targets.",
        "input_schema": {
            "time": "strictly increasing times_h in hours",
            "observations": ["cell_density", *_REPORTER_NAMES],
            "known_imposed_inputs": list(_INPUT_NAMES),
            "observation_units": ["positive proportional cell-density proxy", *["synthetic per-cell reporter units"] * 4],
            "imposed_input_units": "nonnegative normalized stress levels, not chemical concentrations",
            "infer": ["times_h", "observations"],
            "initialize": ["inferred_prefix"],
            "predict_controls": ["latent"],
            "forecast": ["times_h", "inputs", "prefix_times_h", "prefix_observations"],
        },
        "targets": {
            "encoder": list(LATENT_NAMES),
            "dynamics": [f"d_{name}_per_h" for name in LATENT_NAMES],
            "controller": list(CONTROL_NAMES),
            "supervision": "named synthetic latent, latent_derivative, and control labels from training episodes only",
            "definitions": {
                "upr": "synthetic unfolded-protein-response activation coordinate in [0,1]",
                "oxidative": "synthetic oxidative-stress activation coordinate in [0,1]",
                "burden": "synthetic accumulated physiological-burden coordinate in [0,1]",
                "latent_derivative": "instantaneous time derivative of the three synthetic coordinates, per hour",
                "growth_retention": "fraction of reference growth capacity retained, in (0,1]",
                "enzyme_budget_scale": "fraction of reference GEM enzyme budget retained, in (0,1]",
                "ngam_mmol_per_gdcw_h": "positive non-growth-associated ATP maintenance demand, mmol per reference gDCW per hour",
                "allocation_fraction": "fraction allocated to the modeled production program, in (0,1)",
            },
        },
        "learned_models": {
            "encoder": encoder.description(),
            "startup_encoder": startup.description(),
            "dynamics": dynamics.description(),
            "controller": controller.description(),
        },
        "training_errors": {
            "evaluation": "in-sample diagnostics, not held-out evidence",
            "encoder_rmse_per_coordinate": _rmse(state_fit, encoder_targets).tolist(),
            "startup_encoder_rmse_per_coordinate": _rmse(initial_fit, initial_targets).tolist(),
            "dynamics_rmse_per_coordinate": _rmse(dynamics.predict(latent, imposed), derivative).tolist(),
            "controller_teacher_forced_rmse_per_coordinate": _rmse(control_fit, control).tolist(),
            "encoder_bound_excursions": _excursions(np.vstack([state_fit, initial_fit])),
        },
        "hyperparameters": {
            "ridge": float(ridge),
            "seed": int(seed),
            "deterministic_fit": True,
            "regularization_objective": "mean squared residual + ridge * squared standardized coefficients; intercept unpenalized",
            "history_window_h": _HISTORY_WINDOW_H,
            "history_max_polynomial_degree": _HISTORY_POLYNOMIAL_DEGREE,
            "history_degree_rule": "min(3, max(1, (available_past_points - 1) // 2))",
            "encoder_degree": 1,
            "startup_encoder_degree": 2,
            "dynamics_degree": 2,
            "controller_degree": 2,
            "control_link_boundary_epsilon": _LINK_EPSILON,
            "integrator": {
                "method": "DOP853", "rtol": 1e-7, "atol": 1e-9, "max_step_h": 0.1,
                "intervals": "split at every known piecewise-linear input knot after the prefix",
            },
            "hyperparameter_search": "none",
        },
        "initialization": {
            "method": "bounded generalized least squares / Gaussian MAP with a uniform prior on [0,1]^3",
            "objective": "minimize (state - raw_encoder_state)^T covariance^-1 (state - raw_encoder_state) on [0,1]^3",
            "covariance_source": "training encoder residual second moment + ridge * identity; separate startup residuals for one-observation prefixes",
            "encoder_residual_covariance": encoder_covariance.tolist(),
            "startup_residual_covariance": startup_covariance.tolist(),
            "covariance_diagonal_regularization": float(ridge),
            "standardized_residual_limit": _INITIALIZATION_RESIDUAL_LIMIT,
            "reference_chi_square_degrees_of_freedom": 3,
            "reference_tail_probability": _INITIALIZATION_TAIL_PROBABILITY,
            "coverage_claim": "conservative Gaussian reference guard, not calibrated held-out coverage; no validation or product-based threshold tuning",
            "rejection_rule": "forecast refuses a projection exceeding the standardized residual limit; invalid or nonfinite inputs remain errors",
            "scoring": "raw states and raw-state controls remain unprojected; prefix raw/bounded states, projection, residual and acceptance are reported; rejected forecasts have NaN errors and bounds, not zero error",
            "inference_updates_metadata": False,
        },
        "reporter_degradation_per_h": degradation.tolist(),
        "control_links": ["negative_rectified_log_cost", "negative_rectified_log_cost", "log", "logit"],
        "control_target_scale": control_scale.tolist(),
        "assumptions": [
            "Four nonnegative per-cell reporters, positive density, instantaneous maturation, and known declared degradation.",
            "Causal local-polynomial estimates use only samples at or before each query, never centered smoothing or future observations.",
            "Only the four balances dR/dt + (dlog(cell_density)/dt + degradation) * R enter the main encoder; no growth-rate truth is used.",
            "Estimated growth is used only for reporter dilution inversion, never as a separate encoder feature; raw and smoothed reporter levels are excluded from the main encoder.",
            "A separately learned raw-reporter quadratic startup head is a training prior: one observation cannot identify arbitrary hidden initialization.",
            "All centering/scaling, error covariances and learned coefficients are fitted on complete training episodes only; no test-time parameter fitting occurs.",
            "The dynamics hypothesis is a degree-two positive-rate activation/relaxation system, not an unrestricted polynomial.",
            "Nonnegative affine activation/relaxation rates make [0,1]^3 forward invariant for nonnegative imposed inputs.",
            "Infer and state/control scoring retain raw encoder estimates. Forecast initialization alone uses a bounded MAP estimate with training-only residual covariance and explicit projection diagnostics; incompatible prefixes are rejected.",
            "Future known imposed inputs are linearly interpolated on times_h and integrated separately at every knot; no future observations, teacher calls, or latent/control labels are available to forecast.",
            "At query times within the prefix, only the latest observation at or before that time is used; future prefix observations are not interpolated backward.",
            "Growth and enzyme retention use learned nonnegative quadratic log-costs and exp(-cost); NGAM uses a log link and allocation a logit link.",
            "Retention costs are guaranteed nonnegative on the latent cube; costs extrapolating below zero outside that domain are rectified and counted in scoring.",
            "Only allocation labels at fraction one use declared finite log odds; healthy retention one is fitted directly as zero log-cost.",
            "Redundant positive-rate basis weights are regularized representations, not uniquely identified biochemical rate constants.",
            "Whole-episode held-out scoring rejects IDs seen in training; disjoint IDs must identify genuinely distinct episodes.",
        ],
    }
    return InSilicoStudent(
        encoder, startup, dynamics, controller, _readonly(degradation),
        identifiers, control_scale, _readonly(encoder_covariance), _readonly(startup_covariance), metadata,
    )


def score_student(student, episodes, prefix_h=2.0) -> pd.DataFrame:
    prefix_h = float(_array(prefix_h, "prefix_h", ()))
    if prefix_h < 0:
        raise ValueError("prefix_h must be finite and nonnegative")
    episodes = tuple(episodes)
    if not episodes:
        raise ValueError("held-out episodes cannot be empty")
    identifiers = _episode_ids(episodes)
    overlap = set(identifiers).intersection(student.training_episode_ids)
    if overlap:
        raise ValueError(f"held-out episode IDs overlap training episodes: {sorted(overlap)}")
    rows = []
    for episode in episodes:
        times = _times(episode.times_h, minimum=2)
        observations = _observations(episode.observations, times)
        inputs = _inputs(episode.inputs, (len(times), 2))
        truth = _array(episode.latent, "scoring latent", (len(times), 3))
        true_controls = _controls(episode.controls, len(times))
        prefix = times <= times[0] + prefix_h
        future = ~prefix
        if not future.any():
            raise ValueError("prefix_h must leave a strictly positive forecast horizon in every episode")
        states = student.infer(times, observations)
        predicted_controls = student.predict_controls(states)
        control_costs = student.controller.predict(_quadratic(states))[:, :2]
        negative_costs = np.maximum(-control_costs, 0.0)
        initialization = student.initialize(states[prefix])
        if initialization.accepted:
            forecast = student.forecast(times, inputs, times[prefix], observations[prefix])
            forecast_error = _rmse(forecast[future], truth[future])
            forecast_excursion = _excursions(forecast[future])
        else:
            forecast_error = np.full(3, np.nan)
            forecast_excursion = {"count": np.nan, "max_violation": np.nan}
        state_error = _rmse(states, truth)
        control_error = _rmse(predicted_controls, true_controls)
        persistence_error = _rmse(np.broadcast_to(states[prefix][-1], truth[future].shape), truth[future])
        state_excursion = _excursions(states)
        metadata = getattr(episode, "metadata", {})
        row = {
            "episode_id": episode.episode_id,
            "family": metadata.get("family", "unspecified"),
            "protocol": metadata.get("protocol", "unspecified"),
            "prefix_h": float(prefix_h),
            "prefix_start_h": float(times[0]),
            "prefix_end_h": float(times[prefix][-1]),
            "forecast_horizon_h": float(times[-1] - times[prefix][-1]),
            "n_observations": len(times),
            "n_prefix_observations": int(prefix.sum()),
            "n_forecast_observations": int(future.sum()),
            "n_forecast_predictions": int(future.sum()) if initialization.accepted else 0,
            "forecast_rejected": not initialization.accepted,
            "state_rmse": float(np.sqrt(np.mean(state_error ** 2))),
            "control_normalized_rmse": float(np.sqrt(np.mean((control_error / student.control_target_scale) ** 2))),
            "forecast_rmse": float(np.sqrt(np.mean(forecast_error ** 2))),
            "persistence_rmse": float(np.sqrt(np.mean(persistence_error ** 2))),
            "state_bound_excursion_count": state_excursion["count"],
            "state_max_bound_violation": state_excursion["max_violation"],
            "controller_cost_rectification_count": int(np.count_nonzero(negative_costs)),
            "controller_max_negative_cost": float(np.max(negative_costs, initial=0.0)),
            "forecast_bound_excursion_count": forecast_excursion["count"],
            "forecast_max_bound_violation": forecast_excursion["max_violation"],
        }
        row.update(initialization.diagnostics())
        row.update({f"state_rmse_{name}": float(error) for name, error in zip(LATENT_NAMES, state_error)})
        row.update({f"forecast_rmse_{name}": float(error) for name, error in zip(LATENT_NAMES, forecast_error)})
        row.update({f"control_rmse_{name}": float(error) for name, error in zip(CONTROL_NAMES, control_error)})
        rows.append(row)
    return pd.DataFrame(rows)
