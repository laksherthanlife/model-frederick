from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from numbers import Real

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from ..mech.hog import HogModel, HogProtocol
from ..mech.kinetic_sbml import KineticSimulationError

__all__ = ["HogFit", "fit_hog_parameters", "evaluate_hog_fit", "NativeHogFit", "fit_native_hog_parameters"]

_HOG_OBSERVABLE = "Hog1PP_measured"
_DEFAULT_PARAMETERS = ("kv16f_1",)


@dataclass(frozen=True)
class HogFit:
    parameters: tuple[tuple[str, float], ...]
    observation_gain: float
    training_nrmse: float
    diagnostics: dict = field(repr=False)

    @property
    def parameter_values(self):
        return dict(self.parameters)


def _protocol(frame):
    doses = frame.nacl_molar.unique()
    offsets = (frame.time_model_s - frame.time_relative_s).to_numpy(dtype=float)
    if len(doses) != 1 or not len(offsets) or not np.allclose(offsets, offsets[0], rtol=0, atol=1e-8):
        raise ValueError("observations must describe one explicit dose/time protocol")
    return HogProtocol(float(doses[0]), float(offsets[0]))


def _prediction(model, times, protocol, parameters):
    grid = np.unique(np.r_[0.0, np.asarray(times, dtype=float)])
    result = model.simulate(grid, protocol=protocol, parameters=parameters, rtol=1e-8, atol=1e-11)
    values = np.interp(times, grid, result.variables[_HOG_OBSERVABLE])
    if not np.isfinite(values).all():
        raise KineticSimulationError("nonfinite Hog1 model observation")
    return values


def fit_hog_parameters(
    model: HogModel, training: pd.DataFrame, *, parameter_names=_DEFAULT_PARAMETERS,
    fold_bounds=(0.2, 5.0), max_nfev=30,
) -> HogFit:
    required = {"time_model_s", "time_relative_s", "value", "nacl_molar", "observable_id",
                "split", "genotype", "experiment_id", "unit"}
    if required - set(training.columns) or training.empty:
        raise ValueError("complete measured training observations are required")
    if not training["split"].eq("fit").all():
        raise ValueError("only explicitly assigned fit rows may enter parameter estimation")
    if not training.genotype.eq("wild_type").all():
        raise ValueError("this fitting protocol supports wild_type observations only")
    if not training.observable_id.eq(_HOG_OBSERVABLE).all() or training.experiment_id.nunique() != 1:
        raise ValueError("fit one Hog1 observation scale/experiment at a time")
    times = training.time_model_s.to_numpy(dtype=float)
    values = training.value.to_numpy(dtype=float)
    if not np.isfinite(times).all() or not np.isfinite(values).all() or np.any(times < 0):
        raise ValueError("fit times and observed values must be finite and times nonnegative")
    scale = float(np.sqrt(np.mean(values ** 2)))
    names = tuple(parameter_names)
    if not scale > 0 or len(np.unique(times)) < len(names) + 2:
        raise ValueError("insufficient independent times or observation excitation")
    source = model.parameter_values
    if not names or len(set(names)) != len(names) or any(name not in source or source[name] <= 0 for name in names):
        raise ValueError("fit parameter names must identify distinct positive source parameters")
    bounds = np.asarray(fold_bounds, dtype=float)
    if bounds.shape != (2,) or not np.isfinite(bounds).all() or not 0 < bounds[0] < 1 < bounds[1]:
        raise ValueError("fold_bounds must bracket the positive published reference")
    if isinstance(max_nfev, bool) or int(max_nfev) != max_nfev or max_nfev < 1:
        raise ValueError("max_nfev must be a positive integer")
    protocol = _protocol(training)
    cache, failures = {}, []
    lower, upper = np.log(bounds)
    step = 1e-3

    def evaluate(log_folds):
        key = tuple(float(value) for value in log_folds)
        if key not in cache:
            parameters = {name: float(source[name] * np.exp(value)) for name, value in zip(names, key)}
            try:
                predicted = _prediction(model, times, protocol, parameters)
                denominator = float(predicted @ predicted)
                gain = float(predicted @ values / denominator) if denominator > 0 else float("nan")
                if not np.isfinite(gain) or gain <= 0:
                    raise KineticSimulationError("model does not support a positive observation scale")
                residual = (gain * predicted - values) / scale
                cache[key] = (residual, gain, parameters, None)
            except KineticSimulationError as exc:
                failures.append(str(exc))
                cache[key] = (np.full(len(values), 10.0), float("nan"), parameters, str(exc))
        return cache[key]

    def residual(log_folds):
        return evaluate(log_folds)[0]

    def jacobian(log_folds):
        base = residual(log_folds)
        columns = []
        for index in range(len(names)):
            delta = step if log_folds[index] + step <= upper else -step
            shifted = np.array(log_folds, copy=True)
            shifted[index] += delta
            columns.append((residual(shifted) - base) / delta)
        return np.column_stack(columns)

    baseline = evaluate(np.zeros(len(names)))
    if baseline[3] is not None:
        raise RuntimeError(f"published reference cannot simulate the fitting protocol: {baseline[3]}")
    optimized = least_squares(
        residual, np.zeros(len(names)), jac=jacobian, bounds=(lower, upper),
        max_nfev=int(max_nfev), ftol=1e-6, xtol=1e-6, gtol=1e-6)
    residuals, gain, parameters, failure = evaluate(optimized.x)
    if failure is not None:
        raise RuntimeError(f"selected parameter fit is numerically invalid: {failure}")
    singular = np.linalg.svd(optimized.jac, compute_uv=False)
    tolerance = max(optimized.jac.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.count_nonzero(singular > tolerance))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
    diagnostics = {
        "source": "measured training observations, not synthetic state/control labels",
        "estimated_parameter_names": list(names),
        "published_parameters": {name: float(source[name]) for name in names},
        "parameter_multipliers": {name: parameters[name] / source[name] for name in names},
        "observation_gain": gain,
        "gain_source": "positive least-squares observation scale profiled using training rows only",
        "observation_unit": str(training.unit.iloc[0]),
        "native_model_units_relabelled": False,
        "baseline_training_nrmse": float(np.sqrt(np.mean(baseline[0] ** 2))),
        "training_nrmse": float(np.sqrt(np.mean(residuals ** 2))),
        "training_experiments": list(training.experiment_id.unique()),
        "training_rows": len(training), "protocol": {"nacl_molar": protocol.nacl_molar, "shock_time_s": protocol.shock_time_s},
        "optimization_converged": bool(optimized.success), "optimization_message": str(optimized.message),
        "search_fold_bounds": bounds.tolist(), "active_bounds": optimized.active_mask.tolist(),
        "objective_evaluations": len(cache), "failed_simulations": len(failures),
        "local_sensitivity_rank": rank, "local_sensitivity_singular_values": singular.tolist(),
        "local_sensitivity_condition_number": condition, "log_parameter_difference_step": step,
        "identifiability_scope": "local numerical sensitivity only; not global identifiability or a confidence interval",
        "remaining_parameters": "fixed published model priors, historically fitted to broader source data",
        "parameters_fitted_to_measurements": True,
        "mechanistic_structure_learned": False,
    }
    if names == ("kv16f_1",):
        diagnostics.update(
            learned_activation_balance_fold=parameters["kv16f_1"] / source["kv16f_1"],
            fixed_deactivation_parameter=source["kv16r_1"],
            fast_rate_timescale_estimated=False)
    return HogFit(tuple(parameters.items()), gain, diagnostics["training_nrmse"], diagnostics)


def evaluate_hog_fit(model: HogModel, fit: HogFit, observations: pd.DataFrame):
    if observations.empty or observations["split"].eq("fit").any():
        raise ValueError("evaluation requires nonempty held-out observations")
    if not observations.genotype.eq("wild_type").all() or not observations.observable_id.eq(_HOG_OBSERVABLE).all():
        raise ValueError("this evaluation supports wild_type Hog1 observations only")
    if not observations.unit.eq("relative_intensity").all():
        raise ValueError("held-out Westerns require their declared relative-intensity scale")
    rows, predictions = [], []
    for (experiment, replicate), group in observations.groupby(["experiment_id", "replicate_label"], dropna=False, sort=True):
        protocol = _protocol(group)
        times = group.time_model_s.to_numpy(dtype=float)
        values = group.value.to_numpy(dtype=float)
        if not np.isfinite(times).all() or not np.isfinite(values).all():
            raise ValueError("held-out observations and times must be finite")
        fitted = _prediction(model, times, protocol, fit.parameter_values)
        published = _prediction(model, times, protocol, {})
        if max(fitted) <= 0 or max(published) <= 0:
            raise ValueError("peak-normalized predictions require positive model response")
        fitted = fitted / max(fitted)
        published = published / max(published)
        static = (group.time_relative_s.to_numpy(dtype=float) > 0).astype(float)
        rows.append({
            "experiment_id": experiment, "replicate_label": replicate,
            "measurement_block": str(group.measurement_block.iloc[0]),
            "split": str(group["split"].iloc[0]), "nacl_molar": protocol.nacl_molar,
            "observations": len(group),
            "fitted_shape_rmse": float(np.sqrt(np.mean((fitted - values) ** 2))),
            "published_shape_rmse": float(np.sqrt(np.mean((published - values) ** 2))),
            "persistent_shape_rmse": float(np.sqrt(np.mean((static - values) ** 2))),
            "score_scope": "published peak-scaled Western time-course shape, not absolute phosphorylation fraction",
            "test_refitted": False,
        })
        for index, (_, original) in enumerate(group.iterrows()):
            predictions.append({**original.to_dict(), "fitted_prediction": fitted[index],
                                "published_prediction": published[index], "persistent_prediction": static[index]})
    return pd.DataFrame(rows), pd.DataFrame(predictions)


_NATIVE_GROUPS = (
    ("glycerol_synthesis", ("kv6_1", "kv6b_4")),
    ("passive_diffusion", ("kv13a_1",)),
    ("activation_balance", ("kv16f_1",)),
)
_NATIVE_WESTERNS = frozenset({"Hog1PP_measured", "Gpd1_measured"})
_NATIVE_CONCENTRATIONS = frozenset({"glycerol_measured", "glycerol_e"})


@dataclass(frozen=True)
class NativeHogFit:
    parameters: tuple[tuple[str, float], ...]
    gains: tuple[tuple[str, float], ...]
    training_nrmse: float
    diagnostics: dict = field(repr=False)

    @property
    def parameter_values(self):
        return dict(self.parameters)

    @property
    def observation_gains(self):
        return dict(self.gains)


def _native_training(training):
    required = {"time_model_s", "time_relative_s", "value", "nacl_molar", "observable_id",
                "split", "genotype", "medium", "experiment_id", "unit", "supplement_id",
                "measurement_provenance", "sd"}
    if (not isinstance(training, pd.DataFrame) or training.empty
            or not training.columns.is_unique or required - set(training.columns)):
        raise ValueError("complete measured native training observations are required")
    for column, value in (("split", "fit"), ("supplement_id", "s002"),
                          ("genotype", "wild_type"), ("medium", "YPD"), ("unit", "mol/L")):
        if not training[column].isin([value]).all():
            raise ValueError(f"native training requires explicitly assigned {column}={value} rows only")
    if not training.observable_id.isin(_NATIVE_WESTERNS | _NATIVE_CONCENTRATIONS).all():
        raise ValueError("only native measured Hog1, Gpd1 and glycerol observables may be fitted")
    provenance = np.where(training.observable_id.isin(_NATIVE_WESTERNS),
                          "literature_scaled_western_blot", "processed_hplc_measurement")
    if not training.measurement_provenance.eq(provenance).fillna(False).all():
        raise ValueError("native observation provenance must distinguish literature-scaled Westerns from processed HPLC")
    if not all(isinstance(value, str) and value for value in training.experiment_id):
        raise ValueError("explicit native experiment identifiers are required")
    frame = training.reset_index(drop=True).copy()
    numeric = ["time_model_s", "time_relative_s", "value", "nacl_molar"]
    try:
        frame[numeric] = frame[numeric].astype(float)
        sd = frame.sd.to_numpy(dtype=float, na_value=np.nan)
    except (TypeError, ValueError) as exc:
        raise ValueError("native times, observations and SD must be numeric") from exc
    if not np.isfinite(frame[numeric].to_numpy()).all() or (frame.time_model_s < 0).any():
        raise ValueError("native observations and times must be finite and model times nonnegative")
    if np.isinf(sd).any() or (sd < 0).any():
        raise ValueError("reported SD must be nonnegative finite values or missing")
    return frame


def _native_parameter_groups(parameter_groups, source):
    supplied = dict(_NATIVE_GROUPS) if parameter_groups is None else parameter_groups
    if not isinstance(supplied, Mapping) or not supplied:
        raise ValueError("nonempty named parameter groups are required")
    groups = {}
    for name, members in supplied.items():
        if not isinstance(name, str) or not name or isinstance(members, str):
            raise ValueError("parameter groups require names and sequences of source parameters")
        try:
            groups[name] = tuple(members)
        except TypeError as exc:
            raise ValueError("parameter groups require sequences of source parameters") from exc
        if not groups[name]:
            raise ValueError("parameter groups must not be empty")
    names = [name for members in groups.values() for name in members]
    if any(not isinstance(name, str) or name not in source or not np.isfinite(source[name])
           or source[name] <= 0 for name in names):
        raise ValueError("parameter groups must identify positive finite source parameters")
    if len(set(names)) != len(names):
        raise ValueError("parameter groups must contain distinct, nonoverlapping source parameters")
    if set(names) & {"kv16r_1", "parameter_97", "t_stress"}:
        raise ValueError("deactivation remains fixed and experimental inputs are not fitted parameters")
    synthesis = {"kv6_1", "kv6b_4"}
    if synthesis.intersection(names) and not any(synthesis <= set(members) for members in groups.values()):
        raise ValueError("glycerol synthesis requires a shared kv6_1/kv6b_4 multiplier")
    return groups


def _sensitivity_column(column, others):
    norm = float(np.linalg.norm(column))
    independent = column - others @ np.linalg.lstsq(others, column, rcond=1e-3)[0] if others.shape[1] else column
    independent_norm = float(np.linalg.norm(independent))
    return {"column_norm": norm, "orthogonal_column_norm": independent_norm,
            "orthogonal_fraction": independent_norm / norm if norm > 0 else 0.0}


def fit_native_hog_parameters(
    model: HogModel, training: pd.DataFrame, *, parameter_groups=None,
    fold_bounds=(0.2, 5.0), observable_weights=None, max_nfev=30,
) -> NativeHogFit:
    frame = _native_training(training)
    source = model.parameter_values
    groups = _native_parameter_groups(parameter_groups, source)
    bounds = np.asarray(fold_bounds, dtype=float)
    if bounds.shape != (2,) or not np.isfinite(bounds).all() or not 0 < bounds[0] < 1 < bounds[1]:
        raise ValueError("fold_bounds must bracket the positive source reference")
    try:
        budget = int(max_nfev)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("max_nfev must be a positive integer") from exc
    if isinstance(max_nfev, (bool, np.bool_)) or budget != max_nfev or budget < 1:
        raise ValueError("max_nfev must be a positive integer")
    blocks = {name: group.index.to_numpy() for name, group in frame.groupby("observable_id", sort=False)}
    requested_weights = {} if observable_weights is None else observable_weights
    if not isinstance(requested_weights, Mapping) or set(requested_weights) - set(blocks):
        raise ValueError("observable_weights must name only supplied native observables")
    weights = {name: requested_weights.get(name, 1.0) for name in blocks}
    if any(not isinstance(value, Real) or isinstance(value, (bool, np.bool_))
           or not np.isfinite(value) or value <= 0 for value in weights.values()):
        raise ValueError("observable_weights must be finite and positive")
    total_weight = float(sum(weights.values()))
    if not np.isfinite(total_weight):
        raise ValueError("observable_weights must have a finite sum")
    weights = {name: float(value / total_weight) for name, value in weights.items()}
    values = frame.value.to_numpy(dtype=float)
    times = frame.time_model_s.to_numpy(dtype=float)
    observable_ids = frame.observable_id.to_numpy()
    scales, residual_weights = {}, np.empty(len(frame))
    for name, indices in blocks.items():
        scales[name] = float(np.sqrt(np.mean(values[indices] ** 2)))
        if not np.isfinite(scales[name]) or scales[name] <= 0:
            raise ValueError("each native observable needs finite nonzero training excitation")
        residual_weights[indices] = np.sqrt(weights[name] / len(indices)) / scales[name]
    gain_count = len(set(blocks) & _NATIVE_WESTERNS)
    coordinates = len(frame.drop_duplicates(["experiment_id", "observable_id", "time_model_s"]))
    if coordinates <= len(groups) + gain_count or len(np.unique(times)) < 2:
        raise ValueError("insufficient measured times for the effective native parameter count")
    protocol_rows = {}
    for _, experiment in frame.groupby("experiment_id", sort=False):
        protocol_rows.setdefault(_protocol(experiment), []).extend(experiment.index.tolist())
    protocol_rows = {protocol: np.asarray(indices, dtype=int) for protocol, indices in protocol_rows.items()}
    cache, failures, difference_steps = {}, [], {}
    simulations = 0
    penalty = np.full(len(frame), 1e6)
    lower, upper = np.log(bounds)
    step = 1e-3

    def evaluate(log_folds, extra=None):
        nonlocal simulations
        key = (tuple(float(value) for value in log_folds), tuple(sorted((extra or {}).items())))
        if key not in cache:
            parameters = {name: float(source[name] * np.exp(fold))
                          for members, fold in zip(groups.values(), key[0]) for name in members}
            parameters.update(extra or {})
            try:
                predicted = np.empty(len(frame))
                for protocol, indices in protocol_rows.items():
                    grid = np.unique(np.r_[0.0, times[indices]])
                    simulations += 1
                    result = model.simulate(grid, protocol=protocol, parameters=parameters, rtol=1e-8, atol=1e-11)
                    for name in np.unique(observable_ids[indices]):
                        series = np.asarray(result.variables.get(name), dtype=float)
                        if series.shape != grid.shape or not np.isfinite(series).all():
                            raise KineticSimulationError(f"missing, malformed or nonfinite native observable {name}")
                        selected = indices[observable_ids[indices] == name]
                        predicted[selected] = series[np.searchsorted(grid, times[selected])]
                gains = {}
                for name in blocks.keys() & _NATIVE_WESTERNS:
                    indices = blocks[name]
                    response = predicted[indices]
                    denominator = float(response @ response)
                    gain = float(response @ values[indices] / denominator) if denominator > 0 else float("nan")
                    if not np.isfinite(gain) or gain <= 0:
                        raise KineticSimulationError(f"native Western {name} does not support a positive observation gain")
                    gains[name] = gain
                    predicted[indices] *= gain
                residual = (predicted - values) * residual_weights
                if not np.isfinite(residual).all():
                    raise KineticSimulationError("nonfinite normalized native residual")
                cache[key] = (residual, gains, parameters, None)
            except KineticSimulationError as exc:
                failures.append({"parameters": parameters.copy(), "message": str(exc)})
                cache[key] = (penalty.copy(), {}, parameters, str(exc))
        return cache[key]

    def difference(log_folds, index):
        base = evaluate(log_folds)
        if base[3] is not None:
            raise RuntimeError(f"native sensitivity base cannot simulate: {base[3]}")
        deltas = (min(step, upper - log_folds[index]), -min(step, log_folds[index] - lower))
        for delta in sorted(deltas, key=abs, reverse=True):
            if abs(delta) <= 1e-12:
                continue
            shifted = np.array(log_folds, dtype=float, copy=True)
            shifted[index] += delta
            result = evaluate(shifted)
            if result[3] is None:
                difference_steps[list(groups)[index]] = float(delta)
                return (result[0] - base[0]) / delta
        raise RuntimeError("native local sensitivity cannot be estimated with valid in-bound simulations")

    def jacobian(log_folds):
        return np.column_stack([difference(log_folds, index) for index in range(len(groups))])

    baseline = evaluate(np.zeros(len(groups)))
    if baseline[3] is not None:
        raise RuntimeError(f"native source reference cannot simulate the fitting protocols: {baseline[3]}")
    penalty = np.full(len(frame), 1e3 * (1 + np.linalg.norm(baseline[0])) / np.sqrt(len(frame)))
    optimized = least_squares(
        lambda log_folds: evaluate(log_folds)[0], np.zeros(len(groups)), jac=jacobian,
        bounds=(lower, upper), max_nfev=budget, ftol=1e-6, xtol=1e-6, gtol=1e-6)
    residuals, gains, parameters, failure = evaluate(optimized.x)
    if failure is not None:
        raise RuntimeError(f"selected native parameter fit is numerically invalid: {failure}")
    sensitivity = jacobian(optimized.x)
    singular = np.linalg.svd(sensitivity, compute_uv=False)
    numerical_tolerance = max(sensitivity.shape) * np.finfo(float).eps * singular[0]
    practical_tolerance = max(1e-6, 1e-3 * singular[0])
    columns = {name: _sensitivity_column(sensitivity[:, index], np.delete(sensitivity, index, axis=1))
               for index, name in enumerate(groups)}
    weak = [name for name, column in columns.items() if column["orthogonal_column_norm"] <= practical_tolerance]
    fixed_sensitivities = {}
    if "kv13b_1" in source and "kv13b_1" not in parameters and source["kv13b_1"] > 0:
        fixed_sensitivities["kv13b_1"] = {"estimated": False, "available": False}
        for delta in (min(step, upper), max(-step, lower)):
            probe = evaluate(optimized.x, {"kv13b_1": float(source["kv13b_1"] * np.exp(delta))})
            if probe[3] is None:
                fixed_sensitivities["kv13b_1"].update(
                    _sensitivity_column((probe[0] - residuals) / delta, sensitivity),
                    available=True, log_difference_step=float(delta))
                break
    errors = {}
    for name, indices in blocks.items():
        nrmse = float(np.sqrt(residuals[indices] @ residuals[indices] / weights[name]))
        baseline_nrmse = float(np.sqrt(baseline[0][indices] @ baseline[0][indices] / weights[name]))
        errors[name] = {"rows": len(indices), "rmse": nrmse * scales[name], "nrmse": nrmse,
                        "baseline_rmse": baseline_nrmse * scales[name], "baseline_nrmse": baseline_nrmse,
                        "unit": "mol/L", "weight": weights[name], "training_scale": scales[name],
                        "reported_sd_rows": int(frame.sd.iloc[indices].notna().sum())}
    proximity = {name: float(min(value - lower, upper - value) / (upper - lower))
                 for name, value in zip(groups, optimized.x)}
    warnings = []
    if not optimized.success:
        warnings.append("bounded optimization did not converge; returned a numerically valid iterate, not a validated fit")
    if failures:
        warnings.append("failed trial simulations were rejected and recorded")
    if weak:
        warnings.append("weak or confounded local sensitivity: prefer a smaller parameterization")
    if any(value <= 0.05 for value in proximity.values()):
        warnings.append("one or more fitted multipliers are near the search bounds")
    if any(not result["available"] for result in fixed_sensitivities.values()):
        warnings.append("fixed-uptake sensitivity probe failed; no uptake identifiability claim is supported")
    diagnostics = {
        "source": "explicitly supplied s002 wild-type native measured training observations only",
        "parameter_groups": {name: list(members) for name, members in groups.items()},
        "published_parameters": {name: float(source[name]) for name in parameters},
        "parameter_multipliers": {name: float(np.exp(value)) for name, value in zip(groups, optimized.x)},
        "kinetic_parameter_count": len(groups), "profiled_gain_count": gain_count,
        "effective_parameter_count": len(groups) + gain_count,
        "absolute_parameter_override_count": len(parameters),
        "observation_gains": gains, "baseline_observation_gains": baseline[1],
        "fixed_observation_gains": {name: 1.0 for name in blocks.keys() & _NATIVE_CONCENTRATIONS},
        "gain_source": "one positive least-squares gain per literature-scaled Western observable, profiled on training only",
        "observation_transforms": "source measurement variables used directly; no additional volume correction or dilution",
        "native_model_units_relabelled": False, "concentration_observation_gains_fitted": False,
        "normalization_scales": scales, "observable_weights": weights,
        "scale_source": "per-observable RMS of supplied training values; weights normalized across observables, not rows",
        "sd_treatment": "reported SD retained for diagnostics, not inverse-variance weighted",
        "sd_weighting_reason": "comparable SD and independent replicate counts are not available for all native assays",
        "per_observable_errors": errors,
        "baseline_training_nrmse": float(np.linalg.norm(baseline[0])),
        "training_nrmse": float(np.linalg.norm(residuals)),
        "training_experiments": list(frame.experiment_id.unique()), "training_rows": len(frame),
        "distinct_observation_coordinates": coordinates,
        "protocols": [{"nacl_molar": protocol.nacl_molar, "shock_time_s": protocol.shock_time_s,
                       "training_rows": len(indices)} for protocol, indices in protocol_rows.items()],
        "initial_state_policy": "unchanged source initial state at model time zero for each protocol",
        "optimization_converged": bool(optimized.success), "optimization_message": str(optimized.message),
        "optimization_nfev": int(optimized.nfev), "max_nfev": budget, "search_start_count": 1,
        "search_fold_bounds": bounds.tolist(), "active_bounds": optimized.active_mask.tolist(),
        "bound_proximity_fraction": proximity, "near_bound_threshold_fraction": 0.05,
        "objective_evaluations": len(cache), "simulation_evaluations": simulations,
        "failed_simulations": len(failures), "simulation_failures": failures,
        "local_sensitivity_rank": int(np.count_nonzero(singular > practical_tolerance)),
        "local_sensitivity_numerical_rank": int(np.count_nonzero(singular > numerical_tolerance)),
        "local_sensitivity_singular_values": singular.tolist(),
        "local_sensitivity_condition_number": float(singular[0] / singular[-1]) if singular[-1] > 0 else None,
        "local_sensitivity_rank_relative_tolerance": 1e-3, "local_sensitivity_rank_absolute_tolerance": 1e-6,
        "local_sensitivity_columns": columns, "local_sensitivity_steps": difference_steps.copy(),
        "weak_parameter_groups": weak, "fixed_parameter_sensitivities": fixed_sensitivities,
        "log_parameter_difference_step": step,
        "identifiability_scope": "local sensitivity conditional on remaining source priors, with Western gains reprofiled; not global identifiability or confidence intervals",
        "active_uptake_estimated": "kv13b_1" in parameters,
        "active_uptake_policy": "fixed prior by default; a local sensitivity probe alone does not establish independent native-data support",
        "fixed_deactivation_parameter": source.get("kv16r_1"),
        "remaining_parameters": "unchanged source priors; no independent Gpd1/Gpd2 contribution estimate from a shared synthesis multiplier",
        "parameters_fitted_to_measurements": True, "mechanistic_structure_learned": False,
        "warnings": warnings,
    }
    return NativeHogFit(tuple(parameters.items()), tuple(sorted(gains.items())), diagnostics["training_nrmse"], diagnostics)
