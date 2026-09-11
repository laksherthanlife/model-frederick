from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

import numpy as np


TIME_UNITS = {"s": 1.0, "min": 60.0, "h": 3600.0}
DOSE_UNITS = {"M": ("M", 1.0), "mM": ("M", 0.001), "uM": ("M", 0.000001), "g/L": ("g/L", 1.0), "%w/v": ("g/L", 10.0)}
MAX_CANDIDATES = 128
MAX_EXPERIMENTS = 64
MAX_AGGREGATE_TIMES = 10000
MAX_EXPORT_BYTES = 32 * 1024 * 1024
VERIFIER_GAPS = [
    "native_projection_adapter_and_independently_registered_candidate_admission",
    "authenticated_closed_ancestral_lineage_and_current_registry_head",
    "actual_crossed_input_and_absolute_acquisition_clock_strata",
    "verified_source_template_and_input_reassignment_comparators",
    "measurement_uncertainties_and_calibrated_prediction_intervals",
    "independent_final_sham_and_nuclear_marker_observations",
]
EQUATIONS = {
    "observation_calibration": "y_standardized = (measured_proxy - train_center) / train_scale; inverse calibration preserves the declared proxy unit",
    "constant": "y_standardized = intercept",
    "time_only": "y_standardized = intercept + sum(b_k * elapsed_seconds**k); nonmechanistic; elapsed is relative to physical history start",
    "instantaneous_input": "y_standardized = intercept + b_u * u(t)",
    "causal_state": "dz_tau/dt = (u - z_tau) / tau; z_tau(start) = initial_dose; y_standardized = intercept + b_u*u(t) + sum(b_tau*z_tau(t)); terms optional, coefficient signs free",
    "ridge": "minimize sum(equal_experiment_observed_cell_time_weight * residual**2) + ridge * sum(train_standardized_feature_coefficient**2); intercept unpenalized",
    "selection": "minimum equal-experiment development MAE over observed cell-time values; ties within 1e-12 prefer fewer free parameters then model_id; no development refit",
}
LIMITATIONS = [
    "This is an empirical response-family selection, not a universal biological mechanism or a learned biochemical reaction network.",
    "Response states are causal input filters, not measured or teacher-labeled latent biochemical states.",
    "The readout is a nuclear-localization proxy, not phosphorylation, absolute enzyme activity, metabolic flux, or a product outcome.",
    "Support checks are marginal physical-history and filter envelopes, not a certificate of joint support or causal identifiability.",
    "Unknown prehistory, multiple simultaneous input channels, and unpaired mixed readouts are not supported.",
    "Cells in one experiment require a common declared time axis; no cell-time resampling or alignment is inferred.",
    "Train/development comparison is not an independent final test; repeated development use requires verifier governance.",
    "Source and normalization checks are externally verified attestations pinned by SHA256; this learner never audits raw mixed files.",
]


class ContractError(ValueError):
    pass


class FitFailure(ContractError):
    def __init__(self, message, attempts):
        super().__init__(message)
        self.attempts = json.loads(_json_bytes(attempts))


def _require(condition, message):
    if not condition:
        raise ContractError(message)


def _number(value, label):
    _require(not isinstance(value, (bool, np.bool_)) and isinstance(value, (int, float, np.integer, np.floating)), f"{label} must be numeric")
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    return result


def _text(value, label):
    _require(isinstance(value, str) and bool(value.strip()), f"{label} must be explicit")
    return value


def _json_bytes(payload):
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _digest(payload):
    return hashlib.sha256(_json_bytes(payload)).hexdigest()


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


@dataclass(frozen=True)
class Readout:
    observable_id: str
    context_id: str
    unit: str
    quantity_kind: str
    normalization: str
    normalization_reference: str

    def __post_init__(self):
        _text(self.observable_id, "readout observable_id")
        _text(self.context_id, "readout assay context_id")
        _require(self.unit == "dimensionless", "readout unit must remain dimensionless localization proxy")
        _require(self.quantity_kind == "nuclear_localization_proxy", "readout quantity must be nuclear_localization_proxy, not another biological target")
        _require(self.normalization in {"raw_max5_over_median", "audited_fixed_published", "audited_train_only"}, "readout normalization must be explicit and audited")
        _text(self.normalization_reference, "readout normalization evidence")


@dataclass(frozen=True)
class PhysicalHistory:
    quantity: str
    dose_unit: str
    time_unit: str
    start_time: float
    initial_dose: float
    changes: tuple[tuple[float, float], ...]
    interpolation: str = "previous"
    prehistory: str = "constant_equilibrium"

    def __post_init__(self):
        _text(self.quantity, "physical input quantity")
        _require(self.dose_unit in DOSE_UNITS and self.time_unit in TIME_UNITS, "unsupported or ambiguous physical input unit")
        _require(self.interpolation == "previous", "only explicit piecewise-constant physical inputs are supported")
        _require(self.prehistory == "constant_equilibrium", "prehistory must explicitly establish constant equilibrium")
        start = _number(self.start_time, "history start")
        initial = _number(self.initial_dose, "initial dose")
        _require(initial >= 0.0, "concentration cannot be negative")
        events = []
        for pair in self.changes:
            _require(len(pair) == 2, "each physical input change requires time and dose")
            time, dose = _number(pair[0], "event time"), _number(pair[1], "event dose")
            _require(time >= start and (not events or time > events[-1][0]), "physical input events must be strictly ordered after history start")
            _require(dose >= 0.0, "concentration cannot be negative")
            events.append((time, dose))
        object.__setattr__(self, "changes", tuple(events))


def _physical(history, times):
    _require(isinstance(history, PhysicalHistory), "physical history must be explicit")
    time_scale = TIME_UNITS[history.time_unit]
    unit, dose_scale = DOSE_UNITS[history.dose_unit]
    values = np.asarray([_number(t, "observation time") for t in times], dtype=float)
    _require(values.ndim == 1 and values.size > 0, "observation times must be a nonempty vector")
    with np.errstate(over="ignore", invalid="ignore"):
        elapsed = (values - history.start_time) * time_scale
    events = [((t - history.start_time) * time_scale, u * dose_scale) for t, u in history.changes]
    initial = history.initial_dose * dose_scale
    _require(np.isfinite(elapsed).all() and math.isfinite(initial) and all(math.isfinite(t) and math.isfinite(u) for t, u in events), "nonfinite canonical physical history")
    _require(np.all(elapsed >= 0.0) and np.all(np.diff(elapsed) > 0.0), "observation times must increase and cannot precede verified history")
    return elapsed, initial, events, unit


def causal_response_basis(history, times, tau_seconds):
    elapsed, initial, events, _ = _physical(history, times)
    taus = np.asarray([_number(t, "tau") for t in tau_seconds], dtype=float)
    _require(np.all(taus > 0.0) and len(taus) <= 6, "tau basis must be positive and bounded")
    states = np.full(len(taus), initial, dtype=float)
    result = np.empty((len(elapsed), len(taus) + 1), dtype=float)
    dose, previous, index = initial, 0.0, 0
    for row, time in enumerate(elapsed):
        while index < len(events) and events[index][0] <= time:
            event_time, next_dose = events[index]
            states += (dose - states) * (-np.expm1(-(event_time - previous) / taus))
            previous, dose = event_time, next_dose
            index += 1
        states += (dose - states) * (-np.expm1(-(time - previous) / taus))
        previous = time
        result[row] = (dose, *states)
    _require(np.isfinite(result).all(), "nonfinite physical response basis")
    return result


@dataclass(frozen=True)
class NativeExperiment:
    experiment_id: str
    biological_replicate_id: str
    source_groups: tuple[str, ...]
    partition: str
    history: PhysicalHistory
    times: tuple[float, ...]
    values: tuple[tuple[float | None, ...], ...]
    readout: Readout
    mask: tuple[tuple[bool, ...], ...] | None = None

    def __post_init__(self):
        _text(self.experiment_id, "experiment_id")
        _text(self.biological_replicate_id, "biological_replicate_id")
        _require(self.partition in {"train", "development"}, "only train/development partitions are permitted")
        _require(isinstance(self.source_groups, (list, tuple)), "source groups must be a structured sequence")
        groups = tuple(_text(g, "source group") for g in self.source_groups)
        _require(groups and len(set(groups)) == len(groups), "source groups must be nonempty and unique")
        _require(isinstance(self.readout, Readout), "readout contract is required")
        object.__setattr__(self, "source_groups", groups)
        object.__setattr__(self, "times", tuple(self.times))
        object.__setattr__(self, "values", tuple(tuple(row) for row in self.values))
        if self.mask is not None:
            object.__setattr__(self, "mask", tuple(tuple(row) for row in self.mask))
        _physical(self.history, self.times)
        aggregate_observations(self)


def _observed_cells(experiment):
    n_times = len(experiment.times)
    _require(bool(experiment.values) and all(len(row) == n_times for row in experiment.values), "cell-by-time values must match observation times")
    mask = experiment.mask
    if mask is not None:
        _require(len(mask) == len(experiment.values) and all(len(row) == n_times for row in mask), "mask must match the cell-by-time values")
        _require(all(isinstance(v, (bool, np.bool_)) for row in mask for v in row), "mask entries must be boolean")
    values, indices = [], []
    for cell, row in enumerate(experiment.values):
        for time, value in enumerate(row):
            observed = mask[cell][time] if mask is not None else value is not None
            if observed:
                _require(value is not None, "observed values must be finite, not null")
                value = _number(value, "observed value")
                if experiment.readout.normalization == "raw_max5_over_median":
                    _require(value >= 0.0, "raw localization observations must be nonnegative")
                values.append(value)
                indices.append(time)
    _require(bool(values), "experiment has no observed values")
    return np.asarray(values), np.asarray(indices, dtype=int)


def _aggregate_cells(cells, n_times):
    values, indices = cells
    counts = np.bincount(indices, minlength=n_times)
    sums = np.bincount(indices, weights=values, minlength=n_times)
    means = [float(total / count) if count else None for total, count in zip(sums, counts)]
    _require(all(v is None or math.isfinite(v) for v in means), "observed means must be finite")
    return {"mean": means, "count": counts.tolist(), "mask": (counts > 0).tolist()}


def aggregate_observations(experiment):
    return _aggregate_cells(_observed_cells(experiment), len(experiment.times))


@dataclass(frozen=True)
class SelectionRecipe:
    tau_seconds: tuple[float, ...] = (30.0, 120.0, 600.0)
    ridge: tuple[float, ...] = (0.0, 0.01, 0.1)
    max_state_terms: int = 2
    time_only_degrees: tuple[int, ...] = (1, 2)
    minimum_baseline_gain: float = 0.15
    minimum_absolute_scaled_gain: float = 0.05
    support_policy: str = "flag"

    def __post_init__(self):
        taus = tuple(_number(t, "tau") for t in self.tau_seconds)
        penalties = tuple(_number(r, "ridge") for r in self.ridge)
        degrees = tuple(self.time_only_degrees)
        _require(taus and len(taus) <= 6 and all(t > 0.0 for t in taus) and len(set(taus)) == len(taus), "tau candidate grid must be positive, unique and bounded")
        _require(penalties and len(penalties) <= 6 and all(r >= 0.0 for r in penalties) and len(set(penalties)) == len(penalties), "ridge candidate grid must be nonnegative, unique and bounded")
        _require(type(self.max_state_terms) is int and 1 <= self.max_state_terms <= min(2, len(taus)), "state subset search must be bounded to one or two states")
        _require(degrees and all(type(d) is int and d in (1, 2) for d in degrees) and len(set(degrees)) == len(degrees), "time-only degree grid is bounded to degrees one and two")
        count = 1 + len(degrees) * len(penalties) + len(penalties) * (1 + 2 * sum(math.comb(len(taus), k) for k in range(1, self.max_state_terms + 1)))
        _require(count <= MAX_CANDIDATES, "candidate search exceeds the bounded budget")
        _require(_number(self.minimum_absolute_scaled_gain, "absolute scaled gain") >= 0.0, "absolute scaled gain cannot be negative")
        _require(0.0 <= _number(self.minimum_baseline_gain, "baseline gain") < 1.0, "baseline gain must be in [0, 1)")
        _require(self.support_policy in {"flag", "refuse"}, "support policy must flag or refuse, never clip")
        object.__setattr__(self, "tau_seconds", taus)
        object.__setattr__(self, "ridge", penalties)
        object.__setattr__(self, "time_only_degrees", degrees)


def _candidates(recipe):
    forms = [("constant", [])]
    forms.extend(("time_only", [{"kind": "elapsed", "degree": k} for k in range(1, degree + 1)]) for degree in recipe.time_only_degrees)
    forms.append(("instantaneous_input", [{"kind": "input"}]))
    for size in range(1, recipe.max_state_terms + 1):
        for subset in combinations(recipe.tau_seconds, size):
            states = [{"kind": "state", "tau_seconds": tau} for tau in subset]
            forms.extend((("causal_state", states), ("causal_state", [{"kind": "input"}, *states])))
    for family, terms in forms:
        for ridge in (0.0,) if family == "constant" else recipe.ridge:
            yield {
                "id": f"{family}:{_digest([terms, ridge])[:16]}",
                "family": family,
                "terms": terms,
                "ridge": ridge,
                "free_parameter_count": 1 + len(terms) + sum(term["kind"] == "state" for term in terms),
                "interpretation": "nonmechanistic_time_only_comparator" if family == "time_only" else "empirical_response_basis_not_a_biochemical_mechanism",
            }


def candidate_library(recipe):
    _require(isinstance(recipe, SelectionRecipe), "candidate library requires a bounded recipe")
    return json.loads(_json_bytes(list(_candidates(recipe))))


def _validate_candidate(candidate):
    family, terms = candidate["family"], candidate["terms"]
    _require(isinstance(terms, list) and len(terms) <= 3, "candidate terms must be a bounded list")
    kinds = [term["kind"] for term in terms]
    if family == "constant":
        _require(not terms, "constant candidate cannot contain input or clock terms")
    elif family == "time_only":
        _require(terms in ([{"kind": "elapsed", "degree": 1}], [{"kind": "elapsed", "degree": 1}, {"kind": "elapsed", "degree": 2}]), "invalid time-only candidate")
    elif family == "instantaneous_input":
        _require(terms == [{"kind": "input"}], "invalid instantaneous-input candidate")
    else:
        _require(family == "causal_state" and kinds.count("state") in (1, 2) and kinds.count("input") <= 1 and all(kind in {"state", "input"} for kind in kinds), "causal candidate cannot contain clock or undeclared terms")
        taus = []
        for term in terms:
            if term["kind"] == "state":
                _require(set(term) == {"kind", "tau_seconds"} and _number(term["tau_seconds"], "candidate tau") > 0.0, "candidate state requires a positive physical timescale")
                taus.append(term["tau_seconds"])
            else:
                _require(term == {"kind": "input"}, "unexpected candidate input fields")
        _require(len(taus) == len(set(taus)), "duplicate candidate response states")
    _require(_number(candidate["ridge"], "candidate ridge") >= 0.0, "negative candidate ridge")
    expected_interpretation = "nonmechanistic_time_only_comparator" if family == "time_only" else "empirical_response_basis_not_a_biochemical_mechanism"
    _require(candidate["interpretation"] == expected_interpretation, "candidate interpretation cannot assert a mechanism")
    _require(candidate["free_parameter_count"] == 1 + len(terms) + kinds.count("state"), "candidate parameter count mismatch")
    _require(candidate["id"] == f"{family}:{_digest([terms, candidate['ridge']])[:16]}", "candidate id does not bind its equation")


def _design(history, times, terms):
    elapsed, _, _, _ = _physical(history, times)
    taus = tuple(term["tau_seconds"] for term in terms if term["kind"] == "state")
    basis = causal_response_basis(history, times, taus)
    columns, state = [], 1
    for term in terms:
        if term["kind"] == "input":
            columns.append(basis[:, 0])
        elif term["kind"] == "state":
            columns.append(basis[:, state])
            state += 1
        elif term["kind"] == "elapsed":
            with np.errstate(over="ignore", invalid="ignore"):
                columns.append(elapsed ** term["degree"])
        else:
            raise ContractError("unknown candidate term")
    result = np.column_stack(columns) if columns else np.empty((len(times), 0))
    _require(np.isfinite(result).all(), "nonfinite candidate features; refusing prediction")
    return result


def _support_features(history, times, taus):
    elapsed, initial, events, _ = _physical(history, times)
    basis = causal_response_basis(history, times, taus)
    names = ["input_level", "initial_dose", "history_minimum", "history_maximum", "elapsed_seconds", "change_count", "total_variation"]
    names.extend(f"state_tau_{tau:g}_seconds" for tau in taus)
    result = np.empty((len(times), len(names)))
    minimum, maximum, previous, count, variation, index = initial, initial, initial, 0, 0.0, 0
    for row, time in enumerate(elapsed):
        while index < len(events) and events[index][0] <= time:
            value = events[index][1]
            minimum, maximum = min(minimum, value), max(maximum, value)
            if value != previous:
                count += 1
                variation += abs(value - previous)
            previous = value
            index += 1
        result[row] = (basis[row, 0], initial, minimum, maximum, time, count, variation, *basis[row, 1:])
    _require(np.isfinite(result).all(), "nonfinite physical support features")
    return names, result


def _validate_lineage(records, allow_reserved=False):
    seen_ids, seen_replicates, seen_groups = set(), set(), set()
    for row in records:
        partition = row["partition"]
        _require(partition in ({"train", "development", "reserved_test"} if allow_reserved else {"train", "development"}), "unapproved partition")
        identity = _text(row["experiment_id"], "experiment_id")
        replicate = _text(row["biological_replicate_id"], "biological_replicate_id")
        _require(isinstance(row["source_groups"], (list, tuple)), "source groups must be a structured sequence")
        groups = tuple(_text(g, "source group") for g in row["source_groups"])
        _require(groups and len(set(groups)) == len(groups), "source groups must be unique and explicit")
        _require(identity not in seen_ids, "duplicate experiment or split leakage")
        _require(replicate not in seen_replicates, "duplicate biological replicate or split leakage; combine coupled observations into one experiment")
        _require(not (set(groups) & seen_groups), "duplicate source group or split leakage")
        seen_ids.add(identity)
        seen_replicates.add(replicate)
        seen_groups.update(groups)


def _record(experiment):
    return {"experiment_id": experiment.experiment_id, "biological_replicate_id": experiment.biological_replicate_id, "source_groups": list(experiment.source_groups), "partition": experiment.partition}


@dataclass(frozen=True)
class NativeResponseModel:
    readout: Readout
    input_quantity: str
    canonical_dose_unit: str
    candidate: dict
    coefficients: tuple[float, ...]
    intercept: float
    observation_calibration: dict
    support: dict
    training_fingerprint: str

    def to_dict(self):
        return json.loads(_json_bytes({"schema_version": 1, "model_kind": "native_empirical_response", **asdict(self)}))

    @classmethod
    def from_dict(cls, payload):
        try:
            raw = json.loads(_json_bytes(payload))
            _require(raw.pop("schema_version") == 1 and raw.pop("model_kind") == "native_empirical_response", "unsupported model schema")
            raw["readout"] = Readout(**raw["readout"])
            _validate_candidate(raw["candidate"])
            _text(raw["input_quantity"], "model input quantity")
            _require(raw["canonical_dose_unit"] in {unit for unit, _ in DOSE_UNITS.values()}, "model input unit must be canonical")
            _require(_sha(raw["training_fingerprint"]), "model needs a training fingerprint")
            raw["coefficients"] = tuple(_number(c, "coefficient") for c in raw["coefficients"])
            _number(raw["intercept"], "intercept")
            _require(len(raw["coefficients"]) == len(raw["candidate"]["terms"]), "coefficient count must match candidate terms")
            calibration = raw["observation_calibration"]
            _require(calibration["fitted_partition"] == "train" and _number(calibration["scale"], "calibration scale") > 0.0, "calibration must be train-only and invertible")
            _number(calibration["center"], "calibration center")
            _require(calibration["unit"] == raw["readout"].unit, "calibration/readout unit mismatch")
            return cls(**raw)
        except (KeyError, TypeError, ValueError) as error:
            raise ContractError(f"invalid serialized model: {error}") from error

    def predict(self, history, times, readout=None, support_policy="flag"):
        _require(readout is None or readout == self.readout, "prediction readout contract mismatch")
        _require(support_policy in {"flag", "refuse"}, "support policy must flag or refuse")
        elapsed, _, _, unit = _physical(history, times)
        _require(history.quantity == self.input_quantity and unit == self.canonical_dose_unit, "prediction input contract/unit mismatch")
        design = _design(history, times, self.candidate["terms"])
        calibration = self.observation_calibration
        with np.errstate(over="ignore", invalid="ignore"):
            prediction = calibration["center"] + calibration["scale"] * (self.intercept + design @ np.asarray(self.coefficients))
        _require(np.isfinite(prediction).all(), "nonfinite prediction; no output clipping is permitted")
        names, features = _support_features(history, times, self.support["tau_seconds"])
        _require(names == self.support["feature_names"], "support feature contract mismatch")
        low, high = np.asarray(self.support["minimum"]), np.asarray(self.support["maximum"])
        tolerance = 1e-10 * np.maximum(1.0, np.maximum(np.abs(low), np.abs(high)))
        outside = (features < low - tolerance) | (features > high + tolerance)
        flags = [[name for name, flag in zip(names, row) if flag] for row in outside]
        domain_flags = [["negative_localization_proxy"] if value < 0.0 and self.readout.normalization == "raw_max5_over_median" else [] for value in prediction]
        values = [None if domain or (support_policy == "refuse" and flag) else float(value) for value, flag, domain in zip(prediction, flags, domain_flags)]
        any_outside = any(flags)
        return {
            "elapsed_seconds": elapsed.tolist(),
            "values": values,
            "readout": asdict(self.readout),
            "canonical_input_unit": unit,
            "support_flags": flags,
            "domain_flags": domain_flags,
            "support_policy": support_policy,
            "status": "refused_physical_domain" if any(domain_flags) else "partially_refused_out_of_support" if any_outside and support_policy == "refuse" else "flagged_out_of_support" if any_outside else "within_marginal_training_envelopes",
        }


def _training_arrays(train, cells):
    summaries = [_aggregate_cells(c, len(e.times)) for e, c in zip(train, cells)]
    masks = [np.asarray(summary["mask"], dtype=bool) for summary in summaries]
    values = np.concatenate([np.asarray([v for v in summary["mean"] if v is not None]) for summary in summaries])
    weights = np.concatenate([np.asarray(summary["count"])[mask] / (len(train) * len(c[0])) for summary, mask, c in zip(summaries, masks, cells)])
    group_means = [float(np.mean(c[0])) for c in cells]
    center = float(np.mean(group_means))
    variance = float(np.mean([np.mean((c[0] - center) ** 2) for c in cells]))
    scale = math.sqrt(variance) if variance > np.finfo(float).eps * max(1.0, center * center) else 1.0
    calibration = {"center": center, "scale": scale, "unit": train[0].readout.unit, "fitted_partition": "train", "source_groups": sorted(g for e in train for g in e.source_groups), "weighting": "equal_experiment_then_equal_observed_cell_time", "training_variance": variance, "purpose": "invertible_numerical_calibration_not_a_source_observation_transform"}
    quartiles = np.quantile(group_means, (0.25, 0.75), method="linear")
    iqr = float(quartiles[1] - quartiles[0])
    scoring_scale = {"training_group_mean_iqr": iqr if iqr > 0.0 else None, "fitted_partition": "train", "unit": train[0].readout.unit, "status": "usable" if iqr > 0.0 else "unidentifiable_zero_training_group_iqr"}
    fingerprint = _digest([{"lineage": _record(e), "history": asdict(e.history), "times": list(e.times), "readout": asdict(e.readout), "observed_values": c[0].tolist(), "observed_time_indices": c[1].tolist()} for e, c in zip(train, cells)])
    return masks, (values - center) / scale, weights, calibration, scoring_scale, fingerprint


def _fit_coefficients(design, target, weights, ridge):
    if not design.shape[1]:
        return (), float(weights @ target), 0, 1.0
    center = weights @ design
    centered = design - center
    scales = np.sqrt(weights @ (centered * centered))
    if np.any(scales <= np.finfo(float).eps * np.maximum(1.0, np.max(np.abs(design), axis=0))):
        raise ContractError("rank_deficient")
    normalized = centered / scales
    weighted = normalized * np.sqrt(weights[:, None])
    rank = int(np.linalg.matrix_rank(weighted))
    _require(rank == design.shape[1], "rank_deficient")
    condition = float(np.linalg.cond(weighted))
    if ridge:
        beta = np.linalg.solve(weighted.T @ weighted + ridge * np.eye(design.shape[1]), normalized.T @ (weights * target))
    else:
        beta = np.linalg.lstsq(weighted, np.sqrt(weights) * target, rcond=None)[0]
    coefficients = beta / scales
    intercept = float(weights @ target - center @ coefficients)
    _require(np.isfinite(coefficients).all() and math.isfinite(intercept), "nonfinite_fit")
    return tuple(float(c) for c in coefficients), intercept, rank, condition


def _experiment_prediction(fitted, experiment, cells, policy, scoring_iqr):
    prediction = fitted.predict(experiment.history, experiment.times, readout=experiment.readout, support_policy=policy)
    summary = _aggregate_cells(cells, len(experiment.times))
    complete = all(predicted is not None for observed, predicted in zip(summary["mean"], prediction["values"]) if observed is not None)
    mse, mae = None, None
    if complete:
        with np.errstate(over="ignore", invalid="ignore"):
            residual = np.asarray(prediction["values"], dtype=float)[cells[1]] - cells[0]
            mse, mae = float(np.mean(residual ** 2)), float(np.mean(np.abs(residual)))
        _require(math.isfinite(mse) and math.isfinite(mae), "nonfinite observation loss")
    scaled_mae = mae / scoring_iqr if complete and scoring_iqr else None
    _require(scaled_mae is None or math.isfinite(scaled_mae), "nonfinite scaled observation loss")
    return {**_record(experiment), **prediction, "observed_cell_counts": summary["count"], "observation_mask": summary["mask"], "n_observed_times": sum(summary["mask"]), "n_observed_values": len(cells[0]), "mse": mse, "mae": mae, "scaled_mae": scaled_mae, "complete_observed_coverage": complete, "prediction_target": "conditional_population_mean_broadcast_to_each_registered_cell"}


def _select_attempt(attempts):
    best_loss = min(item["row"]["development_mae"] for item in attempts)
    tied = [item for item in attempts if item["row"]["development_mae"] <= best_loss + 1e-12]
    return min(tied, key=lambda item: (item["row"]["candidate"]["free_parameter_count"], item["row"]["candidate_id"]))


def _physical_signature(experiment, observed_mask):
    elapsed, initial, events, unit = _physical(experiment.history, experiment.times)
    last_observed = elapsed[np.flatnonzero(observed_mask)[-1]]
    changes, previous = [], initial
    for time, value in events:
        if time > last_observed:
            break
        if value != previous:
            changes.append((time, value))
        previous = value
    return _digest([experiment.history.quantity, unit, initial, changes]), bool(changes)


def fit_native_response(train, development, recipe=None):
    recipe = SelectionRecipe() if recipe is None else recipe
    _require(isinstance(recipe, SelectionRecipe), "a bounded selection recipe is required")
    train, development = list(train), list(development)
    _require(train and development and len(train) + len(development) <= MAX_EXPERIMENTS, "nonempty, bounded whole-experiment train/development sets are required")
    _require(all(isinstance(e, NativeExperiment) for e in train + development), "native experiment contracts are required")
    _require(all(e.partition == "train" for e in train) and all(e.partition == "development" for e in development), "partition mismatch; final tests are forbidden")
    _validate_lineage([_record(e) for e in train + development])
    _require(sum(len(e.times) for e in train + development) <= MAX_AGGREGATE_TIMES, "aggregate timepoint budget exceeded")
    readout = train[0].readout
    input_contract = (train[0].history.quantity, DOSE_UNITS[train[0].history.dose_unit][0])
    _require(all(e.readout == readout for e in train + development), "mixed readout contracts must be fit separately")
    _require(all((e.history.quantity, DOSE_UNITS[e.history.dose_unit][0]) == input_contract for e in train + development), "mixed physical input contracts must be fit separately")
    cells = [_observed_cells(e) for e in train + development]
    masks, target, weights, calibration, scoring_scale, fingerprint = _training_arrays(train, cells[:len(train)])
    scoring_iqr = scoring_scale["training_group_mean_iqr"]
    support_rows = [_support_features(e.history, e.times, recipe.tau_seconds) for e in train]
    support_array = np.vstack([row[1][mask] for row, mask in zip(support_rows, masks)])
    support = {"tau_seconds": list(recipe.tau_seconds), "feature_names": support_rows[0][0], "minimum": support_array.min(axis=0).tolist(), "maximum": support_array.max(axis=0).tolist(), "interpretation": "marginal_envelopes_not_joint_support"}
    candidates, successes = [], []
    library = candidate_library(recipe)
    for candidate in library:
        row = {"candidate_id": candidate["id"], "family": candidate["family"], "candidate": candidate}
        try:
            design = np.vstack([_design(e.history, e.times, candidate["terms"])[mask] for e, mask in zip(train, masks)])
            coefficients, intercept, rank, condition = _fit_coefficients(design, target, weights, candidate["ridge"])
            row["fitted_parameters"] = {"intercept": intercept, "coefficients": list(coefficients), "training_fingerprint": fingerprint}
            fitted = NativeResponseModel(readout, *input_contract, candidate, coefficients, intercept, calibration, support, fingerprint)
            predictions = [_experiment_prediction(fitted, e, c, recipe.support_policy, scoring_iqr) for e, c in zip(train + development, cells)]
            row["per_experiment_scores"] = [{key: prediction[key] for key in ("experiment_id", "partition", "mae", "mse", "scaled_mae", "n_observed_values", "complete_observed_coverage", "status")} for prediction in predictions]
            _require(all(p["complete_observed_coverage"] for p in predictions), "incomplete_development_coverage")
            train_mse = float(np.mean([p["mse"] for p in predictions[:len(train)]]))
            dev_mse = float(np.mean([p["mse"] for p in predictions[len(train):]]))
            train_mae = float(np.mean([p["mae"] for p in predictions[:len(train)]]))
            dev_mae = float(np.mean([p["mae"] for p in predictions[len(train):]]))
            _require(all(math.isfinite(value) for value in (train_mse, dev_mse, train_mae, dev_mae)), "nonfinite_score")
            row.update(train_mse=train_mse, development_mse=dev_mse, train_mae=train_mae, development_mae=dev_mae, development_scaled_mae=dev_mae / scoring_iqr if scoring_iqr else None, selection_score=dev_mae, design_rank=rank, design_condition=condition)
            successes.append({"row": row, "model": fitted.to_dict(), "predictions": predictions})
        except (ContractError, np.linalg.LinAlgError, FloatingPointError) as error:
            row["failure"] = str(error)
        candidates.append(row)
    if not successes:
        raise FitFailure("no candidates with complete development coverage and identifiable training design", candidates)
    selected = _select_attempt(successes)
    baselines = {}
    for family in ("constant", "time_only", "instantaneous_input"):
        options = [item for item in successes if item["row"]["family"] == family]
        baselines[family] = _select_attempt(options) if options else {"failure": "no_identifiable_complete_coverage_candidate"}
    failures = []
    signatures = [_physical_signature(e, mask) for e, mask in zip(train, masks)]
    if len(train) < 2:
        failures.append("fewer_than_two_training_biological_experiments")
    if len({signature for signature, _ in signatures}) < 2:
        failures.append("single_physical_history_input_vs_time_not_identifiable")
    if not any(excited for _, excited in signatures):
        failures.append("no_temporal_input_excitation")
    if calibration["training_variance"] <= np.finfo(float).eps * max(1.0, calibration["center"] ** 2):
        failures.append("no_training_observation_variation")
    if selected["row"]["design_condition"] > 1e8:
        failures.append("ill_conditioned_selected_response_basis")
    if selected["row"]["family"] in {"constant", "time_only"}:
        failures.append("selected_non_input_baseline")
    if any(flags for row in selected["predictions"][len(train):] for flags, observed in zip(row["support_flags"], row["observation_mask"]) if observed):
        failures.append("development_out_of_support")
    if scoring_iqr is None:
        failures.append("zero_training_group_iqr_no_scaled_scientific_score")
    comparator_scores = [baselines[family]["row"]["development_mae"] for family in ("constant", "time_only") if "row" in baselines[family]]
    improves = scoring_iqr is not None and all(selected["row"]["development_mae"] < score - max(recipe.minimum_baseline_gain * score, recipe.minimum_absolute_scaled_gain * scoring_iqr) for score in comparator_scores)
    if not improves:
        failures.append("no_registered_development_gain_over_non_input_baselines")
    return {
        "schema_version": 1,
        "stage": "training_development_only",
        "data_status": "unverified_or_synthetic",
        "claim": "empirical_response_only_not_a_universal_mechanism",
        "empirical_input_response_supported": False,
        "development_input_response_gain": not failures and selected["row"]["family"] in {"instantaneous_input", "causal_state"},
        "scientific_authorization": "none_no_native_law_claim_is_authorized_by_this_fit",
        "selected_family": selected["row"]["family"],
        "model": selected["model"],
        "selection": selected["row"],
        "selection_recipe": {**asdict(recipe), "refit_on_development": False, "score_weighting": "equal_experiment_then_equal_observed_cell_time", "selection_loss": "equal_experiment_mae", "tie_tolerance": 1e-12, "tie_break": "fewer_free_parameters_then_lexicographic_model_id", "coefficient_sign_constraint": "none", "initial_state_prior": "equilibrium_under_declared_constant_prehistory", "candidate_budget": MAX_CANDIDATES, "optimizer": "deterministic_linear_ridge_no_random_seed"},
        "observation_calibration": calibration,
        "scoring_scale": scoring_scale,
        "uncertainty_status": "not_estimated_no_prediction_intervals_or_scientific_authorization",
        "candidate_library_sha256": _digest(library),
        "candidates": candidates,
        "predictions": selected["predictions"],
        "verifier_gaps": list(VERIFIER_GAPS),
        "baselines": baselines,
        "failures": failures,
        "equations": dict(EQUATIONS),
        "limitations": list(LIMITATIONS),
    }


def _read_json(path, byte_limit):
    path = Path(path)
    _require(path.stat().st_size <= byte_limit, f"file exceeds bounded byte budget: {path.name}")
    raw = path.read_bytes()
    _require(len(raw) <= byte_limit, f"file exceeds bounded byte budget: {path.name}")
    try:
        payload = json.loads(raw, parse_constant=lambda token: (_ for _ in ()).throw(ContractError(f"nonfinite JSON token: {token}")))
    except (ValueError, UnicodeDecodeError) as error:
        raise ContractError(f"invalid JSON contract: {path.name}: {error}") from error
    _require(isinstance(payload, dict), "JSON contract must be an object")
    return payload, hashlib.sha256(raw).hexdigest()


def _metadata_plan(manifest_path, protocol_path):
    manifest_path, protocol_path = Path(manifest_path), Path(protocol_path)
    try:
        protocol, protocol_sha = _read_json(protocol_path, 1024 * 1024)
        manifest, manifest_sha = _read_json(manifest_path, 4 * 1024 * 1024)
        if manifest.get("$schema") == "granados/manifest.schema.json":
            _require(manifest.get("status") != "metadata_schema_verified_ineligible_current_protocol", "custodian marks the native candidate ineligible under the current protocol; fitting, selection, calibration and response exports are forbidden")
            custody = manifest["custody"]
            _require(custody["export_authorized"] is True and bool(custody["training_allowlist"]), "custodian has not authorized any training/development projection; no measurement access")
            raise ContractError("custodian projection-schema adapter and independently registered candidate admission are required; local metadata flags do not authorize fitting")
        _require(protocol["schema_version"] == 1 and protocol["status"] == "frozen_before_learning", "protocol must be frozen before learning")
        _require(protocol["purpose"] == "native_observed_input_response_train_development" and protocol["test_access"] == "forbidden", "protocol must forbid final-test access")
        _require(protocol["split_unit"] == "whole_biological_experiment", "protocol requires whole biological experiment splits")
        recipe = SelectionRecipe(**protocol["selection_recipe"])
        _require(protocol["manifest_sha256"] == manifest_sha, "protocol/manifest digest mismatch")
        _require(manifest["schema_version"] == 1 and manifest["status"] == "verified_train_development_exports", "source manifest is not verified for training/development")
        _require(manifest["observations_origin"] == "experimental_measurements", "only experimentally measured native observables are permitted")
        verification = manifest["verification"]
        for field in ("input_history_and_units", "readout_definition_and_normalization", "whole_experiment_lineage", "reserved_arrays_excluded", "no_teacher_or_kinetic_targets", "positive_raw_denominators_verified", "assay_context_verified"):
            _require(verification[field] is True, f"unverified source requirement: {field}")
        _require(verification["normalization_uses_future_or_development_outcomes"] is False, "source normalization cannot use future or development outcomes")
        _text(verification["evidence"], "source verification evidence")
        _require(verification["normalization_fit_scope"] in {"none", "train_only"}, "normalization fit scope cannot include development or reserved experiments")
        readout = Readout(**manifest["readout"])
        assets = manifest["source_assets"]
        _require(assets and all(_sha(asset["sha256"]) and _text(asset["asset_id"], "source asset id") for asset in assets), "source asset checksums are required; no raw assets are opened")
        records = manifest["experiments"]
        _validate_lineage(records, allow_reserved=True)
        approved_ids = {partition: sorted(r["experiment_id"] for r in records if r["partition"] == partition) for partition in ("train", "development")}
        _require(all(approved_ids.values()), "source registry requires train and development experiments")
        if verification["normalization_fit_scope"] == "train_only":
            _require(readout.normalization == "audited_train_only", "train-fitted source normalization requires its own readout contract")
            _require(sorted(verification["normalization_training_experiment_ids"]) == approved_ids["train"], "normalization training membership must match the registered training set")
        else:
            _require(readout.normalization != "audited_train_only", "train-fitted normalization needs training provenance")
        root = manifest_path.resolve().parent
        paths, export_hashes = {}, {}
        for partition in ("train", "development"):
            entry = manifest["approved_exports"][partition]
            _require(entry["contains_only_partition"] == partition and sorted(entry["experiment_ids"]) == approved_ids[partition], "export partition membership must be approved before opening measurements")
            _require(_sha(entry["sha256"]), "approved export digest is required")
            relative = Path(entry["path"])
            _require(not relative.is_absolute() and ".." not in relative.parts and relative.suffix == ".json", "approved export path must be local, relative JSON")
            target = (root / relative).resolve()
            _require(target.is_relative_to(root) and target not in {manifest_path.resolve(), protocol_path.resolve()}, "export path escapes its approved source directory")
            _require(not any((root.joinpath(*relative.parts[:index])).is_symlink() for index in range(1, len(relative.parts) + 1)), "symlinked measurement exports are not permitted")
            _require(target.is_file() and target.stat().st_size <= MAX_EXPORT_BYTES, "approved measurement export is missing or exceeds the byte budget")
            paths[partition] = target
            export_hashes[partition] = entry["sha256"]
        _require(paths["train"] != paths["development"], "train/development exports must be separate files")
        provenance = {"manifest_sha256": manifest_sha, "protocol_sha256": protocol_sha, "approved_export_sha256": export_hashes, "source_assets": assets, "source_verification": verification, "allowed_partitions": ["train", "development"]}
        return {"manifest": manifest, "recipe": recipe, "readout": readout, "paths": paths, "provenance": provenance}
    except (OSError, KeyError, TypeError, ValueError) as error:
        raise ContractError(f"verified data gate refused: {error}") from error


def check_data_gate(manifest_path, protocol_path):
    try:
        plan = _metadata_plan(manifest_path, protocol_path)
        return {"ready": True, "measurement_arrays_opened": False, "authorization": "structural_projection_checks_only_external_verifier_required", "eligible_partitions": ["train", "development"], "provenance": plan["provenance"], "selection_recipe": asdict(plan["recipe"])}
    except ContractError as error:
        return {"ready": False, "measurement_arrays_opened": False, "failure": str(error)}


def load_verified_train_development(manifest_path, protocol_path):
    plan = _metadata_plan(manifest_path, protocol_path)
    result = {"recipe": plan["recipe"], "provenance": plan["provenance"]}
    registry = {row["experiment_id"]: row for row in plan["manifest"]["experiments"]}
    try:
        for partition in ("train", "development"):
            payload, digest = _read_json(plan["paths"][partition], MAX_EXPORT_BYTES)
            _require(digest == plan["provenance"]["approved_export_sha256"][partition], "approved export digest mismatch")
            _require(set(payload) == {"schema_version", "partition", "experiments"} and payload["schema_version"] == 1 and payload["partition"] == partition, "unexpected export schema or partition")
            rows = []
            for raw in payload["experiments"]:
                parsed = dict(raw)
                parsed["history"] = PhysicalHistory(**parsed["history"])
                parsed["readout"] = Readout(**parsed["readout"])
                experiment = NativeExperiment(**parsed)
                _require(experiment.partition == partition and experiment.readout == plan["readout"], "export readout or partition mismatch")
                _require(_record(experiment) == registry[experiment.experiment_id], "export source lineage does not match approved metadata")
                rows.append(experiment)
            expected = sorted(row["experiment_id"] for row in registry.values() if row["partition"] == partition)
            _require(sorted(e.experiment_id for e in rows) == expected, "export experiment membership mismatch")
            result[partition] = rows
        _validate_lineage([_record(e) for partition in ("train", "development") for e in result[partition]])
        return result
    except (OSError, KeyError, TypeError, ValueError) as error:
        raise ContractError(f"approved training/development export refused: {error}") from error


def _prediction_digest(result):
    return _digest({"selected": result["predictions"], "baselines": {family: row.get("predictions") for family, row in result["baselines"].items()}})


def verify_freeze(envelope):
    try:
        _require(envelope["schema_version"] == 1 and envelope["artifact_kind"] == "native_response_training_development_freeze", "unsupported freeze schema")
        result = envelope["result"]
        _require(result["stage"] == "training_development_only", "freeze cannot include a final test")
        _require(_digest(result) == envelope["result_sha256"], "result digest mismatch")
        _require(_digest(result["model"]) == envelope["model_sha256"], "model digest mismatch")
        _require(_prediction_digest(result) == envelope["predictions_sha256"], "predictions digest mismatch")
        NativeResponseModel.from_dict(result["model"])
        rows = result["predictions"] + [p for baseline in result["baselines"].values() for p in baseline.get("predictions", [])]
        _require(all(row["partition"] in {"train", "development"} for row in rows), "unapproved partition in prediction freeze")
        return True
    except (KeyError, TypeError, ValueError) as error:
        raise ContractError(f"freeze verification failed: {error}") from error


def write_freeze(result, path):
    result = json.loads(_json_bytes(result))
    envelope = {"schema_version": 1, "artifact_kind": "native_response_training_development_freeze", "result": result, "result_sha256": _digest(result), "model_sha256": _digest(result["model"]), "predictions_sha256": _prediction_digest(result)}
    verify_freeze(envelope)
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(envelope, sort_keys=True, indent=2, allow_nan=False) + "\n")
    return envelope
