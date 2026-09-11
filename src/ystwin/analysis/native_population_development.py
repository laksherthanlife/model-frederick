from __future__ import annotations

import hashlib
import json
import math
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares


PROTOCOL_ID = "granados-sfp1-published-cohort-development-01"
ASSET_ID = "1b38366b-7473-4ddf-9cf9-26e625533686"
SOURCE_SHA256 = "1ec5a48955492aade95386d48f15109397a19acb60ef820d3c692ca59d9efc87"
SCIENTIFIC_SCOPE = "retrospective_source_cohort_development"
QUANTITY = "published_prefix_eligible_observed_cell_mean_of_max5_over_median"
PREFIX_OPERATOR = "complete_prefix_raw_ratio_mean_and_mean_frame_sd_o_minus_21_to_o_minus_2"
PARAMETERS = {
    "prefix_persistence": (),
    "training_constant": ("c",),
    "prefix_offset": ("a",),
    "prefix_affine_time": ("a", "b"),
    "first_order_response": ("A", "tau"),
    "transient_response": ("A", "tau_off", "tau_on"),
}
MODEL_IDS = tuple(PARAMETERS)
SEEDS = tuple(range(201801, 201809))
PREFIX_FRAMES = tuple(range(28, 48))
RESPONSE_FRAMES = tuple(range(50, 74))
SETTINGS = {
    "optimizer": "scipy.optimize.least_squares",
    "method": "trf",
    "jacobian": "analytic; derivative of the nonnegative map is zero at and below zero",
    "x_scale": "jac",
    "loss": "linear",
    "ftol": 1e-10,
    "xtol": 1e-10,
    "gtol": 1e-10,
    "max_nfev": 2000,
    "loss_evaluation_accounting": "every residual invocation; analytic Jacobian calls are counted separately",
    "prng": "numpy.random.Generator(numpy.random.PCG64(seed))",
    "starts": 8,
    "seeds": list(SEEDS),
    "reference_time_minutes": 60.0,
    "output_map": "max(0, unbounded_mean) before fitting and scoring",
    "within_family_selection": "finite successful minimum training MSE; numerical ties prefer lexicographic parameters then start index",
    "tie_tolerance": "1e-12 * (1 + abs(best MSE))",
}
LIMITATIONS = [
    "scientific_scope=retrospective_source_cohort_development",
    "The original strong candidates remain ineligible; this protocol never authorizes biological readiness and defines no final test.",
    "All groups have one reported Sfp1/context/glucose schedule. No dose response, causal glucose effect, universal mechanism or intrinsic kinetic parameter is identified.",
    "The historical published cohort and processing may depend on future observations; this is not a certified prospective forecast from untouched acquisitions.",
    "Targets are means of individual raw max5/median ratios among the fixed prefix cohort observed at each frame, not individual trajectories, ratios of mean intensities or phosphorylation fractions.",
    "Finite-cell sampling, within-movie dependence, measurement/calibration error, day/strain effects, attrition and inherited processing-selection bias are distinct; unmeasured components remain unresolved.",
    "Cell sample SD is descriptive, not a biological confidence interval or SD/sqrt(cell count) uncertainty estimate.",
    "Local digests and accountable operational custody are not external registration, encryption, adversarial isolation or a scientific authorization.",
]


class DevelopmentError(ValueError):
    def __init__(self, message, code="dependency_or_partition_violation", groups=(), details=None):
        super().__init__(message)
        self.code = code
        self.groups = list(groups)
        self.details = details


def require(condition, message, code="dependency_or_partition_violation", groups=()):
    if not condition:
        raise DevelopmentError(message, code, groups)


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def digest(value):
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _finite_json(value):
    if isinstance(value, dict):
        for item in value.values():
            _finite_json(item)
    elif isinstance(value, list):
        for item in value:
            _finite_json(item)
    elif isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON numeric value")


def strict_json(content):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    def invalid(token):
        raise DevelopmentError("nonfinite JSON numeric token")

    try:
        result = json.loads(content, object_pairs_hook=pairs, parse_constant=invalid)
        _finite_json(result)
        return result
    except (ValueError, UnicodeDecodeError) as error:
        raise DevelopmentError(f"invalid finite JSON: {error}") from error


def _number(value):
    return not isinstance(value, (bool, np.bool_)) and isinstance(value, (int, float, np.integer, np.floating)) and math.isfinite(float(value))


def _equal(a, b):
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    return a == b


def _schema_validate(value, schema, definitions, location="document"):
    if "$ref" in schema:
        reference = schema["$ref"]
        require(reference.startswith("#/$defs/"), "unsupported schema reference")
        _schema_validate(value, definitions[reference.removeprefix("#/$defs/")], definitions, location)
    if "type" in schema:
        types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        matches = {
            "object": isinstance(value, dict), "array": isinstance(value, list),
            "string": isinstance(value, str), "number": _number(value),
            "integer": type(value) is int, "boolean": type(value) is bool, "null": value is None,
        }
        require(any(matches.get(kind, False) for kind in types), f"{location}: schema type mismatch")
    if "const" in schema:
        require(_equal(value, schema["const"]), f"{location}: schema const mismatch")
    if "enum" in schema:
        require(any(_equal(value, item) for item in schema["enum"]), f"{location}: schema enum mismatch")
    if "oneOf" in schema:
        count = 0
        for choice in schema["oneOf"]:
            try:
                _schema_validate(value, choice, definitions, location)
                count += 1
            except DevelopmentError:
                pass
        require(count == 1, f"{location}: schema oneOf mismatch")
    for choice in schema.get("allOf", []):
        _schema_validate(value, choice, definitions, location)
    if "if" in schema:
        try:
            _schema_validate(value, schema["if"], definitions, location)
            branch = "then"
        except DevelopmentError:
            branch = "else"
        if branch in schema:
            _schema_validate(value, schema[branch], definitions, location)
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        require(set(schema.get("required", ())) <= set(value), f"{location}: required schema fields missing")
        extra = set(value) - set(properties)
        additional = schema.get("additionalProperties", True)
        require(additional is not False or not extra, f"{location}: unexpected schema fields")
        for key, item in value.items():
            if key in properties:
                _schema_validate(item, properties[key], definitions, f"{location}/{key}")
            elif isinstance(additional, dict):
                _schema_validate(item, additional, definitions, f"{location}/{key}")
    if isinstance(value, list):
        require(len(value) >= schema.get("minItems", 0) and len(value) <= schema.get("maxItems", math.inf), f"{location}: schema item count mismatch")
        if schema.get("uniqueItems"):
            require(len({canonical_bytes(item) for item in value}) == len(value), f"{location}: duplicate array entries")
        if "items" in schema:
            for index, item in enumerate(value):
                _schema_validate(item, schema["items"], definitions, f"{location}/{index}")
    if isinstance(value, str):
        require(len(value) >= schema.get("minLength", 0), f"{location}: schema string too short")
        if "pattern" in schema:
            require(re.search(schema["pattern"], value) is not None, f"{location}: schema pattern mismatch")
    if _number(value):
        require(value >= schema.get("minimum", -math.inf) and value <= schema.get("maximum", math.inf), f"{location}: schema numeric bound mismatch")
        if "exclusiveMinimum" in schema:
            require(value > schema["exclusiveMinimum"], f"{location}: schema exclusive lower bound mismatch")


def validate_document(protocol, document, definition):
    definitions = protocol["export_contract"]["schemas"]["$defs"]
    require(definition in definitions, "unknown embedded schema")
    _finite_json(document)
    _schema_validate(document, definitions[definition], definitions)
    return True


def validate_protocol(protocol):
    require(protocol["protocol_id"] == PROTOCOL_ID and protocol["document_type"] == "non_authorizing_operational_development_protocol", "wrong development protocol")
    require(protocol["source"]["asset_id"] == ASSET_ID and protocol["source"]["source_sha256"] == SOURCE_SHA256, "unexpected pinned source")
    require(protocol["estimand"]["quantity"] == QUANTITY and protocol["estimand"]["unit"] == "dimensionless_ratio", "wrong development estimand or unit")
    require(protocol["prefix"]["expected_absolute_columns_1based"] == [28, 47] and protocol["prefix"]["minimum_prefix_cells"] == 2, "unregistered prefix operator")
    require(protocol["response_coordinates"]["expected_absolute_columns_1based"] == [50, 73] and protocol["response_coordinates"]["frame_count"] == 24, "unregistered response window")
    require(protocol["partition"]["source_origin_1based"] == 49 and protocol["partition"]["source_frame_count"] == 97, "source origin is a frame index, not a time offset")
    recipe = protocol["model_recipe"]
    require(tuple(row["model_id"] for row in recipe["canonical_families"]) == MODEL_IDS, "exactly six canonical families are required")
    for row in recipe["canonical_families"]:
        require(tuple(row["parameters"]) == PARAMETERS[row["model_id"]] and row["parameter_count"] == len(row["parameters"]), "canonical parameter contract changed")
    require(recipe["time_unit"] == "minute" and recipe["reference_time_minutes"] == 60.0, "wrong physical time units")
    budget = recipe["budget"]
    require(budget["seeds"] == list(SEEDS) and budget["maximum_starts_per_parameterized_family"] == 8 and budget["maximum_loss_evaluations_per_start"] == 2000 and budget["maximum_parameterized_fit_attempts"] == 40, "unregistered fitting budget or seeds")
    for name in ("tau", "tau_on", "tau_off"):
        require(recipe["parameter_domains"][name] == {"lower": 2.5, "upper": 480.0, "unit": "minute"}, "time-constant bounds changed")
    for name in ("c", "a", "b", "A"):
        expected_lower = "0" if name == "c" else "-10*H_train"
        require(recipe["parameter_domains"][name] == {"lower": expected_lower, "upper": "10*H_train", "unit": "dimensionless_ratio"}, "amplitude bounds changed")
    rows = group_specs(protocol)
    require([row["group_id"] for row in rows] == [f"{ASSET_ID}:/rep{i}" for i in range(1, 6)], "unexpected experiment assignment")
    require([row["source_rows"] for row in rows] == [55, 136, 202, 195, 157], "source row inventory changed")
    require(protocol["context"]["initial_glucose"]["unit"] == "percent_as_reported" and protocol["context"]["final_glucose"]["unit"] == "percent_as_reported", "do not silently convert an unverified glucose-percent convention")
    return protocol


def load_protocol(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink() and path.stat().st_size <= 1024 * 1024, "protocol must be a bounded regular metadata file")
    return validate_protocol(strict_json(path.read_bytes()))


def group_specs(protocol):
    return [{**row, "split": split} for split in ("train", "development") for row in protocol["partition"][split]]


def _frames(window):
    require(window in {"prefix", "response"}, "unexpected source window")
    return PREFIX_FRAMES if window == "prefix" else RESPONSE_FRAMES


def validate_group_records(group, records, windows, times_only=False):
    expected = {(cell, frame) for cell in range(group["source_rows"]) for window in windows for frame in _frames(window)}
    seen = set()
    for record in records:
        require(record["group_id"] == group["group_id"], "cross-group source record", groups=[group["group_id"]])
        cell, frame = record["cell_index"], record["frame_1based"]
        require(type(cell) is int and type(frame) is int and (cell, frame) in expected and (cell, frame) not in seen, "duplicate, missing, renumbered or unexpected source row/frame")
        seen.add((cell, frame))
        parent = group["source_pointer"]
        require(record["time_pointer"] == f"{parent}/general/times/{cell}/{frame - 1}", "source time locator mismatch")
        require(record["relative_time_min"] is None or _number(record["relative_time_min"]), "nonfinite or nonnumeric source time")
        if not times_only:
            require(record["split"] == group["split"] and record["window"] in windows and frame in _frames(record["window"]), "source split/window mismatch")
            for field in ("max5", "median"):
                require(record[f"{field}_pointer"] == f"{parent}/GFP/{field}/{cell}/{frame - 1}", "source intensity locator mismatch")
                require(record[field] is None or _number(record[field]), "nonfinite or nonnumeric source intensity")
    require(seen == expected, "source row/frame inventory is incomplete")
    return True


def _classification(record, prefix):
    missing = [key for key in ("max5", "median", "relative_time_min") if record[key] is None]
    if missing:
        return "source_missing", [f"null_{key}" for key in missing], None
    invalid = []
    if record["max5"] < 0.0:
        invalid.append("negative_max5")
    if record["median"] <= 0.0:
        invalid.append("nonpositive_median")
    if not (record["relative_time_min"] < 0.0 if prefix else record["relative_time_min"] > 0.0):
        invalid.append("nonnegative_prefix_time" if prefix else "nonpositive_response_time")
    if invalid:
        return "invalid", invalid, None
    ratio = record["max5"] / record["median"]
    require(math.isfinite(ratio), "raw ratio is not numerically representable; no amplitude exclusion is permitted", "source_schema_or_alignment")
    return "valid", [], float(ratio)


def derive_prefix(group, cells, projection_sha256):
    records = [row for row in cells if row["window"] == "prefix"]
    validate_group_records(group, records, ("prefix",))
    index = {(row["cell_index"], row["frame_1based"]): row for row in records}
    eligible, ratios, excluded = [], [], []
    for cell in range(group["source_rows"]):
        values, failures = [], []
        for frame in PREFIX_FRAMES:
            status, reasons, value = _classification(index[cell, frame], True)
            if status != "valid":
                failures.append({"frame_1based": frame, "status": status, "reasons": reasons})
            else:
                values.append(value)
        if failures:
            excluded.append({"cell_index": cell, "frames": failures})
        else:
            eligible.append(cell)
            ratios.append(values)
    audit = {"group_id": group["group_id"], "source_rows": group["source_rows"], "eligible_cell_indices": eligible, "excluded_cells": excluded}
    if len(eligible) < 2:
        raise DevelopmentError("prefix cohort has fewer than two eligible cells", "prefix_cohort_unusable", [group["group_id"]], audit)
    matrix = np.asarray(ratios, dtype=float)
    mu = float(np.mean(np.mean(matrix, axis=0)))
    dispersion = float(np.mean(np.std(matrix, axis=0, ddof=1)))
    require(math.isfinite(mu) and math.isfinite(dispersion), "nonfinite prefix summary", "prefix_cohort_unusable", [group["group_id"]])
    return {"group_id": group["group_id"], "origin_1based": 49, "eligible_cell_indices": eligible, "mu": mu, "s": dispersion, "input_projection_sha256": projection_sha256, "operator_id": PREFIX_OPERATOR}, audit


def derive_frames(group, response_cells, response_times, prefix):
    validate_group_records(group, response_times, ("response",), times_only=True)
    if response_cells is not None:
        validate_group_records(group, response_cells, ("response",))
    time_index = {(r["cell_index"], r["frame_1based"]): r for r in response_times}
    intensity_index = None if response_cells is None else {(r["cell_index"], r["frame_1based"]): r for r in response_cells}
    eligible = prefix["eligible_cell_indices"]
    require(prefix["group_id"] == group["group_id"] and len(eligible) >= 2, "invalid fixed prefix cohort")
    frames = []
    for frame in RESPONSE_FRAMES:
        times = [time_index[cell, frame]["relative_time_min"] for cell in eligible]
        valid_times = [value for value in times if value is not None and value > 0.0]
        coordinate = float(np.mean(valid_times)) if valid_times else None
        result = {"group_id": group["group_id"], "frame_1based": frame, "time_min": coordinate}
        if intensity_index is not None:
            values, valid, missing, invalid = [], [], [], []
            for cell in eligible:
                record = intensity_index[cell, frame]
                require(_equal(record["relative_time_min"], time_index[cell, frame]["relative_time_min"]), "response times differ from previously supplied timing metadata")
                status, reasons, value = _classification(record, False)
                if status == "valid":
                    values.append(value)
                    valid.append(cell)
                elif status == "source_missing":
                    missing.append({"cell_index": cell, "reasons": reasons})
                else:
                    invalid.append({"cell_index": cell, "reasons": reasons})
            result.update(observed_mean=float(np.mean(values)) if values else None, prefix_cells=len(eligible), valid_cells=len(values), source_missing_cells=len(missing), invalid_cells=len(invalid), cell_sample_sd=float(np.std(values, ddof=1)) if len(values) >= 2 else None, valid_cell_indices=valid, source_missing_details=missing, invalid_details=invalid, under_resolved_sampling=len(values) < 2)
        frames.append(result)
    return frames


def finalize_prepared(prefixes, training_frames, development_coordinates, protocol):
    specs = group_specs(protocol)
    ids = [row["group_id"] for row in specs]
    require([row["group_id"] for row in prefixes] == ids, "exact five prefix summaries are required")
    for prefix in prefixes:
        validate_document(protocol, prefix, "PrefixSummary")
    for records, expected_ids in ((training_frames, ids[:3]), (development_coordinates, ids[3:])):
        actual = [(row["group_id"], row["frame_1based"]) for row in records]
        expected = [(group, frame) for group in expected_ids for frame in RESPONSE_FRAMES]
        require(actual == expected, "all registered group/frame coordinates must be present in source order")
    require(all(row["observed_mean"] is not None and _number(row["observed_mean"]) and row["observed_mean"] >= 0.0 and row["time_min"] is not None and _number(row["time_min"]) and row["time_min"] > 0.0 for row in training_frames), "complete training means and measured coordinates are required; no subset fit", "timing_or_coverage_unusable", ids[:3])
    training_mu = {row["group_id"]: row["mu"] for row in prefixes[:3]}
    values = np.asarray([row["observed_mean"] for row in training_frames])
    residual = values - np.asarray([training_mu[row["group_id"]] for row in training_frames])
    quartiles = np.quantile(residual, [0.25, 0.75], method="linear")
    return {"prefix_summaries": prefixes, "training_frames": training_frames, "development_coordinates": development_coordinates, "H_train": max(1.0, *training_mu.values(), *values.tolist()), "S_train": float(quartiles[1] - quartiles[0]), "scientific_scope": SCIENTIFIC_SCOPE}


def _unbounded_jacobian(model_id, parameters, mu, time, reference_time):
    require(model_id in PARAMETERS and set(parameters) == set(PARAMETERS[model_id]), "one exact canonical family and parameter set is required")
    require(all(_number(value) for value in parameters.values()), "model parameters must be finite")
    mu, time = np.broadcast_arrays(np.asarray(mu, dtype=float), np.asarray(time, dtype=float))
    require(np.isfinite(mu).all() and np.isfinite(time).all() and np.all(mu >= 0.0) and np.all(time >= 0.0) and reference_time > 0.0, "invalid physical mean-function inputs")
    ones = np.ones(mu.shape)
    if model_id == "prefix_persistence":
        raw, columns = mu, []
    elif model_id == "training_constant":
        raw, columns = parameters["c"] * ones, [ones]
    elif model_id == "prefix_offset":
        raw, columns = mu + parameters["a"], [ones]
    elif model_id == "prefix_affine_time":
        relative = time / reference_time
        raw, columns = mu + parameters["a"] + parameters["b"] * relative, [ones, relative]
    elif model_id == "first_order_response":
        amplitude, tau = parameters["A"], parameters["tau"]
        require(tau > 0.0, "time constants must be positive")
        decay = np.exp(-time / tau)
        rise = -np.expm1(-time / tau)
        raw, columns = mu + amplitude * rise, [rise, -amplitude * decay * time / tau ** 2]
    else:
        amplitude, off, on = parameters["A"], parameters["tau_off"], parameters["tau_on"]
        require(off > 0.0 and on > 0.0, "time constants must be positive")
        decay_off, decay_on = np.exp(-time / off), np.exp(-time / on)
        rise = -np.expm1(-time / on)
        basis = rise * decay_off
        raw, columns = mu + amplitude * basis, [basis, amplitude * basis * time / off ** 2, -amplitude * decay_on * decay_off * time / on ** 2]
    require(np.isfinite(raw).all(), "nonfinite mean-function output", "model_execution_failed")
    jacobian = np.stack(columns, axis=-1) if columns else np.empty((*raw.shape, 0))
    return raw, jacobian * (raw > 0.0)[..., None]


def evaluate_mean(model_id, parameters, mu, time, reference_time=60.0):
    raw, _ = _unbounded_jacobian(model_id, parameters, mu, time, reference_time)
    return raw, np.maximum(0.0, raw), raw <= 0.0


def mean_mse(predicted, observed):
    predicted, observed = np.asarray(predicted, dtype=float), np.asarray(observed, dtype=float)
    require(predicted.shape == observed.shape and predicted.size > 0 and np.isfinite(predicted).all() and np.isfinite(observed).all(), "loss requires complete matching means")
    result = float(np.mean((predicted - observed) ** 2))
    require(math.isfinite(result), "nonfinite mean squared loss", "model_execution_failed")
    return result


def parameter_bounds(model_id, H_train):
    require(model_id in PARAMETERS and _number(H_train) and H_train >= 1.0, "invalid family or training bound scale")
    return {name: [2.5, 480.0] if name in {"tau", "tau_off", "tau_on"} else [0.0 if name == "c" else -10.0 * H_train, 10.0 * H_train] for name in PARAMETERS[model_id]}


def _initial_parameters(model_id, bounds, start, training_level):
    generator = np.random.Generator(np.random.PCG64(SEEDS[start]))
    values = {}
    for name in PARAMETERS[model_id]:
        low, high = bounds[name]
        if start == 0:
            value = math.sqrt(low * high) if name in {"tau", "tau_off", "tau_on"} else training_level if name == "c" else 0.0
        else:
            value = float(np.exp(generator.uniform(np.log(low), np.log(high)))) if name in {"tau", "tau_off", "tau_on"} else float(generator.uniform(low, high))
        values[name] = value
    return values


def _select_start(attempts, names):
    successful = [attempt for attempt in attempts if attempt["status"] == "ok"]
    if not successful:
        return None
    best = min(attempt["training_mean_mse"] for attempt in successful)
    tied = [attempt for attempt in successful if attempt["training_mean_mse"] <= best + 1e-12 * (1.0 + abs(best))]
    return min(tied, key=lambda attempt: (tuple(attempt["parameters"][name] for name in names), attempt["start_index"]))


def fit_family(model_id, mu, time, observed, H_train):
    names = PARAMETERS[model_id]
    bounds = parameter_bounds(model_id, H_train)
    mu, time, observed = np.asarray(mu), np.asarray(time), np.asarray(observed)
    require(mu.shape == time.shape == observed.shape and observed.shape == (72,) and np.isfinite(observed).all(), "exactly three complete 24-frame training groups are required")
    training_level = float(np.mean(observed))
    attempts = []
    for start in range(8 if names else 1):
        initial = _initial_parameters(model_id, bounds, start, training_level)
        attempt = {"model_id": model_id, "start_index": start, "seed": SEEDS[start], "initial_parameters": initial, "parameters": dict(initial), "loss_evaluations": 0, "jacobian_evaluations": 0, "status": "failed", "training_mean_mse": None, "reason": None}

        def residual(vector):
            require(attempt["loss_evaluations"] < 2000, "registered loss evaluation budget exhausted", "model_execution_failed")
            attempt["loss_evaluations"] += 1
            values = dict(zip(names, (float(x) for x in vector)))
            prediction = evaluate_mean(model_id, values, mu, time)[1]
            attempt["parameters"] = values
            return (prediction - observed) / math.sqrt(72.0)

        def jacobian(vector):
            attempt["jacobian_evaluations"] += 1
            values = dict(zip(names, (float(x) for x in vector)))
            return _unbounded_jacobian(model_id, values, mu, time, 60.0)[1] / math.sqrt(72.0)

        try:
            if not names:
                errors = residual(np.empty(0))
                attempt.update(status="ok", training_mean_mse=float(errors @ errors), optimizer_status="parameter_free", optimizer_message="registered persistence family")
            else:
                result = least_squares(residual, [initial[name] for name in names], jac=jacobian, bounds=([bounds[name][0] for name in names], [bounds[name][1] for name in names]), method="trf", x_scale="jac", loss="linear", max_nfev=2000, ftol=1e-10, xtol=1e-10, gtol=1e-10)
                attempt.update(optimizer_status=int(result.status), optimizer_message=str(result.message))
                if np.isfinite(result.x).all():
                    attempt["parameters"] = dict(zip(names, (float(x) for x in result.x)))
                require(result.success and np.isfinite(result.x).all() and np.isfinite(result.fun).all(), "optimizer did not converge to a finite registered fit", "model_execution_failed")
                attempt.update(status="ok", training_mean_mse=float(result.fun @ result.fun))
        except (DevelopmentError, ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
            attempt.update(status="failed", training_mean_mse=None, reason=str(error))
        attempts.append(attempt)
    selected = _select_start(attempts, names)
    return {"model_id": model_id, "parameter_order": list(names), "bounds": bounds, "attempts": attempts, "status": "ok" if selected is not None else "failed", "selected_start_index": selected["start_index"] if selected else None, "parameters": selected["parameters"] if selected else None, "training_mean_mse": selected["training_mean_mse"] if selected else None, "reason": None if selected else "all registered starts failed"}


def fit_training(prepared, protocol):
    validate_protocol(protocol)
    checked = finalize_prepared(prepared["prefix_summaries"], prepared["training_frames"], prepared["development_coordinates"], protocol)
    require(checked["H_train"] == prepared["H_train"] and checked["S_train"] == prepared["S_train"], "training scales must be recomputed from training only")
    prefix = {row["group_id"]: row["mu"] for row in prepared["prefix_summaries"]}
    mu = np.asarray([prefix[row["group_id"]] for row in prepared["training_frames"]])
    time = np.asarray([row["time_min"] for row in prepared["training_frames"]])
    observed = np.asarray([row["observed_mean"] for row in prepared["training_frames"]])
    return {model_id: fit_family(model_id, mu, time, observed, prepared["H_train"]) for model_id in MODEL_IDS}


def predict_development(prepared, fits, protocol):
    validate_protocol(protocol)
    require(set(fits) == set(MODEL_IDS), "all six canonical families must remain present")
    checked = finalize_prepared(prepared["prefix_summaries"], prepared["training_frames"], prepared["development_coordinates"], protocol)
    require(checked["H_train"] == prepared["H_train"], "prediction bounds must use the authoritative training-only H_train")
    prefix = {row["group_id"]: row["mu"] for row in prepared["prefix_summaries"]}
    predictions, audit = [], []
    for model_id in MODEL_IDS:
        fit = fits[model_id]
        parameter_failure = None
        if fit["status"] == "ok":
            bounds = parameter_bounds(model_id, checked["H_train"])
            parameters = fit["parameters"]
            if not isinstance(parameters, dict) or set(parameters) != set(bounds) or not all(_number(value) and bounds[name][0] <= value <= bounds[name][1] for name, value in parameters.items()):
                parameter_failure = "checkpoint parameters are outside the registered protocol domains and training-only bounds"
        for coordinate in prepared["development_coordinates"]:
            row = {"model_id": model_id, **coordinate, "quantity": QUANTITY, "unit": "dimensionless_ratio", "status": "ok", "mean": None, "reason": None}
            detail = {"model_id": model_id, **coordinate, "unbounded_mean": None, "nonnegative_map_binds": None, "map_changed_value": None}
            if coordinate["time_min"] is None:
                row.update(status="missing_coordinate", reason="no positive measured response-time coordinate")
            elif fit["status"] != "ok":
                row.update(status="failed", reason=fit["reason"] or "canonical fitting failure")
            elif parameter_failure is not None:
                row.update(status="failed", reason=parameter_failure)
            else:
                try:
                    raw, mean, binds = evaluate_mean(model_id, fit["parameters"], prefix[coordinate["group_id"]], coordinate["time_min"])
                    row["mean"] = float(mean)
                    detail.update(unbounded_mean=float(raw), nonnegative_map_binds=bool(binds), map_changed_value=bool(raw < 0.0))
                except (DevelopmentError, ValueError, FloatingPointError) as error:
                    row.update(status="failed", reason=str(error))
            validate_document(protocol, row, "FramePrediction")
            predictions.append(row)
            audit.append(detail)
    require(len(predictions) == 288, "full six-by-two-by-24 prediction inventory required")
    return predictions, audit


def select_model(scores):
    require(set(scores) == set(MODEL_IDS), "cannot omit a failed canonical family")
    if not all(row["complete"] and _number(row["mean_mse"]) for row in scores.values()):
        return None
    best = min(row["mean_mse"] for row in scores.values())
    tied = [name for name in MODEL_IDS if scores[name]["mean_mse"] <= best + 1e-12 * (1.0 + abs(best))]
    return min(tied, key=lambda name: (len(PARAMETERS[name]), name))


def environment_record():
    return {"scientific_scope": SCIENTIFIC_SCOPE, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "optimizer_settings": dict(SETTINGS), "observation_transform": "identity raw individual max5/median ratios; prefix mean is a conditioning input", "limitations": list(LIMITATIONS)}


def validate_projection(protocol, projection, phase):
    validate_document(protocol, projection, "Projection")
    require(projection["phase"] == phase, "projection phase mismatch")
    expected_specs = group_specs(protocol)[:3] if phase == "training_values" else group_specs(protocol)[3:]
    require(projection["group_ids"] == [row["group_id"] for row in expected_specs], "projection biological-group assignment mismatch")
    windows = ("prefix", "response") if phase == "training_values" else ("prefix",) if phase == "development_inputs" else ("response",)
    for group in expected_specs:
        records = [row for row in projection["source_cells"] if row["group_id"] == group["group_id"]]
        validate_group_records(group, records, windows)
        times = [row for row in projection["source_times"] if row["group_id"] == group["group_id"]]
        if phase == "development_inputs":
            validate_group_records(group, times, ("response",), times_only=True)
        else:
            require(not times, "unexpected time-only records in value projection")
    return True


def prepare_projections(protocol, training, conditioning, training_sha, conditioning_sha):
    validate_projection(protocol, training, "training_values")
    validate_projection(protocol, conditioning, "development_inputs")
    require(training["approval_sha256"] == conditioning["approval_sha256"], "projection approvals differ")
    prefixes, training_frames, coordinates, audits, failures = [], [], [], [], []
    for group in group_specs(protocol):
        projection, checksum = (training, training_sha) if group["split"] == "train" else (conditioning, conditioning_sha)
        cells = [row for row in projection["source_cells"] if row["group_id"] == group["group_id"]]
        try:
            prefix, audit = derive_prefix(group, cells, checksum)
        except DevelopmentError as error:
            failures.append({"group_id": group["group_id"], "reason": str(error), "reason_code": error.code})
            audits.append(error.details)
            continue
        prefixes.append(prefix)
        audits.append(audit)
        if group["split"] == "train":
            response = [row for row in cells if row["window"] == "response"]
            training_frames.extend(derive_frames(group, response, response, prefix))
        else:
            times = [row for row in projection["source_times"] if row["group_id"] == group["group_id"]]
            coordinates.extend(derive_frames(group, None, times, prefix))
    details = {"scientific_scope": SCIENTIFIC_SCOPE, "prefix_summaries": prefixes, "training_frames": training_frames, "development_coordinates": coordinates, "prefix_audits": audits, "failures": failures}
    if failures:
        raise DevelopmentError("one or more registered prefix cohorts are unusable; no replacement or partial fit", "prefix_cohort_unusable", [row["group_id"] for row in failures], details)
    try:
        prepared = finalize_prepared(prefixes, training_frames, coordinates, protocol)
    except DevelopmentError as error:
        error.details = details
        raise
    prepared["prefix_audits"] = audits
    return prepared


def adapter_checks(protocol, training):
    records = []
    for group in group_specs(protocol)[:3]:
        cells = [dict(row) for row in training["source_cells"] if row["group_id"] == group["group_id"]]
        original, _ = derive_prefix(group, cells, "0" * 64)
        poisoned = [{**row, "max5": 999999.0, "median": 0.0} if row["window"] == "response" else dict(row) for row in cells]
        after, _ = derive_prefix(group, poisoned, "0" * 64)
        require(after == original, "new adapter depends on response intensities when forming the prefix")
        scaled = [{**row, "max5": None if row["max5"] is None else row["max5"] * 2.0, "median": None if row["median"] is None else row["median"] * 2.0} for row in cells]
        scaled_prefix, _ = derive_prefix(group, scaled, "0" * 64)
        require(scaled_prefix["eligible_cell_indices"] == original["eligible_cell_indices"] and np.allclose([scaled_prefix["mu"], scaled_prefix["s"]], [original["mu"], original["s"]], atol=1e-10, rtol=1e-8), "multiplicative intensity gauge changed prefix summaries")
        response = [row for row in cells if row["window"] == "response"]
        scaled_response = [row for row in scaled if row["window"] == "response"]
        a, b = derive_frames(group, response, response, original), derive_frames(group, scaled_response, scaled_response, scaled_prefix)
        for first, second in zip(a, b):
            require(first["valid_cell_indices"] == second["valid_cell_indices"], "multiplicative gauge changed response validity")
            require(np.allclose(first["observed_mean"], second["observed_mean"], atol=1e-10, rtol=1e-8), "multiplicative gauge changed mean-of-ratios targets")
        records.append({"group_id": group["group_id"], "prefix_response_poison_invariant": True, "joint_intensity_scale_invariant": True, "historical_processing_future_independence_established": False})
    return records


def safe_path(root, relative):
    root = Path(root).resolve()
    relative = Path(relative)
    require(not relative.is_absolute() and ".." not in relative.parts and bool(relative.parts), "artifact path must be explicit and root-relative")
    path = root / relative
    forbidden = [root / "data/native_law_v2/granados/sealed", root / "data/native_law_v2/granados/development/private"]
    require(not any(path.is_relative_to(directory) for directory in forbidden), "raw mixed files and private custody are forbidden")
    require(not any(root.joinpath(*relative.parts[:i]).is_symlink() for i in range(1, len(relative.parts) + 1)), "symlinked artifacts are forbidden")
    require(path.resolve().is_relative_to(root), "artifact escapes the project root")
    return path


def _artifact_path_policy(root, path, role):
    root = Path(root).resolve()
    projection_root = root / f"data/native_law_v2/granados/development_exports/{PROTOCOL_ID}"
    if role in {"training_projection", "development_inputs", "development_responses"}:
        require(path.is_relative_to(projection_root), "projection path is outside the registered export namespace")
    else:
        require(not path.is_relative_to(projection_root), "projection bytes cannot be opened as metadata or code")
        allowed = [root / "data/native_law_v2/granados/development", root / "outputs/native_population_development"]
        require(any(path.is_relative_to(directory) for directory in allowed) or path == root / "data/native_law_v2/development_protocol.json" or (role == "execution_code" and path.suffix == ".py"), "artifact is outside the declared development metadata/code scope")


def scoped_protocol(root, protocol_path):
    root = Path(root).resolve()
    path = Path(protocol_path)
    relative = path.relative_to(root) if path.is_absolute() else path
    path = safe_path(root, relative)
    _artifact_path_policy(root, path, "protocol_snapshot")
    return load_protocol(path)


def artifact_reference(root, path, role):
    root = Path(root).resolve()
    path = Path(path)
    relative = path.relative_to(root) if path.is_absolute() else path
    path = safe_path(root, relative)
    _artifact_path_policy(root, path, role)
    require(path.is_file(), "artifact is missing or is not a regular file")
    content = path.read_bytes()
    return {"path": relative.as_posix(), "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content), "role": role}


def read_artifact(root, reference, protocol, role=None, parse=True, byte_limit=64 * 1024 * 1024):
    validate_document(protocol, reference, "ArtifactRef")
    require(role is None or reference["role"] == role, "artifact role mismatch")
    path = safe_path(root, reference["path"])
    _artifact_path_policy(root, path, reference["role"])
    require(path.is_file() and path.stat().st_size == reference["bytes"] and 0 < reference["bytes"] <= byte_limit, "artifact byte inventory mismatch")
    content = path.read_bytes()
    require(hashlib.sha256(content).hexdigest() == reference["sha256"], "artifact digest mismatch", "source_mismatch")
    return strict_json(content) if parse else content


class ArtifactStore:
    def __init__(self, root, directory):
        self.root = Path(root).resolve()
        self.directory = Path(directory)
        require(self.directory.parts[:2] == ("outputs", "native_population_development"), "new development artifacts must stay in their separate output namespace")
        target = safe_path(self.root, self.directory)
        target.mkdir(parents=True, exist_ok=True)
        self.references = []

    def write(self, name, value, role):
        return self.write_bytes(name, canonical_bytes(value), role)

    def write_bytes(self, name, content, role):
        name = Path(name)
        require(not name.is_absolute() and ".." not in name.parts, "artifact name must remain inside the new run directory")
        path = safe_path(self.root, self.directory / name)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
        reference = {"path": path.relative_to(self.root).as_posix(), "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content), "role": role}
        self.references.append(reference)
        return reference


def code_artifacts(root):
    names = ["src/ystwin/analysis/native_population_development.py", "scripts/run_native_population_development.py", "tests/test_native_population_development.py"]
    return [artifact_reference(root, name, "execution_code") for name in names]


def write_implementation_configuration(root, protocol_path, output_directory, configuration_phase="before_training_value_export"):
    protocol = scoped_protocol(root, protocol_path)
    store = ArtifactStore(root, output_directory)
    implementation = code_artifacts(root)
    protocol_ref = artifact_reference(root, protocol_path, "protocol_snapshot")
    snapshots = [store.write_bytes(f"code/{Path(reference['path']).name}", read_artifact(root, reference, protocol, "execution_code", parse=False), "execution_code") for reference in implementation]
    configuration = {"document_type": "population_development_executor_configuration", "protocol_id": protocol["protocol_id"], "protocol_sha256": protocol_ref["sha256"], "code_artifacts": implementation, "code_snapshots": snapshots, "environment": environment_record(), "requested_phase": configuration_phase, "strong_claim_authorized": False}
    configuration_ref = store.write("execution_configuration.json", configuration, "execution_code")
    snapshot = store.write_bytes("development_protocol.json", safe_path(root, protocol_ref["path"]).read_bytes(), "protocol_snapshot")
    return {"scientific_scope": SCIENTIFIC_SCOPE, "status": "implementation_frozen_not_export_permission", "configuration": configuration_ref, "protocol": snapshot, "code_artifacts": implementation, "strong_claim_authorized": False}


def write_release_request(root, protocol_path, output_directory, exporter_paths, requester, configuration_phase="before_training_value_export"):
    root = Path(root).resolve()
    protocol = scoped_protocol(root, protocol_path)
    require(isinstance(requester, str) and bool(requester.strip()), "actual requester identity is required")
    require(exporter_paths, "actual exporter implementation artifacts are required, not a placeholder")
    require(all(Path(path).suffix == ".py" for path in exporter_paths), "exporter implementation must identify actual Python source, not a measurement file")
    implementation = code_artifacts(root)
    existing = safe_path(root, Path(output_directory) / "execution_configuration.json")
    if existing.is_file():
        configuration_ref = artifact_reference(root, existing, "execution_code")
        configuration = read_artifact(root, configuration_ref, protocol, "execution_code")
        snapshot = artifact_reference(root, Path(output_directory) / "development_protocol.json", "protocol_snapshot")
        current_protocol = artifact_reference(root, protocol_path, "protocol_snapshot")
        require(configuration["code_artifacts"] == implementation and configuration["environment"] == environment_record() and configuration["protocol_sha256"] == current_protocol["sha256"] == snapshot["sha256"] and configuration["requested_phase"] == configuration_phase, "existing implementation freeze differs from the current request or actual release phase")
    else:
        frozen = write_implementation_configuration(root, protocol_path, output_directory, configuration_phase)
        configuration_ref, snapshot = frozen["configuration"], frozen["protocol"]
    store = ArtifactStore(root, output_directory)
    request = {"document_type": "development_release_request", "protocol_id": PROTOCOL_ID, "protocol_sha256": snapshot["sha256"], "source_asset_id": ASSET_ID, "source_sha256": SOURCE_SHA256, "train_ids": [row["group_id"] for row in group_specs(protocol)[:3]], "development_ids": [row["group_id"] for row in group_specs(protocol)[3:]], "exporter_code": [artifact_reference(root, path, "execution_code") for path in exporter_paths], "executor_code": [*implementation, configuration_ref], "requested_phases": ["training_and_conditioning_export", "development_response_export"], "requester": requester}
    validate_document(protocol, request, "DevelopmentReleaseRequest")
    request_ref = store.write("release_request.json", request, "operational_event_log")
    return {"scientific_scope": SCIENTIFIC_SCOPE, "status": "request_only_not_export_permission", "request": request_ref, "configuration": configuration_ref, "protocol": snapshot, "code_artifacts": implementation, "strong_claim_authorized": False}


def seal_training_freeze(protocol, prepared, fits, provenance, store):
    predictions, output_audit = predict_development(prepared, fits, protocol)
    official_attempts, selected_refs, checkpoint_refs = [], {}, []
    for model_id in MODEL_IDS:
        fit = fits[model_id]
        for attempt in fit["attempts"]:
            checkpoint = None
            if attempt["status"] == "ok":
                body = {"document_type": "population_development_model_checkpoint", "protocol_id": PROTOCOL_ID, "scientific_scope": SCIENTIFIC_SCOPE, "model_id": model_id, "parameter_order": list(PARAMETERS[model_id]), "parameters": attempt["parameters"], "initial_parameters": attempt["initial_parameters"], "bounds": parameter_bounds(model_id, prepared["H_train"]), "start_index": attempt["start_index"], "seed": attempt["seed"], "training_mean_mse": attempt["training_mean_mse"], "training_projection_sha256": provenance["training_projection_sha256"], "H_train": prepared["H_train"], "S_train": prepared["S_train"], "observation_transform": "identity", "output_map": SETTINGS["output_map"], "selected_for_family": attempt["start_index"] == fit["selected_start_index"], "strong_claim_authorized": False}
                checkpoint = store.write(f"checkpoints/{model_id}_start_{attempt['start_index']}.json", body, "model_checkpoint")
                checkpoint_refs.append(checkpoint)
                if body["selected_for_family"]:
                    selected_refs[model_id] = checkpoint
            official = {key: attempt[key] for key in ("model_id", "start_index", "seed", "loss_evaluations", "status", "parameters", "training_mean_mse", "reason")}
            official["checkpoint"] = checkpoint
            validate_document(protocol, official, "Attempt")
            official_attempts.append(official)
        selected_refs.setdefault(model_id, None)
    attempt_log = store.write("attempt_log.json", {"document_type": "population_development_complete_attempt_log", "protocol_id": PROTOCOL_ID, "scientific_scope": SCIENTIFIC_SCOPE, "families": fits, "environment": environment_record(), "strong_claim_authorized": False}, "attempt_log")
    conditioning = store.write("conditioning_checkpoint.json", prepared, "model_checkpoint")
    training_audit = []
    prefix_means = {row["group_id"]: row["mu"] for row in prepared["prefix_summaries"]}
    for model_id in MODEL_IDS:
        for frame in prepared["training_frames"]:
            detail = {"model_id": model_id, "group_id": frame["group_id"], "frame_1based": frame["frame_1based"], "time_min": frame["time_min"], "unbounded_mean": None, "mean": None, "nonnegative_map_binds": None, "map_changed_value": None}
            if fits[model_id]["status"] == "ok":
                raw, mean, binds = evaluate_mean(model_id, fits[model_id]["parameters"], prefix_means[frame["group_id"]], frame["time_min"])
                detail.update(unbounded_mean=float(raw), mean=float(mean), nonnegative_map_binds=bool(binds), map_changed_value=bool(raw < 0.0))
            training_audit.append(detail)
    map_audit = store.write("nonnegative_map_audit.json", {"scientific_scope": SCIENTIFIC_SCOPE, "output_map": SETTINGS["output_map"], "training": training_audit, "development": output_audit}, "model_checkpoint")
    environment = {**environment_record(), "provenance": provenance, "selected_checkpoints": selected_refs, "attempt_log": attempt_log, "conditioning_checkpoint": conditioning, "nonnegative_map_audit": map_audit}
    freeze = {"document_type": "development_prediction_freeze_not_final_validation", "protocol_id": PROTOCOL_ID, "approval_sha256": provenance["approval_sha256"], "training_projection_sha256": provenance["training_projection_sha256"], "development_inputs_sha256": provenance["development_inputs_sha256"], "code_artifacts": provenance.get("frozen_code_artifacts", provenance["code_artifacts"]), "environment": environment, "prefix_summaries": prepared["prefix_summaries"], "H_train": prepared["H_train"], "S_train": prepared["S_train"], "attempts": official_attempts, "predictions": predictions, "operational_event_reference": provenance["operational_event_reference"], "strong_claim_authorized": False}
    validate_freeze(protocol, freeze)
    reference = store.write("development_prediction_freeze.json", freeze, "development_predictions")
    return {"scientific_scope": SCIENTIFIC_SCOPE, "status": "training_frozen_development_responses_not_accessed", "prediction_freeze": reference, "attempt_log": attempt_log, "conditioning_checkpoint": conditioning, "nonnegative_map_audit": map_audit, "selected_checkpoints": selected_refs, "checkpoint_artifacts": checkpoint_refs, "training_results": {model_id: {"status": fits[model_id]["status"], "parameters": fits[model_id]["parameters"], "training_mean_mse": fits[model_id]["training_mean_mse"], "selected_start_index": fits[model_id]["selected_start_index"]} for model_id in MODEL_IDS}, "development_prediction_rows": len(predictions), "selected_development_model": None, "strong_claim_authorized": False}


def validate_freeze(protocol, freeze):
    validate_document(protocol, freeze, "DevelopmentPredictionFreeze")
    specs = {row["group_id"]: row for row in group_specs(protocol)}
    prefixes = {row["group_id"]: row for row in freeze["prefix_summaries"]}
    require(set(prefixes) == set(specs), "freeze must retain all five unique prefix groups")
    for group_id, prefix in prefixes.items():
        require(all(0 <= cell < specs[group_id]["source_rows"] for cell in prefix["eligible_cell_indices"]), "frozen cohort row is outside its source group")
        expected_sha = freeze["training_projection_sha256"] if specs[group_id]["split"] == "train" else freeze["development_inputs_sha256"]
        require(prefix["input_projection_sha256"] == expected_sha, "prefix summary is not bound to its correct input projection")
    attempts = freeze["attempts"]
    require({row["model_id"] for row in attempts} == set(MODEL_IDS), "freeze omitted a canonical family")
    chosen = {}
    for model_id in MODEL_IDS:
        rows = [row for row in attempts if row["model_id"] == model_id]
        expected_starts = list(range(8 if PARAMETERS[model_id] else 1))
        require(sorted(row["start_index"] for row in rows) == expected_starts, "registered starts are missing, duplicated or expanded")
        bounds = parameter_bounds(model_id, freeze["H_train"])
        for row in rows:
            require(row["seed"] == SEEDS[row["start_index"]] and set(row["parameters"]) == set(PARAMETERS[model_id]), "attempt seed or parameter contract mismatch")
            require(all(bounds[key][0] <= value <= bounds[key][1] for key, value in row["parameters"].items()), "attempt parameters escape frozen training-only bounds")
        chosen[model_id] = _select_start(rows, PARAMETERS[model_id])
    expected = {(model_id, group["group_id"], frame) for model_id in MODEL_IDS for group in group_specs(protocol)[3:] for frame in RESPONSE_FRAMES}
    actual = [(row["model_id"], row["group_id"], row["frame_1based"]) for row in freeze["predictions"]]
    require(set(actual) == expected and len(actual) == len(set(actual)), "freeze prediction inventory is missing or duplicated")
    coordinates = {}
    for row in freeze["predictions"]:
        key = row["group_id"], row["frame_1based"]
        if key in coordinates:
            require(_equal(coordinates[key], row["time_min"]), "canonical families use different response coordinates")
        coordinates[key] = row["time_min"]
        selected = chosen[row["model_id"]]
        if row["time_min"] is None:
            require(row["status"] == "missing_coordinate", "unavailable coordinate must remain explicit")
        elif selected is None:
            require(row["status"] == "failed", "failed family cannot supply a successful prediction")
        elif row["status"] == "ok":
            prediction = float(evaluate_mean(row["model_id"], selected["parameters"], prefixes[row["group_id"]]["mu"], row["time_min"])[1])
            require(math.isclose(prediction, row["mean"], abs_tol=1e-10, rel_tol=1e-8), "prediction differs from the frozen training-selected fit")
    return True


def build_development_report(protocol, freeze, observed_frames, freeze_sha, responses_sha):
    validate_freeze(protocol, freeze)
    development_ids = [row["group_id"] for row in group_specs(protocol)[3:]]
    expected = {(group_id, frame) for group_id in development_ids for frame in RESPONSE_FRAMES}
    observations = {(row["group_id"], row["frame_1based"]): row for row in observed_frames}
    require(set(observations) == expected and len(observed_frames) == 48, "all 48 registered development observations must remain present")
    predictions = {(row["model_id"], row["group_id"], row["frame_1based"]): row for row in freeze["predictions"]}
    evaluations, model_scores = [], {}
    observation_incomplete = False
    for model_id in MODEL_IDS:
        model_groups = []
        for group_id in development_ids:
            frames = []
            for frame in RESPONSE_FRAMES:
                observed, predicted = observations[group_id, frame], predictions[model_id, group_id, frame]
                require(_equal(observed["time_min"], predicted["time_min"]), "observed timing differs from the pre-response freeze")
                reason, status = None, "quantified"
                if predicted["time_min"] is None:
                    status, reason = "missing_coordinate", "registered measured coordinate unavailable"
                    observation_incomplete = True
                elif observed["observed_mean"] is None:
                    status, reason = "missing_observation", "no valid ratio in the fixed prefix cohort at this frame"
                    observation_incomplete = True
                elif predicted["status"] != "ok":
                    status, reason = "failed_prediction", predicted["reason"]
                squared, absolute = None, None
                if status == "quantified":
                    error = predicted["mean"] - observed["observed_mean"]
                    squared, absolute = error * error, abs(error)
                row = {"model_id": model_id, "group_id": group_id, "frame_1based": frame, "time_min": predicted["time_min"], "observed_mean": observed["observed_mean"], "predicted_mean": predicted["mean"], "prefix_cells": observed["prefix_cells"], "valid_cells": observed["valid_cells"], "source_missing_cells": observed["source_missing_cells"], "invalid_cells": observed["invalid_cells"], "cell_sample_sd": observed["cell_sample_sd"], "squared_error": squared, "absolute_error": absolute, "status": status, "reason": reason}
                validate_document(protocol, row, "FrameEvaluation")
                frames.append(row)
            complete = all(row["status"] == "quantified" for row in frames)
            mse = float(np.mean([row["squared_error"] for row in frames])) if complete else None
            evaluation = {"model_id": model_id, "group_id": group_id, "complete": complete, "mean_mse": mse, "rmse": math.sqrt(mse) if complete else None, "mae": float(np.mean([row["absolute_error"] for row in frames])) if complete else None, "frames": frames}
            evaluations.append(evaluation)
            model_groups.append(evaluation)
        complete = all(row["complete"] for row in model_groups)
        model_scores[model_id] = {"complete": complete, "mean_mse": float(np.mean([row["mean_mse"] for row in model_groups])) if complete else None}
    selected = select_model(model_scores)
    status = "incomplete_observation_coverage" if observation_incomplete else "incomplete_model_execution" if selected is None else "complete_development_diagnostics"
    report = {"document_type": "development_diagnostics_not_biological_validation", "protocol_id": PROTOCOL_ID, "status": status, "prediction_freeze_sha256": freeze_sha, "development_responses_sha256": responses_sha, "group_evaluations": evaluations, "selected_development_model": selected, "limitations": list(LIMITATIONS), "strong_claim_authorized": False, "final_test_defined": False}
    validate_document(protocol, report, "DevelopmentReport")
    return report


def admission_contract():
    return {
        "document_type": "population_development_projection_allowlist",
        "required_fields": ["document_type", "protocol_id", "protocol_sha256", "phase", "request", "approval", "configuration", "source_replay", "artifacts", "custodian_identity", "operational_event_reference", "resolved_blocker_sha256", "prediction_freeze_sha256"],
        "artifact_fields": "request, approval, configuration, source_replay and every artifacts entry are exact embedded-protocol ArtifactRef records",
        "initial_phase": "training_and_conditioning_export",
        "initial_artifact_roles": ["training_projection", "development_inputs"],
        "response_phase": "development_response_export",
        "response_artifact_roles": ["development_responses"],
        "projection_root": f"data/native_law_v2/granados/development_exports/{PROTOCOL_ID}",
        "source_replay_document_type": "population_development_projection_source_replay",
        "source_replay_required_fields": ["document_type", "protocol_id", "protocol_sha256", "source_asset_id", "source_sha256", "approval_sha256", "custodian_identity", "phase", "operational_event_reference", "resolved_blocker_sha256", "prediction_freeze_sha256", "artifact_checks"],
        "source_replay_artifact_check_fields": ["artifact", "source_scalar_count", "all_source_scalars_equal", "nulls_preserved", "source_rows_and_pointers_verified", "group_window_inventory_verified"],
        "source_replay_scalar_counts": {"training_projection": 51876, "development_inputs": 29568, "development_responses": 25344},
        "authority_limit": "Only an actual independent custodian may issue this operational source-replay allowlist. The executor never creates approvals, replay receipts or production signatures. These records do not authorize a strong biological claim.",
    }


def _metadata_file(root, path):
    path = Path(path)
    relative = path.relative_to(Path(root).resolve()) if path.is_absolute() else path
    path = safe_path(root, relative)
    _artifact_path_policy(root, path, "operational_event_log")
    require(path.is_file() and path.stat().st_size <= 4 * 1024 * 1024, "required public custody metadata is missing or too large", "custody_not_approved")
    raw = path.read_bytes()
    return strict_json(raw), hashlib.sha256(raw).hexdigest()


def _phase_b_data_v1(root, protocol, manifest, manifest_reference):
    expected_fields = {"document_type", "protocol_id", "protocol_sha256", "operational_phase", "release_request", "custody_decision", "source_asset_id", "source_sha256", "learner_allowlist", "source_binding_artifacts", "frozen_executor_artifacts", "exporter_code_artifacts", "protocol_snapshot", "approved_prefix_cells_by_group", "all_planned_coordinates_and_targets_available_custodian_checked", "coverage_details", "development_response_projection", "development_prediction_freeze_received", "development_responses_released", "further_release_requires", "exposed_biological_group_ids", "unexported_scope", "original_strong_artifact_sha256_unchanged", "strong_claim_authorized", "final_test_defined", "limitations"}
    require(set(manifest) == expected_fields and manifest["document_type"] == "approved_development_phase_b_data_allowlist_not_biological_authorization", "phase-B release-manifest-v1 schema mismatch")
    require(manifest["protocol_id"] == PROTOCOL_ID and manifest["source_asset_id"] == ASSET_ID and manifest["source_sha256"] == SOURCE_SHA256, "phase-B source/protocol identity mismatch", "source_mismatch")
    require(manifest["operational_phase"] == "training_and_conditioning_export_complete" and manifest["development_responses_released"] is False and manifest["development_response_projection"] is None and manifest["strong_claim_authorized"] is False and manifest["final_test_defined"] is False, "phase-B manifest does not authorize the requested limited phase")
    protocol_bytes = read_artifact(root, manifest["protocol_snapshot"], protocol, "protocol_snapshot", parse=False)
    require(validate_protocol(strict_json(protocol_bytes)) == protocol and manifest["protocol_sha256"] == manifest["protocol_snapshot"]["sha256"], "phase-B protocol snapshot differs from the active protocol")
    request_ref = {**manifest["release_request"], "role": "operational_event_log"}
    require(set(manifest["release_request"]) == {"path", "sha256", "bytes"}, "phase-B request reference must use its declared v1 layout")
    request = read_artifact(root, request_ref, protocol, "operational_event_log")
    approval = read_artifact(root, manifest["custody_decision"], protocol, "operational_approval")
    validate_document(protocol, request, "DevelopmentReleaseRequest")
    validate_document(protocol, approval, "DevelopmentCustodyDecision")
    require(request["protocol_sha256"] == manifest["protocol_sha256"] and approval["request_sha256"] == request_ref["sha256"] and approval["decision"] == "approve_development_exports", "phase-B custody decision does not bind the approved request", "custody_not_approved")
    require(approval["custodian_identity"] != request["requester"] and set(approval["permitted_group_ids"]) == {row["group_id"] for row in group_specs(protocol)}, "phase-B independent custody/group assignment mismatch")
    event_path, separator, event_id = approval["operational_event_reference"].partition("#")
    require(separator and event_id == "phase-b-approved", "phase-B v1 approval-event reference mismatch")
    event, event_sha = _metadata_file(root, event_path)
    require(event["document_type"] == "local_operational_custody_event_not_external_registry" and event["event_id"] == event_id and event["request_sha256"] == request_ref["sha256"] and event["approval_sha256"] == manifest["custody_decision"]["sha256"] and event["values_released_before_this_approval"] is False and event["strong_claim_authorized"] is False, "phase-B approval event is not a matching operational authorization", "custody_not_approved")
    blocker_path = safe_path(root, "data/native_law_v2/granados/development/operational_blocker.json")
    if blocker_path.exists():
        blocker, blocker_sha = _metadata_file(root, blocker_path)
        validate_document(protocol, blocker, "DevelopmentBlocker")
        require(event["superseded_blocker_sha256"] == blocker_sha, "phase-B approval did not explicitly supersede the recorded blocker", "custody_not_approved")
    ids = [row["group_id"] for row in group_specs(protocol)]
    require(manifest["exposed_biological_group_ids"] == ids and set(manifest["approved_prefix_cells_by_group"]) == set(ids), "phase-B biological exposure inventory mismatch")
    references = manifest["learner_allowlist"]
    require([ref["role"] for ref in references] == ["training_projection", "development_inputs"] and len(manifest["source_binding_artifacts"]) == 2, "phase-B allowlist must contain exactly the two approved projection roles")
    bindings = {}
    for reference, binding_ref, phase, expected_ids, expected_counts in zip(references, manifest["source_binding_artifacts"], ("training_values", "development_inputs"), (ids[:3], ids[3:]), ((17292, 0, 51876), (7040, 8448, 29568))):
        validate_document(protocol, reference, "ArtifactRef")
        path = safe_path(root, reference["path"])
        _artifact_path_policy(root, path, reference["role"])
        require(path.suffix == ".json" and path.is_file() and path.stat().st_size == reference["bytes"], "phase-B projection is missing or changed", "source_mismatch")
        binding = read_artifact(root, binding_ref, protocol, "source_binding_audit")
        require(binding["document_type"] == "source_scalar_replay_binding_not_biological_validation" and binding["protocol_id"] == PROTOCOL_ID and binding["protocol_sha256"] == manifest["protocol_sha256"] and binding["phase"] == phase and binding["approval_sha256"] == manifest["custody_decision"]["sha256"] and binding["source_asset_id"] == ASSET_ID and binding["source_sha256"] == SOURCE_SHA256 and binding["group_ids"] == expected_ids and binding["development_response_intensities_in_payload"] is False and binding["strong_claim_authorized"] is False, "phase-B source-scalar binding does not match the authorized phase", "source_mismatch")
        checks = binding["verification"]
        require(tuple(checks[key] for key in ("source_cell_records", "source_time_records", "source_scalar_checks")) == expected_counts, "phase-B source replay omitted registered scalar/null records", "source_mismatch")
        require(isinstance(binding["ordered_payload_sha256"], str) and re.fullmatch("[0-9a-f]{64}", binding["ordered_payload_sha256"]) is not None, "phase-B ordered-payload checksum is missing")
        audit = read_artifact(root, binding["source_binding_audit"], protocol, "source_binding_audit")
        require(audit["protocol_sha256"] == manifest["protocol_sha256"] and audit["source_sha256"] == SOURCE_SHA256 and binding["source_binding_audit"]["sha256"] == event["source_binding_audit_sha256"], "phase-B source audit ancestry mismatch", "source_mismatch")
        bindings[reference["role"]] = {"reference": binding_ref, "record": binding}
    require(request["executor_code"] == manifest["frozen_executor_artifacts"] and request["exporter_code"] == manifest["exporter_code_artifacts"], "phase-B pre-value code inventory differs from the request")
    frozen_refs = manifest["frozen_executor_artifacts"]
    require(len(frozen_refs) == 2 and Path(frozen_refs[0]["path"]).suffix == ".py" and Path(frozen_refs[1]["path"]).name == "executor_environment.json", "unsupported phase-B v1 frozen executor inventory")
    read_artifact(root, frozen_refs[0], protocol, "execution_code", parse=False)
    frozen_environment = read_artifact(root, frozen_refs[1], protocol, "execution_code")
    require(frozen_environment["document_type"] == "pre_value_executor_implementation_freeze_not_a_fit" and frozen_environment["protocol_sha256"] == manifest["protocol_sha256"] and frozen_environment["frozen_implementation"] == frozen_refs[0] and frozen_environment["original_implementation_sha256_at_freeze"] == frozen_refs[0]["sha256"] and frozen_environment["optimizer_settings"] == SETTINGS and frozen_environment["canonical_model_ids"] == list(MODEL_IDS), "phase-B frozen scientific recipe or code binding is incompatible")
    require(frozen_environment["versions"] == {key: environment_record()[key] for key in ("python", "numpy", "scipy")}, "phase-B runtime/PRNG versions differ from the frozen recipe")
    for reference in manifest["exporter_code_artifacts"]:
        read_artifact(root, reference, protocol, "execution_code", parse=False)
    return {"adapter_version": "phase_b_release_manifest_v1", "data_ready": True, "executor_ready": False, "phase_b_manifest": manifest_reference, "approval_sha256": manifest["custody_decision"]["sha256"], "approval_artifact": manifest["custody_decision"], "request_artifact": request_ref, "protocol_artifact": manifest["protocol_snapshot"], "historical_executor_artifacts": frozen_refs, "source_binding_artifacts": manifest["source_binding_artifacts"], "source_bindings": bindings, "approved_prefix_cells_by_group": manifest["approved_prefix_cells_by_group"], "approval_event_sha256": event_sha, "operational_event_reference": approval["operational_event_reference"], "source": dict(protocol["source"]), "approved_data_artifacts": references, "refresh_requirement": "POP001/POP003/POP004 fixes require a new code/configuration freeze, a new typed request, and an independently issued executor_refresh_admission_v1; existing source projections remain approved and must not be regenerated"}


def executor_refresh_contract():
    return {"document_type": "population_development_executor_refresh_admission_v1", "schema_version": 1, "required_fields": ["document_type", "schema_version", "protocol_id", "protocol_sha256", "phase_b_manifest", "previous_approval_sha256", "request", "approval", "configuration", "custodian_identity", "operational_event_reference", "reuse_existing_projections", "scientific_recipe_changed", "data_reexport_permitted", "development_responses_released", "strong_claim_authorized"], "authority": "An independent custodian must issue the matching protocol-defined DevelopmentCustodyDecision and this operational refresh record. A request or caller-supplied success flag is not approval."}


def write_executor_refresh_request(root, protocol_path, manifest_path, output_directory, requester):
    protocol, _, provenance = _admission(root, protocol_path, manifest_path, "training_and_conditioning_export")
    require(provenance.get("adapter_version") == "phase_b_release_manifest_v1", "executor refresh requires the validated phase-B release-manifest-v1 adapter")
    manifest = read_artifact(root, provenance["phase_b_manifest"], protocol, "operational_event_log")
    request = write_release_request(root, protocol_path, output_directory, [reference["path"] for reference in manifest["exporter_code_artifacts"]], requester, "after_phase_b_export_before_native_fitting")
    store = ArtifactStore(root, output_directory)
    context = {"document_type": "population_development_executor_refresh_request_v1", "schema_version": 1, "protocol_id": PROTOCOL_ID, "scientific_scope": SCIENTIFIC_SCOPE, "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "phase_b_manifest": provenance["phase_b_manifest"], "previous_request": provenance["request_artifact"], "previous_approval": provenance["approval_artifact"], "previous_executor_artifacts": provenance["historical_executor_artifacts"], "new_request": request["request"], "new_configuration": request["configuration"], "fixes": ["POP001 authoritative protocol-domain validation before inference", "POP003 exact family set with protocol-ordered replay", "POP004 explicit phase-B manifest-v1 administrative adapter"], "reuse_existing_projections": provenance["approved_data_artifacts"], "scientific_recipe_changed": False, "native_fits_performed_by_this_request": False, "development_responses_accessed": False, "data_reexport_requested": False, "required_custodian_admission": executor_refresh_contract(), "strong_claim_authorized": False}
    context_ref = store.write("executor_refresh_request.json", context, "operational_event_log")
    return {**request, "status": "executor_refresh_requested_not_admitted", "refresh_request": context_ref, "data_ready": True, "executor_ready": False, "reuse_existing_projections": provenance["approved_data_artifacts"], "required_custodian_admission": executor_refresh_contract()}


def _apply_executor_refresh_v1(root, protocol, provenance, refresh_path):
    refresh, refresh_sha = _metadata_file(root, refresh_path)
    contract = executor_refresh_contract()
    require(set(refresh) == set(contract["required_fields"]) and refresh["document_type"] == contract["document_type"] and refresh["schema_version"] == 1, "executor-refresh-admission-v1 schema mismatch", "custody_not_approved")
    require(refresh["protocol_id"] == PROTOCOL_ID and refresh["protocol_sha256"] == provenance["protocol_artifact"]["sha256"] and refresh["phase_b_manifest"] == provenance["phase_b_manifest"] and refresh["previous_approval_sha256"] == provenance["approval_sha256"], "executor refresh does not bind the existing phase-B data approval", "custody_not_approved")
    require(refresh["reuse_existing_projections"] == provenance["approved_data_artifacts"] and all(refresh[key] is False for key in ("scientific_recipe_changed", "data_reexport_permitted", "development_responses_released", "strong_claim_authorized")), "executor refresh cannot change the scientific recipe or release additional data")
    request = read_artifact(root, refresh["request"], protocol, "operational_event_log")
    approval = read_artifact(root, refresh["approval"], protocol, "operational_approval")
    validate_document(protocol, request, "DevelopmentReleaseRequest")
    validate_document(protocol, approval, "DevelopmentCustodyDecision")
    require(request["protocol_sha256"] == refresh["protocol_sha256"] and approval["request_sha256"] == refresh["request"]["sha256"] and approval["decision"] == "approve_development_exports" and approval["custodian_identity"] == refresh["custodian_identity"] != request["requester"], "executor refresh lacks a matching independent custody decision", "custody_not_approved")
    require(set(approval["permitted_group_ids"]) == {row["group_id"] for row in group_specs(protocol)} and bool(refresh["operational_event_reference"]), "executor refresh group/event contract mismatch")
    configuration = read_artifact(root, refresh["configuration"], protocol, "execution_code")
    current = code_artifacts(root)
    require(configuration["requested_phase"] == "after_phase_b_export_before_native_fitting" and configuration["code_artifacts"] == current and configuration["environment"] == environment_record() and configuration["protocol_sha256"] == refresh["protocol_sha256"] and request["executor_code"] == [*current, refresh["configuration"]], "current executor/configuration has not been newly frozen and admitted", "executor_admission_refresh_required")
    snapshots = configuration["code_snapshots"]
    require(len(snapshots) == len(current), "new executor source snapshots are incomplete")
    for original, snapshot in zip(current, snapshots):
        require(snapshot["sha256"] == original["sha256"] and snapshot["bytes"] == original["bytes"], "new executor source snapshot mismatch")
        read_artifact(root, snapshot, protocol, "execution_code", parse=False)
    provenance.update(executor_ready=True, executor_admission_sha256=refresh_sha, executor_admission_path=str(refresh_path), executor_request_artifact=refresh["request"], executor_approval_artifact=refresh["approval"], configuration_artifact=refresh["configuration"], code_artifacts=current, frozen_code_artifacts=snapshots, operational_event_reference=refresh["operational_event_reference"])
    return provenance


def _admission(root, protocol_path, allowlist_path, phase, freeze_sha=None, executor_admission_path=None):
    root = Path(root).resolve()
    protocol = scoped_protocol(root, protocol_path)
    protocol_ref = artifact_reference(root, protocol_path, "protocol_snapshot")
    blocker_path = safe_path(root, "data/native_law_v2/granados/development/operational_blocker.json")
    blocker, blocker_sha = (None, None)
    if blocker_path.is_file():
        blocker, blocker_sha = _metadata_file(root, blocker_path)
        validate_document(protocol, blocker, "DevelopmentBlocker")
    if allowlist_path is None:
        released = safe_path(root, "data/native_law_v2/granados/development/release_manifest.json")
        if released.is_file():
            allowlist_path = released
        else:
            raise DevelopmentError(blocker["reason"] if blocker else "missing independently issued source-bound projection allowlist", "custody_not_approved")
    allowlist, allowlist_sha = _metadata_file(root, allowlist_path)
    if allowlist.get("document_type") == "approved_development_phase_b_data_allowlist_not_biological_authorization":
        require(phase == "training_and_conditioning_export" and freeze_sha is None, "phase-B release contains no development-response authorization", "custody_not_approved")
        manifest_ref = artifact_reference(root, allowlist_path, "operational_event_log")
        provenance = _phase_b_data_v1(root, protocol, allowlist, manifest_ref)
        require(allowlist["protocol_sha256"] == protocol_ref["sha256"], "active protocol differs from the source-bound phase-B release")
        if executor_admission_path is not None:
            provenance = _apply_executor_refresh_v1(root, protocol, provenance, executor_admission_path)
        return protocol, {"artifacts": allowlist["learner_allowlist"]}, provenance
    contract = admission_contract()
    require(set(allowlist) == set(contract["required_fields"]) and allowlist["document_type"] == contract["document_type"], "unsupported custodian allowlist contract", "custody_not_approved")
    require(allowlist["protocol_id"] == PROTOCOL_ID and allowlist["protocol_sha256"] == protocol_ref["sha256"] and allowlist["phase"] == phase, "custodian approval does not match protocol bytes or requested phase", "custody_not_approved")
    require(allowlist["prediction_freeze_sha256"] == freeze_sha, "response release is not bound to the requested prediction freeze", "custody_not_approved")
    request = read_artifact(root, allowlist["request"], protocol, parse=True, byte_limit=1024 * 1024)
    approval = read_artifact(root, allowlist["approval"], protocol, role="operational_approval", byte_limit=1024 * 1024)
    validate_document(protocol, request, "DevelopmentReleaseRequest")
    validate_document(protocol, approval, "DevelopmentCustodyDecision")
    require(request["protocol_sha256"] == protocol_ref["sha256"] and approval["request_sha256"] == allowlist["request"]["sha256"], "request/approval/protocol binding mismatch", "custody_not_approved")
    require(approval["decision"] == "approve_development_exports" and approval["custodian_identity"] == allowlist["custodian_identity"] and approval["custodian_identity"] != request["requester"], "independent custodian approval is absent", "custody_not_approved")
    require(set(approval["permitted_group_ids"]) == {row["group_id"] for row in group_specs(protocol)}, "custodian group allowlist differs from all five registered groups")
    configuration = read_artifact(root, allowlist["configuration"], protocol, role="execution_code", byte_limit=1024 * 1024)
    current_code = code_artifacts(root)
    require(configuration["document_type"] == "population_development_executor_configuration" and configuration["protocol_sha256"] == protocol_ref["sha256"] and configuration["code_artifacts"] == current_code and configuration["environment"] == environment_record(), "implementation, PRNG/version or optimizer differs from pre-export configuration", "custody_not_approved")
    require(request["executor_code"] == [*current_code, allowlist["configuration"]], "request did not freeze the actual executor and configuration", "custody_not_approved")
    snapshots = configuration["code_snapshots"]
    require(len(snapshots) == len(current_code), "frozen implementation source copies are incomplete")
    for original, snapshot in zip(current_code, snapshots):
        require(snapshot["sha256"] == original["sha256"] and snapshot["bytes"] == original["bytes"] and Path(snapshot["path"]).parts[:2] == ("outputs", "native_population_development"), "frozen implementation copy differs from the executed code")
        read_artifact(root, snapshot, protocol, "execution_code", parse=False, byte_limit=4 * 1024 * 1024)
    for reference in request["exporter_code"]:
        require(Path(reference["path"]).suffix == ".py", "exporter reference is not implementation code")
        read_artifact(root, reference, protocol, role="execution_code", parse=False, byte_limit=4 * 1024 * 1024)
    replay = read_artifact(root, allowlist["source_replay"], protocol, role="source_binding_audit", byte_limit=4 * 1024 * 1024)
    require(set(replay) == set(contract["source_replay_required_fields"]) and replay["document_type"] == contract["source_replay_document_type"], "projection source-replay receipt is absent; coverage metadata alone is insufficient", "custody_not_approved")
    for key, value in {"protocol_id": PROTOCOL_ID, "protocol_sha256": protocol_ref["sha256"], "source_asset_id": ASSET_ID, "source_sha256": SOURCE_SHA256, "approval_sha256": allowlist["approval"]["sha256"], "custodian_identity": approval["custodian_identity"], "phase": phase, "prediction_freeze_sha256": freeze_sha}.items():
        require(replay[key] == value, "source-replay receipt binding mismatch", "source_mismatch")
    require(bool(replay["operational_event_reference"]) and bool(allowlist["operational_event_reference"]), "custody operational event references are required", "custody_not_approved")
    if blocker is not None:
        require(allowlist["resolved_blocker_sha256"] == blocker_sha and replay["resolved_blocker_sha256"] == blocker_sha, blocker["reason"], "custody_not_approved")
    else:
        require(allowlist["resolved_blocker_sha256"] is None and replay["resolved_blocker_sha256"] is None, "unknown blocker-resolution history", "custody_not_approved")
    roles = ["training_projection", "development_inputs"] if phase == "training_and_conditioning_export" else ["development_responses"]
    require(phase in {"training_and_conditioning_export", "development_response_export"}, "unapproved release phase")
    require([row["role"] for row in allowlist["artifacts"]] == roles and len(replay["artifact_checks"]) == len(roles), "phase must allowlist only its exact projection roles")
    projection_root = root / contract["projection_root"]
    for reference, check in zip(allowlist["artifacts"], replay["artifact_checks"]):
        validate_document(protocol, reference, "ArtifactRef")
        path = safe_path(root, reference["path"])
        require(path.is_relative_to(projection_root) and path.suffix == ".json", "projection is outside the protocol's administrative export root", "custody_not_approved")
        require(set(check) == set(contract["source_replay_artifact_check_fields"]) and check["artifact"] == reference, "source replay does not bind the exact projection bytes and path", "source_mismatch")
        require(check["source_scalar_count"] == contract["source_replay_scalar_counts"][reference["role"]], "source replay did not cover every approved scalar/null", "source_mismatch")
        require(all(check[key] is True for key in ("all_source_scalars_equal", "nulls_preserved", "source_rows_and_pointers_verified", "group_window_inventory_verified")), "source replay or null/row/window preservation failed", "source_mismatch")
        require(path.is_file() and path.stat().st_size == reference["bytes"], "approved projection file is missing or changed", "source_mismatch")
    provenance = {"approval_sha256": allowlist["approval"]["sha256"], "protocol_artifact": protocol_ref, "request_artifact": allowlist["request"], "approval_artifact": allowlist["approval"], "configuration_artifact": allowlist["configuration"], "source_replay_artifact": allowlist["source_replay"], "allowlist_sha256": allowlist_sha, "code_artifacts": current_code, "frozen_code_artifacts": snapshots, "operational_event_reference": allowlist["operational_event_reference"], "source": dict(protocol["source"])}
    return protocol, allowlist, provenance


def check_custody_gate(root, protocol_path, allowlist_path=None, phase="training_and_conditioning_export", freeze_sha=None, executor_admission_path=None):
    try:
        protocol, allowlist, provenance = _admission(root, protocol_path, allowlist_path, phase, freeze_sha, executor_admission_path)
        ready = provenance.get("executor_ready", True)
        result = {"ready": ready, "data_ready": True, "executor_ready": ready, "scientific_scope": SCIENTIFIC_SCOPE, "measurement_arrays_opened": False, "approved_artifacts": allowlist["artifacts"], "provenance": provenance, "strong_claim_authorized": False}
        if not ready:
            result.update(reason_code="executor_admission_refresh_required", reason=provenance["refresh_requirement"])
        return result
    except (DevelopmentError, OSError, KeyError, TypeError, ValueError) as error:
        return {"ready": False, "scientific_scope": SCIENTIFIC_SCOPE, "measurement_arrays_opened": False, "reason_code": getattr(error, "code", "custody_not_approved"), "reason": str(error), "strong_claim_authorized": False}


def execute_training(root, protocol_path, allowlist_path, output_directory, executor_admission_path=None):
    protocol, allowlist, provenance = _admission(root, protocol_path, allowlist_path, "training_and_conditioning_export", executor_admission_path=executor_admission_path)
    require(provenance.get("executor_ready", True), provenance.get("refresh_requirement", "executor admission refresh required"), "executor_admission_refresh_required")
    store = ArtifactStore(root, output_directory)
    intent = store.write("training_access_intent.json", {"scientific_scope": SCIENTIFIC_SCOPE, "provenance": provenance, "artifacts": allowlist["artifacts"], "environment": environment_record(), "strong_claim_authorized": False}, "operational_event_log")
    refs = {row["role"]: row for row in allowlist["artifacts"]}
    training = read_artifact(root, refs["training_projection"], protocol, "training_projection")
    conditioning = read_artifact(root, refs["development_inputs"], protocol, "development_inputs")
    require(training["approval_sha256"] == provenance["approval_sha256"] and conditioning["approval_sha256"] == provenance["approval_sha256"], "projection payload is bound to a different custody approval")
    if provenance.get("adapter_version") == "phase_b_release_manifest_v1":
        for role, projection in (("training_projection", training), ("development_inputs", conditioning)):
            binding = provenance["source_bindings"][role]
            require(projection["source_binding_audit"] == binding["reference"], "projection source-binding reference differs from the admitted phase-B record")
            payload = {key: projection[key] for key in ("group_ids", "source_cells", "source_times")}
            require(digest(payload) == binding["record"]["ordered_payload_sha256"], "phase-B scalar/row/null payload differs from the custodian replay", "source_mismatch")
            null_counts = {key: sum(row[key] is None for row in projection["source_cells"]) for key in ("max5", "median", "relative_time_min")}
            null_counts["relative_time_min"] += sum(row["relative_time_min"] is None for row in projection["source_times"])
            require(null_counts == binding["record"]["verification"]["preserved_null_counts"], "phase-B source null preservation mismatch", "source_mismatch")
    prepared = prepare_projections(protocol, training, conditioning, refs["training_projection"]["sha256"], refs["development_inputs"]["sha256"])
    if provenance.get("adapter_version") == "phase_b_release_manifest_v1":
        require({row["group_id"]: len(row["eligible_cell_indices"]) for row in prepared["prefix_summaries"]} == provenance["approved_prefix_cells_by_group"], "derived cohorts differ from the admitted source coverage audit", "source_mismatch")
    prepared["adapter_checks"] = adapter_checks(protocol, training)
    provenance.update(training_projection_sha256=refs["training_projection"]["sha256"], development_inputs_sha256=refs["development_inputs"]["sha256"], training_projection=refs["training_projection"], development_inputs=refs["development_inputs"], training_access_intent=intent)
    fitted = fit_training(prepared, protocol)
    summary = seal_training_freeze(protocol, prepared, fitted, provenance, store)
    summary["run_summary"] = store.write("training_run_summary.json", summary, "operational_event_log")
    return summary


def execute_scoring(root, protocol_path, allowlist_path, freeze_reference):
    protocol = scoped_protocol(root, protocol_path)
    freeze = read_artifact(root, freeze_reference, protocol, "development_predictions")
    validate_freeze(protocol, freeze)
    protocol, allowlist, provenance = _admission(root, protocol_path, allowlist_path, "development_response_export", freeze_reference["sha256"])
    require(freeze["code_artifacts"] == provenance["frozen_code_artifacts"] and freeze["approval_sha256"] == provenance["approval_sha256"], "freeze code or custody approval changed before scoring")
    original_provenance = freeze["environment"]["provenance"]
    require(original_provenance["protocol_artifact"]["sha256"] == provenance["protocol_artifact"]["sha256"], "frozen protocol differs from the scoring release")
    store = ArtifactStore(root, Path(freeze_reference["path"]).parent)
    store.write("development_response_access_intent.json", {"scientific_scope": SCIENTIFIC_SCOPE, "prediction_freeze": freeze_reference, "provenance": provenance, "artifacts": allowlist["artifacts"], "strong_claim_authorized": False}, "operational_event_log")
    conditioning = read_artifact(root, original_provenance["development_inputs"], protocol, "development_inputs")
    validate_projection(protocol, conditioning, "development_inputs")
    response_ref = allowlist["artifacts"][0]
    response = read_artifact(root, response_ref, protocol, "development_responses")
    validate_projection(protocol, response, "development_responses")
    require(response["approval_sha256"] == freeze["approval_sha256"], "response projection has a different approval")
    prefix = {row["group_id"]: row for row in freeze["prefix_summaries"]}
    frames = []
    for group in group_specs(protocol)[3:]:
        cells = [row for row in response["source_cells"] if row["group_id"] == group["group_id"]]
        times = [row for row in conditioning["source_times"] if row["group_id"] == group["group_id"]]
        frames.extend(derive_frames(group, cells, times, prefix[group["group_id"]]))
    report = build_development_report(protocol, freeze, frames, freeze_reference["sha256"], response_ref["sha256"])
    report_ref = store.write("development_report.json", report, "development_report")
    aggregation_ref = store.write("development_aggregation_audit.json", {"scientific_scope": SCIENTIFIC_SCOPE, "frames": frames, "fixed_prefix_summaries": freeze["prefix_summaries"][3:], "strong_claim_authorized": False}, "development_report")
    models = {}
    for model_id in MODEL_IDS:
        groups = [row for row in report["group_evaluations"] if row["model_id"] == model_id]
        complete = all(row["complete"] for row in groups)
        mse = float(np.mean([row["mean_mse"] for row in groups])) if complete else None
        mae = float(np.mean([row["mae"] for row in groups])) if complete else None
        rmse = math.sqrt(mse) if complete else None
        models[model_id] = {"complete": complete, "mean_mse": mse, "rmse": rmse, "mae": mae, "scaled_rmse": rmse / freeze["S_train"] if complete and freeze["S_train"] > 0.0 else None, "scaled_mae": mae / freeze["S_train"] if complete and freeze["S_train"] > 0.0 else None}
    summary = {"scientific_scope": SCIENTIFIC_SCOPE, "status": report["status"], "prediction_freeze": freeze_reference, "development_responses": response_ref, "report": report_ref, "aggregation_audit": aggregation_ref, "models": models, "selected_development_model": report["selected_development_model"], "limitations": list(LIMITATIONS), "strong_claim_authorized": False, "final_test_defined": False}
    summary["summary_artifact"] = store.write("development_summary.json", summary, "development_report")
    return summary
