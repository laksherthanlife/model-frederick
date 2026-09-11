from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import multiprocessing
import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.analysis.hybrid_stress import fit_sensor_model
from ystwin.analysis.hybrid_validation import score_external_targets, summarize_scores
from ystwin.fba.dynamic_rates import build_network_rates
from ystwin.fba.product_panel import (
    biomass_pool_quota,
    install_product,
    prepare_panel_model,
    product_names,
)
from ystwin.fba.solver import load_model
from ystwin.hybrid import Environment, ScenarioUnsupported, simulate_hybrid_case


ROOT = paths.REPO_ROOT
_WORKER_SENSOR = None
_WORKER_PREPARED = None
_WORKER_ARGS = None


def _initialize_worker(sensor, model_path, args):
    global _WORKER_SENSOR, _WORKER_PREPARED, _WORKER_ARGS
    base, _ = load_model(model_path)
    _WORKER_PREPARED, _ = prepare_panel_model(base)
    _WORKER_SENSOR, _WORKER_ARGS = sensor, args


def _execute_job(job, sensor, prepared, args):
    if job["kind"] == "external":
        frame = _predict_external(prepared, pd.DataFrame([job["target"]]), args)
        return frame.iloc[0].to_dict(), None
    return _case(sensor, prepared, job["product"], job["environment"], args,
                 **job.get("options", {}))


def _worker_job(job):
    return _execute_job(job, _WORKER_SENSOR, _WORKER_PREPARED, _WORKER_ARGS)


def _job_results(jobs, sensor, prepared, model_path, args):
    if args.workers == 1:
        for job in jobs:
            yield _execute_job(job, sensor, prepared, args)
        return
    with ProcessPoolExecutor(
            max_workers=args.workers, mp_context=multiprocessing.get_context("spawn"),
            initializer=_initialize_worker, initargs=(sensor, model_path, args)) as pool:
        yield from pool.map(_worker_job, jobs, chunksize=1)


def environments(hours: float) -> tuple[Environment, ...]:
    baseline = Environment("glucose_aerobic", hours=hours)
    return (
        baseline,
        replace(baseline, name="glucose_low_oxygen", oxygen_uptake=2.0),
        replace(baseline, name="ethanol_aerobic", carbon_source="ethanol", substrate_g_per_l=15.343),
        replace(baseline, name="glucose_37C", temperature_c=37.0),
        replace(baseline, name="glucose_DTT_1mM", stressor="DTT", dose_mM=1.0),
        replace(baseline, name="glucose_H2O2_0p5mM", stressor="H2O2", dose_mM=0.5),
    )


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _write_json(path: pathlib.Path, value) -> None:
    path.write_text(json.dumps(_json_value(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def _hash(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _case(sensor, prepared, product, environment, args, **options):
    installed, task = install_product(prepared, product)
    try:
        result = simulate_hybrid_case(
            sensor, installed, task, environment,
            allocation_fraction=options.pop("allocation_fraction", args.allocation),
            growth_fraction=options.pop("growth_fraction", args.growth_fraction),
            grid_points=args.grid_points, steps=args.steps, **options)
    except ScenarioUnsupported as exc:
        return {"product": product, "environment": environment.name,
                "status": "unsupported", "reason": str(exc)}, None
    except (ValueError, RuntimeError) as exc:
        return {"product": product, "environment": environment.name,
                "status": "failed", "failure_type": type(exc).__name__, "reason": str(exc)}, None
    return result.summary(), result


def _predict_external(prepared, targets: pd.DataFrame, args) -> pd.DataFrame:
    rows, cache = [], {}
    input_columns = [name for name in targets.columns if name not in (
        "observed", "uncertainty", "uncertainty_kind", "observation_kind", "detection_limit")]
    for target in targets[input_columns].to_dict("records"):
        row = {key: target[key] for key in ("study_id", "case_id", "product", "metric")}
        row.update(predicted=np.nan, prediction_units=target["units"], comparable=False,
                   prediction_status="unsupported", prediction_basis="unavailable", reason="")
        if target["study_id"] not in ("hakkaart2020", "vanhoek1998"):
            row["reason"] = (
                "Matched initial dry biomass, product pool, uptake/oxygen-transfer and genotype "
                "response are unavailable; no synthetic initial values substituted")
            rows.append(row)
            continue
        if target["product"] == "glycogen":
            row["reason"] = (
                "Model glucose-equivalent glycogen mass and source hydrolysis-assay mass "
                "have an unresolved conversion; not declared comparable")
            rows.append(row)
            continue
        oxygen = target["oxygen_uptake"]
        if not np.isfinite(oxygen):
            row["reason"] = "Matched oxygen uptake was not determined; not replaced by a fitted oxygen bound"
            rows.append(row)
            continue
        uptake, mu = float(target["glucose_uptake"]), float(target["growth_rate_per_h"])
        key = (target["product"], uptake, mu, float(oxygen))
        try:
            if key not in cache:
                model, task = install_product(prepared, target["product"])
                rates, info = build_network_rates(
                    model, "r_1714", "r_1992", "r_4046",
                    product_reaction_id=task.reaction_id,
                    uptake_vmax_mmol_per_gdcw_h=uptake,
                    oxygen_lower_bound=-float(oxygen), fixed_growth_per_h=mu,
                    allocation_fraction=args.allocation, growth_fraction=args.growth_fraction,
                    grid_points=args.grid_points)
                solved = rates.at_uptake(uptake)
                quota = (biomass_pool_quota(model, target["product"])
                         if target["product"] == "trehalose" else None)
                cache[key] = solved, task, quota, info
            (_, actual_uptake, product_rate), task, quota, info = cache[key]
            if target["product"] == "glycerol":
                prediction = product_rate
                basis = (
                    "physiology-conditioned net secretion; measured growth and uptake/oxygen caps; "
                    "industrial-strain transfer; no measured reporter state")
            else:
                dilution = float(target["dilution_rate_per_h"])
                viable = float(target["viable_fraction"])
                if dilution <= 0 or not 0 < viable <= 1:
                    raise ValueError("Valid dilution and viable fraction are required for total-pool accounting")
                prediction = (quota + product_rate * viable / dilution) * task.molar_mass_g_per_mol
                basis = (
                    "physiology-conditioned total dry-biomass content; frozen GEM basal quota plus "
                    "extra net synthesis times viable fraction divided by dilution; "
                    "assumes equal dry mass per PI-classified cell and no extra-pool loss; "
                    "not an independent prediction of pH effects or latent stress")
            row.update(predicted=prediction, prediction_status="predicted", comparable=True,
                       prediction_basis=basis, reason="",
                       actual_substrate_uptake=actual_uptake,
                       predicted_specific_rate=product_rate,
                       allocation_fraction=args.allocation,
                       product_targets_used_for_fitting=False)
        except ValueError as exc:
            row.update(prediction_status="infeasible", reason=f"{type(exc).__name__}: {exc}")
        except RuntimeError as exc:
            row.update(prediction_status="invalid", reason=f"{type(exc).__name__}: {exc}")
        rows.append(row)
    return pd.DataFrame(rows)


def _external_baselines(prepared, targets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    quota = biomass_pool_quota(prepared, "trehalose")
    _, task = install_product(prepared, "trehalose")
    for target in targets.to_dict("records"):
        row = {key: target[key] for key in ("study_id", "case_id", "product", "metric")}
        row.update(predicted=np.nan, prediction_units=target["units"], comparable=False,
                   prediction_status="unsupported", prediction_basis="unavailable",
                   reason="No declared baseline for this target")
        if target["product"] == "glycerol" and target["metric"] == "specific_production_rate":
            row.update(predicted=0.0, comparable=True, prediction_status="predicted",
                       prediction_basis="zero-product null, not a fitted or metabolic prediction", reason="")
        elif target["product"] == "trehalose" and target["metric"] == "intracellular_content":
            row.update(predicted=quota * task.molar_mass_g_per_mol, comparable=True,
                       prediction_status="predicted", prediction_basis="frozen GEM biomass quota only",
                       reason="")
        rows.append(row)
    return pd.DataFrame(rows)


def _compare_numerics(frame, protocol, reference_dir, destination):
    reference = json.loads((reference_dir / "protocol.json").read_text())
    current = _json_value(protocol)
    for key in ("input_hashes", "training_plates", "heldout_plate", "allocation_fraction",
                "growth_floor_fraction"):
        if reference[key] != current[key]:
            raise ValueError(f"numerical comparison changes {key}, not just numerical resolution")
    previous_seed = json.loads((reference_dir / "sensor_model.json").read_text())["seed"]
    if previous_seed != current["training_seed"]:
        raise ValueError("numerical comparison changes the sensor-training seed")
    previous_cases = {case["name"]: case for case in reference["environments"]}
    if any(previous_cases.get(case["name"]) != case for case in current["environments"]):
        raise ValueError("numerical comparison changes environment inputs")
    for source, fingerprint in current["code_hashes"].items():
        if source != "scripts/run_hybrid_benchmark.py" and reference["code_hashes"].get(source) != fingerprint:
            raise ValueError(f"numerical comparison changes model implementation: {source}")
    keys = ["product", "environment"]
    metrics = ["titre_mg_per_l", "mean_specific_rate_mmol_per_gdcw_h", "yield_g_per_g"]
    previous = pd.read_csv(reference_dir / "scenarios.csv")
    joined = frame[keys + ["status"] + metrics].merge(
        previous[keys + ["status"] + metrics], on=keys, validate="one_to_one",
        suffixes=("", "_reference"))
    if len(joined) != len(frame) or not (
            joined.status.eq("predicted") & joined.status_reference.eq("predicted")).all():
        raise ValueError("numerical comparison requires matching successful scenarios")
    for metric in metrics:
        before = joined[f"{metric}_reference"].to_numpy()
        after = joined[metric].to_numpy()
        change = np.full(len(joined), np.nan)
        np.divide(100.0 * (after - before), before, out=change, where=before != 0)
        change[(before == 0) & (after == 0)] = 0.0
        joined[f"{metric}_change_pct"] = change
        joined[f"{metric}_reference_zero_changed"] = (before == 0) & (after != 0)
    joined.to_csv(destination, index=False)
    print("\nNumerical refinement, percentage changes from the frozen reference:", flush=True)
    print(joined[keys + [f"{metric}_change_pct" for metric in metrics]].to_string(index=False), flush=True)


_ASSESS_TARGET_KEYS = ["study_id", "case_id", "product", "metric"]
_ASSESS_SCENARIO_KEYS = ["product", "environment"]
_ASSESS_SCORE_STATUSES = ("scored", "censored", "infeasible", "unsupported", "invalid", "missing")
_ASSESS_SCENARIO_STATUSES = ("predicted", "unsupported", "failed", "infeasible", "invalid", "pending", "running")
_ASSESS_ERRORS = ["abs_error", "relative_error_pct", "absolute_percent_error", "fold_error"]
_ASSESS_NUMERICAL_METRICS = ["titre_mg_per_l", "mean_specific_rate_mmol_per_gdcw_h", "yield_g_per_g"]
_ASSESS_SUMMARY_COUNTS = ["n_total", "n_scored", "n_censored", "n_unsupported", "n_infeasible",
                          "n_missing", "n_invalid", "n_percent_scored", "n_fold_scored"]
_ASSESS_SUMMARY_METRICS = ["mae", "median_abs_percent_error", "median_fold_error"]


class AssessmentError(ValueError):
    pass


class NumericalCompatibilityError(AssessmentError):
    pass


def _assessment_artifact(directory, name, artifacts, prefix=""):
    path = directory / name
    try:
        raw = path.read_bytes()
    except FileNotFoundError:
        artifacts[prefix + name] = {"status": "missing"}
        return None
    except OSError as exc:
        raise AssessmentError(f"Cannot read {path}: {exc}") from exc
    entry = {"status": "present", "sha256": hashlib.sha256(raw).hexdigest()}
    artifacts[prefix + name] = entry
    if not raw.strip():
        entry["status"] = "empty"
        return None
    try:
        text = raw.decode("utf-8-sig")
        if path.suffix == ".json":
            def pairs(items):
                result = {}
                for key, value in items:
                    if key in result:
                        raise AssessmentError(f"{path}: duplicate JSON key {key}")
                    result[key] = value
                return result

            value = json.loads(text, object_pairs_hook=pairs)
            try:
                json.dumps(value, allow_nan=False)
            except ValueError as exc:
                raise AssessmentError(f"{path}: non-finite JSON number") from exc
            if not isinstance(value, dict):
                raise AssessmentError(f"{path}: expected a JSON object")
            return value
        records = [row for row in csv.reader(io.StringIO(text), strict=True) if row]
        if not records:
            entry["status"] = "empty"
            return None
        if len(set(records[0])) != len(records[0]) or any(not key.strip() for key in records[0]):
            raise AssessmentError(f"{path}: duplicate or empty CSV columns")
        if any(len(row) != len(records[0]) for row in records[1:]):
            raise AssessmentError(f"{path}: incomplete or malformed CSV row")
        frame = pd.read_csv(io.StringIO(text), float_precision="round_trip")
        if frame.empty:
            entry["status"] = "empty"
        return frame
    except (UnicodeError, json.JSONDecodeError, csv.Error, pd.errors.ParserError) as exc:
        raise AssessmentError(f"Malformed artifact {path}: {exc}") from exc


def _assessment_table(frame, required, keys, name):
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise AssessmentError(f"{name}: missing columns {missing}")
    for key in keys:
        if not all(isinstance(value, str) and value.strip() for value in frame[key]):
            raise AssessmentError(f"{name}: {key} must contain nonempty string keys")
    if frame.duplicated(keys).any():
        raise AssessmentError(f"{name}: duplicate keys {keys}")


def _assessment_count(value, name):
    if not _assessment_number(value) or value < 0 or int(value) != value:
        raise AssessmentError(f"{name}: expected a nonnegative integer")
    return int(value)


def _assessment_number(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
        return False
    try:
        return math.isfinite(value)
    except (OverflowError, TypeError, ValueError):
        return False


def _assessment_csv_number(value):
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value
    return value


def _assessment_equal(left, right):
    if pd.isna(left) or pd.isna(right):
        return bool(pd.isna(left) and pd.isna(right))
    if isinstance(left, (bool, np.bool_)) or isinstance(right, (bool, np.bool_)):
        return isinstance(left, (bool, np.bool_)) and isinstance(right, (bool, np.bool_)) and left == right
    if _assessment_number(left) and _assessment_number(right):
        return math.isclose(float(left), float(right), rel_tol=1e-9, abs_tol=1e-12)
    return left == right


def _assessment_protocol(protocol, name):
    if protocol is None:
        return
    if "run_complete" in protocol and not isinstance(protocol["run_complete"], bool):
        raise AssessmentError(f"{name}: run_complete must be a boolean, not a science verdict")
    for field in ("scenario_rows", "successful_scenarios", "planned_jobs", "training_seed", "grid_points", "steps"):
        if field in protocol:
            _assessment_count(protocol[field], f"{name}: {field}")
    for field in ("products", "training_plates"):
        if field in protocol:
            values = protocol[field]
            if (not isinstance(values, list) or not all(isinstance(item, str) and item.strip() for item in values)
                    or len(set(values)) != len(values)):
                raise AssessmentError(f"{name}: {field} must contain distinct nonempty strings")
    if "environments" in protocol:
        cases = protocol["environments"]
        if (not isinstance(cases, list) or any(not isinstance(case, dict)
                or not isinstance(case.get("name"), str) or not case["name"].strip() for case in cases)
                or len({case["name"] for case in cases}) != len(cases)):
            raise AssessmentError(f"{name}: environments must have distinct nonempty names")
    for field in ("input_hashes", "code_hashes"):
        if field in protocol:
            values = protocol[field]
            if not isinstance(values, dict) or not all(isinstance(value, str) and value for value in values.values()):
                raise AssessmentError(f"{name}: {field} must be a mapping of fingerprints")
    for field in ("allocation_fraction", "growth_floor_fraction"):
        if field in protocol:
            value = protocol[field]
            if (not _assessment_number(value) or not 0 <= value <= 1
                    or (field == "growth_floor_fraction" and value == 0)):
                raise AssessmentError(f"{name}: invalid {field}")
    if "product_targets_used_for_fitting" in protocol and not isinstance(protocol["product_targets_used_for_fitting"], bool):
        raise AssessmentError(f"{name}: product_targets_used_for_fitting must be boolean")


def _assessment_scenarios(frame, name):
    if frame is None:
        return pd.DataFrame(columns=[*_ASSESS_SCENARIO_KEYS, "status"])
    _assessment_table(frame, [*_ASSESS_SCENARIO_KEYS, "status"], _ASSESS_SCENARIO_KEYS, name)
    if not frame.status.isin(_ASSESS_SCENARIO_STATUSES).all():
        raise AssessmentError(f"{name}: unknown scenario status")
    return frame


def _assessment_execution(protocol, scenarios, auxiliary, name):
    _assessment_protocol(protocol, f"{name}/protocol.json")
    frame = _assessment_scenarios(scenarios, f"{name}/scenarios.csv")
    declaration = protocol or {}
    products = declaration.get("products")
    environments_declared = declaration.get("environments")
    plan = (None if products is None or environments_declared is None else
            {(product, case["name"]) for product in products for case in environments_declared})
    actual = set(frame[_ASSESS_SCENARIO_KEYS].itertuples(index=False, name=None))
    if plan is not None and actual - plan:
        raise AssessmentError(f"{name}/scenarios.csv: scenarios outside the protocol plan: {sorted(actual - plan)}")
    counts = {status: int(frame.status.eq(status).sum()) for status in _ASSESS_SCENARIO_STATUSES}
    for field, count in (("scenario_rows", len(frame)), ("successful_scenarios", counts["predicted"])):
        if scenarios is not None and field in declaration and declaration[field] != count:
            raise AssessmentError(f"{name}/protocol.json: {field} contradicts scenarios.csv")
    missing = None if plan is None else sorted(plan - actual)
    actual_jobs = len(frame)
    for filename, jobs in auxiliary.items():
        if jobs is not None:
            required = (_ASSESS_TARGET_KEYS + ["prediction_status"] if filename == "external_predictions.csv"
                        else [*_ASSESS_SCENARIO_KEYS, "status"])
            if not set(required).issubset(jobs.columns):
                raise AssessmentError(f"{name}/{filename}: missing job identity or status columns")
            actual_jobs += len(jobs)
    planned_jobs = declaration.get("planned_jobs")
    if planned_jobs is not None and actual_jobs > planned_jobs:
        raise AssessmentError(f"{name}/protocol.json: more recorded jobs than planned_jobs")
    census_complete = plan is not None and bool(plan) and not missing
    jobs_complete = planned_jobs is None or planned_jobs == actual_jobs
    flag = declaration.get("run_complete")
    complete = bool(flag is True and census_complete and jobs_complete
                    and not counts["pending"] and not counts["running"])
    if protocol is None and scenarios is None:
        status = "missing"
    elif plan is None or flag is None:
        status = "unverified"
    else:
        status = "complete" if complete else "incomplete"
    return {
        "status": status, "run_complete": complete, "protocol_run_complete": flag,
        "planned_scenarios": None if plan is None else len(plan), "actual_scenarios": len(frame),
        "scenario_status_counts": counts,
        "missing_scenarios": None if missing is None else [dict(zip(_ASSESS_SCENARIO_KEYS, key)) for key in missing],
        "planned_jobs": planned_jobs, "actual_jobs": actual_jobs,
        "job_census_status": "not_recorded" if planned_jobs is None else ("complete" if jobs_complete else "incomplete"),
        "per_product": {product: {
            "planned": None if plan is None else sum(key[0] == product for key in plan),
            "actual": int(frame["product"].eq(product).sum()),
            "predicted": int((frame["product"].eq(product) & frame.status.eq("predicted")).sum()),
        } for product in sorted(set(products or []) | set(frame["product"]))},
        "meaning": "Execution bookkeeping and scenario census only; not scientific model validation.",
    }


def _assessment_scores(frame, name):
    if frame is None:
        return None
    required = [*_ASSESS_TARGET_KEYS, "units", "observed", "observation_kind", "detection_limit",
                "predicted", "prediction_status", "prediction_units", "comparable", "score_status", *_ASSESS_ERRORS]
    _assessment_table(frame, required, _ASSESS_TARGET_KEYS, name)
    if not frame.score_status.isin(_ASSESS_SCORE_STATUSES).all():
        raise AssessmentError(f"{name}: unknown score_status")
    if not all(isinstance(value, str) and value.strip() for value in frame.units):
        raise AssessmentError(f"{name}: missing measurement units")
    for _, group in frame.groupby(["product", "metric"]):
        if group.units.nunique() != 1:
            raise AssessmentError(f"{name}: mixed units for a product/metric")
    frame = frame.copy()
    for field in ("observed", "predicted", "detection_limit"):
        frame[field] = frame[field].map(_assessment_csv_number)
    for field in _ASSESS_ERRORS:
        try:
            if frame[field].map(lambda value: isinstance(value, (bool, np.bool_))).any():
                raise ValueError("boolean error")
            frame[field] = pd.to_numeric(frame[field], errors="raise")
        except (TypeError, ValueError) as exc:
            raise AssessmentError(f"{name}: {field} must be numeric or missing") from exc
    for row in frame.to_dict("records"):
        scored = row["score_status"] == "scored"
        censored = row["observation_kind"] == "below_detection_limit"
        if "is_censored" in row and (not isinstance(row["is_censored"], (bool, np.bool_))
                                     or row["is_censored"] != censored):
            raise AssessmentError(f"{name}: is_censored contradicts observation_kind")
        if not scored:
            if any(not pd.isna(row[field]) for field in _ASSESS_ERRORS):
                raise AssessmentError(f"{name}: unscored targets must not contain fabricated errors")
            if row["score_status"] == "censored" and not censored:
                raise AssessmentError(f"{name}: censored score contradicts observation_kind")
            continue
        if (row["observation_kind"] != "quantified" or row["prediction_status"] != "predicted"
                or not isinstance(row["comparable"], (bool, np.bool_)) or not row["comparable"]
                or row["prediction_units"] != row["units"]):
            raise AssessmentError(f"{name}: scored row is not a comparable quantified prediction")
        observed, predicted = row["observed"], row["predicted"]
        if not all(_assessment_number(value) and value >= 0 for value in (observed, predicted)):
            raise AssessmentError(f"{name}: scored observed/predicted values must be finite and nonnegative")
        expected = {"abs_error": abs(predicted - observed)}
        if observed > 0:
            expected["relative_error_pct"] = (predicted - observed) / observed * 100.0
            expected["absolute_percent_error"] = abs(expected["relative_error_pct"])
            if predicted > 0:
                expected["fold_error"] = max(predicted / observed, observed / predicted)
        for field in _ASSESS_ERRORS:
            value = row[field]
            if field not in expected:
                if not pd.isna(value):
                    raise AssessmentError(f"{name}: {field} is undefined for this quantified zero")
            elif _assessment_number(value) and not _assessment_equal(value, expected[field]):
                raise AssessmentError(f"{name}: {field} contradicts observed/predicted values")
    return frame


def _assessment_summary(summary, scores, name):
    if summary is not None:
        _assessment_table(summary, ["product", "metric", "units", *_ASSESS_SUMMARY_COUNTS,
                                   *_ASSESS_SUMMARY_METRICS], ["product", "metric"], name)
        for row in summary.to_dict("records"):
            for field in _ASSESS_SUMMARY_COUNTS:
                _assessment_count(row[field], f"{name}: {field}")
            censored_status = row["n_total"] - sum(row[f"n_{status}"] for status in _ASSESS_SCORE_STATUSES if status != "censored")
            if (not 0 <= censored_status <= row["n_censored"] <= row["n_total"]
                    or row["n_scored"] + row["n_censored"] > row["n_total"]
                    or not 0 <= row["n_fold_scored"] <= row["n_percent_scored"] <= row["n_scored"]):
                raise AssessmentError(f"{name}: contradictory coverage counts")
            for field, count in zip(_ASSESS_SUMMARY_METRICS, (row["n_scored"], row["n_percent_scored"], row["n_fold_scored"])):
                value = row[field]
                if not pd.isna(value) and (not isinstance(value, (int, float, np.number)) or value < 0 or count == 0):
                    raise AssessmentError(f"{name}: {field} contradicts its scored count")
    if scores is None:
        return summary
    calculated = summarize_scores(scores)
    if summary is not None:
        saved = {(row["product"], row["metric"]): row for row in summary.to_dict("records")}
        expected = {(row["product"], row["metric"]): row for row in calculated.to_dict("records")}
        if saved.keys() != expected.keys() or any(
                not _assessment_equal(saved[key][field], row[field]) for key, row in expected.items()
                for field in ["units", *_ASSESS_SUMMARY_COUNTS, *_ASSESS_SUMMARY_METRICS]):
            raise AssessmentError(f"{name}: summary contradicts score rows")
    return calculated


def _assessment_score_counts(scores, summary):
    if scores is not None:
        return {status: int(scores.score_status.eq(status).sum()) for status in _ASSESS_SCORE_STATUSES}
    if summary is None:
        return dict.fromkeys(_ASSESS_SCORE_STATUSES)
    counts = {status: int(summary[f"n_{status}"].sum()) for status in _ASSESS_SCORE_STATUSES if status != "censored"}
    counts["censored"] = int(summary.n_total.sum()) - sum(counts.values())
    return counts


def _assessment_observation_counts(scores, summary):
    if scores is None:
        return {"quantified": None, "below_detection_limit": None if summary is None else int(summary.n_censored.sum()),
                "other": None, "unknown_detection_limit": None, "known_detection_limit": None,
                "invalid_detection_limit": None}
    censored = scores.loc[scores.observation_kind.eq("below_detection_limit"), "detection_limit"]
    known = sum(_assessment_number(value) and value > 0 for value in censored)
    unknown = int(censored.isna().sum())
    return {"quantified": int(scores.observation_kind.eq("quantified").sum()),
            "below_detection_limit": len(censored),
            "other": int((~scores.observation_kind.isin(["quantified", "below_detection_limit"])).sum()),
            "unknown_detection_limit": unknown, "known_detection_limit": int(known),
            "invalid_detection_limit": len(censored) - int(known) - unknown}


def _assessment_nonfinite_errors(scores):
    if scores is None:
        return None
    count = 0
    for row in scores.loc[scores.score_status.eq("scored")].to_dict("records"):
        fields = ["abs_error"]
        if row["observed"] > 0:
            fields.extend(["relative_error_pct", "absolute_percent_error"])
            if row["predicted"] > 0:
                fields.append("fold_error")
        count += any(not _assessment_number(row[field]) for field in fields)
    return count


def _assessment_row_map(frame, keys):
    return {} if frame is None else {tuple(row[key] for key in keys): row for row in frame.to_dict("records")}


def _assessment_target_identity(scores, baselines):
    model = _assessment_row_map(scores, _ASSESS_TARGET_KEYS)
    baseline = _assessment_row_map(baselines, _ASSESS_TARGET_KEYS)
    for key in model.keys() & baseline.keys():
        for field in ("observed", "units", "observation_kind", "detection_limit"):
            if not _assessment_equal(model[key][field], baseline[key][field]):
                raise AssessmentError(f"external_baseline_scores.csv: {field} differs for target {key}")


def _assessment_median(values):
    if not values or not all(_assessment_number(value) for value in values):
        return None
    ordered = sorted(float(value) for value in values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else ordered[middle - 1] / 2 + ordered[middle] / 2


def _assessment_comparison(model, baseline, metric, units):
    def scored_rows(frame):
        if frame is None:
            return {}, None
        selected = frame.loc[frame.score_status.eq("scored") & frame.metric.eq(metric)]
        eligible = selected.loc[selected.abs_error.map(_assessment_number)]
        return _assessment_row_map(eligible, _ASSESS_TARGET_KEYS), len(selected)

    left, left_count = scored_rows(model)
    right, right_count = scored_rows(baseline)
    matched = sorted(left.keys() & right.keys())
    report = {"metric": metric, "units": units, "model_scored_rows": left_count,
              "baseline_scored_rows": right_count, "matched_scored_rows": len(matched),
              "model_only_scored_rows": len(left.keys() - right.keys()),
              "baseline_only_scored_rows": len(right.keys() - left.keys()),
              "matched_keys": [dict(zip(_ASSESS_TARGET_KEYS, key)) for key in matched]}
    for name, rows in (("model", left), ("baseline", right)):
        values = [rows[key] for key in matched]
        report[f"{name}_mae"] = math.fsum(row["abs_error"] / len(values) for row in values) if values else None
        report[f"{name}_median_abs_percent_error"] = _assessment_median([
            row["absolute_percent_error"] for row in values if row["observed"] > 0])
        report[f"{name}_median_fold_error"] = _assessment_median([
            row["fold_error"] for row in values if row["observed"] > 0 and row["predicted"] > 0])
    report["result"] = ("insufficient_coverage" if not matched else
                        "beats_baseline_on_scored_rows" if report["model_mae"] < report["baseline_mae"] else
                        "does_not_beat_baseline")
    return report


def _assessment_comparison_result(comparisons):
    results = {row["result"] for row in comparisons if row["matched_scored_rows"]}
    if not results:
        return "insufficient_coverage"
    return results.pop() if len(results) == 1 else "mixed_against_baseline"


def _assessment_external_report(directory, protocol, artifacts):
    frames, summaries = {}, {}
    for prefix in ("external", "external_baseline"):
        frames[prefix] = _assessment_scores(
            _assessment_artifact(directory, f"{prefix}_scores.csv", artifacts), f"{prefix}_scores.csv")
        summaries[prefix] = _assessment_summary(
            _assessment_artifact(directory, f"{prefix}_summary.csv", artifacts), frames[prefix], f"{prefix}_summary.csv")
    scores, baseline = frames["external"], frames["external_baseline"]
    summary, baseline_summary = summaries["external"], summaries["external_baseline"]
    _assessment_target_identity(scores, baseline)
    planned_products = (protocol or {}).get("products", [])
    products = set(planned_products)
    for frame in (summary, baseline_summary):
        if frame is not None:
            products.update(frame["product"])
    per_product = {}
    all_compared, complete_coverage = [], bool(planned_products) and scores is not None and baseline is not None
    for product in sorted(products):
        model_rows, baseline_rows, model_summary, null_summary = [
            None if frame is None else frame.loc[frame["product"].eq(product)]
            for frame in (scores, baseline, summary, baseline_summary)]
        target_rows = len(model_rows) if model_rows is not None else (None if model_summary is None else int(model_summary.n_total.sum()))
        counts = _assessment_score_counts(model_rows, model_summary)
        observations = _assessment_observation_counts(model_rows, model_summary)
        nonfinite = _assessment_nonfinite_errors(model_rows)
        metrics = {} if model_summary is None else {row["metric"]: row["units"] for row in model_summary.to_dict("records")}
        if null_summary is not None:
            for row in null_summary.to_dict("records"):
                metrics.setdefault(row["metric"], row["units"])
        comparisons = [_assessment_comparison(model_rows, baseline_rows, metric, units)
                       for metric, units in sorted(metrics.items())]
        all_compared.extend(comparisons)
        if target_rows is None:
            status = "missing"
        elif target_rows == 0:
            status = "no_targets"
        elif model_rows is None:
            status = "summary_only"
        elif observations["quantified"] == 0:
            status = "no_quantified_targets"
        elif not counts["scored"]:
            status = "no_scored_quantified_targets"
        elif counts["scored"] != observations["quantified"] or nonfinite:
            status = "partial_quantified_coverage"
        else:
            status = "quantified_targets_scored"
        per_product[product] = {
            "status": status, "requested_product": product in planned_products, "target_rows": target_rows,
            "score_status_counts": counts, "observation_counts": observations,
            "nonfinite_scored_error_rows": nonfinite,
            "model_summary": [] if model_summary is None else model_summary.to_dict("records"),
            "baseline_comparisons": comparisons, "comparison_result": _assessment_comparison_result(comparisons),
        }
        if product in planned_products:
            matched = sum(row["matched_scored_rows"] for row in comparisons)
            if (not target_rows or matched != target_rows or nonfinite
                    or baseline_rows is None or len(baseline_rows) != target_rows):
                complete_coverage = False
    return {
        "status": ("assessed" if scores is not None and len(scores) else "empty" if scores is not None else
                   "summary_only" if summary is not None else "missing"),
        "coverage_status": "complete_on_recorded_targets" if complete_coverage else "incomplete",
        "target_rows": len(scores) if scores is not None else (None if summary is None else int(summary.n_total.sum())),
        "score_status_counts": _assessment_score_counts(scores, summary),
        "observation_counts": _assessment_observation_counts(scores, summary),
        "nonfinite_scored_error_rows": _assessment_nonfinite_errors(scores), "per_product": per_product,
        "comparison_result": _assessment_comparison_result(all_compared),
        "summary_consistency": ("verified" if scores is not None and artifacts["external_summary.csv"]["status"] == "present"
                                else "unverified" if scores is None else artifacts["external_summary.csv"]["status"]),
        "comparison_basis": "MAE on the intersection of finite-error, quantified, scored target keys; never unmatched summary averages.",
        "censoring_basis": "Score statuses partition rows; censored observations overlap infeasible/unsupported/invalid statuses. Unknown LOD is not zero.",
    }


def _assessment_seed(protocol, sensor, label):
    declared = (protocol or {}).get("training_seed")
    recorded = (sensor or {}).get("seed")
    for value in (declared, recorded):
        if value is not None:
            _assessment_count(value, f"numerical {label} sensor seed")
    if declared is not None and recorded is not None and declared != recorded:
        raise NumericalCompatibilityError(f"numerical {label}: training_seed contradicts sensor_model.json seed")
    return declared if declared is not None else recorded


def _assessment_numerical_compatibility(reference_dir, numerical_dir, reference, current, sensor, refined_sensor):
    missing, checked = [], []
    reference, current = reference or {}, current or {}
    link = current.get("numerical_reference")
    if link is None:
        missing.append("numerical_reference")
    elif not isinstance(link, str) or not link.strip():
        raise NumericalCompatibilityError("numerical_reference must identify the assessed reference directory")
    else:
        path = pathlib.Path(link)
        candidates = [path] if path.is_absolute() else [ROOT / path, numerical_dir / path, numerical_dir.parent / path]
        if reference_dir not in {candidate.resolve() for candidate in candidates}:
            raise NumericalCompatibilityError("numerical_reference does not link to the assessed run")
        checked.append("numerical_reference")
    for field in ("input_hashes", "training_plates", "heldout_plate", "allocation_fraction",
                  "growth_floor_fraction", "product_targets_used_for_fitting"):
        if field not in reference or field not in current or reference[field] is None or current[field] is None:
            missing.append(field)
        elif reference[field] != current[field]:
            raise NumericalCompatibilityError(f"numerical comparison changes {field}")
        elif field in ("input_hashes", "training_plates", "heldout_plate") and not reference[field]:
            missing.append(field)
        else:
            checked.append(field)
    seeds = [_assessment_seed(reference, sensor, "reference"), _assessment_seed(current, refined_sensor, "refined")]
    if None in seeds:
        missing.append("training_seed")
    elif seeds[0] != seeds[1]:
        raise NumericalCompatibilityError("numerical comparison changes sensor-training seed")
    else:
        checked.append("training_seed")
    if sensor is not None and refined_sensor is not None:
        if sensor != refined_sensor:
            raise NumericalCompatibilityError("numerical comparison changes the fitted sensor metadata or assumptions")
        checked.append("sensor_model")
    for field in ("versions", "solver", "allocation_status", "sensitivity_status", "scope", "external_targets_sha256"):
        if field in reference and field in current:
            if reference[field] != current[field]:
                raise NumericalCompatibilityError(f"numerical comparison changes {field}")
            checked.append(field)
        elif field in ("versions", "solver", "allocation_status", "scope") and (field in reference or field in current):
            missing.append(field)
    fingerprints = [{key: value for key, value in protocol.get("code_hashes", {}).items()
                     if key != "scripts/run_hybrid_benchmark.py"} for protocol in (reference, current)]
    if not all(fingerprints):
        missing.append("model_code_hashes")
    elif fingerprints[0] != fingerprints[1]:
        raise NumericalCompatibilityError("numerical comparison changes model implementation fingerprints")
    else:
        checked.append("model_code_hashes_except_benchmark_driver")
    previous_cases = {case["name"]: case for case in reference.get("environments", [])}
    current_cases = {case["name"]: case for case in current.get("environments", [])}
    if not previous_cases or not current_cases:
        missing.append("environments")
    else:
        for name, case in current_cases.items():
            if name not in previous_cases or case != previous_cases[name]:
                raise NumericalCompatibilityError(f"numerical comparison changes environment inputs: {name}")
            if not set(Environment.__dataclass_fields__).issubset(case):
                missing.append(f"environment_inputs:{name}")
        checked.append("environments")
    if not reference.get("products") or not current.get("products"):
        missing.append("products")
    elif not set(current["products"]).issubset(reference["products"]):
        raise NumericalCompatibilityError("numerical comparison introduces products outside the reference plan")
    else:
        checked.append("products")
    return {"status": "unverified" if missing else "verified", "checked": checked, "missing": missing}


def _assessment_resolution(protocol, frame, label):
    requested = []
    if frame is not None and "grid_points_requested" in frame:
        requested = sorted({_assessment_count(value, f"{label}: grid_points_requested")
                            for value in frame.grid_points_requested if not pd.isna(value)})
    declared = (protocol or {}).get("grid_points")
    if declared is not None:
        if requested and requested != [declared]:
            raise AssessmentError(f"{label}: grid_points contradicts scenario requests")
        requested = [declared]
    return {"grid_points": requested or None, "steps": (protocol or {}).get("steps")}


def _assessment_numerical_change(before, after):
    finite = _assessment_number(before) and _assessment_number(after)
    zero_changed = bool(finite and before == 0 and after != 0)
    if not finite or zero_changed:
        change = None
    elif before == 0:
        change = 0.0
    else:
        change = (float(after) - float(before)) / float(before) * 100.0
    return {"change_pct": change if change is not None and math.isfinite(change) else None,
            "reference_zero_changed": zero_changed, "nonfinite_values": not finite,
            "nonfinite_change": change is not None and not math.isfinite(change)}


def _assessment_numerical_report(reference_dir, numerical_dir, protocol, scenarios, execution, artifacts):
    report = {"status": "not_supplied", "compatibility": "not_checked", "run_dir": None,
              "max_abs_finite_change_pct": None, "any_failed": None,
              "any_zero_reference_changed": None, "any_nonfinite_values": None,
              "compared_scenarios": 0, "per_metric": {}, "cases": [],
              "meaning": "Resolution comparison of allocation-conditioned scenarios only; not external prediction accuracy, a convergence guarantee, or biological validation.",
              "maxima_basis": "Finite percentage changes only; failures, missing cases and undefined zero-reference changes are reported separately. No acceptance threshold is imposed."}
    if numerical_dir is None:
        return report
    numerical_dir = numerical_dir.resolve()
    if numerical_dir == reference_dir:
        raise NumericalCompatibilityError("numerical_run must be a separate refinement run")
    report["run_dir"] = str(numerical_dir)
    current = _assessment_artifact(numerical_dir, "protocol.json", artifacts, "numerical_run/")
    refined = _assessment_artifact(numerical_dir, "scenarios.csv", artifacts, "numerical_run/")
    auxiliary = {name: _assessment_artifact(numerical_dir, name, artifacts, "numerical_run/")
                 for name in ("ablations.csv", "sensitivities.csv", "external_predictions.csv")}
    refined_execution = _assessment_execution(current, refined, auxiliary, str(numerical_dir))
    sensor = _assessment_artifact(reference_dir, "sensor_model.json", artifacts)
    refined_sensor = _assessment_artifact(numerical_dir, "sensor_model.json", artifacts, "numerical_run/")
    saved = _assessment_artifact(numerical_dir, "numerical_refinement.csv", artifacts, "numerical_run/")
    compatibility = _assessment_numerical_compatibility(
        reference_dir, numerical_dir, protocol, current, sensor, refined_sensor)
    report.update(compatibility=compatibility["status"], compatibility_checks=compatibility,
                  execution=refined_execution, planned_scenarios=refined_execution["planned_scenarios"],
                  resolution={"reference": _assessment_resolution(protocol, scenarios, "reference"),
                              "refined": _assessment_resolution(current, refined, "numerical run")},
                  comparison_artifact_status=artifacts["numerical_run/numerical_refinement.csv"]["status"])
    report["resolution"]["status"] = ("recorded" if all(
        report["resolution"][label][field] is not None for label in ("reference", "refined")
        for field in ("grid_points", "steps")) else "partially_recorded")
    if compatibility["status"] != "verified" or scenarios is None or refined is None:
        report["status"] = "unverified" if compatibility["status"] != "verified" else "missing_scenarios"
        return report
    for field in ("grid_points", "steps"):
        before, after = [report["resolution"][label][field] for label in ("reference", "refined")]
        if before is not None and after is not None:
            if (field == "steps" and after < before) or (field == "grid_points" and min(after) < max(before)):
                raise NumericalCompatibilityError(f"numerical run coarsens {field} instead of refining it")
    old_rows = _assessment_row_map(scenarios, _ASSESS_SCENARIO_KEYS)
    new_rows = _assessment_row_map(refined, _ASSESS_SCENARIO_KEYS)
    planned = {(product, case["name"]) for product in current["products"] for case in current["environments"]}
    report["reference_scenarios_not_refined"] = execution["planned_scenarios"] - len(planned)
    report["any_failed"] = False
    report["any_zero_reference_changed"] = False
    report["any_nonfinite_values"] = False
    report["missing_scenarios"] = []
    input_fields = ("allocation_fraction", "growth_fraction", "couple_stress", "product_targets_used_for_fitting",
                    "carbon_source", "initial_substrate_g_per_l", "carbon_uptake_bound_mmol_c_per_gdcw_h",
                    "oxygen_uptake_bound_mmol_per_gdcw_h", "temperature_c", "initial_biomass_g_per_l", "hours",
                    "stressor", "dose_mM", "prediction_kind", "prediction_basis", "assumptions", "state_source",
                    "state_extrapolated", "latent_upr", "latent_oxidative", "learned_growth_retention", "applied_growth_retention")
    for key in sorted(planned):
        before, after = old_rows.get(key), new_rows.get(key)
        case = {**dict(zip(_ASSESS_SCENARIO_KEYS, key)),
                "status_reference": None if before is None else before["status"],
                "status": None if after is None else after["status"], "metrics": {}}
        report["cases"].append(case)
        for row in (before, after):
            if row is not None and row["status"] not in ("predicted", "pending", "running"):
                report["any_failed"] = True
        if before is None or after is None:
            report["missing_scenarios"].append(dict(zip(_ASSESS_SCENARIO_KEYS, key)))
            continue
        if before["status"] != "predicted" or after["status"] != "predicted":
            continue
        for field in input_fields:
            if field in before and field in after and not _assessment_equal(before[field], after[field]):
                raise NumericalCompatibilityError(f"numerical comparison changes scenario {key} input {field}")
        for field, declared in (("allocation_fraction", "allocation_fraction"), ("growth_fraction", "growth_floor_fraction")):
            if field in after and not _assessment_equal(after[field], current[declared]):
                raise NumericalCompatibilityError(f"numerical scenario {key} contradicts protocol {declared}")
        report["compared_scenarios"] += 1
        for metric in _ASSESS_NUMERICAL_METRICS:
            change = _assessment_numerical_change(before.get(metric), after.get(metric))
            case["metrics"][metric] = change
            report["any_zero_reference_changed"] |= change["reference_zero_changed"]
            report["any_nonfinite_values"] |= change["nonfinite_values"] or change["nonfinite_change"]
    if saved is not None:
        required = [*_ASSESS_SCENARIO_KEYS, "status", "status_reference", *_ASSESS_NUMERICAL_METRICS]
        for metric in _ASSESS_NUMERICAL_METRICS:
            required.extend([f"{metric}_reference", f"{metric}_change_pct", f"{metric}_reference_zero_changed"])
        _assessment_table(saved, required, _ASSESS_SCENARIO_KEYS, "numerical_refinement.csv")
        calculated = _assessment_row_map(pd.DataFrame(report["cases"]), _ASSESS_SCENARIO_KEYS)
        for key, row in _assessment_row_map(saved, _ASSESS_SCENARIO_KEYS).items():
            if key not in calculated or key not in old_rows or key not in new_rows:
                raise AssessmentError(f"numerical_refinement.csv: unknown or unavailable scenario {key}")
            expected = {"status": new_rows[key]["status"], "status_reference": old_rows[key]["status"]}
            for metric in _ASSESS_NUMERICAL_METRICS:
                change = calculated[key]["metrics"].get(metric, {})
                expected.update({metric: new_rows[key].get(metric), f"{metric}_reference": old_rows[key].get(metric),
                                 f"{metric}_change_pct": change.get("change_pct"),
                                 f"{metric}_reference_zero_changed": change.get("reference_zero_changed")})
            if any(not _assessment_equal(row[field], value) for field, value in expected.items()):
                raise AssessmentError(f"numerical_refinement.csv: values contradict linked scenarios for {key}")
    maxima = []
    for metric in _ASSESS_NUMERICAL_METRICS:
        values = [case["metrics"][metric] for case in report["cases"] if metric in case["metrics"]]
        finite = [abs(value["change_pct"]) for value in values if value["change_pct"] is not None]
        maximum = max(finite) if finite else None
        report["per_metric"][metric] = {"finite_change_rows": len(finite), "max_abs_finite_change_pct": maximum,
                                        "zero_reference_changes": sum(value["reference_zero_changed"] for value in values),
                                        "nonfinite_value_rows": sum(value["nonfinite_values"] or value["nonfinite_change"] for value in values)}
        if maximum is not None:
            maxima.append(maximum)
    report["max_abs_finite_change_pct"] = max(maxima) if maxima else None
    if (not execution["run_complete"] or not refined_execution["run_complete"]
            or report["compared_scenarios"] != len(planned) or report["any_failed"]):
        report["status"] = "incomplete"
    elif report["any_nonfinite_values"] or report["any_zero_reference_changed"]:
        report["status"] = "assessed_with_undefined_changes"
    else:
        report["status"] = "assessed"
    return report


def assess_run(run_dir: pathlib.Path, *, numerical_run: pathlib.Path | None = None) -> dict:
    directory = run_dir.resolve()
    artifacts = {}
    protocol = _assessment_artifact(directory, "protocol.json", artifacts)
    scenarios = _assessment_artifact(directory, "scenarios.csv", artifacts)
    auxiliary = {name: _assessment_artifact(directory, name, artifacts)
                 for name in ("ablations.csv", "sensitivities.csv", "external_predictions.csv")}
    execution = _assessment_execution(protocol, scenarios, auxiliary, str(directory))
    external = _assessment_external_report(directory, protocol, artifacts)
    numerics = _assessment_numerical_report(directory, numerical_run, protocol, scenarios, execution, artifacts)
    scientific_result = external["comparison_result"]
    if scientific_result == "beats_baseline_on_scored_rows":
        if (external["coverage_status"] != "complete_on_recorded_targets" or not execution["run_complete"]
                or (protocol or {}).get("product_targets_used_for_fitting") is not False):
            scientific_result = "insufficient_coverage"
        else:
            scientific_result = "beats_baseline_on_scored_rows_only"
    return _json_value({
        "schema_version": 1, "run_dir": str(directory), "validated_model": False,
        "scientific_result": scientific_result, "execution": execution, "external_validation": external,
        "numerical_refinement": numerics, "artifacts": artifacts,
        "evidence_limitations": [
            "Scenarios are allocation-conditioned assumptions, not validated realized production forecasts.",
            "The sensor coupling is an empirical growth-cap, not a mechanistic gene-specific regulatory or ATP-cost fit.",
            "There is no paired observed sensor/product validation in these benchmark artifacts; scenario sensors are simulated.",
            "External scores condition on measured physiology and do not validate the complete latent-stress-to-product chain.",
            "Sensitivities are one-at-a-time assumption checks, not confidence intervals.",
            "Baseline comparisons are descriptive on the same scored rows, not statistical significance or a prospective test.",
            "Missing or censored evidence is not a zero error; numerical resolution agreement is not biological validation.",
            "Only recorded artifacts are assessed; missing integration-step metadata is not inferred from defaults.",
        ],
    })


def _plots(frame: pd.DataFrame, destination: pathlib.Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = frame[frame.status == "predicted"]
    if valid.empty:
        return
    metrics = (("titre_mg_per_l", "Extra accumulated titre (mg/L)"),
               ("productivity_mg_per_l_h", "Average productivity (mg/L/h)"),
               ("yield_g_per_g", "Net product / consumed substrate (g/g)"))
    figure, axes = plt.subplots(3, 1, figsize=(13, 12), constrained_layout=True)
    for axis, (metric, label) in zip(axes, metrics):
        matrix = valid.pivot(index="product", columns="environment", values=metric)
        matrix.plot.bar(ax=axis, logy=True, width=0.8)
        axis.set_ylabel(label)
        axis.set_xlabel("")
        axis.tick_params(axis="x", rotation=0)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(fontsize=7, ncol=2)
    figure.suptitle("Five product-label holdouts: allocation-conditioned scenarios, not validated forecasts")
    figure.savefig(destination, dpi=150)
    plt.close(figure)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run or assess the calibrated sensor-to-GEM product challenge")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output-dir", type=pathlib.Path, help="New output directory for model execution")
    mode.add_argument("--assess-run", type=pathlib.Path, help="Read existing artifacts without training or prediction")
    parser.add_argument("--numerical-run", type=pathlib.Path, help="Linked numerical refinement to include in assessment")
    parser.add_argument("--assessment-output", type=pathlib.Path,
                        help="New JSON report path; defaults to ASSESS_RUN/assessment.json; never overwrites")
    parser.add_argument("--sensor-table", type=pathlib.Path,
                        default=ROOT / "outputs" / "sensor_characterisation.csv")
    parser.add_argument("--heldout-plate")
    parser.add_argument("--products", nargs="+", choices=product_names())
    parser.add_argument("--environments", nargs="+", choices=[case.name for case in environments(24.0)])
    parser.add_argument("--compare-to", type=pathlib.Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--grid-points", type=int, default=5)
    parser.add_argument("--steps", type=int, default=120)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--allocation", type=float, default=0.25)
    parser.add_argument("--growth-fraction", type=float, default=0.9)
    parser.add_argument("--skip-sensitivity", action="store_true")
    parser.add_argument("--skip-external", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    if args.assess_run is not None:
        if args.compare_to is not None:
            parser.error("--compare-to is for model execution; use --numerical-run for assessment")
        destination = args.assessment_output or args.assess_run / "assessment.json"
        if destination.exists() or destination.is_symlink():
            raise SystemExit(f"Refusing to overwrite an existing assessment output: {destination}")
        try:
            report = assess_run(args.assess_run, numerical_run=args.numerical_run)
            encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
            with destination.open("x", encoding="utf-8") as stream:
                stream.write(encoded)
        except FileExistsError as exc:
            raise SystemExit(f"Refusing to overwrite an existing assessment output: {destination}") from exc
        except (AssessmentError, OSError) as exc:
            raise SystemExit(f"{type(exc).__name__}: {exc}") from exc
        print(f"Assessment report: {destination.resolve()}\nScientific result: {report['scientific_result']}\n"
              "Report generation completed; validated_model=false. Exit status is not a science pass.", flush=True)
        return 0
    if args.numerical_run is not None or args.assessment_output is not None:
        parser.error("--numerical-run and --assessment-output require --assess-run")
    if not np.isfinite(args.allocation) or not 0 <= args.allocation <= 1:
        parser.error("--allocation must be in [0, 1]")
    if not np.isfinite(args.growth_fraction) or not 0 < args.growth_fraction <= 1:
        parser.error("--growth-fraction must be in (0, 1]")
    if args.grid_points < 3 or args.steps < 2 or args.seed < 0 or args.workers < 1:
        parser.error("grid-points >= 3, steps >= 2, seed >= 0 and workers >= 1 are required")
    cases = tuple(case for case in environments(args.hours)
                  if args.environments is None or case.name in args.environments)
    out = args.output_dir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"Refusing to overwrite an existing run: {out}; choose a new output directory")
    out.mkdir(parents=True, exist_ok=True)
    requested_products = tuple(dict.fromkeys(args.products or product_names()))
    model_path = paths.require(paths.ec_yeast_gem(), "the frozen EC yeast model", "YSTWIN_EC_YEAST_GEM")
    sensor_data = pd.read_csv(args.sensor_table)
    plates = tuple(sorted(sensor_data.plate.unique()))
    if len(plates) < 2:
        raise SystemExit("A held-out plate and at least one different training plate are required")
    heldout = args.heldout_plate or plates[-1]
    if heldout not in plates:
        raise SystemExit(f"Held-out plate {heldout!r} is absent from the sensor table")
    training = tuple(plate for plate in plates if plate != heldout)
    protocol = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "run_complete": False,
        "evaluation": "retrospective product-measurement-withheld challenge; not a blind prospective test",
        "products": requested_products,
        "environments": [asdict(case) for case in cases],
        "training_plates": training,
        "training_seed": args.seed,
        "heldout_plate": heldout,
        "product_targets_used_for_fitting": False,
        "allocation_fraction": args.allocation,
        "growth_floor_fraction": args.growth_fraction,
        "independent_scenario_workers": args.workers,
        "grid_points": args.grid_points,
        "steps": args.steps,
        "allocation_status": "frozen modeling assumption, not fitted to any product target",
        "sensitivity_status": "one-at-a-time assumption checks, not confidence intervals",
        "scope": [
            "Native products with shared EC enzyme budgets and integrated wild-type protein abundances",
            "GEM biomass composition is preexisting product-related prior information, not unseen chemistry",
            "Learned growth-cap coupling is empirical; gene-specific regulatory/ATP costs are not identified",
            "Scenario sensors are generated by the calibrated reporter ODE, not observed on product cultures",
            "External scores condition on measured physiology and do not validate the complete latent-stress chain",
            "No target-derived product calibration, measured product pool, or product flux enters fitting",
        ],
        "input_hashes": {"sensor_table": _hash(args.sensor_table), "ec_model": _hash(model_path),
                         "proteome": _hash(ROOT / "data/proteome/paxdb_scerevisiae_integrated.tsv")},
        "versions": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "pandas", "cobra")},
        "code_hashes": {name: _hash(ROOT / name) for name in (
            "src/ystwin/hybrid.py", "src/ystwin/analysis/hybrid_stress.py",
            "src/ystwin/analysis/hybrid_validation.py", "src/ystwin/analysis/latent.py",
            "src/ystwin/analysis/stress_model.py", "src/ystwin/reporter.py",
            "src/ystwin/growth.py", "src/ystwin/readings.py",
            "src/ystwin/fba/product_panel.py", "src/ystwin/fba/dynamic.py",
            "src/ystwin/fba/dynamic_rates.py", "src/ystwin/fba/physiology.py",
            "src/ystwin/fba/solver.py", "src/ystwin/generator/context.py",
            "src/ystwin/pathway/proteome.py", "scripts/run_hybrid_benchmark.py")},
        "git_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                       capture_output=True, text=True, check=True).stdout.strip(),
    }
    _write_json(out / "protocol.json", protocol)
    print(f"Training only on {training}; held-out sensor plate {heldout}", flush=True)
    sensor = fit_sensor_model(sensor_data, training_plates=training, seed=args.seed)
    _write_json(out / "sensor_model.json", sensor.metadata)
    sensor.validation_table(sensor_data[sensor_data.plate == heldout]).to_csv(
        out / "heldout_sensor_scores.csv", index=False)
    base, solver = load_model(model_path)
    prepared, enzyme_inputs = prepare_panel_model(base)
    enzyme_inputs.to_csv(out / "enzyme_inputs.csv", index=False)
    protocol["solver"] = str(solver)
    jobs = [{"kind": "scenario", "product": product, "environment": environment}
            for product in requested_products for environment in cases]
    stress_case = next((case for case in cases if case.stressor == "DTT"), None)
    if stress_case is not None:
        jobs.extend({"kind": "ablation", "product": product, "environment": stress_case,
                     "options": {"couple_stress": False}} for product in requested_products)
    if not args.skip_sensitivity:
        for product in requested_products:
            for field, values in (("allocation_fraction", (0.1, 0.5)),
                                  ("growth_fraction", (0.8, 0.95))):
                jobs.extend({"kind": "sensitivity", "product": product,
                             "environment": cases[0], "options": {field: value},
                             "sensitivity_parameter": field, "sensitivity_value": value}
                            for value in values)
    if not args.skip_external:
        target_path = ROOT / "data/hybrid_benchmark/external_targets.csv"
        targets = pd.read_csv(target_path)
        covariates = targets.drop(columns=["observed", "uncertainty", "uncertainty_kind",
                                           "observation_kind", "detection_limit"])
        jobs.extend({"kind": "external", "target": target}
                    for target in covariates.to_dict("records"))
        protocol["external_targets_sha256"] = _hash(target_path)
    protocol["planned_jobs"] = len(jobs)
    _write_json(out / "protocol.json", protocol)
    print(f"Executing {len(jobs)} jobs with {args.workers} isolated scenario workers", flush=True)
    rows, provenance, ablations, sensitivities, external = [], [], [], [], []
    for job, (row, result) in zip(jobs, _job_results(jobs, sensor, prepared, model_path, args)):
        if job["kind"] == "external":
            external.append(row)
            continue
        if job["kind"] == "ablation":
            row["ablation"] = "no_learned_stress_growth_cap"
            ablations.append(row)
            continue
        if job["kind"] == "sensitivity":
            row.update(sensitivity_parameter=job["sensitivity_parameter"],
                       sensitivity_value=job["sensitivity_value"])
            sensitivities.append(row)
            continue
        rows.append(row)
        if result is not None:
            provenance.append(result.provenance())
            print(f"{row['product']:12} {row['environment']:24} {row['titre_mg_per_l']:10.4g} mg/L  "
                  f"{row['yield_g_per_g']:.4g} g/g", flush=True)
        else:
            print(f"{row['product']:12} {row['environment']:24} {row['status']}: {row['reason']}", flush=True)
        pd.DataFrame(rows).to_csv(out / "scenarios.csv", index=False)
    frame = pd.DataFrame(rows)
    pd.DataFrame(ablations).to_csv(out / "ablations.csv", index=False)
    if sensitivities:
        pd.DataFrame(sensitivities).to_csv(out / "sensitivities.csv", index=False)
    _write_json(out / "case_provenance.json", provenance)
    if not args.skip_external:
        predictions = pd.DataFrame(external)
        predictions.to_csv(out / "external_predictions.csv", index=False)
        scores = score_external_targets(targets, predictions)
        scores.to_csv(out / "external_scores.csv", index=False)
        summary = summarize_scores(scores)
        summary.to_csv(out / "external_summary.csv", index=False)
        baselines = _external_baselines(prepared, targets)
        baseline_scores = score_external_targets(targets, baselines)
        baseline_scores.to_csv(out / "external_baseline_scores.csv", index=False)
        summarize_scores(baseline_scores).to_csv(out / "external_baseline_summary.csv", index=False)
        matched = scores[scores.score_status == "scored"]
        matched_baseline = baseline_scores.loc[matched.index]
        comparison = matched[["product", "metric", "abs_error"]].copy()
        comparison["baseline_abs_error"] = matched_baseline.abs_error
        comparison.groupby(["product", "metric"], as_index=False).agg(
            n_compared=("abs_error", "count"), model_mae=("abs_error", "mean"),
            baseline_mae=("baseline_abs_error", "mean")
        ).to_csv(out / "external_baseline_comparison.csv", index=False)
        protocol["external_targets_sha256"] = _hash(target_path)
        print("\nExternal challenge (physiology-conditioned, not whole-chain validation):", flush=True)
        print(summary.to_string(index=False), flush=True)
        pd.DataFrame([
            {"product": name, "target_rows": int((targets["product"] == name).sum()),
             "target_conditions": int(targets.loc[targets["product"] == name, "case_id"].nunique())}
            for name in requested_products
        ]).to_csv(out / "validation_coverage.csv", index=False)
    if args.compare_to is not None:
        _compare_numerics(frame, protocol, args.compare_to.resolve(), out / "numerical_refinement.csv")
        protocol["numerical_reference"] = str(args.compare_to.resolve())
    if not args.no_plots:
        _plots(frame, out / "product_scenarios.png")
    protocol.update(run_complete=True, scenario_rows=len(frame),
                    successful_scenarios=int((frame.status == "predicted").sum()))
    _write_json(out / "protocol.json", protocol)
    print(f"\nResults: {out}", flush=True)
    return 0 if protocol["successful_scenarios"] == len(frame) else 1


if __name__ == "__main__":
    raise SystemExit(main())
