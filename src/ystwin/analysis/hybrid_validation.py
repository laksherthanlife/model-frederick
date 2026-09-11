from __future__ import annotations

import math
from numbers import Real

import numpy as np
import pandas as pd


_KEYS = ["study_id", "case_id", "product", "metric"]
_TARGET_FIELDS = ["observed", "units", "observation_kind", "detection_limit"]
_PREDICTION_FIELDS = [
    "predicted",
    "prediction_units",
    "prediction_status",
    "comparable",
    "prediction_basis",
    "reason",
]
_ERROR_COLUMNS = ["abs_error", "relative_error_pct", "absolute_percent_error", "fold_error"]
_SCORE_COLUMNS = ["score_status", "score_reason", *_ERROR_COLUMNS]
_STATUSES = ["scored", "censored", "unsupported", "infeasible", "missing", "invalid"]
_COUNT_COLUMNS = [
    "n_total",
    "n_scored",
    "n_censored",
    "n_unsupported",
    "n_infeasible",
    "n_missing",
    "n_invalid",
    "n_percent_scored",
    "n_fold_scored",
]
_SUMMARY_METRICS = ["mae", "median_abs_percent_error", "median_fold_error"]
_SUMMARY_COLUMNS = ["product", "metric", "units", *_COUNT_COLUMNS, *_SUMMARY_METRICS]


def _has_text(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _validate_frame(frame: pd.DataFrame, required: list[str], name: str) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    if not frame.columns.is_unique:
        raise ValueError(f"{name} contains duplicate column names")
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")
    for key in _KEYS:
        if not all(_has_text(value) for value in frame[key]):
            raise ValueError(f"{name} key {key} must contain nonempty strings")
    if frame.duplicated(_KEYS).any():
        raise ValueError(f"{name} contains duplicate target keys: {_KEYS}")


def _valid_number(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        return False
    try:
        return math.isfinite(value) and value >= 0
    except (OverflowError, TypeError, ValueError):
        return False


def _unscored(status: str, reason: str) -> dict[str, object]:
    return {
        "score_status": status,
        "score_reason": reason,
        **dict.fromkeys(_ERROR_COLUMNS, np.nan),
    }


def _score_row(row: dict[str, object], present: bool) -> dict[str, object]:
    if not present:
        return _unscored("missing", "missing_prediction")
    status = row["prediction_status"]
    if not isinstance(status, str) or status not in {"predicted", "unsupported", "infeasible"}:
        return _unscored("invalid", "invalid_prediction_status")
    if status != "predicted":
        return _unscored(status, f"prediction_{status}")
    comparable = row["comparable"]
    if not isinstance(comparable, (bool, np.bool_)):
        return _unscored("invalid", "comparable_must_be_boolean")
    if not comparable:
        return _unscored("invalid", "not_comparable")
    if (
        not _has_text(row["units"])
        or not isinstance(row["prediction_units"], str)
        or row["prediction_units"] != row["units"]
    ):
        return _unscored("invalid", "units_mismatch")
    if not _has_text(row["prediction_basis"]):
        return _unscored("invalid", "prediction_basis_required")
    if not _valid_number(row["predicted"]):
        return _unscored("invalid", "prediction_must_be_finite_and_nonnegative")
    kind = row["observation_kind"]
    if not isinstance(kind, str) or kind not in {"quantified", "below_detection_limit"}:
        return _unscored("invalid", "invalid_observation_kind")
    if kind == "below_detection_limit":
        limit = row["detection_limit"]
        missing_limit = pd.isna(limit)
        if (
            isinstance(missing_limit, (bool, np.bool_)) and missing_limit
        ) or (isinstance(limit, str) and not limit.strip()):
            return _unscored("censored", "censored_unknown_limit")
        if not _valid_number(limit) or limit <= 0:
            return _unscored("invalid", "invalid_detection_limit")
        return _unscored("censored", "censored_known_limit")
    if not _valid_number(row["observed"]):
        return _unscored("invalid", "invalid_observed")
    observed = float(row["observed"])
    predicted = float(row["predicted"])
    result = _unscored("scored", "quantified_zero" if observed == 0 else "quantified_positive")
    result["abs_error"] = abs(predicted - observed)
    if observed > 0:
        relative_error = (predicted - observed) / observed * 100.0
        result["relative_error_pct"] = relative_error
        result["absolute_percent_error"] = abs(relative_error)
        if predicted > 0:
            result["fold_error"] = max(predicted / observed, observed / predicted)
    return result


def score_external_targets(targets: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    _validate_frame(targets, [*_KEYS, *_TARGET_FIELDS], "targets")
    _validate_frame(predictions, [*_KEYS, *_PREDICTION_FIELDS], "predictions")
    reserved = {*_SCORE_COLUMNS, "is_censored"}
    for name, frame, columns in (
        ("targets", targets, reserved | set(_PREDICTION_FIELDS)),
        ("predictions", predictions, reserved),
    ):
        collisions = sorted(set(frame.columns) & columns)
        if collisions:
            raise ValueError(f"{name} contains reserved scoring columns: {collisions}")
    target_keys = pd.MultiIndex.from_frame(targets[_KEYS])
    prediction_keys = pd.MultiIndex.from_frame(predictions[_KEYS])
    unknown = ~prediction_keys.isin(target_keys)
    if unknown.any():
        raise ValueError(f"predictions contain unknown target keys: {prediction_keys[unknown].tolist()}")
    present = target_keys.isin(prediction_keys)
    scores = targets.merge(
        predictions,
        on=_KEYS,
        how="left",
        sort=False,
        validate="one_to_one",
        suffixes=("", "_prediction"),
    ).copy(deep=True)
    if not scores.columns.is_unique:
        raise ValueError("prediction metadata creates duplicate output columns")
    results = pd.DataFrame(
        [_score_row(row, exists) for row, exists in zip(scores.to_dict("records"), present)],
        columns=_SCORE_COLUMNS,
    )
    results[_ERROR_COLUMNS] = results[_ERROR_COLUMNS].astype(float)
    scores["is_censored"] = scores["observation_kind"].eq("below_detection_limit").fillna(False).astype(bool)
    scores = pd.concat([scores, results], axis=1)
    scores.index = targets.index.copy()
    return scores


def _median_error(errors: pd.Series) -> float:
    values = sorted(float(value) for value in errors)
    if not values:
        return np.nan
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    lower, upper = values[middle - 1:middle + 1]
    if lower == upper:
        return lower
    return lower + (upper - lower) / 2.0


def summarize_scores(scores: pd.DataFrame) -> pd.DataFrame:
    _validate_frame(
        scores,
        [*_KEYS, "units", "observation_kind", "score_status", *_ERROR_COLUMNS],
        "scores",
    )
    if not scores["score_status"].isin(_STATUSES).all():
        raise ValueError("scores contains unknown score_status values")
    summaries = []
    for (product, metric), group in scores.groupby(
        ["product", "metric"], sort=True, observed=True, dropna=False
    ):
        if group["units"].nunique(dropna=False) != 1:
            raise ValueError(f"cannot summarize mixed units for product={product}, metric={metric}")
        scored = group.loc[group["score_status"].eq("scored")]
        absolute = scored["abs_error"].dropna()
        percent = scored["absolute_percent_error"].dropna()
        fold = scored["fold_error"].dropna()
        counts = group["score_status"].value_counts()
        summaries.append(
            {
                "product": product,
                "metric": metric,
                "units": group["units"].iloc[0],
                "n_total": len(group),
                "n_scored": len(scored),
                "n_censored": int(group["observation_kind"].eq("below_detection_limit").sum()),
                **{f"n_{status}": int(counts.get(status, 0)) for status in ("unsupported", "infeasible", "missing", "invalid")},
                "n_percent_scored": len(percent),
                "n_fold_scored": len(fold),
                "mae": math.fsum(value / len(absolute) for value in absolute) if len(absolute) else np.nan,
                "median_abs_percent_error": _median_error(percent),
                "median_fold_error": _median_error(fold),
            }
        )
    return pd.DataFrame(summaries, columns=_SUMMARY_COLUMNS).astype(
        {**dict.fromkeys(_COUNT_COLUMNS, "int64"), **dict.fromkeys(_SUMMARY_METRICS, "float64")}
    )
