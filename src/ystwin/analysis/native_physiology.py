"""Empirical growth-dependent native exchange laws, not molecular mechanism discovery."""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import asdict, dataclass, field
import hashlib
from io import BytesIO
import json
import math
from numbers import Real
from pathlib import Path
from typing import Iterable

import pandas as pd

__all__ = [
    "TRAIN_GROWTH_RATES_PER_H", "HOLDOUT_GROWTH_RATES_PER_H",
    "NativeChemostatData", "NativeExchangeModel", "load_native_chemostat_data",
    "fit_native_exchange_model", "evaluate_native_exchange_model",
]

TRAIN_GROWTH_RATES_PER_H = (.025, .05, .10, .20, .25, .30, .40)
HOLDOUT_GROWTH_RATES_PER_H = (.15, .28, .35)
_SOURCE_PATH = "data/physiology/chemostatData_VanHoek1998.tsv"
_SOURCE_SHA256 = "4c01b3465ea1410c86709568da44a32b1de652c21a88de623bbf566838506117"
_DOI = "10.1128/AEM.64.11.4226-4233.1998"
_RATES = tuple(sorted((*TRAIN_GROWTH_RATES_PER_H, *HOLDOUT_GROWTH_RATES_PER_H)))
_LINES = {rate: line for line, rate in enumerate(_RATES, 2)}
_UNITS = {
    "GlucoseUptake": "mmol/gDW/h", "O2uptake": "mmol/gDW/h",
    "CO2production": "mmol/gDW/h", "Ethanol": "mmol/gDW/h", "Glycerol": "mmol/gDW/h",
    "Acetate": "mmol/gDW/h", "Pyruvate": "mmol/gDW/h",
    "BiomassYield_g_per_g": "gDW/g_glucose", "CarbonRecovery_pct": "%",
}
_PREDICTION_STATUSES = (
    "quantified", "unquantified_censored_support", "unsupported_growth_rate",
)
_COLUMNS = (
    "row_id", "source_line", "source_column", "source_doi", "source_table_sha256",
    "growth_rate_per_h", "split", "observable_id", "unit", "reported_value", "value",
    "observation_status", "detection_limit",
)


def _source() -> dict:
    return {
        "doi": _DOI, "pmid": "9797269", "pmcid": "PMC106631", "year": 1998,
        "authors": ["Pim van Hoek", "Johannes P. van Dijken", "Jack T. Pronk"],
        "title": "Effect of Specific Growth Rate on Fermentative Capacity of Baker's Yeast",
        "journal": "Applied and Environmental Microbiology", "volume": 64, "issue": 11,
        "pages": "4226-4233", "table": "Table 1", "table_page": 4227,
        "verified_primary_url": (
            "https://repository.tudelft.nl/file/File_0981d011-0515-40cd-bc5f-b5004cd9ca42"
        ),
        "verification": "Native TSV rows checked against the published Table 1, page 4227",
        "table_relative_path": _SOURCE_PATH, "table_sha256": _SOURCE_SHA256,
        "conditions": {
            "organism": "Saccharomyces cerevisiae", "strain": "DS28911",
            "cultivation": "aerobic glucose-limited chemostat", "medium": "defined mineral medium",
            "temperature_c": 30.0, "ph": 5.0,
        },
        "growth_basis": "Specific growth rate equals dilution rate at chemostat steady state",
        "exchange_convention": "Positive uptake or production magnitudes, not signed GEM fluxes",
        "biomass_carbon": {
            "mass_fraction": .48, "status": "source_assumption_for_carbon_recovery",
        },
        "censoring": {
            "footnote": "Table 1 footnote c: 0, below detection limit",
            "detection_limit": None, "numeric_detection_limit_reported": False,
        },
        "uncertainty": "Table 1 does not report row-level SD or SEM",
    }


def _split_metadata() -> dict:
    return {
        "train": list(TRAIN_GROWTH_RATES_PER_H),
        "holdout": list(HOLDOUT_GROWTH_RATES_PER_H), "unit": "1/h",
        "scope": "predeclared within-study condition holdout, not independent-laboratory validation",
    }


@dataclass(frozen=True)
class NativeChemostatData:
    observations: pd.DataFrame = field(repr=False)
    metadata: dict = field(repr=False)

    def split(self) -> dict[str, pd.DataFrame]:
        """Split whole conditions by dilution rate, never by row names or outcomes."""
        return {
            name: self.observations.loc[self.observations.growth_rate_per_h.isin(rates)].copy()
            for name, rates in (
                ("train", TRAIN_GROWTH_RATES_PER_H), ("holdout", HOLDOUT_GROWTH_RATES_PER_H),
            )
        }


def _context(growth: float, observable: str, line: int) -> dict:
    return {
        "row_id": f"vanHoek1998:D={growth:g}", "source_line": line,
        "source_column": observable, "source_doi": _DOI, "source_table_sha256": _SOURCE_SHA256,
        "growth_rate_per_h": growth,
        "split": "train" if growth in TRAIN_GROWTH_RATES_PER_H else "holdout",
        "observable_id": observable, "unit": _UNITS[observable], "detection_limit": None,
    }


def load_native_chemostat_data(path: str | Path | None = None) -> NativeChemostatData:
    """Load the verified TSV as 90 tidy observations; an explicit path must have identical bytes."""
    path = Path(path) if path is not None else Path(__file__).resolve().parents[3] / _SOURCE_PATH
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != _SOURCE_SHA256:
        raise ValueError("Native Table 1 SHA-256 mismatch; unverified bytes cannot inherit its provenance")
    table = pd.read_csv(BytesIO(payload), sep="\t")
    if list(table.columns) != ["Drate", *_UNITS] or tuple(table.Drate) != _RATES:
        raise ValueError("Unexpected native Table 1 fields or dilution rates")
    records = []
    for line, row in enumerate(table.to_dict("records"), 2):
        for observable in _UNITS:
            reported = float(row[observable])
            records.append({
                **_context(float(row["Drate"]), observable, line),
                "reported_value": reported, "value": reported if reported != 0 else None,
                "observation_status": "quantified" if reported != 0 else "below_detection_limit",
            })
    observations = pd.DataFrame(records, columns=_COLUMNS)
    observations["value"] = pd.Series([row["value"] for row in records], dtype=object)
    return NativeChemostatData(observations, {
        "source": _source(), "units": dict(_UNITS), "predeclared_split": _split_metadata(),
    })


def _number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def _validated_rows(rows: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(rows, pd.DataFrame) or rows.empty:
        raise ValueError("Native observations must be a nonempty DataFrame")
    if not rows.columns.is_unique or not set(_COLUMNS).issubset(rows.columns):
        raise ValueError("Native observations require explicit per-value provenance, units and censoring")
    rows = rows.loc[:, list(_COLUMNS)].copy()
    for row in rows.to_dict("records"):
        growth = _number(row["growth_rate_per_h"], "growth rate")
        observable = row["observable_id"]
        if growth not in _LINES or observable not in _UNITS:
            raise ValueError("Observation is not a native Table 1 condition/readout")
        context = _context(growth, observable, _LINES[growth])
        for key in ("source_line", "source_column", "source_doi", "source_table_sha256", "unit", "split"):
            if row[key] != context[key]:
                raise ValueError(f"Invalid native observation {key}")
        if not pd.isna(row["detection_limit"]):
            raise ValueError("No numeric detection limit is reported by this source")
        reported = _number(row["reported_value"], "reported value")
        if row["observation_status"] == "below_detection_limit":
            if reported != 0 or not pd.isna(row["value"]):
                raise ValueError("Below-detection observations require reported zero and an unquantified value")
        elif row["observation_status"] == "quantified":
            value = _number(row["value"], "quantified value")
            if value <= 0 or value != reported:
                raise ValueError("A quantified observation must equal its positive reported value")
        else:
            raise ValueError("Unknown native observation status")
    if rows.duplicated(["growth_rate_per_h", "observable_id"]).any():
        raise ValueError("Duplicate native condition/readout")
    if not rows.groupby("growth_rate_per_h").size().eq(len(_UNITS)).all():
        raise ValueError("Retain all nine native readouts for every supplied condition")
    return rows.sort_values(["growth_rate_per_h", "observable_id"])


@dataclass(frozen=True)
class _Knot:
    growth_rate_per_h: float
    value: float | None
    reported_value: float
    observation_status: str
    source_line: int


@dataclass(frozen=True)
class NativeExchangeModel:
    _knots: tuple[tuple[_Knot, ...], ...] = field(repr=False)

    def to_dict(self) -> dict:
        """Detached, JSON-safe freeze payload containing every actual fitted knot."""
        knots = {name: [asdict(knot) for knot in series] for name, series in zip(_UNITS, self._knots)}
        digest = hashlib.sha256(json.dumps(
            knots, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode()).hexdigest()
        return {
            "schema_version": 1, "model_kind": "empirical_growth_dependent_exchange_law",
            "method": "piecewise_linear", "extrapolation": "unsupported",
            "interpretation": "Condition-specific empirical interpolation, not molecular mechanism discovery",
            "input": {"name": "growth_rate_per_h", "unit": "1/h", "role": "imposed_condition_not_scored"},
            "source": _source(), "units": dict(_UNITS), "predeclared_split": _split_metadata(),
            "training_growth_rates_per_h": [k.growth_rate_per_h for k in self._knots[0]],
            "training_rows_sha256": digest, "knots": knots,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> NativeExchangeModel:
        """Validate/reload fitted values without reading the TSV or any omitted observations."""
        try:
            if not isinstance(payload, dict) or payload.get("source") != _source():
                raise ValueError("source provenance mismatch")
            if set(payload["knots"]) != set(_UNITS):
                raise ValueError("incomplete native readouts")
            records = []
            for observable in _UNITS:
                for knot in payload["knots"][observable]:
                    records.append({
                        **_context(knot["growth_rate_per_h"], observable, knot["source_line"]),
                        **knot, "unit": payload["units"][observable],
                    })
            model = fit_native_exchange_model(pd.DataFrame(records))
            if model.to_dict() != payload:
                raise ValueError("noncanonical or inconsistent model payload")
            return model
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid native model payload: {exc}") from exc

    def predict(self, growth_rates: Iterable[float]) -> list[dict]:
        """Return every requested readout; unquantified/unsupported values are explicit None."""
        grid = [k.growth_rate_per_h for k in self._knots[0]]
        predictions = []
        for condition_index, rate in enumerate(growth_rates):
            growth = _number(rate, "growth rate")
            support_indices = []
            if grid[0] <= growth <= grid[-1]:
                upper = bisect_left(grid, growth)
                support_indices = [upper] if grid[upper] == growth else [upper - 1, upper]
            for observable, series in zip(_UNITS, self._knots):
                support = [series[i] for i in support_indices]
                value, method = None, None
                if not support:
                    status = "unsupported_growth_rate"
                elif any(k.observation_status == "below_detection_limit" for k in support):
                    status = "unquantified_censored_support"
                else:
                    status = "quantified"
                    if len(support) == 1:
                        value, method = support[0].value, "observed_knot"
                    else:
                        left, right = support
                        fraction = (growth - left.growth_rate_per_h) / (right.growth_rate_per_h - left.growth_rate_per_h)
                        value = (1 - fraction) * left.value + fraction * right.value
                        method = "linear_interpolation"
                predictions.append({
                    "condition_index": condition_index, "growth_rate_per_h": growth,
                    "observable_id": observable, "unit": _UNITS[observable],
                    "prediction": value, "prediction_status": status, "method": method,
                    "support": [asdict(k) for k in support], "bounds": None, "detection_limit": None,
                })
        return predictions


def fit_native_exchange_model(training_rows: pd.DataFrame) -> NativeExchangeModel:
    """Fit only supplied training conditions; censored values occupy support but are not targets."""
    rows = _validated_rows(training_rows)
    if rows.growth_rate_per_h.isin(HOLDOUT_GROWTH_RATES_PER_H).any():
        raise ValueError("training_rows contain predeclared native holdout conditions")
    if rows.growth_rate_per_h.nunique() < 2:
        raise ValueError("At least two training growth rates are required for a curve")
    series = []
    for observable in _UNITS:
        series.append(tuple(
            _Knot(
                float(row.growth_rate_per_h),
                float(row.value) if row.observation_status == "quantified" else None,
                float(row.reported_value), row.observation_status, int(row.source_line),
            )
            for row in rows.loc[rows.observable_id.eq(observable)].itertuples(index=False)
        ))
    return NativeExchangeModel(tuple(series))


def evaluate_native_exchange_model(model: NativeExchangeModel, validation_rows: pd.DataFrame) -> dict:
    """Score quantified pairs only, preserving every readout and both coverage denominators."""
    rows = _validated_rows(validation_rows)
    growth_rates = sorted(rows.growth_rate_per_h.unique())
    predictions = {
        (p["growth_rate_per_h"], p["observable_id"]): p for p in model.predict(growth_rates)
    }
    readouts = []
    for row in rows.to_dict("records"):
        prediction = predictions[(row["growth_rate_per_h"], row["observable_id"])]
        observed = float(row["value"]) if row["observation_status"] == "quantified" else None
        scored = observed is not None and prediction["prediction_status"] == "quantified"
        readouts.append({
            **prediction, "row_id": str(row["row_id"]), "source_line": int(row["source_line"]),
            "reported_value": float(row["reported_value"]), "observed_value": observed,
            "observation_status": row["observation_status"], "scored": scored,
            "residual": prediction["prediction"] - observed if scored else None,
        })
    metrics = {}
    for observable, unit in _UNITS.items():
        selected = [row for row in readouts if row["observable_id"] == observable]
        errors = [row["residual"] for row in selected if row["scored"]]
        n, scored = len(selected), len(errors)
        quantified = sum(row["observation_status"] == "quantified" for row in selected)
        metrics[observable] = {
            "unit": unit, "n_total": n, "n_observed_quantified": quantified,
            "n_observed_censored": n - quantified, "n_scored": scored,
            "coverage": scored / n,
            "quantified_observation_coverage": scored / quantified if quantified else None,
            "prediction_status_counts": {
                status: sum(row["prediction_status"] == status for row in selected)
                for status in _PREDICTION_STATUSES
            },
            "mae": math.fsum(abs(error) / scored for error in errors) if scored else None,
            "rmse": math.hypot(*errors) / math.sqrt(scored) if scored else None,
        }
    n_scored = sum(row["scored"] for row in readouts)
    return {
        "scope": "native condition evaluation; growth is imposed, not scored",
        "n_conditions": len(growth_rates), "n_readouts": len(readouts), "n_scored": n_scored,
        "coverage": n_scored / len(readouts), "readouts": readouts, "metrics": metrics,
    }
