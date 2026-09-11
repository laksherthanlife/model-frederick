"""G1: is an optical channel quantitative for this well?

Every latent and flux quantity downstream is a function of growth rate, which is a function of the
OD trace. Thresholds are explicit arguments carried on the result, because they are calibration
claims about one reader and plate, and a result is auditable only if they travel with it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..readings import RawOD, RawRFU, require

__all__ = ["OpticalQualityGate", "WellQualityResult", "assess_plate", "assess_well"]


@dataclass(frozen=True)
class OpticalQualityGate:
    """Thresholds defining "quantitative" for one reader and plate format.

    Args:
        od_linear_max: Highest raw OD still on the reader's linear response. In a
            96-well plate the optical path is roughly 0.3-0.5 cm, so absorbance
            departs from linearity well below the cuvette figure of 1.0. This
            default is a placeholder until a dilution series is run: it is the
            single most important number to measure for the platform.
        min_od_above_blank: Smallest blank-corrected OD counted as real biomass.
        max_decline_fraction: Largest share of intervals whose OD may fall before
            the trace is treated as settling, flocculation or lysis rather than
            growth.
        min_od_fold: Smallest robust end/start OD ratio that makes a growth rate
            identifiable.
        min_reporter_above_background: Smallest background-corrected reporter
            signal counted as real.
    """

    od_linear_max: float = 1.0
    min_od_above_blank: float = 0.02
    max_decline_fraction: float = 0.15
    min_od_fold: float = 1.15
    min_reporter_above_background: float = 50.0


@dataclass(frozen=True)
class WellQualityResult:
    """Per-well verdict with the metrics and thresholds that produced it."""

    well: str
    failures: tuple[str, ...]
    metrics: dict[str, float]
    gate: OpticalQualityGate = field(repr=False)

    @property
    def passed(self) -> bool:
        return not self.failures

    def summary(self) -> str:
        state = "PASS" if self.passed else "FAIL: " + ", ".join(self.failures)
        return f"{self.well}: {state}"


def _robust_fold_change(series: np.ndarray, frac: float = 0.1) -> float:
    n = max(3, int(len(series) * frac))
    start = float(np.median(series[:n]))
    end = float(np.median(series[-n:]))
    if start <= 0 or end <= 0:
        return float("nan")
    return end / start


def assess_well(
    times_h: np.ndarray,
    optical_density: RawOD,
    od_blank: float,
    rfu: RawRFU | None = None,
    rfu_background: float = 0.0,
    gate: OpticalQualityGate = OpticalQualityGate(),
    well: str = "well",
) -> WellQualityResult:
    """Check one well's optical channels against the gate.

    Args:
        times_h: Ascending time grid, hours.
        optical_density: Raw OD for the well, blank included.
        od_blank: Media-only OD.
        rfu: Optional reporter trace, raw.
        rfu_background: Media plus autofluorescence baseline.
        gate: Thresholds to apply.
        well: Label carried through to the result.
    """
    density = require(optical_density, RawOD, name="optical_density")
    raw = density.array
    corrected = density.minus_blank(od_blank).array

    decline_fraction = float(np.mean(np.diff(raw) < 0)) if raw.size > 1 else 0.0
    metrics = {
        "max_raw_od": float(np.nanmax(raw)),
        "min_od_above_blank": float(np.nanmin(corrected)),
        "decline_fraction": decline_fraction,
        "od_fold_change": _robust_fold_change(np.clip(corrected, 1e-9, None)),
    }

    failures: list[str] = []
    if metrics["max_raw_od"] > gate.od_linear_max:
        failures.append("linear_range")
    if metrics["min_od_above_blank"] < gate.min_od_above_blank:
        failures.append("above_blank")
    if metrics["decline_fraction"] > gate.max_decline_fraction:
        failures.append("sustained_decline")
    fold = metrics["od_fold_change"]
    if not np.isfinite(fold) or fold < gate.min_od_fold:
        failures.append("dynamic_range")

    if rfu is not None:
        signal = require(rfu, RawRFU, name="rfu").minus_background(rfu_background).array
        metrics["max_reporter_above_background"] = float(np.nanmax(signal))
        if metrics["max_reporter_above_background"] < gate.min_reporter_above_background:
            failures.append("reporter_above_background")

    return WellQualityResult(
        well=well, failures=tuple(failures), metrics=metrics, gate=gate
    )


def assess_plate(
    od: RawOD,
    od_blank: float,
    rfu: RawRFU | None = None,
    rfu_background: float = 0.0,
    gate: OpticalQualityGate = OpticalQualityGate(),
) -> pd.DataFrame:
    """Apply :func:`assess_well` to every column of an OD frame.

    Args:
        od: Time-indexed RAW OD, one column per well, wrapped in
            :class:`~ystwin.readings.RawOD`. The frame keeps its columns and index; the
            reading type carries only whether the blank has come off.
        od_blank: Media-only OD.
        rfu: Optional matching reporter frame.
        rfu_background: Reporter baseline.
        gate: Thresholds to apply.
    """
    frame = require(od, RawOD, name="od").values
    rfu_frame = None if rfu is None else require(rfu, RawRFU, name="rfu").values
    t = frame.index.to_numpy(dtype=float)
    rows = []
    for well in frame.columns:
        result = assess_well(
            t, RawOD(frame[well].to_numpy()), od_blank=od_blank,
            rfu=None if rfu_frame is None else RawRFU(rfu_frame[well].to_numpy()),
            rfu_background=rfu_background, gate=gate, well=str(well),
        )
        rows.append({
            "well": result.well,
            "passed": result.passed,
            "failures": ",".join(result.failures),
            **result.metrics,
        })
    return pd.DataFrame(rows)
