"""S1: fit the generator to an uploaded plate.

Plate terms, then growth, then promoter -- each depends on the one before. What the plate cannot
identify is reported as such: gain and gdcw_per_od appear only as a product, so only the product
is returned unless a measured dry-weight factor is supplied.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

from ..growth import max_specific_growth_rate, specific_growth_rate
from ..observation import ReporterOptics
from ..plate.layout import NEWPROTOCOL_LAYOUT, PlateLayout
from ..qpcr import STRESSOR_FOR_CONSTRUCT
from ..reporter import ReporterKinetics, promoter_activity
from .culture import CultureParameters
from .plate import DEFAULT_DOSES, PlateConditions
from ..readings import CorrectedOD, RawOD, RawRFU, SpecificFluorescence, require

__all__ = ["CalibrationResult", "calibrate_from_plate"]


@dataclass(frozen=True)
class CalibrationResult:
    """A generator fitted to one plate, with its own limits attached."""

    panel: dict[str, CultureParameters]
    conditions: PlateConditions
    gain_times_gdcw: float
    unidentifiable: tuple[str, ...]
    fit_quality: dict[str, float]
    per_condition: pd.DataFrame = field(repr=False, default_factory=pd.DataFrame)

    def summary(self) -> str:
        lines = [
            f"blank OD {self.conditions.od_blank:.4f}, background "
            f"{self.conditions.optics.background:.0f} RFU, reader CV "
            f"{self.conditions.reader_cv:.1%}",
            f"gain x gdcw_per_od = {self.gain_times_gdcw:.3g}"
            + (f"  (unidentifiable separately: {', '.join(self.unidentifiable)})"
               if self.unidentifiable else ""),
        ]
        for construct, params in self.panel.items():
            lines.append(
                f"  {construct:<12} mu_max {params.mu_max:.3f}/h  "
                f"growth IC50 {params.growth_ic50:.2f} mM  "
                f"induction {params.promoter_peak / params.promoter_basal:.2f}x  "
                f"EC50 {params.promoter_ec50:.2f} mM  (R^2 {self.fit_quality[construct]:.2f})"
            )
        return "\n".join(lines)


def _estimate_reader_cv(od: pd.DataFrame, window: int = 7, polyorder: int = 2) -> float:
    """Reader noise from the high-frequency residual of each well's own trace.

    Replicate scatter cannot supply this on its own: it also contains well-to-well
    variation in starting biomass, which is a fixed offset per well rather than an
    independent error per reading. Detrending each well separates them -- what is
    left after a smooth is subtracted is the per-reading noise.
    """
    scatter = []
    for column in od.columns:
        trace = od[column].to_numpy(dtype=float)
        if trace.size < window or not np.all(trace > 0):
            continue
        length = min(window if window % 2 else window + 1, trace.size - (1 - trace.size % 2))
        if length <= polyorder:
            continue
        smooth = savgol_filter(trace, window_length=length, polyorder=polyorder)
        usable = smooth > 0
        if usable.sum() < 3:
            continue
        scatter.append(float(np.std((trace[usable] - smooth[usable]) / smooth[usable])))
    return float(np.median(scatter)) if scatter else 0.01


def _hill_inhibition(dose, top, ic50, hill):
    return top / (1.0 + (np.maximum(dose, 0.0) / ic50) ** hill)


def _hill_induction(dose, basal, peak, ec50, hill):
    d = np.maximum(dose, 0.0) ** hill
    return basal + (peak - basal) * d / (ec50**hill + d)


def _fit(model, x, y, p0, bounds):
    try:
        params, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=40000)
    except (RuntimeError, ValueError):
        return np.asarray(p0, dtype=float), 0.0
    residual = y - model(x, *params)
    total = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - float(np.sum(residual**2)) / total if total > 0 else 1.0
    return params, float(np.clip(r2, 0.0, 1.0))


def calibrate_from_plate(
    od: RawOD,
    rfu: RawRFU,
    blank_wells,
    layout: PlateLayout = NEWPROTOCOL_LAYOUT,
    doses: dict[str, tuple[float, ...]] | None = None,
    k_deg: float = 0.0,
    gdcw_per_od: float | None = None,
) -> CalibrationResult:
    """Fit generator parameters to an uploaded plate.

    Args:
        od: Time-indexed RAW OD, one column per well, as
            :class:`~ystwin.readings.RawOD`. Raw because the blank is what this function
            derives from ``blank_wells``; handing it a corrected frame would subtract the
            blank twice.
        rfu: Matching RAW reporter frame.
        blank_wells: Medium-only wells. Required -- without them the blank, the
            background and every ratio built on them are guesses.
        layout: Which construct and dose each well holds.
        doses: Dose ladder per stressor.
        k_deg: Reporter loss beyond dilution, from a chase if one has been run.
        gdcw_per_od: Measured dry weight per OD unit. Supplying it resolves the
            reporter gain, which is otherwise confounded with it.
    """
    # Unwrapped first, because `blank_wells` is filtered against the frame's columns and
    # a reading carries the correction state rather than the DataFrame's interface.
    od = require(od, RawOD, name="od").values
    rfu = require(rfu, RawRFU, name="rfu").values
    blank_wells = [w for w in blank_wells if w in od.columns]
    if not blank_wells:
        raise ValueError(
            "no blank well supplied or present; the blank, the reporter background "
            "and every ratio built on them cannot be fitted without one"
        )
    doses = doses or dict(DEFAULT_DOSES)
    t = od.index.to_numpy(dtype=float)

    od_blank = float(od[blank_wells].to_numpy().mean())
    background = float(rfu[blank_wells].to_numpy().mean())

    reader_cv = _estimate_reader_cv(od.drop(columns=blank_wells))

    rows, residuals = [], []
    for construct, columns in layout.construct_columns.items():
        stressor = STRESSOR_FOR_CONSTRUCT.get(construct, "DTT")
        ladder = doses.get(stressor, DEFAULT_DOSES["DTT"])
        for index, dose in enumerate(ladder):
            wells = [w for w in layout.wells(construct, index) if w in od.columns]
            if not wells:
                continue
            corrected = (od[wells] - od_blank)
            if not (corrected.to_numpy() > 0).all():
                continue
            mean_od = corrected.mean(axis=1).to_numpy()
            residuals.append(
                (corrected.std(axis=1) / corrected.mean(axis=1)).to_numpy()
            )
            specific = ((rfu[wells] - background) / corrected).mean(axis=1).to_numpy()
            mu = specific_growth_rate(t, CorrectedOD(mean_od))
            activity = promoter_activity(t, SpecificFluorescence(specific), mu, ReporterKinetics(k_deg=k_deg))
            interior = slice(max(2, len(t) // 5), -2)
            rows.append({
                "construct": construct, "stressor": stressor, "dose_mM": float(dose),
                "mu_max": max_specific_growth_rate(t, CorrectedOD(mean_od)),
                "promoter_activity": float(np.nanmedian(activity[interior])),
            })

    per_condition = pd.DataFrame(rows)
    # Replicate scatter mixes reader noise with well-to-well variation.
    combined = float(np.nanmedian(np.concatenate(residuals))) if residuals else reader_cv
    well_cv = float(np.sqrt(max(combined**2 - reader_cv**2, 0.0)))

    panel, quality = {}, {}
    for construct, group in per_condition.groupby("construct"):
        group = group.sort_values("dose_mM")
        dose = group.dose_mM.to_numpy()
        mu = np.nan_to_num(group.mu_max.to_numpy(), nan=0.0)
        activity = np.clip(np.nan_to_num(group.promoter_activity.to_numpy(), nan=0.0), 1e-12, None)

        (top, ic50, growth_hill), growth_r2 = _fit(
            _hill_inhibition, dose, mu,
            p0=[max(mu.max(), 1e-3), max(np.median(dose[dose > 0]), 0.5), 1.5],
            bounds=([1e-6, 1e-3, 0.3], [5.0, 1e3, 8.0]),
        )
        (basal, peak, ec50, induction_hill), promoter_r2 = _fit(
            _hill_induction, dose, activity,
            p0=[activity[0], max(activity.max(), activity[0] * 1.1),
                max(np.median(dose[dose > 0]), 0.5), 2.0],
            bounds=([1e-12, 1e-12, 1e-3, 0.3], [np.inf, np.inf, 1e3, 8.0]),
        )
        panel[construct] = CultureParameters(
            mu_max=float(top), growth_ic50=float(ic50), growth_hill=float(growth_hill),
            promoter_basal=float(basal), promoter_peak=float(peak),
            promoter_ec50=float(ec50), promoter_hill=float(induction_hill),
            carrying_capacity=float(max(od.to_numpy().max() * 6, 1.0)),
            k_deg=k_deg,
        )
        quality[construct] = float(np.mean([growth_r2, promoter_r2]))

    # The fitted "activity" is gdcw*gain*k_synth, so the scale carries the product.
    scale = float(np.nanmedian([p.promoter_basal for p in panel.values()])) or 1.0
    gain_times_gdcw = scale / 1.0e-3  # against the generator's nominal basal activity

    unidentifiable: tuple[str, ...] = ()
    if gdcw_per_od is None:
        unidentifiable = ("gain", "gdcw_per_od")
        effective_gdcw, gain = 0.42, gain_times_gdcw / 0.42
    else:
        effective_gdcw, gain = gdcw_per_od, gain_times_gdcw / gdcw_per_od

    # Express the panel in per-cell units now that a scale has been chosen.
    scaled = {
        construct: replace(
            params,
            promoter_basal=params.promoter_basal / gain_times_gdcw,
            promoter_peak=params.promoter_peak / gain_times_gdcw,
        )
        for construct, params in panel.items()
    }

    conditions = PlateConditions(
        layout=layout, blank_wells=tuple(blank_wells),
        duration_h=float(t[-1] - t[0]), n_timepoints=len(t), doses=doses,
        optics=ReporterOptics(
            gain=float(gain), background=background,
            autofluorescence=0.0, inner_filter_coeff=None,
        ),
        inoculum_od=float(np.median((od.drop(columns=blank_wells).iloc[0] - od_blank))),
        od_blank=od_blank, reader_cv=reader_cv, well_cv=well_cv,
        gdcw_per_od=float(effective_gdcw),
    )
    return CalibrationResult(
        panel=scaled, conditions=conditions, gain_times_gdcw=gain_times_gdcw,
        unidentifiable=unidentifiable, fit_quality=quality, per_condition=per_condition,
    )
