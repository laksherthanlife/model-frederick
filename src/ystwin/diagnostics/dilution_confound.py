"""D2: how much of a reporter's apparent induction is growth dilution?

At quasi-steady state log S = log k_synth - log(mu + k_deg), so regressing log S on
-log(mu + k_deg) measures how much of the reporter's dynamic range is the culture slowing rather
than the promoter turning on. Run before fitting any latent state to a reporter channel.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..growth import specific_growth_rate
from ..plate.layout import recorded_well_roles
from ..reporter import (
    ReporterKinetics,
    default_activity_window_h,
    _smooth_derivative,
    naive_specific_fluorescence,
    promoter_activity,
)
from ..readings import CorrectedOD, CorrectedRFU, RawOD, RawRFU, require


__all__ = [
    "DilutionConfoundReport",
    "PlateDilutionReport",
    "detect_blank_wells",
    "dilution_confound_report",
    "plate_dilution_report",
]

_MIN_GROWTH_FOLD = 1.15
_MU_FLOOR = 1e-3
_MIN_DILUTION_SPREAD = 1e-6
_MODEL_INADEQUATE_FRACTION = 0.25
_MIN_REGRESSION_POINTS = 5


@dataclass(frozen=True)
class DilutionConfoundReport:
    """Per-well verdict on whether a reporter channel is measuring 1/mu."""

    times_h: np.ndarray = field(repr=False)
    growth_rate: np.ndarray = field(repr=False)
    naive_specific: np.ndarray = field(repr=False)
    promoter_activity: np.ndarray = field(repr=False)
    dilution_explained_r2: float
    dilution_slope: float
    naive_fold_change: float
    activity_fold_change: float
    negative_activity_fraction: float
    implied_min_k_deg: float
    assumed_k_deg: float
    qss_fraction: float
    mu_max: float
    mu_min: float

    @property
    def excess_loss_rate(self) -> float:
        """How far the assumed loss rate falls short of what the decline requires (1/h)."""
        return max(0.0, self.implied_min_k_deg - self.assumed_k_deg)

    @property
    def verdict(self) -> str:
        """Coarse label; the numbers, not the label, are the deliverable."""
        # Order matters. negative_activity_fraction comes straight from
        # promoter_activity and reads neither the R^2 nor the QSS fraction, so it is
        # determinable precisely when they are not -- and a plate where most inferred
        # activities are negative has falsified the assumed kinetics whether or not the
        # dilution regression was estimable. Checking it first.
        if self.negative_activity_fraction > _MODEL_INADEQUATE_FRACTION:
            return "model-inadequate"
        # The remaining labels all read the R^2 or the fold change, which need
        # quasi-steady state: mu + k_deg > 2/duration, so mu > 0.042 /h over 48 h. Below
        # that nothing was measured, and "mixed" would report a finding where there is
        # none.
        if self.qss_fraction <= 0.0 or not np.isfinite(self.dilution_explained_r2):
            return "not-estimable"
        if self.activity_fold_change >= 2.0 and self.activity_fold_change >= 0.5 * self.naive_fold_change:
            return "promoter-driven"
        if self.dilution_explained_r2 > 0.75 and self.activity_fold_change < 1.5:
            return "dilution-dominated"
        return "mixed"

    def summary(self) -> str:
        return (
            f"naive RFU/OD fold={self.naive_fold_change:.2f}  "
            f"promoter-activity fold={self.activity_fold_change:.2f}  "
            f"R^2(log S ~ -log(mu+k_deg))={self.dilution_explained_r2:.3f} "
            f"[QSS {self.qss_fraction:.0%}]  "
            f"slope={self.dilution_slope:.2f}  "
            f"neg-activity={self.negative_activity_fraction:.2f}  "
            f"implied min k_deg={self.implied_min_k_deg:.3f}/h "
            f"(excess {self.excess_loss_rate:.3f})  "
            f"mu {self.mu_min:.3f}-{self.mu_max:.3f} /h  -> {self.verdict}"
        )


def _robust_fold_change(series: np.ndarray, frac: float = 0.1) -> float:
    """End-to-start ratio from window medians. Undefined unless both are positive."""
    n = max(3, int(len(series) * frac))
    start = float(np.median(series[:n]))
    end = float(np.median(series[-n:]))
    if start <= 0 or end <= 0 or not np.isfinite(start) or not np.isfinite(end):
        return float("nan")
    return end / start


def dilution_confound_report(
    times_h: np.ndarray,
    optical_density: RawOD,
    rfu: RawRFU,
    kinetics: ReporterKinetics = ReporterKinetics(),
    od_blank: float = 0.0,
    rfu_background: float = 0.0,
    growth_window_h: float = 1.5,
    activity_window_h: float | None = None,
    qss_max_relaxation_h: float | None = None,
) -> DilutionConfoundReport:
    """Quantify the growth-dilution contribution to one well's reporter trace.

    Args:
        times_h: Ascending time grid, hours.
        optical_density: Raw OD600 for the well.
        rfu: Raw fluorescence for the well, same grid.
        kinetics: Reporter constants (a stable FP has ``k_deg`` near 0).
        od_blank: Media-only OD to subtract.
        rfu_background: Media plus autofluorescence baseline to subtract.
        growth_window_h: Smoothing window for the growth-rate estimate.
        activity_window_h: Smoothing window for the promoter-activity inversion.
            Defaults to :func:`~ystwin.reporter.default_activity_window_h`.
        qss_max_relaxation_h: Longest reporter relaxation time still treated as
            quasi-steady. Defaults to half the run duration.
    """
    t = np.asarray(times_h, dtype=float)
    density = require(optical_density, RawOD, name="optical_density").minus_blank(od_blank)
    corrected_rfu = require(rfu, RawRFU, name="rfu").minus_background(rfu_background)
    od = density.array

    if np.any(od <= 0):
        raise ValueError("blank-corrected OD must be positive; check od_blank")
    # Judge growth from robust start/end levels; max/min is noise-dominated.
    od_fold = _robust_fold_change(od)
    if not np.isfinite(od_fold) or od_fold < _MIN_GROWTH_FOLD:
        raise ValueError(
            f"well did not grow (robust OD fold change {od_fold:.3f} < {_MIN_GROWTH_FOLD}); "
            "growth rate and therefore dilution are not identifiable here"
        )

    if activity_window_h is None:
        activity_window_h = default_activity_window_h(t)
    mu = specific_growth_rate(t, density, window_h=growth_window_h)
    # Typed end to end: corrected fluorescence over corrected density gives a
    # SpecificFluorescence, which is exactly what promoter_activity takes. Nothing here
    # unwraps to a bare array until the arithmetic below needs one.
    ratio = naive_specific_fluorescence(corrected_rfu, density)
    activity = promoter_activity(t, ratio, mu, kinetics, window_h=activity_window_h)
    specific = ratio.array

    # Quasi-steady state needs the reporter to relax within the run.
    duration = float(t[-1] - t[0])
    max_relaxation_h = qss_max_relaxation_h if qss_max_relaxation_h is not None else duration / 2.0
    loss = mu + kinetics.k_deg
    qss = loss > (1.0 / max_relaxation_h)
    qss_fraction = float(np.mean(qss))

    usable = qss & (specific > 0) & np.isfinite(loss)
    if usable.sum() < _MIN_REGRESSION_POINTS:
        slope, r2 = float("nan"), float("nan")
    else:
        x = -np.log(loss[usable])
        y = np.log(specific[usable])
        # A flat growth rate leaves the regressor with no variance, so no slope is
        # identifiable.
        if np.std(x) < _MIN_DILUTION_SPREAD:
            slope, r2 = float("nan"), 0.0
        else:
            slope, intercept = np.polyfit(x, y, 1)
            resid = y - (slope * x + intercept)
            ss_tot = float(np.sum((y - y.mean()) ** 2))
            r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else 0.0
            r2 = float(np.clip(r2, 0.0, 1.0))

    # k_synth >= 0 implies k_deg >= max_t[-d(ln R)/dt - mu].
    positive_r = specific > 0
    implied = np.full(t.shape, -np.inf)
    if positive_r.sum() >= _MIN_REGRESSION_POINTS:
        dln = _smooth_derivative(t[positive_r], np.log(specific[positive_r]),
                                 activity_window_h, polyorder=3)
        implied[positive_r] = -dln - mu[positive_r]
    finite = implied[np.isfinite(implied)]
    implied_min_k_deg = float(max(0.0, np.nanmax(finite))) if finite.size else 0.0
    negative_fraction = float(np.mean(activity < 0))

    return DilutionConfoundReport(
        times_h=t,
        growth_rate=mu,
        naive_specific=specific,
        promoter_activity=activity,
        dilution_explained_r2=float(r2),
        dilution_slope=float(slope),
        naive_fold_change=_robust_fold_change(specific),
        activity_fold_change=_robust_fold_change(activity),
        negative_activity_fraction=negative_fraction,
        implied_min_k_deg=implied_min_k_deg,
        assumed_k_deg=float(kinetics.k_deg),
        qss_fraction=qss_fraction,
        mu_max=float(np.nanmax(mu)),
        mu_min=float(np.nanmin(mu)),
    )


@dataclass(frozen=True)
class PlateDilutionReport:
    """Plate-level roll-up of D2, with blanks and exclusions made explicit."""

    per_well: pd.DataFrame
    excluded: pd.DataFrame
    blank_wells: list[str]
    od_blank: float
    rfu_background: float

    @property
    def median_dilution_r2(self) -> float:
        return float(self.per_well["dilution_r2"].median())

    def summary(self) -> str:
        counts = self.per_well["verdict"].value_counts().to_dict()
        return (
            f"{len(self.per_well)} cultures analysed, {len(self.excluded)} excluded, "
            f"{len(self.blank_wells)} media blanks; verdicts {counts}; "
            f"median dilution R^2 {self.median_dilution_r2:.3f}"
        )


def detect_blank_wells(
    od: RawOD | CorrectedOD,
    tolerance: float = 0.02,
    max_growth_fold: float = 1.05,
    rfu: RawRFU | CorrectedRFU | None = None,
    max_reporter_fold: float = 1.15,
) -> list[str]:
    """Media-only wells: near the plate floor *and* not growing.

    Proximity to the floor alone is not sufficient. A culture held back by a high
    stressor dose also sits near the floor, and folding one into the blank raises
    it -- on a real plate that shifted the blank by 0.014 OD and pushed every
    derived growth rate above the physiological maximum for yeast. A media blank
    contains no biomass, so the definitive test is that it does not grow.

    Supplying the reporter channel makes the test stronger still. A medium well's
    OD can drift up through evaporation or settling debris -- on one real plate two
    of them rose 1.10x, which an OD-only test reads as growth -- but with no cells
    the fluorescence stays flat while a real culture's climbs several fold. When
    ``rfu`` is given, a well may drift in OD provided its reporter does not.

    Args:
        od: Time-indexed OD, one column per well. **Either correction state is accepted,
            and that is not laxity.** This function is looking FOR the blank, so demanding
            a corrected frame would require the answer as an input; and it compares wells
            within one frame -- a plate floor and an end-to-start ratio -- both of which a
            constant offset leaves unchanged. It is the one place in the package where the
            state genuinely does not matter, which is why it has to say so out loud rather
            than take a bare array like everything used to.
        tolerance: How far above the plate floor a blank may sit.
        max_growth_fold: Largest OD end-to-start ratio counted as flat.
        rfu: Optional matching reporter frame.
        max_reporter_fold: Largest reporter end-to-start ratio counted as flat.
    """
    frame = require(od, RawOD, CorrectedOD, name="od").values
    rfu_frame = (None if rfu is None
                 else require(rfu, RawRFU, CorrectedRFU, name="rfu").values)
    roles = recorded_well_roles(frame, rfu_frame)
    if roles is not None:
        return sorted(well for well in frame.columns if roles.get(well) == "blank")
    floor = float(frame.to_numpy().min())
    medians = frame.median(axis=0)
    near_floor = medians[medians <= floor + tolerance].index

    blanks = []
    for well in near_floor:
        od_flat = _robust_fold_change(frame[well].to_numpy()) <= max_growth_fold
        if rfu_frame is None or well not in rfu_frame.columns:
            if od_flat:
                blanks.append(well)
            continue
        # With the reporter available it is the decisive channel -- cells make
        # fluorescence and medium does not -- so OD drift alone is allowed.
        if _robust_fold_change(rfu_frame[well].to_numpy()) <= max_reporter_fold:
            blanks.append(well)
    return sorted(blanks)


def plate_dilution_report(
    aligned: pd.DataFrame,
    od_channel: str,
    reporter_channel: str,
    kinetics: ReporterKinetics = ReporterKinetics(),
    blank_wells: list[str] | None = None,
    **report_kwargs,
) -> PlateDilutionReport:
    """Run D2 across every culture on one aligned plate.

    Args:
        aligned: Frame from :meth:`SynergyRun.aligned`, columns ``(channel, well)``.
        od_channel: Channel label carrying optical density.
        reporter_channel: Channel label carrying the reporter.
        kinetics: Assumed reporter constants.
        blank_wells: Media-only wells. Detected from the OD floor when omitted.
        **report_kwargs: Forwarded to :func:`dilution_confound_report`.
    """
    for name in (od_channel, reporter_channel):
        if name not in set(aligned.columns.get_level_values("channel")):
            raise KeyError(
                f"channel {name!r} not on this plate; "
                f"have {sorted(set(aligned.columns.get_level_values('channel')))}"
            )
    t = aligned.index.to_numpy(dtype=float)
    od, rfu = aligned[od_channel], aligned[reporter_channel]

    # Channels may cover different wells; analyse only the shared ones.
    shared = [w for w in od.columns if w in set(rfu.columns)]
    if not shared:
        raise ValueError(
            f"no well appears in both {od_channel!r} and {reporter_channel!r}; "
            "check that the two sheets describe the same plate"
        )
    unpaired = [
        {"well": w, "reason": f"not read in the {reporter_channel!r} channel"}
        for w in od.columns if w not in set(rfu.columns)
    ]
    od, rfu = od[shared], rfu[shared]
    roles = recorded_well_roles(od, rfu)
    if roles is not None and blank_wells is not None:
        conflicts = [well for well in blank_wells if roles.get(well) == "culture"]
        if conflicts:
            raise ValueError(f"supplied blanks are recorded as culture wells: {conflicts}")

    if blank_wells is None:
        blank_wells = detect_blank_wells(RawOD(od), rfu=RawRFU(rfu))
        if not blank_wells:
            raise ValueError(
                "no media-blank well found in both channels. Recorded well identities "
                "are not replaced by an OD heuristic; supply matching blank measurements "
                "or declare blank_wells explicitly for an unrecorded plate."
            )
    blank_wells = [w for w in blank_wells if w in set(shared)]
    if not blank_wells:
        raise ValueError(
            "none of the supplied blank wells are read in both channels; "
            "a blank must be measured on the same channels as the cultures"
        )
    od_blank = float(od[blank_wells].to_numpy().mean())
    rfu_background = float(rfu[blank_wells].to_numpy().mean())

    rows, dropped = [], list(unpaired)
    for well in od.columns:
        if well in blank_wells:
            continue
        if roles is not None and roles.get(well) != "culture":
            dropped.append({"well": well, "reason": "not recorded as a culture"})
            continue
        try:
            r = dilution_confound_report(
                t, RawOD(od[well].to_numpy()), RawRFU(rfu[well].to_numpy()),
                kinetics=kinetics, od_blank=od_blank, rfu_background=rfu_background,
                **report_kwargs,
            )
        except ValueError as exc:
            dropped.append({"well": well, "reason": str(exc)})
            continue
        rows.append({
            "well": well,
            "dilution_r2": r.dilution_explained_r2,
            "dilution_slope": r.dilution_slope,
            "qss_fraction": r.qss_fraction,
            "naive_fold": r.naive_fold_change,
            "activity_fold": r.activity_fold_change,
            "negative_activity_fraction": r.negative_activity_fraction,
            "implied_min_k_deg": r.implied_min_k_deg,
            "excess_loss_rate": r.excess_loss_rate,
            "mu_min": r.mu_min,
            "mu_max": r.mu_max,
            "verdict": r.verdict,
        })
    return PlateDilutionReport(
        per_well=pd.DataFrame(rows),
        excluded=pd.DataFrame(dropped, columns=["well", "reason"]),
        blank_wells=blank_wells,
        od_blank=od_blank,
        rfu_background=rfu_background,
    )
