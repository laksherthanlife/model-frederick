"""Fitting the panel's dose parameters to real plates, and reporting what survives.

The panel's EC50s are literature estimates. Where a plate carries a dose ladder they are
testable, so this fits the induction EC50 and the lethal dose from measured activity and
compares them with what the panel encodes.

Only what is identifiable is fitted. One reporter under one agent gives the shape of the
curve, which fixes both doses; the amplitude is a product of promoter strength, reporter
gain and the target weight, and no single ladder separates them. It is reported as
confounded rather than dressed up as a target weight.

Wells that stopped growing are excluded. Negative blank-corrected activity is an absent
measurement, not a small one, and fitting it would drag the lethal dose toward the noise.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from .stress_panel import STRESSORS, _LETHAL_HILL

__all__ = [
    "CalibrationVerdict",
    "DoseFit",
    "ReplicateVerdict",
    "calibrate_replicates",
    "calibrate_stressor",
    "fit_dose_response",
    "half_maximal_dose",
    "healthy_doses",
    "induction_is_identifiable",
    "lethal_dose_from_growth",
    "usable_points",
]

_MIN_GROWTH = 0.05
_MIN_POINTS = 4
_AGREEMENT_FOLD = 3.0
_MIN_FIT_QUALITY = 0.5
_HEALTHY_FRACTION = 0.7
_SATURATION_FRACTION = 0.15
_MIN_INDUCTION_FOLD = 1.05


def usable_points(doses, activity, growth) -> np.ndarray:
    """Which wells carry a measurement at all.

    A culture that did not grow has no promoter activity to report, and its
    blank-corrected signal can fall below zero. Keeping such a point as though it were a
    low reading would put the lethal dose wherever the noise happened to land.
    """
    doses, activity, growth = (np.asarray(x, dtype=float) for x in (doses, activity, growth))
    if not (len(doses) == len(activity) == len(growth)):
        raise ValueError("doses, activity and growth must be the same length")
    return (activity > 0.0) & (growth > _MIN_GROWTH)


def lethal_dose_from_growth(doses, growth) -> float | None:
    """Dose at which growth falls to half the undosed control's rate.

    Lethality is a growth measurement. Reading it off a reporter instead makes it depend on
    a model of transcription, on the dilution correction, and on whichever curve is being
    fitted -- and on these plates that route gave 7.2 mM for an agent that halves growth
    below 2. Optical density measures it directly, so optical density decides it.
    """
    doses, growth = np.asarray(doses, dtype=float), np.asarray(growth, dtype=float)
    order = np.argsort(doses)
    doses, growth = doses[order], growth[order]
    if doses[0] > 0:
        raise ValueError("need an undosed control well to measure growth against")

    control = float(growth[0])
    if control <= 0:
        return None
    target = 0.5 * control
    for i in range(1, len(growth)):
        if growth[i] <= target:
            span = growth[i - 1] - growth[i]
            frac = (growth[i - 1] - target) / span if span > 0 else 0.0
            return float(doses[i - 1] + frac * (doses[i] - doses[i - 1]))
    return None


def healthy_doses(doses, growth, fraction: float = _HEALTHY_FRACTION) -> np.ndarray:
    """Wells still growing near their control, where induction can be read at all.

    Below this a culture is not merely stressed but failing, and its promoter activity is a
    small number divided by another small number. The first calibration anchored a
    half-maximal estimate on a peak measured at 23% of control growth.
    """
    doses, growth = np.asarray(doses, dtype=float), np.asarray(growth, dtype=float)
    if not (doses == 0).any():
        raise ValueError("need an undosed control well to measure growth against")
    control = float(growth[doses == 0].mean())
    if control <= 0:
        return np.zeros(len(doses), dtype=bool)
    return growth >= fraction * control


def induction_is_identifiable(
    doses, activity, fraction: float = _SATURATION_FRACTION,
    min_fold: float = _MIN_INDUCTION_FOLD,
) -> bool:
    """Whether a real induction flattens before the ladder ends.

    An EC50 is the dose at half of a maximum, so a ladder still climbing at its last healthy
    dose has not shown the maximum and cannot place the half. Fitting one anyway returns
    whatever the model's curvature happens to prefer.

    A flat ladder passes the flattening test trivially, so an induction has to be there
    first. ``min_fold`` is the smallest rise a plate distinguishes from well-to-well
    variation, which these readers put at a few percent.
    """
    doses, activity = np.asarray(doses, dtype=float), np.asarray(activity, dtype=float)
    order = np.argsort(doses)
    doses, activity = doses[order], activity[order]
    if len(doses) < 3:
        return False
    basal = float(activity[0])
    rise = float(activity.max() - basal)
    if rise <= 0 or basal <= 0 or activity.max() / basal < min_fold:
        return False
    final_step = float(activity[-1] - activity[-2])
    return final_step <= fraction * rise


@dataclass(frozen=True)
class DoseFit:
    """A biphasic dose response fitted to one ladder."""

    ec50: float
    lethal_dose: float
    amplitude: float
    basal: float
    r_squared: float
    n_points: int
    identifiable: bool
    unidentified_reason: str = ""
    amplitude_is_confounded: bool = True


def _model(dose, ec50, lethal, amplitude, basal):
    """Basal plus induction, the whole of it scaled by viability.

    An uninduced promoter is not silent, so a curve through the origin cannot describe a
    real ladder. Viability multiplies the sum rather than the induced part alone, because a
    cell that has stopped transcribing has stopped transcribing the basal too.
    """
    d = np.asarray(dose, dtype=float)
    induction = d / (ec50 + d)
    return (basal + amplitude * induction) / (1.0 + (d / lethal) ** _LETHAL_HILL)


def fit_dose_response(doses, activity) -> DoseFit:
    """Fit induction EC50, lethal dose and a confounded amplitude to one ladder."""
    doses, activity = np.asarray(doses, dtype=float), np.asarray(activity, dtype=float)
    if len(doses) < _MIN_POINTS:
        raise ValueError(f"need at least {_MIN_POINTS} usable doses, got {len(doses)}")

    positive = doses[doses > 0]
    highest = float(positive.max())
    span = float(np.max(activity) - np.min(activity)) or 1.0
    floor = float(activity[doses == 0].mean()) if (doses == 0).any() else float(np.min(activity))
    guess = [float(np.median(positive)), float(positive.max()), span, max(floor, 0.0)]
    ec50_ceiling, lethal_ceiling = highest * 50, highest * 500
    bounds = ([1e-4, 1e-3, 1e-9, 0.0], [ec50_ceiling, lethal_ceiling, np.inf, np.inf])
    try:
        (ec50, lethal, amplitude, basal), _ = curve_fit(
            _model, doses, activity, p0=guess, bounds=bounds, maxfev=40000)
    except RuntimeError:
        ec50, lethal, amplitude, basal = guess

    residual = activity - _model(doses, ec50, lethal, amplitude, basal)
    centred = activity - activity.mean()
    total = float(centred @ centred)
    r2 = 1.0 - float(residual @ residual) / total if total > 0 else float("nan")

    reason = ""
    if lethal >= lethal_ceiling * 0.99 or ec50 >= ec50_ceiling * 0.99:
        reason = ("the ladder has no falling limb, so the lethal dose is unconstrained and "
                  "the optimiser ran to its bound; dose above the peak to settle it")
    elif lethal <= ec50:
        reason = (f"the fitted lethal dose {lethal:.3g} is below the inducing {ec50:.3g}, "
                  "which no cell does; the ladder does not separate the two limbs")
    elif float(np.max(doses)) < lethal * 0.5:
        reason = ("the highest dose is well below the fitted lethal dose, so the falling "
                  "limb was never observed")
    return DoseFit(float(ec50), float(lethal), float(amplitude), float(basal), r2,
                   len(doses), not reason, reason)


def half_maximal_dose(doses, activity) -> float | None:
    """Dose where activity first crosses halfway from basal to its observed peak.

    A LOWER BOUND on the EC50, not an estimate of it. The peak it measures against is
    already pulled down by the toxicity that turns the curve over, so half of that peak is
    reached earlier than half of the induction the promoter would have reached unimpeded.
    Against known parameters it returns between 0.89 and 0.13 of the truth, worsening as
    toxicity bites nearer the inducing dose.

    It earns its place by assuming nothing about the shape, which matters when the fitted
    model is misspecified. It does not earn the right to be called an EC50.
    """
    doses, activity = np.asarray(doses, dtype=float), np.asarray(activity, dtype=float)
    order = np.argsort(doses)
    doses, activity = doses[order], activity[order]
    peak = int(np.argmax(activity))
    if peak == 0:
        return None
    basal, top = float(activity[0]), float(activity[peak])
    if top <= basal:
        return None
    target = basal + 0.5 * (top - basal)
    rise_d, rise_a = doses[: peak + 1], activity[: peak + 1]
    for i in range(1, len(rise_a)):
        if rise_a[i] >= target:
            span = rise_a[i] - rise_a[i - 1]
            frac = (target - rise_a[i - 1]) / span if span > 0 else 0.0
            return float(rise_d[i - 1] + frac * (rise_d[i] - rise_d[i - 1]))
    return None


@dataclass(frozen=True)
class CalibrationVerdict:
    """What a plate says about a dose parameter the panel encodes."""

    stressor: str
    verdict: str
    encoded_ec50: float
    fitted_ec50: float | None
    fitted_lethal_dose: float | None
    r_squared: float | None
    n_usable: int
    note: str
    half_maximal_dose: float | None = None


def calibrate_stressor(stressor: str, doses, activity, growth) -> CalibrationVerdict:
    """Test the panel's encoded EC50 for one agent against a measured ladder.

    Args:
        stressor: Agent whose encoded parameters are being tested.
        doses: Applied concentrations, in the agent's own units.
        activity: Dilution-corrected promoter activity per dose.
        growth: Specific growth rate per dose, used to drop dead wells.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    encoded = STRESSORS[stressor].ec50

    doses_all, activity_all = np.asarray(doses, dtype=float), np.asarray(activity, dtype=float)
    keep = usable_points(doses, activity, growth) & healthy_doses(doses, growth)
    doses, activity = doses_all[keep], activity_all[keep]
    if len(doses) < _MIN_POINTS:
        return CalibrationVerdict(
            stressor, "INCONCLUSIVE", encoded, None, None, None, int(keep.sum()),
            f"only {int(keep.sum())} doses left a healthy culture, below the {_MIN_POINTS} "
            "a curve needs; the rest were too slowed for activity to be read")
    crossing = half_maximal_dose(doses, activity)
    if not induction_is_identifiable(doses, activity):
        return CalibrationVerdict(
            stressor, "INCONCLUSIVE", encoded, None, None, None, int(keep.sum()),
            "induction is still climbing at the last dose the culture tolerates, so no "
            "maximum was reached and no half of one can be placed; the crossing at "
            f"{crossing:.3g} is a lower bound, not an estimate"
            if crossing else "induction never rose above plate variation", crossing)

    fit = fit_dose_response(doses, activity)
    if not fit.identifiable:
        return CalibrationVerdict(
            stressor, "INCONCLUSIVE", encoded, fit.ec50, fit.lethal_dose, fit.r_squared,
            int(keep.sum()), fit.unidentified_reason, crossing)
    if not np.isfinite(fit.r_squared) or fit.r_squared < _MIN_FIT_QUALITY:
        return CalibrationVerdict(
            stressor, "INCONCLUSIVE", encoded, fit.ec50, fit.lethal_dose, fit.r_squared,
            int(keep.sum()),
            f"the curve fit poorly (R2 {fit.r_squared:.2f}), so it says nothing about the "
            f"encoded EC50 whatever value it landed on", crossing)

    fold = max(fit.ec50 / encoded, encoded / fit.ec50)
    if fold <= _AGREEMENT_FOLD:
        verdict = "CONFIRMED"
        note = f"fitted EC50 {fit.ec50:.3g} is within {fold:.1f}-fold of the encoded {encoded:.3g}"
    else:
        verdict = "REFUTED"
        note = f"fitted EC50 {fit.ec50:.3g} is {fold:.1f}-fold from the encoded {encoded:.3g}"
    return CalibrationVerdict(
        stressor, verdict, encoded, fit.ec50, fit.lethal_dose, fit.r_squared,
        int(keep.sum()), note, crossing)


@dataclass(frozen=True)
class ReplicateVerdict:
    """What independent replicates jointly say about an encoded parameter."""

    stressor: str
    verdict: str
    encoded: float
    estimates: list[float]
    low: float
    high: float
    note: str


def calibrate_replicates(
    stressor: str, doses, activity_by_replicate: dict, growth_by_replicate: dict,
    encoded: float | None = None, spread_limit: float = 3.0,
) -> ReplicateVerdict:
    """Fit each replicate separately and let their spread decide.

    A point estimate from one pooled curve cannot say whether an encoded value is excluded,
    and pooling hides replicates disagreeing with each other. Fitting them separately gives
    an interval, and the encoded value is refuted only when it falls outside it.

    Args:
        stressor: Agent whose encoded parameter is being tested.
        doses: Applied concentrations, shared across replicates.
        activity_by_replicate: Dilution-corrected activity per replicate name.
        growth_by_replicate: Matching growth rates, used to drop dead wells.
        encoded: Value under test; defaults to what the panel carries.
        spread_limit: Largest fold-range across replicates that still permits a verdict.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    encoded = STRESSORS[stressor].ec50 if encoded is None else encoded

    estimates = []
    for name, activity in activity_by_replicate.items():
        growth = growth_by_replicate[name]
        keep = usable_points(doses, activity, growth)
        d, a = np.asarray(doses, float)[keep], np.asarray(activity, float)[keep]
        if len(d) < _MIN_POINTS:
            continue
        fit = fit_dose_response(d, a)
        if fit.identifiable and np.isfinite(fit.r_squared) and fit.r_squared >= _MIN_FIT_QUALITY:
            estimates.append(float(fit.ec50))

    if len(estimates) < 2:
        return ReplicateVerdict(
            stressor, "INCONCLUSIVE", encoded, estimates, float("nan"), float("nan"),
            f"{len(estimates)} replicate(s) gave an identifiable fit; one replicate cannot "
            "bound anything")

    low, high = float(np.min(estimates)), float(np.max(estimates))
    if high / max(low, 1e-12) > spread_limit:
        return ReplicateVerdict(
            stressor, "INCONCLUSIVE", encoded, estimates, low, high,
            f"replicates span {low:.3g} to {high:.3g}, a {high / low:.1f}-fold "
            "disagreement among themselves, so they cannot exclude a third value")
    if low <= encoded <= high:
        return ReplicateVerdict(
            stressor, "CONFIRMED", encoded, estimates, low, high,
            f"the encoded {encoded:.3g} lies inside the replicate range {low:.3g}-{high:.3g}")
    return ReplicateVerdict(
        stressor, "REFUTED", encoded, estimates, low, high,
        f"every replicate fell {'below' if high < encoded else 'above'} the encoded "
        f"{encoded:.3g}; they span {low:.3g}-{high:.3g}")
