"""Acceptance gates for simulated panel data, which had none.

Three grounds for refusing a simulated plate: it could not have been measured, its summary
statistics sit outside a real plate's spread, or a batch of them is a set of near-copies.

The posterior-predictive gate is the one that earns its keep. It is the only check here
that compares simulated data against measured data rather than against the generator's own
assumptions, and a generator can be perfectly self-consistent while producing activity
distributions no culture ever made.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "GateReport",
    "accept_panel",
    "diversity_ok",
    "physically_possible",
    "resembles_reference",
]

_CEILING = 1e6
_MIN_INDUCTION = 1e-9
_MIN_BATCH_CV = 0.01


def physically_possible(readings, reporters=None, ceiling: float = _CEILING) -> str | None:
    """Why a set of readings could not have come off a plate, or None if it could.

    No channel may read below zero. A promoter cannot produce less than nothing, and a
    ratiometric sensor reports a ratio that moves either side of a set point it never
    reaches -- both sit on a constitutive floor. A negative reading therefore means either
    a blank subtraction gone wrong or a perturbation the generator let exceed the floor,
    and both are worth refusing.
    """
    values = np.asarray(readings, dtype=float)
    if values.size == 0:
        return "physical: no readings at all"
    if not np.all(np.isfinite(values)):
        return "physical: non-finite readings"

    if values.min() < 0:
        return "physical: a reading below zero, which is under the constitutive floor"
    if np.abs(values).max() > ceiling:
        return f"physical: a reading above the reader ceiling of {ceiling:g}"
    if np.abs(values).max() <= _MIN_INDUCTION:
        return "physical: no induction anywhere on the plate"
    return None


def _shape_features(values: np.ndarray) -> np.ndarray:
    """Scale-free summary, because reporter gain is arbitrary between builds."""
    values = np.asarray(values, dtype=float)
    centre = float(np.nanmedian(values))
    if abs(centre) < 1e-12:
        return np.array([0.0, 0.0, 0.0])
    scaled = values / centre
    return np.array([float(np.nanstd(scaled)), float(np.nanmax(scaled)),
                     float(np.nanpercentile(scaled, 90) - np.nanpercentile(scaled, 10))])


def resembles_reference(readings, reference, tolerance: float = 0.5) -> bool:
    """Whether simulated readings have the spread of a real plate.

    Compared on scale-free features, since a reporter's gain says nothing about whether the
    biology behind it is right. What is being asked is whether the distribution has the
    shape measured data has, not whether it lands on the same numbers.
    """
    target = _shape_features(reference)
    got = _shape_features(readings)
    scale = np.where(np.abs(target) > 1e-12, np.abs(target), 1.0)
    return bool(np.all(np.abs(got - target) / scale <= tolerance))


def diversity_ok(batch, min_cv: float = _MIN_BATCH_CV) -> bool:
    """Whether a batch carries more information than one of its members.

    A batch of near-copies passes every other check while being worth a single plate.
    """
    batch = [np.asarray(b, dtype=float) for b in batch]
    if len(batch) < 2:
        return False
    summaries = np.array([np.nanmean(b) for b in batch])
    centre = float(np.abs(summaries).mean())
    if centre < 1e-12:
        return False
    return bool(float(np.std(summaries)) / centre >= min_cv)


@dataclass(frozen=True)
class GateReport:
    """Which gates ran, which were skipped, and which one refused."""

    accepted: bool
    refused_by: str | None = None
    reason: str = ""
    checked: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)


def accept_panel(readings, reference=None, batch=None, reporters=None,
                 tolerance: float = 0.5) -> GateReport:
    """Run every applicable gate and report the first refusal.

    Args:
        readings: Simulated panel readings under test.
        reference: Measured readings to compare the spread against. Omitted means the
            posterior-predictive gate is skipped, and the report says so rather than
            quietly passing.
        batch: Sibling draws, for the diversity check.
        reporters: Channel names. Accepted and currently unused: both sensor kinds sit
            on a constitutive floor, so :func:`physically_possible` applies one sign
            rule to every channel.
        tolerance: Allowed relative departure of each summary feature.
    """
    checked, skipped = [], []

    checked.append("physical")
    failure = physically_possible(readings, reporters)
    if failure:
        return GateReport(False, "physical", failure, checked, skipped)

    if reference is None:
        skipped.append("posterior-predictive")
    else:
        checked.append("posterior-predictive")
        if not resembles_reference(readings, reference, tolerance):
            return GateReport(
                False, "posterior-predictive",
                "summary spread sits outside the reference plate's own", checked, skipped)

    if batch is None:
        skipped.append("diversity")
    else:
        checked.append("diversity")
        if not diversity_ok(batch):
            return GateReport(False, "diversity", "the batch is a set of near-copies",
                              checked, skipped)

    return GateReport(True, None, "", checked, skipped)
