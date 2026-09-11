"""Where the reader stops responding linearly, measured from a dilution series.

`OpticalQualityGate.od_linear_max` is a placeholder and the gate says so: everything
downstream divides by optical density, so a growth rate taken above the linear range is
wrong by whatever the reader is compressing, and no later modelling recovers it.

A dilution series answers it directly. The dilution factor says what should have been read
and the reader says what was; they agree until absorbance begins to saturate. The answer is
the highest density still inside a stated tolerance, and it is never extrapolated past the
densest well actually run.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..readings import CorrectedOD, require

__all__ = ["LinearityFit", "fit_linear_range"]

_MIN_POINTS = 4
_ANCHOR_POINTS = 3


@dataclass(frozen=True)
class LinearityFit:
    """The linear range of one reader in one plate format."""

    od_linear_max: float
    slope: float
    tolerance: float
    saturation_observed: bool
    n_points: int

    def summary(self) -> str:
        edge = "saturation seen" if self.saturation_observed else "no saturation reached"
        return (f"linear to OD {self.od_linear_max:.3f} within {self.tolerance:.0%} "
                f"({edge}, {self.n_points} dilutions)")


def fit_linear_range(expected_od, measured_od: CorrectedOD,
                     tolerance: float = 0.05) -> LinearityFit:
    """Highest optical density the reader still reports within ``tolerance``.

    Args:
        expected_od: Density each well should read, from its dilution factor.
        measured_od: Density the reader actually reported, blank corrected.
        tolerance: Largest relative departure from the low-density line still called linear.
    """
    expected = np.asarray(expected_od, dtype=float)
    measured = require(measured_od, CorrectedOD, name="measured_od").array
    if expected.shape != measured.shape:
        raise ValueError("expected_od and measured_od must be the same length")
    if len(expected) < _MIN_POINTS:
        raise ValueError(f"need at least {_MIN_POINTS} dilutions, got {len(expected)}")

    order = np.argsort(expected)
    expected, measured = expected[order], measured[order]

    anchor = slice(0, max(_ANCHOR_POINTS, len(expected) // 3))
    slope = float(np.sum(expected[anchor] * measured[anchor])
                  / max(float(np.sum(expected[anchor] ** 2)), 1e-12))

    predicted = slope * expected
    departure = np.abs(measured - predicted) / np.maximum(np.abs(predicted), 1e-12)
    within = departure <= tolerance

    saturated = np.flatnonzero(~within)
    if saturated.size == 0:
        return LinearityFit(float(expected[-1]), slope, tolerance, False, len(expected))
    first = int(saturated[0])
    highest_linear = float(expected[first - 1]) if first > 0 else float(expected[0])
    return LinearityFit(highest_linear, slope, tolerance, True, len(expected))
