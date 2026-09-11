"""Compare cultures at the same optical density rather than the same clock time.

Conditions that grow at different rates are in different physiological states at any
fixed timepoint, so a plate read at t = 12 h compares exponential cells in one well with
stationary cells in another. The difference that produces looks exactly like a dose
response and is not one.

Matching on density instead asks the only fair question: at the same amount of biomass,
does this culture carry more reporter than that one? What survives is response; what
disappears was growth phase.
"""

from __future__ import annotations

import numpy as np

from ..readings import CorrectedOD, CorrectedRFU, require

__all__ = ["fold_range", "specific_fluorescence_at_density", "time_to_density"]


def _monotone(optical_density: CorrectedOD) -> np.ndarray:
    """Densities dip with settling and refocus; the running maximum is what was reached."""
    return np.maximum.accumulate(
        require(optical_density, CorrectedOD, name="optical_density").array)


def time_to_density(times_h, optical_density: CorrectedOD, target: float) -> float:
    """Hour at which the culture first reaches ``target``, or nan if it never does."""
    od = _monotone(optical_density)
    t = np.asarray(times_h, dtype=float)
    i = int(np.searchsorted(od, target))
    return float(t[i]) if i < len(t) else float("nan")


def specific_fluorescence_at_density(times_h, optical_density: CorrectedOD,
                                     rfu: CorrectedRFU, target: float) -> float:
    """RFU per OD at the moment the culture first reaches ``target`` density."""
    od = _monotone(optical_density)
    signal = require(rfu, CorrectedRFU, name="rfu").array
    i = int(np.searchsorted(od, target))
    if i >= len(od):
        return float("nan")
    return float(signal[i] / od[i])


def fold_range(values) -> float:
    """Largest over smallest, ignoring conditions that never reached the density."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size < 2 or v.min() <= 0:
        return float("nan")
    return float(v.max() / v.min())
