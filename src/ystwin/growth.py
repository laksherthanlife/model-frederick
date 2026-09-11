"""Specific growth rate from an optical-density trace: mu = d(ln X)/dt.

Differentiating noise amplifies it, so the log trace is smoothed and the derivative taken from the
fitted polynomial. Not a convenience statistic: mu divides the reporter by dilution, so a biased
mu propagates straight into a biased promoter-activity estimate.

**Everything here takes a** :class:`~ystwin.readings.CorrectedOD`, and that is now a type
rather than a sentence in an Args block. A raw trace still carries the medium's absorbance
-- about 0.09 on this project's plates, against cultures spanning 0.09 to 0.61 -- so
``d ln(OD)/dt`` computed on one understates the rate early and converges late. It returns a
number, not an error, and that number divides the reporter.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter

from .readings import CorrectedOD, require

__all__ = [
    "growth_rate_uncertainty",
    "growth_window",
    "max_specific_growth_rate",
    "specific_growth_rate",
]

_MIN_WINDOW_POINTS = 5


def _regular_grid(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    """Savitzky-Golay assumes uniform spacing; resample when the clock drifts."""
    dt = np.diff(t)
    if dt.size and np.allclose(dt, dt[0], rtol=1e-6, atol=1e-9):
        return t, y, False
    grid = np.linspace(t[0], t[-1], t.size)
    return grid, np.interp(grid, t, y), True


def specific_growth_rate(
    times_h: np.ndarray,
    optical_density: CorrectedOD,
    window_h: float = 1.5,
    polyorder: int = 2,
) -> np.ndarray:
    """Instantaneous specific growth rate (1/h) at each supplied timepoint.

    Args:
        times_h: Elapsed time in hours, ascending.
        optical_density: Blank-corrected OD, strictly positive. ``RawOD.minus_blank`` is
            how one is made.
        window_h: Width of the smoothing window in hours.
        polyorder: Order of the local polynomial fit.

    Raises:
        TypeError: if handed a raw density or a bare array.
    """
    t = np.asarray(times_h, dtype=float)
    od = require(optical_density, CorrectedOD, name="optical_density").array
    if t.shape != od.shape:
        raise ValueError("times and optical density must have the same shape")
    if t.size < _MIN_WINDOW_POINTS:
        raise ValueError(f"need at least {_MIN_WINDOW_POINTS} timepoints")
    if np.any(np.diff(t) <= 0):
        raise ValueError("times must be strictly increasing")
    if not np.all(od > 0):
        raise ValueError(
            "optical density must be strictly positive to take a logarithm; "
            "blank-correct first and drop or floor non-positive wells"
        )

    grid, log_od, resampled = _regular_grid(t, np.log(od))
    dt = float(grid[1] - grid[0])

    window = int(round(window_h / dt))
    if window % 2 == 0:
        window += 1
    window = max(window, polyorder + 2 + (polyorder % 2 == 0))
    if window % 2 == 0:
        window += 1
    window = min(window, grid.size if grid.size % 2 == 1 else grid.size - 1)
    if window <= polyorder:
        polyorder = max(1, window - 1)

    mu = savgol_filter(log_od, window_length=window, polyorder=polyorder, deriv=1, delta=dt)
    if resampled:
        mu = np.interp(t, grid, mu)
    return mu


_DEFAULT_DETECTION_OD = 0.02
_DEFAULT_RATE_QUANTILE = 0.9


def growth_rate_uncertainty(times_h, od: CorrectedOD) -> float:
    """Standard error of the specific growth rate estimated from one optical density trace.

    A growth rate is the slope of log optical density against time, so this is the standard
    error of that slope: the residual spread about the fit, over the spread of the time
    points and the root of their number. Every term is measurable on the trace, so the
    quantity is derived rather than fitted -- which matters because it is the error that
    propagates into every dilution-corrected activity, and fitting it alongside the
    measurement noise leaves both unidentifiable.

    It is an absolute error, not a relative one. A slowed culture is measured about as
    accurately in absolute terms as a fast one, which is exactly why the relative error
    grows as growth falls, and why the confound has dose structure.
    """
    t = np.asarray(times_h, dtype=float)
    y = require(od, CorrectedOD, name="od").array
    if t.shape != y.shape:
        raise ValueError("times_h and od must be the same length")
    if len(t) < 3:
        raise ValueError(f"need at least 3 points to fit a slope, got {len(t)}")

    usable = np.isfinite(t) & np.isfinite(y) & (y > 0)
    if usable.sum() < 3:
        return float("inf")
    t, logy = t[usable], np.log(y[usable])

    centred = t - t.mean()
    denominator = float(centred @ centred)
    if denominator <= 0:
        return float("inf")
    slope = float(centred @ (logy - logy.mean())) / denominator
    residual = logy - (logy.mean() + slope * centred)
    dof = len(t) - 2
    if dof <= 0:
        return float("inf")
    variance = float(residual @ residual) / dof
    return float(np.sqrt(variance / denominator))


def growth_window(
    times_h: np.ndarray,
    optical_density: CorrectedOD,
    min_od: float = _DEFAULT_DETECTION_OD,
) -> np.ndarray:
    """Mask of the timepoints where density is above the detection floor.

    These plates are inoculated at raw OD ~0.11 against a blank of ~0.11, so the
    blank-corrected density starts near zero. The log-derivative there is dominated
    by reader noise, and any statistic computed over it measures the reader.

    Args:
        times_h: Ascending time grid, hours.
        optical_density: Blank-corrected OD.
        min_od: Smallest density treated as a real measurement.
    """
    od = require(optical_density, CorrectedOD, name="optical_density").array
    return np.asarray(od >= min_od)


def max_specific_growth_rate(
    times_h: np.ndarray,
    optical_density: CorrectedOD,
    min_od: float = _DEFAULT_DETECTION_OD,
    window_h: float = 1.5,
    quantile: float = _DEFAULT_RATE_QUANTILE,
) -> float:
    """Peak specific growth rate, taken robustly.

    Two differences from ``max(specific_growth_rate(...))``: points below the
    detection floor are excluded, and a high quantile is taken rather than the
    maximum, so one spike in one well does not become the reported rate.

    Args:
        times_h: Ascending time grid, hours.
        optical_density: Blank-corrected OD.
        min_od: Detection floor.
        window_h: Smoothing window for the underlying rate estimate.
        quantile: Quantile of the in-window rates to report.

    Returns:
        The rate in 1/h, or NaN if the culture never clears the detection floor.
    """
    t = np.asarray(times_h, dtype=float)
    od = require(optical_density, CorrectedOD, name="optical_density").array
    mask = growth_window(t, CorrectedOD(od), min_od)
    if mask.sum() < _MIN_WINDOW_POINTS:
        return float("nan")
    mu = specific_growth_rate(t[mask], CorrectedOD(od[mask]), window_h=window_h)
    return float(np.nanquantile(mu, quantile))
