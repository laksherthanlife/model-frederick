"""Reporter loss-rate calibration.

k_deg is the difference between a reporter that reports and one that reports 1/mu. Three routes:
a translation-shutoff chase (strongest), splitting total loss across two read intervals to
separate degradation from photobleaching, and a free lower bound from ordinary growth data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..growth import specific_growth_rate
from ..reporter import _smooth_derivative, default_activity_window_h
from ..readings import RawOD, RawRFU, require

__all__ = [
    "ChaseFit",
    "DecayFit",
    "LossPartition",
    "fit_decay_rate",
    "fit_kdeg_from_chase",
    "minimum_consistent_kdeg",
    "partition_loss",
]


@dataclass(frozen=True)
class DecayFit:
    """Exponential decay rate from a log-linear regression."""

    rate: float
    stderr: float
    r_squared: float
    n_points: int

    @property
    def half_life_h(self) -> float:
        return float(np.log(2) / self.rate) if self.rate > 0 else float("inf")


@dataclass(frozen=True)
class ChaseFit:
    """Reporter degradation rate from a translation-shutoff chase."""

    k_deg: float
    stderr: float
    r_squared: float
    total_loss_rate: float
    mean_growth_rate: float
    bleaching_rate: float

    @property
    def half_life_h(self) -> float:
        return float(np.log(2) / self.k_deg) if self.k_deg > 0 else float("inf")

    def summary(self) -> str:
        return (
            f"k_deg = {self.k_deg:.4f} +/- {self.stderr:.4f} /h "
            f"(half-life {self.half_life_h:.1f} h); "
            f"total loss {self.total_loss_rate:.4f}, growth {self.mean_growth_rate:.4f}, "
            f"bleaching {self.bleaching_rate:.4f}; R^2 {self.r_squared:.3f}"
        )


@dataclass(frozen=True)
class LossPartition:
    """Total loss split into time- and read-proportional parts."""

    k_deg: float
    bleach_per_read: float
    n_conditions: int

    def total_at_interval(self, read_interval_h: float) -> float:
        """Predicted loss rate when reading every ``read_interval_h`` hours."""
        return self.k_deg + self.bleach_per_read / read_interval_h


def fit_decay_rate(times_h: np.ndarray, signal: np.ndarray) -> DecayFit:
    """Fit ``signal ~ exp(-rate * t)`` by regression on the logarithm.

    A positive rate means decay. A negative rate is reported as such: on a chase it
    means translation was not actually blocked, which the experimenter needs to know.

    Args:
        times_h: Ascending time grid, hours.
        signal: Strictly positive signal on the same grid.
    """
    t = np.asarray(times_h, dtype=float)
    y = np.asarray(signal, dtype=float)
    if t.shape != y.shape:
        raise ValueError("times and signal must have the same shape")
    if t.size < 3:
        raise ValueError("need at least 3 timepoints to fit a decay rate")
    if not np.all(y > 0):
        raise ValueError("signal must be strictly positive to take a logarithm")

    log_y = np.log(y)
    slope, intercept = np.polyfit(t, log_y, 1)
    predicted = slope * t + intercept
    resid = log_y - predicted
    ss_tot = float(np.sum((log_y - log_y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else 1.0

    dof = max(t.size - 2, 1)
    s_xx = float(np.sum((t - t.mean()) ** 2))
    stderr = float(np.sqrt(np.sum(resid**2) / dof / s_xx)) if s_xx > 0 else float("inf")
    return DecayFit(rate=float(-slope), stderr=stderr, r_squared=float(r2), n_points=int(t.size))


def fit_kdeg_from_chase(
    times_h: np.ndarray,
    rfu: RawRFU,
    optical_density: RawOD,
    od_blank: float,
    rfu_background: float = 0.0,
    bleaching_rate: float = 0.0,
) -> ChaseFit:
    """Reporter degradation rate from a translation-shutoff chase.

    After shutoff ``k_synth = 0``, so the per-cell signal obeys
    ``d(ln R)/dt = -(mu + k_deg + bleaching)``. Residual growth and any measured
    photobleaching are subtracted to leave degradation.

    Args:
        times_h: Ascending time grid from the moment of shutoff, hours.
        rfu: Raw reporter signal.
        optical_density: Raw OD on the same grid, used for the residual growth rate.
        od_blank: Media-only OD.
        rfu_background: Media plus autofluorescence baseline.
        bleaching_rate: Loss per hour attributable to illumination, from a
            no-chase control read on the same schedule.
    """
    t = np.asarray(times_h, dtype=float)
    density = require(optical_density, RawOD, name="optical_density").minus_blank(od_blank)
    signal = require(rfu, RawRFU, name="rfu").minus_background(rfu_background)
    od = density.array
    if not np.all(od > 0):
        raise ValueError("blank-corrected OD must be positive")

    per_cell = signal.per(density).array
    decay = fit_decay_rate(t, per_cell)
    mu = float(np.mean(specific_growth_rate(t, density))) if od.max() / od.min() > 1.02 else 0.0
    k_deg = decay.rate - mu - bleaching_rate
    return ChaseFit(
        k_deg=float(k_deg),
        stderr=decay.stderr,
        r_squared=decay.r_squared,
        total_loss_rate=decay.rate,
        mean_growth_rate=mu,
        bleaching_rate=float(bleaching_rate),
    )


def partition_loss(observations: list[tuple[float, float]]) -> LossPartition:
    """Split total loss into degradation and photobleaching.

    Degradation runs with the clock; photobleaching runs with the number of reads.
    At read interval ``dt`` the observed rate is ``k_deg + bleach_per_read / dt``,
    so the two components separate only across different read intervals.

    Args:
        observations: ``(read_interval_h, observed_loss_rate)`` pairs, one per plate.

    Raises:
        ValueError: unless at least two distinct read intervals are present. One
            interval leaves the split unidentifiable, however many plates were run.
    """
    if len(observations) < 2:
        raise ValueError(
            "need at least two distinct read intervals to separate photobleaching "
            "from degradation; one interval gives only the total loss rate"
        )
    intervals = np.array([o[0] for o in observations], dtype=float)
    rates = np.array([o[1] for o in observations], dtype=float)
    if not np.all(intervals > 0):
        raise ValueError("read intervals must be positive")
    inverse = 1.0 / intervals
    if np.ptp(inverse) < 1e-12:
        raise ValueError(
            "need at least two distinct read intervals to separate photobleaching "
            "from degradation; all supplied plates share one interval"
        )
    bleach, k_deg = np.polyfit(inverse, rates, 1)
    return LossPartition(
        k_deg=float(k_deg), bleach_per_read=float(bleach), n_conditions=len(observations)
    )


def minimum_consistent_kdeg(
    times_h: np.ndarray,
    optical_density: RawOD,
    rfu: RawRFU,
    od_blank: float = 0.0,
    rfu_background: float = 0.0,
    window_h: float | None = None,
) -> float:
    """Smallest loss rate under which inferred promoter activity is never negative.

    ``k_synth = dR/dt + (mu + k_deg) R >= 0`` implies
    ``k_deg >= max_t[ -d(ln R)/dt - mu ]``. A lower bound available from ordinary
    growth data with no extra experiment, and the right prior until a chase is run.
    """
    t = np.asarray(times_h, dtype=float)
    density = require(optical_density, RawOD, name="optical_density").minus_blank(od_blank)
    signal = require(rfu, RawRFU, name="rfu").minus_background(rfu_background)
    od = density.array
    if not np.all(od > 0):
        raise ValueError("blank-corrected OD must be positive")

    if window_h is None:
        window_h = default_activity_window_h(t)
    mu = specific_growth_rate(t, density)
    per_cell = signal.per(density).array
    positive = per_cell > 0
    if positive.sum() < 5:
        raise ValueError("too few positive reporter readings to bound the loss rate")
    dln = _smooth_derivative(t[positive], np.log(per_cell[positive]), window_h, polyorder=3)
    implied = -dln - mu[positive]
    return float(max(0.0, np.nanmax(implied)))
