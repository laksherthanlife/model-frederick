"""Joint deconvolution of N co-expressed reporters, treating maturation as information.

Steady state with a maturation step is

    R_i = k_i * m_i / [(mu + m_i + kdeg_i)(mu + kdeg_i)]

A *matched* pair cancels mu, so its ratio reads relative promoter activity. A *mismatched*
pair does not: its ratio is monotonic in mu, so it measures growth. Both are useful, and a
panel wants both -- a slow fluorophore is the best growth sensor available, not a liability.

Maturation enters as a known parameter, so it is calibrated out rather than avoided.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "growth_from_pair",
    "growth_from_reference",
    "growth_sensitivity",
    "solve_activities",
    "steady_state_gain",
]

_MIN_MISMATCH = 1e-6
_MU_BOUNDS = (1e-4, 5.0)


def steady_state_gain(growth_rate: float, spec) -> float:
    """Mature reporter per unit promoter activity at a given growth rate."""
    mu = max(float(growth_rate), 1e-12)
    maturation = 1.0 / max(spec.maturation_h, 1e-12)
    return maturation / ((mu + maturation + spec.k_deg) * (mu + spec.k_deg))


def _gain_ratio(growth_rate: float, slow, fast) -> float:
    return steady_state_gain(growth_rate, slow) / steady_state_gain(growth_rate, fast)


def growth_sensitivity(growth_rate: float, slow, fast, step: float = 1e-4) -> float:
    """|d log(ratio) / d log(mu)| -- how much a pair's ratio moves with growth."""
    mu = max(float(growth_rate), 1e-6)
    low, high = _gain_ratio(mu * (1 - step), slow, fast), _gain_ratio(mu * (1 + step), slow, fast)
    if low <= 0 or high <= 0:
        return 0.0
    return float(abs((np.log(high) - np.log(low)) / (2 * step)))


def growth_from_pair(observed_ratio: float, slow, fast, relative_activity: float = 1.0) -> float:
    """Recover growth rate from the ratio of two constitutive reporters.

    Args:
        observed_ratio: Measured ``R_slow / R_fast``.
        slow, fast: Their specs; maturation must differ or mu does not enter.
        relative_activity: Ratio of their promoter activities, calibrated once.

    Raises:
        ValueError: if maturation matches, in which case the ratio carries no mu.
    """
    if abs(slow.maturation_h - fast.maturation_h) < _MIN_MISMATCH:
        raise ValueError(
            "the two reporters share a maturation rate, so mu cancels from their ratio "
            "and cannot be recovered from it; pair a slow fluorophore with a fast one"
        )

    def residual(mu):
        return relative_activity * _gain_ratio(mu, slow, fast) - observed_ratio

    low, high = _MU_BOUNDS
    if residual(low) * residual(high) > 0:
        raise ValueError(
            f"ratio {observed_ratio:.4g} is outside the range this pair can produce "
            f"over mu in [{low}, {high}]"
        )
    return float(brentq(residual, low, high, xtol=1e-12, rtol=1e-12))


def growth_from_reference(observed: float, activity: float, spec) -> float:
    """Recover growth from one reference channel whose absolute activity is known.

    Needs absolute calibration of that channel, and works whatever the maturation is,
    because the gain itself is monotonic in mu. The pair route below needs only a
    relative calibration but does require a maturation mismatch.
    """
    target = float(observed) / float(activity)
    low, high = _MU_BOUNDS
    if (steady_state_gain(low, spec) - target) * (steady_state_gain(high, spec) - target) > 0:
        raise ValueError(
            f"gain {target:.4g} is outside the range this reporter can produce over "
            f"mu in [{low}, {high}]"
        )
    return float(brentq(lambda mu: steady_state_gain(mu, spec) - target,
                        low, high, xtol=1e-14, rtol=1e-13))


def solve_activities(
    observed: np.ndarray,
    specs,
    growth_rate: float,
) -> np.ndarray:
    """Promoter activity per channel, given the shared growth rate.

    Maturation enters as a known per-channel parameter and is divided out, so a slow
    fluorophore needs no special handling once it has been calibrated.
    """
    observed = np.asarray(observed, dtype=float)
    if observed.size != len(specs):
        raise ValueError("one reading is required per channel")
    if observed.size < 1:
        raise ValueError("need at least one channel")
    return np.array([o / steady_state_gain(growth_rate, s) for o, s in zip(observed, specs)])
