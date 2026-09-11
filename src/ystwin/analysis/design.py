"""Sensor-set design on Fisher information, with an explicit noise model.

Three earlier criteria were ad hoc: a spanning score that deduplicated and normalised
rows, a projection-weight threshold of 0.35, and forcing a reporter in by hand. All three
were symptoms of one gap -- with no noise model, "identifiable" had no statistical meaning.

With observations y = L x + e, e ~ N(0, S), everything follows from one object. The
information is L' S^-1 L; a module's standard error is the root of the matching diagonal
of its pseudo-inverse; identifiability is that error being small against an effect worth
detecting; and a target state is served by weighting it, not by overriding the answer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..generator.panel_experiment import OBSERVED_ACTIVITY_CV
from ..generator.stress_panel import MODULES

__all__ = [
    "NoiseModel",
    "design_score",
    "design_value",
    "fisher_information",
    "identifiable_by_precision",
    "module_standard_errors",
]

_EIGENVALUE_FLOOR = 1e-9


@dataclass(frozen=True)
class NoiseModel:
    """Measurement noise on one reporter channel.

    The default is the noise this project's own reader produces, not a round number. It was
    0.02 -- seven times tighter than measured, and a default is what a caller who does not
    know the reader silently gets, so the optimism propagated into every precision claim
    made without an explicit argument. Standard errors scale linearly in it, so the number
    does not change which reporter set wins; it changes what that set can honestly claim to
    resolve, and by a factor of seven.

    Args:
        relative_cv: Multiplicative noise per reading, the dominant term on a plate.
            Defaults to `panel_experiment.OBSERVED_ACTIVITY_CV`, measured across 28 matched
            conditions on two real plates.
        floor: Additive noise, which matters for a dim channel near background.
    """

    relative_cv: float = OBSERVED_ACTIVITY_CV
    floor: float = 0.0

    def variance(self, expected_signal: float = 1.0) -> float:
        return (self.relative_cv * expected_signal) ** 2 + self.floor**2


def fisher_information(loadings: np.ndarray, noise: NoiseModel, extra_cv=None) -> np.ndarray:
    """``L' S^-1 L``, with each channel's noise scaled to the signal it carries.

    Plate noise is multiplicative, so a reporter that reads more loudly carries
    proportionally more noise and is no better for it. Holding the noise constant instead
    would reward a promoter for picking up crosstalk, which is precisely what makes its
    module hard to separate. With relative noise the two cancel and only orthogonality
    survives -- an additive floor then puts genuine brightness back in, where it belongs.

    ``extra_cv`` adds relative noise to individual channels. Its use here is the growth
    correction: a promoter fusion is diluted by growth, so reading activity off it means
    dividing out a growth rate estimated from noisy optical density, and that variance
    lands on that channel alone. A ratiometric sensor carries none of it.
    """
    matrix = np.asarray(loadings, dtype=float)
    n_modules = len(MODULES)
    if matrix.size == 0 or matrix.shape[0] == 0:
        return np.zeros((n_modules, n_modules))
    signal = np.linalg.norm(matrix, axis=1)
    variance = np.array([noise.variance(s) for s in signal])
    if extra_cv is not None:
        extra = np.asarray(extra_cv, dtype=float)
        if extra.shape != signal.shape:
            raise ValueError(f"extra_cv needs one entry per channel, got {extra.shape}")
        variance = variance + (extra * signal) ** 2
    precision = np.where(signal > 0, 1.0 / np.maximum(variance, 1e-300), 0.0)
    return matrix.T @ np.diag(precision) @ matrix


def module_standard_errors(loadings: np.ndarray, noise: NoiseModel, extra_cv=None) -> dict[str, float]:
    """Standard error on each module's activity, infinite where it is unobserved."""
    info = fisher_information(loadings, noise, extra_cv)
    order = list(MODULES)
    covariance = np.linalg.pinv(info, rcond=1e-10)
    rank_deficient = np.abs(np.diag(info)) < 1e-12
    errors = {}
    for i, name in enumerate(order):
        if rank_deficient[i]:
            errors[name] = float("inf")
            continue
        variance = float(covariance[i, i])
        errors[name] = float(np.sqrt(variance)) if variance > 0 else float("inf")
    # A module inside the null space is unobservable however large its diagonal looks.
    null_dim = len(order) - int(np.linalg.matrix_rank(info, tol=1e-10))
    if null_dim:
        _, _, vt = np.linalg.svd(info)
        null_space = vt[len(order) - null_dim:]
        leakage = np.sum(null_space**2, axis=0)
        for i, name in enumerate(order):
            if leakage[i] > 0.5:
                errors[name] = float("inf")
    return errors


def identifiable_by_precision(
    loadings: np.ndarray, noise: NoiseModel, min_effect: float, confidence: float = 2.0,
    extra_cv=None,
) -> list[str]:
    """Modules whose standard error is small enough to resolve ``min_effect``.

    Replaces a projection-weight threshold with a statement about what can actually be
    detected: an effect of ``min_effect`` must exceed ``confidence`` standard errors.
    """
    errors = module_standard_errors(loadings, noise, extra_cv)
    return [name for name, se in errors.items() if np.isfinite(se) and confidence * se <= min_effect]


def _weighted_information(
    loadings: np.ndarray, noise: NoiseModel, priority: dict[str, float] | None, extra_cv=None
) -> np.ndarray:
    order = list(MODULES)
    if priority:
        unknown = set(priority) - set(order)
        if unknown:
            raise KeyError(f"no such module(s) to prioritise: {sorted(unknown)}")
    info = fisher_information(loadings, noise, extra_cv)
    if priority:
        scaling = np.diag(np.sqrt([priority.get(name, 1.0) for name in order]))
        info = scaling @ info @ scaling
    return info


def design_score(
    loadings: np.ndarray, noise: NoiseModel, priority: dict[str, float] | None = None,
    extra_cv=None,
) -> float:
    """Information a reporter set carries: log pseudo-determinant of ``L' S^-1 L``.

    Fewer reporters than modules leaves the matrix singular, so a plain log-determinant
    needs a ridge -- and the size of that ridge silently prices how much one more estimable
    module is worth against better conditioning. Summing over the non-null eigenvalues
    instead prices nothing, and comparing sets of unequal rank is left to `design_value`,
    which says outright that estimable modules come first.
    """
    info = _weighted_information(loadings, noise, priority, extra_cv)
    eigenvalues = np.linalg.eigvalsh(info)
    live = eigenvalues[eigenvalues > _EIGENVALUE_FLOOR * max(eigenvalues.max(initial=0.0), 1.0)]
    if live.size == 0:
        return float("-inf")
    return float(np.sum(np.log(live)))


def design_value(
    loadings: np.ndarray,
    noise: NoiseModel,
    min_effect: float,
    priority: dict[str, float] | None = None,
    confidence: float = 2.0,
    extra_cv=None,
) -> tuple[float, float]:
    """What a reporter set is worth: modules it can actually measure, then information.

    Ranking on information alone rewards a set that spans widely but resolves nothing, so
    the primary term counts modules whose standard error clears ``min_effect``, weighted by
    what the build is for. Information breaks ties among sets that measure the same things.
    """
    if priority:
        unknown = set(priority) - set(MODULES)
        if unknown:
            raise KeyError(f"no such module(s) to prioritise: {sorted(unknown)}")
    resolved = identifiable_by_precision(loadings, noise, min_effect, confidence, extra_cv)
    weight = sum((priority or {}).get(name, 1.0) for name in resolved)
    return float(weight), design_score(loadings, noise, extra_cv=extra_cv)
