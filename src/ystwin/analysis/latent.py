"""Multi-channel latent stress state, with its dimension chosen by evidence.

A latent state must not be a reporter rename, so the dimension is selected by BIC over
the channels rather than set to the number of sensors. Inputs are dilution-corrected
promoter activities: the growth confound is removed before anything is called latent.

Missing readings are masked, never zero-filled -- a zero asserts a dark cell rather than
an absent measurement, and it drags the recovered state toward it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["LatentFit", "fit_latent", "select_dimension"]

_MAX_ITERATIONS = 200
_TOLERANCE = 1e-9


@dataclass(frozen=True)
class LatentFit:
    """Recovered states, how each channel loads on them, and the evidence."""

    states: np.ndarray = field(repr=False)
    loadings: np.ndarray = field(repr=False)
    explained_per_channel: np.ndarray = field(repr=False)
    bic: float
    heldout_error: float
    n_states: int
    n_observations: int
    n_channels: int

    def summary(self) -> str:
        worst = float(np.min(self.explained_per_channel))
        return (
            f"{self.n_states} state(s) over {self.n_channels} channels, "
            f"{self.n_observations} observations; held-out error "
            f"{self.heldout_error:.4f}; least-explained channel {worst:.0%}"
        )


def _centre(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centre each channel on its observed mean; do not rescale.

    Scaling to unit variance would inflate a low-variance channel until it looked like a
    real axis, and a pure-noise sensor would earn its own latent state. Inputs are
    dilution-corrected activities, already in comparable units, so centring is enough.
    Channels in genuinely different units must be rescaled by the caller.
    """
    means = np.nanmean(data, axis=0)
    return data - means, means


def _factorise(scaled, observed, n_states):
    """Iterated SVD, re-imputing missing entries from the current fit rather than zero."""
    filled = np.where(observed, scaled, 0.0)
    reconstruction = np.zeros_like(filled)
    previous = None
    for _ in range(_MAX_ITERATIONS):
        u, sv, vt = np.linalg.svd(filled, full_matrices=False)
        states = u[:, :n_states] * sv[:n_states]
        loadings = vt[:n_states]
        reconstruction = states @ loadings
        filled = np.where(observed, scaled, reconstruction)
        if previous is not None and np.max(np.abs(reconstruction - previous)) < _TOLERANCE:
            break
        previous = reconstruction
    return states, loadings, reconstruction


def _heldout_error(scaled, observed, n_states, folds=4, fraction=0.2, seed=0):
    """Mean error on entries hidden from the fit; immune to the perfect-fit degeneracy."""
    rng = np.random.default_rng(seed)
    rows, cols = np.where(observed)
    if len(rows) < 10:
        return float("inf")
    errors = []
    for fold in range(folds):
        pick = rng.choice(len(rows), size=max(1, int(len(rows) * fraction)), replace=False)
        train = observed.copy()
        train[rows[pick], cols[pick]] = False
        if not train.any(axis=0).all():
            continue
        _, means = _centre(np.where(train, scaled, np.nan))
        centred = scaled - means
        _, _, reconstruction = _factorise(centred, train, n_states)
        held = centred[rows[pick], cols[pick]] - reconstruction[rows[pick], cols[pick]]
        errors.append(float(np.mean(held**2)))
    return float(np.mean(errors)) if errors else float("inf")


def fit_latent(activities: np.ndarray, n_states: int, seed: int = 0) -> LatentFit:
    """Fit ``n_states`` shared factors across channels, masking missing readings.

    Args:
        activities: ``(observations, channels)`` of dilution-corrected activity. NaN
            marks a channel not measured for that observation.
        n_states: Number of latent states.
    """
    data = np.asarray(activities, dtype=float)
    if data.ndim != 2:
        raise ValueError("activities must be a 2-D observations-by-channels array")
    n_observations, n_channels = data.shape
    if not isinstance(n_states, (int, np.integer)) or isinstance(n_states, bool) or n_states < 1:
        raise ValueError("n_states must be a positive integer")
    if n_states > n_channels:
        raise ValueError(f"n_states ({n_states}) cannot exceed channels ({n_channels})")
    if n_observations < 2 or n_states > n_observations:
        raise ValueError("need at least two observations and no more states than observations")
    if np.isinf(data).any():
        raise ValueError("activities must be finite or NaN for missing readings")
    observed = np.isfinite(data)
    empty = [c for c in range(n_channels) if not observed[:, c].any()]
    if empty:
        raise ValueError(f"channels {empty} have no observations at all")

    scaled, _ = _centre(data)
    states, loadings, reconstruction = _factorise(scaled, observed, n_states)

    residual = np.where(observed, scaled - reconstruction, 0.0)
    n_points = int(observed.sum())
    sse = float((residual**2).sum())
    # Floor the variance: an exact fit would otherwise send the log to -inf and win.
    variance = max(sse / max(n_points, 1), 1e-12)
    n_parameters = n_states * (n_observations + n_channels)
    bic = n_points * np.log(variance) + n_parameters * np.log(max(n_points, 2))
    heldout = _heldout_error(scaled, observed, n_states, seed=seed)

    explained = []
    for c in range(n_channels):
        mask = observed[:, c]
        total = float((scaled[mask, c] ** 2).sum())
        left = float((residual[mask, c] ** 2).sum())
        explained.append(0.0 if total <= 0 else max(0.0, 1.0 - left / total))

    return LatentFit(
        states=states, loadings=loadings,
        explained_per_channel=np.asarray(explained, dtype=float),
        bic=float(bic), heldout_error=heldout, n_states=n_states,
        n_observations=n_observations, n_channels=n_channels,
    )


def select_dimension(activities: np.ndarray, max_states: int = 4, seed: int = 0) -> int:
    """States minimising held-out reconstruction error, capped by the channel count.

    Cross-validation rather than BIC: with as many factors as channels the fit is exact,
    the residual variance goes to zero and any in-sample criterion picks the largest model.
    """
    data = np.asarray(activities, dtype=float)
    if data.ndim != 2 or min(data.shape) < 1:
        raise ValueError("activities must be a nonempty 2-D observations-by-channels array")
    if not isinstance(max_states, (int, np.integer)) or isinstance(max_states, bool) or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    ceiling = min(max_states, *data.shape)
    scores = {n: fit_latent(data, n, seed=seed).heldout_error for n in range(1, ceiling + 1)}
    if not any(np.isfinite(score) for score in scores.values()):
        raise ValueError("too few observed entries to select a dimension by held-out error")
    return min(scores, key=scores.get)
