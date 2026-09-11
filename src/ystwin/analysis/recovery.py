"""Whether a fitted latent state carries the module activities, or only compresses them.

Predicting held-out channels well proves the states span the data; it does not prove they
mean anything. The states are identified only up to a rotation, so a model can impute
perfectly while no single state corresponds to any regulon.

In simulation the true activities are known, so the question is answerable directly:
regress each module on the fitted states and report the variance recovered. Rotation does
not matter to that, because a linear readout absorbs it -- what is being asked is whether
the module lives in the fitted subspace at all.
"""

from __future__ import annotations

import numpy as np

from ..generator.stress_panel import MODULES
from .latent import fit_latent

__all__ = ["module_recovery"]


def module_recovery(dataset, n_states: int, seed: int = 0, shuffle: bool = False) -> dict[str, float]:
    """Fraction of each true module activity recoverable by a linear readout of the states.

    Args:
        dataset: Panel readings carrying the true module activities behind them.
        n_states: Latent dimension to fit.
        seed: Seeds the shuffle **only**. It does not perturb the fit, and it cannot:
            :func:`~ystwin.analysis.latent.fit_latent` recovers its states by SVD,
            which is deterministic, so the seed reaches only the cross-validated
            held-out error -- a quantity this function does not use. With
            ``shuffle=False`` the seed therefore has no effect whatsoever.

            Worth stating because the failure mode is silent and inviting: computing a
            spread of unshuffled recoveries across seeds returns exactly zero, which
            reads as a perfectly precise estimate and is in fact no estimate at all.
            Resample the data, or vary the dataset seed, to get a spread.
        shuffle: Permute the states across samples, breaking their link to the truth. A
            metric that still scores under this is measuring the regression, not the model.
    """
    readings = dataset.readings
    fit = fit_latent(readings - readings.mean(axis=0), n_states=n_states, seed=seed)
    states = fit.states
    if shuffle:
        states = states[np.random.default_rng(seed).permutation(len(states))]
    design = np.column_stack([states, np.ones(len(states))])

    scores = {}
    for index, name in enumerate(MODULES):
        truth = dataset.modules[:, index]
        centred = truth - truth.mean()
        total = float(centred @ centred)
        if total <= 1e-12:
            scores[name] = 0.0
            continue
        coefficients, *_ = np.linalg.lstsq(design, truth, rcond=None)
        residual = truth - design @ coefficients
        scores[name] = float(np.clip(1.0 - float(residual @ residual) / total, 0.0, 1.0))
    return scores
