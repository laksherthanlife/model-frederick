"""Does a stress state learned from some stressors describe one it never saw?

A latent state that only re-describes its training stressors is a lookup table with extra
steps. The test that separates the two is to withhold a stressor entirely, fit on the
rest, then reveal only some of the held-out stressor's channels and predict the others.

The result is reported next to the geometry that explains it. Each stressor is a ray in
module space, so transfer can only reach a held-out stressor to the extent its direction
lies in the span of the training ones -- and the suite proves that by showing transfer
collapse for a stressor whose direction is private.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .latent import fit_latent
from .validation import stressor_fold, subset_panel

__all__ = [
    "TransferResult",
    "leave_one_stressor_out",
    "select_dimension_by_transfer",
    "subspace_alignment",
    "transfer_test",
]


@dataclass(frozen=True)
class TransferResult:
    """Prediction of unseen channels under an unseen stressor."""

    held_out: str
    trained_on: list[str]
    n_states: int
    observed: tuple[int, ...]
    predicted: tuple[int, ...]
    r2: float
    oracle_r2: float
    alignment: float
    n_train: int = 0
    n_test: int = 0
    n_excluded: int = 0
    selection: str = "fixed width supplied before outer scoring"


def subspace_alignment(dataset, held_out: str, n_states: int) -> float:
    """Fraction of a held-out stressor's variation living in the fitted latent subspace.

    This is what bounds transfer. Measuring against the full span of the training
    directions reads 1.0 for almost everything, because seven stressors already span seven
    modules; the model only ever carries ``n_states`` directions, so the question that
    predicts anything is how much of the held-out stressor those directions reach.
    """
    fold = stressor_fold(dataset, held_out)
    train, held = dataset.readings[fold.train], dataset.readings[fold.test]
    centre = train.mean(axis=0)
    _, _, vt = np.linalg.svd(train - centre, full_matrices=False)
    basis = vt[:n_states]
    centred = held - centre
    total = float(np.sum(centred**2))
    if total <= 0:
        return 0.0
    captured = float(np.sum((centred @ basis.T) ** 2))
    return float(np.clip(captured / total, 0.0, 1.0))


def _fit_model(train, n_states, seed):
    """Centre, factorise, and keep what prediction needs: loadings, prior scale, noise."""
    centre = train.mean(axis=0)
    fit = fit_latent(train - centre, n_states=n_states, seed=seed)
    return centre, fit.loadings, np.var(fit.states, axis=0), max(fit.heldout_error, 1e-12)


def _predict_from(model, held, observed, predicted):
    """Infer the held-out state by MAP under the prior the training states define.

    Least squares would be free to pick an enormous state that matches the revealed
    channels while the predicted ones diverge, which is what an unseen stressor provokes.
    Both the prior scale and the noise level come from the training fit, so shrinkage
    costs no new constant and an unreachable stressor lands on the training mean instead
    of somewhere absurd. The noise is the cross-validated error rather than the in-sample
    residual, which collapses to zero once the states saturate the channels and would hand
    the problem straight back to least squares.
    """
    observed, predicted = np.asarray(observed), np.asarray(predicted)
    centre, loadings, prior, noise = model
    design = loadings[:, observed]
    precision = (design @ design.T) / noise + np.diag(1.0 / np.maximum(prior, 1e-12))
    target = design @ (held - centre)[:, observed].T / noise
    states = np.linalg.solve(precision, target)
    return (states.T @ loadings[:, predicted]) + centre[predicted]


def _r2(truth: np.ndarray, prediction: np.ndarray, baseline: np.ndarray) -> float:
    residual = float(np.sum((truth - prediction) ** 2))
    total = float(np.sum((truth - baseline) ** 2))
    return 1.0 - residual / total if total > 0 else float("nan")


def transfer_test(
    dataset,
    held_out: str,
    n_states: int | None,
    observed,
    seed: int = 0,
    oracle_model=None,
    max_states: int = 4,
    include_oracle: bool = True,
) -> TransferResult:
    """Fit without one stressor, then predict its unrevealed channels.

    Args:
        dataset: Panel readings with a stressor label per sample.
        held_out: Stressor withheld from fitting.
        n_states: Latent dimension of the fitted state.
        observed: Channel indices revealed for the held-out samples.
        seed: Seed passed to the latent fit.
        oracle_model: A fit over the whole dataset, which is the same whichever stressor is
            held out. A sweep computes it once and passes it in rather than refitting it
            per stressor; the result is identical either way.
    """
    observed = tuple(sorted(set(observed)))
    n_channels = dataset.readings.shape[1]
    if not observed or any(not isinstance(i, (int, np.integer)) or i < 0 or i >= n_channels
                           for i in observed):
        raise ValueError("observed must contain valid non-negative channel indices")
    predicted = tuple(i for i in range(n_channels) if i not in observed)
    if not predicted:
        raise ValueError("every channel was revealed, so there is nothing left to predict")
    fold = stressor_fold(dataset, held_out)
    selection = "fixed width supplied before outer scoring"
    if n_states is None:
        n_states = _select_nested_width(subset_panel(dataset, fold.train), max_states,
                                        observed, seed)
        selection = "nested leave-stressor-out on outer training only"
    if not isinstance(n_states, (int, np.integer)) or n_states < 1:
        raise ValueError("n_states must be a positive integer or None for nested selection")
    if n_states > len(observed):
        raise ValueError(
            f"cannot determine {n_states} states from {len(observed)} observed channels")

    train, held = dataset.readings[fold.train], dataset.readings[fold.test]
    if not np.isfinite(train).all() or not np.isfinite(held).all():
        raise ValueError("channel transfer requires finite readings")
    trained_on = sorted(set(dataset.labels[fold.train]))

    index = np.asarray(predicted)
    truth = held[:, index]
    baseline = np.repeat(train[:, index].mean(axis=0)[None, :], len(held), axis=0)
    prediction = _predict_from(_fit_model(train, n_states, seed), held, observed, predicted)
    oracle_r2 = float("nan")
    if include_oracle:
        oracle_model = oracle_model or _fit_model(dataset.readings, n_states, seed)
        oracle = _predict_from(oracle_model, held, observed, predicted)
        oracle_r2 = _r2(truth, oracle, baseline)

    return TransferResult(
        held_out=held_out,
        trained_on=trained_on,
        n_states=n_states,
        observed=observed,
        predicted=predicted,
        r2=_r2(truth, prediction, baseline),
        oracle_r2=oracle_r2,
        alignment=subspace_alignment(dataset, held_out, n_states),
        n_train=len(fold.train), n_test=len(fold.test), n_excluded=len(fold.excluded),
        selection=selection,
    )


def leave_one_stressor_out(dataset, n_states: int | None, observed, seed: int = 0,
                          only=None, max_states: int = 4, include_oracle: bool = True) -> pd.DataFrame:
    """Transfer score for every stressor in turn, beside the alignment that explains it.

    Args:
        dataset: Panel readings with a stressor label per sample.
        n_states: Latent dimension of the fitted state.
        observed: Channel indices revealed for the held-out samples.
        seed: Seed passed to the latent fit.
        only: Stressors to hold out; defaults to every label present. A sweep costs two
            latent fits per stressor, so a broad landscape wants this narrowed.
    """
    wanted = set(only) if only is not None else None
    oracle_models = {}
    rows = []
    for name in sorted(set(dataset.labels)):
        if wanted is not None and name not in wanted:
            continue
        chosen = n_states
        if chosen is None:
            fold = stressor_fold(dataset, name)
            chosen = _select_nested_width(subset_panel(dataset, fold.train),
                                          max_states, observed, seed)
        if include_oracle and chosen not in oracle_models:
            oracle_models[chosen] = _fit_model(dataset.readings, chosen, seed)
        result = transfer_test(dataset, name, chosen, observed, seed, oracle_models.get(chosen),
                               include_oracle=include_oracle)
        rows.append({
            "held_out": name,
            "alignment": result.alignment,
            "r2": result.r2,
            "oracle_r2": result.oracle_r2,
            "n_states": chosen,
            "n_train": result.n_train,
            "n_test": result.n_test,
            "n_excluded": result.n_excluded,
            "selection": ("nested leave-stressor-out on outer training only"
                          if n_states is None else result.selection),
            "baseline": "outer training channel mean",
            "oracle_role": ("in-sample diagnostic, not a held-out prediction"
                            if include_oracle else "not computed during selection"),
        })
    return pd.DataFrame(rows)


def _select_nested_width(training, max_states, observed, seed):
    labels = [str(label) for label in sorted(set(training.labels)) if "+" not in str(label)]
    if len(labels) < 2:
        raise ValueError("nested dimension selection needs at least two training stressors")
    return select_dimension_by_transfer(training, max_states, observed, seed, only=labels)


def select_dimension_by_transfer(dataset, max_states: int, observed, seed: int = 0, only=None) -> int:
    """Latent dimension that best predicts a stressor the model never saw.

    Entry-wise cross-validation picks the rank that best fills a gap in a well whose other
    channels are known, which on a smoothly decaying spectrum is one. Predicting an unseen
    stressor is a different question and generally wants more states, so a model built for
    transfer has to be selected on transfer.
    """
    if not isinstance(max_states, (int, np.integer)) or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    if len(set(dataset.labels)) < 2:
        raise ValueError("dimension selection needs at least two training treatments")
    best, best_score = None, -np.inf
    for n_states in range(1, min(max_states, len(observed)) + 1):
        table = leave_one_stressor_out(dataset, n_states, observed, seed, only=only, include_oracle=False)
        score = float(table["r2"].median())
        if np.isfinite(score) and score > best_score:
            best, best_score = n_states, score
    if best is None:
        raise ValueError("no dimension has a finite inner transfer score")
    return best
