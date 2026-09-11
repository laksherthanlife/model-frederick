from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ValidationFold:
    held_out: str
    train: np.ndarray
    test: np.ndarray
    excluded: np.ndarray

    def __post_init__(self):
        for name in ("train", "test", "excluded"):
            indices = np.asarray(getattr(self, name))
            if (indices.ndim != 1 or indices.dtype.kind not in "iu" or np.any(indices < 0)
                    or len(np.unique(indices)) != len(indices)):
                raise ValueError("fold indices must be unique non-negative integer arrays")
            indices = indices.copy()
            indices.setflags(write=False)
            object.__setattr__(self, name, indices)
        if not len(self.train) or not len(self.test):
            raise ValueError("a validation fold needs both training and test observations")
        if any(np.intersect1d(a, b).size for a, b in (
                (self.train, self.test), (self.train, self.excluded),
                (self.test, self.excluded))):
            raise ValueError("training, test and excluded observations must be disjoint")


def group_folds(groups):
    groups = np.asarray(groups)
    if groups.ndim != 1 or pd.isna(groups).any():
        raise ValueError("groups must be a one-dimensional array of nonmissing identifiers")
    if len(pd.unique(groups)) < 2:
        raise ValueError("need at least two independent groups for held-out validation")
    return tuple(ValidationFold(str(group), np.flatnonzero(groups != group),
                                np.flatnonzero(groups == group), np.array([], dtype=int))
                 for group in sorted(pd.unique(groups)))


def stressor_fold(dataset, held_out):
    labels = np.asarray(dataset.labels)
    test = np.asarray(dataset.mask(held_out), dtype=bool)
    if test.shape != labels.shape or not test.any():
        raise KeyError(f"{held_out!r} is not in this dataset")
    members = set(str(held_out).split("+"))
    excluded = np.array([members <= set(str(label).split("+")) for label in labels]) & ~test
    return ValidationFold(str(held_out), np.flatnonzero(~(test | excluded)),
                          np.flatnonzero(test), np.flatnonzero(excluded))


def subset_panel(dataset, indices):
    indices = np.asarray(indices)
    if indices.dtype == bool:
        if indices.shape != (len(dataset.labels),):
            raise ValueError("panel mask must match the observation count")
        indices = np.flatnonzero(indices)
    if (indices.ndim != 1 or indices.dtype.kind not in "iu" or np.any(indices < 0)
            or np.any(indices >= len(dataset.labels)) or len(np.unique(indices)) != len(indices)):
        raise ValueError("panel indices must be unique non-negative observation indices")
    updates = {name: np.asarray(getattr(dataset, name))[indices].copy()
               for name in ("readings", "labels", "doses", "modules")}
    for name in ("wells", "contexts", "read_times", "growth_rates"):
        values = getattr(dataset, name, None)
        if values is not None and len(values):
            if len(values) != len(dataset.labels):
                raise ValueError(f"{name} does not match the panel observation count")
            updates[name] = [values[int(i)] for i in indices]
    return replace(dataset, **updates)


def training_mean_skill(truth, prediction, training):
    truth, prediction, training = (np.asarray(a, dtype=float)
                                   for a in (truth, prediction, training))
    if truth.shape != prediction.shape or truth.ndim != 2 or training.ndim != 2:
        raise ValueError("truth and prediction must be matching matrices with a training matrix")
    if training.shape[1] != truth.shape[1] or not np.isfinite(training).all():
        raise ValueError("training targets must be finite with the same columns as truth")
    baseline = training.mean(axis=0)
    total = np.sum((truth - baseline) ** 2, axis=0)
    error = np.sum((truth - prediction) ** 2, axis=0)
    score = np.full(truth.shape[1], np.nan)
    usable = (total > 1e-12) & (np.ptp(truth, axis=0) > 1e-3)
    np.divide(error, total, out=score, where=usable)
    score[usable] = 1.0 - score[usable]
    return score


@dataclass(frozen=True)
class ModuleValidation:
    fold: ValidationFold
    n_states: int
    prediction: np.ndarray
    baseline: np.ndarray
    scores: dict[str, float]
    inner_scores: dict[int, float]
    selection: str


def validate_modules(dataset, held_out, max_states=6, width=None, seed=0):
    from .stress_model import train_stress_model

    outer = stressor_fold(dataset, held_out)
    training = subset_panel(dataset, outer.train)
    if not isinstance(max_states, (int, np.integer)) or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    inner_scores = {}
    if width is None:
        labels = [str(label) for label in sorted(set(training.labels)) if "+" not in str(label)]
        if len(labels) < 2:
            raise ValueError("nested width selection needs at least two training stressors")
        folds = [stressor_fold(training, label) for label in labels]
        ceiling = min(max_states, training.readings.shape[1],
                      min(len(fold.train) for fold in folds))
        for candidate in range(1, ceiling + 1):
            scores = []
            for fold in folds:
                inner = subset_panel(training, fold.train)
                model = train_stress_model(inner, candidate, seed)
                predicted = model.predict_modules(training.readings[fold.test]).to_numpy()
                values = training_mean_skill(training.modules[fold.test], predicted, inner.modules)
                finite = values[np.isfinite(values)]
                if finite.size:
                    scores.append(float(np.mean(finite)))
            inner_scores[candidate] = float(np.mean(scores)) if len(scores) == len(folds) else -np.inf
        if not inner_scores or not any(np.isfinite(value) for value in inner_scores.values()):
            raise ValueError("no width has a defined score on every inner training fold")
        width = max(inner_scores, key=inner_scores.get)
        selection = "nested leave-stressor-out on outer training only"
    else:
        if not isinstance(width, (int, np.integer)) or width < 1:
            raise ValueError("width must be a positive integer")
        selection = "fixed width supplied before outer scoring"
    model = train_stress_model(training, width, seed)
    prediction = model.predict_modules(dataset.readings[outer.test]).to_numpy()
    values = training_mean_skill(dataset.modules[outer.test], prediction, training.modules)
    return ModuleValidation(outer, width, prediction, training.modules.mean(axis=0),
                            dict(zip(model.modules, map(float, values))), inner_scores, selection)


def matched_product_branches(spec, states, selected_predictions):
    from ..pathway.flux import ProductCandidate, fit_product_candidate, predict_product_candidate

    if selected_predictions.state.duplicated().any() or states.condition.duplicated().any():
        raise ValueError("matched branch comparison needs unique condition identifiers")
    selected = selected_predictions.set_index("state")
    if set(selected.index) != set(states.condition):
        raise ValueError("selected predictions must account for every condition exactly once")
    rows = []
    for fold in group_folds(states.strain.to_numpy()):
        train, test = states.iloc[fold.train], states.iloc[fold.test]
        selection = selected.loc[test.condition]
        if selection.selected_model.nunique() != 1:
            raise ValueError("all held-out conditions of a strain must share one training-selected model")
        gene = selection.selected_gene.iloc[0]
        law = selection.selected_entry_law.iloc[0]
        for branch in ("partition", "saturating"):
            try:
                if pd.isna(gene):
                    raise ValueError("training selected a no-expression baseline; no entry branch to ablate")
                fitted = fit_product_candidate(spec, train, ProductCandidate(gene, law, branch))
                fit_failure = ""
            except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                fitted, fit_failure = None, str(exc)
            truth = test.reindex(columns=["q_betacarotene", "q_lycopene"]).to_numpy(dtype=float)
            # The fit is shared, but a prediction refusal belongs only to its condition.
            for i, row in enumerate(test.itertuples()):
                try:
                    if fitted is None:
                        raise ValueError(fit_failure)
                    query = test.iloc[[i]][["condition", "strain", "mu_per_h", gene]].copy()
                    predicted = predict_product_candidate(spec, fitted, query)[0]
                    status, reason, alpha = "ok", "", fitted.alpha
                except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                    predicted = np.full(2, np.nan)
                    status, reason, alpha = "failed", str(exc), np.nan
                valid = np.isfinite(predicted).all() and (predicted > 0).all()
                valid &= np.isfinite(truth[i]).all() and (truth[i] > 0).all()
                errors = np.log(predicted / truth[i]) if valid else np.full(2, np.inf)
                rows.append({
                    "condition": row.condition, "strain": row.strain,
                    "branch": branch, "selected_gene": gene, "selected_entry_law": law,
                    "alpha": alpha, "flux_predicted": float(predicted.sum()),
                    "product_predicted": predicted[0], "lycopene_rate_predicted": predicted[1],
                    "product_measured": truth[i, 0], "lycopene_rate_measured": truth[i, 1],
                    "joint_log_mse": float(np.mean(errors ** 2)),
                    "status": status if valid else "failed",
                    "failure_reason": reason or ("invalid observation or prediction" if not valid else ""),
                    "comparison": "same training-selected gene and entry law; only partition versus saturation changes",
                })
    return pd.DataFrame(rows)
