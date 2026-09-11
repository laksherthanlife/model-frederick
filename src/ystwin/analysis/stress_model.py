"""The trained artefact: reporter readings in, named module activities out.

Two pieces are fitted, and they are fitted differently. The latent basis comes from the
readings alone and is unsupervised, so it transfers to any plate read on the same channels.
The map from that basis to named modules is supervised, and simulation is the only place
the labels exist to fit it -- which is what a mechanistic generator is for. Once fitted the
readout applies to real readings, which carry no labels at all.

Inference is MAP under the prior the training states define, not least squares. A plate
outside the training range is then pulled toward the training mean rather than producing
an enormous state that happens to match the channels it was shown.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..generator.stress_panel import MODULES
from .latent import fit_latent

# Span of held-out truth below which a module is not exercised and cannot be recovered.
# Module activities are fractions of full induction, so this is 0.1% of the scale.
_RECOVERABLE_SPAN = 1e-3

__all__ = [
    "StressModel",
    "module_transfer",
    "select_width_for_modules",
    "train_stress_model",
]


@dataclass(frozen=True)
class StressModel:
    """A fitted stress state, with the channels it expects and how well it recovers each module."""

    reporters: list[str]
    modules: list[str]
    centre: np.ndarray = field(repr=False)
    loadings: np.ndarray = field(repr=False)
    prior: np.ndarray = field(repr=False)
    noise: float
    readout: np.ndarray = field(repr=False)
    recovery: dict[str, float] = field(repr=False)

    def _check(self, readings: np.ndarray) -> np.ndarray:
        data = np.atleast_2d(np.asarray(readings, dtype=float))
        if data.shape[1] != len(self.reporters):
            raise ValueError(
                f"model expects {len(self.reporters)} channels {self.reporters}, "
                f"got {data.shape[1]}")
        return data

    def infer(self, readings) -> np.ndarray:
        """MAP latent state for each reading, shrunk toward the training mean."""
        data = self._check(readings)
        precision = (self.loadings @ self.loadings.T) / self.noise + np.diag(
            1.0 / np.maximum(self.prior, 1e-12))
        target = self.loadings @ (data - self.centre).T / self.noise
        return np.linalg.solve(precision, target).T

    def predict_modules(self, readings) -> pd.DataFrame:
        """Estimated activity of every module, one row per reading."""
        states = self.infer(readings)
        design = np.column_stack([states, np.ones(len(states))])
        return pd.DataFrame(design @ self.readout, columns=self.modules)

    def save(self, path) -> None:
        np.savez(
            path, reporters=np.array(self.reporters), modules=np.array(self.modules),
            centre=self.centre, loadings=self.loadings, prior=self.prior,
            noise=np.array([self.noise]), readout=self.readout,
            recovery=np.array([self.recovery[m] for m in self.modules]))

    @classmethod
    def load(cls, path) -> StressModel:
        blob = np.load(path, allow_pickle=False)
        modules = [str(m) for m in blob["modules"]]
        return cls(
            reporters=[str(r) for r in blob["reporters"]], modules=modules,
            centre=blob["centre"], loadings=blob["loadings"], prior=blob["prior"],
            noise=float(blob["noise"][0]), readout=blob["readout"],
            recovery=dict(zip(modules, (float(v) for v in blob["recovery"]))))


def train_stress_model(dataset, n_states: int, seed: int = 0) -> StressModel:
    """Fit a latent stress state and the readout that names its modules.

    Args:
        dataset: Panel readings carrying the true module activities behind them.
        n_states: Width of the latent state.
        seed: Seed for the fit.
    """
    readings = np.asarray(dataset.readings, dtype=float)
    if n_states > readings.shape[1]:
        raise ValueError(f"n_states ({n_states}) cannot exceed channels ({readings.shape[1]})")

    # nanmean, not mean: fit_latent below masks missing readings and documents that it
    # does, so a plain mean here would throw that away -- one NaN anywhere in a channel
    # makes its centre NaN and the channel useless. Measured on the Gasch arrays, this
    # single character was 81 arrays x 17 channels against 96 x 19.
    centre = np.nanmean(readings, axis=0)
    fit = fit_latent(readings - centre, n_states=n_states, seed=seed)
    design = np.column_stack([fit.states, np.ones(len(fit.states))])
    readout, *_ = np.linalg.lstsq(design, dataset.modules, rcond=None)

    residual = dataset.modules - design @ readout
    recovery = {}
    for index, name in enumerate(MODULES):
        truth = dataset.modules[:, index] - dataset.modules[:, index].mean()
        total = float(truth @ truth)
        left = float(residual[:, index] @ residual[:, index])
        recovery[name] = 0.0 if total <= 1e-12 else float(np.clip(1.0 - left / total, 0.0, 1.0))

    return StressModel(
        reporters=list(dataset.reporters), modules=list(MODULES), centre=centre,
        loadings=fit.loadings, prior=np.var(fit.states, axis=0),
        noise=max(fit.heldout_error, 1e-12), readout=readout, recovery=recovery)


def module_transfer(dataset, held_out: str, n_states: int, seed: int = 0) -> dict[str, float]:
    """Train without one stressor, then score the module activities it reports for it.

    Stronger than filling in a channel, and in the units that matter: trained on other
    stressors, does the model report the right biology when an unseen one arrives? Both
    stages are blind to it -- the latent basis and the readout are fitted without it.

    Scored against predicting the training mean of each module, so a module the training
    stressors never excited scores zero rather than borrowing credit from the intercept.

    A module the held-out stressor leaves untouched scores zero too, and that guard is
    load-bearing. Its truth is a constant, so any prediction nearer that constant than the
    training mean beats the baseline without tracking anything -- and the training mean is
    far from it precisely when OTHER stressors drive the module hard. Held out caffeine,
    which touches cell_wall, nitrogen and ESR, the unguarded score read 0.64 for peroxide,
    0.62 for redox and 0.60 for oxidative, all three identically zero in every held-out
    well. Shrinking toward zero is free; paying for it inverts the metric into a measure of
    how strong the rest of the panel is.
    """
    mask = dataset.mask(held_out)
    if not mask.any():
        raise KeyError(f"{held_out!r} is not in this dataset; have {sorted(set(dataset.labels))}")

    training = _subset(dataset, ~mask)
    model = train_stress_model(training, n_states=n_states, seed=seed)
    predicted = model.predict_modules(dataset.readings[mask]).to_numpy()
    truth = dataset.modules[mask]
    baseline = training.modules.mean(axis=0)

    scores = {}
    for index, name in enumerate(model.modules):
        total = float(np.sum((truth[:, index] - baseline[index]) ** 2))
        left = float(np.sum((truth[:, index] - predicted[:, index]) ** 2))
        if total <= 1e-12 or float(np.ptp(truth[:, index])) <= _RECOVERABLE_SPAN:
            scores[name] = 0.0
            continue
        scores[name] = float(np.clip(1.0 - left / total, 0.0, 1.0))
    return scores


def _subset(dataset, mask):
    from ..generator.panel_experiment import PanelDataset

    return PanelDataset(
        readings=dataset.readings[mask], labels=dataset.labels[mask],
        doses=dataset.doses[mask], modules=dataset.modules[mask],
        reporters=list(dataset.reporters))


def select_width_for_modules(dataset, max_states: int, sample: int = 4, seed: int = 0) -> int:
    """Latent width that best reports module activities for a stressor never seen.

    Width has to be chosen on the task the model is deployed for. Selecting on channel
    prediction picks one state for a narrow build, and one state cannot hold a general
    stress response and a metabolite pool at the same time, so the pool is reported as
    nothing however good the sensor reading it is.

    Args:
        dataset: Panel readings carrying the true module activities behind them.
        max_states: Widest state to consider; capped at the channels available.
        sample: Stressors held out to score each width. A full sweep is unnecessary here
            and costs a model fit per stressor per width.
        seed: Seed for the fits.
    """
    labels = sorted(set(dataset.labels))
    held = labels[:: max(1, len(labels) // max(sample, 1))][:sample] or labels[:1]
    best, best_score = 1, -np.inf
    for n_states in range(1, min(max_states, dataset.readings.shape[1]) + 1):
        scores = [v for name in held
                  for v in module_transfer(dataset, name, n_states, seed).values()]
        mean = float(np.mean(scores)) if scores else -np.inf
        if mean > best_score:
            best, best_score = n_states, mean
    return best
