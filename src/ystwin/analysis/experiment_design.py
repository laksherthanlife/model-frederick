"""Choosing which stressor combinations to run, scored on transfer rather than a proxy.

Analytical measures of design coverage failed repeatedly here. Eight single-agent rays
already span all seven modules, so any full-span score ranks every design identically;
restricting to the fitted subspace ranks them but is not monotonic in adding treatments,
because a new treatment can rotate that subspace away from a stressor it used to reach.

Transfer is simulable, so it is measured instead of approximated: run the candidate
design, hold each stressor out in turn, take the median R2. Each candidate costs a full
leave-one-out sweep, which is the price of scoring the thing actually wanted.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations as _pairs

import numpy as np

from ..generator.panel_experiment import panel_dataset
from ..generator.stress_panel import REPORTERS, STRESSORS
from .transfer import leave_one_stressor_out

__all__ = ["RECOMMENDED_DESIGN", "design_transfer", "select_combinations"]

_OBSERVED = (0, 1, 2, 3)


def design_transfer(
    treatments,
    n_states: int = 3,
    reporters=None,
    doses=(0.25, 0.5, 1.0, 2.0),
    replicates: int = 2,
    noise_cv: float = 0.02,
    seed: int = 0,
    stressors=None,
    observed=_OBSERVED,
) -> float:
    """Median leave-one-stressor-out transfer R2 a design delivers.

    Args:
        treatments: Stressor tuples co-dosed alongside the single-agent series.
        n_states: Latent dimension the fitted model carries.
        reporters: Reporters read; defaults to the whole panel.
        doses: Multiples of each stressor's EC50 applied.
        replicates: Wells per treatment and dose.
        noise_cv: Measurement noise, held fixed so designs are compared on one draw.
        seed: Seed for the noise.
        stressors: Single agents dosed; defaults to all of them. Each one costs two latent
            fits in the sweep, so a broad landscape wants this narrowed.
        observed: Channel indices revealed for the held-out stressor.
    """
    treatments = [tuple(t) for t in treatments]
    unknown = sorted({s for t in treatments for s in t if s not in STRESSORS})
    if unknown:
        raise KeyError(f"no stressor(s) {unknown}; have {sorted(STRESSORS)}")
    single = list(stressors or STRESSORS)
    data = panel_dataset(
        reporters=list(reporters or REPORTERS), stressors=single, doses=doses,
        replicates=replicates, noise_cv=noise_cv, seed=seed, combinations=treatments)
    table = leave_one_stressor_out(data, n_states=n_states, observed=observed, only=single)
    return float(table["r2"].median())


def _single_threaded():
    """Pin each worker to one BLAS thread.

    Workers are spawned per core and every one of them would otherwise start its own pool
    of linear algebra threads, oversubscribing the machine several times over. That is
    slower at best and deadlocks at worst, which is what it did here.
    """
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"


def _score_candidate(job):
    treatments, scoring = job
    return design_transfer(treatments, **scoring)


def select_combinations(n_pairs: int, library=None, workers: int = 1, **scoring) -> list[tuple[str, str]]:
    """Greedily add the pair that most improves simulated transfer.

    Args:
        n_pairs: How many co-dosed pairs the plate can hold.
        library: Stressors available to pair; defaults to all of them.
        workers: Processes scoring candidates at once. Candidates within a round are
            independent, so this changes only how long the search takes. Worker processes
            are spawned, which re-imports the calling module, so a script using this must
            guard its entry point with ``if __name__ == "__main__"`` or the search restarts
            inside every worker.
        **scoring: Passed to `design_transfer`, so the design is chosen under the same
            reporters, doses and noise the experiment will actually run.
    """
    names = list(library or scoring.get("stressors") or STRESSORS)
    available = list(_pairs(names, 2))
    if n_pairs > len(available):
        raise ValueError(f"asked for {n_pairs} pairs but only {len(available)} exist")

    chosen: list[tuple[str, str]] = []
    for _ in range(n_pairs):
        candidates = [p for p in available if p not in chosen]
        jobs = [(chosen + [pair], scoring) for pair in candidates]
        if workers > 1:
            with ProcessPoolExecutor(max_workers=workers,
                                     initializer=_single_threaded) as pool:
                scores = list(pool.map(_score_candidate, jobs))
        else:
            scores = [_score_candidate(job) for job in jobs]
        chosen.append(candidates[int(np.argmax(scores))])
    return chosen


RECOMMENDED_DESIGN = [
    ("menadione", "diamide"),
    ("H2O2", "diamide"),
    ("H2O2", "menadione"),
    ("DTT", "MMS"),
    ("calcium_chloride", "rapamycin"),
    ("diamide", "copper_sulfate"),
    ("diamide", "rapamycin"),
    ("menadione", "NaCl"),
]
"""Eight co-dosed pairs, searched over all 25 stressors and 300 candidates.

Derived after the generator was corrected: per-module potency, a dose ladder that stays
inside the healthy range, the constitutive floor every promoter sits on, and a noise split
calibrated so a simulated plate reproduces the measured activity spread. The previous set
was derived before any of that and no longer holds.

Scored on the full panel they lift median transfer from 0.038 to 0.429. Crossing all 300
pairs reaches only 0.161 -- less than half, at nearly forty times the wells -- because the
fitted subspace follows whatever the design samples most, and redundant mixtures crowd it
with directions already covered.

Oxidants dominate the list because per-module potency gave them the most to separate: a
peroxide moves its pool, its regulon and the proteasome at three different concentrations,
so crossing two of them visits directions no single agent reaches. Ranking was done on a
representative ten held out; the numbers above are scored on all twenty-five and on fresh
noise seeds, so none of them is the number the search optimised.
"""
