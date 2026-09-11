"""Transfer: train on some stressors, predict one never seen -- the test of "general".

A stress state that only re-describes its training stressors is a lookup table. The claim
worth making is that a state learned from some stressors predicts reporters under another,
and it is checked by withholding a stressor, revealing only some of its channels and
predicting the rest, against an oracle that was allowed to see it.

Usage: python scripts/run_transfer.py [--scale]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import dataclasses
import itertools
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
from ystwin.analysis.latent import select_dimension
from ystwin.analysis.recovery import module_recovery
from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD
from ystwin.analysis.nulls import (
    _MIN_DRAWS,
    compare_to_null,
    matched_marginals,
    rotated_subspace,
    skill_score_if_defined,
)
from ystwin.analysis.transfer import leave_one_stressor_out, select_dimension_by_transfer
from ystwin.generator.panel_experiment import (
    MEASURED_ACTIVITY_CV,
    MEASURED_GROWTH_RATE_SE,
    panel_dataset,
)
from ystwin.generator.stress_panel import (
    MODULES,
    REPORTERS,
    STRESSORS,
)

OUT = paths.outputs_dir()
PANEL = ["DTT", "H2O2", "heat", "NaCl", "glucose_starvation", "BPS", "MG132",
         "MMS", "congo_red", "calcium_chloride", "rapamycin", "copper_sulfate"]
READERS = ["STRE-general", "UPRE-ER", "TRX2-oxidative", "HSE-heat",
           "STRE-osmotic", "PACE-proteasome", "FeRE-iron", "Xbox-dna"]
ALL_PAIRS = list(itertools.combinations(PANEL, 2))

NULL_DRAWS = 40
"""Surrogate draws per null. Each costs a full leave-one-out sweep, so the
p-value floor here is 1/41 = 0.024 rather than anything smaller."""


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def _validate_workers(workers) -> None:
    if (not isinstance(workers, (int, np.integer)) or isinstance(workers, (bool, np.bool_))
            or workers < 1):
        raise ValueError("workers must be a positive integer")


def transfer_table(data, n_states, observed, *, workers: int = 1):
    _validate_workers(workers)
    if workers == 1:
        return leave_one_stressor_out(data, n_states=n_states, observed=observed, only=PANEL)
    held_out = [name for name in sorted(set(data.labels)) if name in PANEL]
    if len(held_out) < 2:
        return leave_one_stressor_out(data, n_states=n_states, observed=observed, only=PANEL)

    def fold_table(name):
        return leave_one_stressor_out(data, n_states=n_states, observed=observed, only=[name])

    with ThreadPoolExecutor(max_workers=min(workers, len(held_out))) as executor:
        return pd.concat(executor.map(fold_table, held_out), ignore_index=True)


def main(scale: bool, output_dir=None, draws: int = NULL_DRAWS, transfer_only: bool = False,
         *, workers: int = 1) -> None:
    if not transfer_only and (not isinstance(draws, (int, np.integer)) or draws < _MIN_DRAWS):
        raise ValueError(f"null comparisons require at least {_MIN_DRAWS} draws")
    _validate_workers(workers)
    out = pathlib.Path(output_dir) if output_dir is not None else OUT
    out.mkdir(parents=True, exist_ok=True)
    doses = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0) if scale else (0.25, 0.5, 1.0, 2.0)
    replicates = 6 if scale else 3
    seeds = range(5) if scale else range(2)
    full = READERS
    build = RECOMMENDED_BUILD.stress_reporters

    rule("how many dimensions does the stress landscape actually have?")
    data = panel_dataset(reporters=full, stressors=PANEL, doses=doses,
                         replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV, growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0,
                         combinations=RECOMMENDED_DESIGN)
    centred = data.readings - data.readings.mean(axis=0)
    share = np.linalg.svd(centred, compute_uv=False) ** 2
    share = share / share.sum()
    print(f"  variance per component   : {np.round(share, 3)}")
    print(f"  best rank for imputing a gap : {select_dimension(centred, max_states=6)}")
    k = select_dimension_by_transfer(data, max_states=4, observed=(0, 1, 2, 3))
    print(f"  full-data exploratory rank   : {k} (not used to score outer folds)")
    print(f"  {len(data.labels)} wells over {len(set(data.labels))} treatments")
    print(f"  the landscape has {len(MODULES)} state variables and {len(STRESSORS)} stressors")
    print("  Every scored outer fold reselects its width using only its training treatments.")
    print("  Co-doses containing the held-out agent are excluded from its training set.")
    if transfer_only:
        table = transfer_table(data, None, (0, 1, 2, 3), workers=workers)
        table.to_csv(out / "transfer_by_stressor.csv", index=False)
        print(table.to_string(index=False))
        print("Transfer-only run: no design comparison, null test or descriptive recovery table was requested.")
        print(f"wrote {out / 'transfer_by_stressor.csv'}")
        return

    rule("does the experimental design decide whether transfer works?")
    # Sizes are read off the designs, not written into the labels: PANEL and
    # RECOMMENDED_DESIGN have both grown since these were last typed by hand.
    designs = {"blocked, one at a time": [],
               f"all {len(ALL_PAIRS)} pairs": ALL_PAIRS,
               f"recommended {len(RECOMMENDED_DESIGN)} pairs": RECOMMENDED_DESIGN}
    rows = []
    for tag, design in designs.items():
        scores = [float(transfer_table(panel_dataset(
            reporters=full, stressors=PANEL, doses=doses, replicates=2,
            noise_cv=MEASURED_ACTIVITY_CV, seed=s, combinations=design),
            None, (0, 1, 2, 3), workers=workers).r2.median()) for s in seeds]
        wells = len(panel_dataset(reporters=full, stressors=PANEL, doses=doses,
                                  replicates=replicates, noise_cv=0.0,
                                  combinations=design).labels)
        rows.append({"design": tag, "wells": wells, "transfer_r2": np.mean(scores),
                     "sd": np.std(scores)})
    print(pd.DataFrame(rows).round(3).to_string(index=False))

    rule("per-stressor transfer under the recommended design, nested width selection")
    table = transfer_table(data, None, (0, 1, 2, 3), workers=workers)
    scored_transfer = table.copy()
    print(table.round(3).to_string(index=False))

    rule("and is that better than no structure at all?")
    # A transfer R2 alone is not evidence. Any rank-k projection of correlated channels
    # captures shared variance whether or not its axes mean anything, so the number to
    # beat is what arbitrary axes of the same rank score on the same data -- not zero.
    observed_median = float(np.nanmedian(table.r2))

    def median_transfer(dataset) -> float:
        return float(np.nanmedian(transfer_table(dataset, None, (0, 1, 2, 3), workers=workers).r2))

    def rotate(dataset, rng):
        return dataclasses.replace(dataset, readings=rotated_subspace(dataset.readings, rng))

    def scramble(dataset, rng):
        return dataclasses.replace(dataset, readings=matched_marginals(dataset.readings, rng))

    rotated_null = None
    for surrogate, label in ((rotate, "random axes, nested rank"),
                             (scramble, "no cross-channel structure")):
        result = compare_to_null(median_transfer, data, surrogate,
                                 n_draws=draws, greater_is_better=True,
                                 label=label, seed=0)
        if surrogate is rotate:
            rotated_null = result.null_median
        verdict = "BEATS" if result.beats_null else "does NOT beat"
        print(f"  {label:32s} null median {result.null_median:+.4f}  "
              f"observed {result.observed:+.4f}  p={result.p_value:.3f}  {verdict}")
    skill = skill_score_if_defined(observed_median, rotated_null, perfect=1.0)
    print(f"  skill over the rotated null: "
          f"{'undefined -- the null already predicts perfectly' if skill is None else f'{skill:+.3f}'} "
          "(1.0 = closes the whole gap to perfect, negative = worse than no structure)")

    rule("does it survive the five-channel build that can actually be made?")
    built = panel_dataset(reporters=build, stressors=PANEL, doses=doses,
                          replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV, growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0,
                          combinations=RECOMMENDED_DESIGN)
    kinds = {r: REPORTERS[r].kind.value for r in build}
    print(f"  build: {kinds} + {RECOMMENDED_BUILD.reference_fluorophore} reference")
    print(f"  spectrally buildable: {RECOMMENDED_BUILD.buildable}")
    print("  every channel held out in turn, because which one is predicted decides the")
    print("  score: a channel no stressor in the panel drives cannot be predicted by anyone")
    rows = []
    for hidden in range(len(build)):
        seen = tuple(i for i in range(len(build)) if i != hidden)
        table = transfer_table(built, None, seen, workers=workers)
        rows.append({"predicted": build[hidden], "transfer_r2": table["r2"].median(),
                     "oracle_r2": table["oracle_r2"].median()})
    frame = pd.DataFrame(rows)
    print(frame.round(3).to_string(index=False))
    print(f"  averaged over channels: transfer {frame['transfer_r2'].mean():.3f}"
          f"   oracle {frame['oracle_r2'].mean():.3f}")

    rule("do the latent states carry the modules, or only compress them?")
    recovery = pd.DataFrame({
        "module": list(MODULES),
        "full panel (7)": [module_recovery(data, n_states=5)[m] for m in MODULES],
        "built panel (4)": [module_recovery(built, n_states=3)[m] for m in MODULES],
        "shuffled control": [module_recovery(data, n_states=5, shuffle=True)[m] for m in MODULES],
    })
    print(recovery.round(3).to_string(index=False))

    scored_transfer.to_csv(out / "transfer_by_stressor.csv", index=False)
    recovery.to_csv(out / "module_recovery.csv", index=False)
    print(f"\nwrote {out/'transfer_by_stressor.csv'} and {out/'module_recovery.csv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", action="store_true")
    parser.add_argument("--output-dir", type=pathlib.Path, default=None)
    parser.add_argument("--draws", type=int, default=NULL_DRAWS)
    parser.add_argument("--transfer-only", action="store_true")
    parser.add_argument("--workers", type=int, default=1,
                        help="maximum threads for independent outer folds (default: 1)")
    args = parser.parse_args()
    main(scale=args.scale, output_dir=args.output_dir, draws=args.draws,
         transfer_only=args.transfer_only, workers=args.workers)
