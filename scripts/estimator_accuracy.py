"""Measure the activity inversion's error at the geometry it actually runs at.

    python scripts/estimator_accuracy.py

Writes ``outputs/estimator_accuracy.csv``. Needs no wet-lab data: the plates supply the
geometry, growth and noise as constants in :mod:`ystwin.analysis.estimator_accuracy`, and
everything below is simulated from them.

This script exists because the number it produces did not have one. `reporter.py` chose the
smoothing window for every dilution-corrected result in this repository, and justified the
choice with a five-cell table in its own docstring that no script wrote and no test read.
The table was also taken over a 24-hour run at 145 points; the committed exports are 25
points over 4.00 hours, where the old rule returned 0.667 h -- below the smallest cell in
the table it was resting on.

Three findings, in the order they matter:

1. **The error is set by the window's duration in hours, not by its share of the run and
   not by its point count.** Two geometries six times apart in density agree at matched
   duration and disagree at matched points. A rule that sets the window as a fraction of
   the run is therefore right only at the run length it was tuned on.
2. **A fold change is accurate at every window tried** -- under 0.6% low across a four-fold
   range of them -- because a fold is a ratio of two wells on one plate and the smoother
   treats numerator and denominator alike. What the window buys is precision, and precision
   is what reaches an interval.
3. **Ignoring chromophore maturation costs about 10% on a single well's activity and almost
   nothing on a fold**, for the same reason. That assumption had never been costed;
   `promoter_activity_from_total` refuses maturation outright, so no committed number could
   ever have carried a correction for it.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.estimator_accuracy import (
    DOCSTRING_REGIME,
    MATURATION_HALF_LIFE_H,
    PLATE,
    Geometry,
    accuracy_curve,
    fold_recovery,
    maturation_bias,
)

# A sparse-long and a dense-short run, included only to separate two explanations that the
# two headline geometries cannot: they differ in BOTH density and duration, so agreement
# between them is consistent with either variable governing the error. These two cross the
# design -- same durations, densities swapped -- and settle it.
CONTROLS = (
    Geometry("sparse and long", 25, 24.0),
    Geometry("dense and short", 145, 4.0),
)

REPLICATES = 200


def main() -> int:
    rows: list[dict] = []

    for geometry in (PLATE, DOCSTRING_REGIME, *CONTROLS):
        print(f"\n{geometry.summary()}")
        for point in accuracy_curve(geometry, n_replicates=REPLICATES):
            print(f"    {point.summary()}")
            rows.append({"measurement": "accuracy_curve", **point.__dict__})

    print("\nWhat a window buys on a FOLD, at the committed-plate geometry.")
    print("The fold, not the activity, is what every table in outputs/ holds.")
    for window in (0.667, 1.0, 1.333, 2.0):
        got = fold_recovery(PLATE, window_h=window, n_replicates=REPLICATES)
        print(f"    {window:6.3f} h  target {got['target_fold']:.4f}  "
              f"recovered {got['recovered_fold']:.4f} ({got['relative_error']:+.2%})  "
              f"sd {got['spread_of_recovered_fold']:.4f}  naive {got['naive_fold']:.3f}")
        rows.append({"measurement": "fold_recovery", **got})

    print("\nWhat ignoring chromophore maturation costs.")
    for geometry in (PLATE, DOCSTRING_REGIME):
        got = maturation_bias(geometry, n_replicates=REPLICATES)
        print(f"    {geometry.label}: a {MATURATION_HALF_LIFE_H * 60:.0f}-minute maturation is "
              f"{got['half_life_as_fraction_of_run']:.1%} of the run and adds "
              f"{got['added_by_maturation']:+.1%} to one well's late activity")
        rows.append({"measurement": "maturation_bias", **got})
        paired = fold_recovery(geometry, half_life_h=MATURATION_HALF_LIFE_H,
                               n_replicates=REPLICATES)
        clean = fold_recovery(geometry, n_replicates=REPLICATES)
        print(f"      and {paired['relative_error'] - clean['relative_error']:+.2%} to a FOLD, "
              f"because a ratio cancels what both wells share")
        rows.append({"measurement": "fold_recovery_with_maturation", **paired})

    out = paths.outputs_dir() / "estimator_accuracy.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
