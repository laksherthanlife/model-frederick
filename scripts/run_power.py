"""How many biological replicates each design needs, and what the design buys.

Usage: python scripts/run_power.py
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.power import (
    DEFAULT_WELL_CV,
    discrimination_power,
    replicates_needed,
)
from ystwin.generator.design import collinearity, decoupling_grid, stressor_only_series
from ystwin.generator.literature import parameters_from_literature

OUT = paths.outputs_dir()
DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0)
FOLDS = (1.2, 1.5, 2.0, 3.0)
SIMS = 120


def main() -> None:
    params = parameters_from_literature("UPRE2")
    designs = {
        "dose only (as run)": stressor_only_series(params, doses=DOSES),
        "dose x 2 nutrients": decoupling_grid(params, doses=DOSES, nutrient_factors=(1.0, 0.5)),
        "dose x 3 nutrients": decoupling_grid(params, doses=DOSES, nutrient_factors=(1.0, 0.6, 0.35)),
    }

    print("=" * 84)
    print("Design collinearity  (1.0 = stress and growth are the same axis)")
    print("=" * 84)
    for name, conditions in designs.items():
        print(f"  {name:22s} {len(conditions):>3} conditions   collinearity {collinearity(conditions):.3f}")

    print("\n" + "=" * 84)
    print("Replicates to DETECT a dose slope at 80% power  (the easy question)")
    print("=" * 84)
    rows = []
    print(f"  {'design':22s}" + "".join(f"{f'{f}x':>10}" for f in FOLDS))
    for name, conditions in designs.items():
        cells = []
        for fold in FOLDS:
            n = replicates_needed(params, conditions, fold, target_power=0.8,
                                  max_replicates=16, n_simulations=SIMS, seed=0)
            cells.append("none" if n is None else str(n))
            rows.append({"design": name, "question": "detect", "fold": fold,
                         "replicates": n, "collinearity": collinearity(conditions)})
        print(f"  {name:22s}" + "".join(f"{c:>10}" for c in cells))

    print("\n" + "=" * 84)
    print("Power to ATTRIBUTE it to dose with growth also in the model  (the real question)")
    print("=" * 84)
    print(f"  {'design':22s}" + "".join(f"{f'{f}x':>10}" for f in FOLDS) + f"{'  (n=6 plates)':>16}")
    for name, conditions in designs.items():
        cells = []
        for fold in FOLDS:
            p = discrimination_power(params, conditions, fold, n_replicates=6,
                                     n_simulations=SIMS, seed=0)
            cells.append(f"{p:.0%}")
            rows.append({"design": name, "question": "attribute", "fold": fold,
                         "power_n6": float(p), "power_n6_low": p.low,
                         "power_n6_high": p.high, "n_simulations": p.n_simulations,
                         "well_cv": DEFAULT_WELL_CV,
                         "collinearity": collinearity(conditions)})
        print(f"  {name:22s}" + "".join(f"{c:>10}" for c in cells))

    print("\n" + "=" * 84)
    print("What a chemostat buys  (mu known exactly rather than estimated from OD)")
    print("=" * 84)
    grid = designs["dose x 3 nutrients"]
    print(f"  {'effect':>8}{'plate':>10}{'chemostat':>12}")
    for fold in FOLDS:
        plate = discrimination_power(params, grid, fold, n_replicates=4,
                                     n_simulations=SIMS, seed=1)
        chem = discrimination_power(params, grid, fold, n_replicates=4,
                                    n_simulations=SIMS, seed=1, growth_known=True)
        print(f"  {f'{fold}x':>8}{plate:>10.0%}{chem:>12.0%}")
        rows.append({"design": "dose x 3 nutrients", "question": "chemostat", "fold": fold,
                     "power_plate_n4": float(plate), "power_plate_n4_low": plate.low,
                     "power_plate_n4_high": plate.high,
                     "power_chemostat_n4": float(chem),
                     "power_chemostat_n4_low": chem.low,
                     "power_chemostat_n4_high": chem.high,
                     "n_simulations": plate.n_simulations,
                     "well_cv": DEFAULT_WELL_CV})

    pd.DataFrame(rows).to_csv(OUT / "power_analysis.csv", index=False)
    print(f"\nwrote {OUT / 'power_analysis.csv'}")


if __name__ == "__main__":
    main()
