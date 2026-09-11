"""Test the panel's dose parameters against the uploaded plates.

The panel's EC50s are literature estimates. Two of them have measured dose ladders behind
them, so those two are testable: fit the biphasic response, cross-check with a shape-free
half-maximal crossing, and say whether the encoded value survived.

Usage: python scripts/run_calibration.py
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.generator.panel_calibration import calibrate_stressor
from ystwin.generator.stress_panel import STRESSORS

OUT = paths.outputs_dir()
MEASURED = OUT / "sensor_characterisation.csv"


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def main() -> None:
    if not MEASURED.exists():
        raise SystemExit(f"run scripts/run_sensor_characterisation.py first: {MEASURED} missing")
    measured = pd.read_csv(MEASURED)

    rule("what the plates say about each encoded EC50")
    rows = []
    for (construct, stressor), group in measured.groupby(["construct", "stressor"]):
        ladder = group.groupby("dose_mM")[["activity_late", "mu_late"]].mean().reset_index()
        verdict = calibrate_stressor(stressor, ladder.dose_mM, ladder.activity_late,
                                     ladder.mu_late)
        rows.append({
            "construct": construct, "stressor": stressor, "verdict": verdict.verdict,
            "encoded_ec50": verdict.encoded_ec50, "fitted_ec50": verdict.fitted_ec50,
            "half_maximal": verdict.half_maximal_dose, "r_squared": verdict.r_squared,
            "n_usable": verdict.n_usable, "n_total": len(ladder),
        })
        print(f"  {construct:<12} {stressor:<5} {verdict.verdict}")
        print(f"  {'':<12} {verdict.note}")
    table = pd.DataFrame(rows)

    rule("fitted against encoded, with the shape-free crossing beside it")
    print(table.drop(columns=["stressor"]).round(3).to_string(index=False))

    rule("what the panel now carries")
    for name in sorted({r["stressor"] for r in rows}):
        spec = STRESSORS[name]
        print(f"  {name:<6} EC50 {spec.ec50:<6.3g} lethal {spec.lethal_dose:<6.3g} "
              f"({spec.lethal_dose / spec.ec50:.0f}x window)")
        print(f"  {'':<6} {spec.source.splitlines()[0]}")

    rule("what would settle what is still open")
    dead = table[table.n_usable < table.n_total]
    for _, row in dead.iterrows():
        print(f"  {row.construct}: {row.n_total - row.n_usable} of {row.n_total} doses had no "
              f"growing culture, so the falling limb was never measured")
    print("  peroxide needs doses between 1 and 2 mM, where the plates jump straight from")
    print("  a live culture to a dead one and the lethal dose hides in the gap")

    table.to_csv(OUT / "panel_calibration.csv", index=False)
    print(f"\nwrote {OUT / 'panel_calibration.csv'}")


if __name__ == "__main__":
    main()
