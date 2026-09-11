"""Do the ER sensors see oxidative stress, or the oxidative sensors see ER stress?

    python3 scripts/score_sensor_crosstalk.py

The reciprocal experiment, and the answer is **no in both directions** -- which is the
useful answer, because a sensor that responds to the wrong stressor cannot be read.

**The raw fluorescence says otherwise, and that is the trap this repository exists for.**
UPRE1 under H2O2 falls from 1774 to 194 RFU, a 9.1-fold drop that reads as a large response.
Over the same doses the culture's OD600 falls from 0.320 to 0.145 -- **45% of control**. Per
cell the signal is flat to 1 mM and then follows the death, and `docs/FINDINGS.md` section 1
is the same finding on the dose-response plates.

**It is three biological replicates, and it was one until 2026-08-29.** The table this reads
held the workbook's first `Summary` block only, while `SOURCE.md` described it as a mean over
three. It now holds all three, extracted by `scripts/extract_crosstalk_workbook.py` and
checked value-for-value against the raw instrument export. That matters here more than it
usually would: on replicate 1 UPRE1 is at 90% viability at 1 mM H2O2, and on the other two it
is at 40% and 31% -- the viability gate below passed a dose on the one plate where the cells
were alive.

**The on-target comparison is the point, and it was missing entirely.** "The sensor stayed
quiet for the wrong stressor" is weak on its own, because a dead sensor also stays quiet.
`data/crosstalk/erox_2026-08_on_vs_off_target.tsv` puts each sensor's response to its OWN
stressor beside its response to the other one, over the same three replicates.

**What this cannot do.** It is an ENDPOINT read, so the dilution correction that produced
the numbers in `outputs/late_window_sensitivity.csv` cannot be applied -- that inverts the
reporter ODE using a growth rate from an OD trace, and there is one timepoint per well here.
The per-cell ratio removes how many cells there are; it does not remove how long a slowed
culture accumulated reporter for. So this supports a negative -- a flat ratio cannot hide a
response -- and would not support a positive.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths

# Below this fraction of the undosed OD, a culture is dying and a signal ratio describes the
# dying rather than the sensing. Not a fitted threshold -- a reading aid, and stated so.
VIABILITY_FLOOR = 0.6


def load() -> pd.DataFrame:
    path = paths.data_dir() / "crosstalk" / "erox_2026-08_endpoint.tsv"
    if not path.exists():
        raise SystemExit(f"{path} is tracked and missing; bad checkout")
    return pd.read_csv(path, sep="\t")


def scored(states: pd.DataFrame) -> pd.DataFrame:
    """Fold and viability per dose, formed WITHIN each plate and then combined.

    The same discipline `analysis/uncertainty.py::_combine` applies to the kinetic plates,
    and for the same reason: a plate-level shift -- inoculum, medium batch, reader gain --
    multiplies a plate's dosed and control wells alike and cancels in a within-plate ratio,
    while a ratio of values pooled across plates imports it. Here it is not hypothetical.
    The three control ODs are 0.320, 0.517 and 0.488, so pooling would let the plate with
    the thinnest culture set the denominator.

    `n_replicates` and the spread are carried on every row because the previous version of
    this table had one plate and no way to say so.
    """
    rows = []
    for (construct, stressor, dose), group in states.groupby(
            ["construct", "stressor", "dose_mM"]):
        controls = states[(states.construct == construct) & (states.dose_mM == 0.0)]
        per_plate = []
        for r in group.itertuples():
            control = controls[controls.replicate == r.replicate]
            if control.empty:
                continue
            control = control.iloc[0]
            per_plate.append({
                "viability": r.od600_mean / control.od600_mean,
                "raw_fold": r.rfu_mean / control.rfu_mean,
                "per_cell_fold": r.per_cell / control.per_cell,
            })
        if not per_plate:
            continue
        frame = pd.DataFrame(per_plate)
        rows.append({
            "construct": construct, "stressor": stressor, "dose_mM": dose,
            "n_replicates": len(frame),
            "od600": float(group.od600_mean.mean()),
            "viability": float(frame.viability.mean()),
            "viability_min": float(frame.viability.min()),
            "raw_fold": float(frame.raw_fold.mean()),
            "per_cell_fold": float(frame.per_cell_fold.mean()),
            "per_cell_fold_half_sd": float(frame.per_cell_fold.std(ddof=0) / 2),
            # The WORST plate decides, not the mean. A dose where one plate's culture has
            # collapsed is a dose where the reading is partly about death, and averaging
            # that away is how the single-replicate version passed 1 mM H2O2 for UPRE1.
            "culture_healthy": bool(frame.viability.min() >= VIABILITY_FLOOR),
        })
    return pd.DataFrame(rows).sort_values(["construct", "dose_mM"]).reset_index(drop=True)


def _on_versus_off() -> None:
    """The comparison that makes the negative worth anything: right stressor beside wrong."""
    path = paths.data_dir() / "crosstalk" / "erox_2026-08_on_vs_off_target.tsv"
    if not path.exists():
        print("\n  (no on-vs-off-target table; run scripts/extract_crosstalk_workbook.py)")
        return
    table = pd.read_csv(path, sep="\t")
    print("\n" + "=" * 90)
    print("And the half that makes it a result: each sensor against its OWN stressor")
    print("=" * 90)
    print("\n  A sensor that stays quiet for the wrong stressor has said nothing until you")
    print("  know it shouts for the right one. Same three plates, same days.\n")
    for construct, group in table.groupby("construct", sort=False):
        group = group.sort_values("on_target_dose_mM")
        on = " ".join(f"{v:6.2f}" for v in group.on_target_fold)
        off = " ".join(f"{v:6.2f}" for v in group.off_target_fold_derived)
        print(f"    {construct:12s} {group.on_target_stressor.iloc[0]:>4s}, its own  {on}")
        print(f"    {'':12s} {group.off_target_stressor.iloc[0]:>4s}, the wrong {off}")
    print("\n  The negative folds are not repression. Above 1 mM the oxidative sensors'")
    print("  cultures are dying and a blank-subtracted signal goes below the blank.")


def main() -> int:
    table = scored(load())

    print("=" * 90)
    print("Sensor crosstalk: ER sensors vs H2O2, oxidative sensors vs DTT")
    print("=" * 90)
    for (construct, stressor), group in table.groupby(["construct", "stressor"]):
        print(f"\n  {construct} + {stressor}")
        print(f"    {'dose':>6} {'OD600':>7} {'viable':>7} {'worst':>7} {'raw fold':>9} "
              f"{'per-cell':>9} {'+-':>7}   reading")
        for r in group.itertuples():
            note = "" if r.culture_healthy else "  <- a culture is dying, not responding"
            print(f"    {r.dose_mM:6g} {r.od600:7.3f} {r.viability:6.0%} "
                  f"{r.viability_min:6.0%} {r.raw_fold:9.2f}x {r.per_cell_fold:8.2f}x "
                  f"{r.per_cell_fold_half_sd:6.3f}{note}")

    healthy = table[table.culture_healthy]
    span = healthy.per_cell_fold.max() - healthy.per_cell_fold.min()
    print("\n" + "=" * 90)
    print("The verdict")
    print("=" * 90)
    print(f"\n  Across every dose where EVERY plate's culture stays above "
          f"{VIABILITY_FLOOR:.0%} of control viability,")
    print(f"  the per-cell signal spans {healthy.per_cell_fold.min():.2f}x to "
          f"{healthy.per_cell_fold.max():.2f}x -- a range of {span:.2f}, over "
          f"{int(table.n_replicates.max())} biological replicates.")
    print("  NO CROSSTALK WORTH READING IN EITHER DIRECTION, and the wording is deliberate:")
    print(f"  at one replicate this range was 0.15 and now it is {span:.2f}. The off-target")
    print("  response is small, not absent, and the comparison below is what makes that")
    print("  distinction safe to draw.")
    _on_versus_off()

    worst = table.loc[table.raw_fold.idxmin()]
    print(f"\n  The trap: {worst.construct} at {worst.dose_mM:g} mM {worst.stressor} shows a")
    print(f"  raw fold of {worst.raw_fold:.2f}x, which reads as a large response. Its culture")
    print(f"  is at {worst.viability:.0%} of control. Per cell it is {worst.per_cell_fold:.2f}x,")
    print("  and what the raw number is measuring is the cells that are no longer there.")

    out = paths.outputs_dir() / "sensor_crosstalk.csv"
    table.to_csv(out, index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
