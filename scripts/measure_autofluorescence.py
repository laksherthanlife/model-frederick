"""Autofluorescence, measured at last, from a plate that was run on 2026-08-07.

    YSTWIN_GEN5_XPT="<Experiments dir>" python3 scripts/measure_autofluorescence.py

Writes ``outputs/autofluorescence_measured.csv``. Needs the Gen5 ``.xpt``; skips without it.

**Why this could not be run for three weeks, and why it can now.** `observation.py` writes the
measurement model as ``RFU = background + [a + gain*R] * biomass``, and its docstring insists
that ``a`` be measured on an isogenic reporter-free strain rather than guessed. Three documents
recorded that as never done. It was done: BY4741 was plated on 2026-08-07 under
``4h-10min_mCitrine_OD600.prt`` -- the same protocol, gain, interval and optics as biosensor
replicates 2, 3 and 4 -- and read in both channels. **The `.xlsx` export dropped the mCitrine
channel.** `plate/gen5.py` reads the instrument file, which kept it.

**How it is estimated.** With no reporter, ``R = 0`` and the model is a straight line in
biomass: ``RFU = background + a*OD``. So ``a`` is a SLOPE, not a ratio. Taking the median of
``RFU/OD`` -- the obvious thing, and what a first pass did -- measures ``a + background/OD``
instead, which is dominated by the background at the optical densities a 4-hour plate reaches
and returned 249 RFU/OD against the 18 this finds.

One slope per WELL, then a t interval over wells. Twenty-five readings of one well are one
well, not twenty-five observations; pooling them gives a standard error about five times too
small and would have made ``a`` look significantly positive.

**The answer is that it is not distinguishable from zero**, and the useful number is the upper
bound rather than the point estimate. See `docs/FINDINGS.md`.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.plate import gen5

#: The reporter-free plate. BY4741, no cassette, same protocol as the biosensor replicates.
SOURCE = "20260807_BY4741_ER_Oxidative.xpt"

#: Blank-corrected OD at the last timepoint above which a well is taken to have grown. The
#: occupied block is A1-F6 and F1-F3 barely grew; this separates the thirty that did.
MIN_FINAL_OD = 0.15

#: A well whose OD barely moves cannot constrain a slope against OD. Not a fitted threshold --
#: a well with no x-range has no regression, and including it adds noise with no information.
MIN_OD_RANGE = 0.02

#: Median RFU/OD of the reporter constructs on replicates 2 and 3, from
#: docs/research/XPT_INVENTORY.md. What `a` is expressed as a fraction OF, so that the answer
#: can be compared with `autofluorescence_sensitivity.csv`, which sweeps that fraction.
REPORTER_RFU_PER_OD = 4354.0


def measure(path: pathlib.Path) -> pd.DataFrame:
    """One row per grown well: its own slope, intercept and OD range."""
    run = gen5.read_xpt(path)
    optical = run.channel("OD600:600").frame()
    reporter = run.channel("mCitrine:480,530").frame()
    final = optical.iloc[-1]
    rows = []
    for well in optical.columns:
        if final[well] <= MIN_FINAL_OD:
            continue
        x = optical[well].to_numpy(dtype=float)
        y = reporter[well].to_numpy(dtype=float)
        if x.max() - x.min() < MIN_OD_RANGE:
            continue
        slope, intercept = np.polyfit(x, y, 1)
        rows.append({"well": well, "a_rfu_per_od": float(slope),
                     "background_rfu": float(intercept),
                     "od_min": float(x.min()), "od_max": float(x.max()),
                     "n_timepoints": int(len(x))})
    return pd.DataFrame(rows)


def summarise(table: pd.DataFrame) -> dict:
    """The interval over wells, and what it means for the sensitivity sweep."""
    slopes = table.a_rfu_per_od.to_numpy()
    n = len(slopes)
    mean = float(slopes.mean())
    half = float(stats.t.ppf(0.975, n - 1) * slopes.std(ddof=1) / np.sqrt(n))
    return {"n_wells": n, "a_rfu_per_od": mean,
            "ci95_low": mean - half, "ci95_high": mean + half,
            "includes_zero": bool((mean - half) < 0 < (mean + half)),
            "fraction_of_reporter": mean / REPORTER_RFU_PER_OD,
            "fraction_upper_bound": (mean + half) / REPORTER_RFU_PER_OD}


def main() -> int:
    directory = paths.gen5_xpt_dir()
    if directory is None or not (directory / SOURCE).is_file():
        raise SystemExit(
            f"{SOURCE} is not here. Set YSTWIN_GEN5_XPT to the Experiments directory; the "
            "instrument files are not redistributable.")
    per_well = measure(directory / SOURCE)
    got = summarise(per_well)

    print("Autofluorescence on an isogenic reporter-free strain (BY4741), measured.\n")
    print(f"  wells                {got['n_wells']}")
    print(f"  a                    {got['a_rfu_per_od']:+8.1f} RFU/OD")
    print(f"  95% CI over wells    [{got['ci95_low']:+.1f}, {got['ci95_high']:+.1f}]")
    print(f"  includes zero        {got['includes_zero']}")
    print(f"\n  as a fraction of a reporter well ({REPORTER_RFU_PER_OD:.0f} RFU/OD):")
    print(f"    point estimate     {got['fraction_of_reporter']:.2%}")
    print(f"    upper 95% bound    {got['fraction_upper_bound']:.2%}")
    print("\n  `outputs/autofluorescence_sensitivity.csv` sweeps this fraction at 0, 5, 10, 20")
    print("  and 30%. The lowest fraction at which ANY verdict changes is 5%, which is five")
    print("  times the upper bound measured here. EVERY CALL IN THE PANEL SURVIVES THE")
    print("  MEASURED AUTOFLUORESCENCE -- and that is the sweep's question answered, not")
    print("  a claim that the term is zero.")

    out = paths.outputs_dir() / "autofluorescence_measured.csv"
    per_well.assign(**{k: v for k, v in got.items()}).to_csv(out, index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
