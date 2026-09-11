"""The auto feedback loop plates, read the way this repository says a reporter must be read.

    YSTWIN_GEN5_XPT="<Experiments dir>" python3 scripts/score_afl_circuit.py

Writes ``outputs/afl_circuit.csv``. Needs the Gen5 ``.xpt`` files; skips cleanly without them.

**Two experiments, not four replicates.** The lab ran the auto negative feedback loop twice,
against different comparators, and the notebook keeps them apart:

``AFL-vs-control``
    2026-07-05, the first biological replicate of "AFL vs control". Two circuits only --
    the loop and its open-loop control. ``No-Gal`` and ``No-LacI`` did not exist yet; the
    notebook builds them in the second half of July.
``AFL-debugging``
    2026-08-02, 08-03, 08-04. Four strains, after the team went back to work out what was
    wrong with the loop.

They are scored side by side here and never pooled. Two plates from different weeks with
different strain sets are not replicates of each other.

**Why this dataset and not another.** Every biosensor plate this project analyses is four
hours long, and `analysis/estimator_accuracy.py` established that the activity inversion's
error is governed by the smoothing window's duration in hours -- so on those plates the rule
lands at its floor and ~8% single-well error. These runs are 18-24 hours, where the rule picks
a 3.00-4.00 h window. The 2026-07-05 plate is the closest match in the project to that table:
145 points over 24 h is exactly the geometry `outputs/estimator_accuracy.csv` labels "the
regime the table was measured in", whose auto window is 4.0 h at 1.31% median single-well
error. These are the only plates in the project measured in the regime the estimator is good
at.

**What the experiment is.** An auto negative feedback loop: a galactose-inducible circuit
carrying its own repressor, so that raising the inducer raises the repressor which holds the
output down. Both plate maps put the inducer down the rows and the strains across the columns.

**Read the plates at the same time or do not compare them.** The four runs were stopped at
18.1, 20.1 and 24.1 hours. `RFU/OD` at the end of a run is not a property of the strain: it is
how much stable reporter has piled up since the culture stopped dividing, so a plate read for
six hours longer answers a different question. Every plate is therefore scored twice, at its
own length and truncated to the 109 timepoints all four share, and the two answers differ --
see ``docs/research/AUTO_FEEDBACK_LOOP.md``, which is mostly about that difference.

**What this script does not do.** It does not say whether the circuit is working. That is the
team's call against what they expected, and three of the four `.xpt` files are named
"debugging". It reports the numbers and their caveats.
"""
from __future__ import annotations

import pathlib
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.plate import gen5
from ystwin.reporter import (ReporterKinetics, default_activity_window_h,
                             promoter_activity_from_total)
from ystwin.readings import CorrectedOD, CorrectedRFU

#: Plate row -> galactose percentage. Identical in both notebook maps, and identical to the
#: stock ladder they describe: 20 uL of a 20% sugar stock into 180 uL of culture, with
#: glucose made up to 2% total, so the carbon load is constant down the plate and only the
#: inducer changes. Confirmed against the readings on every plate: time to reach blank-
#: corrected OD 1.0 lengthens down the rows, ~8 h in row A against ~11-12 h in row G.
GALACTOSE_PERCENT = dict(zip("ABCDEFG", (0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8)))

#: The medium blanks, row H of both maps. H4-H12 are empty, and the difference is what
#: confirms the map: on 2026-07-05 H1-H3 read 123 RFU against H4-H12's 55, and on the
#: August plates 172 against 59. That is medium against air.
BLANK_WELLS = ("H1", "H2", "H3")

#: Every plate truncated to its first 109 timepoints for the cross-plate comparison. 109 is
#: the length of the two shortest runs (2026-08-03 and 08-04), and it is a real truncation
#: rather than a resampling: all four files record the same 10-minute grid from 445 s to
#: 65,245 s, so the first 109 elapsed values are equal to the millisecond across all four.
#: Checked in ``tests/test_afl_circuit.py``.
MATCHED_TIMEPOINTS = 109

#: mCitrine is a stable YFP: loss is dilution, not degradation. The same assumption the
#: sensor pipeline makes, and the same one `analysis/estimator_accuracy.py` prices.
KINETICS = ReporterKinetics(k_deg=0.0)

#: Where the late window starts, as a fraction of the run. Matches the sensor pipeline.
LATE_FRACTION = 0.75

#: Below this blank-corrected optical density a well did not grow and its per-cell ratio is
#: a division by noise. Not a fitted threshold -- a reading aid, and stated as one.
MIN_BLANK_CORRECTED_OD = 0.02


@dataclass(frozen=True)
class Series:
    """One experimental series: plates that share a strain layout and a comparator set."""

    name: str
    plates: tuple[str, ...]
    """``.xpt`` file names, in date order."""

    strain_columns: dict[str, tuple[int, ...]]
    """Strain -> well columns."""

    strain_source: str
    """How :attr:`strain_columns` was established. ``'notebook'`` when the notebook's plate
    map carries a strain header row; ``'inferred-from-signal'`` when it does not and the
    assignment rests on the readings, which the class docstring of the series must then
    justify and ``tests/test_afl_circuit.py`` must check."""


SERIES = (
    Series(
        name="AFL-vs-control",
        plates=("20260705_AFL&control_1st.xpt",),
        # The notebook's map for 2026-07-05 gives the galactose ladder and the blanks and
        # STOPS: unlike the August maps it carries no strain header row, so which half of
        # the plate held which circuit is not written down. Three things fix it anyway.
        #   1. Only two circuits existed. The entry says "compare our autonegative feedback
        #      loop with the open-loop circuit" and prepares "5mL master stocks for each
        #      circuit (so 2 falcon tubes)". No-Gal and No-LacI are built weeks later.
        #   2. The plate splits 6/6 and nowhere else. Blank-corrected reporter at the first
        #      read is ~42-50 RFU in columns 1-6 and ~70-94 in columns 7-12; late RFU/OD is
        #      106-165 on the left and 194-442 on the right, at every dose.
        #   3. Those two ranges are the August AFL and Control ranges. AFL runs 84-146 there
        #      and Control 136-336, while No-Gal and No-LacI run 1,250-3,700 -- more than
        #      ten times anything on this plate, which is what "the strain is absent" looks
        #      like. The dim half is the closed loop because closing the loop is what makes
        #      it dim, which is the circuit's entire design.
        # If the assignment were swapped the only thing that changes is which block is
        # called the loop; every number below is computed per block and stands either way.
        strain_columns={"AFL": (1, 2, 3, 4, 5, 6), "Control": (7, 8, 9, 10, 11, 12)},
        strain_source="inferred-from-signal",
    ),
    Series(
        name="AFL-debugging",
        plates=("20260802_AFL-debugging_1st.xpt",
                "20260803_AFL-debugging_2nd.xpt",
                "20260804_AFL-debugging_3rd.xpt"),
        # Straight off the notebook maps, which head the column blocks
        # "AFL | Control | No-Gal | No-LacI".
        strain_columns={"AFL": (1, 2, 3), "Control": (4, 5, 6),
                        "No-Gal": (7, 8, 9), "No-LacI": (10, 11, 12)},
        strain_source="notebook",
    ),
)


def score_plate(path: pathlib.Path, series: Series, *,
                timepoints: int | None = None) -> pd.DataFrame:
    """One plate: naive RFU/OD and corrected activity, per strain and inducer dose.

    Args:
        path: The ``.xpt``.
        series: Which strain layout to read it with.
        timepoints: Keep only the first *n* reads, or all of them when ``None``. Used for
            the matched-length comparison; see :data:`MATCHED_TIMEPOINTS`.

    Raises:
        ValueError: if the file is shorter than *timepoints*, rather than silently
            returning a run of a different length under a label that says otherwise.
    """
    run = gen5.read_xpt(path)
    optical = run.channel("OD600:600").frame()
    reporter = run.channel("mCitrine:480,530").frame()
    hours = run.channel("OD600:600").elapsed_hours()
    if timepoints is not None:
        if len(hours) < timepoints:
            raise ValueError(
                f"{path.name} has {len(hours)} timepoints, fewer than the {timepoints} "
                "the matched comparison needs")
        # .iloc, not []: these frames are indexed by elapsed hours, a float index, on which
        # a plain slice is label-based and would silently keep the wrong rows.
        optical = optical.iloc[:timepoints]
        reporter = reporter.iloc[:timepoints]
        hours = hours[:timepoints]
    blank_od = optical[list(BLANK_WELLS)].to_numpy().mean(axis=1)
    blank_rfu = reporter[list(BLANK_WELLS)].to_numpy().mean(axis=1)
    late = slice(int(len(hours) * LATE_FRACTION), None)
    shared = {
        "experiment": series.name,
        "plate": path.stem,
        "truncation": "as-run" if timepoints is None else f"first-{timepoints}",
        "n_timepoints": len(hours),
        "run_hours": round(float(hours[-1]), 3),
        "activity_window_h": round(default_activity_window_h(hours), 3),
        "late_window_start_h": round(float(hours[late.start]), 3),
    }

    rows = []
    for strain, columns in series.strain_columns.items():
        for row, galactose in GALACTOSE_PERCENT.items():
            naive, corrected, grown, growth = [], [], 0, []
            for column in columns:
                well = f"{row}{column}"
                biomass = optical[well].to_numpy() - blank_od
                signal = reporter[well].to_numpy() - blank_rfu
                if biomass[late].mean() < MIN_BLANK_CORRECTED_OD:
                    continue
                grown += 1
                naive.append(float(np.mean((signal / biomass)[late])))
                activity = promoter_activity_from_total(
                    hours, CorrectedRFU(signal), CorrectedOD(np.maximum(biomass, 1e-6)), KINETICS)
                corrected.append(float(np.mean(activity[late])))
                growth.append(_late_growth_rate(hours, biomass, late))
            if not naive:
                continue
            rows.append({
                **shared, "strain": strain, "strain_source": series.strain_source,
                "well_columns": "-".join((str(columns[0]), str(columns[-1]))),
                "galactose_percent": galactose, "n_wells": grown,
                "naive_rfu_per_od": float(np.mean(naive)),
                "corrected_activity": float(np.mean(corrected)),
                "mu_late": float(np.mean(growth)),
            })
    return pd.DataFrame(rows)


def _late_growth_rate(hours: np.ndarray, biomass: np.ndarray, late: slice) -> float:
    """Specific growth rate over the late window, /h, as ``d(ln OD)/dt`` by least squares.

    **Reported because it is the confound this dataset cannot rule out by itself.** The
    2026-08-30 audit measured Spearman(corrected_activity, mu_late) at a median of 0.964
    across the strain x plate blocks, against Spearman(galactose, mu_late) of 0.25 to 1.00.
    A galactose ladder is also a growth ladder here, so "activity rises with dose" and
    "activity rises with whatever is still growing" make the same prediction, and nothing in
    this script could previously tell them apart -- mu was not even a column.

    Emitting it does not resolve the confound. It makes it checkable, which is the most this
    experiment supports: see `docs/research/AUTO_FEEDBACK_LOOP.md`.
    """
    window = np.asarray(hours)[late]
    values = np.asarray(biomass)[late]
    usable = values > 0
    if usable.sum() < 2:
        return float("nan")
    slope, _ = np.polyfit(window[usable], np.log(values[usable]), 1)
    return float(slope)


def _ladder(block: pd.DataFrame, column: str) -> pd.Series:
    return block.groupby("galactose_percent")[column].mean()


def _report(table: pd.DataFrame, truncation: str, heading: str) -> None:
    print(heading)
    for series in SERIES:
        for strain in series.strain_columns:
            block = table[(table.experiment == series.name) & (table.strain == strain)
                          & (table.truncation == truncation)]
            if block.empty:
                continue
            naive, corrected = _ladder(block, "naive_rfu_per_od"), _ladder(block, "corrected_activity")
            print(f"  {series.name:15s} {strain:8s} n={block.plate.nunique()}")
            print("    galactose      " + "".join(f"{g:8.1f}%" for g in naive.index))
            print("    naive RFU/OD   " + "".join(f"{v:9.0f}" for v in naive))
            print("    corrected      " + "".join(f"{v:9.1f}" for v in corrected))
            print(f"    -> naive {naive.iloc[-1] / naive.iloc[0]:.2f}x, "
                  f"corrected {corrected.iloc[-1] - corrected.iloc[0]:+.1f} over the ladder")
    print()


def main() -> int:
    directory = paths.gen5_xpt_dir()
    if directory is None:
        raise SystemExit(
            "the Gen5 .xpt files are not here. Set YSTWIN_GEN5_XPT to the Experiments "
            "directory. They are the raw instrument files and are not redistributable.")
    wanted = [name for series in SERIES for name in series.plates]
    missing = [n for n in wanted if not (directory / n).is_file()]
    if missing:
        raise SystemExit(f"missing from {directory}: {missing}")

    frames = []
    for series in SERIES:
        for name in series.plates:
            frames.append(score_plate(directory / name, series))
            frames.append(score_plate(directory / name, series,
                                      timepoints=MATCHED_TIMEPOINTS))
    table = pd.concat(frames, ignore_index=True)

    print("Auto feedback loop: galactose ladder, two experiments.\n")
    _report(table, "as-run",
            "AS RUN -- each plate to its own end, 18.1 h / 20.1 h / 24.1 h.\n"
            "The naive column's DIRECTION is not stable across these; that is the finding.\n")
    _report(table, f"first-{MATCHED_TIMEPOINTS}",
            f"MATCHED -- every plate truncated to the first {MATCHED_TIMEPOINTS} reads, "
            "18.12 h, an identical\ntime grid on all four. Now they can be compared.\n")

    out = paths.outputs_dir() / "afl_circuit.csv"
    table.to_csv(out, index=False)
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
