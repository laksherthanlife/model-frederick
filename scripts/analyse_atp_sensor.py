"""Reproduce the ATP-sensor analysis: growth phase versus response, on the real plates.

Reads the two Synergy exports in data/atp_sensor, compares every condition first at a fixed
timepoint and then at matched optical density, and prints the diauxic contrast that is left
once growth phase is taken out. See docs/ATP_SENSOR.md.

Usage:
    python scripts/analyse_atp_sensor.py
    python scripts/analyse_atp_sensor.py --figures            # into outputs/
    python scripts/analyse_atp_sensor.py --figures DIR

Why ``--figures`` exists
------------------------
Three ``outputs/atp_sensor_*.png`` panels were committed in August 2026 and no script
produced them; ``docs/AUDIT_2026_08.md`` recorded that as the pass's only new liability and
they were moved to ``archive/outputs/``. This flag is the writer they never had. It draws
**the panels of those renders whose numbers come from the plates**, from the same reading
path the printed report uses, and it deliberately draws nothing else:

* the ICL plate's four glucose traces, each with the hour it peaked -- archived panel
  ``out_of_sync`` A and ``reasoning_chain`` 2;
* the fold spread at one clock time against the fold spread at matched density --
  archived ``reasoning_chain`` 3;
* each condition's peak hour against the hour its culture stopped growing -- archived
  ``out_of_sync`` B.

Not drawn, and not because it was hard. The archived renders also carry a literature bar
(ICL1 >200-fold on non-fermentable carbon), a three-bar model-against-data comparison, and
an analytic OD-saturation curve. The first is a citation rather than a measurement on these
plates, and it is one ``docs/SOURCE_AUDIT.md`` section 3.1 records as **broken**: the
identifier the panel and ``generator/context.py`` both carry resolves to a bacteriophage
genome paper, and no Schuller/Entian 1996 paper could be found to replace it. The second is
an output of
``generator/context.py`` whose input has since moved -- ``docs/ATP_SENSOR.md`` states that
``_CARBON["galactose"]["carbon"]`` "has been put back to 0.45" and the table today reads
0.05 -- so the archived bars cannot be recomputed at today's constants and a figure that
redrew them at a value nobody can pin would be a guess. The third has a free parameter the
plates do not measure. ``docs/research/XPT_INVENTORY.md`` records the split panel by panel.

PNG only, no SVG, which is a departure from ``viz/figures.py`` and is deliberate twice
over: these are diagnostics rather than the tracked story set (``outputs/*.svg`` and
``outputs/*.png`` are both gitignored except ``fig[0-9][0-9]_*``), and the writer-map
cross-check in ``audit_reproducibility.py`` recognises ``csv|npz|png|xlsx`` in the
``docs/REPRODUCING.md`` table, so a documented SVG could not be cross-checked against the
source at all.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.matched_density import (
    fold_range,
    specific_fluorescence_at_density,
    time_to_density,
)
from ystwin.generator.synergy import read_export

# ``viz/figures.py`` forces the Agg backend on import, which is what makes this script safe
# to run over ssh, and it owns the page style and the caption block for every figure in the
# repository. ``_RC`` and ``_note`` are private to that module and are imported here on
# purpose: a second copy of the palette, the rcParams and the caption geometry is exactly
# how two figure sets drift apart, and this script is inside the same package.
from ystwin.viz.figures import _RC, _note, series_colours  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402  - must follow that import, for the backend
from ystwin.readings import CorrectedOD, CorrectedRFU

RATIOS = {1: "glu 10:0", 2: "gal 3:7", 3: "gal 6:4", 4: "gal 9:1"}
# Export filenames, not paths: which directory holds them is resolved at run time.
PLATES = [("ICL 2026-08-05", "20260805_ICL.xlsx", None, dict(
               A="high", B="mid-high", C="mid-low", D="low")),
          ("ACS/ICL sensor 1", "20260809_ACS_ICL.xlsx", set("ABC"), dict(
               A="2%", B="1.5%", C="0.5%")),
          ("ACS/ICL sensor 2", "20260809_ACS_ICL.xlsx", set("DEF"), dict(
               D="2%", E="1.5%", F="0.5%"))]

CLOCK_HOURS = (10.0, 12.0)
"""The two single timepoints this dataset has been compared at, and both are reported.

10 h is the hour the team's own bar chart used and the hour the archived panels marked;
12 h is the hour ``docs/ATP_SENSOR.md`` quotes for the plate-wide spread ("2.9-4.7x").
They are not the same number -- on the ICL plate the spread is 3.66x at 10 h and 4.66x at
12 h -- and which one a reader meets should not depend on which document they opened."""

MATCHED_DENSITIES = (0.5, 0.7, 0.9, 1.1, 1.3)
"""Optical densities to compare the conditions at, spanning what every well reached (0.5)
up to what only the fastest reached (1.3). The count of conditions that got there is
printed and drawn beside each one, because the fold range at 1.3 is over a subset."""

GROWTH_STOP_FRACTION = 0.9
"""Fraction of a well's final optical density that counts as "growth has stopped".

The cultures approach their plateau asymptotically, so the last reading is not the moment
growth ended; 90% of it is the conventional shorthand and it is stated here rather than
buried in the panel that uses it."""

PEAK_DOUBLINGS = 2.0
"""A peak in RFU/OD600 counts only if it arrives after the culture has doubled.

Every well's specific fluorescence starts high and falls through the first two to three
hours -- the inoculum's carried-over signal diluting into a growing population -- so the
unconstrained maximum of the series sits at t = 0 for eleven of the ICL plate's sixteen
conditions. Requiring the maximum to fall after the well has reached twice its starting
density is what separates a rise from that settling transient, and it reproduces exactly
the peak/no-peak split ``docs/ATP_SENSOR.md`` tabulates by hand."""

_OD_FLOOR = 0.02
"""Denominator floor for RFU/OD600, as in the printed report. A blanked OD at the start of
a run can be a few thousandths, and dividing by it turns instrument noise into a spike."""


def traces(path, keep_rows):
    """Every condition on one plate block, with its replicate wells kept apart.

    Args:
        path: The Synergy export to read.
        keep_rows: Plate rows to keep, or a falsy value for all of them. The 2026-08-09
            export carries two different builds, rows A-C and rows D-F.

    Returns:
        ``{(row, ratio): (times_h, od_frame, fluorescence_frame)}``, each frame one column
        per replicate well in column order.
    """
    layout, tables = read_export(path)
    od, fl = tables["Blank OD600:600"], tables["Blank mCitrine:480,530"]
    grouped = {}
    for well, sample in layout.items():
        if sample != "BLK" and (not keep_rows or well[0] in keep_rows):
            grouped.setdefault(sample, []).append(well)
    out = {}
    for wells in grouped.values():
        ordered = sorted(wells, key=lambda w: int(w[1:]))
        first = ordered[0]
        key = (first[0], RATIOS[(int(first[1:]) - 1) // 3 + 1])
        out[key] = (od.index.values, od[ordered], fl[ordered])
    return out


def conditions(path, keep_rows):
    """The same blocks as :func:`traces`, reduced to the mean over replicate wells."""
    return {key: (times, od.mean(axis=1).values, fl.mean(axis=1).values)
            for key, (times, od, fl) in traces(path, keep_rows).items()}


def specific(od, fluorescence):
    """RFU per OD600, with the floor the printed report uses."""
    return np.asarray(fluorescence, dtype=float) / np.maximum(np.asarray(od, dtype=float),
                                                              _OD_FLOOR)


def growth_stops_at(times, od):
    """Hour the culture first reached :data:`GROWTH_STOP_FRACTION` of its final density."""
    reached = np.maximum.accumulate(np.asarray(od, dtype=float))
    index = int(np.searchsorted(reached, GROWTH_STOP_FRACTION * reached[-1]))
    return float(times[index]) if index < len(times) else float("nan")


def peak_index(od, values):
    """Index of the RFU/OD600 peak, or ``None`` when the series never rises to one.

    ``None`` is the point of this function. Eleven of the ICL plate's sixteen conditions
    have their maximum in the first three hours, before the culture has doubled, and
    reporting that as a peak would put a settling transient on the same axis as a diauxic
    rise. See :data:`PEAK_DOUBLINGS`.
    """
    reached = np.maximum.accumulate(np.asarray(od, dtype=float))
    doubled = int(np.searchsorted(reached, PEAK_DOUBLINGS * reached[0]))
    index = int(np.argmax(values))
    return index if doubled < index < len(values) - 1 else None


def clock_folds(wells, hour):
    """Fold range across every condition read at one clock time."""
    values = []
    for times, od, fl in wells.values():
        index = int(np.abs(np.asarray(times) - hour).argmin())
        values.append(specific(od, fl)[index])
    return fold_range(values)


def matched_folds(wells):
    """``{density: (fold_range, n_reached, n_conditions)}`` across :data:`MATCHED_DENSITIES`."""
    out = {}
    for target in MATCHED_DENSITIES:
        got = [specific_fluorescence_at_density(t, CorrectedOD(o), CorrectedRFU(f), target)
               for t, o, f in wells.values()]
        out[target] = (fold_range(got), int(sum(np.isfinite(x) for x in got)), len(got))
    return out


# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------

_LINESTYLES = ("-", "--", "-.", ":")
"""Doubled encoding, so a greyscale print keeps every distinction the colour made."""

_SLUGS = {"ICL 2026-08-05": "icl_20260805",
          "ACS/ICL sensor 1": "acs_icl_sensor1_20260809",
          "ACS/ICL sensor 2": "acs_icl_sensor2_20260809"}
"""Filename stem per plate block. Written out rather than derived from the block name, so
that renaming a heading does not silently rename a file a document points at."""


def _peak_panel(axis, name, block, sugar, colours):
    """Panel A: the pure-glucose column, one trace per sugar concentration."""
    rows = sorted(sugar)
    for position, row in enumerate(rows):
        key = (row, "glu 10:0")
        if key not in block:
            continue
        times, od_frame, fl_frame = block[key]
        per_well = specific(od_frame.values, fl_frame.values)
        mean = per_well.mean(axis=1)
        axis.fill_between(times, per_well.min(axis=1), per_well.max(axis=1),
                          color=colours[position], alpha=0.18, linewidth=0)
        axis.plot(times, mean, color=colours[position], linewidth=1.6,
                  linestyle=_LINESTYLES[position % len(_LINESTYLES)],
                  label=f"{sugar[row]} sugar")
        index = peak_index(od_frame.mean(axis=1).values, mean)
        if index is not None:
            axis.plot(times[index], mean[index], marker="o", markersize=6,
                      color=colours[position], markeredgecolor="white")
            axis.annotate(f"h{times[index]:.1f}", (times[index], mean[index]),
                          textcoords="offset points", xytext=(4, 6), fontsize=8,
                          color=colours[position])
    for hour in CLOCK_HOURS:
        axis.axvline(hour, color="#555555", linewidth=0.9, linestyle=(0, (4, 3)))
        axis.annotate(f"read at {hour:.0f} h", (hour, 1.0), xycoords=("data", "axes fraction"),
                      textcoords="offset points", xytext=(3, -10), fontsize=7.5,
                      color="#555555", rotation=90, va="top")
    axis.set_title(f"{name}: each condition peaks at its own time")
    axis.set_xlabel("time (h)")
    axis.set_ylabel("RFU / OD600  (glucose-only column)")
    axis.legend(loc="upper left")


def _collapse_panel(axis, wells, clock, matched):
    """Panel B: one clock time against matched density."""
    densities = list(matched)
    folds = [matched[d][0] for d in densities]
    axis.plot(densities, folds, color="#0072B2", marker="o", linewidth=1.6,
              label="matched on OD600")
    for density in densities:
        fold, reached, total = matched[density]
        axis.annotate(f"{reached}/{total}", (density, fold), textcoords="offset points",
                      xytext=(0, -13), fontsize=7.5, ha="center", color="#333333")
    for position, hour in enumerate(sorted(clock)):
        axis.axhline(clock[hour], color="#D55E00", linewidth=1.1,
                     linestyle=_LINESTYLES[1 + position % 2])
        axis.annotate(f"read at {hour:.0f} h: {clock[hour]:.2f}x",
                      (0.02, clock[hour]), xycoords=("axes fraction", "data"),
                      textcoords="offset points", xytext=(0, 4), fontsize=8,
                      color="#D55E00")
    axis.set_title("Matched on density, the spread collapses")
    axis.set_xlabel("optical density compared at")
    axis.set_ylabel("fold range across all conditions")
    axis.set_ylim(bottom=1.0)
    axis.legend(loc="lower right")


def _timing_panel(axis, block, sugar):
    """Panel C: peak hour against the hour growth stopped, for conditions that peaked."""
    points = []
    for (row, ratio), (times, od_frame, fl_frame) in block.items():
        od = od_frame.mean(axis=1).values
        values = specific(od, fl_frame.mean(axis=1).values)
        index = peak_index(od, values)
        if index is not None:
            points.append((growth_stops_at(times, od), float(times[index]),
                           f"{sugar[row]} {ratio}", ratio == "glu 10:0"))
    if not points:
        axis.text(0.5, 0.5, "no condition on this block rises to a peak\nafter its culture "
                            "has doubled, so there is nothing to time",
                  ha="center", va="center", fontsize=9, color="#555555",
                  transform=axis.transAxes)
        axis.set_title("The peak arrives after growth stops")
        axis.set_xlabel(f"hour the culture reached {GROWTH_STOP_FRACTION:.0%} of its final OD600")
        axis.set_ylabel("hour of the RFU / OD600 peak")
        return
    span = [min(min(p[0], p[1]) for p in points) - 1.0,
            max(max(p[0], p[1]) for p in points) + 1.0]
    axis.plot(span, span, color="#999999", linewidth=1.0, linestyle=(0, (5, 4)),
              label="peak = growth arrest")
    for stop, peak, label, is_glucose in points:
        axis.plot(stop, peak, marker="o" if is_glucose else "^", markersize=7,
                  color="#009E73" if peak >= stop else "#CC79A7", markeredgecolor="white")
        axis.annotate(label, (stop, peak), textcoords="offset points", xytext=(6, -3),
                      fontsize=7.5, color="#333333")
    early = [label for stop, peak, label, _ in points if peak < stop]
    if early:
        # Pink is the only thing separating a peak that precedes growth arrest from one
        # that follows it, and a colour with nothing to read it against is decoration.
        axis.plot([], [], marker="s", linestyle="none", color="#CC79A7",
                  label=f"peak before growth arrest ({len(early)})")
    axis.set_title("The peak arrives after growth stops"
                   + (" -- mostly" if early else ""))
    axis.set_xlabel(f"hour the culture reached {GROWTH_STOP_FRACTION:.0%} of its final OD600")
    axis.set_ylabel("hour of the RFU / OD600 peak")
    axis.legend(loc="upper left")


def build_figure(name, filename, block, sugar, clock, matched):
    """One three-panel figure for one plate block."""
    labelled = {row: label for row, label in sugar.items() if any(k[0] == row for k in block)}
    with plt.rc_context(_RC):
        fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), layout="constrained")
        _peak_panel(axes[0], name, block, labelled, series_colours(max(len(labelled), 1)))
        _collapse_panel(axes[1], {k: (t, o.mean(axis=1).values, f.mean(axis=1).values)
                                  for k, (t, o, f) in block.items()}, clock, matched)
        _timing_panel(axes[2], block, labelled)
        below = [f"{d}: {matched[d][2] - matched[d][1]} of {matched[d][2]}"
                 for d in matched if matched[d][1] < matched[d][2]]
        _note(fig, "\n".join([
            f"Source: data/atp_sensor/{filename}, rows {''.join(sorted(labelled))}, three "
            f"wells per condition, read through generator/synergy.py's blank-subtracted "
            f"channels. Rebuild: python scripts/analyse_atp_sensor.py --figures",
            "Left: the shaded band is the range across the three replicate wells, not an "
            "interval on anything. A marked hour is a maximum reached after the culture "
            "doubled; a trace with no marker had its maximum in the settling transient and "
            "is reported as no peak rather than as a peak at t = 0.",
            "Middle: the fold range at a density only counts conditions that reached it -- "
            + ("the counts under each point say how many did; conditions short at "
               + "; ".join(below) + "." if below
               else "every condition reached every density drawn here."),
            "Right: nothing here says the rise is the promoter. Neither plate carried a "
            "promoterless strain, and flavin autofluorescence rises with mitochondrial mass "
            "at the same moment and scales with glucose the same way. See docs/ATP_SENSOR.md.",
        ]))
    return fig


def write_figures(out_dir, blocks):
    """Write one PNG per plate block and return the paths, newest content overwriting."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, filename, block, sugar, clock, matched in blocks:
        fig = build_figure(name, filename, block, sugar, clock, matched)
        png = out_dir / f"atp_sensor_{_SLUGS[name]}.png"
        with plt.rc_context(_RC):
            fig.savefig(png, format="png", metadata={"Software": "ystwin.viz"})
        plt.close(fig)
        written.append(png)
    return written


# ---------------------------------------------------------------------------
# the printed report
# ---------------------------------------------------------------------------


def rule(title):
    print(f"\n{title}\n{'-' * len(title)}")


def main(argv: list[str] | None = None) -> None:
    """Print the report, and write the figures when asked.

    ``argv`` defaults to no arguments rather than to ``sys.argv[1:]``: this function is
    called directly by ``tests/test_analyse_atp_sensor.py``, and a default that read the
    process's own arguments would make the script's behaviour depend on whatever launched
    it -- under pytest, on the test paths. The entry point below passes them explicitly.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures", nargs="?", const="", default=None, metavar="DIR",
                        help="also write one PNG per plate block; DIR defaults to outputs/")
    args = parser.parse_args(list(argv or ()))

    data = paths.resolve_or_exit(paths.atp_sensor_plates(), "the ATP-sensor plate exports",
                                 "YSTWIN_ATP_SENSOR")
    missing = sorted({f for _, f, _, _ in PLATES if not (data / f).exists()})
    if missing:
        raise SystemExit(f"{data} is missing {', '.join(missing)}; "
                         "point YSTWIN_ATP_SENSOR at the directory that holds them")
    summary = []
    blocks = []
    for name, filename, rows, sugar in PLATES:
        path = data / filename
        rule(name)
        block = traces(path, rows)
        wells = {key: (times, od.mean(axis=1).values, fl.mean(axis=1).values)
                 for key, (times, od, fl) in block.items()}
        clock = {hour: clock_folds(wells, hour) for hour in CLOCK_HOURS}
        matched = matched_folds(wells)
        for hour in sorted(clock):
            print(f"  read at t={hour:.0f}h        fold range {clock[hour]:.2f}")
        for target in MATCHED_DENSITIES:
            fold, reached, total = matched[target]
            print(f"  matched at OD {target}    fold range {fold:5.2f}"
                  f"   ({reached}/{total} conditions reached it)")
        hours = [time_to_density(t, CorrectedOD(o), 1.0) for t, o, _ in wells.values()]
        hours = [h for h in hours if np.isfinite(h)]
        print(f"  conditions are {max(hours) - min(hours):.1f} h apart at OD 1.0, "
              f"which is where the fixed-time spread comes from")

        top = sorted({r for r, _ in wells})[0]
        for ratio in ("glu 10:0", "gal 9:1"):
            t, o, f = wells[(top, ratio)]
            low = specific_fluorescence_at_density(t, CorrectedOD(o), CorrectedRFU(f), 0.5)
            high = specific_fluorescence_at_density(t, CorrectedOD(o), CorrectedRFU(f), 1.3)
            when = time_to_density(t, CorrectedOD(o), 1.3)
            print(f"  {ratio:9s} OD 0.5 -> 1.3: {low:7.0f} -> {high:7.0f} "
                  f"({high / low:.2f}x), reached h {when:.0f} of {t[-1]:.0f}")
        summary.append(dict(plate=name, at_clock=clock[CLOCK_HOURS[-1]],
                            at_density=matched[MATCHED_DENSITIES[0]][0]))
        blocks.append((name, filename, block, sugar, clock, matched))

    rule("what the plates measured")
    frame = pd.DataFrame(summary)
    print(frame.round(2).to_string(index=False))
    print("\nThe spread at one timepoint is growth phase. What is left at matched density is\n"
          "the diauxic shift, in glucose only -- which is what a CSRE reporter should read.\n"
          "See docs/ATP_SENSOR.md for what to run next.")

    if args.figures is not None:
        out_dir = pathlib.Path(args.figures) if args.figures else paths.outputs_dir()
        rule("figures")
        for written in write_figures(out_dir, blocks):
            print(f"  wrote {written}")


if __name__ == "__main__":
    main(sys.argv[1:])
