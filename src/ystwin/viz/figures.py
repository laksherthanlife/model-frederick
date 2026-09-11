"""The five figures the findings actually rest on, and what each one refuses to draw.

A figure is an argument, and the arguments here are all about the gap between a number
and the evidence behind it. So the rendering rules are not decoration:

**No error bar without an interval.** ``analysis/uncertainty.py::fold_change`` refuses a
cluster bootstrap below three biological replicates, and the sensor-characterisation
table has two. Every corrected fold change in figure 02 is therefore a point estimate
with no interval, and the figure says so in words on its face. Drawing a whisker from the
spread of two plates would be inventing the very quantity the analysis declined to
compute.

**No point where the estimator broke.** At 2 and 4 mM H2O2 the recovered promoter
activity goes negative -- the culture is dying and the ODE inversion has nothing to
invert. A fold change from a negative denominator is not a small fold change; it is not a
fold change. Those doses are drawn as marked gaps.

**No number that is not in the file that was read.** The G1 panel plots the plates for
which a G1 table exists in ``outputs/`` and names, in a footnote, the recorded plates for
which one does not. The README quotes pass rates for those plates; this module does not,
because it has not seen them.

**Colour never carries a distinction on its own.** The palette is Okabe-Ito, which
survives every common form of colour vision deficiency, and every categorical encoding is
doubled by marker shape, line style or hatch so that a greyscale print still reads.

Two figures are generated rather than read, because their inputs are cheap and
deterministic simulations rather than measurements: the decoupling grid (figure 03) comes
from ``generator/design.py`` and the attribution power (figure 04) from
``analysis/power.py``. Both record the parameters they were run with on the figure
itself, since a simulated number with no visible settings is unreproducible.
"""

from __future__ import annotations

import pathlib
import textwrap
from collections.abc import Sequence

import matplotlib

# This module exists to write files. It never opens a window, and an interactive backend
# picked up from the environment is the difference between a figure build that works over
# ssh and one that dies importing Tk. Forced here rather than left to the caller so that
# tests, scripts and a wiki build all agree.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  - must follow the backend selection
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..analysis.power import (  # noqa: E402
    DEFAULT_WELL_CV,
    discrimination_power,
    ratio_discrimination_power,
)
from ..analysis.uncertainty import MIN_PLATES_FOR_INTERVAL, fold_change  # noqa: E402
from ..generator.design import collinearity, decoupling_grid, stressor_only_series  # noqa: E402
from ..generator.literature import parameters_from_literature  # noqa: E402
from ..generator.panel_calibration import fit_dose_response  # noqa: E402
from ..generator.unmixing import PANEL  # noqa: E402
from ..qpcr import REFERENCE_TARGET, TARGET_FOR_CONSTRUCT, gdna_share  # noqa: E402

__all__ = [
    "ARTEFACT_ZONES",
    "CORRECTABLE_GDNA_CEILING",
    "DESIGN_CONSTRUCT",
    "DESIGN_DOSES",
    "HEADLINE_NULLS",
    "NUTRIENT_LADDERS",
    "PALETTE",
    "PASS_MARGIN_CYCLES",
    "PEER_BINS",
    "POWER_FOLDS",
    "POWER_REPLICATES",
    "RATIO_PAIR",
    "ROBUSTNESS_SWEEPS",
    "VERDICT_ABOVE",
    "VERDICT_BELOW",
    "VERDICT_INCONCLUSIVE",
    "VERDICT_NOT_ESTIMABLE",
    "VERDICT_NO_INTERVAL",
    "attribution_power",
    "attribution_power_table",
    "decoupling_designs",
    "decoupling_grid_scatter",
    "dose_response_correction",
    "dose_response_table",
    "g1_pass_rates",
    "g1_summary",
    "g4_contamination",
    "g4_contamination_table",
    "heldout_skill_table",
    "heldout_skill_vs_saturation",
    "identifiable_ec50",
    "interval_robustness",
    "interval_robustness_table",
    "null_comparison_table",
    "read_table",
    "redundancy_recovery_table",
    "save_figure",
    "series_colours",
    "transfer_honesty",
]


# ---------------------------------------------------------------------------
# palette and page style
# ---------------------------------------------------------------------------

PALETTE = (
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
    "#000000",  # black
)
"""Okabe-Ito, in a fixed order, and eight is the whole of it.

Chosen because it is the palette with published separation under deuteranopia,
protanopia and tritanopia rather than under one of them, and because red and green never
appear as the only difference between two series here. It is deliberately not extensible:
:func:`series_colours` refuses a ninth series rather than cycling, because a cycled
palette silently gives two categories the same colour, which is worse than no figure.
"""

_HATCHES = ("", "//", "xx", "\\\\", "..")
"""Secondary encoding, so a greyscale print keeps every distinction the colour made."""

_RC = {
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 10.0,
    "axes.titlesize": 11.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    # matplotlib derives SVG element ids from a hash salted with a per-process value
    # unless one is pinned, which makes two byte-identical figures differ on disk.
    "svg.hashsalt": "ystwin",
    "svg.fonttype": "path",
}


def series_colours(n: int) -> tuple[str, ...]:
    """The first ``n`` palette entries, or a refusal.

    Args:
        n: Number of categories that need to be told apart.

    Returns:
        ``n`` hex colours from :data:`PALETTE`, in its fixed order.

    Raises:
        ValueError: if ``n`` is not positive, or exceeds the palette. Cycling would
            give two categories one colour; splitting the figure, or grouping
            categories before plotting, is the fix the message asks for.
    """
    if n <= 0:
        raise ValueError(f"a figure needs at least one series; got {n}")
    if n > len(PALETTE):
        raise ValueError(
            f"{n} series requested but the colour-blind-safe palette has "
            f"{len(PALETTE)}; group the categories or split the figure rather than "
            f"cycling colours, which would give two categories the same colour"
        )
    return PALETTE[:n]


def _hatch(index: int) -> str:
    return _HATCHES[index % len(_HATCHES)]


# ---------------------------------------------------------------------------
# reading and writing
# ---------------------------------------------------------------------------


def read_table(path: pathlib.Path, what: str, produced_by: str) -> pd.DataFrame:
    """A committed table, or an error naming the script that writes it.

    A figure built from a table that is not there is a figure with no series in it,
    which renders perfectly and says nothing. That failure is silent at exactly the
    moment it matters -- a regenerated wiki -- so absence is an error here and not a
    skip.

    Args:
        path: Location of the CSV.
        what: Human name of the table, for the message.
        produced_by: Command that writes it.

    Returns:
        The table, guaranteed non-empty.

    Raises:
        FileNotFoundError: if the file is absent, naming ``produced_by``.
        ValueError: if the file is empty or has no data rows.
    """
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{what} not found at {path}. It is a committed output; regenerate it with "
            f"`{produced_by}`."
        )
    try:
        frame = pd.read_csv(path)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(
            f"{what} at {path} is empty. An empty gate report reads like a plate that "
            f"passed nothing rather than like a table that was never written; "
            f"regenerate it with `{produced_by}`."
        ) from exc
    if frame.empty:
        raise ValueError(
            f"{what} at {path} has a header and no rows; regenerate it with "
            f"`{produced_by}`."
        )
    return frame


def save_figure(fig, out_dir: pathlib.Path, stem: str) -> tuple[pathlib.Path, pathlib.Path]:
    """Write one figure as both SVG and PNG under a deterministic stem.

    SVG because the report embeds it inline and a wiki reader may zoom; PNG because a
    slide deck and a poster cannot take an SVG. The filenames are a pure function of
    ``stem``, so a rebuild overwrites rather than accumulating a second copy under a
    timestamp.

    The SVG's ``Date`` metadata is suppressed: matplotlib stamps the current time into
    it otherwise, and a file whose bytes change every build cannot be reviewed in a diff.

    Args:
        fig: The figure to write.
        out_dir: Directory to write into. Created if absent.
        stem: Filename without extension, e.g. ``01_g1_pass_rates``.

    Returns:
        ``(svg_path, png_path)``.
    """
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    svg = out_dir / f"{stem}.svg"
    png = out_dir / f"{stem}.png"
    # The rc context has to cover the write, not only the drawing: svg.hashsalt is read
    # by the SVG writer, and without it here every rebuild produces different element ids
    # and therefore different bytes for an identical figure.
    with plt.rc_context(_RC):
        fig.savefig(svg, format="svg", metadata={"Date": None})
        fig.savefig(png, format="png", metadata={"Software": "ystwin.viz"})
    plt.close(fig)
    return svg, png


def _figure(*args, **kwargs):
    """A figure under this module's rcParams, with constrained layout."""
    kwargs.setdefault("layout", "constrained")
    with plt.rc_context(_RC):
        return plt.subplots(*args, **kwargs)


_NOTE_FONTSIZE = 7.8
_NOTE_LINE_IN = 0.148
_NOTE_CHAR_IN = 0.0565
"""Width of one character of note text, in inches, at :data:`_NOTE_FONTSIZE`.

Measured rather than derived, and used only to wrap the note to the figure. Matplotlib
will not wrap for us without also refusing to report how tall the result is, and a note
that runs off the right edge loses exactly the caveat it was written to carry.
"""


def _note(fig, text: str) -> None:
    """A provenance or refusal block along the bottom of a figure.

    The block is not decoration and must never be cropped, so it is given room in the
    layout rather than written past the edge of the canvas and left to ``bbox_inches``
    to rescue. Constrained layout is told to keep the bottom strip clear, and the text
    is placed inside it.

    Args:
        fig: Figure to annotate. Must have been created with constrained layout.
        text: The note. Newlines are respected; nothing is re-wrapped, so keep lines
            short enough for the figure's width.
    """
    width = max(40, int(fig.get_figwidth() / _NOTE_CHAR_IN))
    paragraphs = []
    for line in text.split("\n"):
        indent = line[: len(line) - len(line.lstrip(" "))]
        paragraphs.append(textwrap.fill(
            line.strip(), width=width, initial_indent=indent,
            subsequent_indent=indent + "  ") if line.strip() else "")
    text = "\n".join(paragraphs)
    n_lines = text.count("\n") + 1
    reserved = min((_NOTE_LINE_IN * n_lines + 0.12) / fig.get_figheight(), 0.45)
    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(rect=(0.004, reserved, 0.992, 0.996 - reserved))
    fig.text(0.006, reserved * 0.90, text, fontsize=_NOTE_FONTSIZE, va="top",
             ha="left", color="#333333", linespacing=1.35)


def _readable_on(colour: str) -> str:
    """Black or white, whichever the eye can read on ``colour``.

    A fixed white label disappears on the pale half of any palette, and the palette here
    deliberately spans light yellow to black. Relative luminance decides it.
    """
    rgb = [int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    luminance = 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
    return "#000000" if luminance > 0.28 else "#FFFFFF"


# ---------------------------------------------------------------------------
# 01 - G1 optical quality, before and after the NewProtocol
# ---------------------------------------------------------------------------


def g1_summary(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per-plate outcome counts, partitioned by exact failure signature.

    Counting failure *reasons* would double-count: a well that fails both
    ``linear_range`` and ``sustained_decline`` appears under each, and the reasons then
    sum past the number of wells, so a stacked bar of them is taller than the plate. The
    signature -- the full comma-joined set of criteria a well failed -- partitions the
    wells exactly, which is what a stacked bar requires to be honest.

    Args:
        tables: Plate label -> the frame read from one ``outputs/g1_*.csv``, which must
            carry ``well``, ``passed`` and ``failures``.

    Returns:
        Long frame with ``plate``, ``outcome``, ``n_wells_with_outcome``, ``n_wells``,
        ``n_pass``, ``pass_rate``, ``n_criteria``. ``outcome`` is ``"passed"`` or a
        failure signature.

    Raises:
        ValueError: if ``tables`` is empty, or a frame lacks the gate's columns.
    """
    if not tables:
        raise ValueError(
            "no G1 tables to summarise; pass at least one outputs/g1_*.csv frame"
        )
    rows = []
    for label, frame in tables.items():
        missing = {"passed", "failures"} - set(frame.columns)
        if missing:
            raise ValueError(
                f"G1 table {label!r} is missing {sorted(missing)}; this is not a "
                f"gates/g1_optical.py assess_plate export"
            )
        if frame.empty:
            raise ValueError(f"G1 table {label!r} has no wells in it")
        passed = frame.passed.astype(bool)
        n_wells = int(len(frame))
        n_pass = int(passed.sum())
        counts = {"passed": n_pass}
        signatures = frame.loc[~passed, "failures"].fillna("").astype(str)
        for signature, n in signatures.value_counts().items():
            counts[signature] = int(n)
        for outcome, n in counts.items():
            rows.append({
                "plate": label,
                "outcome": outcome,
                "n_wells_with_outcome": int(n),
                "n_wells": n_wells,
                "n_pass": n_pass,
                "pass_rate": n_pass / n_wells,
                "n_criteria": 0 if outcome == "passed" else outcome.count(",") + 1,
            })
    return pd.DataFrame(rows)


def _outcome_order(summary: pd.DataFrame) -> list[str]:
    """Passed first, then failure signatures by how many criteria they broke."""
    failures = summary[summary.outcome != "passed"]
    ordered = sorted(set(failures.outcome), key=lambda s: (s.count(","), s))
    return (["passed"] if (summary.outcome == "passed").any() else []) + ordered


def g1_pass_rates(
    summary: pd.DataFrame,
    out_dir: pathlib.Path,
    plates_without_a_table: Sequence[str] = (),
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure 01: what G1 rejected on each plate, and why.

    Args:
        summary: Output of :func:`g1_summary`.
        out_dir: Directory for the SVG and PNG.
        plates_without_a_table: Recorded plates with no G1 export in ``outputs/``.
            Named in a footnote rather than left out, because the interesting
            comparison is with plates whose numbers this figure has not seen.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``summary`` is empty or carries more outcomes than the palette.
    """
    if summary.empty:
        raise ValueError("nothing to plot: the G1 summary has no rows")
    plates = list(dict.fromkeys(summary.plate))
    outcomes = _outcome_order(summary)
    colours = series_colours(len(outcomes))

    fig, ax = _figure(figsize=(9.6, 5.6))
    x = np.arange(len(plates), dtype=float)
    bottom = np.zeros(len(plates))
    handles = []
    for index, (outcome, colour) in enumerate(zip(outcomes, colours)):
        heights = np.array([
            float(summary[(summary.plate == p) & (summary.outcome == outcome)]
                  .n_wells_with_outcome.sum())
            for p in plates
        ])
        bars = ax.bar(
            x, heights, bottom=bottom, width=0.56, color=colour,
            edgecolor="white", linewidth=0.8,
            hatch=_hatch(0 if outcome == "passed" else outcome.count(",") + 1),
        )
        handles.append(bars)
        for xi, h, b in zip(x, heights, bottom):
            if h >= 3:
                ax.text(xi, b + h / 2, f"{int(h)}", ha="center", va="center",
                        fontsize=8.5, color=_readable_on(colour), fontweight="bold")
        bottom = bottom + heights

    for xi, plate in zip(x, plates):
        block = summary[summary.plate == plate].iloc[0]
        ax.text(xi, block.n_wells + max(bottom) * 0.02,
                f"{int(block.n_pass)}/{int(block.n_wells)} pass  ({block.pass_rate:.0%})",
                ha="center", va="bottom", fontsize=9.5, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(plates, fontsize=9)
    ax.set_ylabel("cultures on the plate (wells)")
    ax.set_ylim(0, max(bottom) * 1.18)
    ax.set_title("G1 optical quality: every well, and the exact criteria it failed")
    ax.legend(
        handles,
        [("passed all four criteria" if o == "passed" else o.replace(",", " + "))
         for o in outcomes],
        title="outcome (signatures partition the wells; no well is counted twice)",
        loc="upper left", bbox_to_anchor=(1.01, 1.0), title_fontsize=8.5,
    )
    ax.margins(x=0.12)

    note = (
        "G1 asks whether an optical channel is quantitative for a well: linear range, "
        "blank separation, sustained decline, dynamic range.\n"
        "Read from outputs/g1_*.csv. Gate: od_linear_max=1.0, max_decline_fraction=0.15 "
        "(scripts/run_gates.py)."
    )
    if plates_without_a_table:
        note += (
            "\nNo G1 table in outputs/ for the recorded plates "
            + ", ".join(plates_without_a_table) + " -- the NewProtocol replicates."
            "\ndocs/FINDINGS.md quotes their pass rates; this figure does not, because it "
            "has not read them. To produce them:"
            "\n    YSTWIN_IGEM_RESULTS=<plate directory> python3 scripts/run_gates.py"
        )
    _note(fig, note)
    return save_figure(fig, out_dir, "01_g1_pass_rates")


# ---------------------------------------------------------------------------
# 02 - naive versus dilution-corrected dose response
# ---------------------------------------------------------------------------

ARTEFACT_ZONES = {
    "H2O2": (0.08, 0.25),
    "DTT": (1.7, 6.0),
}
"""Dose windows where the naive readout and the corrected one disagree about the biology.

Not derived from the plotted numbers -- they are an annotation, and deriving them from
the same data they annotate would make the highlight circular. They record the two
claims in ``docs/FINDINGS.md`` that this figure exists to check: that the oxidative sensors do not
detect 0.1-0.2 mM H2O2 but do register the growth slowdown it causes, and that the ER
sensors' largest apparent inductions at 2-5 mM DTT are dilution rather than induction.
A stressor absent from this map is simply not shaded.
"""

_STRESSOR_UNITS = {"H2O2": "mM H2O2", "DTT": "mM DTT"}

VERDICT_ABOVE = "resolved above 1.0"
VERDICT_BELOW = "resolved below 1.0"
VERDICT_INCONCLUSIVE = "inconclusive"
VERDICT_NO_INTERVAL = "no interval"
VERDICT_NOT_ESTIMABLE = "not estimable"
"""What a single dose's corrected fold change is allowed to be read as.

The distinction the figures draw, and the reason the column exists rather than being
recomputed at draw time: a fold of 1.4 whose interval runs [1.2, 1.8] and a fold of 1.4
whose interval runs [0.7, 2.9] are the same point estimate and different results, and a
plot that draws them identically has thrown away the finding. ``VERDICT_NOT_ESTIMABLE``
is stronger than ``VERDICT_INCONCLUSIVE``: it means the recovered promoter activity went
negative and there is no ratio to take at all, which is not a wide interval but the
absence of one.
"""


def _fold_verdict(row) -> str:
    """Classify one row of :func:`dose_response_table`.

    ``n_plates_usable_corrected`` and not ``n_plates`` decides estimability, for the
    reason :func:`_usable_plate_count` exists: ``fold_change`` drops the plates a
    geometric mean cannot take, so a point can rest on fewer plates than it claims.
    """
    if not np.isfinite(row.corrected_fold) or row.n_plates_usable_corrected < 2:
        return VERDICT_NOT_ESTIMABLE
    if not row.interval_estimable:
        return VERDICT_NO_INTERVAL
    if row.interval_low > 1.0:
        return VERDICT_ABOVE
    if row.interval_high < 1.0:
        return VERDICT_BELOW
    return VERDICT_INCONCLUSIVE


def _usable_plate_count(
    readings: pd.DataFrame, construct: str, dose: float, value: str,
    control_dose: float = 0.0,
) -> int:
    """Plates whose own dosed/control ratio is finite and positive.

    ``uncertainty.fold_change`` combines per-plate folds and silently drops the ones a
    geometric mean cannot take -- which is correct, and which means a reported point can
    rest on fewer plates than its ``n_plates`` says. That distinction decides whether a
    dose is drawable at all, and the public result object does not carry it, so it is
    recomputed here under the same rule.
    """
    subset = readings[readings.construct == construct]
    usable = 0
    for plate, group in subset.groupby("plate", sort=True):
        dosed = group[np.isclose(group.dose_mM, dose)][value].to_numpy(dtype=float)
        control = group[np.isclose(group.dose_mM, control_dose)][value].to_numpy(dtype=float)
        if dosed.size == 0 or control.size == 0:
            continue
        denominator = float(np.mean(control))
        if not np.isfinite(denominator) or abs(denominator) < 1e-12:
            continue
        ratio = float(np.mean(dosed)) / denominator
        if np.isfinite(ratio) and ratio > 0:
            usable += 1
    return usable


def dose_response_table(readings: pd.DataFrame, control_dose: float = 0.0) -> pd.DataFrame:
    """Naive and dilution-corrected fold change at every dose of every construct.

    Both quantities come from ``analysis/uncertainty.py::fold_change``, which computes
    the fold within each plate and combines afterwards, so a plate-level shift that
    multiplies dosed and control wells alike cancels instead of moving the estimate.

    Args:
        readings: The ``outputs/sensor_characterisation.csv`` frame, or one shaped like
            it: ``plate``, ``construct``, ``stressor``, ``dose_mM``, ``naive_late``,
            ``activity_late``.
        control_dose: The reference condition each fold is taken against.

    Returns:
        One row per construct and dose above the control, carrying both folds, the
        number of biological replicates behind them, the number of plates on which each
        fold was actually computable, whether an interval was estimable, the reason it
        was not, and a ``verdict`` from :data:`VERDICT_ABOVE` and its neighbours saying
        what the interval permits the point to be read as.

    Raises:
        ValueError: if a required column is absent or the frame carries no dosed rows.
    """
    required = {"plate", "construct", "stressor", "dose_mM", "naive_late", "activity_late"}
    missing = required - set(readings.columns)
    if missing:
        raise ValueError(
            f"sensor readings are missing {sorted(missing)}; expected the columns "
            f"scripts/run_sensor_characterisation.py writes"
        )
    if readings.empty:
        raise ValueError("sensor readings frame has no rows")

    rows = []
    for construct in sorted(set(readings.construct)):
        subset = readings[readings.construct == construct]
        stressor = str(subset.stressor.iloc[0])
        doses = sorted(float(d) for d in set(subset.dose_mM) if not np.isclose(d, control_dose))
        if not doses:
            raise ValueError(
                f"construct {construct!r} has only its control dose; nothing to ratio"
            )
        for dose in doses:
            naive = fold_change(readings, construct, dose, value="naive_late",
                                control_dose=control_dose)
            corrected = fold_change(readings, construct, dose, value="activity_late",
                                    control_dose=control_dose)
            rows.append({
                "construct": construct,
                "stressor": stressor,
                "dose_mM": dose,
                "naive_fold": naive.point,
                "corrected_fold": corrected.point,
                "n_plates": corrected.n_plates,
                "n_plates_usable_naive": _usable_plate_count(
                    readings, construct, dose, "naive_late", control_dose),
                "n_plates_usable_corrected": _usable_plate_count(
                    readings, construct, dose, "activity_late", control_dose),
                "interval_estimable": bool(corrected.estimable),
                "interval_low": corrected.low,
                "interval_high": corrected.high,
                "interval_method": corrected.method,
            })
    table = pd.DataFrame(rows)
    table["verdict"] = [_fold_verdict(row) for row in table.itertuples()]
    return table


def dose_response_correction(
    folds: pd.DataFrame, out_dir: pathlib.Path,
    stem: str = "02_dose_response_correction",
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure 02: what the dilution correction does to every dose of every construct.

    The headline figure, and it has to carry three things at once. The open dashed series
    is the raw per-cell readout, which says "induction". The filled solid series is the
    same data after the reporter ODE is inverted for growth dilution, which says much
    less. The whisker is the 95% interval, which says whether we can tell at all -- and a
    filled marker means the interval clears 1.0 while a hollow one means it does not.

    The interval is drawn only where ``fold_change`` estimated one. Below
    ``MIN_PLATES_FOR_INTERVAL`` biological replicates it refuses a cluster bootstrap and
    there is nothing to draw; a whisker taken from two plates would describe which of two
    plates was drawn rather than the population of plates. The note says which regime the
    figure is in, because "no error bars" and "error bars that cross 1.0" are different
    results and a reader must not have to guess which one they are looking at.

    Args:
        folds: Output of :func:`dose_response_table`.
        out_dir: Directory for the SVG and PNG.
        stem: Filename without extension. Defaults to the report's own numbering; the
            story set in ``outputs/`` passes its own.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``folds`` is empty or lacks the columns the table builder writes.
    """
    required = {"construct", "stressor", "dose_mM", "naive_fold", "corrected_fold",
                "n_plates", "n_plates_usable_corrected", "n_plates_usable_naive"}
    missing = required - set(folds.columns)
    if missing:
        raise ValueError(f"fold table is missing {sorted(missing)}")
    if folds.empty:
        raise ValueError("nothing to plot: the fold-change table has no rows")
    folds = folds.copy()
    if "verdict" not in folds.columns:
        folds["verdict"] = [_fold_verdict(row) for row in folds.itertuples()]
    has_intervals = {"interval_estimable", "interval_low", "interval_high"} <= set(folds.columns)

    constructs = sorted(set(folds.construct),
                        key=lambda c: (folds[folds.construct == c].stressor.iloc[0], c))
    naive_colour, corrected_colour = PALETTE[1], PALETTE[0]
    n_cols = 2
    n_rows = int(np.ceil(len(constructs) / n_cols))
    fig, axes = _figure(n_rows, n_cols, figsize=(10.4, 4.4 * n_rows), squeeze=False)
    flat = [ax for row in axes for ax in row]

    # The interval bounds set the extent of the drawing, not the point estimates: an
    # axis scaled to the points alone clips the very whisker that carries the result.
    spans = [folds.naive_fold.to_numpy(dtype=float),
             folds.corrected_fold.to_numpy(dtype=float)]
    if has_intervals:
        drawn = folds.interval_estimable.astype(bool) & (folds.n_plates_usable_corrected >= 2)
        spans += [folds.loc[drawn, "interval_low"].to_numpy(dtype=float),
                  folds.loc[drawn, "interval_high"].to_numpy(dtype=float)]
    finite = np.concatenate(spans)
    finite = finite[np.isfinite(finite)]
    top = float(np.nanmax(finite)) * 1.16 if finite.size else 2.0
    bottom = min(0.85, float(np.nanmin(finite)) * 0.92) if finite.size else 0.8

    for ax, construct in zip(flat, constructs):
        panel = folds[folds.construct == construct].sort_values("dose_mM")
        stressor = str(panel.stressor.iloc[0])
        zone = ARTEFACT_ZONES.get(stressor)
        if zone is not None:
            ax.axvspan(zone[0], zone[1], color="#000000", alpha=0.055, lw=0, zorder=0)
        ax.axhline(1.0, color="#555555", lw=1.0, ls=(0, (4, 3)), zorder=1)

        naive = panel[np.isfinite(panel.naive_fold) & (panel.n_plates_usable_naive >= 2)]
        ax.plot(naive.dose_mM, naive.naive_fold, marker="o", ls=(0, (5, 2)),
                color=naive_colour, markerfacecolor="none", markeredgecolor=naive_colour,
                markersize=7, lw=1.8, zorder=3)

        drawable = panel[np.isfinite(panel.corrected_fold)
                         & (panel.n_plates_usable_corrected >= 2)]
        ax.plot(drawable.dose_mM, drawable.corrected_fold, ls="-", color=corrected_colour,
                lw=1.8, zorder=3)
        for row in drawable.itertuples():
            resolved = row.verdict in (VERDICT_ABOVE, VERDICT_BELOW)
            if has_intervals and bool(row.interval_estimable):
                ax.errorbar(
                    row.dose_mM, row.corrected_fold,
                    yerr=[[row.corrected_fold - row.interval_low],
                          [row.interval_high - row.corrected_fold]],
                    ecolor=corrected_colour, elinewidth=1.6 if resolved else 1.0,
                    capsize=4.0, capthick=1.6 if resolved else 1.0, ls="none", zorder=4,
                )
            ax.plot([row.dose_mM], [row.corrected_fold], marker="s", ls="none",
                    color=corrected_colour, markersize=7.5,
                    markerfacecolor=corrected_colour if resolved else "#FFFFFF",
                    markeredgecolor=corrected_colour, markeredgewidth=1.6, zorder=5)

        refused = panel[panel.verdict == VERDICT_NOT_ESTIMABLE]
        for dose in refused.dose_mM:
            ax.axvline(dose, color="#777777", lw=8, alpha=0.18, zorder=0)
            ax.text(dose, bottom + (top - bottom) * 0.03, "not estimable",
                    rotation=90, ha="center", va="bottom", fontsize=7.5,
                    color="#444444")

        ax.set_xscale("log")
        ax.set_xticks(list(panel.dose_mM))
        ax.set_xticklabels([f"{d:g}" for d in panel.dose_mM])
        ax.minorticks_off()
        ax.set_xlim(float(panel.dose_mM.min()) * 0.7, float(panel.dose_mM.max()) * 1.5)
        ax.set_ylim(bottom, top)
        ax.set_xlabel(f"dose ({_STRESSOR_UNITS.get(stressor, 'mM ' + stressor)}, log scale)")
        ax.set_ylabel("fold change vs 0 mM control")
        ax.set_title(f"{construct}  ({stressor})")

    for ax in flat[len(constructs):]:
        ax.set_visible(False)

    resolved = folds[folds.verdict.isin((VERDICT_ABOVE, VERDICT_BELOW))]
    estimable = folds[folds.verdict.isin(
        (VERDICT_ABOVE, VERDICT_BELOW, VERDICT_INCONCLUSIVE))]
    n_plates = int(folds.n_plates.max())

    handles = [
        Line2D([], [], color=naive_colour, ls=(0, (5, 2)), marker="o", lw=1.8,
               markerfacecolor="none", markersize=7, label="naive RFU/OD"),
        Line2D([], [], color=corrected_colour, ls="-", marker="s", lw=1.8, markersize=7.5,
               markerfacecolor=corrected_colour, markeredgewidth=1.6,
               label="dilution-corrected, interval excludes 1.0"),
        Line2D([], [], color=corrected_colour, ls="-", marker="s", lw=1.8, markersize=7.5,
               markerfacecolor="#FFFFFF", markeredgewidth=1.6,
               label="dilution-corrected, INCONCLUSIVE (interval spans 1.0)"),
        Line2D([], [], color="#555555", ls=(0, (4, 3)), lw=1.0, label="no change (1.0)"),
    ]
    # The key goes inside the first panel rather than in a figure-level strip. A figure
    # legend has to be given vertical room by the layout engine, and `_note` has already
    # taken that room for the caption; the two then land on top of each other and on the
    # title. The band above the naive curve in the first panel is empty by construction --
    # the y limit is set by the widest interval anywhere, which is in another panel.
    flat[0].legend(handles=handles, loc="upper left", fontsize=8.0, ncols=1,
                   handlelength=2.4, labelspacing=0.35, borderaxespad=0.4)
    fig.suptitle(
        "Naive says induction; the correction says much less; the interval says whether "
        "we can tell", fontsize=12.5, fontweight="bold",
    )

    if len(estimable):
        calls = ", ".join(
            f"{row.construct} at {row.dose_mM:g} mM ({row.corrected_fold:.2f}, "
            f"[{row.interval_low:.2f}, {row.interval_high:.2f}])"
            for row in resolved.sort_values(["construct", "dose_mM"]).itertuples()
        ) if len(resolved) else "none"
        interval_line = (
            f"n = {n_plates} biological replicates (3 plates x 3 wells per condition). "
            f"Whiskers are 95% cluster-bootstrap intervals on the log ratio, plate as "
            f"the resampling unit, t-inflated for the replicate count -- "
            f"analysis/uncertainty.py::fold_change, which refuses one below "
            f"{MIN_PLATES_FOR_INTERVAL} plates.\n"
            f"{len(resolved)} of {len(estimable)} corrected calls clear 1.0: {calls}. "
            f"Every hollow marker is INCONCLUSIVE -- the point estimate is not the "
            f"result, the interval is. Clearing 1.0 here is not the end of it: "
            f"fig04_interval_robustness sweeps the two analysis knobs behind these "
            f"intervals, and two of the resolved calls lose their lower bound at an "
            f"extreme late window. Read alongside it before quoting a count."
        )
    else:
        interval_line = (
            f"n = {n_plates} biological replicates. NO INTERVAL IS ESTIMABLE: "
            f"analysis/uncertainty.py::fold_change refuses a cluster bootstrap below "
            f"{MIN_PLATES_FOR_INTERVAL} plates, so these are point estimates and the "
            f"absence of error bars is the result, not an omission."
        )

    _note(fig, (
        f"{interval_line}\n"
        "Shaded bands mark the doses docs/FINDINGS.md calls artefacts: 0.1-0.2 mM H2O2, "
        "where the naive readout rises and the corrected one does not, and 2-5 mM DTT, "
        "where the largest apparent inductions vanish under correction and the folds "
        "collapse to about 1.0 as growth collapses.\n"
        "Grey bars mark doses where recovered promoter activity went negative on at "
        "least one plate -- the culture is dying and there is no ratio to take. "
        "Source: outputs/sensor_characterisation.csv "
        "(python3 scripts/run_sensor_characterisation.py)."
    ))
    return save_figure(fig, out_dir, stem)


# ---------------------------------------------------------------------------
# 03 - the decoupling grid
# ---------------------------------------------------------------------------

DESIGN_CONSTRUCT = "UPRE2"
DESIGN_DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0)
NUTRIENT_LADDERS = {
    "dose only (as run)": (1.0,),
    "dose x 2 nutrients": (1.0, 0.5),
    "dose x 3 nutrients": (1.0, 0.6, 0.35),
}
"""The three designs, with exactly the doses and nutrient levels ``scripts/run_power.py``
scored -- so the collinearity drawn here is the collinearity in
``outputs/power_analysis.csv`` and not a second, differently-parameterised number."""


def decoupling_designs(construct: str = DESIGN_CONSTRUCT) -> pd.DataFrame:
    """Every planned condition in the three designs, with its two axes.

    Generated rather than read: ``generator/design.py`` is deterministic and costs
    nothing, and the whole point of a planned design is that it exists before any plate
    does.

    Args:
        construct: Which literature parameter set to plan against.

    Returns:
        Long frame with ``design``, ``dose_mM``, ``nutrient_factor``, ``growth_rate``,
        ``promoter_activity`` and the design's ``collinearity``.
    """
    params = parameters_from_literature(construct)
    rows = []
    for design, nutrients in NUTRIENT_LADDERS.items():
        conditions = (
            stressor_only_series(params, doses=DESIGN_DOSES)
            if nutrients == (1.0,)
            else decoupling_grid(params, doses=DESIGN_DOSES, nutrient_factors=nutrients)
        )
        score = collinearity(conditions)
        for condition in conditions:
            rows.append({
                "design": design,
                "dose_mM": condition.dose_mM,
                "nutrient_factor": condition.nutrient_factor,
                "growth_rate": condition.growth_rate,
                "promoter_activity": condition.promoter_activity,
                "collinearity": score,
                "construct": construct,
            })
    return pd.DataFrame(rows)


def decoupling_grid_scatter(
    designs: pd.DataFrame, out_dir: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure 03: why the collinearity number is the design's whole problem.

    The dose-only design puts every condition on one curve: knowing the dose is knowing
    the growth rate, so no fit can say which of them moved the reporter. Crossing the
    ladder with nutrient limitation lifts the same dose to several growth rates and
    fills the plane.

    Args:
        designs: Output of :func:`decoupling_designs`.
        out_dir: Directory for the SVG and PNG.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``designs`` is empty or lacks a required column.
    """
    required = {"design", "dose_mM", "nutrient_factor", "growth_rate", "collinearity"}
    missing = required - set(designs.columns)
    if missing:
        raise ValueError(f"design table is missing {sorted(missing)}")
    if designs.empty:
        raise ValueError("nothing to plot: the design table has no rows")

    panels = [d for d in NUTRIENT_LADDERS if d in set(designs.design)]
    panels += [d for d in sorted(set(designs.design)) if d not in panels]
    nutrients = sorted(set(designs.nutrient_factor), reverse=True)
    colours = series_colours(len(nutrients))
    markers = ("o", "s", "^", "D", "v", "P", "X", "*")

    fig, axes = _figure(1, len(panels), figsize=(4.2 * len(panels), 4.6), squeeze=False)
    flat = list(axes[0])
    x_lo = float(designs.growth_rate.min())
    x_hi = float(designs.growth_rate.max())
    pad = (x_hi - x_lo) * 0.12 or 0.05

    for ax, design in zip(flat, panels):
        panel = designs[designs.design == design]
        for index, nutrient in enumerate(nutrients):
            points = panel[np.isclose(panel.nutrient_factor, nutrient)]
            if points.empty:
                continue
            ax.plot(points.growth_rate, points.dose_mM, ls="-", lw=1.0, alpha=0.45,
                    color=colours[index], zorder=2)
            ax.scatter(points.growth_rate, points.dose_mM, s=64,
                       facecolor=colours[index], edgecolor="white", linewidth=0.8,
                       marker=markers[index % len(markers)], zorder=3,
                       label=f"nutrients x{nutrient:g}")
        score = float(panel.collinearity.iloc[0])
        ax.set_title(f"{design}\ncollinearity {score:.3f}")
        ax.set_xlabel("specific growth rate (1/h)")
        ax.set_xlim(x_lo - pad, x_hi + pad)
        ax.set_ylim(-0.12, float(designs.dose_mM.max()) * 1.14)
        if ax is flat[0]:
            ax.set_ylabel("planned stressor dose (mM)")
        # Every design's conditions run from slow-and-dosed to fast-and-undosed, so the
        # upper-right corner is empty in all three panels.
        ax.legend(loc="upper right", fontsize=8.5)

    fig.suptitle(
        "Dose and growth rate are one axis until the design separates them",
        fontsize=12.5, fontweight="bold",
    )
    _note(fig, (
        "Generated by ystwin.generator.design from literature parameters for "
        f"{designs.construct.iloc[0] if 'construct' in designs else DESIGN_CONSTRUCT}; "
        f"doses {', '.join(f'{d:g}' for d in DESIGN_DOSES)} mM. Deterministic -- no "
        "random draw, no plate.\n"
        "collinearity() is |correlation| between promoter activity and growth rate "
        "across the planned conditions; 1.000 means no fit can tell them apart at any "
        "sample size. The same three designs and the same scores appear in "
        "outputs/power_analysis.csv."
    ))
    return save_figure(fig, out_dir, "03_decoupling_grid")


# ---------------------------------------------------------------------------
# 04 - attribution power: one reporter, a nutrient axis, a co-expressed ratio
# ---------------------------------------------------------------------------

POWER_FOLDS = (1.2, 1.5, 2.0)
POWER_REPLICATES = 6
RATIO_PAIR = ("YFP", "eGFP")
"""Sensor and reference fluorophore for the co-expressed arm.

Maturation must match or growth does not cancel in the ratio. YFP against eGFP departs
from invariance by 0.017; mOrange2 against YFP departs by 0.226, which is larger than the
1.2-1.5x effects being measured, so the pairing is a design choice and not a detail.
"""


def attribution_power_table(
    n_simulations: int = 240,
    seed: int = 0,
    folds: Sequence[float] = POWER_FOLDS,
    n_replicates: int = POWER_REPLICATES,
    construct: str = DESIGN_CONSTRUCT,
) -> pd.DataFrame:
    """Power to attribute a dose response to dose, for three ways of running the assay.

    Regenerated rather than read. ``outputs/power_analysis.csv`` carries point estimates
    only -- no Wilson bounds, and no co-expressed-ratio arm -- and a power quoted without
    its Monte Carlo interval is a number whose last digit is noise.

    The single-reporter arms use ``discrimination_power``, which puts growth rate in the
    model alongside dose; that is the question that matters, because detecting a dose
    slope is easy on any design. The ratio arm uses ``ratio_discrimination_power`` with
    growth left out, since growth cancels in a ratio of two reporters in one cell and
    carrying the term back would inflate the dose standard error on a collinear design.

    Args:
        n_simulations: Simulated campaigns per cell. The Wilson interval reported is the
            precision this bought.
        seed: Random seed, so the table is reproducible.
        folds: True inductions to simulate.
        n_replicates: Biological replicates per condition.
        construct: Literature parameter set to plan against.

    Returns:
        One row per arm and effect size, with ``power``, ``power_low``, ``power_high``,
        ``n_simulations``, ``collinearity`` and the ``well_cv`` in force.

    Raises:
        ValueError: if ``n_simulations`` is not positive.
    """
    if n_simulations <= 0:
        raise ValueError("n_simulations must be positive")
    params = parameters_from_literature(construct)
    well_cv = float(DEFAULT_WELL_CV)
    dose_only = stressor_only_series(params, doses=DESIGN_DOSES)
    crossed = decoupling_grid(params, doses=DESIGN_DOSES,
                              nutrient_factors=NUTRIENT_LADDERS["dose x 3 nutrients"])
    stress_spec, reference_spec = PANEL[RATIO_PAIR[0]], PANEL[RATIO_PAIR[1]]

    arms = (
        ("one reporter,\ndose only", dose_only,
         lambda conds, fold: discrimination_power(
             params, conds, fold, n_replicates, n_simulations=n_simulations, seed=seed)),
        ("one reporter,\n+ nutrient axis", crossed,
         lambda conds, fold: discrimination_power(
             params, conds, fold, n_replicates, n_simulations=n_simulations, seed=seed)),
        (f"co-expressed ratio,\ndose only ({RATIO_PAIR[0]}/{RATIO_PAIR[1]})", dose_only,
         lambda conds, fold: ratio_discrimination_power(
             params, conds, fold, n_replicates, stress_spec, reference_spec,
             n_simulations=n_simulations, seed=seed)),
    )

    rows = []
    for arm, conditions, estimate_for in arms:
        score = collinearity(conditions)
        for fold in folds:
            estimate = estimate_for(conditions, fold)
            rows.append({
                "arm": arm.replace("\n", " "),
                "arm_label": arm,
                "induction_fold": float(fold),
                "n_replicates": int(n_replicates),
                "power": float(estimate),
                "power_low": estimate.low,
                "power_high": estimate.high,
                "n_simulations": estimate.n_simulations,
                "collinearity": score,
                "well_cv": well_cv,
            })
    return pd.DataFrame(rows)


def attribution_power(
    power: pd.DataFrame, out_dir: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure 04: what the design and the strain each buy, with the simulation's error.

    Args:
        power: Output of :func:`attribution_power_table`.
        out_dir: Directory for the SVG and PNG.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``power`` is empty, lacks a column, or has more arms than the
            palette can separate.
    """
    required = {"arm", "induction_fold", "power", "power_low", "power_high",
                "n_simulations"}
    missing = required - set(power.columns)
    if missing:
        raise ValueError(f"power table is missing {sorted(missing)}")
    if power.empty:
        raise ValueError("nothing to plot: the power table has no rows")

    arms = list(dict.fromkeys(power.arm))
    labels = {a: (power[power.arm == a].arm_label.iloc[0]
                  if "arm_label" in power.columns else a) for a in arms}
    folds = sorted(set(power.induction_fold))
    colours = series_colours(len(arms))
    width = 0.8 / len(arms)

    fig, ax = _figure(figsize=(11.4, 5.6))
    x = np.arange(len(folds), dtype=float)
    for index, (arm, colour) in enumerate(zip(arms, colours)):
        offset = (index - (len(arms) - 1) / 2) * width
        block = power[power.arm == arm].set_index("induction_fold").reindex(folds)
        heights = block.power.to_numpy(dtype=float) * 100
        low = np.clip(heights - block.power_low.to_numpy(dtype=float) * 100, 0, None)
        high = np.clip(block.power_high.to_numpy(dtype=float) * 100 - heights, 0, None)
        ax.bar(x + offset, heights, width=width * 0.9, color=colour,
               edgecolor="white", linewidth=0.8, hatch=_hatch(index), zorder=2,
               label=labels[arm].replace("\n", " "))
        ax.errorbar(x + offset, heights, yerr=np.vstack([low, high]), fmt="none",
                    ecolor="#222222", elinewidth=1.3, capsize=4, zorder=4)
        # Above the whisker, not above the bar: a value label struck through by its own
        # confidence interval is the one thing this figure must not do.
        tops = block.power_high.to_numpy(dtype=float) * 100
        for xi, h, cap in zip(x + offset, heights, tops):
            ax.text(xi, cap + 3.0, f"{h:.0f}%", ha="center", va="bottom", fontsize=8.5,
                    fontweight="bold", color=colour, zorder=5)

    ax.axhline(80, color="#555555", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.text(-0.55, 81.5, "80% target", fontsize=8.5, va="bottom", ha="left",
            color="#555555")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{f:g}x" for f in folds])
    ax.set_xlabel("true induction at saturating dose")
    ax.set_ylabel("power to attribute the response to dose (%)")
    ax.set_ylim(0, 118)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xlim(-0.6, len(folds) - 0.4)
    ax.set_title("Detecting a dose slope is easy; attributing it to dose is not")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), ncols=1, fontsize=9)

    n_sims = int(power.n_simulations.max())
    n_reps = int(power.n_replicates.iloc[0]) if "n_replicates" in power.columns else 6
    well_cv = float(power.well_cv.iloc[0]) if "well_cv" in power.columns else float("nan")
    _note(fig, (
        f"Simulated: {n_sims} campaigns per bar, {n_reps} biological replicates, "
        "ystwin.analysis.power, seed 0. Whiskers are Wilson score intervals on the\n"
        "simulated proportion -- the precision the simulation bought, not biological "
        "spread. The single-reporter arms put growth rate in the model alongside dose;\n"
        "the ratio arm leaves it out, because growth cancels between two reporters in "
        "one cell and putting it back on a collinear design throws that away again.\n"
        "outputs/power_analysis.csv carries the two single-reporter arms as bare point "
        "estimates, with no interval and no ratio arm, which is why this figure\n"
        f"regenerates rather than reads. Its numbers (0%, 99%) also predate "
        f"analysis/power.py's well-level variance term, here well_cv={well_cv:g}, which "
        "that\nmodule's own docstring calls the term that limits power; well_cv=0 "
        "reproduces the committed table exactly."
    ))
    return save_figure(fig, out_dir, "04_attribution_power")


# ---------------------------------------------------------------------------
# 05 - G4: how much of the anchor signal is genomic DNA
# ---------------------------------------------------------------------------

CORRECTABLE_GDNA_CEILING = 0.60
"""Share of a qPCR signal above which no published gDNA correction holds.

Mirrors ``qpcr._MAX_CORRECTABLE_GDNA_SHARE``. Restated here rather than imported because
it is a private name there, and a figure that silently tracked a private constant would
move without notice.
"""

PASS_MARGIN_CYCLES = 3.0
"""The G4 no-RT margin, in cycles. Mirrors ``qpcr._MIN_RT_MINUS_MARGIN``."""


def g4_contamination_table(
    rt_minus: pd.DataFrame, verdicts: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Every no-RT reading as a contamination share rather than a cycle count.

    The instrument prints cycles, but cycles are a logarithmic proxy for the quantity
    that decides whether a reading means anything. ``qpcr.gdna_share`` converts exactly,
    and the conversion is worth making visible: 0.6 cycles sounds only a little worse
    than 1.0 and is 66% contamination against 50%, which straddles the limit of what any
    correction can undo.

    Args:
        rt_minus: The ``outputs/g4_rt_minus_qc.csv`` frame: ``construct``, ``target``,
            ``margin_cycles`` at minimum.
        verdicts: Optional ``outputs/g4_verdicts.csv`` frame, merged on ``construct`` so
            each row carries its gate verdict.

    Returns:
        The readings with ``gdna_share``, a ``role`` naming what the target is for that
        construct, and ``correctable``.

    Raises:
        ValueError: if a required column is absent or the frame is empty.
    """
    required = {"construct", "target", "margin_cycles"}
    missing = required - set(rt_minus.columns)
    if missing:
        raise ValueError(f"RT-minus QC table is missing {sorted(missing)}")
    if rt_minus.empty:
        raise ValueError("RT-minus QC table has no readings in it")

    out = rt_minus.copy()
    out["gdna_share"] = gdna_share(out.margin_cycles.to_numpy(dtype=float))
    out["role"] = [
        "anchor" if t == TARGET_FOR_CONSTRUCT.get(c)
        else ("reference" if t == REFERENCE_TARGET else "other")
        for c, t in zip(out.construct, out.target)
    ]
    out["correctable"] = out.gdna_share.notna() & (out.gdna_share <= CORRECTABLE_GDNA_CEILING)
    if verdicts is not None and not verdicts.empty:
        if "construct" not in verdicts.columns or "verdict" not in verdicts.columns:
            raise ValueError("verdict table needs 'construct' and 'verdict' columns")
        out = out.merge(verdicts[["construct", "verdict"]].drop_duplicates(),
                        on="construct", how="left")
    return out


def g4_contamination(
    readings: pd.DataFrame, out_dir: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure 05: the G4 verdict, on the axis that explains it.

    Four constructs, four INCONCLUSIVE verdicts, two different reasons. The oxidative
    anchor is clean and its verdict turns on replicate count; the ER anchor is not
    measuring transcript at all, and sits past the share above which no correction
    applies.

    Args:
        readings: Output of :func:`g4_contamination_table`.
        out_dir: Directory for the SVG and PNG.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``readings`` is empty or lacks a required column.
    """
    required = {"construct", "target", "gdna_share", "role"}
    missing = required - set(readings.columns)
    if missing:
        raise ValueError(f"contamination table is missing {sorted(missing)}")
    if readings.empty:
        raise ValueError("nothing to plot: the contamination table has no readings")

    plotted = readings[readings.role.isin(("anchor", "reference"))].copy()
    if plotted.empty:
        raise ValueError(
            "no anchor or reference readings in the table; qpcr.TARGET_FOR_CONSTRUCT "
            "did not recognise any construct/target pair"
        )
    constructs = sorted(set(plotted.construct),
                        key=lambda c: (TARGET_FOR_CONSTRUCT.get(c, ""), c))
    roles = [r for r in ("anchor", "reference") if r in set(plotted.role)]
    colours = series_colours(len(roles))
    markers = {"anchor": "o", "reference": "^"}
    offsets = {"anchor": 0.17, "reference": -0.17}

    fig, ax = _figure(figsize=(10.6, 6.2))
    # Room on the left for the direct target labels, on the right for the medians, and a
    # clear band at the top for the two reference lines -- none of which may sit on data.
    lo = max(float(plotted.gdna_share.min()) * 0.12, 1e-5)
    hi = float(plotted.gdna_share.max()) * 1.7

    ax.axvspan(1.0, hi, color="#000000", alpha=0.07, lw=0, zorder=0)
    ax.axvline(CORRECTABLE_GDNA_CEILING, color=PALETTE[1], lw=2.0, zorder=1)
    ax.axvline(2.0 ** -PASS_MARGIN_CYCLES, color="#555555", lw=1.4, ls=(0, (4, 3)),
               zorder=1)

    for y, construct in enumerate(constructs):
        for role, colour in zip(roles, colours):
            points = plotted[(plotted.construct == construct) & (plotted.role == role)]
            if points.empty:
                continue
            target = str(points.target.iloc[0])
            ys = np.full(len(points), y + offsets[role], dtype=float)
            ax.scatter(points.gdna_share, ys, s=42, marker=markers[role],
                       facecolor=colour, edgecolor="white", linewidth=0.7,
                       alpha=0.85, zorder=3)
            median = float(points.gdna_share.median())
            ax.scatter([median], [y + offsets[role]], s=170, marker="|",
                       color="#000000", linewidth=2.2, zorder=4)
            ax.text(lo * 1.05, y + offsets[role], f"{target}  ({role})", fontsize=8.5,
                    va="center", ha="left", color=colour, fontweight="bold", zorder=5)
            # Outside the axes, in axes fractions, so a reading near the top of the
            # range can never be overwritten by its own construct's median label.
            ax.annotate(f"median {median:.1%}", (1.012, y + offsets[role]),
                        xycoords=("axes fraction", "data"), fontsize=8.5,
                        va="center", ha="left", color=colour, annotation_clip=False)

    if "verdict" in plotted.columns:
        for y, construct in enumerate(constructs):
            verdict = plotted[plotted.construct == construct].verdict.dropna()
            if len(verdict):
                ax.text(lo * 1.05, y + 0.42, f"G4: {verdict.iloc[0]}", fontsize=8.0,
                        va="center", ha="left", color="#444444", style="italic")

    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.62, len(constructs) - 0.5 + 1.4)
    ax.set_yticks(range(len(constructs)))
    ax.set_yticklabels(constructs, fontweight="bold")
    ax.set_xlabel("share of the +RT signal that is genomic DNA "
                  "(log scale; from the no-RT margin, at 100% efficiency)")
    ax.set_title("G4: the ER anchor is not measuring transcript, in all three replicates")
    ax.grid(axis="y", visible=False)

    # Rotated labels on the two reference lines would lie across the readings, so the
    # lines are keyed in the clear band above the topmost construct instead.
    ax.legend(
        handles=[
            Line2D([], [], color=PALETTE[1], lw=2.0,
                   label=f"{CORRECTABLE_GDNA_CEILING:.0%} -- past this, no published "
                         f"gDNA correction holds"),
            Line2D([], [], color="#555555", lw=1.4, ls=(0, (4, 3)),
                   label=f"G4 pass: no-RT margin >= {PASS_MARGIN_CYCLES:g} cycles "
                         f"({2.0 ** -PASS_MARGIN_CYCLES:.1%})"),
            Patch(facecolor="#000000", alpha=0.07,
                  label="above 100%: the no-RT control amplified before the +RT sample"),
            Line2D([], [], color="#000000", lw=2.2, marker="|", ls="none", markersize=11,
                   label="median of the readings"),
        ],
        loc="upper left", ncols=2, fontsize=8.5, handlelength=1.8,
    )

    _note(fig, (
        "Every no-RT reading, converted with ystwin.qpcr.gdna_share; black bar is the "
        "median. Read from outputs/g4_rt_minus_qc.csv and outputs/g4_verdicts.csv.\n"
        "The two ER constructs share the Hac1 anchor and the logbook names the primer "
        "pair 'Hac1 Exon' -- exon primers amplify the genomic locus readily, which is "
        "both failures at once.\n"
        "All four verdicts are INCONCLUSIVE, for two different reasons: the ER anchor is "
        "past the correction ceiling, and the oxidative anchor is clean but too noisy at "
        "the replicate count available."
    ))
    return save_figure(fig, out_dir, "05_g4_contamination")


# ---------------------------------------------------------------------------
# 06 - the held-out prediction, and the dose range where it fails
# ---------------------------------------------------------------------------


def identifiable_ec50(
    readings: pd.DataFrame, value: str = "activity_late",
) -> pd.DataFrame:
    """Refit the repo's own biphasic dose response per construct, and say which it vouches for.

    The held-out score is reported against ``dose / EC50`` rather than against dose,
    because that is the axis on which its sign separates. No committed table carries the
    fitted EC50 -- ``outputs/panel_calibration.csv`` leaves it empty, because that module
    refuses a fit it cannot identify -- so it is refitted here from the committed
    characterisation table with ``generator/panel_calibration.fit_dose_response``, the
    same function ``scripts/run_heldout_score.py`` scores with. Refitted rather than
    transcribed from ``docs/WHY_INTERPOLATION_WORKS.md``: a number copied out of prose
    stops tracking the data the moment a plate is added, which is exactly what happened
    to the table this figure plots.

    The fit here sees every plate, where the scorer's fits each see one partition's
    training rows, so these EC50s are close to but not identical with the ones the doc
    quotes. That is stated on the figure rather than smoothed over.

    Args:
        readings: The ``outputs/sensor_characterisation.csv`` frame: ``construct``,
            ``dose_mM`` and ``value`` at minimum.
        value: Column to fit. The default is the dilution-corrected activity, which is
            what the held-out score predicts.

    Returns:
        One row per construct with ``ec50``, ``lethal_dose``, ``r_squared``,
        ``identifiable`` and, when it is not, ``unidentified_reason``.

    Raises:
        ValueError: if a required column is absent or the frame is empty.
    """
    required = {"construct", "dose_mM", value}
    missing = required - set(readings.columns)
    if missing:
        raise ValueError(
            f"readings are missing {sorted(missing)}; expected the columns "
            f"scripts/run_sensor_characterisation.py writes"
        )
    if readings.empty:
        raise ValueError("readings frame has no rows")

    rows = []
    for construct, group in readings.groupby("construct", sort=True):
        doses = group.dose_mM.to_numpy(dtype=float)
        activity = group[value].to_numpy(dtype=float)
        usable = np.isfinite(doses) & np.isfinite(activity)
        row = {
            "construct": str(construct),
            "stressor": str(group.stressor.iloc[0]) if "stressor" in group else "",
            "n_readings": int(usable.sum()),
            "ec50": np.nan, "lethal_dose": np.nan, "r_squared": np.nan,
            "identifiable": False, "unidentified_reason": "",
        }
        if usable.sum() < 5:
            row["unidentified_reason"] = (
                f"{int(usable.sum())} usable readings; the four-parameter fit needs 5"
            )
            rows.append(row)
            continue
        try:
            fit = fit_dose_response(doses[usable], activity[usable])
        except (ValueError, RuntimeError) as exc:
            row["unidentified_reason"] = f"fit failed: {exc}"
            rows.append(row)
            continue
        row.update({
            "ec50": float(fit.ec50), "lethal_dose": float(fit.lethal_dose),
            "r_squared": float(fit.r_squared), "identifiable": bool(fit.identifiable),
            "unidentified_reason": str(fit.unidentified_reason),
        })
        rows.append(row)
    return pd.DataFrame(rows)


def heldout_skill_table(
    by_dose: pd.DataFrame, ec50: pd.DataFrame, readings: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """The held-out skill per dose, placed on the EC50 axis, and checked for staleness.

    Two things are added to the committed table. The first is ``dose_over_ec50``, the
    axis on which the sign of the skill separates. The second is a freshness check, and
    it is the reason ``readings`` is accepted at all: this table is written by a script
    that reads the characterisation table, so if the characterisation table has grown
    since, every number here is out of date and the figure would be a lie told with real
    data. The check is exact rather than a guess at a date -- the scorer records how many
    fits it vouched for, and the readings table says how many wells each construct has at
    that dose, so the product is how many test rows a current run would have scored.

    Args:
        by_dose: The ``outputs/heldout_interpolation_by_dose.csv`` frame.
        ec50: Output of :func:`identifiable_ec50`. Only rows it vouches for are used.
        readings: Optional ``outputs/sensor_characterisation.csv`` frame. Without it the
            staleness columns are NaN and the figure says the check was not run.

    Returns:
        ``by_dose`` with ``scoreable``, ``beats_baseline``, ``reference_ec50_mM``,
        ``dose_over_ec50``, and -- when ``readings`` is given -- ``n_test_expected`` and
        ``is_stale``.

    Raises:
        ValueError: if a required column is absent, the frame is empty, or no fit in
            ``ec50`` is identifiable, in which case there is no axis to draw on.
    """
    required = {"held_out_dose", "skill", "n_test", "n_fits_vouched"}
    missing = required - set(by_dose.columns)
    if missing:
        raise ValueError(
            f"held-out table is missing {sorted(missing)}; expected the columns "
            f"scripts/run_heldout_score.py writes"
        )
    if by_dose.empty:
        raise ValueError("held-out table has no rows")
    vouched = ec50[ec50.identifiable.astype(bool) & np.isfinite(ec50.ec50)]
    if vouched.empty:
        raise ValueError(
            "no construct has an identifiable dose-response fit, so there is no EC50 to "
            "put the held-out doses on. Every scored prediction comes from a vouched "
            "fit; if there are none there is also nothing to score."
        )

    out = by_dose.copy()
    reference = float(np.median(vouched.ec50.to_numpy(dtype=float)))
    out["reference_ec50_mM"] = reference
    out["dose_over_ec50"] = out.held_out_dose.astype(float) / reference
    out["scoreable"] = np.isfinite(out.skill.to_numpy(dtype=float))
    out["beats_baseline"] = out.scoreable & (out.skill.to_numpy(dtype=float) > 0.0)

    if readings is None:
        out["n_test_expected"] = np.nan
        out["is_stale"] = pd.NA
        return out

    expected = []
    for row in out.itertuples():
        at_dose = readings[np.isclose(readings.dose_mM.astype(float), row.held_out_dose)]
        if at_dose.empty or not np.isfinite(row.n_fits_vouched):
            expected.append(np.nan)
            continue
        per_construct = at_dose.groupby("construct").size().max()
        expected.append(float(row.n_fits_vouched) * float(per_construct))
    out["n_test_expected"] = expected
    out["is_stale"] = [
        pd.NA if not np.isfinite(e) or not np.isfinite(n) else bool(n < e)
        for n, e in zip(out.n_test.astype(float), out.n_test_expected)
    ]
    return out


def _skill_sign_crossover(table: pd.DataFrame) -> tuple[float, float] | None:
    """The bracket in ``dose / EC50`` between the last win and the first loss.

    Returns None when the sign does not separate cleanly, which is the case worth
    knowing about: the claim is that it does.
    """
    scored = table[table.scoreable].sort_values("dose_over_ec50")
    if scored.empty:
        return None
    wins = scored[scored.beats_baseline]
    losses = scored[~scored.beats_baseline]
    if wins.empty or losses.empty:
        return None
    last_win = float(wins.dose_over_ec50.max())
    first_loss = float(losses.dose_over_ec50.min())
    if last_win >= first_loss:
        return None
    return last_win, first_loss


def heldout_skill_vs_saturation(
    table: pd.DataFrame, out_dir: pathlib.Path,
    stem: str = "fig02_heldout_skill_vs_saturation",
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure: the held-out prediction wins near the EC50 and loses in saturation.

    The zero line is the whole figure. Below it the fitted dose response is worse than
    carrying the nearest measured dose forward, which is the baseline any future
    prediction here has to beat, so a point below zero is a negative result and is drawn
    as one.

    The lower panel is the mechanism rather than more data. The induction term is
    Michaelis-Menten, ``d / (K + d)``, whose sensitivity to dose is ``K / (K + d)^2``.
    In units of ``u = d / K`` that sensitivity relative to its value at ``u = 1`` is
    ``4 / (1 + u)^2``, with no free parameter at all -- so the curve is a property of the
    functional form and not a fit to these five points. Deep in saturation the model
    cannot move its prediction with dose, and the residual is whatever the plateau height
    happens to be.

    Args:
        table: Output of :func:`heldout_skill_table`.
        out_dir: Directory for the SVG and PNG.
        stem: Filename without extension.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if ``table`` lacks the builder's columns or carries no scoreable row.
    """
    required = {"held_out_dose", "skill", "n_test", "dose_over_ec50", "scoreable",
                "beats_baseline", "reference_ec50_mM"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"held-out skill table is missing {sorted(missing)}")
    scored = table[table.scoreable].sort_values("dose_over_ec50")
    if scored.empty:
        raise ValueError(
            "nothing to plot: no dose in the held-out table has a skill score. The "
            "unscoreable rows carry the reason in their `note` column."
        )

    win_colour, loss_colour = PALETTE[2], PALETTE[1]
    fig, axes = _figure(2, 1, figsize=(9.6, 7.6), sharex=True,
                        height_ratios=(1.55, 1.0))
    top_ax, bottom_ax = axes

    crossover = _skill_sign_crossover(scored)
    if crossover is not None:
        top_ax.axvspan(crossover[1] * 0.999, float(scored.dose_over_ec50.max()) * 4.0,
                       color="#000000", alpha=0.055, lw=0, zorder=0)
        bottom_ax.axvspan(crossover[1] * 0.999,
                          float(scored.dose_over_ec50.max()) * 4.0,
                          color="#000000", alpha=0.055, lw=0, zorder=0)

    top_ax.plot(scored.dose_over_ec50, scored.skill, ls="-", lw=1.2, color="#888888",
                zorder=2)
    top_ax.axhline(0.0, color="#000000", lw=1.9, zorder=3)
    for row in scored.itertuples():
        wins = bool(row.beats_baseline)
        top_ax.plot([row.dose_over_ec50], [row.skill],
                    marker="o" if wins else "v", ls="none", markersize=10,
                    color=win_colour if wins else loss_colour,
                    markerfacecolor=win_colour if wins else "#FFFFFF",
                    markeredgecolor=win_colour if wins else loss_colour,
                    markeredgewidth=1.8, zorder=4)
        top_ax.annotate(
            f"{row.held_out_dose:g} mM\nn = {int(row.n_test)}",
            (row.dose_over_ec50, row.skill), textcoords="offset points",
            xytext=(0, 15 if row.skill > 0 else -30), ha="center", fontsize=8.2,
            color="#333333",
        )

    reference = float(scored.reference_ec50_mM.iloc[0])
    lo = float(scored.dose_over_ec50.min()) * 0.55
    hi = float(scored.dose_over_ec50.max()) * 1.9
    span = float(np.nanmax(np.abs(scored.skill.to_numpy(dtype=float)))) * 1.55
    top_ax.set_xscale("log")
    top_ax.set_xlim(lo, hi)
    top_ax.set_ylim(-span, span)
    top_ax.set_ylabel("skill against the nearest-dose baseline\n"
                      "(1 = perfect, 0 = no better than the baseline)")
    top_ax.set_title("Below the zero line the model is worse than carrying the nearest "
                     "measured dose forward")
    top_ax.legend(handles=[
        Line2D([], [], color=win_colour, marker="o", ls="none", markersize=9,
               markerfacecolor=win_colour, label="beats the baseline"),
        Line2D([], [], color=loss_colour, marker="v", ls="none", markersize=9,
               markerfacecolor="#FFFFFF", markeredgewidth=1.8,
               label="LOSES to the baseline"),
        Line2D([], [], color="#000000", lw=1.9, label="no skill (0)"),
    ], loc="lower left", fontsize=8.5, ncols=1)

    u = np.geomspace(lo, hi, 400)
    bottom_ax.plot(u, 4.0 / (1.0 + u) ** 2, ls="-", lw=2.0, color=PALETTE[0], zorder=2)
    bottom_ax.plot(scored.dose_over_ec50, 4.0 / (1.0 + scored.dose_over_ec50) ** 2,
                   marker="D", ls="none", markersize=7, color=PALETTE[0],
                   markerfacecolor="#FFFFFF", markeredgewidth=1.6, zorder=3)
    bottom_ax.axhline(1.0, color="#555555", lw=1.0, ls=(0, (4, 3)), zorder=1)
    bottom_ax.set_xscale("log")
    bottom_ax.set_yscale("log")
    bottom_ax.set_xlim(lo, hi)
    bottom_ax.set_xlabel("held-out dose, in multiples of the fitted EC50 "
                         f"(EC50 = {reference:.3f} mM, log scale)")
    bottom_ax.set_ylabel("d/dd [ d/(K+d) ], relative\nto its value at d = K")
    bottom_ax.set_title("Why: the model's sensitivity to dose falls as 1/dose^2 past its "
                        "own EC50")

    ticks = [float(v) for v in scored.dose_over_ec50]
    bottom_ax.set_xticks(ticks)
    bottom_ax.set_xticklabels([f"{t:.1f}" for t in ticks])
    bottom_ax.minorticks_off()

    stale = [bool(s) for s in table.get("is_stale", []) if s is not pd.NA and s is not None]
    if not stale:
        freshness = ("Staleness check not run: the characterisation table was not passed "
                     "to the table builder.")
    elif any(stale):
        freshness = (
            "STALE: the scored test rows are fewer than the current "
            "outputs/sensor_characterisation.csv would produce, so this table predates a "
            "plate. Regenerate with `python3 scripts/run_heldout_score.py` before "
            "quoting these numbers."
        )
    else:
        freshness = (
            "Freshness checked: every partition's scored test rows equal the vouched fits "
            "times the wells per construct in the current "
            "outputs/sensor_characterisation.csv, so the table is not stale."
        )

    crossover_line = (
        f"The sign separates on dose/EC50, with the crossover between "
        f"{crossover[0]:.1f} and {crossover[1]:.1f} multiples of the fitted EC50."
        if crossover is not None else
        "The sign does NOT separate cleanly on dose/EC50 in this table, which is the "
        "claim this figure exists to check."
    )
    fig.suptitle("The fit beats the baseline near its own EC50, and loses in saturation",
                 fontsize=12.5, fontweight="bold")
    _note(fig, (
        f"Skill and test-row counts from outputs/heldout_interpolation_by_dose.csv "
        f"(python3 scripts/run_heldout_score.py). {freshness}\n"
        f"{crossover_line} The EC50 is refitted per construct from "
        f"outputs/sensor_characterisation.csv with "
        f"generator/panel_calibration.fit_dose_response and the median over the fits it "
        f"vouches for is used; because that fit sees every plate rather than one "
        f"partition's training rows, it is close to but not identical with the per-"
        f"partition EC50s docs/WHY_INTERPOLATION_WORKS.md quotes.\n"
        "What to conclude: on an interior dose at or below the crossover the fitted "
        "dose response beats persistence, MODESTLY -- and above it, it is decisively "
        "worse, so a prediction from this model outside its EC50 range should not be "
        "used. The lower panel is analytic and has no free parameter."
    ))
    return save_figure(fig, out_dir, stem)


# ---------------------------------------------------------------------------
# 07 - what transfer does and does not recover
# ---------------------------------------------------------------------------

PEER_BINS = ("0", "1-2", "3-5", "6+")
"""The peer-count bins ``run_external_validation.py`` writes, in the order they mean.

Read off an axis they must be ordered by peer count and not alphabetically, which would
put ``6+`` before ``1-2``. An empty bin is drawn as an explicit gap rather than skipped,
because "no pairs fell here" and "recovery was zero here" are opposite readings of the
same blank space.
"""

HEADLINE_NULLS = (
    ("transfer (median own-target score)", "rotated_subspace"),
    ("transfer (median own-target score)", "matched_marginals"),
    ("redundancy law, peers from panel target lists", "rotated_subspace"),
    ("known-loadings readout (median Spearman)", "shuffled_loadings"),
)
"""The (claim, null) pairs the transfer figure draws, and why these four.

The first is the load-bearing failure: arbitrary axes of the same rank carry as much
transfer information as the fitted ones. The second is the null the same claim does clear,
included because leaving it out would make the failure look like a total absence of
structure when it is not -- it says only that the channels are correlated at all. The
third asks whether the redundancy law is evidence about the fitted basis, and the fourth
is the one claim in the set that beats its null, which is what makes the others legible as
a specific failure rather than a broken pipeline.
"""


def redundancy_recovery_table(redundancy: pd.DataFrame) -> pd.DataFrame:
    """Module recovery against how many other stressors drive the same module.

    Args:
        redundancy: An ``outputs/external_redundancy*.csv`` frame. Its first column is
            the peer bin and is unnamed in the file, so it arrives as ``Unnamed: 0``.

    Returns:
        ``peers``, ``peer_definition``, ``pairs``, ``mean_score``, ``recovered``, with
        ``peers`` ordered by :data:`PEER_BINS` and empty bins kept.

    Raises:
        ValueError: if the frame is empty, lacks a column, or carries a bin this module
            does not know how to order.
    """
    if redundancy.empty:
        raise ValueError("redundancy table has no rows")
    out = redundancy.rename(columns={redundancy.columns[0]: "peers"}).copy()
    required = {"peers", "pairs", "recovered", "peer_definition"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(
            f"redundancy table is missing {sorted(missing)}; expected the columns "
            f"scripts/run_external_validation.py writes"
        )
    out["peers"] = out.peers.astype(str)
    unknown = sorted(set(out.peers) - set(PEER_BINS))
    if unknown:
        raise ValueError(
            f"peer bins {unknown} are not in PEER_BINS {list(PEER_BINS)}, so they cannot "
            f"be ordered on an axis; add them there rather than sorting alphabetically"
        )
    out["peer_index"] = [PEER_BINS.index(p) for p in out.peers]
    return out.sort_values(["peer_definition", "peer_index"]).reset_index(drop=True)


def null_comparison_table(
    panels: Sequence[tuple[str, pd.DataFrame]],
    simulated: pd.DataFrame | None = None,
    simulated_family: str = "baseline",
    simulated_label: str = "simulation, default configuration",
) -> pd.DataFrame:
    """Every claim beside the null it was tested against, on one comparable axis.

    The axis is ``observed - null_median``, the gap, because the claims are scored in
    different units -- a transfer R2 in simulation and a median own-target score on real
    arrays -- and putting those on one score axis would invite a comparison neither
    supports. The gap is a difference of like with like within each row, and zero means
    the same thing on every row: the null did as well as the model.

    Args:
        panels: ``(label, frame)`` pairs, each frame shaped like
            ``outputs/external_nulls*.csv``: ``test``, ``null``, ``observed``,
            ``null_median``, ``p_value``, ``beats_null`` and optionally ``draws``.
        simulated: Optional ``outputs/cross_family_transfer.csv`` frame, whose one
            ``simulated_family`` row carries the simulated transfer number and its
            rotated-subspace null.
        simulated_family: Which family row to take from ``simulated``.
        simulated_label: Dataset label for that row.

    Returns:
        ``claim``, ``dataset``, ``null``, ``observed``, ``null_median``, ``gap``,
        ``p_value``, ``beats_null``, ``draws``, ``headline``.

    Raises:
        ValueError: if no panel is given, a frame lacks a column, or ``simulated`` does
            not carry ``simulated_family``.
    """
    if not panels:
        raise ValueError("no null tables to compare; pass at least one (label, frame)")
    required = {"test", "null", "observed", "null_median", "p_value", "beats_null"}
    rows = []
    for label, frame in panels:
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(
                f"null table {label!r} is missing {sorted(missing)}; expected the "
                f"columns scripts/run_external_validation.py writes"
            )
        if frame.empty:
            raise ValueError(f"null table {label!r} has no rows")
        for row in frame.itertuples():
            rows.append({
                "claim": str(row.test), "dataset": str(label), "null": str(row.null),
                "observed": float(row.observed), "null_median": float(row.null_median),
                "gap": float(row.observed) - float(row.null_median),
                "p_value": float(row.p_value), "beats_null": bool(row.beats_null),
                "draws": float(getattr(row, "draws", np.nan)),
            })

    if simulated is not None and not simulated.empty:
        needed = {"family", "transfer_r2", "null_median", "p_value", "beats_null"}
        missing = needed - set(simulated.columns)
        if missing:
            raise ValueError(f"simulated transfer table is missing {sorted(missing)}")
        chosen = simulated[simulated.family == simulated_family]
        if chosen.empty:
            raise ValueError(
                f"no family {simulated_family!r} in the simulated transfer table; it "
                f"carries {sorted(set(simulated.family))}"
            )
        row = chosen.iloc[0]
        rows.append({
            "claim": "transfer (median own-target score)",
            "dataset": str(simulated_label), "null": "rotated_subspace",
            "observed": float(row.transfer_r2), "null_median": float(row.null_median),
            "gap": float(row.transfer_r2) - float(row.null_median),
            "p_value": float(row.p_value), "beats_null": bool(row.beats_null),
            "draws": np.nan,
        })

    out = pd.DataFrame(rows)
    out["headline"] = [(c, n) in HEADLINE_NULLS for c, n in zip(out.claim, out.null)]
    return out


def transfer_honesty(
    recovery: pd.DataFrame, nulls: pd.DataFrame, out_dir: pathlib.Path,
    stem: str = "fig03_transfer_honesty",
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure: transfer recovers redundant coverage, and fails its rotated-subspace null.

    Two panels, one negative result each. The left one is what transfer does: recovery
    tracks how many *other* stressors in the panel drive the same module, so a module has
    to appear twice in the design to be recoverable at all. The right one is what it does
    not do: the observed score sits on top of a null that keeps the rank and destroys the
    axes, so the identity of the fitted latent basis contributes nothing.

    The right panel is drawn on the gap rather than on the scores. Zero is then the
    result, and a bar that reaches zero is a claim that did not survive.

    Args:
        recovery: Output of :func:`redundancy_recovery_table`.
        nulls: Output of :func:`null_comparison_table`. Only ``headline`` rows are drawn.
        out_dir: Directory for the SVG and PNG.
        stem: Filename without extension.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if either frame is empty, lacks a column, or leaves nothing to draw.
    """
    for name, frame, needed in (
        ("recovery", recovery, {"peers", "peer_definition", "pairs", "recovered"}),
        ("nulls", nulls, {"claim", "dataset", "null", "gap", "p_value", "beats_null",
                          "observed", "headline"}),
    ):
        missing = needed - set(frame.columns)
        if missing:
            raise ValueError(f"{name} table is missing {sorted(missing)}")
        if frame.empty:
            raise ValueError(f"nothing to plot: the {name} table has no rows")
    drawn = nulls[nulls.headline.astype(bool)]
    if drawn.empty:
        raise ValueError(
            "no headline rows in the null table; HEADLINE_NULLS names the (claim, null) "
            "pairs this figure draws and none of them appear"
        )

    fig, axes = _figure(1, 2, figsize=(15.2, 6.8), width_ratios=(1.0, 1.5))
    left, right = axes

    definitions = list(dict.fromkeys(recovery.peer_definition))
    colours = series_colours(len(definitions))
    markers = ("o", "s", "^", "D")
    styles = ("-", (0, (5, 2)), (0, (1, 1.6)), (0, (7, 2, 1, 2)))
    x = np.arange(len(PEER_BINS), dtype=float)
    for index, (definition, colour) in enumerate(zip(definitions, colours)):
        series = recovery[recovery.peer_definition == definition]
        present = series[series.pairs.astype(float) > 0]
        xs = [PEER_BINS.index(p) for p in present.peers]
        left.plot(xs, present.recovered.astype(float), marker=markers[index % 4],
                  ls=styles[index % 4], color=colour, markersize=9, lw=1.9,
                  markerfacecolor=colour, markeredgecolor=colour, zorder=3,
                  label=f"peers from {definition}")
        for row in present.itertuples():
            left.annotate(
                f"{int(row.pairs)} pairs",
                (PEER_BINS.index(row.peers), float(row.recovered)),
                textcoords="offset points", xytext=(0, 13 + 13 * index), ha="center",
                fontsize=8.2, color=colour,
            )
        for row in series[series.pairs.astype(float) == 0].itertuples():
            left.annotate("no pairs", (PEER_BINS.index(row.peers), 0.02),
                          textcoords="offset points", xytext=(0, 13 * index), ha="center",
                          fontsize=8.0, color=colour, style="italic")
    left.set_xticks(x)
    left.set_xticklabels(PEER_BINS)
    left.set_xlim(-0.45, len(PEER_BINS) - 0.55)
    left.set_ylim(-0.06, 1.26)
    left.set_xlabel("peers: other stressors in the panel\ndriving the same module (count)")
    left.set_ylabel("modules recovered (fraction of own-target pairs)")
    left.set_title("A. Transfer recovers what the panel covers twice")
    left.legend(loc="upper left", fontsize=8.5)

    # Failures first. A reader who stops after the top block has read the result, and a
    # reader who goes on sees which nulls the same machinery does clear.
    order = drawn.sort_values(["beats_null", "claim", "dataset"]).reset_index(drop=True)
    labels = []
    for y, row in enumerate(order.itertuples()):
        colour = PALETTE[2] if row.beats_null else PALETTE[1]
        right.plot([0.0, row.gap], [y, y], ls="-", lw=2.2, color=colour, zorder=2)
        right.plot([row.gap], [y], marker="o" if row.beats_null else "X", ls="none",
                   markersize=11, color=colour,
                   markerfacecolor=colour if row.beats_null else "#FFFFFF",
                   markeredgecolor=colour, markeredgewidth=2.0, zorder=3)
        verdict = "beats" if row.beats_null else "DOES NOT beat"
        # Outside the axes, in axes fractions, so the readouts form one tidy column
        # instead of colliding with the tick labels whenever a gap comes out negative.
        right.annotate(
            f"{row.observed:+.4f} vs {row.null_median:+.4f}   p = {row.p_value:.3f}   "
            f"{verdict}",
            (1.014, y), xycoords=("axes fraction", "data"), va="center", ha="left",
            fontsize=8.2, color="#333333", annotation_clip=False,
        )
        labels.append(f"{textwrap.fill(row.claim, 44)}\nvs {row.null}\n{row.dataset}")
    right.axvline(0.0, color="#000000", lw=1.9, zorder=4)
    gaps = order.gap.to_numpy(dtype=float)
    pad = max(float(np.nanmax(np.abs(gaps))) * 0.12, 0.03)
    right.set_xlim(float(np.nanmin(gaps)) - pad, float(np.nanmax(gaps)) + pad)
    right.set_ylim(-0.7, len(order) - 0.3)
    right.set_yticks(range(len(order)))
    right.set_yticklabels(labels, fontsize=7.6)
    right.invert_yaxis()
    right.set_xlabel("observed score minus its null's median score\n"
                     "(same metric within each row; 0 = the null did as well)")
    right.set_title("B. The fitted basis fails its null; the given loading matrix does not")
    right.grid(axis="y", visible=False)
    right.legend(handles=[
        Line2D([], [], color=PALETTE[1], marker="X", ls="-", lw=2.2, markersize=10,
               markerfacecolor="#FFFFFF", markeredgewidth=2.0,
               label="DOES NOT beat its null"),
        Line2D([], [], color=PALETTE[2], marker="o", ls="-", lw=2.2, markersize=10,
               markerfacecolor=PALETTE[2], label="beats its null"),
    ], loc="upper right", fontsize=8.5)

    fig.suptitle("Transfer is interpolation across redundant coverage, not generalisation",
                 fontsize=12.5, fontweight="bold")
    _note(fig, (
        "A: outputs/external_redundancy.csv. B: outputs/external_nulls.csv, "
        "outputs/external_nulls_independent.csv and outputs/cross_family_transfer.csv "
        "(python3 scripts/run_external_validation.py; python3 scripts/run_cross_family.py). "
        "Both peer definitions in A are plotted because the panel's own target lists are "
        "incomplete -- see docs/EXTERNAL_VALIDATION.md.\n"
        "What to conclude, and it is negative: a module is recoverable when some other "
        "stressor in the panel drives it, and not otherwise, so what looks like "
        "generalisation to an unseen stressor is interpolation across redundant coverage. "
        "And the transfer score itself does not beat a random rotation of the same rank, "
        "in simulation or on real arrays, so the identity of the fitted latent basis is "
        "not what produces the score.\n"
        "The one claim that does beat its null is the readout through the literature "
        "loading matrix, which is given rather than fitted. That is the asymmetry: what "
        "is known works, what is estimated does not."
    ))
    return save_figure(fig, out_dir, stem)


# ---------------------------------------------------------------------------
# 08 - do the resolved calls survive the analyst's choices
# ---------------------------------------------------------------------------

ROBUSTNESS_SWEEPS = {
    "late_fraction": "late window: fraction of the trace treated as late",
    "autofluorescence_fraction":
        "autofluorescence assumed: share of the 0 mM reporter signal",
}
"""Knob columns that identify a sensitivity sweep, and how to label each on an axis.

Which sweep a frame is comes from the column it carries rather than from a label passed
alongside it, so a frame cannot be plotted under the wrong axis title.
"""


def interval_robustness_table(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Fold changes and their intervals across every analysis knob that was swept.

    Args:
        frames: Frames shaped like ``outputs/late_window_sensitivity.csv`` and
            ``outputs/autofluorescence_sensitivity.csv``. Each must carry exactly one
            column named in :data:`ROBUSTNESS_SWEEPS`.

    Returns:
        Long frame with ``sweep``, ``sweep_label``, ``knob_value``, ``construct``,
        ``dose_mM``, ``fold``, ``low``, ``high``, ``estimable``, ``excludes_unity``.

    Raises:
        ValueError: if a frame carries no knob column, more than one, or is missing the
            fold columns.
    """
    if not frames:
        raise ValueError("no sensitivity sweeps to tabulate")
    rows = []
    for frame in frames:
        knobs = [c for c in ROBUSTNESS_SWEEPS if c in frame.columns]
        if len(knobs) != 1:
            raise ValueError(
                f"a sweep frame must carry exactly one of "
                f"{sorted(ROBUSTNESS_SWEEPS)}; this one carries {knobs}"
            )
        needed = {"construct", "dose_mM", "fold", "low", "high", "estimable"}
        missing = needed - set(frame.columns)
        if missing:
            raise ValueError(f"sweep frame for {knobs[0]} is missing {sorted(missing)}")
        knob = knobs[0]
        part = frame.rename(columns={knob: "knob_value"}).copy()
        part["sweep"] = knob
        part["sweep_label"] = ROBUSTNESS_SWEEPS[knob]
        part["excludes_unity"] = (part.estimable.astype(bool)
                                  & ((part.low.astype(float) > 1.0)
                                     | (part.high.astype(float) < 1.0)))
        rows.append(part[["sweep", "sweep_label", "knob_value", "construct", "dose_mM",
                          "fold", "low", "high", "estimable", "excludes_unity"]])
    return pd.concat(rows, ignore_index=True)


def interval_robustness(
    table: pd.DataFrame, calls: pd.DataFrame, out_dir: pathlib.Path,
    stem: str = "fig04_interval_robustness",
) -> tuple[pathlib.Path, pathlib.Path]:
    """Figure: how far the resolved calls survive two arbitrary analysis choices.

    The calls in the headline figure clear 1.0 at one late window and one assumed
    autofluorescence. Neither was chosen by the biology, so the question this figure
    settles is whether the calls are properties of the plates or of those two choices. The
    quantity plotted is the *lower* bound of each interval, because that is what the claim
    rests on: a point estimate that stays above 1.0 while its lower bound crosses has not
    survived anything.

    Args:
        table: Output of :func:`interval_robustness_table`.
        calls: The rows of :func:`dose_response_table` to track -- in practice the ones
            whose ``verdict`` resolved above or below 1.0. Passed in rather than
            recomputed so the figure tracks exactly the calls the headline figure made.
        out_dir: Directory for the SVG and PNG.
        stem: Filename without extension.

    Returns:
        ``(svg_path, png_path)``.

    Raises:
        ValueError: if a frame is empty or lacks a column, or if no call in ``calls``
            appears in the sweeps.
    """
    missing = {"sweep", "sweep_label", "knob_value", "construct", "dose_mM", "low",
               "estimable"} - set(table.columns)
    if missing:
        raise ValueError(f"robustness table is missing {sorted(missing)}")
    if table.empty:
        raise ValueError("nothing to plot: the robustness table has no rows")
    missing = {"construct", "dose_mM"} - set(calls.columns)
    if missing:
        raise ValueError(f"call table is missing {sorted(missing)}")
    if calls.empty:
        raise ValueError(
            "no calls to track: dose_response_table resolved nothing away from 1.0, so "
            "there is no claim whose robustness this figure could test"
        )

    tracked = [(str(c), float(d)) for c, d in
               sorted(set(zip(calls.construct, calls.dose_mM.astype(float))))]
    doses = sorted({d for _, d in tracked})
    constructs = sorted({c for c, _ in tracked})
    if len(doses) > len(PALETTE):
        raise ValueError(
            f"{len(doses)} distinct doses need {len(doses)} colours and the palette has "
            f"{len(PALETTE)}; split the figure rather than cycling"
        )
    dose_colour = dict(zip(doses, series_colours(len(doses))))
    dose_marker = dict(zip(doses, ("o", "s", "^", "D", "v", "P", "X", "*")))
    construct_style = dict(zip(constructs, ("-", (0, (5, 2)), (0, (1, 1.6)),
                                            (0, (7, 2, 1, 2)))))

    sweeps = list(dict.fromkeys(table.sweep))
    # Shared y: the comparison the figure is for is between the two sweeps, and on
    # separate scales a knob that moves nothing looks exactly like one that moves a lot.
    fig, axes = _figure(1, len(sweeps), figsize=(6.4 * len(sweeps), 5.8), squeeze=False,
                        sharey=True)
    flat = [ax for row in axes for ax in row]
    drawn_any = False
    lost = []
    lost_calls = set()

    for ax, sweep in zip(flat, sweeps):
        panel = table[table.sweep == sweep]
        ax.axhline(1.0, color="#000000", lw=1.9, zorder=4)
        for construct, dose in tracked:
            series = panel[(panel.construct == construct)
                           & np.isclose(panel.dose_mM.astype(float), dose)]
            series = series[series.estimable.astype(bool)].sort_values("knob_value")
            if series.empty:
                continue
            drawn_any = True
            ax.plot(series.knob_value.astype(float), series.low.astype(float),
                    marker=dose_marker[dose], ls=construct_style.get(construct, "-"),
                    color=dose_colour[dose], markersize=7.5, lw=1.7,
                    markerfacecolor=dose_colour[dose], zorder=3)
            below = series[series.low.astype(float) <= 1.0]
            for row in below.itertuples():
                ax.plot([float(row.knob_value)], [float(row.low)], marker="o", ls="none",
                        markersize=14, markerfacecolor="none",
                        markeredgecolor="#000000", markeredgewidth=1.5, zorder=5)
                lost.append(f"{construct} at {dose:g} mM when "
                            f"{sweep.replace('_', ' ')} = {float(row.knob_value):g}")
                lost_calls.add((construct, dose))
        ax.set_xlabel(str(panel.sweep_label.iloc[0]))
        ax.set_title(sweep.replace("_", " "))
    flat[0].set_ylabel("lower bound of the 95% interval\non the corrected fold change")
    if not drawn_any:
        raise ValueError(
            "none of the tracked calls appears in the sweeps; the call table and the "
            "sensitivity tables do not describe the same constructs and doses"
        )

    handles = [
        Line2D([], [], color="#000000", lw=1.9, label="no change (1.0)"),
        Line2D([], [], color="#000000", marker="o", ls="none", markersize=12,
               markerfacecolor="none", markeredgewidth=1.5,
               label="ringed: the interval no longer clears 1.0"),
    ]
    handles += [
        Line2D([], [], color=dose_colour[d], marker=dose_marker[d], ls="none",
               markersize=8, label=f"{d:g} mM")
        for d in doses
    ]
    handles += [
        Line2D([], [], color="#555555", ls=construct_style.get(c, "-"), lw=1.7, label=c)
        for c in constructs
    ]
    flat[0].legend(handles=handles, loc="upper right", fontsize=8.0, ncols=2,
                   handlelength=2.4, labelspacing=0.32)

    kept = len(tracked) - len(lost_calls)
    fig.suptitle(
        f"Two analysis choices, {len(tracked)} resolved calls: {kept} keep their "
        f"interval everywhere swept", fontsize=12.5, fontweight="bold",
    )
    _note(fig, (
        "Sources: outputs/late_window_sensitivity.csv and "
        "outputs/autofluorescence_sensitivity.csv "
        "(python3 scripts/run_sensor_characterisation.py). Tracked calls are the rows "
        "the headline dose-response table resolved away from 1.0, so this figure and that "
        "one always describe the same claims.\n"
        + ("What to conclude: every tracked call keeps a lower bound above 1.0 across "
           "both sweeps, so none of them is an artefact of where the late window was cut "
           "or of how much autofluorescence was assumed."
           if not lost else
           "What to conclude: the calls are robust to the assumed autofluorescence over "
           "the whole range swept, and not uniformly robust to the late window. These "
           "lose their lower bound at an extreme of it: " + "; ".join(lost)
           + ". The point estimate stays above 1.0 throughout; the interval does not, "
             "and the interval is the claim.")
    ))
    return save_figure(fig, out_dir, stem)
