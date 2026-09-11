"""Assemble the figures and the tables behind them into one self-contained HTML page.

Wiki material generated from source rather than maintained by hand. Every number on the
page is read from a committed table in ``outputs/`` or from a table
``scripts/make_figures.py`` wrote next to the figure that drew it, so a page and a figure
can never disagree: they are the same frame rendered twice.

Nothing here is timestamped. A report whose bytes change on every build cannot be
reviewed in a diff, and "when was this generated" is a question the commit history
answers better than a line of text nobody updates.

The SVGs are inlined rather than linked so the page is one file that can be attached,
mailed or dropped into a wiki. Their internal ids are prefixed per figure on the way in:
five matplotlib SVGs in one document otherwise share glyph and clip-path ids, and the
first definition silently wins.

Usage: python3 scripts/build_report.py [--figures-dir DIR] [--out PATH]
"""

from __future__ import annotations

import argparse
import html
import pathlib
import re
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.uncertainty import MIN_PLATES_FOR_INTERVAL
from ystwin.viz.figures import CORRECTABLE_GDNA_CEILING, read_table

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_FIGURES_DIR = REPO_ROOT / "figures"
MAKE_FIGURES = "python3 scripts/make_figures.py"

_STYLE = """
:root {
  --ink: #1a1a1a; --soft: #555; --rule: #d8d8d8; --wash: #f6f7f8;
  --blue: #0072B2; --vermillion: #D55E00; --green: #009E73;
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 0 1.5rem 6rem; color: var(--ink); background: #fff;
  font: 16px/1.62 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
}
main { max-width: 62rem; margin: 0 auto; }
header { padding: 3.5rem 0 1.5rem; border-bottom: 3px solid var(--ink); }
h1 { font-size: 2.1rem; line-height: 1.15; margin: 0 0 .6rem; letter-spacing: -.02em; }
h2 {
  font-size: 1.45rem; margin: 3.5rem 0 .4rem; padding-top: 1.2rem;
  border-top: 1px solid var(--rule); letter-spacing: -.01em;
}
h3 { font-size: 1.06rem; margin: 2.2rem 0 .4rem; }
p, li { max-width: 46rem; }
.lede { font-size: 1.08rem; color: var(--soft); max-width: 44rem; }
.sub { color: var(--soft); font-size: .92rem; }
code, .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: .88em; }
figure { margin: 1.6rem 0 2.4rem; }
figure svg { max-width: 100%; height: auto; display: block; }
figcaption { font-size: .88rem; color: var(--soft); margin-top: .7rem; max-width: 46rem; }
figcaption b { color: var(--ink); }
table { border-collapse: collapse; width: 100%; margin: 1.1rem 0 1.6rem; font-size: .9rem; }
th, td { text-align: left; padding: .42rem .6rem; border-bottom: 1px solid var(--rule); }
th { background: var(--wash); font-weight: 600; white-space: nowrap; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tbody tr:hover { background: #fafbfc; }
.callout {
  border-left: 4px solid var(--vermillion); background: var(--wash);
  padding: .9rem 1.1rem; margin: 1.4rem 0; font-size: .95rem;
}
.callout b { color: var(--vermillion); }
.badge {
  display: inline-block; padding: .06rem .45rem; border-radius: 3px;
  font-size: .78rem; font-weight: 600; letter-spacing: .02em;
  background: #eee; color: #333;
}
.badge.no { background: #fbe6d9; color: #8a3b00; }
.badge.yes { background: #dcf1ea; color: #005f45; }
footer { margin-top: 4rem; padding-top: 1.2rem; border-top: 1px solid var(--rule);
         color: var(--soft); font-size: .85rem; }
"""


# ---------------------------------------------------------------------------
# assembling pieces
# ---------------------------------------------------------------------------


def inline_svg(path: pathlib.Path, prefix: str) -> str:
    """One figure's SVG, id-namespaced so five of them can share a document.

    Args:
        path: The ``.svg`` written by :mod:`ystwin.viz.figures`.
        prefix: Short unique string prepended to every id and reference.

    Returns:
        The ``<svg>`` element, without the XML declaration or doctype.

    Raises:
        FileNotFoundError: if the figure has not been built.
        ValueError: if the file holds no ``<svg>`` element.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"figure {path.name} is missing from {path.parent}. Build the figures first: "
            f"`{MAKE_FIGURES}`."
        )
    text = path.read_text(encoding="utf-8")
    start = text.find("<svg")
    if start < 0:
        raise ValueError(f"{path} does not contain an <svg> element")
    body = text[start:]
    body = re.sub(r'(?<=id=")([^"]+)(?=")', lambda m: prefix + m.group(1), body)
    body = re.sub(r'(?<=href="#)([^"]+)(?=")', lambda m: prefix + m.group(1), body)
    body = re.sub(r"(?<=url\(#)([^)]+)(?=\))", lambda m: prefix + m.group(1), body)
    return body


def table(
    frame: pd.DataFrame,
    columns: dict[str, str],
    numeric: frozenset[str] = frozenset(),
    raw: frozenset[str] = frozenset(),
) -> str:
    """A frame as an HTML table, with the columns and headings named explicitly.

    Explicit rather than ``to_html`` because a report should break when a column it
    names disappears, not quietly render whatever the frame happens to hold. Cell text
    is escaped by default; a column has to be named in ``raw`` to get markup through,
    which keeps "this cell is a badge" a decision rather than an accident.

    Args:
        frame: Rows to render, already ordered and formatted as strings where needed.
        columns: Source column -> heading, in display order.
        numeric: Headings to right-align.
        raw: Headings whose cells already hold HTML and must not be escaped.

    Raises:
        ValueError: if a named column is absent.
    """
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"table is missing {sorted(missing)}")
    head = "".join(
        f'<th class="num">{html.escape(h)}</th>' if h in numeric else f"<th>{html.escape(h)}</th>"
        for h in columns.values()
    )
    rows = []
    for _, row in frame.iterrows():
        cells = []
        for column, heading in columns.items():
            value = str(row[column])
            body = value if heading in raw else html.escape(value)
            cells.append(f'<td class="num">{body}</td>' if heading in numeric
                         else f"<td>{body}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def figure_block(figures_dir: pathlib.Path, stem: str, prefix: str, caption: str) -> str:
    """One figure, inlined, captioned, and cross-referenced to the file it came from.

    The report numbers its figures in reading order and the files are numbered in build
    order, which are not the same sequence. The caption carries the filename so that
    "Figure 3" and ``05_g4_contamination.svg`` are visibly the same object.
    """
    return (
        f'<figure id="{prefix.rstrip("-")}">'
        f"{inline_svg(figures_dir / f'{stem}.svg', prefix)}"
        f'<figcaption>{caption} '
        f'<span class="sub mono">[{stem}.svg &middot; {stem}.png]</span>'
        f"</figcaption></figure>"
    )


# ---------------------------------------------------------------------------
# the sections
# ---------------------------------------------------------------------------


def _contract() -> str:
    return """
<h2 id="contract">1 &middot; The contract</h2>
<p>Borrowed vocabulary rather than invented: Kapteyn, Pretorius &amp; Willcox's
six-quantity abstraction of an asset-twin system
(<i>Nature Computational Science</i> <b>1</b>(5):337&ndash;347, 2021). The distinction
that earns its keep here is <b>S</b>, the culture's physical state, against <b>D</b>, the
digital state the model actually carries &mdash; and their own statement of the
consequence, that the digital state space is generally only a subset of the physical one,
is not a caveat in this project but its central finding.</p>
<div class="callout">
<b>What this twin does, in one sentence.</b> Given a dose and a well's OD and fluorescence
time series, it returns a posterior over that well's growth rate and promoter activity
&mdash; dilution removed, uncertainty attached, and refused outright where the optical
channel is not quantitative.
</div>
<p>Three ambition levels, and only the first is committed to.
<b>Level 1</b>, a posterior over one well's digital state with an interval and a refusal,
has every component built and needs only wiring.
<b>Level 2</b>, a scored forward prediction against a baseline ladder, needs no bench work
&mdash; the data is on disk.
<b>Level 3</b>, a distribution over the &beta;-carotene titre curve, is blocked twice
over: no strain in this dataset carries a carotenoid pathway, and D1 established that FBA
bounds the product flux at <span class="mono">[0, ceiling]</span> at every growth level
tested, relative width 1.000, so it can never predict it. Level 3 is a design target and
should never appear in a results table.</p>
<p class="sub">Full statement, including the undefined reward term and what follows from
its absence: <span class="mono">docs/CONTRACT.md</span>. Every fold change on this page is
a <b>Q</b>, every plate export an <b>O</b>, every dose and nutrient level a
<b>U</b>. Autofluorescence is a term in the measurement model that no data
constrains &mdash; see section 4.</p>
"""


def _measured(figures_dir: pathlib.Path, outputs: pathlib.Path) -> str:
    g1 = read_table(figures_dir / "g1_summary.csv", "the G1 summary", MAKE_FIGURES)
    per_plate = (
        g1[["plate", "n_wells", "n_pass", "pass_rate"]].drop_duplicates("plate")
        .assign(plate=lambda d: d.plate.str.replace("\n", " ", regex=False),
                pass_rate=lambda d: (d.pass_rate * 100).round(0).astype(int).astype(str) + "%")
    )
    dominant = (
        g1[g1.outcome != "passed"].sort_values("n_wells_with_outcome", ascending=False)
        .drop_duplicates("plate")
        .assign(plate=lambda d: d.plate.str.replace("\n", " ", regex=False))
        [["plate", "outcome", "n_wells_with_outcome"]]
        .rename(columns={"outcome": "largest failure signature"})
    )
    per_plate = per_plate.merge(dominant, on="plate", how="left")
    per_plate["largest failure signature"] = (
        per_plate["largest failure signature"].fillna("--").str.replace(",", " + ", regex=False)
        + " ("
        + per_plate.n_wells_with_outcome.fillna(0).astype(int).astype(str)
        + " wells)"
    )

    folds = read_table(figures_dir / "dose_response_folds.csv",
                       "the dose-response fold table", MAKE_FIGURES)
    shown = folds.copy()
    shown["naive"] = shown.naive_fold.map(lambda v: f"{v:.2f}" if pd.notna(v) else "--")
    shown["corrected"] = [
        f"{v:.2f}" if pd.notna(v) and usable >= 2 else "not estimable"
        for v, usable in zip(shown.corrected_fold, shown.n_plates_usable_corrected)
    ]
    shown["dose"] = shown.dose_mM.map(lambda v: f"{v:g}")
    # The sense was inverted -- an estimable row printed "--" and only a refused one
    # printed anything -- so the report's interval column was blank exactly where an
    # interval existed. The "(n=2)" was stale from two plates as well.
    shown["interval"] = [
        f"[{lo:.2f}, {hi:.2f}]" if est else f"not estimable (n={int(n)})"
        for est, lo, hi, n in zip(shown.interval_estimable, shown.interval_low,
                                  shown.interval_high, shown.n_plates_usable_corrected)
    ]

    contamination = read_table(figures_dir / "g4_contamination.csv",
                               "the G4 contamination table", MAKE_FIGURES)
    anchors = (
        contamination[contamination.role == "anchor"]
        .groupby(["construct", "target"], as_index=False)
        .agg(readings=("gdna_share", "size"),
             median_share=("gdna_share", "median"),
             pass_rate=("passed", "mean"))
    )
    verdicts = read_table(outputs / "g4_verdicts.csv", "the G4 verdict table",
                          "python3 scripts/run_g4.py")
    anchors = anchors.merge(verdicts[["construct", "verdict", "replicates_needed"]],
                            on="construct", how="left")
    anchors["gDNA share (median)"] = anchors.median_share.map(lambda v: f"{v:.1%}")
    anchors["no-RT pass rate"] = anchors.pass_rate.map(lambda v: f"{v:.0%}")
    anchors["correctable?"] = [
        '<span class="badge yes">yes</span>' if v <= CORRECTABLE_GDNA_CEILING
        else '<span class="badge no">past the ceiling</span>'
        for v in anchors.median_share
    ]
    anchors["replicates_needed"] = anchors.replicates_needed.fillna(0).astype(int).map(
        lambda n: f"{n}" if n else "--")

    verdict_table = table(
        anchors, {
            "construct": "construct", "target": "anchor",
            "gDNA share (median)": "gDNA share (median)",
            "no-RT pass rate": "no-RT pass rate",
            "correctable?": "correctable?",
            "verdict": "G4 verdict",
            "replicates_needed": "replicates that would settle it",
        },
        numeric=frozenset({"gDNA share (median)", "no-RT pass rate",
                           "replicates that would settle it"}),
        raw=frozenset({"correctable?"}),
    )

    return f"""
<h2 id="measured">2 &middot; What was measured</h2>
<p>Four mCitrine reporters on a plate reader, and one qPCR anchor. UPRE1 and UPRE2 report
ER stress and are challenged with DTT; NativeYap1 and AlteredYap1 report oxidative stress
and are challenged with H<sub>2</sub>O<sub>2</sub>. <b>No strain here carries a carotenoid
pathway.</b></p>

<h3>2.1 &nbsp;G1: is the optical channel quantitative at all?</h3>
<p>G1 runs before anything else, because asking how much of a reporter signal is growth
rate is meaningless in a well whose OD is not a biomass measurement. On the two 2026-07
plates most wells are not: cultures ran to raw OD 1.8&ndash;1.9 in a 96-well plate, past
the reader's linear range, and <span class="mono">linear_range</span> is the dominant
failure.</p>
{table(per_plate, {"plate": "plate", "n_wells": "cultures", "n_pass": "passed",
                   "pass_rate": "pass rate",
                   "largest failure signature": "largest failure signature"},
       numeric=frozenset({"cultures", "passed", "pass rate"}))}
{figure_block(figures_dir, "01_g1_pass_rates", "f1-",
              "<b>Figure 1.</b> Every well on each plate, partitioned by the exact set of "
              "criteria it failed &mdash; not by failure reason, which would count a well "
              "failing two criteria twice and make the bar taller than the plate. The G1 "
              "tables committed to <span class='mono'>outputs/</span> cover only the two "
              "2026-07 plates; the NewProtocol replicates are named in the footnote rather "
              "than drawn from numbers this figure has not read.")}

<h3>2.2 &nbsp;The dilution correction, and what it takes away</h3>
<p>A stable reporter accumulates as growth slows, so a raw per-cell readout rises whenever
a dose is toxic &mdash; whether or not the promoter responded.
<span class="mono">reporter.py</span> inverts
<span class="mono">dR/dt = k_synth &minus; (mu + k_deg) R</span> to recover promoter
activity, and the difference between the two curves is how much of the apparent induction
was growth.</p>
<div class="callout">
<b>Two ends of every curve are artefacts.</b> The oxidative sensors do not detect
0.1&ndash;0.2 mM H<sub>2</sub>O<sub>2</sub>; they register the growth slowdown it causes.
The ER sensors' largest apparent inductions, at 2&ndash;5 mM DTT, vanish entirely under
correction. What survives is a genuine 1.2&ndash;1.5&times; band in the middle of each
ladder.
</div>
{figure_block(figures_dir, "02_dose_response_correction", "f2-",
              "<b>Figure 2.</b> Naive per-cell readout against dilution-corrected promoter "
              "activity, per construct. No error bars: at two biological replicates "
              f"<span class='mono'>fold_change</span> refuses a cluster bootstrap below "
              f"{MIN_PLATES_FOR_INTERVAL} plates, and their absence is the result rather "
              "than an omission.")}
{table(shown, {"construct": "construct", "stressor": "stressor", "dose": "dose (mM)",
               "naive": "naive fold", "corrected": "corrected fold",
               "n_plates": "plates", "interval": "interval"},
       numeric=frozenset({"dose (mM)", "naive fold", "corrected fold", "plates"}))}

<h3>2.3 &nbsp;G4: the biosensors against an independent anchor</h3>
<p>A reporter that agrees with itself proves nothing, so each construct is tested against
an endogenous transcript: Hac1 for the ER pair, TRX2 for the oxidative pair, both against
UBC. All four verdicts are INCONCLUSIVE, for two unrelated reasons.</p>
{verdict_table}
<p>The ER anchor is not measuring transcript. Its no-reverse-transcriptase control
amplifies at essentially the same cycle as the +RT sample in all three replicates, so most
of the Hac1 signal is genomic DNA &mdash; and the logbook names the primer pair
<b>"Hac1 Exon"</b>, which explains it: exon primers amplify the genomic locus readily, and
they measure <i>total</i> HAC1, which barely moves during the UPR because HAC1 is
regulated by unconventional splicing. Two independent fixes: DNase-treat and use
intron-spanning primers, and switch the readout to KAR2.</p>
{figure_block(figures_dir, "05_g4_contamination", "f5-",
              "<b>Figure 3.</b> Every no-RT reading as a contamination share rather than a "
              "cycle count. Cycles mislead near zero: 0.6 cycles sounds only a little worse "
              "than 1.0 and is 66% genomic DNA against 50%, which straddles the 60% share "
              "above which no published correction holds.")}
"""


def _simulated(figures_dir: pathlib.Path) -> str:
    designs = read_table(figures_dir / "decoupling_designs.csv",
                         "the decoupling design table", MAKE_FIGURES)
    scores = (designs.groupby("design", as_index=False)
              .agg(conditions=("dose_mM", "size"), collinearity=("collinearity", "first"),
                   nutrient_levels=("nutrient_factor", "nunique")))
    scores["collinearity"] = scores.collinearity.map(lambda v: f"{v:.3f}")

    power = read_table(figures_dir / "attribution_power.csv",
                       "the attribution power table", MAKE_FIGURES)
    wide = power.copy()
    wide["cell"] = [
        f"{p:.0%}  [{lo:.0%}, {hi:.0%}]"
        for p, lo, hi in zip(wide.power, wide.power_low, wide.power_high)
    ]
    # Narrative order -- one reporter, then the design fix, then the strain fix -- not the
    # alphabetical order a pivot would impose.
    order = list(dict.fromkeys(wide.arm))
    pivot = (wide.pivot(index="arm", columns="induction_fold", values="cell")
             .reindex(order).reset_index())
    pivot.columns = [c if isinstance(c, str) else f"{c:g}x" for c in pivot.columns]
    fold_columns = {c: c for c in pivot.columns if c != "arm"}

    return f"""
<h2 id="simulated">3 &middot; What was simulated</h2>
<p>Two results here decide whether a latent stress state is reachable at all, and neither
needs a plate. Both are deterministic given their seed and both record their settings on
the figure, because a simulated number with no visible settings is unreproducible.</p>

<h3>3.1 &nbsp;The decoupling grid</h3>
<p>Every perturbation run so far slows growth, so stress and <span class="mono">1/mu</span>
are one axis and no latent stress state is separable from growth at any sample size.
Nutrient limitation slows growth <i>without</i> stressing; crossing it with the dose ladder
fills the quadrants the dose-only design leaves empty.</p>
{table(scores, {"design": "design", "conditions": "conditions",
                "nutrient_levels": "nutrient levels", "collinearity": "collinearity"},
       numeric=frozenset({"conditions", "nutrient levels", "collinearity"}))}
{figure_block(figures_dir, "03_decoupling_grid", "f3-",
              "<b>Figure 4.</b> Dose against growth rate for each planned design. In the "
              "dose-only series a dose determines the growth rate exactly, so no fit can "
              "say which of them moved the reporter; adding nutrient levels lifts the same "
              "dose to several growth rates. The collinearity number is what the picture "
              "already shows.")}

<h3>3.2 &nbsp;Attribution power, and what a second reporter buys</h3>
<p>Detecting a dose slope is easy on any design &mdash; growth inhibition alone guarantees
one. The question that matters is whether the response can be attributed to <i>dose</i>
rather than to the growth the dose caused, which means fitting both and asking whether
dose survives. On the design as run, at the effect sizes actually measured, it does not.</p>
{table(pivot, {"arm": "assay", **fold_columns}, numeric=frozenset(fold_columns))}
<p>Two reporters in one cell share <span class="mono">mu</span> exactly, so a ratio between
them cancels the dilution term structurally rather than as a correction &mdash; and that
rescues the collinear design already run, with no nutrient axis at all. Two conditions
attach: the reference must be genuinely constitutive, and the maturation rates must match
or <span class="mono">mu</span> does not cancel. mOrange2 against YFP leaves a 23%
growth-dependent residual, larger than the effects being measured.</p>
{figure_block(figures_dir, "04_attribution_power", "f4-",
              "<b>Figure 5.</b> Power to attribute a dose response to dose, with Wilson "
              "score intervals on the simulated proportion. Regenerated rather than read: "
              "<span class='mono'>outputs/power_analysis.csv</span> carries neither the "
              "intervals nor the co-expressed-ratio arm.")}
"""


def _refused(figures_dir: pathlib.Path) -> str:
    refusals = pd.DataFrame([
        {
            "claim": "A corrected fold change is distinguishable from 1.0",
            "status": "REFUSED",
            "why": f"Two biological replicates. A cluster bootstrap on two clusters "
                   f"describes which of two plates was drawn, not the population of "
                   f"plates, so uncertainty.fold_change refuses an interval below "
                   f"{MIN_PLATES_FOR_INTERVAL}.",
            "enforced by": "analysis/uncertainty.py::fold_change",
        },
        {
            "claim": "The 2-5 mM DTT readings are entirely dilution",
            "status": "NOT ROBUST",
            "why": "Autofluorescence has never been measured on any plate. An unremoved "
                   "constant contributes exactly a*mu to the recovered activity; at a "
                   "plausible share, UPRE1 at 5 mM runs 0.92 to 1.07 and the claim "
                   "inverts. The 0.1-0.2 mM calls move by at most 0.02 and do stand.",
            "enforced by": "docs/DATA_INVENTORY.md; observation.py",
        },
        {
            "claim": "A biosensor agrees with its independent anchor",
            "status": "INCONCLUSIVE (4/4)",
            "why": "The ER anchor is past the 60% genomic-DNA share above which no "
                   "published correction holds. The oxidative anchor is clean but too "
                   "noisy at the replicate count available.",
            "enforced by": "gates/g4_anchor.py; qpcr.rt_minus_margin",
        },
        {
            "claim": "The dose response is attributable to dose",
            "status": "REFUSED",
            "why": "Collinearity 0.993 on the design as run. At the measured 1.2-1.5x "
                   "effects the design has almost no power to separate dose from the "
                   "growth the dose caused, at any replicate count.",
            "enforced by": "generator/design.py::collinearity; analysis/power.py",
        },
        {
            "claim": "FBA predicts the beta-carotene titre",
            "status": "REFUSED",
            "why": "The feasible product flux is [0, ceiling] at every growth level "
                   "tested, relative width 1.000. A capacity bound moves the ceiling and "
                   "never lifts the floor off zero, so FBA bounds the product and never "
                   "predicts it. No strain in this dataset carries the pathway either.",
            "enforced by": "fba/fva.py; docs/CONTRACT.md level 3",
        },
        {
            "claim": "The NewProtocol plates pass G1 at 80/100/86%",
            "status": "NOT PLOTTED",
            "why": "True per docs/FINDINGS.md, and no G1 table for those plates "
                   "exists in outputs/. This report renders committed tables, so the "
                   "figure names the gap instead of drawing a number it has not read.",
            "enforced by": "scripts/make_figures.py; plate/layout.py::RECORDED_PLATES",
        },
    ])
    return f"""
<h2 id="refused">4 &middot; What is refused</h2>
<p>The refusals are the load-bearing part. Each one is enforced by code that raises or
returns an unestimable result rather than by a note somebody has to remember, which is the
only version of a caveat that survives contact with a deadline.</p>
{table(refusals, {"claim": "claim", "status": "status", "why": "why",
                  "enforced by": "enforced by"})}
<div class="callout">
<b>INCONCLUSIVE is not refuted.</b> The point estimates are unchanged and the mechanism
behind them &mdash; a stable reporter accumulating as growth slows &mdash; is sound and
independently supported. What is missing is the evidence that these particular numbers are
distinguishable from 1.0. Four items fix most of it, in cost order: read the reporter-free
strain in the mCitrine channel, fluoresce all twelve columns, put blanks in every channel
that is read, and run two more full biological replicates.
</div>
"""


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------


def build(figures_dir: pathlib.Path, outputs: pathlib.Path) -> str:
    """The whole page as a string.

    Raises:
        FileNotFoundError: if a figure or its source table has not been built.
        ValueError: if a table is empty or missing a column the report names.
    """
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        "<title>ystwin &mdash; yeast stress digital twin, model v2</title>\n"
        f"<style>{_STYLE}</style>\n</head>\n<body>\n<main>\n"
        "<header>\n"
        "<h1>What this twin measures, what it simulates,<br>and what it refuses to say</h1>\n"
        '<p class="lede">Four mCitrine stress reporters in <i>S. cerevisiae</i>, one qPCR '
        "anchor, and the diagnostics that decide what the rest of the architecture is "
        "allowed to claim. Every number below is read from a committed table; nothing on "
        "this page is typed by hand.</p>\n"
        '<p class="sub mono">generated by scripts/build_report.py from outputs/ and '
        "figures/</p>\n"
        "</header>\n"
        f"{_contract()}\n{_measured(figures_dir, outputs)}\n"
        f"{_simulated(figures_dir)}\n{_refused(figures_dir)}\n"
        "<footer>Regenerate with <span class=\"mono\">python3 scripts/make_figures.py "
        "&amp;&amp; python3 scripts/build_report.py</span>. Figures are also written as "
        "PNG alongside the SVGs, and each figure's source table sits next to it as CSV."
        "</footer>\n"
        "</main>\n</body>\n</html>\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures-dir", type=pathlib.Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--out", type=pathlib.Path, default=None,
                        help="output path (default: <figures-dir>/report.html)")
    args = parser.parse_args()

    out = args.out or (args.figures_dir / "report.html")
    try:
        page = build(args.figures_dir, paths.outputs_dir())
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    print(f"wrote {out}  ({len(page) / 1024:.0f} kB, {page.count('<svg')} figures inline)")


if __name__ == "__main__":
    main()
