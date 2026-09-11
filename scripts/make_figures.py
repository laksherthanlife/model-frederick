"""Build every figure in the project from the committed tables, and nothing else.

The figures are a rendering of ``outputs/``, not a second analysis of it. This script is
the only place that knows which table feeds which figure, so a table that moves or
disappears fails here, loudly, with the name of the script that writes it -- rather than
producing a figure with no series in it, which renders perfectly and says nothing.

Two figures are simulated rather than measured. The decoupling grid comes from
``generator/design.py`` and costs nothing; the attribution power comes from
``analysis/power.py`` and costs about a minute, because ``outputs/power_analysis.csv``
carries neither the Monte Carlo intervals nor the co-expressed-ratio arm. Both write the
frame they simulated next to the figure so that ``scripts/build_report.py`` reports the
same numbers the figure drew.

There are two sets, and they are not the same set at two sizes. The **report set**
(``figures/``, numbered ``01``-``05``) is what ``build_report.py`` inlines, one figure per
section, including the ones that exist to show a refusal. The **story set**
(``outputs/``, numbered ``fig01``-``fig04``) is the four figures the project's argument
actually rests on -- the dilution correction and its intervals, the held-out prediction
and where it fails, what transfer does not recover, and whether the resolved calls survive
the analysis choices behind them. Every one of the four traces to a committed CSV, and
``outputs/figures_manifest.csv`` records which, alongside the command that regenerates it,
so a reader who distrusts a figure can find the number behind it without reading this
script.

Usage: python3 scripts/make_figures.py [--figures-dir DIR] [--story-dir DIR]
                                       [--power-simulations N]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import textwrap
from collections.abc import Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.plate.layout import RECORDED_PLATES
from ystwin.viz import figures

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_FIGURES_DIR = REPO_ROOT / "figures"


def _default_figures_dir() -> pathlib.Path:
    """``figures/`` in the checkout, or ``figures/`` beside a redirected ``outputs/``.

    ``--figures-dir`` defaulted to the tracked ``figures/`` whatever the environment said,
    so ``YSTWIN_OUTPUTS=/tmp/check python3 scripts/make_figures.py`` -- the dry run
    docs/REPRODUCING.md prescribes for exactly this -- rewrote four tracked SVGs and two
    tracked tables while presenting itself as a check. Read off ``paths.outputs_dir()``
    rather than the variable, so there is one answer to where a run is allowed to write.
    """
    outputs = paths.outputs_dir()
    if outputs == REPO_ROOT / "outputs":
        return DEFAULT_FIGURES_DIR
    return outputs / "figures"


def _plate_label(tag: str) -> str:
    """A readable plate name from a ``g1_<tag>.csv`` stem.

    The tag is the export's filename with spaces and ampersands substituted, so it
    carries the date the plate was run in its first eight characters. That date is the
    thing a reader needs on an axis; the rest is the lab's own file naming, and it is
    wrapped narrow because unwrapped it runs into the neighbouring plate's label -- two
    plate names overlapping on an axis is worse than either being on three lines.
    """
    if len(tag) >= 8 and tag[:8].isdigit():
        date = f"{tag[:4]}-{tag[4:6]}-{tag[6:8]}"
        rest = tag[8:].strip("_").replace("_", " ").replace("(RAW)", "").strip()
        wrapped = textwrap.fill(rest, 16, break_long_words=False)
        return f"{date}\n{wrapped}" if rest else date
    return tag


def _g1_tables(outputs: pathlib.Path) -> dict[str, pd.DataFrame]:
    """Every committed G1 export, keyed by a readable plate label."""
    found = sorted(outputs.glob("g1_*.csv"))
    if not found:
        raise SystemExit(
            f"no G1 tables under {outputs}. Figure 01 is the gate report; without one "
            f"there is nothing to draw. Regenerate with `python3 scripts/run_gates.py`."
        )
    return {
        _plate_label(path.stem[len("g1_"):]):
            figures.read_table(path, f"G1 gate report {path.name}",
                               "python3 scripts/run_gates.py")
        for path in found
    }


def _plates_without_a_table(outputs: pathlib.Path) -> list[str]:
    """Recorded plates that no committed G1 export covers.

    ``plate/layout.py::RECORDED_PLATES`` is the logbook-confirmed list of plates that
    were actually run, so it -- and not a list written into this script -- decides what
    is missing.
    """
    have = {path.stem[len("g1_"):][:8] for path in outputs.glob("g1_*.csv")}
    return [f"{d[:4]}-{d[4:6]}-{d[6:8]}" for d in sorted(RECORDED_PLATES) if d not in have]


def _pooled_to_palette(summary: pd.DataFrame) -> pd.DataFrame:
    """Pool the rarest failure signatures until the palette can separate what is left.

    ``figures.g1_pass_rates`` refuses more outcomes than the colour-blind-safe palette
    has, and its message tells the caller to group the categories rather than cycle
    colours -- because a cycled colour is two failure modes wearing one colour, which is
    worse than no figure. This is that grouping, and it belongs here rather than in the
    drawing code: which signatures are worth their own colour is a property of the plates
    that happen to be in ``outputs/``, not of the figure.

    The pooling keeps the partition exact. Every well still appears in exactly one
    outcome, so the stacked bars still reach the plate's well count, and the pooled
    outcome names how many signatures went into it so the figure does not imply there
    was only one.

    Args:
        summary: Output of :func:`figures.g1_summary`.

    Returns:
        ``summary`` unchanged when it already fits the palette; otherwise the same frame
        with its rarest failure signatures collapsed into one outcome.
    """
    budget = len(figures.PALETTE)
    signatures = (summary[summary.outcome != "passed"]
                  .groupby("outcome").n_wells_with_outcome.sum()
                  .sort_values(ascending=False))
    has_passed = bool((summary.outcome == "passed").any())
    room = budget - (1 if has_passed else 0)
    if len(signatures) <= room:
        return summary

    keep = set(signatures.index[:room - 1])
    pooled_count = len(signatures) - len(keep)
    label = f"other ({pooled_count} rarer signatures)"
    rows = []
    for plate, group in summary.groupby("plate", sort=False):
        kept = group[(group.outcome == "passed") | group.outcome.isin(keep)]
        rows.extend(kept.to_dict("records"))
        rest = group[(group.outcome != "passed") & ~group.outcome.isin(keep)]
        if rest.empty:
            continue
        first = rest.iloc[0]
        rows.append({
            "plate": plate, "outcome": label,
            "n_wells_with_outcome": int(rest.n_wells_with_outcome.sum()),
            "n_wells": int(first.n_wells), "n_pass": int(first.n_pass),
            "pass_rate": float(first.pass_rate),
            "n_criteria": int(rest.n_criteria.max()),
        })
    return pd.DataFrame(rows)


_CI_CLEAR = "CI excludes 1.0"
_CI_SPANS = "CI spans 1.0"
_CI_REFUSED = "CI REFUSED (finite fold)"
_FOLD_MISSING = "CI REFUSED (fold missing)"


def _valid_bounds(estimable, low, high):
    return estimable.eq(True) & np.isfinite(low) & np.isfinite(high) & (low <= high)


def _fold_status(folds: pd.DataFrame) -> pd.Series:
    """Presentation categories from interval availability, never from a finite point alone."""
    finite = np.isfinite(folds.corrected_fold)
    valid = finite & _valid_bounds(folds.interval_estimable,
                                   folds.interval_low, folds.interval_high)
    status = pd.Series(_FOLD_MISSING, index=folds.index)
    status.loc[finite] = _CI_REFUSED
    status.loc[valid] = _CI_SPANS
    status.loc[valid & ((folds.interval_low > 1) | (folds.interval_high < 1))] = _CI_CLEAR
    return status


def _dose_coverage(sensor: pd.DataFrame) -> pd.DataFrame:
    """Count recorded plate/well identities, separately for every construct and dose."""
    keys = ["construct", "stressor", "dose_mM"]
    wells = sensor.groupby([*keys, "plate"]).well.nunique()
    return wells.groupby(level=keys).agg(
        plates="size", wells="sum", min_wells="min", max_wells="max").reset_index()


def _count_range(values) -> str:
    lo, hi = int(min(values)), int(max(values))
    return str(lo) if lo == hi else f"{lo}-{hi}"


def _caption(fig, text: str) -> None:
    """Reserve the measured caption height rather than estimating its line height."""
    figures._note(fig, text)
    caption = fig.texts[-1]
    caption.set_position((.006, .02))
    caption.set_verticalalignment("bottom")
    height = caption.get_window_extent(fig.canvas.get_renderer()).height / fig.bbox.height
    reserved = height + .04
    fig.get_layout_engine().set(rect=(.004, reserved, .992, .996 - reserved))


def _dose_response_figure(folds, sensor, out_dir, stem="02_dose_response_correction"):
    """Render the unchanged fold table with coverage and all interval refusal states."""
    table = folds.assign(status=_fold_status(folds)).merge(
        _dose_coverage(sensor), on=["construct", "stressor", "dose_mM"],
        how="left", validate="one_to_one")
    if table.empty or table.plates.isna().any():
        raise ValueError("fold conditions need recorded plate/well coverage")
    valid = table.status.isin((_CI_CLEAR, _CI_SPANS))
    groups = list(table.groupby(["stressor", "construct"], sort=True))
    fig, axes = figures._figure((len(groups) + 1) // 2, 2, squeeze=False,
                                figsize=(10.8, 4.5 * ((len(groups) + 1) // 2)))
    spans = np.concatenate([table.naive_fold, table.corrected_fold,
                            table.loc[valid, "interval_low"], table.loc[valid, "interval_high"]])
    finite = spans[np.isfinite(spans)]
    bottom = min(-0.08, float(finite.min()) - 0.08) if finite.size else -0.08
    top = max(1.0, float(finite.max())) * 1.32 if finite.size else 2.0
    for ax, ((stressor, construct), panel) in zip(axes.flat, groups):
        panel = panel.sort_values("dose_mM")
        has_ci = panel.status.isin((_CI_CLEAR, _CI_SPANS))
        ax.axhline(1, color="#555555", ls="--", lw=1, label="reference (1.0)")
        ax.plot(panel.dose_mM, panel.naive_fold, "o--", color=figures.PALETTE[1],
                mfc="none", lw=1.5, label="naive RFU/OD (point only)")
        ax.plot(panel.dose_mM, panel.corrected_fold.where(has_ci),
                color=figures.PALETTE[0], lw=1.5)
        for row in panel[has_ci].itertuples():
            ax.vlines(row.dose_mM, row.interval_low, row.interval_high,
                      color=figures.PALETTE[0], lw=1.3)
            ax.plot([row.dose_mM] * 2, [row.interval_low, row.interval_high],
                    marker="_", ls="none", color=figures.PALETTE[0])
        for status in (_CI_CLEAR, _CI_SPANS, _CI_REFUSED, _FOLD_MISSING):
            part = panel[panel.status == status]
            if status == _FOLD_MISSING:
                ax.scatter(part.dose_mM, np.full(len(part), 0.045), marker="X", s=45,
                           color="#666666", transform=ax.get_xaxis_transform(),
                           label=status, zorder=5)
                for dose in part.dose_mM:
                    ax.axvline(dose, color="#777777", alpha=0.15, lw=7)
                    ax.text(dose, 0.09, "fold missing", rotation=90, ha="center",
                            va="bottom", transform=ax.get_xaxis_transform(), fontsize=7)
            elif status == _CI_REFUSED:
                ax.scatter(part.dose_mM, part.corrected_fold, marker="x", s=65,
                           color=figures.PALETTE[1], label=status, zorder=5)
            else:
                ax.scatter(part.dose_mM, part.corrected_fold, marker="s", s=44,
                           facecolors=figures.PALETTE[0] if status == _CI_CLEAR else "white",
                           edgecolors=figures.PALETTE[0], linewidths=1.5,
                           label=status, zorder=5)
        ax.set_xscale("log")
        ax.set_xticks(panel.dose_mM)
        ax.set_xticklabels([
            f"{r.dose_mM:g}\n{int(r.plates)}/{int(r.wells)}\n"
            f"{r.n_plates_usable_naive}/{r.n_plates_usable_corrected}"
            for r in panel.itertuples()], fontsize=7.5)
        ax.minorticks_off()
        ax.set_xlim(panel.dose_mM.min() * .7, panel.dose_mM.max() * 1.4)
        ax.set_ylim(bottom, top)
        ax.set_xlabel(f"dose (mM {stressor}, log scale); P/W; N/C")
        ax.set_ylabel("fold change vs 0 mM control")
        ax.set_title(f"{construct} ({stressor})")
    for ax in list(axes.flat)[len(groups):]:
        ax.set_visible(False)
    axes.flat[0].legend(loc="upper left", fontsize=6.8, ncols=2, columnspacing=.8)
    counts = table.status.value_counts()
    fig.suptitle("Dilution correction: pointwise intervals and explicit refusals",
                 fontsize=12.5, fontweight="bold")
    _caption(fig, (
        f"{len(table)} dose conditions: {int(valid.sum())} estimable "
        f"({counts.get(_CI_CLEAR, 0)} exclude 1.0; {counts.get(_CI_SPANS, 0)} span 1.0); "
        f"{int((~valid).sum())} CI REFUSED ({counts.get(_CI_REFUSED, 0)} finite folds; "
        f"{counts.get(_FOLD_MISSING, 0)} missing folds).\n"
        f"Coverage at plotted doses: {_count_range(table.plates)} recorded plates and "
        f"{_count_range(table.wells)} recorded wells per condition; "
        f"{_count_range([*table.min_wells, *table.max_wells])} wells per plate. "
        "Below each dose, P/W: recorded plates/wells; N/C: usable naive/corrected plate folds. "
        "Wells within a plate are not independent biological replicates.\n"
        "Whiskers: nominal 95% two-stage cluster-bootstrap intervals on the log ratio, "
        "with approximate t inflation; coverage is not guaranteed. These are pointwise, "
        "not family-wise acceptance. Hollow squares require an estimable CI spanning 1.0; "
        "crosses mark refusals. Missing-fold crosses sit on an axis rail, not at zero.\n"
        "A finite fold is not a valid interval. Negative fluorescence or recovered activity "
        "is not a viability diagnosis. Neither CI overlap nor refusal establishes absence "
        "of induction. Fig04 follows the selected pointwise calls through the declared sweeps.\n"
        "Source: outputs/sensor_characterisation.csv "
        "(python3 scripts/run_sensor_characterisation.py)."
    ))
    return figures.save_figure(fig, out_dir, stem)


def _heldout_figure(table, sensor, out_dir, stem="fig02_heldout_skill_vs_saturation"):
    """Keep the fixed-form diagnostic and identify its actual baseline and scored cohort."""
    table = table.sort_values("held_out_dose")
    scored = table[np.isfinite(table.skill)]
    refused = table[~np.isfinite(table.skill)]
    count_columns = ["n_test", "n_test_available", "n_unscoreable_test", "n_fits_vouched"]
    if not np.isfinite(scored[count_columns].to_numpy(dtype=float)).all():
        raise ValueError("scored held-out doses need their reported cohort counts")
    fig, (top, bottom) = figures._figure(2, 1, figsize=(11.4, 9.6), sharex=True,
                                        height_ratios=(1.6, 1))
    top.plot(table.dose_over_ec50, table.skill, color="#888888", lw=1.2,
             label="fixed biphasic skill")
    top.axhline(0, color="black", lw=1.5)
    for row in scored.itertuples():
        colour = figures.PALETTE[2] if row.skill > 0 else figures.PALETTE[1]
        top.plot(row.dose_over_ec50, row.skill, marker="o" if row.skill > 0 else "v",
                 color=colour, mfc=colour if row.skill > 0 else "white", mew=1.5, ms=8)
        top.annotate(f"{row.held_out_dose:g} mM\n{int(row.n_test)}/{int(row.n_test_available)} "
                     f"scored; {int(row.n_unscoreable_test)} excluded\n"
                     f"{int(row.n_fits_vouched)} vouched fits",
                     (row.dose_over_ec50, row.skill), xytext=(0, 12 if row.skill > 0 else -42),
                     textcoords="offset points", ha="center", fontsize=7.5)
    top.scatter(refused.dose_over_ec50, np.full(len(refused), .06), marker="X", s=60,
                color="#666666", transform=top.get_xaxis_transform(),
                label="REFUSED (not a score)", zorder=4)
    for row in refused.itertuples():
        top.text(row.dose_over_ec50, .12, f"{row.held_out_dose:g} mM\nREFUSED",
                 ha="center", transform=top.get_xaxis_transform(), fontsize=8)
    span = max(.25, float(scored.skill.abs().max())) * 1.9 if len(scored) else 1.0
    top.set_ylim(-span, span)
    top.set_ylabel("fixed biphasic RMSE skill\nvs training-selected baseline (0 = equal RMSE)")
    top.set_title("A. Fixed-form diagnostic, not selected-pipeline performance")
    reference = float(table.reference_ec50_mM.iloc[0])
    lo, hi = table.dose_over_ec50.min() * .55, table.dose_over_ec50.max() * 1.5
    u = np.geomspace(lo, hi, 400)
    bottom.plot(u, 4 / (1 + u) ** 2, color=figures.PALETTE[0], lw=2)
    bottom.axhline(1, color="#555555", ls="--", lw=1)
    bottom.set_yscale("log")
    bottom.set_ylabel("Michaelis-term sensitivity\nrelative to its value at d = K")
    bottom.set_title("B. Functional-form illustration, conditional on full-data K")
    bottom.set_xlabel(f"dose / K (K = {reference:.3f} mM; full-data reference, log scale)")
    bottom.set_xscale("log")
    bottom.set_xlim(lo, hi)
    bottom.set_xticks(table.dose_over_ec50)
    bottom.set_xticklabels([f"{u:.1f}" for u in table.dose_over_ec50])
    bottom.minorticks_off()
    inventory = _dose_coverage(sensor).groupby("dose_mM").wells.sum()
    requested_wells = sum(int(inventory.get(d, 0)) for d in table.held_out_dose)
    refusals = " ".join(
        f"{row.held_out_dose:g} mM: {int(inventory.get(row.held_out_dose, 0))} recorded wells; "
        f"REFUSED, split counts not reported (not zero). CSV note: {row.note}."
        for row in refused.itertuples())
    baselines = []
    for row in scored.itertuples():
        mapping = json.loads(row.selected_baselines)
        identities = "; ".join(
            f"{name} = {', '.join(sorted(c for c, b in mapping.items() if b == name))}"
            for name in sorted(set(mapping.values())))
        baselines.append(f"{row.held_out_dose:g} mM: {identities}.")
    cohort = "; ".join(dict.fromkeys(scored.comparison_cohort))
    selection = "; ".join(dict.fromkeys(scored.baseline_selection))
    fig.suptitle("Fixed-form biphasic diagnostic: held-out skill varies by dose",
                 fontsize=12.5, fontweight="bold")
    _caption(fig, (
        "Plotted CSV column: skill = 1 - rmse_biphasic / rmse_selected_baseline, "
        f"not skill_selected_model. Baseline selection: {selection}.\n"
        f"{len(table)} requested doses: {len(scored)} scored, {len(refused)} REFUSED; "
        f"{requested_wells} recorded wells at these doses. On scored dose cohorts: "
        f"{int(scored.n_test.sum())}/{int(scored.n_test_available.sum())} test wells scored; "
        f"{int(scored.n_unscoreable_test.sum())} excluded. Cohort: {cohort}.\n"
        f"{refusals}\n"
        "Selected baseline identities (all fitted constructs; scoring uses only the vouched cohort):\n"
        + "\n".join(baselines) + "\n"
        "The lower curve is 4/(1 + d/K)^2 for the Michaelis induction term. K is the median "
        "of identifiable full-data biphasic EC50 fits, not an outer-training estimate. "
        "This is not out-of-sample evidence or causal proof; the diagnostic does not establish "
        "a crossover or an extrapolation rule.\n"
        "Sources: outputs/heldout_interpolation_by_dose.csv; outputs/sensor_characterisation.csv "
        "(python3 scripts/run_heldout_score.py; python3 scripts/run_sensor_characterisation.py)."
    ))
    return figures.save_figure(fig, out_dir, stem)


def _transfer_figure(recovery, nulls, out_dir, stem="fig03_transfer_honesty"):
    """Display the same headline comparisons without treating non-rejection as equivalence."""
    drawn = nulls[nulls.headline | ~nulls.beats_null].sort_values(
        ["beats_null", "claim", "dataset"]).reset_index(drop=True)
    fig, (left, right) = figures._figure(1, 2, figsize=(16.4, 8.4), width_ratios=(1, 1.5))
    for index, (definition, series) in enumerate(recovery.groupby("peer_definition", sort=False)):
        colour = figures.PALETTE[index]
        series = series.set_index("peers").reindex(figures.PEER_BINS)
        present = series.pairs > 0
        left.plot(range(len(series)), series.recovered.where(present),
                  marker=("o", "s")[index], ls=("-", "--")[index],
                  color=colour, lw=1.8, ms=8, label=definition)
        for x, row in enumerate(series.itertuples()):
            left.annotate(f"{int(row.pairs)} pairs" if row.pairs > 0 else "no pairs",
                          (x, row.recovered if row.pairs > 0 else .01),
                          xytext=(0, 12 + 12 * index), textcoords="offset points",
                          ha="center", fontsize=8, color=colour)
    left.set_xticks(range(len(figures.PEER_BINS)), figures.PEER_BINS)
    left.set_ylim(-.06, 1.25)
    left.set_xlabel("other stressors driving the same module (peers)")
    left.set_ylabel("recovered own-target pairs (fraction)")
    left.set_title("A. Coverage and recovery: descriptive association")
    left.legend(loc="lower right", fontsize=8)
    claim_labels = {
        "transfer (median own-target score)": "fitted-basis transfer",
        "redundancy law, peers from panel target lists": "peer-count association (Spearman)",
        "known-loadings readout (median Spearman)": "known-loadings readout (Spearman)",
    }
    labels = []
    for y, row in enumerate(drawn.itertuples()):
        colour = figures.PALETTE[2] if row.beats_null else figures.PALETTE[1]
        right.plot([0, row.gap], [y, y], color=colour, lw=2)
        right.plot(row.gap, y, marker="o" if row.beats_null else "X", color=colour,
                   mfc=colour if row.beats_null else "white", mew=1.8, ms=9)
        verdict = "rejects null" if row.beats_null else "does not reject"
        right.annotate(f"{row.observed:+.4f} vs {row.null_median:+.4f}   "
                       f"p = {row.p_value:.3f}\n{verdict}", (1.015, y),
                       xycoords=("axes fraction", "data"), ha="left", va="center",
                       fontsize=8, annotation_clip=False)
        claim = claim_labels.get(row.claim, row.claim)
        labels.append(f"{textwrap.fill(claim, 44)}\nvs {row.null}\n{row.dataset}")
    right.axvline(0, color="black", lw=1.5)
    right.set_yticks(range(len(drawn)), labels, fontsize=7.6)
    right.set_ylim(len(drawn) - .3, -.7)
    right.set_xlabel("observed score minus its null median\n(compare like with like within each row)")
    right.set_title("B. Fitted-basis tests and known-loadings control")
    right.grid(axis="y", visible=False)
    zero_peer = "; ".join(
        f"{r.recovered:.0%} recovery with 0 peers ({int(r.pairs)} pairs, {r.peer_definition})"
        for r in recovery[recovery.peers.eq("0") & (recovery.pairs > 0)].itertuples())
    fig.suptitle("Transfer: coverage associations and the limits of the null tests",
                 fontsize=12.5, fontweight="bold")
    _caption(fig, (
        f"A is descriptive, not a necessary-condition result: {zero_peer}. "
        "Empty bins are not zero recovery. Both recorded peer definitions are retained.\n"
        "B retains the negative rotated-subspace outcomes for transfer and peer association. "
        "Non-rejection is not equivalence or proof of no generalisation. "
        "Rejection of matched-marginals nulls does not resolve the rotated-subspace test.\n"
        "Known-loadings control and fitted-basis transfer ask different questions: "
        "a supplied loading matrix versus shuffled loadings is not a test of learning a "
        "transfer basis. Transfer uses own-target scores on arrays and held-out channel R2 "
        "in simulation; peer association and the known-loadings control use Spearman correlation. "
        "Compare each observed score only with its own null.\n"
        "Sources: outputs/external_redundancy.csv, outputs/external_nulls.csv, "
        "outputs/external_nulls_independent.csv, outputs/cross_family_transfer.csv "
        "(python3 scripts/run_external_validation.py; python3 scripts/run_cross_family.py)."
    ))
    return figures.save_figure(fig, out_dir, stem)


def _robustness_figure(table, folds, out_dir, stem="fig04_interval_robustness"):
    """Track the originally selected pointwise calls, including missing sweep bounds."""
    status = _fold_status(folds)
    calls = folds[status.eq(_CI_CLEAR)].sort_values(["construct", "dose_mM"])
    if calls.empty:
        raise ValueError("no pointwise calls with valid bounds to track")
    sweeps = list(dict.fromkeys(table.sweep))
    fig, axes = figures._figure(1, len(sweeps), figsize=(6.4 * len(sweeps), 6.4),
                                squeeze=False, sharey=True)
    doses = sorted(calls.dose_mM.unique())
    colours = dict(zip(doses, figures.series_colours(len(doses))))
    styles = dict(zip(sorted(calls.construct.unique()), ("-", "--", ":", "-.")))
    retained = np.ones(len(calls), dtype=bool)
    kept_bounds, expected_bounds = 0, 0
    settings = []
    for ax, sweep in zip(axes.flat, sweeps):
        panel = table[table.sweep == sweep]
        grid = sorted(panel.knob_value.unique())
        settings.append(f"{sweep}: " + ", ".join(f"{v:g}" for v in grid))
        ax.axhline(1, color="black", lw=1.5, label="reference (1.0)")
        for index, call in enumerate(calls.itertuples()):
            series = panel[(panel.construct == call.construct) & (panel.dose_mM == call.dose_mM)]
            series = series.set_index("knob_value").reindex(grid)
            valid = np.isfinite(series.fold) & _valid_bounds(series.estimable, series.low, series.high)
            above = call.interval_low > 1
            bound = series.low if above else series.high
            keeps = valid & ((bound > 1) if above else (bound < 1))
            retained[index] &= bool(keeps.all())
            kept_bounds += int(keeps.sum())
            expected_bounds += len(grid)
            ax.plot(grid, bound.where(valid), color=colours[call.dose_mM], marker="o", ms=5,
                    ls=styles[call.construct], lw=1.5,
                    label=f"{call.construct} {call.dose_mM:g} mM")
            ax.scatter(np.asarray(grid)[~valid], np.full(int((~valid).sum()), .04),
                       transform=ax.get_xaxis_transform(), marker="X", color=colours[call.dose_mM],
                       label="CI REFUSED / missing sweep row", s=45, zorder=5)
            ax.scatter(np.asarray(grid)[valid & ~keeps], bound[valid & ~keeps],
                       facecolors="none", edgecolors="black", s=110, zorder=5)
        ax.set_xlabel(str(panel.sweep_label.iloc[0]))
        ax.set_title(sweep.replace("_", " "))
        ax.set_xticks(grid)
    axes.flat[0].set_ylabel("95% pointwise CI bound toward 1.0\n(lower for above-1 calls; upper for below-1 calls)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    keep = [(h, label) for h, label in zip(handles, labels)
            if label != "CI REFUSED / missing sweep row"]
    axes.flat[0].legend([h for h, _ in keep], [label for _, label in keep],
                        loc="upper right", fontsize=7.3, ncols=2)
    ymin, ymax = axes.flat[0].get_ylim()
    axes.flat[0].set_ylim(ymin, ymax + (ymax - ymin) * .35)
    fig.suptitle(f"{int(retained.sum())}/{len(calls)} selected pointwise calls retain bounds "
                 "across the declared sweeps", fontsize=12.5, fontweight="bold")
    _caption(fig, (
        f"This figure follows selected {len(calls)} of {int(status.isin((_CI_CLEAR, _CI_SPANS)).sum())} "
        f"estimable dose conditions ({len(folds)} total), selected by their original pointwise "
        f"CI excluding 1.0. {kept_bounds}/{expected_bounds} bounds remain estimable and on their "
        "original side of 1.0. Crosses on the axis rail mean CI REFUSED / missing sweep row; "
        "rings mark bounds that no longer clear 1.0.\n"
        + "; ".join(settings) + ".\n"
        "This supports robustness at only the declared sweep settings, not family-wise "
        "acceptance, proof that window/autofluorescence artefacts are absent, or validity "
        "outside those settings. The two sweeps are separate, not a joint sensitivity grid.\n"
        "Sources: outputs/late_window_sensitivity.csv, outputs/autofluorescence_sensitivity.csv, "
        "outputs/sensor_characterisation.csv (python3 scripts/run_sensor_characterisation.py)."
    ))
    return figures.save_figure(fig, out_dir, stem)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures-dir", type=pathlib.Path, default=_default_figures_dir(),
                        help="where to write the SVG/PNG pairs and their source tables")
    parser.add_argument("--story-dir", type=pathlib.Path, default=None,
                        help="where to write the four story figures and their manifest; "
                             "defaults to outputs/, beside the tables they trace to")
    parser.add_argument("--power-simulations", type=int, default=240,
                        help="simulated campaigns per bar in figure 04; the Wilson "
                             "interval drawn is the precision this buys")
    args = parser.parse_args()

    outputs = paths.outputs_dir()
    figures_dir = args.figures_dir
    figures_dir.mkdir(parents=True, exist_ok=True)
    story_dir = args.story_dir if args.story_dir is not None else outputs
    story_dir.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []
    manifest: list[dict[str, str]] = []

    def record(paths_pair: tuple[pathlib.Path, pathlib.Path]) -> None:
        written.extend(paths_pair)
        for path in paths_pair:
            print(f"  wrote {paths.display_path(path)}  "
                  f"({path.stat().st_size / 1024:.0f} kB)")

    def sidecar(frame: pd.DataFrame, name: str) -> None:
        path = figures_dir / name
        frame.to_csv(path, index=False)
        written.append(path)
        print(f"  wrote {paths.display_path(path)}  ({len(frame)} rows)")

    def entry(paths_pair: tuple[pathlib.Path, pathlib.Path], question: str,
              sources: Sequence[str], produced_by: Sequence[str]) -> None:
        """Record one figure, what it settles, and the trail back to the numbers.

        The sources are named here rather than inside the drawing code because this
        script is the only place that knows which table feeds which figure, and a
        provenance line that lives anywhere else drifts from the call that reads the file.
        """
        record(paths_pair)
        svg, png = paths_pair
        manifest.append({
            "figure": svg.stem,
            "svg": paths.display_path(svg),
            "png": paths.display_path(png),
            "question": question,
            "source_tables": "; ".join(sources),
            "sources_produced_by": "; ".join(produced_by),
            "figure_produced_by": "python3 scripts/make_figures.py",
        })

    print("01  G1 optical quality, per plate and per failure signature")
    g1_summary = _pooled_to_palette(figures.g1_summary(_g1_tables(outputs)))
    sidecar(g1_summary, "g1_summary.csv")
    record(figures.g1_pass_rates(g1_summary, figures_dir,
                                 plates_without_a_table=_plates_without_a_table(outputs)))

    print("02  naive versus dilution-corrected dose response")
    sensor = figures.read_table(
        outputs / "sensor_characterisation.csv", "the sensor characterisation table",
        "python3 scripts/run_sensor_characterisation.py")
    folds = figures.dose_response_table(sensor)
    sidecar(folds, "dose_response_folds.csv")
    record(_dose_response_figure(folds, sensor, figures_dir))

    print("03  the decoupling grid")
    designs = figures.decoupling_designs()
    sidecar(designs, "decoupling_designs.csv")
    record(figures.decoupling_grid_scatter(designs, figures_dir))

    print(f"04  attribution power ({args.power_simulations} campaigns per bar; "
          f"this is the slow one)")
    power = figures.attribution_power_table(n_simulations=args.power_simulations)
    sidecar(power, "attribution_power.csv")
    record(figures.attribution_power(power, figures_dir))

    print("05  G4 anchor contamination")
    rt_minus = figures.read_table(
        outputs / "g4_rt_minus_qc.csv", "the G4 no-RT QC table",
        "python3 scripts/run_g4.py")
    verdicts = figures.read_table(
        outputs / "g4_verdicts.csv", "the G4 verdict table",
        "python3 scripts/run_g4.py")
    contamination = figures.g4_contamination_table(rt_minus, verdicts)
    sidecar(contamination, "g4_contamination.csv")
    record(figures.g4_contamination(contamination, figures_dir))

    print(f"\n{len(written)} files in {paths.display_path(figures_dir)}/")

    # ------------------------------------------------------------------
    # the story set: the four figures the argument rests on, in outputs/
    # ------------------------------------------------------------------
    print(f"\nstory set, into {paths.display_path(story_dir)}/")

    print("fig01  the dilution correction, with intervals")
    entry(
        _dose_response_figure(
            folds, sensor, story_dir, stem="fig01_dose_response_correction"),
        "Which corrected pointwise intervals exclude or span 1.0, and which of all "
        "recorded dose conditions are refused?",
        ["outputs/sensor_characterisation.csv"],
        ["python3 scripts/run_sensor_characterisation.py"],
    )

    print("fig02  the held-out prediction, and where it fails")
    by_dose = figures.read_table(
        outputs / "heldout_interpolation_by_dose.csv",
        "the per-dose held-out interpolation scores",
        "python3 scripts/run_heldout_score.py")
    ec50 = figures.identifiable_ec50(sensor)
    skill = figures.heldout_skill_table(by_dose, ec50, readings=sensor)
    entry(
        _heldout_figure(skill, sensor, story_dir),
        "How does fixed-form biphasic skill compare with training-selected baselines "
        "on the reported held-out cohort, including refused doses?",
        ["outputs/heldout_interpolation_by_dose.csv",
         "outputs/sensor_characterisation.csv (for the fitted EC50)"],
        ["python3 scripts/run_heldout_score.py",
         "python3 scripts/run_sensor_characterisation.py"],
    )

    print("fig03  what transfer does and does not recover")
    redundancy = figures.redundancy_recovery_table(figures.read_table(
        outputs / "external_redundancy.csv", "the external redundancy table",
        "python3 scripts/run_external_validation.py"))
    nulls = figures.null_comparison_table(
        [("Gasch 2000, 8 stressors", figures.read_table(
            outputs / "external_nulls.csv", "the external null table",
            "python3 scripts/run_external_validation.py")),
         ("Gasch 2000, 6 stressors (independent)", figures.read_table(
             outputs / "external_nulls_independent.csv",
             "the independent-panel external null table",
             "python3 scripts/run_external_validation.py"))],
        simulated=figures.read_table(
            outputs / "cross_family_transfer.csv", "the cross-family transfer table",
            "python3 scripts/run_cross_family.py"),
    )
    entry(
        _transfer_figure(redundancy, nulls, story_dir),
        "What descriptive coverage associations and null-test outcomes are reported "
        "for fitted-basis transfer and the known-loadings control?",
        ["outputs/external_redundancy.csv", "outputs/external_nulls.csv",
         "outputs/external_nulls_independent.csv", "outputs/cross_family_transfer.csv"],
        ["python3 scripts/run_external_validation.py",
         "python3 scripts/run_cross_family.py"],
    )

    print("fig04  do the resolved calls survive the analysis choices")
    robustness = figures.interval_robustness_table([
        figures.read_table(outputs / "late_window_sensitivity.csv",
                           "the late-window sensitivity sweep",
                           "python3 scripts/run_sensor_characterisation.py"),
        figures.read_table(outputs / "autofluorescence_sensitivity.csv",
                           "the autofluorescence sensitivity sweep",
                           "python3 scripts/run_sensor_characterisation.py"),
    ])
    entry(
        _robustness_figure(robustness, folds, story_dir),
        "Do the selected pointwise calls retain their bounds at every declared "
        "late-window and autofluorescence sweep setting?",
        ["outputs/late_window_sensitivity.csv",
         "outputs/autofluorescence_sensitivity.csv",
         "outputs/sensor_characterisation.csv (for which calls to track)"],
        ["python3 scripts/run_sensor_characterisation.py"],
    )

    manifest_path = story_dir / "figures_manifest.csv"
    pd.DataFrame(manifest).to_csv(manifest_path, index=False)
    print(f"  wrote {paths.display_path(manifest_path)}  ({len(manifest)} figures)")
    print("\nAssemble the report with: python3 scripts/build_report.py")


if __name__ == "__main__":
    main()
