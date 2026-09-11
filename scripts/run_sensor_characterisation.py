"""Sensor characterisation: is each reporter's induction promoter activity or growth?

DTT and H2O2 both slow growth. A stable reporter accumulates as growth slows, purely
by reduced dilution, so an apparent dose response can be a growth response. With the
plate map recovered, growth rate is available per dose and the two can be separated.

Every fold change here carries what is known about its spread, or says that nothing is.
The dose response used to be printed as bare ratios of pooled means -- ``0.96``, ``1.01``,
``0.97`` -- which read as though they were distinguishable from 1.0 when no one had
computed whether they are. ``analysis/uncertainty.py`` supplies the cluster bootstrap that
answers that, and refuses an interval below three biological replicates rather than
manufacturing a tight-looking one. Most rows below are refused, and the refusal is the
result: it is what two plates buy.

Two consequences of wiring it in are visible in the printed tables and are not bugs.

The point estimates move, because the fold is now formed **within** each plate and only
then combined across plates, rather than as a ratio of means pooled over plates. A
plate-level shift -- inoculum density, medium batch, reader gain -- multiplies a plate's
dosed and control wells alike and therefore cancels in the first and not in the second.
The rows that move are listed after the table rather than left for a reader to notice.

And the table is swept against the one term the plates never measured. No reporter-free
strain was ever read in the fluorescence channel, so per-biomass autofluorescence is still
inside the per-cell signal; the media blank that is subtracted is the *background* term and
does not touch it. Because the ODE inversion is exact, an unremoved constant ``a`` in the
per-cell signal contributes exactly ``a * mu`` to the recovered activity -- a constant has
zero derivative -- so the whole table can be re-read at any assumed ``a`` without
re-inverting anything. That sweep separates the conclusions that survive the missing
measurement from the ones that depend on it.

Usage: python scripts/run_sensor_characterisation.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.uncertainty import (
    MIN_PLATES_FOR_INTERVAL,
    activity_uncertainty,
    equivalence,
    fold_change,
)
from ystwin.diagnostics.dilution_confound import detect_blank_wells
from ystwin.growth import (
    growth_rate_uncertainty,
    growth_window,
    max_specific_growth_rate,
    specific_growth_rate,
)
from ystwin.plate.dose_response import read_dose_response
from ystwin.plate.layout import (
    NEWPROTOCOL_LAYOUT,
    blanks_for_export,
    plate_key,
    recorded_for_export,
    recorded_well_roles,
    recover_layout,
)
from ystwin.plate.synergy import read_synergy_kinetic
from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT
from ystwin.reporter import (
    ReporterKinetics,
    promoter_activity,
    promoter_activity_from_total,
)
from ystwin.readings import CorrectedOD, CorrectedRFU, RawOD, RawRFU, SpecificFluorescence

OUT = paths.outputs_dir()
KINETICS = ReporterKinetics(k_deg=0.0)  # D2: mCitrine loss is pure dilution

NO_INDUCTION_MARGIN = (0.9, 1.1)
"""What counts as "no induction", named before the numbers were seen.

An equivalence claim needs its margin fixed in advance or it is a description of the
data rather than a test of it. Ten percent either way is the smallest fold this assay
could plausibly resolve at three technical replicates.
"""

AUTOFLUORESCENCE_FRACTIONS = (0.0, 0.05, 0.10, 0.20, 0.30)

# Where the "late" window starts, as a fraction of the points that clear the density
# floor. 0.75 is a choice, not a measurement, so it is swept rather than asserted.
LATE_WINDOW_FRACTION = 0.75
LATE_WINDOW_SWEEP = (0.5, 0.6, 0.7, 0.75, 0.8, 0.9)
"""Assumed autofluorescence, as a fraction of each construct's own 0 mM per-cell signal.

Expressed as a fraction rather than an absolute level because the constructs sit at very
different signal scales, so no single absolute figure would mean the same thing to all
four. Zero is what the analysis currently assumes without saying so; the upper end is
chosen to be larger than anyone would defend, since the point of the sweep is to find
where a conclusion breaks rather than to guess the true value. Nothing here picks a value
out of the range: the measurement that would is missing -- see ``docs/DATA_INVENTORY.md``.
"""

_SWING_TOLERANCE = 0.05
"""How far a fold may move across the sweep before the sweep is worth reporting.

Half the half-width of :data:`NO_INDUCTION_MARGIN`. A fold that moves by less than this
under every assumed autofluorescence is being reported to a precision the assumption does
not disturb; one that moves by more is partly a statement about the assumption, whether or
not it happens to cross the margin.
"""


def prepare(path, recorded_blanks=None):
    run = read_synergy_kinetic(path)
    odc = run.raw_channel("OD600").channel
    flc = run.raw_channel("mCitrine").channel
    al = run.aligned(odc)
    od, fl = al[odc], al[flc]
    shared = [w for w in od.columns if w in set(fl.columns)]
    roles = recorded_well_roles(od, fl)
    # Two blanks: layout matching needs the team's, physiology needs the non-growing one.
    permissive = ([w for w in shared if od[w].median() < 0.15]
                  if roles is None and recorded_blanks is None else [])

    # The logbook records blank positions per plate, and they moved between plates, so
    # a recorded position is used when the export still has it and detection is fallback.
    strict = []
    if recorded_blanks is not None:
        conflicts = ([] if roles is None else
                     [well for well in recorded_blanks if roles.get(well) == "culture"])
        if conflicts:
            raise ValueError(f"supplied blanks are recorded as culture wells: {conflicts}")
        strict = [w for w in recorded_blanks if w in shared]
        if not strict:
            print(f"  logbook records blanks at {recorded_blanks}, none present in "
                  "both channels; recorded cultures cannot replace the missing blanks")
            return None
    else:
        strict = [w for w in detect_blank_wells(RawOD(od[shared]), rfu=RawRFU(fl[shared])) if w in shared]
    if not strict:
        return None
    if roles is not None or recorded_blanks is not None:
        permissive = list(strict)
    ob_match = float(od[permissive].to_numpy().mean()) if permissive else 0.0
    rb_match = float(fl[permissive].to_numpy().mean()) if permissive else 0.0
    ob = float(od[strict].to_numpy().mean())
    rb = float(fl[strict].to_numpy().mean())
    cult = ([w for w in shared if roles.get(w) == "culture"] if roles is not None
            else [w for w in shared if w not in set(permissive) | set(strict)])
    for_matching = (fl[cult] - rb_match) / (od[cult] - ob_match)
    normalised = (fl[cult] - rb) / (od[cult] - ob)
    if abs(ob_match - ob) > 1e-6:
        print(f"  blanks: physiology {sorted(strict)} -> OD {ob:.4f}; "
              f"layout-matching {sorted(permissive)} -> OD {ob_match:.4f}")
    return od[cult], fl[cult], normalised, for_matching, ob, rb


def _plate_level(readings: pd.DataFrame) -> pd.DataFrame:
    """Collapse the per-well rows back to one row per plate x construct x dose.

    The well rows are a decomposition of the plate aggregate, not a different estimate:
    with the growth rate taken from the condition's pooled optical density, the activity
    inversion is a *linear* operator on the per-cell signal, so averaging the wells'
    activities and inverting the wells' average signal give the same number to machine
    precision. That is why this collapse is lossless and why the printed table is
    unchanged by the switch to per-well rows.
    """
    keys = ["plate", "construct", "stressor", "dose_mM"]
    numeric = [c for c in readings.columns if c not in {*keys, "well"}]
    return readings.groupby(keys, as_index=False)[numeric].mean()


def _pooled_ratio(readings: pd.DataFrame, construct: str, dose: float,
                  value: str = "activity_late", control_dose: float = 0.0) -> float:
    """The fold as it used to be reported: a ratio of means pooled across plates.

    Kept only so the size of the change can be shown rather than asserted. It is the
    wrong estimator -- it lets a plate-level gain difference into a quantity that is
    supposed to be free of it -- but it is the one behind the numbers in ``README.md``.
    """
    subset = readings[readings.construct == construct]
    per_plate = subset.groupby(["plate", "dose_mM"])[value].mean()
    dosed = per_plate[np.isclose(per_plate.index.get_level_values("dose_mM"), dose)]
    control = per_plate[np.isclose(per_plate.index.get_level_values("dose_mM"), control_dose)]
    if dosed.empty or control.empty or abs(control.mean()) < 1e-12:
        return float("nan")
    return float(dosed.mean() / control.mean())


def _interval_cell(fold) -> str:
    """Render an interval, or render its absence so it cannot be mistaken for one."""
    if not fold.estimable:
        return f"not estimable (n={fold.n_plates})"
    return f"[{fold.low:.2f}, {fold.high:.2f}]"


def _folds_for(readings: pd.DataFrame, construct: str, doses) -> dict:
    """Every non-control dose's fold change, skipping doses this plate set cannot pair."""
    out = {}
    for dose in doses:
        if np.isclose(dose, 0.0):
            continue
        try:
            out[dose] = fold_change(readings, construct, dose)
        except ValueError:
            continue
    return out


def one_file_per_plate(exports) -> tuple[list, list[str]]:
    """Keep one file per logbook plate, and say which were set aside.

    A biological replicate is a culture grown and read once, and ``n_plates`` is what every
    interval in this analysis is computed from -- so counting a plate twice does not add
    evidence, it manufactures it. Two files can describe one read: ``20260804`` reaches
    ``data/plates`` both as Gen5's own ``.xpt`` and as the blank-subtracted ``.xlsx`` it
    exported from that file, and both are committed on purpose (the export is the evidence
    of what the transform removed).

    When that happens the instrument file wins, and this is an ordering rather than a
    preference: a ``.xpt`` is what Gen5 wrote and an ``.xlsx`` is an export of one, so the
    archive cannot be a lossy copy of the workbook while the workbook demonstrably can be a
    lossy copy of the archive. On ``20260804`` it is: the export went through a
    blank-subtraction transform, and with the background gone seven high-dose culture wells
    fall under the density floor this analysis screens on and are lost.

    Two files for one plate with **no** instrument file among them is a different situation
    and is refused rather than guessed. Nothing distinguishes them, and picking either would
    decide which measurement is the replicate by list order.

    Args:
        exports: Plate export paths, in manifest order.

    Returns:
        ``(kept, set_aside)`` -- the paths to analyse, and the names of the files dropped
        because another file already covers their plate.

    Raises:
        ValueError: for two non-instrument files describing one logbook plate.
    """
    by_plate: dict[str, list] = {}
    kept, set_aside = [], []
    for path in exports:
        key = plate_key(path.name)
        if key is None:  # not a logbook plate: the July exports, the BY4741 controls
            kept.append(path)
            continue
        by_plate.setdefault(key, []).append(path)
    # `copies`, not `paths`: `ystwin.paths` is imported at module scope and is how every
    # real plate in this script is located, and a loop variable of that name shadowed it.
    # Function-local, so nothing was ever wrong -- but the day this loop moves out of a
    # function it would resolve plate locations against a list of files.
    for key, copies in by_plate.items():
        if len(copies) == 1:
            kept.append(copies[0])
            continue
        archives = [p for p in copies if p.suffix.lower() == ".xpt"]
        if len(archives) != 1:
            raise ValueError(
                f"plate {key} has {len(copies)} committed files "
                f"({', '.join(p.name for p in copies)}) and "
                f"{len(archives)} of them are instrument files. One plate is one "
                "biological replicate, and nothing here says which of these is it.")
        kept.append(archives[0])
        for other in copies:
            if other is not archives[0]:
                set_aside.append(other.name)
                print(f"{other.name[:48]}: plate {key} is also committed as "
                      f"{archives[0].name}, the instrument's own file; using that one")
    keep = {path.name for path in kept}
    return [path for path in exports if path.name in keep], set_aside


def collect(exports, late_fraction: float = LATE_WINDOW_FRACTION
            ) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """One row per well, with the condition-level quantities alongside.

    The committed table used to hold plate means only, which makes a within-plate
    bootstrap impossible to run from the artefact: the resampling stage that needs the
    wells had nothing to resample. Emitting the wells costs three lines and is the
    difference between a table that supports the analysis it feeds and one that does not.

    Returns the readings and, separately, every export that contributed nothing and the
    reason. A plate dropped in silence is indistinguishable from a plate that was never
    run, and the count of biological replicates is the number the whole analysis turns on.

    Args:
        exports: Plate export paths.
        late_fraction: Where the late window starts, as a fraction of the points that
            clear the density floor. Swept by :func:`report_late_window_sensitivity`.
    """
    rows: list[dict] = []
    skipped: dict[str, list[str]] = {"no dose ladder": [], "no fluorescence blank": [],
                                     "a second file for a plate already counted": []}
    exports, duplicates = one_file_per_plate(exports)
    skipped["a second file for a plate already counted"].extend(duplicates)
    for path in exports:
        recorded = recorded_for_export(path.name)
        try:
            derived = read_dose_response(path)
        except ValueError:
            # No per-construct sheets. The layout and the dose ladder are both in the
            # logbook registry, so a recorded plate still places every well; only an
            # unrecorded one has nothing to go on. 20260728 is recovered this way, and
            # its recorded blanks are the wells detection independently finds.
            if recorded is None:
                skipped["no dose ladder"].append(path.name)
                continue
            derived = None
        prepared = prepare(path, recorded_blanks=blanks_for_export(path.name))
        if prepared is None:
            print(f"{path.name[:48]}: no blanks in the shared wells; skipped")
            skipped["no fluorescence blank"].append(path.name)
            continue
        od, fl, normalised, for_matching, ob, rb = prepared

        print(f"\n{path.name[:60]}")
        if derived is None:
            layout = recorded.layout
            print("  layout: no per-construct sheets in this export; using the recorded "
                  "layout and dose ladder, which cannot be cross-checked against it")
        else:
            found = recover_layout(for_matching, derived)
            print(f"  layout: {found.summary()}")
            layout = found.layout if found.margin_over_runner_up > 1.5 else NEWPROTOCOL_LAYOUT
            if found.margin_over_runner_up <= 1.5:
                print("  margin too small to trust this plate's own recovery; "
                      "using the recorded layout")

        t = normalised.index.to_numpy(dtype=float)
        constructs = (sorted(derived.construct.unique()) if derived is not None
                      else sorted(layout.construct_columns))
        for construct in constructs:
            if derived is not None:
                doses = sorted(derived[derived.construct == construct].dose_mM.unique())
            else:
                doses = sorted(recorded.doses_mM[STRESSOR_FOR_CONSTRUCT[construct]])
            for i, dose in enumerate(doses):
                try:
                    wells = layout.wells(construct, i)
                except (IndexError, KeyError):
                    continue
                wells = [w for w in wells if w in normalised.columns]
                if not wells:
                    continue
                od_mean = (od[wells].mean(axis=1) - ob).to_numpy()
                if not np.all(od_mean > 0):
                    continue
                # Restrict to where density clears the detection floor.
                keep = growth_window(t, CorrectedOD(od_mean))
                if keep.sum() < 10:
                    continue
                tw, odw = t[keep], od_mean[keep]
                mu = specific_growth_rate(tw, CorrectedOD(odw))
                late = slice(int(len(tw) * late_fraction), None)

                # The growth rate is a property of the condition, not of a single well:
                # the three replicates share one culture, and a rate fitted to one noisy
                # trace is a worse estimate than one fitted to their mean. Sharing it also
                # keeps the well rows an exact decomposition of the plate aggregate.
                pooled = normalised[wells].mean(axis=1).to_numpy()[keep]
                spread = activity_uncertainty(tw, SpecificFluorescence(pooled), CorrectedOD(odw), mu, KINETICS)
                sigma_late = float(np.nanmean(spread.sigma[late]))
                # Fitted over the late window rather than the whole trace, because the
                # reported rate is the late-window one and a log-linear fit across the
                # full curve would charge the growth rate for the curve's own curvature.
                mu_se = growth_rate_uncertainty(tw[late], CorrectedOD(odw[late]))

                condition = {
                    "plate": path.stem[:22], "construct": construct,
                    "stressor": STRESSOR_FOR_CONSTRUCT[construct], "dose_mM": dose,
                    "mu_max": max_specific_growth_rate(t, CorrectedOD(od_mean)),
                    "n_points": int(keep.sum()),
                    "mu_late": float(np.nanmean(mu[late])),
                    "mu_late_se": float(mu_se),
                    "activity_late_sigma": sigma_late,
                    "reporter_cv": float(spread.reporter_cv),
                }
                for well in wells:
                    signal = normalised[well].to_numpy()[keep]
                    own_od = (od[well] - ob).to_numpy()
                    # The total-signal route: mu cancels algebraically, so the growth
                    # estimate's error does not propagate. 40% smoother on these wells.
                    # Kept beside the per-cell route rather than replacing it, so the
                    # shift is measured on every run instead of asserted once.
                    total = (fl[well] - rb).to_numpy()[keep]
                    usable_od = own_od[keep]
                    if np.all(usable_od > 0):
                        activity = promoter_activity_from_total(
                            tw, CorrectedRFU(total), CorrectedOD(usable_od), KINETICS)
                    else:
                        activity = promoter_activity(tw, SpecificFluorescence(signal), mu, KINETICS)
                    per_cell_activity = promoter_activity(tw, SpecificFluorescence(signal), mu, KINETICS)
                    # One well's own growth rate is not what the activity was corrected
                    # by, but its error is what says whether the condition's pooled rate
                    # was worth pooling, so it is carried rather than recomputed later.
                    own_se = (
                        float(growth_rate_uncertainty(tw[late], CorrectedOD(own_od[keep][late])))
                        if np.all(own_od[keep][late] > 0) else float("nan")
                    )
                    rows.append({
                        **condition, "well": well,
                        "naive_late": float(np.nanmean(signal[late])),
                        "activity_late": float(np.nanmean(activity[late])),
                        "activity_late_percell": float(np.nanmean(per_cell_activity[late])),
                        "mu_max_well": max_specific_growth_rate(t, CorrectedOD(own_od)),
                        "mu_late_se_well": own_se,
                    })
    return pd.DataFrame(rows), skipped


def report_dose_response(df: pd.DataFrame) -> None:
    """The dose-response table, with an interval on every fold or a refusal in its place."""
    print("\n" + "=" * 92)
    print("Does the induction survive dilution correction, and is the answer resolvable?")
    print("=" * 92)

    moved: list[tuple[str, float, float, float]] = []
    for (construct, stressor), g in df.groupby(["construct", "stressor"]):
        plates = _plate_level(g)
        agg = plates.groupby("dose_mM").agg(
            mu_max=("mu_max", "mean"), naive=("naive_late", "mean"),
            activity=("activity_late", "mean"), mu_se=("mu_late_se", "mean"),
            sigma=("activity_late_sigma", "mean"), cv=("reporter_cv", "mean"),
            n=("plate", "nunique"),
        )
        base = agg.iloc[0]
        folds = _folds_for(g, construct, agg.index)
        wells_per_plate = int(round(len(g) / max(g.plate.nunique() * len(agg), 1)))

        print(f"\n  {construct}  ({stressor}, {int(agg.n.max())} plates x "
              f"{wells_per_plate} wells)")
        print(f"    growth-rate SE median {agg.mu_se.median():.3f} /h; reader CV median "
              f"{agg.cv.median():.2%}; activity sigma median "
              f"{(agg.sigma / agg.activity.abs()).median():.1%} of activity")
        print(f"   {'dose':>5} {'mu_max':>7} {'growth%':>8} {'naive':>7} {'corrected':>10} "
              f"  {'95% CI on corrected':<21} {'verdict':<12}")
        for dose, r in agg.iterrows():
            nf = r.naive / base.naive if base.naive else np.nan
            af = r.activity / base.activity if base.activity else np.nan
            if np.isclose(dose, 0.0):
                cell, verdict = "(control)", "--"
            elif dose not in folds:
                cell, verdict = "no paired control", "--"
            else:
                fold = folds[dose]
                af = fold.point
                if not np.isfinite(af):
                    # Recovered activity went negative: the culture is below background
                    # here, so a ratio is not a small fold, it is undefined.
                    cell, verdict = "activity <= 0", "--"
                else:
                    cell = _interval_cell(fold)
                    verdict = equivalence(fold, NO_INDUCTION_MARGIN).verdict
                    old = _pooled_ratio(g, construct, float(dose))
                    if old > 0 and abs(old - af) >= 0.02:
                        moved.append((construct, float(dose), old, af))
            shown = f"{af:.2f}" if np.isfinite(af) else "n/a"
            print(f"   {dose:>5.1f} {r.mu_max:>7.3f} {r.mu_max / base.mu_max:>7.0%} "
                  f"{nf:>7.2f} {shown:>10}   {cell:<21} {verdict:<12}")

    print("\n  Reading the table. A refused interval is not a small one: at two biological")
    print("  replicates a cluster bootstrap has three distinct resamples, two of them")
    print("  degenerate, so any interval it produced would describe which of two plates was")
    print("  drawn. INCONCLUSIVE therefore means unresolved, not refuted -- the point")
    print("  estimates and the mechanism behind them are unchanged.")

    if moved:
        print("\n  Folds that move against the pooled-ratio estimator behind README.md's")
        print("  table, which ratios means pooled across plates rather than forming the")
        print("  fold within each plate first:")
        for construct, dose, old, new in sorted(moved, key=lambda m: -abs(m[2] - m[3])):
            print(f"    {construct:<12} @ {dose:>4.1f} mM   {old:.2f} -> {new:.2f}"
                  f"   ({new - old:+.2f})")


def report_late_window_sensitivity(exports) -> pd.DataFrame:
    """Re-derive every fold with the late window starting at 0.5 through 0.9.

    ``int(len(tw) * 0.75)`` is a choice about where a trace has settled, not something
    measured, and every reported fold and growth rate is computed inside it. Sweeping it
    is the difference between "the folds are near 1.0" and "the folds are near 1.0 for
    one arbitrary window". Reports the span each fold moves across the sweep, so a call
    that only holds at 0.75 cannot pass as a stable one.
    """
    swept: list[dict] = []
    for fraction in LATE_WINDOW_SWEEP:
        frame, _ = collect(exports, late_fraction=fraction)
        if frame.empty:
            continue
        for (construct, stressor), g in frame.groupby(["construct", "stressor"]):
            # The uncorrected fold, carried alongside the corrected one. It is what the
            # plate reader shows before anything is done to it, so it is half of the
            # headline claim -- "looks like X, actually is Y" -- and until it was written
            # here only Y had a cell an audit could point at. It moves with the window for
            # the same reason Y does, so it belongs in the sweep and not in a constant.
            naive = _plate_level(g).groupby("dose_mM").naive_late.mean()
            control = naive.iloc[0] if len(naive) else float("nan")
            for dose, fold in _folds_for(frame, construct, sorted(g.dose_mM.unique())).items():
                swept.append({"late_fraction": fraction, "construct": construct,
                              "stressor": stressor, "dose_mM": dose,
                              "naive_fold": float(naive.get(dose, float("nan")) / control)
                              if control else float("nan"),
                              "fold": float(fold.point), "low": float(fold.low),
                              "high": float(fold.high), "n_plates": int(fold.n_plates),
                              "estimable": bool(fold.estimable)})
    sweep = pd.DataFrame(swept)
    if sweep.empty:
        print("\n  no fold was estimable at any late window; nothing to sweep")
        return sweep

    print("\n" + "=" * 92)
    print("Late-window sensitivity: does the call survive moving the window?")
    print("=" * 92)
    print(f"  {'construct':14s} {'dose':>6s} {'at 0.75':>8s} {'min':>7s} {'max':>7s} "
          f"{'span':>7s}  verdict")
    for (construct, dose), g in sweep.groupby(["construct", "dose_mM"]):
        at_default = g[np.isclose(g.late_fraction, LATE_WINDOW_FRACTION)].fold
        base = float(at_default.iloc[0]) if len(at_default) else float("nan")
        lo, hi, verdict = window_verdict(g.fold)
        print(f"  {construct:14s} {dose:6g} {base:8.2f} {lo:7.2f} {hi:7.2f} "
              f"{hi - lo:7.2f}  {verdict}")
    return sweep


def window_verdict(folds) -> tuple[float, float, str]:
    """Range of a fold across the window sweep, and what that range permits saying.

    A fold whose sweep straddles 1.0 cannot be called induction or repression at all --
    the direction is a property of the window, not of the biology. One that stays on a
    side is at least consistent across it.
    """
    values = np.asarray(folds, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan"), "no window gave an estimable fold"
    lo, hi = float(values.min()), float(values.max())
    straddles = lo < 1.0 < hi
    return lo, hi, "straddles 1.0 across the sweep" if straddles else "keeps its sign"


def report_autofluorescence(df: pd.DataFrame) -> pd.DataFrame:
    """Sweep the never-measured autofluorescence and print which calls depend on it."""
    swept: list[dict] = []
    for fraction in AUTOFLUORESCENCE_FRACTIONS:
        frame = df.copy()
        # a is per construct and per plate, because the 0 mM per-cell signal it is scaled
        # against is: the constructs sit at different expression levels and the plates at
        # different gains, and a single global constant would mix the two.
        control = frame[np.isclose(frame.dose_mM, 0.0)]
        level = control.groupby(["plate", "construct"]).naive_late.mean().rename("control_signal")
        frame = frame.merge(level, on=["plate", "construct"], how="left")
        frame["activity_corrected"] = (
            frame.activity_late - fraction * frame.control_signal * frame.mu_late
        )
        for (construct, stressor), g in frame.groupby(["construct", "stressor"]):
            for dose in sorted(g.dose_mM.unique()):
                if np.isclose(dose, 0.0):
                    continue
                try:
                    fold = fold_change(g, construct, float(dose), value="activity_corrected")
                except ValueError:
                    continue
                swept.append({
                    "construct": construct, "stressor": stressor, "dose_mM": float(dose),
                    "autofluorescence_fraction": fraction, "fold": fold.point,
                    "low": fold.low, "high": fold.high, "estimable": fold.estimable,
                    "n_plates": fold.n_plates,
                    "verdict": equivalence(fold, NO_INDUCTION_MARGIN).verdict,
                })
    sweep = pd.DataFrame(swept)

    print("\n" + "=" * 92)
    print("Autofluorescence was never measured on any plate: how much does that cost?")
    print("=" * 92)
    print("  Both BY4741 reporter-free controls carry an OD600 channel and no fluorescence")
    print("  channel, so the per-biomass term was never read. Subtracting the media blank")
    print("  removes the background term and leaves this one untouched. An unremoved")
    print("  constant a in the per-cell signal adds exactly a*mu to the recovered activity,")
    print("  a constant having zero derivative, so the table can be re-read at any a.")
    print("  a is a fraction of each construct's own 0 mM per-cell signal on that plate.")
    header = " ".join(f"{f:>7.0%}" for f in AUTOFLUORESCENCE_FRACTIONS)
    print(f"\n  {'construct':<12} {'dose':>5} {header} {'swing':>7}  depends on it?")
    wide = sweep.pivot_table(index=["construct", "dose_mM"],
                             columns="autofluorescence_fraction", values="fold")
    for (construct, dose), r in wide.iterrows():
        values = [r.get(f, np.nan) for f in AUTOFLUORESCENCE_FRACTIONS]
        if not np.all(np.isfinite(values)):
            print(f"  {construct:<12} {dose:>5.1f} " + " ".join(f"{'n/a':>7}" for _ in values)
                  + f" {'':>7}  activity below background at this dose")
            continue
        swing = float(np.nanmax(values) - np.nanmin(values))
        low, high = NO_INDUCTION_MARGIN
        inside = [low <= v <= high for v in values]
        if not (all(inside) or not any(inside)):
            note = "CHANGES THE CALL"
        elif swing > _SWING_TOLERANCE:
            note = "shifts"
        else:
            note = "robust"
        print(f"  {construct:<12} {dose:>5.1f} " + " ".join(f"{v:>7.2f}" for v in values)
              + f" {swing:>7.2f}  {note}")
    print(f"\n  CHANGES THE CALL: the row crosses the no-induction margin "
          f"{NO_INDUCTION_MARGIN} inside the")
    print("  sweep, so at one end of the range it reads as no induction and at the other it")
    print("  does not. The conclusion is a property of the unmeasured constant.")
    print(f"  shifts: the fold moves by more than {_SWING_TOLERANCE:.2f} without crossing the")
    print("  margin. The direction of the call survives; its size does not, and since the")
    print("  correction only ever raises a fold it enlarges the genuine inductions too.")
    print("  robust: the row does not depend on the missing measurement.")
    print("  One plate of BY4741 across the same ladder in the mCitrine channel, same gain,")
    print("  same protocol, replaces the whole sweep with a number.")
    return sweep


def main() -> None:
    plates = paths.resolve_or_exit(paths.biosensor_plates(),
                                   "the NewProtocol biosensor plate exports", "YSTWIN_PLATES")
    exports = sorted(plates.glob("*.xlsx"))
    if not exports:
        raise SystemExit(f"no .xlsx exports found in {plates}")

    df, skipped = collect(exports)
    if df.empty:
        # Overwriting the table with an empty one would read as "no induction anywhere".
        raise SystemExit(f"no export in {plates} yielded a usable dose ladder; "
                         f"{OUT / 'sensor_characterisation.csv'} left untouched")
    df.to_csv(OUT / "sensor_characterisation.csv", index=False)

    # A consumer that averages the well rows weights a plate by its well count, where the
    # old plate rows weighted every plate equally. The two agree exactly on a balanced
    # condition, so say which conditions are not balanced rather than leave the difference
    # to be discovered as a discrepancy.
    counts = df.groupby(["plate", "construct", "dose_mM"]).well.count()
    unbalanced = counts[counts != counts.mode().iloc[0]]
    for (plate, construct, dose), n in unbalanced.items():
        print(f"\n  {construct} @ {dose:g} on {plate} has {n} wells, not "
              f"{counts.mode().iloc[0]}: a well fell into the blank set. Averaging the "
              f"well rows across plates weights this one unequally.")

    report_dose_response(df)
    sweep = report_autofluorescence(df)
    sweep.to_csv(OUT / "autofluorescence_sensitivity.csv", index=False)

    window = report_late_window_sensitivity(exports)
    if not window.empty:
        window.to_csv(OUT / "late_window_sensitivity.csv", index=False)

    realised = df.groupby("construct").plate.nunique()
    print("\n" + "=" * 92)
    print("Replication actually behind the table")
    print("=" * 92)
    for construct, n in realised.items():
        marker = "" if n >= MIN_PLATES_FOR_INTERVAL else "  <- below the interval floor"
        print(f"  {construct:<12} {n} biological replicate(s){marker}")
    for reason, names in skipped.items():
        for name in names:
            print(f"  dropped: {name[:56]:<58} ({reason})")
    if int(realised.max()) < MIN_PLATES_FOR_INTERVAL:
        print(f"\n  No construct reaches {MIN_PLATES_FOR_INTERVAL} plates, so every fold above is")
        print("  refused an interval. That is a statement about the dataset, not about the")
        print("  sensors: nothing was wrong with the dropped cultures, the measurement was")
        print("  simply not taken. See docs/DATA_INVENTORY.md for what each plate is missing.")
    print(f"\n  wrote {OUT / 'sensor_characterisation.csv'} ({len(df)} well rows)")
    print(f"  wrote {OUT / 'autofluorescence_sensitivity.csv'} ({len(sweep)} rows)")


if __name__ == "__main__":
    main()
