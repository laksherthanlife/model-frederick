"""D2: how much of each reporter channel is growth rate rather than promoter activity?

Every number printed here is a median over the wells of one plate, and a median over wells
used to be printed bare. It now carries the spread of the wells it was taken over, with the
one qualification that makes the spread readable: it is a within-plate quantity and not a
biological one. A plate is a single biological replicate however many wells it holds, so
nothing in this script can say how a second plate would come out. See
:func:`well_interval`.

Usage: python scripts/run_d2.py
Writes per-well tables to outputs/ and prints the plate-level decision summary.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.diagnostics.dilution_confound import plate_dilution_report
from ystwin.plate import replay
from ystwin.plate.synergy import read_synergy_kinetic

OUT = paths.outputs_dir()
# mCitrine is a stable YFP: the textbook assumption is that loss is pure dilution.
ASSUMED = {"k_deg": 0.0}

MIN_WELLS_FOR_INTERVAL = 5
"""Fewest wells over which a within-plate percentile interval is worth printing.

Below this the bootstrap of a median is enumerating a handful of orderings rather than
approximating a sampling distribution, and the interval it draws is an artefact of the
resample grid. Refused rather than shown, on the same principle as
``analysis/uncertainty.MIN_PLATES_FOR_INTERVAL``.
"""


def well_interval(values, fmt: str = "{:.3f}", statistic=np.median,
                  confidence: float = 0.95, seed: int = 0) -> str:
    """A statistic over one plate's wells, with the spread of those wells around it.

    This is deliberately *not* the interval that ``analysis/uncertainty.py`` computes, and
    the difference is the whole reason it is labelled. That module resamples plates,
    because plates are the unit of biological replication; wells on one plate share an
    inoculum, a medium batch, a reader and a position in an incubator, so resampling them
    reproduces only the within-plate variance component. These two survey scripts read one
    plate at a time and therefore have exactly one biological replicate each -- there is no
    between-plate component available to estimate, and none is claimed.

    What the interval does say is worth having: how firmly this plate's own wells pin the
    number down. A median dilution R^2 of 0.9 over wells that all agree is a different
    statement from the same median over wells that range from 0.2 to 1.0, and the bare
    median cannot tell those apart.

    Args:
        values: One value per well; NaNs are dropped and counted out.
        fmt: Format applied to the point estimate and both bounds.
        statistic: Summary to bootstrap. Anything ``scipy.stats.bootstrap`` accepts.
        confidence: Nominal coverage of the percentile interval.
        seed: Random seed, so a rerun prints the same bounds.
    """
    finite = np.asarray(pd.Series(values, dtype="float64").dropna(), dtype=float)
    if finite.size == 0:
        return "not reportable (no well carried this quantity)"
    point = fmt.format(float(statistic(finite)))
    if finite.size < MIN_WELLS_FOR_INTERVAL:
        return (f"{point} (no interval: {finite.size} well(s), "
                f"{MIN_WELLS_FOR_INTERVAL} needed)")
    if float(np.ptp(finite)) == 0.0:
        return f"{point} (no interval: all {finite.size} wells identical)"
    from scipy import stats

    ci = stats.bootstrap(
        (finite,), statistic, confidence_level=confidence, n_resamples=2000,
        method="percentile", random_state=seed,
    ).confidence_interval
    return (f"{point} [{fmt.format(ci.low)}, {fmt.format(ci.high)}] "
            f"over {finite.size} wells, within-plate")


def analyse(path: pathlib.Path) -> None:
    run = read_synergy_kinetic(path)
    reporters = [b.channel for b in run.blocks if b.fluorophore != "OD600"]
    if "OD600" not in run.channel_names or not reporters:
        return
    aligned = run.aligned("OD600")
    for reporter in reporters:
        from ystwin.reporter import ReporterKinetics

        rep = plate_dilution_report(
            aligned, od_channel="OD600", reporter_channel=reporter,
            kinetics=ReporterKinetics(**ASSUMED),
        )
        pw = rep.per_well
        tag = f"{path.stem}__{reporter}".replace(" ", "_").replace("&", "and")
        pw.to_csv(OUT / f"d2_{tag}.csv", index=False)

        clean = pw[(pw.mu_min > -1e-9) & (pw.qss_fraction > 0.5)]
        print(f"\n{'=' * 78}\n{path.name}  |  reporter {reporter}\n{'=' * 78}")
        print(f"  {rep.summary()}")
        print(f"  media blanks {rep.blank_wells} -> OD {rep.od_blank:.3f}, RFU {rep.rfu_background:.0f}")
        print(f"  verdicts: {pw.verdict.value_counts().to_dict()}")
        # Verdict counts and the mu<0 tally are censuses of this plate's wells rather than
        # estimates from them, so they carry no sampling interval; the medians below do.
        print(f"  negative-activity fraction  median "
              f"{well_interval(pw.negative_activity_fraction, '{:.2f}')}")
        k_deg = pw.implied_min_k_deg.median()
        print(f"  implied min k_deg (1/h)     median {well_interval(pw.implied_min_k_deg, '{:.4f}')}"
              f"  -> apparent half-life {np.log(2) / max(k_deg, 1e-9):.1f} h")
        print(f"  QSS-valid fraction          median {well_interval(pw.qss_fraction, '{:.2f}')}")
        print(f"  wells with declining OD (mu<0): {int((pw.mu_min < 0).sum())}/{len(pw)}")
        r2 = pw.dilution_r2.dropna()
        if len(r2):
            print(f"  dilution R^2 (reportable in {len(r2)}/{len(pw)} wells): "
                  f"median {well_interval(r2)}  IQR {r2.quantile(.25):.3f}-{r2.quantile(.75):.3f}")
        if len(clean):
            c2 = clean.dilution_r2.dropna()
            print(f"  -- restricted to wells with monotone OD and QSS>50% (n={len(clean)}):")
            if len(c2):
                print(f"     dilution R^2 median {well_interval(c2)}")
                print(f"     slope median {well_interval(clean.dilution_slope, '{:+.3f}')}")
            print(f"     negative-activity median "
                  f"{well_interval(clean.negative_activity_fraction, '{:.2f}')}")
            print(f"     implied k_deg (1/h) median "
                  f"{well_interval(clean.implied_min_k_deg, '{:.4f}')}")
        print("  'within-plate' bounds resample wells, which share an inoculum and a reader:")
        print("  they say how firmly this plate pins the number down, not how a second plate")
        print("  would come out. One plate is one biological replicate.")


if __name__ == "__main__":
    # This used to `resolve_or_exit` on YSTWIN_IGEM_RESULTS, so the script ran on exactly
    # one machine -- the one with the July workbooks in a sibling directory. Their numbers
    # are now committed as text under data/plates, verified value-for-value, so the replay
    # is the ordinary path and the workbooks are the optional one.
    results = paths.igem_results()
    exports = sorted(results.rglob("*.xlsx")) if results is not None else []
    if not exports:
        # The July set only, which is what this script read when it resolved one directory.
        # The replay covers every committed export, so an unfiltered call would quietly
        # widen this survey from two plates to seven and write five tracked tables nobody
        # asked for.
        replay.install(globals())
        exports = replay.exports(source_set="july_2026_07")
        print(f"replaying the committed text in {replay.committed_dir()} "
              f"({len(exports)} July exports); set YSTWIN_IGEM_RESULTS to read the workbooks")
    if not exports:
        # An empty survey table is indistinguishable from a plate that had nothing in it.
        raise SystemExit("no .xlsx exports found and no committed plate text either")
    for f in exports:
        try:
            analyse(f)
        except Exception as exc:  # noqa: BLE001 - survey script, report and continue
            print(f"\n{f.name}: skipped ({type(exc).__name__}: {str(exc)[:70]})")
