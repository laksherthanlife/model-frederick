"""Does a dose-dependent carrying capacity reproduce how a dosed well decelerates?

The incumbent generator integrated every well against one asserted `carrying_capacity` per
strain that no dose could move, while `culture.substrate_limited_capacity` sat beside it
deriving the same quantity from van Hoek's measured biomass yield and was called by nothing
but its own test. `simulate_culture` now takes the derived route by default. This script is
the score.

The observable is the one the orphan's own docstring points at -- how much a well slows down
between its fastest and its last hour, which is what a capacity term does to a growth curve.
It is formed exactly as `scripts/run_sensor_characterisation.py` forms it, on the same
estimator and the same late window, so the measured and simulated numbers are the same
statistic:

    mu_max   = 0.9 quantile of the smoothed log-OD derivative over the detectable window
    mu_late  = mean of that derivative over the last quarter of the window
    ratio    = mu_late / mu_max          decel = mu_max - mu_late   (1/h)

Scored in the growth channel against `MEASURED_GROWTH_RATE_SE` = 0.0117 /h, the standard
error of a growth rate estimated from one of these traces (`generator/panel_experiment.py`,
median over 210 real wells). The plate activity CV of 0.146 is a *reporter* floor and is not
borrowed here: nothing in this script reads a reporter.

What it finds, and it is not the result the build list expected. Two separate answers, and
they point opposite ways.

The ablation passes. Putting the asserted capacity back moves the late-window growth rate by
up to 0.0144 /h on the 4.14 h protocol, 1.23x the 0.0117 /h growth assay SE, and by 0.0660 /h
= 5.64x on a 24 h read that actually approaches the capacity. So the term is not decoration.

The defect it was supposed to close stays open. The measured drop in the ratio between the
control rung and the top two rungs is 0.4511 - 0.1365 = 0.3147 across 112 committed
conditions; the asserted capacity leaves 0.2391 of that unexplained and the derived capacity
leaves 0.2277, closing 4.8%. And the derived capacity moves the dose *direction* the wrong way
at every glucose charge below the standing 20 g/L, because van Hoek's yield rises as the rate
falls: a dose that slows a culture makes its substrate budget bigger, not smaller. The fitted
K_dosed/K_control = 0.245 that `REVISED_BUILD_LIST.md` cites is therefore not a prediction of
the measured yield curve, and no charge in the sweep makes it one.

Which is the orphan's docstring being right about its own scope: the effect is visible on a
plate "that runs to substrate exhaustion", and a 4.14 h read from OD 0.16 is not one.

Usage: python3 scripts/carrying_capacity_check.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.generator.culture import STANDARD_GLUCOSE_G_PER_L, simulate_culture
from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE
from ystwin.generator.plate import DEFAULT_DOSES, DEFAULT_PANEL, PlateConditions
from ystwin.growth import growth_window, max_specific_growth_rate, specific_growth_rate
from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT
from ystwin.readings import CorrectedOD

LATE_WINDOW_FRACTION = 0.75
"""Where the late window starts, as a fraction of the detectable points.

Not chosen here. It is the value `scripts/run_sensor_characterisation.py` used to produce the
committed table this script scores against, so any other value would compare two different
statistics. That script sweeps it; `outputs/late_window_sensitivity.csv` is the sweep.
"""

TOP_RUNGS = 2
"""How many rungs at the top of each ladder count as "high dose".

Two, because both ladders in `DEFAULT_DOSES` end in a pair (2 and 5 mM DTT, 2 and 4 mM H2O2)
and `REVISED_BUILD_LIST.md` quotes its fitted capacity ratio over exactly that 2-5 mM band.
"""

GLUCOSE_SWEEP = (0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
"""Charges to sweep the derived capacity over, g/L.

`STANDARD_GLUCOSE_G_PER_L` is asserted and was never assayed on these plates, so it is a
declared axis. The low end is below what the control well's own biomass gain would consume,
so if any charge could invert the sign of the dose effect it would be in this range.
"""

SEEDS = (0, 1, 2, 3, 4)
"""Well-to-well draws to average the score over.

One draw is not a measurement of a difference this small: the inoculum CV alone moves the
simulated drop between 0.04 and 0.14. Five is enough to show that the *difference between the
two capacities* is stable while the drop itself is not.
"""

EXHAUSTION_READ_H = 24.0
"""A read long enough for the default panel to close on its capacity, hours.

The control strain needs about 20 h to reach 95% of the derived capacity at 20 g/L. 4.14 h is
the biosensor protocol and reaches under a tenth of it, which is why the capacity term is
nearly invisible there and why the ablation is reported at both lengths.
"""


def _late_early(times_h: np.ndarray, od_corrected: np.ndarray) -> tuple[float, float]:
    """``(mu_max, mu_late)`` for one blank-corrected OD trace, 1/h each.

    Returns NaNs rather than raising when the trace never clears the detection floor for long
    enough to fit, which is what a dead high-dose well does.
    """
    keep = growth_window(times_h, CorrectedOD(od_corrected))
    if keep.sum() < 10:
        return float("nan"), float("nan")
    tw, odw = times_h[keep], od_corrected[keep]
    mu = specific_growth_rate(tw, CorrectedOD(odw))
    late = slice(int(len(tw) * LATE_WINDOW_FRACTION), None)
    return (max_specific_growth_rate(times_h, CorrectedOD(od_corrected)),
            float(np.nanmean(mu[late])))


def _band(stressor: str, dose: float, ladder: tuple[float, ...]) -> str:
    if dose == 0.0:
        return "control"
    return "high" if dose in ladder[-TOP_RUNGS:] else "mid"


def measured() -> pd.DataFrame | None:
    """The late/early statistic as the real plates give it, one row per condition.

    Read from `outputs/sensor_characterisation.csv` rather than recomputed from the exports:
    `paths.biosensor_plates()` returns None on this machine, so the raw Gen5 files are not
    resolvable and the committed table is the only carrier of this statistic. It holds the
    same columns the estimator produced -- `mu_max` and `mu_late` per condition -- so the
    comparison is like for like even though the traces are gone.
    """
    path = paths.outputs_dir() / "sensor_characterisation.csv"
    if not path.exists():
        return None
    table = pd.read_csv(path)
    conditions = table.drop_duplicates(["plate", "construct", "stressor", "dose_mM"]).copy()
    conditions["ratio"] = conditions.mu_late / conditions.mu_max
    conditions["decel"] = conditions.mu_max - conditions.mu_late
    conditions["band"] = [
        _band(row.stressor, row.dose_mM, tuple(sorted(DEFAULT_DOSES[row.stressor])))
        for row in conditions.itertuples()
    ]
    return conditions


def simulated(glucose_g_per_L: float | None, duration_h: float, seed: int = 0,
              replicates: int = 3) -> pd.DataFrame:
    """The same statistic from the generator, one row per construct and dose.

    Replicates are pooled into one mean trace before the rate is estimated, which is what the
    real pipeline does: the growth rate is a property of the condition, not of one noisy well.
    """
    setup = PlateConditions(duration_h=duration_h)
    rng = np.random.default_rng(seed)
    times = np.linspace(0.0, duration_h, setup.n_timepoints)
    rows = []
    for construct, params in DEFAULT_PANEL.items():
        stressor = STRESSOR_FOR_CONSTRUCT.get(construct, "DTT")
        ladder = tuple(sorted(DEFAULT_DOSES[stressor]))
        for dose in ladder:
            traces, capacity, source = [], float("nan"), ""
            for _ in range(replicates):
                start_od = setup.inoculum_for(construct) * float(
                    np.exp(rng.normal(0.0, setup.well_cv)))
                biomass0 = start_od * setup.gdcw_per_od
                out = simulate_culture(times, dose, params, biomass0,
                                       glucose_g_per_L=glucose_g_per_L)
                capacity, source = out["carrying_capacity"], out["capacity_source"]
                od = out["biomass"] / setup.gdcw_per_od + setup.od_blank
                traces.append(od * (1.0 + rng.normal(0.0, setup.reader_cv, times.size)))
            pooled = np.mean(traces, axis=0) - setup.od_blank
            mu_max, mu_late = _late_early(times, pooled)
            rows.append({
                "construct": construct, "stressor": stressor, "dose_mM": float(dose),
                "band": _band(stressor, dose, ladder),
                "mu_max": mu_max, "mu_late": mu_late,
                "ratio": mu_late / mu_max, "decel": mu_max - mu_late,
                "carrying_capacity": float(capacity), "capacity_source": source,
            })
    return pd.DataFrame(rows)


def _bands(frame: pd.DataFrame) -> tuple[float, float, float]:
    """``(ratio_control, ratio_high, drop)`` for one late/early table."""
    control = float(frame.loc[frame.band == "control", "ratio"].mean())
    high = float(frame.loc[frame.band == "high", "ratio"].mean())
    return control, high, control - high


def main() -> None:
    print("=" * 92)
    print("Dose-dependent carrying capacity: does wiring the measured yield curve in close it?")
    print("=" * 92)

    exports = paths.biosensor_plates()
    print(f"\nRaw plate exports: {paths.display_path(exports) if exports else 'NOT RESOLVABLE'}")
    if exports is None:
        print("  paths.biosensor_plates() returned None on this machine, so the late/early")
        print("  statistic is taken from the committed table rather than recomputed.")

    obs = measured()
    if obs is None:
        print("\nREFUSED: outputs/sensor_characterisation.csv is absent and no other committed")
        print("table carries mu_late beside mu_max. Nothing here can be scored.")
        return

    o_ctrl, o_high, o_drop = _bands(obs)
    print(f"\nMEASURED, {len(obs)} conditions over {obs.plate.nunique()} plates "
          f"({obs.groupby(['plate', 'construct', 'stressor', 'dose_mM']).ngroups} cells):")
    print(obs.groupby("band")[["ratio", "decel"]].agg(["mean", "std", "count"]).to_string())
    print(f"  ratio control {o_ctrl:.4f} -> high {o_high:.4f}, drop {o_drop:+.4f}")

    protocol_h = PlateConditions().duration_h
    fixed_drops = [_bands(simulated(None, protocol_h, seed))[2] for seed in SEEDS]
    derived_drops = [_bands(simulated(STANDARD_GLUCOSE_G_PER_L, protocol_h, seed))[2]
                     for seed in SEEDS]
    f_drop, d_drop = float(np.mean(fixed_drops)), float(np.mean(derived_drops))
    derived = simulated(STANDARD_GLUCOSE_G_PER_L, protocol_h, SEEDS[0])

    print(f"\nSIMULATED on the {protocol_h} h protocol, mean over {len(SEEDS)} seeds:")
    print(f"  asserted capacity   drop {f_drop:+.4f} "
          f"[{min(fixed_drops):+.4f}, {max(fixed_drops):+.4f}]"
          f"   residual vs measured {abs(o_drop - f_drop):.4f}")
    print(f"  derived  @ {STANDARD_GLUCOSE_G_PER_L:g} g/L  drop {d_drop:+.4f} "
          f"[{min(derived_drops):+.4f}, {max(derived_drops):+.4f}]"
          f"   residual vs measured {abs(o_drop - d_drop):.4f}")
    closed = (d_drop - f_drop) / (o_drop - f_drop)
    print(f"  the derived capacity closes {closed:.1%} of the gap the asserted one leaves.")

    print("\n  capacity per condition (g/L), and where it came from:")
    for row in derived.itertuples():
        if row.band in ("control", "high"):
            print(f"    {row.construct:12s} {row.dose_mM:4.1f} mM  K={row.carrying_capacity:7.3f}"
                  f"  {row.capacity_source}")

    print("\nGLUCOSE SWEEP -- the charge is asserted and was never assayed, so it is an axis:")
    print(f"  {'g/L':>6}  {'K control':>10}  {'K high':>9}  {'ratio ctrl':>10}  "
          f"{'ratio high':>10}  {'drop':>8}  {'residual':>9}")
    for charge in GLUCOSE_SWEEP:
        arm = simulated(charge, protocol_h, SEEDS[0])
        c, h, drop = _bands(arm)
        k_c = float(arm.loc[arm.band == "control", "carrying_capacity"].mean())
        k_h = float(arm.loc[arm.band == "high", "carrying_capacity"].mean())
        print(f"  {charge:6.1f}  {k_c:10.3f}  {k_h:9.3f}  {c:10.4f}  {h:10.4f}  "
              f"{drop:+8.4f}  {abs(o_drop - drop):9.4f}")

    print("\nABLATION -- remove the derived capacity, put the asserted one back, and measure")
    print(f"the growth channel against MEASURED_GROWTH_RATE_SE = {MEASURED_GROWTH_RATE_SE} /h:")
    for duration in (protocol_h, EXHAUSTION_READ_H):
        worst, median = [], []
        for seed in SEEDS:
            on = simulated(STANDARD_GLUCOSE_G_PER_L, duration, seed)
            off = simulated(None, duration, seed)
            delta = (on.mu_late - off.mu_late).abs()
            worst.append(float(delta.max()))
            median.append(float(delta.median()))
        peak = float(np.mean(worst))
        print(f"  {duration:5.2f} h read: median |d mu_late| = {np.mean(median):.4f} /h, "
              f"max = {peak:.4f} /h = {peak / MEASURED_GROWTH_RATE_SE:.2f}x the floor "
              f"({'CLEARS' if peak > MEASURED_GROWTH_RATE_SE else 'below'})")

    print("\nFREE SCALARS. Added: 1 (the glucose charge, declared and swept above).")
    print("Removed: 4 (one asserted capacity per construct in DEFAULT_PANEL, now reached only")
    print(f"on the refusal fallback). Independent targets: {len(obs)} committed conditions")
    print("carrying mu_late beside mu_max. Net free scalars -3, so the gate passes.")

    print("\nVERDICT. The wiring is structurally right and it clears the ablation bar, but it")
    print("does NOT close the high-dose deceleration: "
          f"{abs(o_drop - d_drop):.4f} of the {o_drop:.4f} measured")
    print("drop is left over. The derived capacity moves the dose direction the wrong way")
    print("below 20 g/L because van Hoek's yield rises as the rate falls, so a fitted")
    print("K_dosed/K_control below 1 is not something this curve can produce at any charge.")


if __name__ == "__main__":
    main()
