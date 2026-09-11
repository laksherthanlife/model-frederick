"""End-to-end scenario predictions across conditions, attributed by branch.

Chain: mechanistic generator -> observed OD/RFU -> promoter activity ->
measured physiology -> FBA constraints -> flux. That third step is a smoothed point
estimate from ``reporter.promoter_activity``, not a filter: ``estimator.py`` is not on
this path. Branch B (latent) is solved beside
Branch A (physiology) so the difference is visible, and every row says which branch
produced it.

Usage: python scripts/run_scenarios.py
"""
from __future__ import annotations

import pathlib
import sys

import cobra
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.bridge.latent_bridge import LatentState, compare_branches, latent_constraints
from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints
from ystwin.generator.plate import DEFAULT_PANEL, DEFAULT_DOSES, PlateConditions, generate_plate
from ystwin.growth import max_specific_growth_rate, specific_growth_rate
from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT
from ystwin.reporter import ReporterKinetics, promoter_activity
from ystwin.readings import CorrectedOD, SpecificFluorescence

OUT = paths.outputs_dir()
REACTIONS = ["r_2111", "r_1714", "r_1992", "r_1672", "r_1761"]
NAMES = {"r_2111": "growth", "r_1714": "glucose", "r_1992": "oxygen",
         "r_1672": "CO2", "r_1761": "ethanol"}

# The ATP cost of the latent branch is NOT set here any more. This script used to pin
# `maintenance_per_activity = 1200.0` -- "a plausible scale, not a fitted one" -- and that
# parameter was removed from `latent_constraints` when it was found to be dimensionally
# incoherent: it multiplied a promoter activity in RFU/OD/h and called the product
# mmol ATP/gDW/h, demanding an NGAM around 3e5 against a resting value of 0.7, when Yeast9
# goes infeasible at 19.19. The replacement is a bounded scale (`_MAX_STRESS_MAINTENANCE`)
# over a measured activity range (`_ACTIVITY_FULL_SCALE`), and this script takes both as
# the library's defaults so there is one place they can be wrong.


def _corrected_activity(plate, conditions, t, construct: str, index: int) -> float:
    """Dilution-corrected promoter activity at one rung, RFU/OD/h, autofluorescence removed.

    Pulled out because the unstressed rung is now read twice -- once as the baseline the
    excess is measured from and once as a scenario in its own right -- and two copies of
    this arithmetic would be two chances to correct different things.
    """
    wells = [w for w in conditions.layout.wells(construct, index) if w in plate.od.columns]
    if not wells:
        return 0.0
    od = (plate.od[wells].mean(axis=1) - conditions.od_blank).to_numpy()
    rfu = plate.rfu[wells].mean(axis=1).to_numpy()
    specific = (rfu - conditions.optics.background) / od
    mu_series = specific_growth_rate(t, CorrectedOD(od))
    activity = promoter_activity(t, SpecificFluorescence(specific), mu_series, ReporterKinetics(k_deg=0.0))
    late = slice(int(len(t) * 0.6), None)
    autofluorescence_term = (
        float(np.nanmedian(mu_series[late])) * conditions.gdcw_per_od
        * conditions.optics.autofluorescence
    )
    return max(float(np.nanmedian(activity[late])) - autofluorescence_term, 0.0)


def main() -> None:
    gem = paths.resolve_or_exit(paths.yeast_gem(), "the Yeast9 GSMM (yeast-GEM.xml)",
                                "YSTWIN_YEAST_GEM")
    conditions = PlateConditions(seed=42)
    plate = generate_plate(DEFAULT_PANEL, conditions)
    model = cobra.io.read_sbml_model(str(gem))
    t = plate.times_h

    rows = []
    for construct in DEFAULT_PANEL:
        stressor = STRESSOR_FOR_CONSTRUCT[construct]
        # The unstressed rung's own corrected activity is what "excess" is measured from,
        # so it has to be in hand before the first constrained solve. Every ladder starts
        # at 0 mM -- there is a test that says so -- and this is the one place that fact
        # is load-bearing rather than decorative.
        basal = _corrected_activity(plate, conditions, t, construct, 0)
        for index, dose in enumerate(DEFAULT_DOSES[stressor]):
            wells = [w for w in conditions.layout.wells(construct, index) if w in plate.od.columns]
            if not wells:
                continue
            od = (plate.od[wells].mean(axis=1) - conditions.od_blank).to_numpy()
            rfu = plate.rfu[wells].mean(axis=1).to_numpy()
            specific = (rfu - conditions.optics.background) / od

            mu_series = specific_growth_rate(t, CorrectedOD(od))
            mu = max_specific_growth_rate(t, CorrectedOD(od))
            activity = promoter_activity(t, SpecificFluorescence(specific), mu_series, ReporterKinetics(k_deg=0.0))
            late = slice(int(len(t) * 0.6), None)
            activity_late = float(np.nanmedian(activity[late]))
            biomass = float(od[-1] * conditions.gdcw_per_od)

            measured = MeasuredState(
                biomass_gl=biomass, growth_rate=float(np.nanmedian(mu_series[late])),
                glucose_mM=20.0, time_h=float(t[-1]),
            )
            physical = physiology_constraints(measured)
            # Recovered activity is gdcw*gain*k_synth + mu*gdcw*autofluorescence. Both
            # terms scale with biomass, so the autofluorescence one is subtracted here,
            # and the result is LEFT IN RFU/OD/h.
            #
            # It used to be divided through by gdcw_per_od * gain as well, which turned it
            # into a synthesis rate around 1e-3. `latent_constraints` normalises by
            # `_ACTIVITY_FULL_SCALE = 1305 RFU/OD/h`, so a 1e-3 input made the stress
            # fraction about 3e-7 and the latent NGAM 0.7000002 against a resting 0.7 --
            # the branch was arithmetically incapable of differing from the physiology
            # branch, and the script's closing line duly reported "shifts by -0.0% to
            # +0.0%" as though that were a finding about biology. It was a unit mismatch,
            # left behind when `maintenance_per_activity` was replaced by a normalised
            # scale and the caller was not moved with it.
            autofluorescence_term = (
                measured.growth_rate * conditions.gdcw_per_od
                * conditions.optics.autofluorescence
            )
            corrected = max(activity_late - autofluorescence_term, 0.0)
            latent = latent_constraints(
                LatentState(corrected, construct, float(t[-1])),
                physical, acknowledge_unvalidated=True, basal_activity=basal,
            )
            try:
                table = compare_branches(model, physical, latent, reactions=REACTIONS)
            except ValueError as exc:
                rows.append({"construct": construct, "stressor": stressor, "dose_mM": dose,
                             "branch": "physiology", "note": str(exc)[:60]})
                continue
            for _, r in table.iterrows():
                rows.append({
                    "construct": construct, "stressor": stressor, "dose_mM": dose,
                    "branch": r.branch, "load_bearing": r.load_bearing,
                    "mu_max": mu, "growth_rate": measured.growth_rate,
                    "biomass_gl": biomass, "promoter_activity": activity_late,
                    "atp_maintenance": r.atp_maintenance,
                    **{NAMES[k]: r[f"flux_{k}"] for k in REACTIONS},
                })

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "scenario_predictions.csv", index=False)

    print("=" * 96)
    print("Scenario predictions.  Branch A = measured physiology (load-bearing).")
    print("                       Branch B = latent (research only; G4 refuses all modules).")
    print("=" * 96)
    for construct, group in frame.groupby("construct"):
        stressor = group.stressor.iloc[0]
        print(f"\n{construct}  ({stressor})")
        print(f"  {'dose':>5} {'mu':>7} {'promoter':>10} | "
              f"{'glucose A':>10} {'glucose B':>10} | {'oxygen A':>9} {'oxygen B':>9} | {'NGAM B':>7}")
        for dose, sub in group.groupby("dose_mM"):
            a = sub[sub.branch == "physiology"]
            b = sub[sub.branch == "latent"]
            if a.empty or b.empty:
                continue
            a, b = a.iloc[0], b.iloc[0]
            print(f"  {dose:>5.1f} {a.mu_max:>7.3f} {a.promoter_activity:>10.0f} | "
                  f"{a.glucose:>10.3f} {b.glucose:>10.3f} | "
                  f"{a.oxygen:>9.3f} {b.oxygen:>9.3f} | {b.atp_maintenance:>7.2f}")

    print("\n" + "=" * 96)
    print("Difference the latent branch makes, if it were ever admitted")
    print("=" * 96)
    wide = frame.pivot_table(index=["construct", "dose_mM"], columns="branch",
                             values=["glucose", "oxygen"])
    delta = ((wide[("glucose", "latent")] - wide[("glucose", "physiology")])
             / wide[("glucose", "physiology")].abs())
    print(f"  glucose uptake shifts by {delta.min():+.1%} to {delta.max():+.1%} across scenarios")
    print(f"  median shift {delta.median():+.1%}")
    print("\n  Every latent number above is REFUSED by G4 and must not be quoted as a result.")


if __name__ == "__main__":
    main()
