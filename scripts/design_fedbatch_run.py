"""Plan the fed-batch runs that replace this project's chemostat states.

    python3 scripts/design_fedbatch_run.py

Writes ``outputs/fedbatch_design.csv``. Needs no wet-lab data and no GEM -- the whole design
is closed-form, which is itself part of the case for the vessel.

WHAT THIS ANSWERS. `pathway/calibrations.py` fits its entry-flux scalar on six Elizondo
chemostat states: three strains at two growth rates. Replacing the chemostat means asking,
per state, whether an exponential fed-batch could produce it -- and the answer is not uniform.
The lower rate transfers. The upper rate does not, because it sits ON the strain's measured
uptake ceiling, and a feed cannot set a growth rate the strain is already maxed out at.

The refusals are the point of the table, so they are rows in it rather than omissions from it.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.fba.fedbatch import SetpointAboveCapacity, design_fedbatch

#: The two Elizondo growth rates, and the strain physiology measured at each. mu_max is taken
#: as the strain's own measured uptake ceiling expressed as a growth rate -- which is exactly
#: why the upper state has no headroom: at mu = 0.254 /h the culture IS at its ceiling.
ELIZONDO_STATES = [
    {"state": "lower", "mu_per_h": 0.101, "yield_g_per_g": 0.098,
     "max_uptake_g_per_gdcw_h": 2.59, "mu_max_per_h": 0.254},
    {"state": "upper", "mu_per_h": 0.254, "yield_g_per_g": 0.098,
     "max_uptake_g_per_gdcw_h": 2.59, "mu_max_per_h": 0.254},
]

#: A 2 L bench vessel, the size this project would actually buy. Sf 500 g/L is a standard
#: concentrated glucose feed.
VESSEL = {"initial_volume_l": 1.0, "max_volume_l": 2.0, "feed_substrate_g_per_l": 500.0}

#: Two inocula, because the inoculum is the only free lever on the generation budget:
#: generations = log2(1 + (Vmax-V0)*Y*Sf/(X0*V0)), and X0 is the sole term a protocol sets.
INOCULA = (0.5, 0.05)


def main() -> int:
    rows = []
    for state in ELIZONDO_STATES:
        for inoculum in INOCULA:
            row = {
                "state": state["state"],
                "target_mu_per_h": state["mu_per_h"],
                "inoculum_g_per_l": inoculum,
                "assumed_yield_g_per_g": state["yield_g_per_g"],
                "mu_max_per_h": state["mu_max_per_h"],
                "max_uptake_g_per_gdcw_h": state["max_uptake_g_per_gdcw_h"],
            }
            try:
                design = design_fedbatch(
                    state["mu_per_h"], state["mu_max_per_h"],
                    initial_biomass_g_per_l=inoculum,
                    assumed_yield_g_per_g=state["yield_g_per_g"],
                    max_uptake_g_per_gdcw_h=state["max_uptake_g_per_gdcw_h"],
                    **VESSEL)
            except SetpointAboveCapacity as refusal:
                row.update({
                    "feasible": False,
                    "refusal": str(refusal).split(".")[0],
                    "opening_uptake_g_per_gdcw_h": round(
                        state["mu_per_h"] / state["yield_g_per_g"], 3),
                })
                rows.append(row)
                continue

            budget = design.generation_budget_at(state["yield_g_per_g"])
            settling = design.settling_generations(2.0)
            row.update({
                "feasible": True,
                "refusal": "",
                "opening_uptake_g_per_gdcw_h": round(
                    design.initial_uptake_g_per_gdcw_h, 3),
                "initial_feed_ml_per_h": round(
                    1000 * design.initial_feed_rate_l_per_h, 3),
                "capacity_margin": round(design.capacity_margin, 3),
                "generation_budget": round(budget, 2),
                "settling_generations_at_2x_yield_error": round(settling, 2),
                "usable_generations": round(budget - settling, 2),
                "final_biomass_g_per_l": round(design.final_biomass_g_per_l, 1),
                "lowest_safe_assumed_yield_g_per_g": round(
                    design.worst_tolerable_yield_overestimate(
                        state["max_uptake_g_per_gdcw_h"]), 4),
            })
            rows.append(row)

    frame = pd.DataFrame(rows)
    out = paths.outputs_dir() / "fedbatch_design.csv"
    frame.to_csv(out, index=False)

    for _, row in frame.iterrows():
        if row["feasible"]:
            print(f"{row['state']:>6} mu={row['target_mu_per_h']:.3f} X0={row['inoculum_g_per_l']:.2f} "
                  f"-> F0 {row['initial_feed_ml_per_h']:6.3f} mL/h, "
                  f"{row['generation_budget']:.2f} gen budget, "
                  f"{row['usable_generations']:.2f} usable, "
                  f"X_final {row['final_biomass_g_per_l']:.1f} g/L")
        else:
            print(f"{row['state']:>6} mu={row['target_mu_per_h']:.3f} X0={row['inoculum_g_per_l']:.2f} "
                  f"-> REFUSED: {row['refusal'][:70]}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
