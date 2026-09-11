"""How much product a stress response costs, at a held growth rate. The channel that is big.

    python3 scripts/stress_diversion_scan.py

Writes ``outputs/stress_diversion.csv``. Needs the Yeast9 GSMM and nothing else.

Six routes from a stress state into the GEM were measured and all failed, the largest moving
titre 0.17-1.73% against a 22.2% noise floor. Every one was a SCALAR TAX -- an ATP demand, a
shrunken protein budget. This scans a different shape: flux REQUIRED into the named metabolic
effectors of the stress response, which is what a regulator actually does and what the GEM
has no reason to do on its own.

Yeast9 carries those effectors in full and none of the signalling that drives them, so the
missing decision is precisely the learned latent state. This table is the transfer function
from that decision to the product -- a property of the host, and therefore the half that
generalises. The other half, how much effector flux a given environment causes, is unmeasured
here and `bridge/stress_diversion.py` refuses to invent it.

Growth is HELD throughout. Every earlier stress route reached the product through mu and was
indistinguishable from the growth-rate response `generator/stress_panel.py` already measures
directly; holding mu makes diversion the only thing that can move the answer.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.bridge.stress_diversion import STRESS_EFFECTORS, diversion_cost
from ystwin.pathway.spec import available_pathways, load_pathway

warnings.filterwarnings("ignore")

FLUXES = (0.25, 0.5, 1.0, 2.0, 5.0)
GROWTH_RATE = 0.18


def main() -> int:
    path = paths.yeast_gem()
    if path is None:
        print("yeast-GEM not present", file=sys.stderr)
        return 1

    import cobra

    model = cobra.io.read_sbml_model(str(path))
    rows = []
    for product in available_pathways():
        spec = load_pathway(product)
        for effector in STRESS_EFFECTORS:
            for flux in FLUXES:
                result = diversion_cost(model, spec, effector, flux, GROWTH_RATE)
                rows.append({
                    "product": product,
                    "precursor": spec.precursor_metabolite,
                    "effector": effector,
                    "flux_mmol_per_gdcw_h": flux,
                    "growth_rate_per_h": GROWTH_RATE,
                    "precursor_supply": round(result.precursor_supply_mmol_per_gdcw_h, 6)
                    if result.precursor_supply_mmol_per_gdcw_h ==
                    result.precursor_supply_mmol_per_gdcw_h else None,
                    "relative_change": round(result.relative_change, 4)
                    if result.relative_change == result.relative_change else None,
                    "clears_22pct_floor": (bool(result.clears_the_calibration_floor)
                                           if result.relative_change == result.relative_change
                                           else None),
                    "conditions": result.conditions})

    frame = pd.DataFrame(rows)
    out = paths.outputs_dir() / "stress_diversion.csv"
    frame.to_csv(out, index=False)

    pivot = frame[frame["product"] == "beta_carotene"].pivot(
        index="effector", columns="flux_mmol_per_gdcw_h", values="relative_change")
    print(f"beta-carotene precursor supply, mu HELD at {GROWTH_RATE} /h "
          f"(relative change; NaN = infeasible)")
    print(pivot.to_string())
    cleared = frame.dropna(subset=["clears_22pct_floor"])
    print(f"\n{int(cleared.clears_22pct_floor.sum())} of {len(cleared)} "
          f"(product, effector, flux) combinations clear the 22.2% calibration floor")
    print("for comparison, the NGAM channel moved titre 0.17-1.73%")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
