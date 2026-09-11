"""The degradation outlet, calibrated on real yeast — and the model it falsifies.

    python3 scripts/calibrate_glycogen.py

`solve_pathway` gained a third outlet on 2026-08-31: a node whose fate is `DEGRADED` loses
its pool to a native enzyme as well as to the next step and to growth. Glycogen is the case
that prompted it — GPH1 and SGA1 phosphorylase and glucosidase — and the solver refuses it
without a measured turnover rate. This script measures one.

**The falsification comes first, because it is the part that does not need a fit.** Boender
2009 (PMID 19592533) holds anaerobic glucose-limited *S. cerevisiae* at two growth rates: a
chemostat at mu = 0.025 /h and a retentostat where mu falls below 0.001 /h. If growth
dilution were the only outlet, `content = flux/mu` would put the retentostat pool **25x**
above the chemostat's. It is measured at **2.12x**. The two-outlet model is wrong here by an
order of magnitude, which is what the load-time refusal was protecting against and why
moving it to solve time had to keep refusing rather than start guessing.

**What the fit is and is not.** Two states and two unknowns is a saturated system: it
identifies `flux` and `k_deg` exactly and leaves no residual, so it cannot be wrong on its
own data. The only out-of-sample point available is a nitrogen-starved shake flask, which is
a different limitation entirely, and it is reported as a comparison rather than a test.
A third glucose-limited growth rate would make this a fit that could fail.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.pathway.solve import NodeKinetics, solve_pathway
from ystwin.pathway.spec import load_pathway

CONTENT = "glycogen_mmol_glucosyl_per_gdcw_as_glucose"


def states() -> pd.DataFrame:
    return pd.read_csv(paths.data_dir() / "glycogen" / "boender2009_storage_carbohydrates.tsv",
                       sep="\t").set_index("state")


def calibrate(frame: pd.DataFrame) -> tuple[float, float]:
    """Solve ``content = flux / (mu + k_deg)`` on the two glucose-limited states.

    Returns ``(flux, k_deg)``. Saturated by construction -- two equations, two unknowns.
    """
    chem_mu = float(frame.loc["chemostat_D0025", "mu_per_h"])
    ret_mu = float(frame.loc["retentostat_22d", "mu_upper_bound_per_h"])
    chem_c = float(frame.loc["chemostat_D0025", CONTENT])
    ret_c = float(frame.loc["retentostat_22d", CONTENT])
    ratio = ret_c / chem_c
    k_deg = (chem_mu - ratio * ret_mu) / (ratio - 1.0)
    return chem_c * (chem_mu + k_deg), k_deg


def main() -> int:
    frame = states()
    chem_mu = float(frame.loc["chemostat_D0025", "mu_per_h"])
    ret_mu = float(frame.loc["retentostat_22d", "mu_upper_bound_per_h"])
    chem_c = float(frame.loc["chemostat_D0025", CONTENT])
    ret_c = float(frame.loc["retentostat_22d", CONTENT])

    print("\nBoender 2009, anaerobic glucose-limited S. cerevisiae\n")
    print(f"  chemostat    mu = {chem_mu:<8g} glycogen {chem_c:.4f} mmol glucosyl/gDCW")
    print(f"  retentostat  mu <= {ret_mu:<6g} glycogen {ret_c:.4f}")
    print(f"\n  measured rise as growth stops        {ret_c / chem_c:6.2f}x")
    print(f"  what growth dilution alone predicts  {chem_mu / ret_mu:6.0f}x")
    print(f"  -> the two-outlet model is wrong here by {(chem_mu / ret_mu) / (ret_c / chem_c):.0f}x."
          "\n     A storage pool does not diverge when growth stops; an enzyme eats it.")

    flux, k_deg = calibrate(frame)
    print(f"\n  fitted third outlet:  k_deg = {k_deg:.4f} /h   flux = {flux:.5f} mmol/gDCW/h")
    print(f"  plateau at mu -> 0:   flux/k_deg = {flux / k_deg:.4f} mmol/gDCW")
    print("  (two states, two unknowns -- this identifies the pair, it cannot fail on them)")

    # Does the solver, given that rate, reproduce the states it was fitted on?
    spec = load_pathway("glycogen")
    kinetics = {"glycogen": NodeKinetics(degradation_rate_per_h=k_deg)}
    rows = []
    for label, mu, measured in (("chemostat", chem_mu, chem_c),
                                ("retentostat", ret_mu, ret_c)):
        solved = solve_pathway(spec, flux, mu, kinetics).terminal.content_mmol_per_gdcw
        rows.append({"state": label, "mu_per_h": mu, "measured": measured,
                     "solved": solved, "fold": solved / measured})
    n_state = float(frame.loc["n_starved_shake_flask", CONTENT])
    rows.append({"state": "N-starved flask", "mu_per_h": float("nan"),
                 "measured": n_state, "solved": flux / k_deg,
                 "fold": (flux / k_deg) / n_state})
    table = pd.DataFrame(rows)
    print("\n  solver against the data:")
    print(table.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\n  The first two are the fit restating itself. The third is a DIFFERENT")
    print("  limitation -- nitrogen starvation, not glucose -- so read it as a comparison,")
    print("  not a held-out test. A third glucose-limited growth rate would make this a")
    print("  fit that could fail.")

    out = paths.outputs_dir() / "glycogen_calibration.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
