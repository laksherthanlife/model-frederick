"""D1: is beta-carotene flux determined by the GSMM, or merely bounded by it?

Run in a regime that has been checked against a measured phenotype, because a
product ceiling computed in an unvalidated regime describes a different organism.

Usage: python scripts/parked/run_d1.py
"""
from __future__ import annotations

import pathlib
import sys

import cobra

# parents[2], not parents[1]: this script lives one directory deeper than the others.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))

from ystwin import paths
from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID, add_beta_carotene_pathway
from ystwin.fba.fva import capacity_sweep, product_flux_range
from ystwin.fba.physiology import validate_physiology
from ystwin.fba.solver import growth_or_none

OUT = paths.outputs_dir()
MW = 536.87  # g/mol beta-carotene
FERMENTATIVE = ("growth_rate", "glucose", "ethanol")


def flux_from_content(mg_per_gdcw: float, mu: float) -> float:
    return (mg_per_gdcw * 1e-3 / MW) * 1e3 * mu


def banner(text: str) -> None:
    print("\n" + "=" * 78 + f"\n{text}\n" + "=" * 78)


def main() -> None:
    gem = paths.resolve_or_exit(paths.yeast_gem(), "the Yeast9 GSMM (yeast-GEM.xml)",
                                "YSTWIN_YEAST_GEM")
    ec_gem = paths.resolve_or_exit(paths.ec_yeast_gem(),
                                   "the GECKO model (ecYeastGEM_batch.xml)",
                                   "YSTWIN_EC_YEAST_GEM")
    plain = cobra.io.read_sbml_model(str(gem))
    ec = cobra.io.read_sbml_model(str(ec_gem))

    banner("D1.0  Which model may be used to bound anything?")
    print("Yeast9 v9.0.2, measured uptakes applied:")
    print(validate_physiology(plain).summary())
    print("\necYeastGEM v8.3.4 batch, no hand-set uptakes:")
    print(validate_physiology(ec, apply_constraints=False).summary())
    print("\n-> fermentative subset only, ec model:")
    print(validate_physiology(ec, apply_constraints=False, keys=FERMENTATIVE).summary())

    model = add_beta_carotene_pathway(ec)
    banner("D1.a  Feasible beta-carotene range vs enforced growth (ec model)")
    print(f"{'growth %':>9} {'mu (1/h)':>10} {'min':>12} {'max':>12} {'rel. width':>11}")
    for frac in (1.0, 0.99, 0.95, 0.90, 0.75, 0.50, 0.25):
        r = product_flux_range(model, PRODUCT_DEMAND_ID, glucose_uptake=None, growth_fraction=frac)
        print(f"{frac:>8.0%} {r.enforced_growth:>10.4f} {r.minimum:>12.3e} "
              f"{r.maximum:>12.3e} {r.relative_width:>11.3f}")

    r90 = product_flux_range(model, PRODUCT_DEMAND_ID, glucose_uptake=None, growth_fraction=0.90)
    mu = r90.enforced_growth
    banner("D1.b  Stoichiometric ceiling vs reported engineered titres")
    print(f"  ec ceiling at 90% growth (mu={mu:.4f}/h): {r90.maximum:.4e} mmol/gDW/h")
    print(f"  ec floor:                                 {r90.minimum:.4e} mmol/gDW/h")
    for mg in (5, 10, 20, 50):
        obs = flux_from_content(mg, mu)
        print(f"  {mg:>3} mg/gDCW -> {obs:.4e} mmol/gDW/h   ceiling/observed = {r90.maximum / obs:>6.1f}x")

    banner("D1.c  Effective pathway capacity is what carries the prediction")
    grid = [1, 5, 10, 20, 50, 100]
    sweep = capacity_sweep(
        model, PRODUCT_DEMAND_ID,
        capacities=[flux_from_content(mg, mu) for mg in grid],
        glucose_uptake=None, growth_fraction=0.90,
    )
    sweep.insert(0, "mg_per_gdcw", grid)
    print(sweep.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    sweep.to_csv(OUT / "d1_capacity_sweep_ec.csv", index=False)

    banner("D1.d  Does the pathway compete with growth at realistic rates?")
    base_growth = r90.max_growth
    for mg in (5, 20, 100, 500, 2000):
        q = flux_from_content(mg, mu)
        with model as m:
            m.reactions.get_by_id(PRODUCT_DEMAND_ID).lower_bound = q
            m.objective = "r_2111"
            g = growth_or_none(m)
        drop = "infeasible" if g is None else f"{g:.5f} /h  ({g / base_growth:6.2%} of max)"
        print(f"  forcing {mg:>5} mg/gDCW-equivalent ({q:.3e}): max growth {drop}")


if __name__ == "__main__":
    main()
