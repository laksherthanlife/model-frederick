"""Can the chain predict product amount for many products across many environments? Measure it.

    python3 scripts/product_environment_sweep.py

Writes ``outputs/product_environment_sweep.csv`` and
``outputs/product_environment_coverage.csv``. Needs no wet-lab data and no GEM.

This script exists because "environment conditions in, titre out, for any product" is the
project's stated goal, and nothing in the repository had ever run that sweep end to end and
reported what came back. Three findings, and none of them is a bug -- each is the
architecture working exactly as documented, made visible:

1.  ONE PRODUCT OF FOUR HAS A FITTED FLUX CALIBRATION. `data/pathways/` declares
    beta_carotene, gadusol, glycogen and phb. Only beta_carotene has a fitted flux calibration.
    Gadusol and PHB have no measured intermediate pool to constrain branch kinetics;
    glycogen additionally requires a measured terminal turnover rate. Loading those specs
    reports the missing measurements rather than borrowing carotenoid parameters.

2.  THE FIXED EMPIRICAL COMPARISON IS mu-ONLY AMONG REACHABLE HELD STATES. The setpoint
    must stay below the context's growth capacity even without a stressor, and every solve
    retains the kinetics' measured growth-rate window. Unreachable or out-of-range states
    return refusals, not content values. Equal content among the answered states describes
    the comparison's construction, not a validated lack of environmental effects.

3.  AND WIRING THE STRESS ROUTE WOULD NOT CHANGE THAT ENOUGH TO SEE. The one candidate
    channel is stress -> maintenance -> growth -> content, and its size is now measured
    (`outputs/maintenance_scale.csv`). Carried through to content it moves the answer by
    0.17-1.73% at the measured maintenance scale, against the flux calibration's own
    leave-one-strain-out error of 22.2%. **The signal is one to two orders of magnitude
    below the noise floor of the scalar it would have to pass through.** That is the finding
    that matters for planning, and it holds whichever way the unresolved sign goes.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.generator.context import CultureContext
from ystwin.pathway import calibrations
from ystwin.pathway.flux import FluxCalibration
from ystwin.pathway.spec import available_pathways, load_pathway
from ystwin.predict import Environment, Genotype, predict_product

#: The shipped lycopene branch parameters, including their measured growth-rate window.
#: `scripts/predict_product.py` refits this fixed empirical comparison leave-one-strain-out.
KINETICS = calibrations.BETA_CAROTENE_KINETICS

CARBON = ("glucose", "galactose", "ethanol")
TEMPERATURE_C = (25.0, 30.0, 37.0)
OXYGEN = (0.21, 0.05, 0.01)
STRESSORS = ((None, 0.0), ("DTT", 1.0), ("NaCl", 0.4), ("H2O2", 0.5))
HELD_RATE_PER_H = 0.18


def _coverage() -> pd.DataFrame:
    """Which declared pathways can be answered for at all, and why not where not."""
    rows = []
    for name in available_pathways():
        spec = load_pathway(name)
        matches = [k for k in dir(calibrations)
                   if k.isupper()
                   and isinstance(getattr(calibrations, k), FluxCalibration)
                   and getattr(calibrations, k).entry_enzyme == spec.entry_enzyme]
        rows.append({
            "product": name,
            "entry_enzyme": spec.entry_enzyme,
            "nodes": len(spec.nodes),
            "flux_calibration": matches[0] if matches else "",
            "answerable": bool(matches),
        })
    return pd.DataFrame(rows)


def main() -> int:
    coverage = _coverage()
    out_dir = paths.outputs_dir()
    coverage.to_csv(out_dir / "product_environment_coverage.csv", index=False)
    print(coverage.to_string(index=False))

    answerable = coverage[coverage.answerable]
    if answerable.empty:
        print("\nno product has a flux calibration; nothing to sweep", file=sys.stderr)
        return 1

    rows = []
    for product in answerable["product"]:
        spec = load_pathway(product)
        calibration = getattr(calibrations, coverage.set_index("product")
                              .loc[product, "flux_calibration"])
        genotype = Genotype(entry_expression=1.0)
        for held in (True, False):
            for carbon in CARBON:
                for temperature in TEMPERATURE_C:
                    for oxygen in OXYGEN:
                        for stressor, dose in STRESSORS:
                            context = CultureContext(carbon_source=carbon,
                                                     temperature_c=temperature,
                                                     oxygen=oxygen)
                            environment = Environment(
                                context=context, stressor=stressor, dose=dose,
                                growth_rate_setpoint_per_h=(HELD_RATE_PER_H if held
                                                            else None))
                            row = {"product": product, "mode": "held" if held else "batch",
                                   "prediction_mode": "empirical",
                                   "prediction_scope": "fixed-gene comparison, not environment validation",
                                   "carbon_source": carbon, "temperature_c": temperature,
                                   "oxygen": oxygen, "stressor": stressor or "",
                                   "dose": dose}
                            try:
                                got = predict_product(spec, genotype, environment,
                                                      calibration, KINETICS, mode="empirical")
                            except ValueError as refusal:
                                row.update({"refused": type(refusal).__name__,
                                            "refusal_reason": str(refusal),
                                            "growth_rate_per_h": None,
                                            "content_mg_per_gdcw": None})
                            else:
                                row.update({
                                    "refused": "", "refusal_reason": "",
                                    "growth_rate_per_h": round(got.growth_rate_per_h, 6),
                                    "content_mg_per_gdcw": round(
                                        got.content_mg_per_gdcw, 9)})
                            rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "product_environment_sweep.csv", index=False)

    for mode in ("held", "batch"):
        attempted = frame[frame["mode"] == mode]
        part = attempted[attempted.refused == ""]
        refusals = attempted[attempted.refused != ""]
        distinct = part.content_mg_per_gdcw.nunique()
        print(f"\n{mode:>5}: {len(part)} answered, {len(refusals)} refused out of "
              f"{len(attempted)} attempted environments -> {distinct} distinct content value(s)"
              f"; mu spans {part.growth_rate_per_h.nunique()} value(s)")
        for reason, count in refusals.refused.value_counts().items():
            print(f"       {reason}: {count}")
        if distinct == 1:
            print("       mu-only among answered fixed-gene comparisons; not environment validation")

    print(f"\nwrote {out_dir / 'product_environment_sweep.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
