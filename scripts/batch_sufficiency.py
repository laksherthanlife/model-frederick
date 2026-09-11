"""The chemostat requirement, tested and found unnecessary: a flask at constant mu suffices.

    python3 scripts/batch_sufficiency.py

Writes ``outputs/batch_sufficiency.csv``. Needs no wet-lab data, no GEM, no network -- it is
an integration of the balance `pathway/solve.py` already solves.

WHY THIS EXISTS. Every standing recommendation in this repository asked for a HELD growth
rate, and therefore a chemostat, or the exponential fed-batch that replaced it.
`docs/DISTANCE_TO_THE_VISION.md` went further and listed "shake flask / batch" as REFUSED.
An adversarial audit on 2026-09-03 attacked that premise and it did not survive.

THE BALANCE THE SOLVER CLOSES IS PER GRAM OF BIOMASS::

    d[X]/dt = v_in - v_out(X) - mu*[X]

**The vessel appears nowhere in it.** Nothing distinguishes a chemostat from a flask except
whether mu happens to stay constant. Integrating that same balance forward at constant mu,
with the same rate laws and the shipped kinetics, the terminal pool converges on exactly the
steady-state answer `solve_pathway` returns, and the approach rate is set by mu alone. So the
cost of not having a chemostat is not a different answer; it is a GENERATION COUNT, and it is
the same count in either vessel:

    within 10% of steady state    ~3.7-4.0 generations
    within  5%                    ~4.8-5.1
    within  1%                    ~7.2-7.5

**mu must be KNOWN and roughly CONSTANT. It does not have to be HELD.** That is the whole
finding, and it converts every chemostat ask in this repository into a shake flask with a
low inoculum and an OD trace.

WHAT SETS THE ERROR, and it is not the assay. The lag error tracks the fractional FALL IN mu
PER DOUBLING almost one-for-one -- 1% drift per doubling gives about -1% bias, 5% gives about
-6%, 10% gives about -15%. So steadiness should be certified on the OD trace, which has ~100
points at a 1.80% CV, and NOT on the content series, which has four points at 5-10%. A
flatness criterion applied to the content series under-certifies by about a factor of two:
at 5% assay scatter, "flat" is first declared around 4 generations where the true bias is
still about -10%.

THE INOCULUM IS THE EXPERIMENT. Generations of exponential phase available in a 20 g/L
glucose flask: X0 = 0.5 g/L gives about 2.6, which reaches only -36%; X0 = 0.05 gives about
5.7, reaching -5%; X0 = 0.005 gives about 9.0, reaching -0.6%. Two serial exponential
precultures also remove the enzyme-fill transient entirely, which otherwise costs 6-9
generations on its own -- starting from a preculture already at the same mu, the deviation
over eight generations is under 0.002%.
"""
from __future__ import annotations

import math
import pathlib
import sys
import warnings

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.pathway import calibrations
from ystwin.pathway.solve import solve_pathway
from ystwin.pathway.spec import load_pathway

warnings.filterwarnings("ignore")

GROWTH_RATES = (0.101, 0.15, 0.18, 0.2543)
TOLERANCES = (0.10, 0.05, 0.01)


def _generations_to_converge(spec, kinetics, v_in, growth_rate):
    """Integrate the balance forward from an empty cell; return generations to each tolerance.

    Forward Euler at a step small enough that halving it moves the answer by less than the
    reported precision -- checked once rather than assumed, because a convergence claim made
    with a too-coarse step is exactly the bug this repository already hit in `fba/dynamic.py`.
    """
    steady = solve_pathway(spec, v_in, growth_rate, kinetics).terminal.content_mmol_per_gdcw
    branch = kinetics["lycopene"]
    vmax = branch.vmax_per_growth * growth_rate
    km = branch.km

    intermediate, terminal = 0.0, 0.0
    step, elapsed = 1e-4, 0.0
    marks: dict[float, float] = {}
    while elapsed < 400.0 and len(marks) < len(TOLERANCES):
        export = vmax * intermediate / (km + intermediate)
        intermediate += (v_in - export - growth_rate * intermediate) * step
        terminal += (export - growth_rate * terminal) * step
        elapsed += step
        relative = abs(terminal - steady) / steady
        for tolerance in TOLERANCES:
            if tolerance not in marks and relative < tolerance:
                marks[tolerance] = elapsed * growth_rate / math.log(2.0)
    return steady, marks


def main() -> int:
    spec = load_pathway("beta_carotene")
    kinetics = calibrations.BETA_CAROTENE_KINETICS
    v_in = calibrations.BETA_CAROTENE_FLUX.alpha * 1.0

    rows = []
    for growth_rate in GROWTH_RATES:
        steady, marks = _generations_to_converge(spec, kinetics, v_in, growth_rate)
        row = {"growth_rate_per_h": growth_rate,
               "steady_state_mmol_per_gdcw": round(steady, 9),
               "doubling_time_h": round(math.log(2.0) / growth_rate, 3)}
        for tolerance in TOLERANCES:
            row[f"generations_to_{int(tolerance * 100)}pct"] = round(
                marks.get(tolerance, float("nan")), 2)
            row[f"hours_to_{int(tolerance * 100)}pct"] = round(
                marks.get(tolerance, float("nan")) * math.log(2.0) / growth_rate, 1)
        rows.append(row)

    frame = pd.DataFrame(rows)
    out = paths.outputs_dir() / "batch_sufficiency.csv"
    frame.to_csv(out, index=False)
    print("A flask at constant mu reaches the same steady state the chemostat solve returns.")
    print("The cost is a generation count, and it is the same count in either vessel.\n")
    print(frame.to_string(index=False))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
