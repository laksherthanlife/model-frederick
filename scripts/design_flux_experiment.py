"""How many strains does the chemostat run need? Simulation can answer this. It cannot
answer whether the law is true.

    python3 scripts/design_flux_experiment.py

**The distinction this script exists to hold.** `docs/SECOND_DATASET_HUNT.md` ends with a
wet-lab ask -- entry-enzyme expression on >= 5 producing strains at one dilution rate in
chemostat -- and the obvious question is why that cannot be simulated instead. It can be,
and doing so buys exactly one thing: the **design**. It buys none of the evidence.

The reason is the tier system this repository already grades every claim against. Simulating
a strain panel means drawing expression, pushing it through ``flux = alpha * expression``,
adding noise, and fitting ``flux = alpha * expression`` back. What comes out is a measurement
of the fitter, because the generator and the model are the same equation. That is **Tier 1**
-- "works in simulation, tested against our own simulator" -- which the README already
awards to every power and transfer table here. **Tier 3 is "predicted first, measured
after"**, and no amount of simulation crosses that line, because the thing being asked is
whether the *organism* obeys the law, and a simulator that was told it does cannot report
otherwise.

So this script assumes the law and asks the only question left: **if it is true, how big
does the experiment have to be to show it?** That is a real answer, it de-risks a chemostat
run, and it is the honest use of a simulator.

Two bars have to clear, and they are different:

* **the resolution bar** -- with n strains there are n! label assignments, so the smallest
  exact-permutation p is 1/n!. At three strains that is 0.167: this particular test cannot
  clear 0.05 however clean the data. For large panels, B sampled assignments instead give
  finite Monte Carlo resolution 1/(B + 1), reported separately from the theoretical floor.
* **the power bar** -- given the residual scatter actually observed, how often would a
  panel of n strains land a leave-one-strain-out skill that beats the permutation null.

The noise level is not invented: it is the leave-one-strain-out RMSE the real CrtE fit
produced, so the simulated panel is as noisy as the one measurement anybody has.
"""

from __future__ import annotations

import argparse
import itertools
import math
import pathlib
import sys
from numbers import Integral

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX
from ystwin.pathway.flux import fit_flux_law


#: Above this many label assignments, draw this many independent uniform permutations
#: with replacement instead of enumerating. 720 is 6!, so panels through six remain exact.
#: A sampled p includes the observed assignment and has resolution 1/(draws + 1), not 1/n!;
#: Monte Carlo error can still change a verdict near the chosen significance threshold.
_EXACT_PERMUTATION_LIMIT = 720


def permutation_floor(n_strains: int) -> float:
    """Theoretical exact-permutation floor, 1/n!, not the sampled-test resolution.

    Three strains floor at 0.167 for this permutation test. Log-factorial evaluation avoids
    constructing a huge integer; sufficiently small theoretical floors underflow to zero.
    The finite sampled-test resolution is reported separately and never underflows.
    """
    if isinstance(n_strains, bool) or not isinstance(n_strains, Integral) or n_strains < 1:
        raise ValueError("n_strains must be a positive integer")
    return math.exp(-math.lgamma(n_strains + 1))


def _panel(n_strains, alpha, sigma_log, spread, rng):
    """One simulated panel: expression spread over `spread`-fold, flux from the law."""
    expression = np.exp(rng.uniform(0.0, math.log(spread), n_strains))
    flux = alpha * expression * np.exp(rng.normal(0.0, sigma_log, n_strains))
    return expression, flux


def _skill(expression, flux, strains) -> float:
    fit = fit_flux_law(list(expression), list(flux), list(strains),
                       entry_enzyme="CrtE", source="simulated")
    return fit.loso_skill


def power_at(n_strains: int, alpha: float, sigma_log: float, spread: float,
             n_trials: int, rng, alpha_level: float = 0.05) -> dict:
    """Fraction of simulated panels whose permutation p clears ``alpha_level``.

    Small panels enumerate all label assignments, including the observed assignment once.
    Larger panels draw independent uniform permutations with replacement and use
    ``(1 + exceedances)/(1 + draws)``. The extra observation prevents zero Monte Carlo p
    values. A design whose actual p-value resolution exceeds the level is unreachable,
    even if its theoretical 1/n! floor is lower.
    """
    floor = permutation_floor(n_strains)
    if n_strains < 2:
        raise ValueError("power requires at least two strains")
    if isinstance(n_trials, bool) or not isinstance(n_trials, Integral) or n_trials < 1:
        raise ValueError("n_trials must be a positive integer")
    if not math.isfinite(alpha_level) or not 0 < alpha_level <= 1:
        raise ValueError("alpha_level must be finite and in (0, 1]")

    n_permutations, exact = 1, True
    for size in range(2, n_strains + 1):
        n_permutations *= size
        if n_permutations > _EXACT_PERMUTATION_LIMIT:
            n_permutations, exact = _EXACT_PERMUTATION_LIMIT, False
            break
    resolution = 1.0 / (n_permutations + (0 if exact else 1))
    if exact:
        floor = resolution
    method = "exact" if exact else "monte_carlo_with_replacement"
    result = {"n_strains": n_strains, "permutation_floor": floor,
              "permutation_method": method, "permutations_per_trial": n_permutations,
              "p_value_resolution": resolution}
    if resolution > alpha_level:
        return {**result, "power": 0.0, "median_skill": float("nan"),
                "verdict": f"UNREACHABLE: {method} resolution {resolution:.6g} > {alpha_level:g}"}

    strains = [f"s{i}" for i in range(n_strains)]
    cleared, skills = 0, []
    for _ in range(n_trials):
        expression, flux = _panel(n_strains, alpha, sigma_log, spread, rng)
        observed = _skill(expression, flux, strains)
        skills.append(observed)
        # Permute expression against flux: the null is "expression carries nothing".
        # Exact enumeration already includes the observed assignment. Above the bounded
        # count, sample directly without ever constructing the factorial assignment space.
        # The Monte Carlo count adds the observation to both numerator and denominator.
        orders = (itertools.permutations(range(n_strains)) if exact else
                  (rng.permutation(n_strains) for _ in range(n_permutations)))
        exceedances = sum(_skill(expression[list(order)], flux, strains) >= observed
                          for order in orders)
        p = (exceedances / n_permutations if exact else
             (exceedances + 1) / (n_permutations + 1))
        cleared += p <= alpha_level
    return {**result, "power": cleared / n_trials,
            "median_skill": float(np.median(skills)), "verdict": ""}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--spread", type=float, default=4.0,
                        help="fold range of entry-enzyme expression across the panel")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    alpha = BETA_CAROTENE_FLUX.alpha
    # The scatter the real fit produced, not a guess: this is how noisy the one measured
    # panel actually is, carried straight into the simulated one.
    sigma = BETA_CAROTENE_FLUX.loso_rmse_log
    rng = np.random.default_rng(args.seed)

    print(f"\nassuming the law is TRUE, at the measured scatter "
          f"(rmse_log = {sigma:.4f}, {math.exp(sigma):.2f}x typical)")
    print(f"expression spread across the panel: {args.spread:g}-fold, "
          f"{args.trials} simulated panels each\n")

    rows = [power_at(n, alpha, sigma, args.spread, args.trials, rng) for n in range(3, 9)]
    table = pd.DataFrame(rows)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}",
                          formatters={"permutation_floor": "{:.6g}".format,
                                      "p_value_resolution": "{:.6g}".format}))

    usable = table[(table.power >= 0.8) & (table.p_value_resolution <= 0.05)]
    print()
    if usable.empty:
        print("No panel size in this range reaches 80% power.")
    else:
        smallest = int(usable.n_strains.min())
        print(f"Smallest panel with >= 80% power at alpha = 0.05: **{smallest} strains**.")
        print("Three strains -- what exists -- cannot clear 0.05 at all: the permutation")
        print(f"floor is {permutation_floor(3):.3f}. Four is the first size that can "
              f"({permutation_floor(4):.3f}).")

    print("\nWhat this does NOT establish: that the organism obeys the law. The generator")
    print("and the model are the same equation, so this is Tier 1 by the README's own")
    print("grading. It sizes the experiment; it cannot replace it.")

    out = paths.outputs_dir() / "flux_experiment_design.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
