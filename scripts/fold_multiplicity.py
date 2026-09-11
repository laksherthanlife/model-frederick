"""Report per-dose folds and sensor-level trends with their distinct inferential methods.

    python3 scripts/fold_multiplicity.py --output-dir /tmp/fold-multiplicity

By default the input is ``outputs/sensor_characterisation.csv``; ``--readings`` selects
another table independently of the output destination. Every observed positive-dose design
cell stays in the family, including cells whose estimate is refused.

Intervals use ``analysis.uncertainty.fold_change``: positive contributing plates, one summary
per well, and two-stage plate/well resampling when well identifiers exist. Bonferroni bounds
are nominal, approximate simultaneous intervals, not a coverage guarantee. The cluster-only
atom policy does not apply to two-stage resampling; finite Monte Carlo tail resolution and
the shared bootstrap's refusal rules still do. Resolution is checked for each cell.

An empirical bootstrap tail around the fitted distribution is not a confirmatory p-value.
It is neither exported as one nor fed into a multiple-testing correction. Per-dose Holm
and BH columns instead adjust two-sided exact sign tests of a median fold of one, using
non-tied positive per-plate folds. That median null differs from the geometric-mean target
of the interval, and exclusions condition both analyses on usable plates. The sign-test
floor bounds that test only, not all inference at the same plate count. Holm controls FWER
for valid marginal p-values; BH additionally needs independence or suitable positive
dependence, which this shared-control panel does not establish.

The separate sensor family uses within-plate dose-label permutation p-values and actual
Holm and BH corrections across sensors. Those tests require exchangeable dose assignments
under the null and a prespecified dose cutoff; they do not establish a particular fold size.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.multiplicity import (
    benjamini_hochberg,
    dose_response_permutation,
    finest_resolvable_alpha,
    holm,
    sign_test_floor,
)
from ystwin.analysis.uncertainty import FoldChange, fold_change

#: Above this the stressor is killing the culture and activity falls, so a monotone trend
#: test over the whole ladder would be looking for a hump. 1.0 mM is the top of the window
#: README section 1 already identifies as the one where the sensors work.
RESOLVABLE_MAX_DOSE_MM = 1.0

N_RESAMPLES = 4000
N_PERMUTATIONS = 20000


def _interval(readings, construct, dose, alpha, n_resamples, seed) -> FoldChange:
    """Use the shared estimator, with method-specific reporting-resolution policies."""
    estimate = fold_change(readings, construct, dose, n_resamples=n_resamples,
                           confidence=1.0 - alpha, seed=seed)
    if not estimate.estimable:
        return estimate
    reason = ""
    if (estimate.resampling == "cluster_only"
            and alpha < finest_resolvable_alpha(estimate.n_plates_contributing)):
        reason = "refused: requested alpha is finer than the cluster-only tail-resolution policy"
    elif estimate.n_resamples * alpha / 2.0 < 1:
        reason = "refused: too few Monte Carlo draws to resolve both requested tails"
    if reason:
        return replace(estimate, low=float("nan"), high=float("nan"),
                       estimable=False, method=reason)
    return estimate


def _sign_summary(readings, construct, dose) -> dict:
    """Obtain each plate's point from the shared estimator; exact ties carry no sign."""
    subset = readings[readings.construct == construct]
    dosed = subset[np.isclose(subset.dose_mM, dose)]
    control = subset[np.isclose(subset.dose_mM, 0.0)]
    shared = sorted(set(dosed.plate) & set(control.plate))
    points = np.array([
        fold_change(subset[subset.plate == plate], construct, dose, n_resamples=2).point
        for plate in shared
    ])
    usable = points[np.isfinite(points) & (points > 0)]
    above, below = int((usable > 1.0).sum()), int((usable < 1.0).sum())
    n_sign = above + below
    p_value = (float(binomtest(above, n_sign, p=0.5).pvalue) if n_sign
               else 1.0 if usable.size else float("nan"))
    floor = sign_test_floor(n_sign) if n_sign else 1.0 if usable.size else float("nan")
    return {"plates_above_one": above, "plates_below_one": below,
            "plates_tied_one": int(usable.size - n_sign), "n_plates_sign": n_sign,
            "p_exact_sign": p_value, "exact_sign_floor": floor}


def analyse(readings: pd.DataFrame, *, n_resamples: int = N_RESAMPLES,
            n_permutations: int = N_PERMUTATIONS, seed: int = 0) -> pd.DataFrame:
    """Keep each planned cell and correct its family without substituting bootstrap tails.

    ``estimable`` denotes nominal interval availability for fold rows and statistic
    availability for permutation rows; family-wise interval flags apply to fold rows only.
    """
    required = {"construct", "dose_mM", "plate", "activity_late"}
    if required - set(readings):
        raise ValueError(f"readings is missing required columns: {sorted(required - set(readings))}")
    if readings[["construct", "dose_mM", "plate"]].isna().any().any():
        raise ValueError("construct, dose and plate identifiers must not be missing")
    if not np.isfinite(readings.dose_mM).all():
        raise ValueError("dose_mM must be finite")
    if "well" in readings and readings.well.isna().any():
        raise ValueError("well identifiers must not be missing when wells are resampled")
    for name, value, minimum in (("n_resamples", n_resamples, 2),
                                 ("n_permutations", n_permutations, 1)):
        if (not isinstance(value, (int, np.integer)) or isinstance(value, bool)
                or value < minimum):
            raise ValueError(f"{name} must be an integer of at least {minimum}")
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, bool) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    rows = []

    # The FAMILY is every construct x dose cell the panel asks about, not the subset that
    # happened to return an interval. Until 2026-08-30 the Bonferroni divisor was
    # `len(table)` -- the estimable folds -- so a fold that could not be estimated made the
    # correction on all the others WEAKER, and the headline "one fold survives its whole
    # family" was produced by dividing alpha by 20 against a family of 24. Found by the
    # 2026-08-30 audit. A family size must be a property of the design, not of the data.
    cells = readings[readings.dose_mM > 0].groupby(["construct", "dose_mM"])
    family = cells.ngroups
    if not family:
        raise ValueError("readings must contain at least one positive-dose design cell")
    family_alpha = 0.05 / family

    for (construct, dose), dosed_rows in cells:
        try:
            estimate = _interval(readings, construct, dose, 0.05, n_resamples, seed)
            simultaneous = (_interval(readings, construct, dose, family_alpha, n_resamples, seed)
                            if family > 1 else estimate)
            signs = _sign_summary(readings, construct, dose)
        except ValueError as exc:
            estimate = FoldChange(construct, float(dose), float("nan"), float("nan"),
                                  float("nan"), f"refused: {exc}", 0, 0, 0, 0, False,
                                  resampling="two_stage" if "well" in readings else "cluster_only")
            simultaneous = replace(estimate, confidence=1.0 - family_alpha)
            signs = {"plates_above_one": 0, "plates_below_one": 0, "plates_tied_one": 0,
                     "n_plates_sign": 0, "p_exact_sign": float("nan"),
                     "exact_sign_floor": float("nan")}
        rows.append({
            "construct": construct, "dose_mM": dose,
            "stressor": dosed_rows.stressor.iloc[0] if "stressor" in dosed_rows else "",
            "measurement": "per_dose_fold", "n_plates": estimate.n_plates,
            "n_plates_contributing": estimate.n_plates_contributing,
            "n_wells": estimate.n_wells, "fold": estimate.point,
            # The empirical bootstrap tail is not a confirmatory p-value and is not corrected.
            # Exact sign-test p-values use only finite positive per-plate folds,
            # excluding exact ties; their null is a plate-level median fold of one.
            **signs,
            "p_value_method": "two-sided exact sign test of median usable plate fold = 1; ties excluded",
            "low": estimate.low, "high": estimate.high, "estimable": estimate.estimable,
            "clears": estimate.excludes_unity, "confidence": estimate.confidence,
            "interval_method": estimate.method, "resampling": estimate.resampling,
            "n_resamples": estimate.n_resamples, "target": estimate.target,
            "coverage_status": estimate.coverage_status,
            "assumptions": "; ".join(estimate.assumptions),
            "fw_low": simultaneous.low, "fw_high": simultaneous.high,
            "family_wise_resolvable": simultaneous.estimable,
            "clears_family_wise": simultaneous.excludes_unity,
            "family_alpha": family_alpha, "family_wise_method": simultaneous.method,
            "fw_n_resamples": simultaneous.n_resamples,
        })

    table = pd.DataFrame(rows)
    table["family_size"] = family  # set by the planned cells, including refusals
    table["p_holm"] = holm(list(table.p_exact_sign), family_size=family)
    table["p_bh"] = benjamini_hochberg(list(table.p_exact_sign), family_size=family)
    table = table.sort_values("fold", ascending=False, kind="stable")

    trends = []
    constructs = sorted(readings.construct.unique())
    for construct in constructs:
        subset = readings[(readings.construct == construct)
                          & (readings.dose_mM <= RESOLVABLE_MAX_DOSE_MM)]
        row = {"construct": construct, "measurement": "dose_response_permutation",
               "stressor": subset.stressor.iloc[0] if "stressor" in subset and len(subset) else "",
               "family_size": len(constructs), "n_doses": subset.dose_mM.nunique(),
               "n_plates": subset.plate.nunique(), "n_informative_plates": 0,
               "statistic": float("nan"), "p_permutation": float("nan"),
               "n_permutations": n_permutations, "estimable": False,
               "clears": False, "clears_family_wise": False, "family_wise_resolvable": False}
        try:
            got = dose_response_permutation(readings, construct,
                                            max_dose=RESOLVABLE_MAX_DOSE_MM,
                                            n_permutations=n_permutations, seed=seed)
            row.update({"stressor": got.stressor, "n_doses": got.n_doses,
                        "n_plates": got.n_plates, "n_informative_plates": got.n_informative_plates,
                        "statistic": got.statistic, "p_permutation": got.p_value,
                        "estimable": True,
                        "n_permutations": got.n_permutations, "p_value_method": got.method,
                        "assumptions": got.exchangeability})
        except ValueError as exc:
            row["p_value_method"] = f"refused: {exc}"
        trends.append(row)
    trends = pd.DataFrame(trends)
    trends["p_holm"] = holm(list(trends.p_permutation), family_size=len(constructs))
    trends["p_bh"] = benjamini_hochberg(list(trends.p_permutation), family_size=len(constructs))
    return pd.DataFrame([*table.to_dict("records"), *trends.to_dict("records")])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readings", type=pathlib.Path,
                        help="input sensor-characterisation CSV; independent of --output-dir")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="destination for fold_multiplicity.csv (default: YSTWIN_OUTPUTS or outputs)")
    parser.add_argument("--resamples", type=int, default=N_RESAMPLES)
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    source = args.readings if args.readings is not None else paths.outputs_dir() / "sensor_characterisation.csv"
    table = analyse(pd.read_csv(source), n_resamples=args.resamples,
                    n_permutations=args.permutations, seed=args.seed)
    folds = table[table.measurement == "per_dose_fold"]
    family = int(folds.family_size.iloc[0])
    print(f"{int(folds.estimable.sum())} estimable intervals among {family} planned per-dose cells.")
    print(f"Family-wise alpha is 0.05/{family} = {0.05 / family:.5f}; refused cells remain in the family.")
    print(folds[["construct", "dose_mM", "fold", "n_plates_contributing", "n_plates_sign",
                 "p_exact_sign", "p_holm", "p_bh", "low", "high", "fw_low", "fw_high"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nIntervals clearing 1.0: unadjusted {int(folds.clears.sum())}; "
          f"family-wise {int(folds.clears_family_wise.sum())}.")
    print("Bonferroni simultaneous bounds are nominal and approximate; coverage is not guaranteed.")
    for row in folds[~folds.family_wise_resolvable].itertuples():
        # A missing family-wise interval is not an absence-of-response result.
        # The reason may be too few usable plates, unstable ratios, insufficient
        # Monte Carlo tail draws, or the cluster-only policy when wells are absent.
        # No such refusal licenses a universal resolution bound for other methods.
        print(f"  {row.construct} @ {row.dose_mM:g}: {row.family_wise_method}")
    print("The empirical bootstrap tail is not a confirmatory p-value; it is not exported or corrected.")
    print("Exact sign-test floors use each row's non-tied contributing plates and are "
          "not a resolution bound for bootstrap, parametric or within-plate permutation inference.")
    print("The sign test concerns a median usable-plate fold; the interval targets a geometric mean.")
    print(f"Exact-sign Holm/BH calls at 0.05: {int((folds.p_holm <= 0.05).sum())}/"
          f"{int((folds.p_bh <= 0.05).sum())}. BH requires an additional dependence assumption.")
    print(f"\nOne permutation test per sensor at cutoff <= {RESOLVABLE_MAX_DOSE_MM:g} mM. "
          "Confirmatory interpretation requires exchangeable within-plate dose labels "
          "and a prespecified cutoff.")
    print(table[table.measurement == "dose_response_permutation"]
          [["construct", "statistic", "p_permutation", "p_holm", "p_bh", "p_value_method"]]
          .to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    directory = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    directory.mkdir(parents=True, exist_ok=True)
    out = directory / "fold_multiplicity.csv"
    table.to_csv(out, index=False)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
