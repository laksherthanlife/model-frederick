"""Descriptive error budget for the fixed-CrtE comparison.

    python3 scripts/error_budget.py

Rate, yield, content and titre share a relative error only when derived from the SAME
predicted flux using fixed measured denominators. Separate target refits need not share
that error. This script decomposes the fixed-CrtE leave-one-strain-out flux residuals;
it does not validate selecting CrtE or identify what future data will improve.

Three descriptive cuts, not guarantees about attainable accuracy:

1. **Between strain means vs within a strain.** This partitions mean squared log residuals,
   including bias. An independent, unobserved strain offset can impose an irreducible floor
   in a simulation, but that is an assumed fixture property, not a consequence of
   leave-one-strain-out validation. An approximately flat simulated curve is conditional
   on that noise model; more strains can still reduce estimation error. This decomposition
   alone cannot estimate either a learning curve or a biological noise floor.

2. **By dilution rate.** The law is fitted over mu in [0.101, 0.254] and the residuals are
   not uniform across it.

3. **Against a measurement-precision proxy.** Summing component interval endpoints is not
   a joint confidence interval for total flux. Comparing its mean fold width to a typical
   residual fold is descriptive, not a coverage test or evidence of measurement-limited
   accuracy.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.pathway.flux import fit_flux_law

ENTRY_GENE = "CrtE"


def states() -> pd.DataFrame:
    directory = paths.data_dir() / "carotenoid"
    frame = pd.read_csv(directory / "elizondo2025_steady_states.tsv", sep="\t").merge(
        pd.read_csv(directory / "elizondo2025_relative_mrna.tsv", sep="\t")
        .query("gene == @ENTRY_GENE")[["condition", "rel_expression"]], on="condition")
    frame["q_total"] = frame.q_lycopene + frame.q_betacarotene
    frame["rate"] = np.where(frame.mu_per_h > 0.2, "D=0.25", "D=0.10")
    # Sum-of-endpoints precision proxy, not a joint confidence interval on total flux.
    frame["ci_fold"] = np.sqrt(
        (frame.q_lycopene_hi95 + frame.q_betacarotene_hi95)
        / (frame.q_lycopene_lo95 + frame.q_betacarotene_lo95))
    return frame


def residuals(frame: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-strain-out log residual per state."""
    rows = []
    for held in sorted(frame.strain.unique()):
        train, test = frame[frame.strain != held], frame[frame.strain == held]
        fit = fit_flux_law(list(train.rel_expression), list(train.q_total),
                           list(train.strain), ENTRY_GENE, "loso")
        for r in test.itertuples():
            rows.append({"strain": r.strain, "rate": r.rate, "ci_fold": r.ci_fold,
                         "biomass_g_per_l": r.biomass_g_per_l,
                         "log_residual": float(np.log(fit.alpha * r.rel_expression
                                                      / r.q_total))})
    return pd.DataFrame(rows)


def budget(res: pd.DataFrame) -> dict:
    total = float(np.mean(res.log_residual ** 2))
    per_strain = res.groupby("strain").log_residual.mean()
    between = float(np.mean(per_strain.reindex(res.strain).to_numpy() ** 2))
    return {"total": total, "between_strain": between, "within_strain": total - between}


def main() -> int:
    frame = states()
    res = residuals(frame)
    parts = budget(res)

    print(f"\nfixed-CrtE comparison: typical leave-one-strain-out error "
          f"{np.exp(np.sqrt(parts['total'])):.3f}x\n")
    print("1. mean-squared log-residual decomposition")
    for label, key in (("between strain mean residuals", "between_strain"),
                       ("within strain (same strain, two rates)", "within_strain")):
        share = parts[key] / parts["total"]
        print(f"     {label:40s} {share:5.0%}")
    print("     -> descriptive split, not a learning curve. An independent, unobserved")
    print("        strain offset is a conditional simulation assumption, not an earned floor.")
    print("        More strains can still reduce estimation error under that assumption.")

    print("\n2. by dilution rate")
    for rate, group in res.groupby("rate"):
        print(f"     {rate}   model {np.exp(np.sqrt(np.mean(group.log_residual**2))):.3f}x"
              f"   precision proxy {group.ci_fold.mean():.3f}x"
              f"   biomass {group.biomass_g_per_l.mean():.3f} g/L")

    print("\n3. typical residual fold against the mean sum-of-endpoints precision proxy")
    for rate, group in res.groupby("rate"):
        model = float(np.exp(np.sqrt(np.mean(group.log_residual ** 2))))
        ci = float(group.ci_fold.mean())
        verdict = ("typical residual fold is below the mean precision proxy" if model <= ci else
                   "typical residual fold exceeds the mean precision proxy")
        print(f"     {rate}   {verdict}")
    print("     Descriptive only: not a coverage test, not a measurement-limit guarantee.")

    out = paths.outputs_dir() / "error_budget.csv"
    res.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
