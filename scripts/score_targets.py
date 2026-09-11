"""Compare separate fixed-CrtE target refits, not the error cost of converting a rate.

    python3 scripts/score_targets.py

Each row fits ``target = alpha * expression`` independently. A content-proportional law
and a rate-proportional law are different hypotheses, not two units for the same fitted
prediction. Their leave-one-strain-out skills use training-only geometric-mean baselines.
This is an exploratory fixed-gene comparison on the Elizondo states: neither CrtE nor the
winning target was prespecified independently of these data, and ranking them here is not
nested validation or evidence for a causal cassette-dosage mechanism.

Derive every target from the SAME rate prediction using the SAME measured growth rate,
uptake and biomass on both sides of the comparison, and relative and log errors are
identical: ``(q_hat / mu) / (q_measured / mu) = q_hat / q_measured``. Common factors cancel
exactly. This identity does not remove uncertainty in a new, unmeasured denominator or
biomass estimate, nor does it describe separately refitted target laws.

The permutation null reassigns expression profiles across strains, keeping their paired
dilution-rate states together. Renaming strain labels alone would leave the validation
partition unchanged. The exact floor is set by the number of profile assignments; the
printed ranks are conditional on the fixed gene and target, not on selecting them.
"""

from __future__ import annotations

import itertools
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.pathway.flux import fit_flux_law

ENTRY_GENE = "CrtE"


def _states() -> pd.DataFrame:
    directory = paths.data_dir() / "carotenoid"
    states = pd.read_csv(directory / "elizondo2025_steady_states.tsv", sep="\t")
    mrna = pd.read_csv(directory / "elizondo2025_relative_mrna.tsv", sep="\t")
    entry = mrna[mrna.gene == ENTRY_GENE][["condition", "rel_expression"]]
    frame = states.merge(entry, on="condition")
    frame["q_total"] = frame.q_lycopene + frame.q_betacarotene
    # Two dilution rates per strain; the tag pairs a strain's states with another's when
    # the expression profiles are permuted below.
    frame["rate_tag"] = (frame.mu_per_h > 0.2).map({True: "high", False: "low"})
    return frame


def candidate_targets(frame: pd.DataFrame) -> dict[str, pd.Series]:
    """Every quantity a user might ask for, built from the same six steady states."""
    return {
        "q  (specific rate, mmol/gDCW/h)": frame.q_total,
        "yield on glucose (mol/mol)": frame.q_total / frame.q_glucose.abs(),
        "content (q/mu, mmol/gDCW)": frame.q_total / frame.mu_per_h,
        "volumetric rate (q x biomass, mmol/L/h)": frame.q_total * frame.biomass_g_per_l,
        "titre (content x biomass, mmol/L)":
            frame.q_total / frame.mu_per_h * frame.biomass_g_per_l,
    }


def score(frame: pd.DataFrame) -> pd.DataFrame:
    strains = sorted(frame.strain.unique())
    profile = {(r.strain, r.rate_tag): r.rel_expression for r in frame.itertuples()}
    rows = []
    for name, target in candidate_targets(frame).items():
        observed = fit_flux_law(list(frame.rel_expression), list(target),
                                list(frame.strain), ENTRY_GENE, "elizondo2025")
        null = []
        for order in itertools.permutations(strains):
            mapping = dict(zip(strains, order))
            shuffled = [profile[(mapping[r.strain], r.rate_tag)] for r in frame.itertuples()]
            null.append(fit_flux_law(shuffled, list(target), list(frame.strain),
                                     ENTRY_GENE, "permuted").loso_skill)
        beaten = sum(1 for value in null if value >= observed.loso_skill)
        rows.append({
            "target": name,
            "loso_skill": observed.loso_skill,
            "typical_fold_error": observed.typical_fold_error,
            "permutation_rank": f"{beaten}/{len(null)}",
            "permutation_p": beaten / len(null),
            "entry_gene": ENTRY_GENE,
            "validation_scope": "exploratory fixed-gene comparison; separate target refit; not nested selection",
        })
    return pd.DataFrame(rows).sort_values("loso_skill", ascending=False)


def main() -> int:
    frame = _states()
    table = score(frame)
    n_strains = frame.strain.nunique()
    print(f"\nRelative {ENTRY_GENE} expression, separate target refits.\n"
          f"{len(frame)} states in {n_strains} independent strain groups, leave-one-STRAIN-out.\n"
          "Exploratory fixed-gene comparison against training-only baselines, not nested selection.\n")
    print(table.drop(columns=["entry_gene", "validation_scope"])
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    best = table.iloc[0]
    worst = table.iloc[-1]
    print(f"\nHighest fixed-gene skill: {best.target.strip()} -- {best.loso_skill:+.3f}, "
          f"{best.typical_fold_error:.2f}x typical")
    print(f"Lowest fixed-gene skill: {worst.target.strip()} -- {worst.loso_skill:+.3f}, "
          f"{worst.typical_fold_error:.2f}x typical")
    print(f"\nThe {worst.typical_fold_error / best.typical_fold_error:.1f}x spread across "
          "these separate target refits is not the cost of deriving content or titre.")
    print("Using the SAME rate prediction and the SAME measured conversion factors, relative "
          "and log errors are identical: mu, uptake and biomass cancel exactly.")
    print("That identity does not eliminate uncertainty in a new or unmeasured conversion factor.")
    floor = 1.0 / math.factorial(n_strains)
    print(f"\nThe expression-profile permutation floor is {floor:.3g} at {n_strains} strains.")
    if table.permutation_p.eq(floor).all():
        print("Every row attains that floor; this does not validate selecting the gene or target.")
    if floor > 0.05:
        print("The fixed-gene permutation test cannot reach 0.05 at this strain count.")

    out = paths.outputs_dir() / "target_scoreboard.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
