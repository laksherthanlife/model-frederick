"""Reproduce the fixed flux calibration and separately validate nested gene/model selection.

    python3 scripts/fit_pathway_flux.py
    python3 scripts/fit_pathway_flux.py --gene CrtYB
    python3 scripts/fit_pathway_flux.py --output-dir /tmp/pathway-flux-validation

The fourteen-gene ranking is exploratory: CrtE was selected using these same outcomes.
It reproduces BETA_CAROTENE_FLUX, including its training-only baseline, but it does not
validate choosing the winner. The second table selects both the gene and rate/content
entry plus partition/saturating branch family inside the outer training folds. Both
carotenoid channels are scored on three strain groups and six condition predictions.
Unsupported genes and failed conditions remain in their respective denominators.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.pathway.flux import (
    CAROTENOID_GENES,
    carotenoid_measurements,
    fit_flux_law,
    score_product_validation,
    summarize_product_validation,
)
from ystwin.pathway.spec import load_pathway

STATES = "elizondo2025_steady_states.tsv"
MRNA = "elizondo2025_relative_mrna.tsv"
SOURCE = "PMID 40891387 (doi 10.1021/acssynbio.5c00256), Elizondo 2025"

# Genes native to the host, labelled explicitly in the exploratory comparison. They also
# remain in the declared nested candidate set; no outer outcome excludes a gene. Relative
# expression shares a normalisation, so this comparison does not establish causal dosage.
NATIVE = ("ERG10", "ERG12", "ERG13", "ERG20", "ERG8", "ERG9",
          "HMG1", "HMG2", "BTS1", "IDI1", "MVD1")


def load() -> pd.DataFrame:
    directory = paths.data_dir() / "carotenoid"
    for name in (STATES, MRNA):
        if not (directory / name).exists():
            raise SystemExit(f"{directory / name} is tracked and missing; bad checkout")
    return carotenoid_measurements(directory)


def rank_genes(frame: pd.DataFrame, genes=CAROTENOID_GENES) -> pd.DataFrame:
    """Exploratory fixed-gene LOSO scores, retaining every requested gene on failure."""
    rows = []
    for gene in genes:
        kind = ("native" if gene in NATIVE else "heterologous"
                if gene in {"CrtE", "CrtI", "CrtYB"} else "undeclared")
        row = {"gene": gene, "kind": kind,
               "n_independent_strains": frame.strain.nunique(),
               "n_condition_predictions": len(frame)}
        try:
            if gene not in frame:
                raise ValueError(f"missing expression for {gene}")
            fit = fit_flux_law(list(frame[gene]), list(frame.flux), list(frame.strain), gene, SOURCE)
            row.update(alpha=fit.alpha, loso_rmse_log=fit.loso_rmse_log,
                       loso_skill=fit.loso_skill, typical_fold_error=fit.typical_fold_error,
                       status="ok", failure_reason="")
        except ValueError as exc:
            row.update(alpha=float("nan"), loso_rmse_log=float("inf"),
                       loso_skill=float("-inf"), typical_fold_error=float("inf"),
                       status="unsupported", failure_reason=str(exc))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("loso_skill", ascending=False, kind="stable")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", default=None, help="restrict the declared gene candidates")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="directory for the validation tables; default paths.outputs_dir()")
    args = parser.parse_args()

    frame = load()
    genes = (args.gene,) if args.gene else CAROTENOID_GENES
    ranking = rank_genes(frame, genes)
    print("Exploratory fixed-gene flux LOSO; ranking and evaluating its winner is not nested validation")
    print(ranking.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("Skill = 1 - log RMSE / training-only geometric-mean log RMSE.")
    natives = ranking[ranking.kind == "native"]
    if not natives.empty:
        best_native = natives.iloc[0]
        print(f"Best native {best_native.gene}: skill {best_native.loso_skill:+.6g}.")
    print("This exploratory separation does not establish causal cassette dosage.")

    scored = score_product_validation(load_pathway("beta_carotene"), frame, genes=genes)
    summary = summarize_product_validation(scored)
    print(f"\nNested gene AND model selection: {frame.strain.nunique()} independent strains, "
          f"{len(frame)} condition predictions, both beta-carotene and lycopene")
    print(scored[["state", "strain", "selected_model", "inner_rmse_log", "status"]].to_string(index=False))
    columns = ["model", "n_independent_strains", "n_condition_predictions", "n_failed_predictions",
               "product_rmse_log", "lycopene_rmse_log", "joint_rmse_log",
               "rmse_skill_vs_constant_rate", "rmse_skill_vs_constant_content"]
    print(summary[columns].to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("Fixed CrtE is a selected-on-this-dataset comparison, not a prespecified model.")
    # Defaults to the outputs directory like `predict_product.py` and
    # `fit_carotenoid_kinetics.py`, so the committed tables have a no-argv generator.
    directory = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    directory.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(directory / "pathway_flux_ranking.csv", index=False)
    scored.to_csv(directory / "pathway_flux_nested.csv", index=False)
    summary.to_csv(directory / "pathway_flux_scores.csv", index=False)
    # The selection scores are the same object predict_product writes; one owner, one file.
    # See outputs/product_prediction_selection.csv.
    pd.DataFrame(scored.attrs["comparison_predictions"]).to_csv(
        directory / "pathway_flux_comparisons.csv", index=False)
    print(f"wrote validation tables to {directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
