"""Compare mRNA and protein assays under an explicitly assumed translation generator.

    python3 scripts/simulate_proteomics_gain.py

Strain offsets and growth-rate slopes are motivated by residuals of the named fixed-CrtE
empirical comparison, not measured translation efficiencies. A variance decomposition of
those residuals describes where error occurred; it does not establish an irreducible
accuracy ceiling or identify its cause. Protein/mRNA studies motivate alternatives (Cell
Systems 2017 S2405-4712(17)30088-1; eLife 65722; Nat Commun 2022 s41467-022-30513-2),
but do not establish the mechanism in the Elizondo strains.

The generator explicitly makes flux proportional to latent protein, with independent
product-assay noise. One estimator sees mRNA and the other a noisy protein measurement.
With a perfect protein assay the protein model is the generator's own law, so its advantage
is conditional on that assumption, not an independent structural validation. This is Tier 1
simulation evidence about estimators, not biological evidence for a causal expression law.

More strains can reduce estimation error in the shared coefficient even though an unseen
strain's independent offset and slope remain unpredictable from its mRNA alone. Neither
leave-one-strain-out validation nor the residual decomposition implies a flat finite-sample
accuracy curve. A lower permutation floor makes significance possible, not guaranteed.

The historical column ``protein_assay_cv_log`` is the standard deviation of Gaussian log
measurement error, not the exact raw-scale coefficient of variation. The simulated noise
grid is illustrative, not a measured assay specification. The first sampled noise level
with no median advantage is a grid crossing, not a precisely estimated breakeven or an
economic decision. Numerical summaries are computed from the current seeded run rather
than repeated here from a committed artifact.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.pathway.flux import fit_flux_law

#: Motivated by fixed-CrtE leave-one-strain-out residuals, not measured translation. The
#: slope of residual on log(mu) is -0.205 / +0.452 / -0.356 across the three Elizondo
#: strains, and the strain offsets spread by 0.116 in log space.
SIGMA_SLOPE = 0.351
SIGMA_OFFSET = 0.116
#: Product-assay noise. The authors' own 95% interval is about 1.11x at low dilution rate
#: and 1.23x at high, so 0.10 in log space is the optimistic end of what they achieved.
SIGMA_PRODUCT = 0.10


def one_panel(n_strains: int, rates: tuple[float, ...], sigma_protein: float,
              rng: np.random.Generator) -> tuple[float, float]:
    """Simulate one panel and return (mRNA-model error, protein-model error), as folds.

    Each strain has an mRNA level and an assumed translation efficiency with its own
    growth-rate slope. Flux follows latent protein. The protein estimator therefore
    matches the generator when its assay is noiseless; the mRNA estimator omits that term.
    """
    strains = [f"s{i}" for i in range(n_strains)]
    slope = {s: rng.normal(0.0, SIGMA_SLOPE) for s in strains}
    offset = {s: rng.normal(0.0, SIGMA_OFFSET) for s in strains}
    mrna = {s: float(np.exp(rng.uniform(np.log(0.25), np.log(1.0)))) for s in strains}
    centre = float(np.mean(np.log(rates)))

    rows = []
    for strain in strains:
        for mu in rates:
            log_te = slope[strain] * (np.log(mu) - centre) + offset[strain]
            protein = mrna[strain] * np.exp(log_te)
            flux = protein * np.exp(rng.normal(0.0, SIGMA_PRODUCT))
            rows.append({
                "strain": strain, "mrna": mrna[strain], "flux": flux,
                # What a targeted MS assay would report, with its own error.
                "protein": protein * np.exp(rng.normal(0.0, sigma_protein)),
            })
    frame = pd.DataFrame(rows)
    out = []
    for column in ("mrna", "protein"):
        fit = fit_flux_law(list(frame[column]), list(frame.flux), list(frame.strain),
                           "entry", "simulated")
        out.append(fit.typical_fold_error)
    return out[0], out[1]


def sweep(n_strains: int, rates, trials: int, rng) -> pd.DataFrame:
    rows = []
    for sigma in (0.00, 0.05, 0.10, 0.15, 0.20, 0.30):
        pairs = [one_panel(n_strains, rates, sigma, rng) for _ in range(trials)]
        mrna = float(np.median([p[0] for p in pairs]))
        protein = float(np.median([p[1] for p in pairs]))
        rows.append({
            "protein_assay_cv_log": sigma,
            "mrna_model_fold": mrna,
            "protein_model_fold": protein,
            "error_removed_pct": 100.0 * (1.0 - (protein - 1.0) / max(mrna - 1.0, 1e-9)),
        })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strains", type=int, default=3)
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    rates = (0.101, 0.254)
    table = sweep(args.strains, rates, args.trials, rng)

    print("\nassumed generator: flux follows PROTEIN; translation efficiency has a per-strain slope")
    print(f"on log(mu) (sigma {SIGMA_SLOPE}) and a per-strain offset (sigma {SIGMA_OFFSET}),")
    print("motivated by fixed-CrtE residuals, not measured translation efficiencies.")
    print(f"Assumed product-assay log SD: {SIGMA_PRODUCT}; protein_assay_cv_log also denotes log SD.")
    print(f"{args.strains} strains at D = {rates}, {args.trials} panels per row.\n")
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    ideal = table.iloc[0]
    illustrative = table[table.protein_assay_cv_log == 0.15].iloc[0]
    print(f"\n  under this generator, a perfect protein assay removes {ideal.error_removed_pct:.0f}% "
          "of excess fold error")
    print(f"  illustrative assay log SD 0.15 removes {illustrative.error_removed_pct:.0f}%")
    crossing = table[table.protein_model_fold >= table.mrna_model_fold]
    if not crossing.empty:
        print(f"  first simulated log SD with no median advantage: "
              f"{crossing.iloc[0].protein_assay_cv_log:.2f}; a grid point, not a fitted breakeven")
    else:
        print("  a median advantage remains at each simulated log SD; no crossing was observed")

    print("\n  strain count, against the SAME mRNA model:")
    for count in (3, 6, 12):
        folds = [one_panel(count, rates, 0.10, rng)[0] for _ in range(200)]
        print(f"    {count:2d} strains   typical {float(np.median(folds)):.4f}x")
    print("    More strains can reduce estimation error in the shared coefficient; the")
    print("    new strain's independent offset and slope remain unobserved. This does not")
    print("    imply a flat accuracy curve. A lower permutation floor permits significance,")
    print("    but does not guarantee it. A residual decomposition alone establishes neither.")

    print("\n  This is Tier 1 simulation evidence about two ESTIMATORS under an assumed mechanism.")
    print("  With a noiseless protein assay that estimator is the generator's own law.")
    print("  The comparison does not establish a biological translation or cassette-dosage mechanism.")

    out = paths.outputs_dir() / "proteomics_gain.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
