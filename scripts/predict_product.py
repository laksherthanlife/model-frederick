"""Nested strain-group validation of expression-to-product model selection.

    python3 scripts/predict_product.py --output-dir /tmp/product-validation
    python3 scripts/predict_product.py --mode fixed_crte_comparison --draws 500
    python3 scripts/predict_product.py --legacy-demo --audit

Default validation selects the entry gene AND entry/branch model family in inner
strain folds, then refits every scale and branch parameter on the outer training
strains. Both beta-carotene and lycopene are predicted without held-out outcomes.
Three strains generate six condition predictions; those are not six independent units.

The constant-content, constant-rate, entry-plus-training-partition and saturating
families are compared on the same outer folds. The historical fixed CrtE model is an
explicit data-selected comparison, not a prespecified validation. Conditional bands
resample paired training strains and released expression/observation intervals, refitting
entry and kinetics together. Nominal and actual channel/joint coverage are distinct.
Optional public-chain demonstrations use the actual legacy diagnostic mode, never
silently converting a refused prior into a supported prediction.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.bridge.latent_bridge import LatentState
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.flux import (
    CAROTENOID_GENES,
    INTERVAL_ASSUMPTIONS,
    carotenoid_measurements,
    fit_saturating_branch,
    score_product_validation,
    summarize_product_validation,
)
from ystwin.pathway.solve import NodeKinetics
from ystwin.pathway.spec import PathwaySpec, load_pathway
from ystwin.predict import Environment, Genotype, SetpointUnreachable, predict_product

PATHWAY = "beta_carotene"
ENTRY_GENE = "CrtE"

# Cofactor and precursor levels the solve does not pin. Order-of-magnitude placeholders,
# and nothing printed here rests on their values: the carotenoid steps have no formation
# energy in the vendored tables either way, so the gate reports the coverage gap rather
# than a verdict. They are supplied so the gate runs and says so, instead of sitting out.
CAROTENOID_BACKGROUND_M = {
    "s_0189": 1e-5,   # GGPP, the precursor
    "s_0633": 1e-4,   # diphosphate
    "s_0687": 1e-4,   # FAD
    "s_0689": 1e-4,   # FADH2
    "s_0794": 1e-7,   # proton, pH 7
}


def _thermodynamics():
    """The vendored formation energies, or ``None`` where they are not on this machine.

    Returned rather than raised: the gate is one of five layers and the other four have
    answers without it, so a missing table degrades the report by one line instead of
    stopping the script. `predict_product` records it as ``not-run`` and says what to pass.
    """
    from ystwin.bridge.thermodynamic import ThermodynamicData

    try:
        return ThermodynamicData.load()
    except (FileNotFoundError, OSError) as exc:
        print(f"  (thermodynamic tables unavailable: {exc}; that layer will read not-run)")
        return None


def rule(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def measured() -> pd.DataFrame:
    return carotenoid_measurements()


def fit_branch(spec: PathwaySpec, train: pd.DataFrame) -> dict[str, NodeKinetics]:
    """Refit both training pools through the shared solver, including the training window."""
    # The window these parameters were actually fitted over, read off the training
    # subset rather than copied from the shipped calibration -- a leave-one-strain-out
    # fit sees a different set of rates from the full one, and hard-coding the full
    # window would claim support this fit does not have.
    return fit_saturating_branch(spec, train)


def score_forward(
    spec: PathwaySpec, states: pd.DataFrame, *, mode: str = "nested",
    candidate_genes: tuple[str, ...] = CAROTENOID_GENES,
    entry_laws: tuple[str, ...] = ("rate", "content"),
    branch_families: tuple[str, ...] = ("partition", "saturating"),
    draws: int = 500, seed: int = 0, nominal_coverage: float = 0.95,
) -> pd.DataFrame:
    """Select genes and model families inside outer strain folds; never rank globally.

    ``fixed_crte_comparison`` reproduces the historically selected model as a named
    comparison. It must not be interpreted as a prespecified score. Interval assumptions
    and full failed-condition/candidate denominators are returned with the predictions.
    """
    # The held-out bands refit entry and kinetics on the same paired-strain draws.
    # Source expression and product intervals supply the second, observation stage.
    # This includes kinetic, entry-expression and observation uncertainty conditionally,
    # not just a monotone remapping of the flux law's held-out spread. Nominal coverage
    # and actual channel/joint coverage are reported separately on three strain groups.
    return score_product_validation(
        spec, states, genes=candidate_genes, entry_laws=entry_laws, branches=branch_families,
        mode=mode, draws=draws, seed=seed, nominal_coverage=nominal_coverage)


def legacy_demo(spec: PathwaySpec, audit: bool = False) -> None:
    """Explicit public-API legacy diagnostics, separate from the validated model comparison."""
    rule("Legacy diagnostic comparison, NOT a supported product prediction")
    genotype = Genotype(1.0, "b-car4")
    thermo = _thermodynamics()
    diagnostic = predict_product(spec, genotype, Environment(growth_rate_setpoint_per_h=0.18),
                                 BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
                                 thermo=thermo, thermo_background_m=CAROTENOID_BACKGROUND_M,
                                 latent=LatentState(0.02, "UPRE2", 4.0), mode="legacy")
    got = diagnostic.calculation
    print(f"  in:  {ENTRY_GENE} relative expression {genotype.entry_expression:g}, "
          f"chemostat D = 0.18 /h")
    print(f"  out: {diagnostic.summary()}")
    print("\n  every layer the goal names, and what each one did:")
    for layer in got.layers:
        print(f"       {layer.summary()}")
    for name, content in got.intermediates().items():
        print(f"       {name:<14} {content:.4g} mmol/gDCW  (predicted, not supplied)")
    for note in got.notes:
        print(f"       note: {note}")

    rule("A stressor in a chemostat does not change mu -- it decides whether one exists")
    for dose in (0.2, 1.0, 5.0):
        try:
            s = predict_product(spec, genotype,
                                Environment(growth_rate_setpoint_per_h=0.15, stressor="DTT", dose=dose),
                                BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, mode="legacy")
            print(f"  DTT {dose:>4g} mM  {s.summary()}")
        except SetpointUnreachable as exc:
            print(f"  DTT {dose:>4g} mM  REFUSED: {str(exc)[:74]}")

    if audit:
        rule("What the GEM says about a flux it did not supply")
        import warnings

        from ystwin.fba.audit import audit_predicted_flux
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID, add_beta_carotene_pathway
        from ystwin.fba.solver import load_model

        gem = paths.yeast_gem()
        if gem is None or not gem.is_file():
            print("  Yeast9 not present; set YSTWIN_YEAST_GEM. Nothing audited.")
        else:
            warnings.simplefilter("ignore")
            model, settings = load_model(gem)
            model = add_beta_carotene_pathway(model)
            print(f"  solver {settings}")
            for uptake in (10.0, 1.5):
                verdict = audit_predicted_flux(
                    model, PRODUCT_DEMAND_ID, got.flux.flux_mmol_per_gdcw_h,
                    glucose_uptake=uptake,
                    precursor_metabolite=spec.precursor_metabolite)
                print(f"  {verdict.summary()}")
            print("\n  The feasible minimum is 0 at every uptake, which is why this layer")
            print("  audits rather than supplies: an upper bound cannot make a flux happen.")



def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("nested", "fixed_crte_comparison"), default="nested")
    parser.add_argument("--genes", nargs="+", default=list(CAROTENOID_GENES))
    parser.add_argument("--entry-laws", nargs="+", choices=("rate", "content"), default=["rate", "content"])
    parser.add_argument("--branches", nargs="+", choices=("partition", "saturating"),
                        default=["partition", "saturating"])
    parser.add_argument("--draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--nominal-coverage", type=float, default=0.95)
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument("--legacy-demo", action="store_true",
                        help="also run explicitly unsupported legacy public-chain diagnostics")
    parser.add_argument("--audit", action="store_true", help="audit the optional legacy diagnostic against Yeast9")
    args = parser.parse_args()
    spec, states = load_pathway(PATHWAY), measured()
    if args.legacy_demo or args.audit:
        legacy_demo(spec, args.audit)

    rule("Nested gene AND model selection, outer leave-one-strain-out" if args.mode == "nested"
         else "Fixed CrtE comparison: selected on these data, not prespecified")
    scored = score_forward(spec, states, mode=args.mode, candidate_genes=tuple(args.genes),
                           entry_laws=tuple(args.entry_laws), branch_families=tuple(args.branches),
                           draws=args.draws, seed=args.seed, nominal_coverage=args.nominal_coverage)
    summary = summarize_product_validation(scored)
    print(f"  {states.strain.nunique()} independent strains; {len(scored)} condition predictions; "
          "beta-carotene AND lycopene. Failures stay in the denominator.")
    print(scored[["state", "strain", "selected_model", "status", "product_measured",
                  "product_predicted", "lycopene_rate_measured", "lycopene_rate_predicted"]]
          .to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    columns = ["model", "n_independent_strains", "n_condition_predictions", "n_failed_predictions",
               "product_rmse_log", "lycopene_rmse_log", "joint_rmse_log",
               "rmse_skill_vs_constant_rate", "rmse_skill_vs_constant_content"]
    print("\n  Both channels are rates, mmol/gDCW/h. Content = rate / measured mu exactly.")
    print(summary[columns].to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("  Fixed CrtE was selected on this dataset; it is not a prespecified validation.")
    print("\n  Conditional uncertainty: " + INTERVAL_ASSUMPTIONS)
    own = summary[summary.model == args.mode].iloc[0]
    if args.draws:
        print(f"  nominal marginal coverage {args.nominal_coverage:.1%}; "
              f"beta-carotene {own.product_covered_conditions}/{len(scored)}, "
              f"lycopene {own.lycopene_covered_conditions}/{len(scored)}")
        print(f"  nominal simultaneous strain coverage {args.nominal_coverage:.1%}; "
              f"both channels {own.joint_covered_conditions}/{len(scored)} conditions, "
              f"all conditions {own.joint_covered_strains}/{states.strain.nunique()} strains")
        print(f"  intervals available: {own.n_intervals_available}/{len(scored)}; "
              "these counts do not establish nominal coverage.")
    else:
        print("  intervals not requested (--draws 0); empirical coverage is unavailable.")

    directory = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    directory.mkdir(parents=True, exist_ok=True)
    scored.to_csv(directory / "product_prediction.csv", index=False)
    summary.to_csv(directory / "product_prediction_scores.csv", index=False)
    pd.DataFrame(scored.attrs["selection_scores"]).to_csv(
        directory / "product_prediction_selection.csv", index=False)
    pd.DataFrame(scored.attrs["comparison_predictions"]).to_csv(
        directory / "product_prediction_comparisons.csv", index=False)
    print(f"\n  wrote validation tables to {directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
