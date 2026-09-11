"""Report fixed-configuration transfer and a separate candidate optimality-gap bound.

Per-family transfer and rotated-subspace null comparisons use three fixed latent states
and the recommended co-doses, not `run_transfer.py`'s nested width selection. Scores are
median stressor-held-out channel R2 relative to outer-training channel means. Oracle scores
are in-sample diagnostics. Null non-rejection establishes neither equivalence nor absence.

A dimension/design configuration is selected on named training families and evaluated on
different named families. Their mean-score difference is descriptive: the fixed split
changes family difficulty as well as selection exposure, so it is not selection optimism.

SPOTA fits its own candidate on sampled families, independently of the reference sets.
Its approximate bound targets that fixed candidate's expected configuration optimality gap
under the uniform synthetic-family sampler, not each future loss or the named-split
configuration's optimism. Full selected-pipeline optimism is pending independent evaluation;
this script does not transport fitted latent weights or repeat the original co-dose search.
These synthetic comparisons are not biological validation.

Usage: YSTWIN_OUTPUTS=/external/output/dir python scripts/run_cross_family.py [--scale]
"""
from __future__ import annotations

import dataclasses
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
from ystwin.analysis.nulls import compare_to_null, rotated_subspace, skill_score_if_defined
from ystwin.analysis.optimism import OptimismProblem, estimate_optimism
from ystwin.analysis.transfer import leave_one_stressor_out
from ystwin.generator.families import build_family, family_dataset, family_names, sample_family
from ystwin.generator.panel_experiment import MEASURED_ACTIVITY_CV, MEASURED_GROWTH_RATE_SE

OUT = paths.outputs_dir()

PANEL = ["DTT", "H2O2", "heat", "NaCl", "glucose_starvation", "BPS", "MG132",
         "MMS", "congo_red", "calcium_chloride", "rapamycin", "copper_sulfate"]
READERS = ["STRE-general", "UPRE-ER", "TRX2-oxidative", "HSE-heat",
           "STRE-osmotic", "PACE-proteasome", "FeRE-iron", "Xbox-dna"]
OBSERVED = (0, 1, 2, 3)
"""Panel, readers and revealed channels shared with `run_transfer.py`.

Matching these does not match the selection procedure: this script uses fixed widths,
whereas that script now selects widths inside each outer training fold.
"""

TRAINING = ("baseline", "edge_dropped", "edge_added", "biphasic")
HELD_OUT = tuple(name for name in family_names() if name not in TRAINING)
"""Fit the configuration on some families, evaluate on the rest -- the M-open split of
CIRCULARITY.md section 8. The held-out families are the ones whose perturbation the selected
configuration has never been allowed to see."""

K_GRID = (2, 3, 4)
DESIGNS = {"blocked": (), "recommended pairs": tuple(RECOMMENDED_DESIGN)}
"""The dimension/design grid searched separately for the split and SPOTA candidates.

The designs are two declared alternatives, not a repetition of the search that produced
RECOMMENDED_DESIGN. Latent loadings are refitted per held-out stressor, not carried across
families. The fixed-candidate bound does not evaluate the full selected pipeline.
"""

MODEL = f"{leave_one_stressor_out.__module__}.{leave_one_stressor_out.__qualname__}"
SCORE_TARGET = "median stressor-held-out channel R2 relative to outer-training channel means"
PIPELINE_OPTIMISM_STATUS = "pending; independent evaluation of the full selection pipeline required"


def compact_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def configuration(solution, doses, replicates: int) -> str:
    n_states, design = solution
    return compact_json({
        "n_states": n_states, "combinations": DESIGNS[design],
        "doses_ec50": doses, "replicates": replicates,
        "observed": OBSERVED, "reporters": READERS, "stressors": PANEL,
        "noise_cv": MEASURED_ACTIVITY_CV, "growth_rate_se_per_h": MEASURED_GROWTH_RATE_SE,
    })


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def dataset(family, design, seed: int, doses, replicates: int):
    return family_dataset(family, reporters=READERS, stressors=PANEL, doses=doses,
                          replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV,
                          growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=seed,
                          combinations=design)


def median_transfer(data, n_states: int) -> float:
    table = leave_one_stressor_out(data, n_states=n_states, observed=OBSERVED, only=PANEL)
    return float(np.nanmedian(table.r2))


_SCORES: dict = {}


def family_score(family, design: str, seed: int, n_states: int, doses, replicates) -> float:
    """One family's transfer score, memoised.

    SPOTA scores every candidate on a reference's configurations while fitting it, then scores
    the winner and the candidate on those same configurations again, so a quarter of the
    sweeps are repeats. The key is a complete identity -- `build_family` is deterministic in
    (name, draw) and the rest of the key is every argument that reaches the generator -- so
    the cache cannot return a number from a different question.
    """
    key = (family.name, family.draw, design, seed, n_states, tuple(doses), replicates)
    if key not in _SCORES:
        _SCORES[key] = median_transfer(
            dataset(family, DESIGNS[design], seed, doses, replicates), n_states)
    return _SCORES[key]


def main(scale: bool) -> None:
    doses = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0) if scale else (0.25, 0.5, 1.0, 2.0)
    replicates = 6 if scale else 3
    null_draws = 40 if scale else 20
    n_reference = 20 if scale else 8
    n_candidate = 8 if scale else 4
    fixed = (3, "recommended pairs")
    family_seed = 0
    data_seed = 0
    spota_seed = 0

    rule("what would count as success, and what as failure")
    print("  Stated before the numbers; these are outcomes of the synthetic null tests.")
    print("  SUCCESS on the null -- detected superiority includes a non-baseline family.")
    print("  FAILURE on the null -- only the baseline shows detected superiority; this")
    print("    does not establish inferiority or absence of an effect in other families.")
    print("  THIRD OUTCOME -- no family beats its own rotated-subspace null: no detected")
    print("    superiority, not equivalence or absence of transfer information.")
    print("  The named-family split and the separate SPOTA candidate answer different")
    print("    questions. These comparisons are not biological validation.")
    print()
    print(f"  budget: {'--scale' if scale else 'default'}, {len(doses)} doses, "
          f"{replicates} replicates, {null_draws} null draws (p floor "
          f"{1 / (null_draws + 1):.3f}), {n_reference} reference solutions")

    rule(f"transfer per family, {len(family_names())} families, and the null for each")
    print(f"  Fixed configuration: {fixed[0]} states, {len(DESIGNS[fixed[1]])} co-dosed pairs.")
    print("  This is not run_transfer.py's nested-width pipeline or its outer score.")
    print("  Each fold excludes the held-out stressor and co-doses containing it; the")
    print("  R2 baseline is the outer-training channel mean. Oracle R2 is an in-sample")
    print("  diagnostic, not additional held-out evidence. Stressor folds share training data.")
    rows = []
    for name in family_names():
        family = build_family(name, seed=family_seed)
        data = dataset(family, DESIGNS[fixed[1]], data_seed, doses, replicates)
        table = leave_one_stressor_out(data, n_states=fixed[0], observed=OBSERVED, only=PANEL)
        observed = float(np.nanmedian(table.r2))

        def rotate(candidate, rng):
            return dataclasses.replace(
                candidate, readings=rotated_subspace(candidate.readings, rng))

        null = compare_to_null(lambda d: median_transfer(d, fixed[0]), data, rotate,
                               n_draws=null_draws, greater_is_better=True, label=name,
                               seed=data_seed)
        rows.append({
            "family": name,
            "perturbs": family.perturbs,
            "transfer_r2": observed,
            "oracle_r2": float(np.nanmedian(table.oracle_r2)),
            "alignment": float(np.nanmedian(table.alignment)),
            "null_median": null.null_median,
            "p_value": null.p_value,
            "beats_null": null.beats_null,
            "skill_over_null": skill_score_if_defined(observed, null.null_median, perfect=1.0),
        })
    families = pd.DataFrame(rows)
    baseline_r2 = float(families.loc[families.family == "baseline", "transfer_r2"].iloc[0])
    families["delta_vs_baseline"] = families.transfer_r2 - baseline_r2
    print(families.drop(columns=["perturbs"]).round(4).to_string(index=False))
    print()
    for _, row in families.iterrows():
        print(f"  {row.family:14s} {row.perturbs}")

    rule("does a configuration chosen on some families hold on families it never saw")
    print("  The optimisation is the choice of latent dimension and co-dosing design.")
    print(f"  chosen on : {', '.join(TRAINING)}")
    print(f"  applied to: {', '.join(HELD_OUT)}")
    training = [build_family(name, seed=family_seed) for name in TRAINING]
    chosen = best_solution(training, data_seed, doses, replicates, initial=None)
    print(f"  chosen    : {chosen[0]} states, {chosen[1]} design")
    rows = []
    for name in family_names():
        family = build_family(name, seed=family_seed)
        score = objective(chosen, [family], data_seed, doses, replicates)
        rows.append({"family": name,
                     "split": "chosen on" if name in TRAINING else "held out",
                     "transfer_r2": score})
    held = pd.DataFrame(rows)
    baseline_chosen = float(held.loc[held.family == "baseline", "transfer_r2"].iloc[0])
    held["delta_vs_baseline"] = held.transfer_r2 - baseline_chosen
    print(held.round(4).to_string(index=False))
    train_mean = float(held.loc[held.split == "chosen on", "transfer_r2"].mean())
    held_mean = float(held.loc[held.split == "held out", "transfer_r2"].mean())
    print(f"  mean on the families it was chosen on : {train_mean:+.4f}")
    print(f"  mean on the families it never saw     : {held_mean:+.4f}")
    print(f"  descriptive chosen-on minus held-out : {train_mean - held_mean:+.4f}")
    print("  Units are named synthetic families, not independent wells. This fixed split")
    print("  changes family difficulty as well as selection exposure; it is not an")
    print("  estimate of selection optimism.")

    rule("separate fixed-candidate expected optimality gap (UCBOG)")
    print("  SPOTA, Muratore et al. CoRL 2018: fit a separate candidate on sampled")
    print("  families, then hold it fixed independently of fresh reference family sets.")
    print("  Fit references and compare on their own sets under common random numbers;")
    print("  bootstrap one gap per independent reference family set, not per well.")
    problem = OptimismProblem(
        sample=sample_family,
        fit=lambda domains, seed, initial: best_solution(
            domains, seed, doses, replicates, initial),
        score=lambda solution, domains, seed: objective(
            solution, domains, seed, doses, replicates),
    )
    bound = estimate_optimism(problem, n_candidate_domains=n_candidate,
                              n_reference_domains=2, n_reference=n_reference,
                              alpha=0.05, n_bootstrap=1000, seed=spota_seed)
    print(f"  candidate solution : {bound.candidate[0]} states, {bound.candidate[1]} design")
    print(f"  {bound.summary()}")
    print(f"  raw gap samples    : {bound.raw_gaps.tolist()}")
    print(f"  clipped gap samples: {bound.gaps.tolist()}")
    if n_reference < 20:
        print("  n_G is below the published 20 to keep the default runnable; --scale uses 20.")
    print("  This approximate bound targets the separate candidate's expected optimality")
    print("  gap under iid uniform draws from the synthetic-family registry, not each")
    print("  future loss. It does not bound the named-split configuration's optimism or")
    print("  transfer to an unknown biological distribution.")

    rule("verdict")
    for line in verdict(families, train_mean, held_mean, bound.bound, baseline_r2):
        print(f"  {line}")

    fixed_configuration = configuration(fixed, doses, replicates)
    families.assign(
        model=MODEL, configuration=fixed_configuration, target=SCORE_TARGET,
        evaluation_unit="one synthetic family dataset; stressor folds share training data",
        selection_scope="fixed dimension/co-doses; latent weights refit per outer training fold",
        oracle_role="in-sample diagnostic, not a held-out prediction",
    ).to_csv(OUT / "cross_family_transfer.csv", index=False)
    pd.DataFrame([{
        "ucbog": bound.bound, "mean_gap": bound.mean_gap, "alpha": bound.alpha,
        "n_bootstrap": bound.n_bootstrap, "n_reference": bound.n_reference,
        "n_candidate_domains": bound.n_candidate_domains,
        "n_reference_domains": bound.n_reference_domains,
        "clipped_fraction": bound.clipped_fraction,
        "candidate_states": bound.candidate[0], "candidate_design": bound.candidate[1],
        "baseline_transfer_r2": baseline_r2,
        "chosen_on_mean": train_mean, "held_out_mean": held_mean,
        "model": MODEL,
        "transfer_configuration": fixed_configuration,
        "selected_configuration": configuration(chosen, doses, replicates),
        "candidate_configuration": configuration(bound.candidate, doses, replicates),
        "score_target": SCORE_TARGET,
        "configuration_grid": compact_json({"n_states": K_GRID, "designs": DESIGNS}),
        "split_target": "mean family transfer on named training/held-out sets; descriptive difference",
        "split_evaluation_unit": "named synthetic family; fixed split, not iid family draws",
        "split_selection_scope": compact_json({
            "training_families": [family.name for family in training],
            "held_out_families": HELD_OUT, "family_seed": family_seed, "data_seed": data_seed,
        }),
        "selected_pipeline_optimism_status": PIPELINE_OPTIMISM_STATUS,
        "bound_target": bound.target,
        "bound_resampling_unit": "independent reference family set; one gap per reference solution",
        "bound_selection_scope": compact_json({
            "sampler": f"{problem.sample.__module__}.{problem.sample.__qualname__}",
            "family_registry": family_names(), "family_weights": "uniform",
            "seed": spota_seed, "candidate_fixed_for_references": True,
        }),
        "bound_method": bound.method,
        "bound_coverage_status": bound.coverage_warning,
    }]).to_csv(OUT / "cross_family_optimism.csv", index=False)
    print(f"\nwrote {OUT/'cross_family_transfer.csv'} and "
          f"{OUT/'cross_family_optimism.csv'}")


def objective(solution, families, seed: int, doses, replicates: int) -> float:
    """Empirical mean of per-family median transfer scores, not a worst-case score.

    With iid sampled families this targets the risk-neutral expected score under the
    sampler. A mean over fixed named families describes only that set. Neither average
    supplies a worst-case or biological robustness guarantee.
    """
    n_states, design = solution
    return float(np.mean([family_score(family, design, seed, n_states, doses, replicates)
                          for family in families]))


def best_solution(families, seed: int, doses, replicates: int, initial):
    """Exhaustive argmax over the dimension and design grid.

    `initial` enters only as the tie-break. SPOTA passes the candidate in so a reference
    starts from it and cannot report a negative gap by landing in a worse local optimum; a
    discrete exhaustive search has no local optima to land in, so the requirement is met
    trivially and the initialisation does nothing except settle exact ties toward the
    candidate. Worth saying rather than leaving the argument unused and unexplained.
    """
    scored = {(k, design): objective((k, design), families, seed, doses, replicates)
              for k in K_GRID for design in DESIGNS}
    best = max(scored.values())
    winners = [key for key, value in scored.items() if value == best]
    if initial in winners:
        return initial
    return winners[0]


def verdict(families: pd.DataFrame, train_mean: float, held_mean: float,
            ucbog: float, baseline_r2: float) -> list[str]:
    """Report null-test decisions without equating non-rejection, equality and absence."""
    lines = []
    winners = families.loc[families.beats_null, "family"].tolist()
    variants = families.loc[families.family != "baseline"]
    worst = variants.loc[variants.transfer_r2.idxmin()]
    if not winners:
        lines.append("THIRD OUTCOME: no family beats its own rotated-subspace null.")
        lines.append("No detected superiority of fitted axes in these comparisons;")
        lines.append("non-rejection does not establish equivalence or absence of transfer information.")
    elif winners == ["baseline"]:
        lines.append("FAILURE on the null: only the baseline shows detected superiority.")
        lines.append("This does not establish inferiority or absence of an effect elsewhere.")
    else:
        lines.append(f"SUCCESS on the null: {len(winners)} of {len(families)} families beat "
                     "their own")
        lines.append(f"rotated-subspace null ({', '.join(winners)}).")
    lines.append("These synthetic comparisons are not biological validation.")
    lines.append("")
    lines.append(f"Spread across families: baseline {baseline_r2:+.4f}, worst variant "
                 f"{worst.transfer_r2:+.4f} ({worst.family}),")
    lines.append(f"range {families.transfer_r2.max() - families.transfer_r2.min():.4f} on a "
                 f"baseline score of {baseline_r2:.4f}.")
    lines.append(f"Descriptive split difference, chosen-on minus held-out: {train_mean - held_mean:+.4f}.")
    lines.append("The named-family split changes difficulty as well as selection exposure.")
    lines.append(f"Selected-pipeline optimism: {PIPELINE_OPTIMISM_STATUS}.")
    lines.append(f"Separate fixed-candidate UCBOG {ucbog:+.4f}: an approximate bound on the")
    lines.append("expected configuration optimality gap under the uniform synthetic-family sampler,")
    lines.append("not each future loss or the named-split configuration's optimism.")
    if baseline_r2:
        lines.append(f"For scale only: {abs(ucbog) / abs(baseline_r2):.1f}x the fixed-configuration "
                     "baseline score; these are distinct targets.")
    if ucbog > abs(baseline_r2):
        lines.append("This upper bound alone does not establish that the true gap is large")
        lines.append("or that any fresh family will incur excess loss.")
    return lines


if __name__ == "__main__":
    main(scale="--scale" in sys.argv)
