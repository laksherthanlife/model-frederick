"""Random-design demonstration and descriptive PHB candidate scores, not significance bars.

    python3 scripts/flux_null_baseline.py

``flux_null_baseline.csv`` retains five model-size rows. Each draw uses new Gaussian
covariates AND shuffled log flux, so its percentiles describe that artificial design
ensemble, not the null distribution of any real candidate. The old universal-threshold
interpretation, and its hard-coded numerical claims, were incorrect.

``flux_candidate_scores.csv`` reports the four fixed candidates separately. Both the model
and the constant comparator are fitted only on each training fold. State-LOO skill is
1 - RMSE(model) / RMSE(training-mean comparator), pooled across all held-out states; it is
not R-squared or evidence of independent biological replication. Undefined scores and
rank-deficient training fits are retained as failures, not pseudoinverse forecasts.

The tracked Kocharin 2013 table contains eleven condition summaries: four glucose, three
ethanol and four mixed-feed states, with four dilution rates overall. It has no cultivation
or biological-replicate identifiers, or randomization scheme. ``data/phb/SOURCE.md`` does
not supply that provenance for 2013. Neither eleven rows nor three treatment labels
establishes independent experimental units or exchangeability. Real-data inference is
therefore PENDING, with no candidate p-value or null percentile issued.

``scripts/environment_channel.py`` instead scores log CONTENT using an ethanol-carbon
fraction descriptor, state/feed holdouts and descriptor relabellings across three feeds.
That is a different target, candidate family and permutation operation; its six
relabellings do not justify response shuffles for these log-FLUX candidates.

``candidate_null_distribution`` holds a supplied design and folds fixed and refits every
explicit response permutation, including its training-only comparator. It is computational
machinery, not a certificate of exchangeability: the caller must justify the permutation
scheme for the null being tested. No default scheme is inferred from this dataset. If
candidate selection is part of a future claim, that selection must also be validated;
a fixed-candidate distribution alone is not a selection-adjusted test.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.analysis.validation import ValidationFold, group_folds

SEED = 0
SHUFFLES = 2000
VALIDATION = "leave_one_state_out_descriptive"
BASELINE = "fold_training_mean_log_flux"
INFERENCE_REASON = (
    "The condition-summary TSV lacks cultivation/biological-replicate identifiers and "
    "a justified randomization or exchangeability scheme; carbon-source labels do not "
    "establish independent units. No candidate null or p-value is issued."
)


def load() -> pd.DataFrame:
    path = paths.data_dir() / "phb" / "kocharin2013_chemostat_states.tsv"
    if not path.exists():
        raise SystemExit(f"{path} is tracked and missing; bad checkout")
    return pd.read_csv(path, sep="\t")


def _log_positive(values, name: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError(f"{name} must be finite and positive before log scoring")
    return np.log(values)


def _skill(observed: np.ndarray, design: np.ndarray, *,
           folds: tuple[ValidationFold, ...] | None = None) -> float:
    """Pooled RMSE skill against each fold's TRAINING mean, state-LOO by default.

    Folds are operational partitions, not assertions of biological independence. Every
    state must be tested exactly once. A missing/rank-deficient fit or zero baseline
    error refuses the whole score rather than dropping an inconvenient held-out state.
    """
    observed, design = np.asarray(observed, dtype=float), np.asarray(design, dtype=float)
    if (observed.ndim != 1 or len(observed) < 2 or design.ndim != 2
            or design.shape[0] != len(observed) or design.shape[1] < 1):
        raise ValueError("need at least two observations and a matching nonempty design matrix")
    if not np.isfinite(observed).all() or not np.isfinite(design).all():
        raise ValueError("observations and design must be finite")
    n = len(observed)
    folds = group_folds(np.arange(n)) if folds is None else tuple(folds)
    coverage = np.zeros(n, dtype=int)
    for fold in folds:
        if any(np.any(indices >= n) for indices in (fold.train, fold.test, fold.excluded)):
            raise ValueError("validation fold index out of range")
        coverage[fold.test] += 1
    if not np.all(coverage == 1):
        raise ValueError("validation folds must test every observation exactly once")
    residuals = np.empty(n)
    baseline_residuals = np.empty(n)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        for fold in folds:
            coefficients, _, rank, _ = np.linalg.lstsq(
                design[fold.train], observed[fold.train], rcond=None)
            if rank < design.shape[1]:
                raise ValueError(
                    f"rank-deficient training design in fold {fold.held_out!r}: "
                    f"rank {rank} < {design.shape[1]}")
            residuals[fold.test] = design[fold.test] @ coefficients - observed[fold.test]
            baseline_residuals[fold.test] = observed[fold.train].mean() - observed[fold.test]
        model = float(np.sqrt(np.mean(residuals ** 2)))
        baseline = float(np.sqrt(np.mean(baseline_residuals ** 2)))
    if not np.isfinite(model) or not np.isfinite(baseline):
        raise ValueError("model and baseline RMSE must be finite")
    if baseline == 0.0 or np.all(observed == observed[0]):
        raise ValueError("training-mean baseline RMSE is zero; skill is undefined")
    skill = 1.0 - model / baseline
    if not np.isfinite(skill):
        raise ValueError("RMSE ratio is nonfinite; skill is undefined")
    return skill


def _try_skill(observed, design, *, folds=None) -> tuple[float, str]:
    """Retain failed draws/candidates as undefined, with a reason for candidate reports."""
    try:
        return _skill(observed, design, folds=folds), ""
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        return float("nan"), str(exc)


def candidate_null_distribution(observed: np.ndarray, design: np.ndarray,
                                permutations: np.ndarray | None = None, *,
                                folds: tuple[ValidationFold, ...] | None = None) -> np.ndarray:
    """Refit a FIXED candidate on explicitly supplied response-index permutations.

    No permutation scheme is inferred: raw-response exchangeability must be justified
    for the intended null, including any nuisance structure. The caller also owns whether
    permutations are sampled or exhaustive; this helper computes no inferential p-value.
    Each permutation keeps the design/folds fixed and refits the training comparator.
    Failed fits occupy their original positions as NaN, never a shortened denominator.
    """
    if permutations is None:
        raise ValueError(
            "explicit response permutations are required; exchangeability is not established")
    observed = np.asarray(observed, dtype=float)
    orders = np.asarray(permutations)
    if (observed.ndim != 1 or orders.ndim != 2 or orders.shape[0] == 0
            or orders.shape[1] != observed.size or orders.dtype.kind not in "iu"
            or not np.all(np.sort(orders, axis=1) == np.arange(observed.size))):
        raise ValueError("each permutation must contain every observation index exactly once")
    folds = group_folds(np.arange(len(observed))) if folds is None else tuple(folds)
    return np.array([_try_skill(observed[order], design, folds=folds)[0] for order in orders])


def null_distribution(states: pd.DataFrame, n_parameters: int,
                      rng: np.random.Generator) -> np.ndarray:
    """Legacy helper: RANDOM-DESIGN demonstration, not a fixed-candidate null.

    Gaussian covariates are redrawn on every response shuffle. This artificial ensemble
    does not establish exchangeability of the biological observations. Failed draws are
    NaN and must not be discarded to manufacture a finite percentile.
    """
    if (not isinstance(n_parameters, (int, np.integer)) or isinstance(n_parameters, bool)
            or n_parameters < 1):
        raise ValueError("n_parameters must be a positive integer")
    observed = _log_positive(states.q_phb_mmol_per_gdcw_h, "PHB flux")
    n = len(observed)
    folds = group_folds(np.arange(n))
    # An intercept plus (k-1) random covariates: the same degrees of freedom a real model of
    # this size would spend, with none of the meaning.
    skills = np.empty(SHUFFLES)
    for i in range(SHUFFLES):
        design = np.column_stack(
            [np.ones(n)] + [rng.normal(size=n) for _ in range(n_parameters - 1)])
        skills[i] = _try_skill(rng.permutation(observed), design, folds=folds)[0]
    return skills


def main() -> int:
    states = load()
    rng = np.random.default_rng(SEED)
    observed = _log_positive(states.q_phb_mmol_per_gdcw_h, "PHB flux")

    print("=" * 84)
    print(f"Kocharin 2013: {len(states)} condition summaries, "
          f"{states.carbon_source.nunique()} carbon feeds, "
          f"{states.dilution_rate_per_h.nunique()} dilution rates; one reported genotype/study.")
    print(f"log(q_PHB) spans {observed.max() - observed.min():.2f}; "
          "biological independence is not established by the row count.")
    print("=" * 84)
    print(f"\nRandom-design demonstration: {SHUFFLES} draws per size, seed {SEED}.")
    print("State-LOO RMSE skill versus the fold's training mean of log flux.")
    print("Redraws Gaussian covariates AND shuffles responses; not a candidate-specific threshold.\n")
    print(f"  {'parameters':>11} {'median':>9} {'90th pct':>10} {'95th pct':>10} "
          f"{'99th pct':>10} {'max':>8} {'failed':>8}")
    demonstrations = []
    for k in (1, 2, 3, 4, 5):
        skills = null_distribution(states, k, rng)
        n_failed = int((~np.isfinite(skills)).sum())
        quantiles = (np.percentile(skills, [50, 90, 95, 99, 100]) if n_failed == 0
                     else np.full(5, np.nan))
        median, p90, p95, p99, maximum = quantiles
        demonstrations.append({
            "parameters": k,
            "null_skill_p95": p95,
            "analysis_kind": "random_design_demonstration_not_inference",
            "validation": VALIDATION,
            "baseline": BASELINE,
            "n_states": len(states),
            "shuffles": SHUFFLES,
            "seed": SEED,
            "n_failed": n_failed,
            "score_status": "failed" if n_failed else "ok",
        })
        print(f"  {k:>11} {median:9.3f} {p90:10.3f} {p95:10.3f} "
              f"{p99:10.3f} {maximum:8.3f} {n_failed:8d}")
    print("\nPercentiles are withheld if any draw fails; no finite-only filtering.")

    print("\n" + "=" * 84)
    print("Fixed candidates: descriptive state-LOO scores, not independent-unit validation")
    print("=" * 84)
    n = len(states)
    feeds = pd.get_dummies(states.carbon_source, drop_first=True).to_numpy(dtype=float)
    log_mu = _log_positive(states.dilution_rate_per_h, "dilution rate")
    candidates = {
        "constant only": np.ones((n, 1)),
        "log(mu)": np.column_stack([np.ones(n), log_mu]),
        "carbon source": np.column_stack([np.ones(n), feeds]),
        "carbon source + log(mu)": np.column_stack([np.ones(n), feeds, log_mu]),
    }
    rows = []
    print(f"\n  {'model':<26} {'k':>3} {'LOO skill':>10}  score status")
    for name, design in candidates.items():
        skill, failure = _try_skill(observed, design)
        status = "failed" if failure else "ok"
        rows.append({
            "model": name,
            "parameters": design.shape[1],
            "loo_skill": skill,
            "score_status": status,
            "failure_reason": failure,
            "inference_status": "pending_exchangeability",
            "inference_reason": INFERENCE_REASON,
            "validation": VALIDATION,
            "baseline": BASELINE,
            "n_states": n,
            "n_carbon_sources": states.carbon_source.nunique(),
        })
        print(f"  {name:<26} {design.shape[1]:>3} {skill:10.3f}  {status} {failure}")
    print(f"\nInference PENDING: {INFERENCE_REASON}")

    out_dir = paths.outputs_dir()
    pd.DataFrame(demonstrations).to_csv(out_dir / "flux_null_baseline.csv", index=False)
    pd.DataFrame(rows).to_csv(out_dir / "flux_candidate_scores.csv", index=False)
    print(f"\n  wrote {out_dir / 'flux_null_baseline.csv'} and "
          f"{out_dir / 'flux_candidate_scores.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
