"""Separate per-dose fold inference from whole-sensor dose association.

``scripts/fold_multiplicity.py`` retains all 24 planned positive-dose cells, including
refusals. Per-dose ``p_exact_sign`` tests a median usable-plate fold of one, with exact
ties excluded. Holm and BH adjust those sign-test p-values, not fitted bootstrap tails.
The sign-test floor at four non-tied plates is 0.125; it bounds this test only.

The shared ``analysis/uncertainty.py`` estimator targets a geometric-mean fold. With
well identifiers it uses two-stage plate/well resampling; otherwise it uses cluster-only
resampling. The cluster-only empirical atom policy is not a universal resolution limit
and does not apply to two-stage draws. Finite Monte Carlo tails and estimator refusals
still matter. Bonferroni intervals are approximate simultaneous bounds, not a guarantee
of nominal coverage, and a refused interval does not establish no response.

A separate family of four sensor tests permutes dose labels within plates, using
plate-averaged Spearman correlation and a plus-one Monte Carlo p-value. The configured
cutoff is 1 mM, with 20000 permutations and seed 0. Under that configuration the reported
sensor summaries are:

    construct      plate-averaged rho   Holm-adjusted p
    UPRE1          <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=UPRE1;measurement=dose_response_permutation" -->+0.9465   <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=UPRE1;measurement=dose_response_permutation" -->0.000200
    UPRE2          <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=UPRE2;measurement=dose_response_permutation" -->+0.9165   <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=UPRE2;measurement=dose_response_permutation" -->0.000200
    NativeYap1     <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=NativeYap1;measurement=dose_response_permutation" -->+0.4692   <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=NativeYap1;measurement=dose_response_permutation" -->0.000500
    AlteredYap1    <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=AlteredYap1;measurement=dose_response_permutation" -->+0.4337   <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=AlteredYap1;measurement=dose_response_permutation" -->0.001000

AlteredYap1's unadjusted value is
<!-- audit:value table=outputs/fold_multiplicity.csv column=p_permutation row="construct=AlteredYap1;measurement=dose_response_permutation" -->0.001000 in ``p_permutation``. These selectors require the current artifact receipt;
cell agreement alone does not establish provenance or the test's scientific assumptions.

Permutation inference requires exchangeable dose assignments within plates under the null
and a prespecified cutoff. It concerns dose association, not a causal mechanism, a
particular fold size, or practical equivalence. Holm needs valid marginal p-values; BH
also needs independence or suitable positive dependence, not established by shared controls.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

__all__ = [
    "DoseResponse",
    "MultipleTests",
    "adjusted_interval",
    "benjamini_hochberg",
    "dose_response_permutation",
    "holm",
    "finest_resolvable_alpha",
    "sign_test_floor",
]


# --------------------------------------------------------------------------- corrections


def _family(p_values, family_size):
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or np.isinf(values).any() or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be probabilities in [0, 1], or NaN for unestimated tests")
    size = len(values) if family_size is None else family_size
    if not isinstance(size, (int, np.integer)) or isinstance(size, bool) or size < len(values):
        raise ValueError("family_size must be an integer at least as large as the submitted family")
    return values, int(size)


def holm(p_values: list[float], family_size: int | None = None) -> list[float]:
    """Holm-Bonferroni adjustment of valid marginal p-values for family-wise error control.

    ``family_size`` includes planned tests that were not estimable. Missing submitted
    p-values remain missing and do not reduce that denominator. This correction cannot
    turn a fitted empirical bootstrap tail into a valid confirmatory p-value.
    """
    values, size = _family(p_values, family_size)
    order = np.argsort(np.where(np.isnan(values), 1.0, values))
    adjusted = np.full(len(values), np.nan)
    running = 0.0
    for rank, index in enumerate(order):
        if np.isnan(values[index]):
            continue
        running = max(running, (size - rank) * values[index])
        adjusted[index] = min(running, 1.0)
    return [float(value) for value in adjusted]


def benjamini_hochberg(p_values: list[float], family_size: int | None = None) -> list[float]:
    """Benjamini-Hochberg adjustment, retaining the full planned family denominator.

    False-discovery-rate control requires valid marginal p-values and independence or
    suitable positive dependence. Shared-control dose tests do not automatically satisfy
    that condition. These columns accompany Holm, not replace it or certify its inputs.
    """
    values, size = _family(p_values, family_size)
    n = len(values)
    order = np.argsort(np.where(np.isnan(values), 1.0, values))[::-1]
    adjusted = np.full(n, np.nan)
    running = 1.0
    for rank, index in enumerate(order):
        if np.isnan(values[index]):
            continue
        running = min(running, (size / (n - rank)) * values[index])
        adjusted[index] = min(running, 1.0)
    return [float(value) for value in adjusted]


def sign_test_floor(n_clusters: int) -> float:
    """Smallest two-sided exact sign-test p for this many non-tied independent clusters.

    When all signs agree the floor is ``2 * (1/2)**n``: four usable, non-tied plates
    give 0.125; six can first reach 0.05. Ten can first reach a Bonferroni level of
    ``0.05 / 24``. Exclusions and ties reduce the usable count. This bounds the sign
    test, not parametric tests, bootstrap intervals or within-plate dose permutations.
    """
    if not isinstance(n_clusters, (int, np.integer)) or isinstance(n_clusters, bool) or n_clusters < 1:
        raise ValueError(f"need at least one cluster and an integer count, got {n_clusters}")
    return min(2.0 * 0.5 ** n_clusters, 1.0)


def finest_resolvable_alpha(n_clusters: int) -> float:
    """Conservative empirical atom policy for cluster-only percentile reporting.

    Resampling n whole plates n times with replacement assigns probability ``(1/n)**n``
    to an ordered draw; this policy uses twice that probability as a two-sided threshold.
    It is not a universal limit on valid inference at n plates and does not describe
    two-stage well resampling. Monte Carlo tail counts are checked separately. Neither
    passing this policy nor increasing draw count establishes nominal interval coverage.
    """
    if not isinstance(n_clusters, (int, np.integer)) or isinstance(n_clusters, bool) or n_clusters < 1:
        raise ValueError("n_clusters must be a positive integer")
    return float(min(1.0, 2.0 * (1.0 / n_clusters) ** n_clusters))


def adjusted_interval(log_draws: np.ndarray, point: float, n_clusters: int,
                      alpha: float, *, resampling: str = "cluster_only") -> tuple[float, float]:
    """Approximate t-inflated percentile bounds under an explicit resampling method.

    This helper accepts supplied log-fold draws; the producer uses the shared
    ``fold_change`` estimator directly. A t/normal quantile ratio is a small-sample
    adjustment, not proof of calibrated coverage.

    Args:
        log_draws: At least two finite log-fold draws; their logs need not be positive.
        point: Finite, positive fold estimate used as the log-space centre.
        n_clusters: Independent contributing plates, setting the degrees of freedom.
        alpha: Two-sided level; ``0.05 / m`` requests nominal Bonferroni bounds.
        resampling: ``cluster_only`` applies the empirical atom policy;
            ``two_stage`` does not. Both require enough draws in each requested tail.

    Raises:
        ValueError: For invalid inputs, fewer than two clusters, insufficient Monte
            Carlo tail draws, or a cluster-only alpha below the reporting policy.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    finest = finest_resolvable_alpha(n_clusters)
    if resampling not in ("cluster_only", "two_stage"):
        raise ValueError("resampling must be 'cluster_only' or 'two_stage'")
    if n_clusters < 2:
        raise ValueError("an interval needs at least two independent clusters")
    if resampling == "cluster_only" and alpha < finest:
        raise ValueError(
            f"alpha {alpha:.5g} is finer than the cluster-only tail-resolution policy "
            f"{finest:.5g} at {n_clusters} clusters. This empirical atom rule does not "
            "apply to two-stage well resampling, a parametric model, or other tests; "
            "none of those automatically establishes nominal coverage either")
    log_draws = np.asarray(log_draws, dtype=float)
    if log_draws.ndim != 1 or log_draws.size < 2 or not np.isfinite(log_draws).all():
        raise ValueError("log_draws must contain at least two finite log ratios")
    if not np.isfinite(point) or point <= 0:
        raise ValueError("point must be finite and positive")
    if log_draws.size * alpha / 2 < 1:
        raise ValueError("too few Monte Carlo draws to resolve both requested tails")
    low, high = np.quantile(log_draws, [alpha / 2.0, 1.0 - alpha / 2.0])
    inflation = float(stats.t.ppf(1.0 - alpha / 2.0, df=max(n_clusters - 1, 1))
                      / stats.norm.ppf(1.0 - alpha / 2.0))
    centre = float(np.log(point))
    return (float(np.exp(centre + (low - centre) * inflation)),
            float(np.exp(centre + (high - centre) * inflation)))


@dataclass(frozen=True)
class MultipleTests:
    """A fold family with labelled interval calls and its exact sign-test floor."""

    n_tests: int
    n_clusters: int
    unadjusted_calls: int
    family_wise_calls: int
    exact_p_floor: float
    table: pd.DataFrame
    exact_p_floor_method: str = "two-sided exact sign test on non-tied independent clusters"
    family_wise_method: str = "Bonferroni simultaneous intervals"

    def summary(self) -> str:
        return (f"{self.n_tests} tests on {self.n_clusters} plates: "
                f"{self.unadjusted_calls} clear 1.0 unadjusted, "
                f"{self.family_wise_calls} after {self.family_wise_method}; the exact sign test "
                f"on {self.n_clusters} non-tied clusters cannot go below p = {self.exact_p_floor:.3g}. "
                "This is not a resolution bound for parametric or within-plate permutation tests")


# --------------------------------------------------------------------------- the right test


@dataclass(frozen=True)
class DoseResponse:
    """One sensor's dose response, tested once rather than a dose at a time."""

    construct: str
    stressor: str
    n_doses: int
    n_plates: int
    statistic: float
    p_value: float
    n_permutations: int
    method: str = "two-sided within-plate dose-label permutation; plus-one Monte Carlo p-value"
    exchangeability: str = "dose assignments exchangeable within each plate under the null; cutoff prespecified"
    n_informative_plates: int = 0

    def summary(self) -> str:
        return (f"{self.construct} / {self.stressor}: plate-averaged Spearman "
                f"{self.statistic:+.3f} over {self.n_doses} doses, "
                f"permutation p = {self.p_value:.2g}")


def dose_response_permutation(
    readings: pd.DataFrame,
    construct: str,
    *,
    value: str = "activity_late",
    max_dose: float | None = None,
    n_permutations: int = 20000,
    seed: int = 0,
) -> DoseResponse:
    """Test within-plate monotone dose association, not the fold at an individual dose.

    Spearman's rho is computed within each informative plate and averaged, avoiding
    pooled plate offsets. Dose labels are shuffled within plates and the two-sided
    Monte Carlo p-value includes a plus-one correction. Its assignment space differs
    from the exact sign test, so :func:`sign_test_floor` does not bound it.

    Interpretation requires exchangeable dose assignments within each plate under the
    null and a prespecified dose cutoff. This association test does not identify a causal
    mechanism or establish practical equivalence, response size, or generalization to
    new plates. Constant-activity plates carry no rank-association information.

    Args:
        readings: Long-form, with ``plate``, ``construct``, ``dose_mM`` and ``value``.
        construct: Which reporter.
        value: Column to correlate against dose.
        max_dose: Prespecified upper dose cutoff, not one selected for significance.
            A nonmonotone response may be poorly described by the rank statistic.
        n_permutations: Monte Carlo draws.
        seed: Random seed.

    Raises:
        ValueError: For missing columns, invalid draws, an absent construct, fewer
            than two retained doses, or no informative plate.
    """
    required = {"plate", "construct", "dose_mM", value}
    if required - set(readings):
        raise ValueError(f"readings is missing required columns: {sorted(required - set(readings))}")
    if not isinstance(n_permutations, (int, np.integer)) or n_permutations < 1:
        raise ValueError("n_permutations must be a positive integer")
    subset = readings[readings.construct == construct]
    if subset.empty:
        raise ValueError(f"no rows for construct {construct!r}")
    if max_dose is not None:
        subset = subset[subset.dose_mM <= max_dose]
    subset = subset.reset_index(drop=True)
    if subset.dose_mM.nunique() < 2:
        raise ValueError(
            f"construct {construct!r} has {subset.dose_mM.nunique()} distinct dose(s) at or "
            f"below {max_dose}; a trend needs at least two")

    if subset.plate.isna().any() or not np.isfinite(subset[["dose_mM", value]].to_numpy(dtype=float)).all():
        raise ValueError("within-plate permutation requires nonmissing plates and finite doses and responses")
    if "stressor" in subset and subset.stressor.nunique() > 1:
        raise ValueError("one permutation test cannot pool different stressors")
    blocks = [group.index.to_numpy() for _, group in subset.groupby("plate")]
    observed_values = subset[value].to_numpy(dtype=float)
    doses = subset.dose_mM.to_numpy(dtype=float)

    def statistic(assignment: np.ndarray) -> float:
        rhos = [(stats.spearmanr(assignment[block], observed_values[block]).statistic
                 if np.ptp(observed_values[block]) > 0 else 0.0)
                for block in blocks if len(np.unique(assignment[block])) > 1]
        return float(np.mean(rhos)) if rhos else float("nan")

    observed = statistic(doses)
    if not np.isfinite(observed):
        # This returned p = 1/(n+1) -- the SMALLEST p the test can produce -- until
        # 2026-08-30. `statistic` yields NaN when no plate carries more than one dose, and
        # `np.abs(null) >= np.nan` is False for every draw, so the count was zero and the
        # most significant possible answer came back for a question that was never asked.
        # Found by the 2026-08-30 audit. A test whose failure mode is "maximally
        # significant" is worse than one that has no failure mode at all.
        raise ValueError(
            f"construct {construct!r}: the within-plate rank statistic is undefined -- no "
            f"plate of the {len(blocks)} carries more than one distinct dose at or below "
            f"{max_dose}, so there is nothing to correlate within a plate. Permuting dose "
            "labels inside a plate cannot test a ladder that does not exist inside one")
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for draw in range(n_permutations):
        shuffled = doses.copy()
        for block in blocks:
            shuffled[block] = rng.permutation(doses[block])
        null[draw] = statistic(shuffled)

    # The +1 is not a rounding convenience. A permutation p of exactly zero claims the
    # observed statistic is unreachable under the null, which no finite set of draws can
    # establish -- the observed assignment is itself one of the permutations.
    if not np.isfinite(null).all():
        raise ValueError(
            f"construct {construct!r}: {int((~np.isfinite(null)).sum())} of "
            f"{n_permutations} permutation draws gave an undefined statistic. Counting only "
            "the finite ones would divide by a denominator the null does not have")
    p_value = float((1 + np.sum(np.abs(null) >= abs(observed))) / (n_permutations + 1))
    return DoseResponse(
        construct=construct, stressor=str(subset.stressor.iloc[0]) if "stressor" in subset else "unspecified",
        n_doses=int(subset.dose_mM.nunique()), n_plates=len(blocks),
        statistic=observed, p_value=p_value, n_permutations=n_permutations,
        n_informative_plates=sum(np.unique(doses[block]).size > 1 for block in blocks))
