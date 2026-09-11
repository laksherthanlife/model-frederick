"""An approximate optimality-gap bound under an assumed simulation distribution.

Muratore, Gienger & Peters (IEEE TPAMI 43(4):1172-1183, 2021) define Simulation Optimization
Bias as E[max J_n] - max E[J]. It is non-negative under unbiased sampling and well-defined
optima, IN EXPECTATION over samples, not for every fitted score. The expected sample optimum
is monotone with iid sample size under the assumptions of Mak, Morton & Wood (Oper. Res.
Lett. 24(1-2):47-56, 1999). Neither result fixes the sign of a simulator-versus-biology gap
or guarantees monotonicity of this finite bootstrap estimate.

For a candidate independent of the reference sample, E[estimated gap] = population gap +
optimization bias. Thus the gap is not the bias alone. This module implements the SPOTA-style
procedure of Muratore, Treede, Gienger & Peters (CoRL 2018, PMLR 87:700-713): fit a candidate,
fit references on fresh configurations, compare under common random numbers, and bootstrap
a one-sided bound. A local reference fit is not necessarily the global sample optimum;
clipping negative gaps does not restore that guarantee.

The target is the candidate's expected optimality gap under the declared distribution, not
its loss on each fresh draw and not transfer to unknown biology. In `scripts/run_cross_family.py`
the candidate is a dimension/design choice, not the fitted latent weights. The sampler there
is uniform over `generator.families.sample_family`'s registry.

The basic bootstrap used here is an approximation. Sparse nonzero gaps, clipping and an atom
at zero do not themselves establish confidence coverage. Report the objective, distribution,
raw gaps and sampling budget beside the bound; a small value does not validate a real plate.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "ALPHA",
    "N_BOOTSTRAP",
    "N_REFERENCE",
    "OptimismBound",
    "OptimismProblem",
    "basic_bootstrap_upper",
    "estimate_optimism",
]

ALPHA = 0.05
"""One-sided tail probability for the nominal confidence bound, as published."""

N_BOOTSTRAP = 1000
"""Bootstrap replications, as published."""

N_REFERENCE = 20
"""Reference solutions, as published. Each costs a full fit, so a caller on a budget will
lower it -- and should say that it did, because the quantile is what pays for it."""


@dataclass(frozen=True)
class OptimismProblem:
    """The three callables the estimator needs, and the contract each has to keep.

    Kept abstract on purpose. The estimator is a statement about sample-based optimisation
    and holds whatever is being optimised over whatever is being sampled, so tying it to
    this package's latent fit would make it untestable except through the thing under test.

    Args:
        sample: ``(rng) -> domain``. One draw from the assumed distribution over
            configurations. Must consume `rng` for anything random, so the whole procedure
            is reproducible from one seed.
        fit: ``(domains, seed, initial) -> solution``. Optimise on those configurations.
            `initial` is the candidate, passed so a reference starts where the candidate
            already is: SPOTA requires it because a reference that lands in a worse local
            optimum reports a negative gap, which is noise about the optimiser rather than
            information about the bias. An exhaustive search over a small discrete set can
            ignore it beyond tie-breaking, and should say so.
        score: ``(solution, domains, seed) -> float``, larger is better. Called once for
            the reference and once for the candidate with the *same* seed and the *same*
            configurations, which is the common-random-numbers step: same wells, same noise
            draws, so the difference is the solutions and not the draw. A `score` that
            ignores its seed silently discards that variance reduction.
    """

    sample: Callable[[np.random.Generator], Any]
    fit: Callable[[Sequence[Any], int, Any], Any]
    score: Callable[[Any, Sequence[Any], int], float]


@dataclass(frozen=True)
class OptimismBound:
    """An approximate bound on a candidate's expected optimality gap, not each future loss.

    Args:
        bound: The UCBOG estimate. Nominal coverage is ``1 - alpha`` under the stated
            sampling, reference-optimization and bootstrap assumptions.
        mean_gap: Mean of the clipped gap samples, the point estimate the bound is built on.
        gaps: Clipped gap samples, one per reference.
        raw_gaps: The same samples before clipping, kept because their negative fraction is
            the diagnostic that says whether the references were actually better optimisers
            than the candidate or just differently seeded.
        candidate: The solution the bound is about.
        alpha, n_bootstrap, n_candidate_domains, n_reference_domains: Settings, carried so a
            reported bound cannot be separated from the budget that produced it.
    """

    bound: float
    mean_gap: float
    gaps: np.ndarray = field(repr=False)
    raw_gaps: np.ndarray = field(repr=False)
    candidate: Any = None
    alpha: float = ALPHA
    n_bootstrap: int = N_BOOTSTRAP
    n_candidate_domains: int = 0
    n_reference_domains: int = 0
    target: str = "candidate population optimality gap under the declared simulator distribution"
    method: str = "basic bootstrap upper bound on clipped reference sample-optimum gaps"
    coverage_status: str = "nominal approximation, not empirically coverage-calibrated or distribution-free"
    assumptions: tuple[str, ...] = (
        "independent identically distributed reference configuration sets",
        "fixed candidate independent of all reference samples",
        "globally optimal reference fits needed for the population-gap bound interpretation",
        "sample gaps include reference optimization bias; clipping does not establish coverage",
    )

    @property
    def coverage_warning(self) -> str:
        if self.clipped_fraction > 0:
            return "reference fits scored worse than the candidate; clipping is not a certificate of optimality"
        if self.gaps.size and np.all(self.gaps == 0):
            return "all observed gaps are zero; rare unseen gaps can make bootstrap coverage arbitrarily poor"
        return self.coverage_status

    @property
    def n_reference(self) -> int:
        return int(self.gaps.size)

    @property
    def clipped_fraction(self) -> float:
        """Share of references that scored *below* the candidate on their own domains.

        A high fraction may reflect weak local reference optimization, not a good candidate.
        Clipping hides those negative differences from the bootstrap but preserves them here
        for diagnosis; it does not establish the nominal confidence coverage.
        """
        return float(np.mean(self.raw_gaps < 0.0)) if self.raw_gaps.size else float("nan")

    def summary(self) -> str:
        return (
            f"UCBOG {self.bound:+.4f} at alpha={self.alpha} "
            f"(mean gap {self.mean_gap:+.4f} over {self.n_reference} references, "
            f"{self.clipped_fraction:.0%} clipped; {self.n_candidate_domains} candidate and "
            f"{self.n_reference_domains} reference configurations, "
            f"{self.n_bootstrap} bootstrap draws; nominal {1 - self.alpha:.1%}). "
            f"Target: {self.target}. {self.coverage_warning}"
        )


def _bootstrap_budget(alpha, n_bootstrap):
    if not 0.0 < alpha < 0.5:
        raise ValueError(f"alpha={alpha} is not a one-sided level in (0, 0.5)")
    if not isinstance(n_bootstrap, (int, np.integer)) or isinstance(n_bootstrap, bool) or n_bootstrap < 1:
        raise ValueError("n_bootstrap must be a positive integer")
    if n_bootstrap * alpha < 1.0:
        raise ValueError(f"{n_bootstrap} replications cannot resolve the {alpha} quantile; "
                         f"use at least {int(np.ceil(1.0 / alpha))}")


def basic_bootstrap_upper(
    gaps,
    alpha: float = ALPHA,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = 0,
) -> float:
    """One-sided upper bound on the mean gap: ``2 * mean - Q_alpha[bootstrap means]``.

    This implements the basic bootstrap form used by the SPOTA-style procedure, not a
    distribution-free guarantee. Like other bootstrap intervals it needs a representative
    sampling distribution; clipping and an atom at zero do not automatically validate it.

    The result is the unmodified reflection of the sampled quantile about the observed mean.
    A finite resampling budget does not guarantee that a lower-tail empirical quantile is
    below that mean, so this helper does not promise a bound above the observed mean.

    Args:
        gaps: Gap samples, one per reference solution.
        alpha: One-sided tail probability; nominal coverage is ``1 - alpha``, subject
            to the sampling and bootstrap assumptions.
        n_bootstrap: Resampling replications.
        seed: Seed for the resampling.

    Raises:
        ValueError: if there is nothing to bootstrap, or if `alpha` and `n_bootstrap` cannot
            between them resolve the quantile being asked for.
    """
    samples = np.asarray(gaps, dtype=float)
    if samples.ndim != 1 or samples.size < 2:
        raise ValueError(
            f"need at least 2 gap samples to bootstrap, got {samples.size}; one reference "
            "gives a point and not a distribution")
    if not np.all(np.isfinite(samples)):
        raise ValueError("gap samples must all be finite; a failed evaluation is not a gap")
    _bootstrap_budget(alpha, n_bootstrap)

    rng = np.random.default_rng(seed)
    draws = rng.choice(samples, size=(n_bootstrap, samples.size), replace=True)
    quantile = float(np.quantile(draws.mean(axis=1), alpha))
    return float(2.0 * samples.mean() - quantile)


def estimate_optimism(
    problem: OptimismProblem,
    n_candidate_domains: int = 4,
    n_reference_domains: int = 2,
    n_reference: int = N_REFERENCE,
    alpha: float = ALPHA,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = 0,
) -> OptimismBound:
    """Run SPOTA: fit a candidate, fit references, bootstrap the gap they leave.

    The four published steps, in order. A candidate is fitted on `n_candidate_domains` sampled
    configurations. Then `n_reference` reference solutions are each fitted on their own fresh
    `n_reference_domains`, each initialised from the candidate. Both are scored on each
    reference's own configurations under the same seed, and the difference is that reference's
    gap sample. Negative samples are clipped to zero, because a reference is only locally
    optimal and a negative gap is a statement about the optimiser rather than about the bias.
    The clipped samples are bootstrapped into a one-sided upper bound.

    Scoring a reference on its fitting configurations estimates a sample optimality gap.
    For a globally optimized reference and an independent fixed candidate its expectation
    is the candidate's population optimality gap PLUS the reference sample's optimization
    bias, not the bias alone. Local fits and clipping qualify this interpretation.
    Evaluating the reference out of sample would estimate a different quantity.

    Args:
        problem: The sampler, fitter and scorer; see `OptimismProblem` for the contract.
        n_candidate_domains: Configurations used to fit the candidate. Increasing this
            changes its training budget, not a guarantee of smaller gap or bootstrap bound.
        n_reference_domains: Configurations behind each reference. This sample size controls
            the reference sample-optimum bias under the theorem's assumptions.
        n_reference: Reference solutions, hence gap samples.
        alpha: One-sided level for the bound.
        n_bootstrap: Bootstrap replications.
        seed: Seed for sampling configurations and for every evaluation seed drawn from it.

    Raises:
        ValueError: on a budget that cannot produce a bound.
    """
    _bootstrap_budget(alpha, n_bootstrap)
    if any(not isinstance(value, (int, np.integer)) or isinstance(value, bool)
           for value in (n_candidate_domains, n_reference_domains, n_reference)):
        raise ValueError("configuration and reference counts must be integers")
    if n_candidate_domains < 1:
        raise ValueError("the candidate needs at least one configuration to be fitted on")
    if n_reference_domains < 1:
        raise ValueError("a reference needs at least one configuration to be fitted on")
    if n_reference < 2:
        raise ValueError(
            f"n_reference={n_reference} gives fewer than 2 gap samples, which cannot be "
            "bootstrapped")

    rng = np.random.default_rng(seed)
    candidate_domains = [problem.sample(rng) for _ in range(n_candidate_domains)]
    candidate = problem.fit(candidate_domains, int(rng.integers(0, 2**31 - 1)), None)

    raw = np.empty(n_reference, dtype=float)
    for i in range(n_reference):
        domains = [problem.sample(rng) for _ in range(n_reference_domains)]
        # One seed per reference, used for its fit and for both evaluations on it: the
        # candidate and the reference then meet the same wells and the same noise.
        common = int(rng.integers(0, 2**31 - 1))
        reference = problem.fit(deepcopy(domains), common, deepcopy(candidate))
        raw[i] = (float(problem.score(deepcopy(reference), deepcopy(domains), common))
                  - float(problem.score(deepcopy(candidate), deepcopy(domains), common)))
    if not np.all(np.isfinite(raw)):
        raise ValueError(
            "some reference produced a non-finite score, so its gap is undefined; a scorer "
            "must return a number or raise")

    gaps = np.clip(raw, 0.0, None)
    return OptimismBound(
        bound=basic_bootstrap_upper(gaps, alpha, n_bootstrap, seed=seed),
        mean_gap=float(gaps.mean()),
        gaps=gaps,
        raw_gaps=raw,
        candidate=candidate,
        alpha=alpha,
        n_bootstrap=n_bootstrap,
        n_candidate_domains=n_candidate_domains,
        n_reference_domains=n_reference_domains,
    )
