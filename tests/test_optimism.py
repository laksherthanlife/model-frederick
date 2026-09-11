"""The bound has to be the right number, not merely a non-negative one.

A bias estimator is easy to write and hard to check, because on a real problem there is
nothing to check it against. So the suite checks it against a problem where the answer is
known in closed form: configurations drawn from a standard normal, a solution that is a single
number, and a fit that returns the mean of the configurations it saw. There the expected
optimality gap is exactly ``1/n_candidate + 1/n_reference_domains`` -- the candidate's own
sampling variance plus the variance of the reference it is being compared against -- so the
estimator can be scored rather than merely exercised.

Three properties beyond that value are load-bearing and each has a test that would fail if it
broke. The gap samples must be paired by seed, or the difference between two solutions is
swamped by the difference between two draws. Negative gaps must be clipped and counted, not
averaged in, because a reference that landed worse than the candidate says something about the
optimiser and nothing about the bias. And the bound must shrink with more configurations *in
expectation only* -- there is a test that it does on average and a test that it does not always,
because a report that reads a single-run rise as a regression has misread the theory.
"""

from __future__ import annotations

import numpy as np
import pytest

from ystwin.analysis.optimism import (
    OptimismBound,
    OptimismProblem,
    basic_bootstrap_upper,
    estimate_optimism,
)


def quadratic(spread: float = 1.0) -> OptimismProblem:
    """Configurations on a line, and a fitter that is exactly optimal on what it saw.

    Chosen because the optimality gap has a closed form: the mean minimises squared error on
    its own configurations, so the reference is the true in-sample optimum and the gap is
    ``(candidate - reference)**2``, whose expectation is the sum of the two sampling
    variances. That makes the estimator falsifiable rather than merely runnable.
    """
    return OptimismProblem(
        sample=lambda rng: float(spread * rng.normal()),
        fit=lambda domains, seed, initial: float(np.mean(domains)),
        score=lambda solution, domains, seed: float(
            -np.mean([(solution - d) ** 2 for d in domains])),
    )


def seeded_only() -> OptimismProblem:
    """A scorer that reads its seed and ignores its solution.

    The instrument for the pairing test: with common random numbers every gap must be exactly
    zero, and without them the gaps would be the difference between two random draws.
    """
    return OptimismProblem(
        sample=lambda rng: float(rng.normal()),
        fit=lambda domains, seed, initial: float(np.mean(domains)),
        score=lambda solution, domains, seed: float(
            np.random.default_rng(seed).normal()),
    )


class TestTheBootstrapBound:
    def test_constant_gaps_bound_at_their_own_value(self):
        """Every resample has the same mean, so reflection returns the mean untouched."""
        assert basic_bootstrap_upper(np.full(20, 0.25)) == pytest.approx(0.25)

    def test_it_lies_at_or_above_the_mean_gap(self):
        rng = np.random.default_rng(0)
        for _ in range(40):
            gaps = np.abs(rng.normal(size=int(rng.integers(2, 60))))

            assert basic_bootstrap_upper(gaps, seed=1) >= gaps.mean() - 1e-9

    def test_non_negative_gaps_give_a_non_negative_bound(self):
        """The sign the theory demands: the bias is non-negative by Jensen's inequality."""
        rng = np.random.default_rng(1)
        for _ in range(40):
            gaps = np.clip(rng.normal(size=30), 0.0, None)

            assert basic_bootstrap_upper(gaps, seed=2) >= 0.0

    def test_a_skewed_sample_bounds_well_above_its_mean(self):
        """Mostly zeros with one large gap: the atom at zero is exactly the case a
        standard-error interval would mishandle, and the basic form still widens."""
        gaps = np.array([0.0] * 9 + [10.0])

        assert basic_bootstrap_upper(gaps, seed=0) > 1.5 * gaps.mean()

    def test_it_is_reproducible_from_its_seed(self):
        gaps = np.abs(np.random.default_rng(3).normal(size=25))

        assert basic_bootstrap_upper(gaps, seed=5) == basic_bootstrap_upper(gaps, seed=5)
        assert basic_bootstrap_upper(gaps, seed=5) != basic_bootstrap_upper(gaps, seed=6)

    def test_one_sample_is_refused(self):
        with pytest.raises(ValueError, match="at least 2 gap samples"):
            basic_bootstrap_upper([0.4])

    def test_a_non_finite_gap_is_refused(self):
        with pytest.raises(ValueError, match="finite"):
            basic_bootstrap_upper([0.1, np.nan, 0.3])

    def test_a_two_sided_level_is_refused(self):
        with pytest.raises(ValueError, match="one-sided"):
            basic_bootstrap_upper([0.1, 0.2, 0.3], alpha=0.5)

    def test_too_few_replications_for_the_quantile_are_refused(self):
        """At ten replications the 5% quantile is the smallest resample mean, not a quantile."""
        with pytest.raises(ValueError, match="cannot resolve"):
            basic_bootstrap_upper([0.1, 0.2, 0.3], alpha=0.05, n_bootstrap=10)


class TestTheEstimatorRecoversTheKnownGap:
    def test_it_matches_the_closed_form_optimality_gap(self):
        """Mean gap should be ``1/n_candidate + 1/n_reference_domains`` on the toy problem.

        The candidate's own sampling variance plus the reference's. Averaged over seeds
        because a single run of twenty references is a small sample of a heavy-tailed
        quantity, which is the reason the bound exists in the first place.
        """
        gaps = [estimate_optimism(quadratic(), n_candidate_domains=8,
                                  n_reference_domains=8, n_reference=40, seed=s).mean_gap
                for s in range(30)]

        assert float(np.mean(gaps)) == pytest.approx(1 / 8 + 1 / 8, rel=0.25)

    def test_the_bound_is_non_negative_and_covers_the_mean_gap(self):
        result = estimate_optimism(quadratic(), n_candidate_domains=4, n_reference=20, seed=0)

        assert result.bound >= result.mean_gap >= 0.0

    def test_it_is_reproducible_from_its_seed(self):
        first = estimate_optimism(quadratic(), n_candidate_domains=3, n_reference=10, seed=4)
        again = estimate_optimism(quadratic(), n_candidate_domains=3, n_reference=10, seed=4)

        assert first.bound == again.bound
        assert np.array_equal(first.raw_gaps, again.raw_gaps)

    def test_the_seed_actually_reaches_the_sampling(self):
        """The other half of what `scripts/audit_determinism.py` asks of a seeded entry
        point: a seed that is accepted and never drawn from reproduces perfectly and makes
        any spread computed across seeds zero by construction."""
        bounds = {estimate_optimism(quadratic(), n_candidate_domains=3, n_reference=10,
                                    seed=s).bound for s in range(4)}

        assert len(bounds) == 4

    def test_it_carries_the_budget_that_produced_it(self):
        result = estimate_optimism(quadratic(), n_candidate_domains=5,
                                   n_reference_domains=3, n_reference=11, seed=0)

        assert (result.n_candidate_domains, result.n_reference_domains) == (5, 3)
        assert result.n_reference == 11
        assert "UCBOG" in result.summary() and "11 references" in result.summary()


class TestItShrinksInExpectationAndNotOtherwise:
    def test_more_candidate_configurations_lower_the_bound_on_average(self):
        """Mak, Morton & Wood Theorem 2, which is the reason to sample families at all."""
        few = [estimate_optimism(quadratic(), n_candidate_domains=1, n_reference=20,
                                 seed=s).bound for s in range(12)]
        many = [estimate_optimism(quadratic(), n_candidate_domains=16, n_reference=20,
                                  seed=s).bound for s in range(12)]

        assert float(np.mean(many)) < 0.8 * float(np.mean(few))

    def test_a_single_run_can_still_rise(self):
        """Pinned because the docstring promises it and a report must not read a rise as a
        regression: the monotonicity is in expectation, not per run."""
        few = [estimate_optimism(quadratic(), n_candidate_domains=1, n_reference=20,
                                 seed=s).bound for s in range(12)]
        many = [estimate_optimism(quadratic(), n_candidate_domains=16, n_reference=20,
                                  seed=s).bound for s in range(12)]

        assert any(m > f for m, f in zip(many, few))


class TestCommonRandomNumbers:
    def test_the_candidate_and_the_reference_are_scored_on_the_same_draw(self):
        """The pairing, tested directly: same seed and the same configurations, both calls."""
        calls = []
        problem = OptimismProblem(
            sample=lambda rng: float(rng.normal()),
            fit=lambda domains, seed, initial: float(np.mean(domains)),
            score=lambda solution, domains, seed: (
                calls.append((tuple(domains), seed)) or 0.0),
        )
        estimate_optimism(problem, n_candidate_domains=2, n_reference_domains=3,
                          n_reference=5, seed=0)

        assert len(calls) == 10
        for reference, candidate in zip(calls[0::2], calls[1::2]):
            assert reference == candidate

    def test_a_scorer_that_only_reads_its_seed_reports_no_gap_at_all(self):
        """If the seeds were not synchronised these gaps would be two unrelated draws."""
        result = estimate_optimism(seeded_only(), n_candidate_domains=3, n_reference=20,
                                   seed=0)

        assert np.all(result.raw_gaps == 0.0)
        assert result.bound == 0.0

    def test_a_reference_is_initialised_from_the_candidate(self):
        """SPOTA requires it: a reference that starts elsewhere can land in a worse local
        optimum and report a negative gap, which is noise about the optimiser."""
        seen = []
        problem = OptimismProblem(
            sample=lambda rng: float(rng.normal()),
            fit=lambda domains, seed, initial: (seen.append(initial)
                                                or float(np.mean(domains))),
            score=lambda solution, domains, seed: float(
                -np.mean([(solution - d) ** 2 for d in domains])),
        )
        estimate_optimism(problem, n_candidate_domains=2, n_reference=4, seed=0)

        assert seen[0] is None
        assert all(start is not None for start in seen[1:])
        assert len(set(seen[1:])) == 1

    def test_a_reference_is_fitted_and_scored_on_the_same_configurations(self):
        """The in-sample optimum is the object the theory differences; scoring a reference
        out of sample would estimate something else."""
        fitted, scored = [], []
        problem = OptimismProblem(
            sample=lambda rng: float(rng.normal()),
            fit=lambda domains, seed, initial: (fitted.append(tuple(domains))
                                                or float(np.mean(domains))),
            score=lambda solution, domains, seed: (scored.append(tuple(domains)) or 0.0),
        )
        estimate_optimism(problem, n_candidate_domains=2, n_reference_domains=2,
                          n_reference=4, seed=0)

        assert fitted[1:] == scored[0::2]


class TestClipping:
    def test_a_reference_that_lands_worse_is_clipped_and_counted(self):
        """A deliberately bad reference fitter: every raw gap is negative, so the bound is
        zero and the clipped fraction says why rather than hiding it."""
        problem = OptimismProblem(
            sample=lambda rng: float(rng.normal()),
            fit=lambda domains, seed, initial: (0.0 if initial is None else 50.0),
            score=lambda solution, domains, seed: float(
                -np.mean([(solution - d) ** 2 for d in domains])),
        )
        result = estimate_optimism(problem, n_candidate_domains=4, n_reference=20, seed=0)

        assert np.all(result.raw_gaps < 0.0)
        assert np.all(result.gaps == 0.0)
        assert result.bound == 0.0
        assert result.clipped_fraction == 1.0

    def test_the_raw_gaps_survive_the_clipping(self):
        """Kept because their negative fraction is the diagnostic, not an implementation
        detail: it distinguishes "the candidate is already good" from "the fits are noise"."""
        result = estimate_optimism(quadratic(), n_candidate_domains=4, n_reference=20, seed=0)

        assert result.raw_gaps.shape == result.gaps.shape
        assert np.all(result.gaps >= 0.0)


class TestItRefusesRatherThanGuesses:
    def test_one_reference_cannot_be_bootstrapped(self):
        with pytest.raises(ValueError, match="fewer than 2 gap samples"):
            estimate_optimism(quadratic(), n_reference=1)

    def test_a_candidate_needs_something_to_be_fitted_on(self):
        with pytest.raises(ValueError, match="at least one configuration"):
            estimate_optimism(quadratic(), n_candidate_domains=0)

    def test_a_reference_needs_something_to_be_fitted_on(self):
        with pytest.raises(ValueError, match="at least one configuration"):
            estimate_optimism(quadratic(), n_reference_domains=0)

    def test_a_non_finite_score_is_refused_rather_than_dropped(self):
        """Silently dropping it would bias the gap toward whichever references happened to
        evaluate, which is the failure this whole module is about."""
        problem = OptimismProblem(
            sample=lambda rng: float(rng.normal()),
            fit=lambda domains, seed, initial: 0.0,
            score=lambda solution, domains, seed: float("nan"),
        )
        with pytest.raises(ValueError, match="non-finite score"):
            estimate_optimism(problem, n_reference=4)


def test_the_result_reports_its_own_clipping_on_an_empty_sample():
    """A guard on the property rather than on the estimator: `nan` and not a crash."""
    empty = OptimismBound(bound=0.0, mean_gap=0.0, gaps=np.array([]), raw_gaps=np.array([]))

    assert np.isnan(empty.clipped_fraction)
    assert empty.n_reference == 0


def test_mutating_a_reference_start_cannot_modify_the_frozen_candidate():
    def fit(domains, seed, initial):
        if initial is None:
            return {"value": 0.0}
        initial["value"] += 1.0
        return initial

    problem = OptimismProblem(sample=lambda rng: float(rng.normal()), fit=fit,
                              score=lambda solution, domains, seed: solution["value"])
    result = estimate_optimism(problem, n_reference=5, seed=1)
    assert result.candidate == {"value": 0.0}
    assert result.raw_gaps == pytest.approx(np.ones(5))


def test_invalid_bootstrap_budget_is_refused_before_any_fits_or_sampling():
    sampled = []
    problem = OptimismProblem(sample=lambda rng: sampled.append(1) or 0.0,
                              fit=lambda domains, seed, initial: 0.0,
                              score=lambda solution, domains, seed: 0.0)
    with pytest.raises(ValueError, match="cannot resolve"):
        estimate_optimism(problem, n_bootstrap=5)
    assert sampled == []


def test_zero_gaps_carry_a_coverage_warning_not_a_certificate():
    result = estimate_optimism(seeded_only(), n_reference=20)
    assert result.bound == 0.0
    assert "rare unseen gaps" in result.coverage_warning
    assert "not empirically coverage-calibrated" in result.coverage_status
    assert "population optimality gap" in result.target


def test_basic_bootstrap_can_severely_undercover_rare_nonzero_gaps():
    rng = np.random.default_rng(2)
    trials, hits = 200, 0
    for trial in range(trials):
        samples = rng.binomial(1, 0.01, 20).astype(float)
        hits += int(basic_bootstrap_upper(samples, n_bootstrap=200, seed=trial) >= 0.01)
    actual_coverage = hits / trials
    assert actual_coverage < 0.4
    assert actual_coverage < 0.95
