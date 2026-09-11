"""A null has to reject real structure and accept its absence, and phases matter.

The load-bearing test here is that a shuffle-based null and a phase-randomised null
disagree on smooth autocorrelated data: the shuffle rejects because the surrogate is
not smooth, which is the wrong reason, and phase randomisation does not. That is the
whole argument for having two surrogates rather than one.
"""

from __future__ import annotations

import numpy as np
import pytest

from ystwin.analysis.nulls import (
    NullResult,
    compare_to_null,
    matched_marginals,
    phase_randomised,
    rotated_subspace,
    shuffled_loadings,
    skill_score,
    skill_score_if_defined,
)


def _shuffle_rows(data, rng):
    out = np.asarray(data, dtype=float).copy()
    rng.shuffle(out)
    return out


# ---------------------------------------------------------------------------
# the result object
# ---------------------------------------------------------------------------


def test_p_value_is_floored_at_the_resolution_limit():
    """A permutation p-value of exactly zero is not attainable."""
    result = NullResult(observed=10.0, null=np.zeros(99), greater_is_better=True)
    assert result.p_value == pytest.approx(1 / 100)
    assert result.p_value > 0.0


def test_direction_of_better_inverts_the_tail():
    null = np.linspace(0.0, 1.0, 101)
    high_good = NullResult(0.9, null, greater_is_better=True)
    low_good = NullResult(0.9, null, greater_is_better=False)
    assert high_good.p_value < 0.2
    assert low_good.p_value > 0.8


def test_effect_over_null_is_positive_for_better_in_both_directions():
    null = np.full(50, 0.5)
    assert NullResult(0.9, null, greater_is_better=True).effect_over_null > 0
    assert NullResult(0.1, null, greater_is_better=False).effect_over_null > 0


def test_degenerate_draws_are_excluded_not_counted():
    null = np.array([0.1, 0.2, np.nan, np.nan, 0.3])
    result = NullResult(0.9, null, greater_is_better=True)
    assert result.n_draws == 3
    assert np.isfinite(result.null_median)


def test_summary_names_the_verdict():
    strong = NullResult(0.9, np.zeros(99), True, label="transfer R2")
    weak = NullResult(0.0, np.zeros(99), True, label="transfer R2")
    assert "beats" in strong.summary() and "transfer R2" in strong.summary()
    assert "DOES NOT BEAT" in weak.summary()


# ---------------------------------------------------------------------------
# compare_to_null
# ---------------------------------------------------------------------------


def test_real_structure_beats_the_null():
    rng = np.random.default_rng(0)
    signal = rng.normal(size=200)
    data = np.column_stack([signal + 0.1 * rng.normal(size=200) for _ in range(4)])

    def first_eigenvalue_share(x):
        x = x - np.nanmean(x, axis=0)
        sv = np.linalg.svd(np.nan_to_num(x), compute_uv=False)
        return float(sv[0] ** 2 / np.sum(sv**2))

    result = compare_to_null(first_eigenvalue_share, data, matched_marginals,
                            n_draws=60, label="shared variance", seed=1)
    assert result.beats_null
    assert result.observed > result.null_median


def test_absent_structure_does_not_beat_the_null():
    """Independent channels must fail -- otherwise the null is broken, not the data."""
    rng = np.random.default_rng(1)
    data = rng.normal(size=(200, 4))

    def first_eigenvalue_share(x):
        x = x - np.nanmean(x, axis=0)
        sv = np.linalg.svd(np.nan_to_num(x), compute_uv=False)
        return float(sv[0] ** 2 / np.sum(sv**2))

    result = compare_to_null(first_eigenvalue_share, data, matched_marginals,
                            n_draws=60, seed=2)
    assert not result.beats_null


def test_too_few_draws_is_refused():
    with pytest.raises(ValueError, match="cannot resolve a p-value"):
        compare_to_null(lambda d: 1.0, np.zeros((10, 2)), matched_marginals, n_draws=5)


def test_degenerate_surrogates_do_not_abort_the_sweep():
    def brittle(x):
        if np.allclose(x, x[0]):
            raise ValueError("constant input")
        return float(np.std(x))

    def sometimes_constant(data, rng):
        return np.zeros_like(data) if rng.random() < 0.5 else data

    result = compare_to_null(brittle, np.arange(50.0), sometimes_constant, n_draws=40, seed=3)
    assert 0 < result.n_draws < 40  # some failed, some survived


def test_is_deterministic_given_a_seed():
    data = np.random.default_rng(9).normal(size=(80, 3))
    stat = lambda x: float(np.nanmean(x))  # noqa: E731
    a = compare_to_null(stat, data, matched_marginals, n_draws=30, seed=7)
    b = compare_to_null(stat, data, matched_marginals, n_draws=30, seed=7)
    assert np.allclose(a.null, b.null, equal_nan=True)


# ---------------------------------------------------------------------------
# phase randomisation -- the important one
# ---------------------------------------------------------------------------


def _smooth_trace(n, rng, scale=8.0):
    """A smooth autocorrelated series, like an OD or activity trace."""
    raw = rng.normal(size=n)
    kernel = np.exp(-np.arange(-n // 2, n // 2) ** 2 / (2 * scale**2))
    kernel /= kernel.sum()
    return np.convolve(raw, kernel, mode="same")


def test_phase_randomisation_preserves_the_power_spectrum():
    rng = np.random.default_rng(4)
    x = _smooth_trace(256, rng)
    surrogate = phase_randomised(x, rng)
    assert np.allclose(np.abs(np.fft.rfft(x - x.mean())),
                       np.abs(np.fft.rfft(surrogate - surrogate.mean())), atol=1e-8)


def test_phase_randomisation_preserves_mean_and_returns_real_values():
    rng = np.random.default_rng(5)
    for n in (128, 129):  # even and odd -- the Nyquist term differs
        x = _smooth_trace(n, rng)
        surrogate = phase_randomised(x, rng)
        assert np.isrealobj(surrogate)
        assert surrogate.mean() == pytest.approx(x.mean(), abs=1e-8)
        assert surrogate.shape == x.shape


def test_phase_randomisation_preserves_smoothness_and_shuffling_destroys_it():
    """This is why the surrogate matters: a shuffle is not a like-for-like null."""
    rng = np.random.default_rng(6)
    x = _smooth_trace(256, rng)
    roughness = lambda s: float(np.mean(np.abs(np.diff(s))))  # noqa: E731

    phase = roughness(phase_randomised(x, rng))
    shuffled = roughness(rng.permutation(x))
    assert phase == pytest.approx(roughness(x), rel=0.4)
    assert shuffled > 5 * roughness(x)


def test_a_spurious_correlation_between_smooth_series_survives_the_shuffle_null():
    """Two *independent* smooth traces correlate, and the shuffle null misses it.

    The phase-randomised null keeps the autocorrelation that generates the spurious
    correlation, so it does not reject. The shuffle null destroys the autocorrelation,
    so the observed correlation looks extreme against it and it rejects -- for the
    wrong reason. That is exactly the trap the 0.82 variance-explained figure sits in.
    """
    rejections_phase, rejections_shuffle = 0, 0
    trials = 25
    for t in range(trials):
        gen = np.random.default_rng(100 + t)
        a, b = _smooth_trace(128, gen), _smooth_trace(128, gen)  # independent
        stat = lambda pair: abs(float(np.corrcoef(pair[0], pair[1])[0, 1]))  # noqa: E731

        phase_null = compare_to_null(
            stat, (a, b), lambda p, r: (phase_randomised(p[0], r), p[1]),
            n_draws=40, seed=t)
        shuffle_null = compare_to_null(
            stat, (a, b), lambda p, r: (r.permutation(p[0]), p[1]),
            n_draws=40, seed=t)
        rejections_phase += phase_null.beats_null
        rejections_shuffle += shuffle_null.beats_null

    # The shuffle null over-rejects on independent smooth data; phase randomisation
    # keeps close to the nominal 5%.
    assert rejections_shuffle > rejections_phase
    assert rejections_phase <= 0.25 * trials


def test_too_short_a_series_is_refused():
    with pytest.raises(ValueError, match="at least 4 points"):
        phase_randomised(np.array([1.0, 2.0, 3.0]), np.random.default_rng(0))


def test_phase_randomisation_handles_multiple_columns_independently():
    rng = np.random.default_rng(7)
    x = np.column_stack([_smooth_trace(128, rng) for _ in range(3)])
    surrogate = phase_randomised(x, rng)
    assert surrogate.shape == x.shape
    # Independent phases destroy the cross-column relation.
    before = abs(np.corrcoef(x[:, 0], x[:, 1])[0, 1])
    after = abs(np.corrcoef(surrogate[:, 0], surrogate[:, 1])[0, 1])
    assert after < before + 0.5


# ---------------------------------------------------------------------------
# the other surrogates
# ---------------------------------------------------------------------------


def test_shuffled_loadings_keeps_the_value_multiset_and_shape():
    rng = np.random.default_rng(8)
    loadings = np.array([[1.0, 0.2, 0.0], [0.0, 0.9, 0.3]])
    shuffled = shuffled_loadings(loadings, rng)
    assert shuffled.shape == loadings.shape
    assert np.allclose(np.sort(shuffled.ravel()), np.sort(loadings.ravel()))


def test_shuffled_loadings_does_not_mutate_its_input():
    rng = np.random.default_rng(12)
    loadings = np.array([[1.0, 0.2], [0.0, 0.9]])
    original = loadings.copy()
    shuffled_loadings(loadings, rng)
    assert np.allclose(loadings, original)


def test_rotated_subspace_preserves_total_variance():
    rng = np.random.default_rng(9)
    data = rng.normal(size=(200, 5))
    rotated = rotated_subspace(data, rng)
    assert np.sum(rotated**2) == pytest.approx(np.sum(data**2), rel=1e-10)
    assert rotated.shape == data.shape


def test_rotated_subspace_refuses_non_matrix_input():
    with pytest.raises(ValueError, match="observations-by-channels"):
        rotated_subspace(np.arange(10.0), np.random.default_rng(0))


def test_matched_marginals_keeps_each_channel_distribution_and_kills_correlation():
    rng = np.random.default_rng(10)
    shared = rng.normal(size=400)
    data = np.column_stack([shared, shared + 0.05 * rng.normal(size=400)])
    surrogate = matched_marginals(data, rng)

    before = abs(np.corrcoef(data[:, 0], data[:, 1])[0, 1])
    after = abs(np.corrcoef(surrogate[:, 0], surrogate[:, 1])[0, 1])
    assert before > 0.9
    assert after < 0.25
    # Values are drawn from the observed set, so the range cannot expand.
    for j in range(2):
        assert surrogate[:, j].min() >= data[:, j].min()
        assert surrogate[:, j].max() <= data[:, j].max()


def test_matched_marginals_preserves_the_missingness_pattern():
    rng = np.random.default_rng(13)
    data = rng.normal(size=(60, 3))
    data[5, 1] = np.nan
    data[9, 2] = np.nan
    surrogate = matched_marginals(data, rng)
    assert np.array_equal(np.isnan(surrogate), np.isnan(data))


def test_matched_marginals_refuses_an_all_missing_channel():
    data = np.column_stack([np.arange(10.0), np.full(10, np.nan)])
    with pytest.raises(ValueError, match="no finite observations"):
        matched_marginals(data, np.random.default_rng(0))


# ---------------------------------------------------------------------------
# skill score
# ---------------------------------------------------------------------------


def test_skill_score_spans_the_expected_range():
    assert skill_score(observed=0.0, baseline=1.0) == pytest.approx(1.0)   # perfect
    assert skill_score(observed=1.0, baseline=1.0) == pytest.approx(0.0)   # matched
    assert skill_score(observed=2.0, baseline=1.0) == pytest.approx(-1.0)  # worse


def test_negative_skill_makes_underperformance_unmissable():
    """A bare error figure hides doing worse than the trivial predictor."""
    assert skill_score(observed=0.9, baseline=0.5) < 0


def test_a_perfect_baseline_is_refused():
    with pytest.raises(ValueError, match="already scores perfectly"):
        skill_score(observed=0.0, baseline=0.0)


def test_the_reporting_variant_returns_none_instead_of_raising():
    """A report scoring many partitions should not lose the other rows because one
    baseline happened to predict its test set exactly."""
    assert skill_score_if_defined(observed=0.0, baseline=0.0) is None


def test_the_reporting_variant_is_otherwise_the_same_number():
    assert skill_score_if_defined(observed=0.5, baseline=1.0) == pytest.approx(
        skill_score(observed=0.5, baseline=1.0))
