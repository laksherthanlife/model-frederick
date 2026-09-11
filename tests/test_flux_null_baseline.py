"""Descriptive scoring is not evidence that PHB conditions are exchangeable units."""

from __future__ import annotations

import importlib.util
import itertools
import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.validation import ValidationFold, group_folds


@pytest.fixture(scope="module")
def runner():
    script = pathlib.Path(__file__).resolve().parents[1] / "scripts/flux_null_baseline.py"
    spec = importlib.util.spec_from_file_location("flux_null_baseline_regressions", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _reference_skill(observed, design, folds):
    """Independent normal-equation calculation, with one training mean per test fold."""
    prediction = np.full(len(observed), np.nan)
    baseline = np.full(len(observed), np.nan)
    for fold in folds:
        train = design[fold.train]
        beta = np.linalg.solve(train.T @ train, train.T @ observed[fold.train])
        prediction[fold.test] = design[fold.test] @ beta
        baseline[fold.test] = observed[fold.train].mean()
    return 1.0 - np.linalg.norm(observed - prediction) / np.linalg.norm(observed - baseline)


def test_intercept_only_matches_training_mean_baseline(runner):
    observed = np.array([-9.0, -5.0, -4.0, -3.0, 40.0])
    assert runner._skill(observed, np.ones((len(observed), 1))) == pytest.approx(0.0, abs=1e-14)


def test_skill_uses_training_only_comparator_for_every_holdout(runner):
    observed = np.array([0.0, 1.0, 1.5, 3.0, 12.0])
    design = np.column_stack([np.ones(len(observed)), np.arange(len(observed))])
    expected = _reference_skill(observed, design, group_folds(np.arange(len(observed))))
    assert runner._skill(observed, design) == pytest.approx(expected, abs=1e-12)


def test_grouped_scoring_reuses_project_folds_without_claiming_independence(runner):
    observed = np.array([-1.0, 0.5, 2.0, 3.0, 4.5, 12.0])
    design = np.column_stack([np.ones(len(observed)), np.arange(len(observed))])
    folds = group_folds(np.repeat(["a", "b", "c"], 2))
    assert runner._skill(observed, design, folds=folds) == pytest.approx(
        _reference_skill(observed, design, folds), abs=1e-12)


@pytest.mark.parametrize("groups", [np.arange(4), np.array([0, 0, 1, 1])])
def test_fixed_design_candidate_null_refits_each_explicit_permutation(runner, monkeypatch, groups):
    observed = np.array([0.0, 1.0, 3.0, 9.0])
    design = np.column_stack([np.ones(len(observed)), [0.0, 1.0, 3.0, 8.0]])
    permutations = np.array(list(itertools.permutations(range(len(observed)))))
    folds = group_folds(groups)
    expected = [_reference_skill(observed[order], design, folds) for order in permutations]
    calls = []
    original = runner._skill

    def record(target, matrix, *, folds=None):
        calls.append((target.copy(), matrix.copy(), folds))
        return original(target, matrix, folds=folds)

    monkeypatch.setattr(runner, "_skill", record)
    result = runner.candidate_null_distribution(observed, design, permutations, folds=folds)
    assert result == pytest.approx(expected, abs=1e-12)
    assert len(calls) == len(permutations)
    for call, order in zip(calls, permutations):
        np.testing.assert_array_equal(call[0], observed[order])
        np.testing.assert_array_equal(call[1], design)
        assert call[2] == folds
    np.testing.assert_array_equal(observed, [0.0, 1.0, 3.0, 9.0])
    np.testing.assert_array_equal(design[:, 1], [0.0, 1.0, 3.0, 8.0])


def test_fixed_design_null_depends_on_candidate_not_just_parameter_count(runner):
    observed = np.array([0.0, 1.0, 3.0, 9.0])
    first = np.column_stack([np.ones(len(observed)), [0.0, 1.0, 2.0, 3.0]])
    second = np.column_stack([np.ones(len(observed)), [0.0, 1.0, 2.0, 20.0]])
    permutations = np.array(list(itertools.permutations(range(len(observed)))))
    a = runner.candidate_null_distribution(observed, first, permutations)
    b = runner.candidate_null_distribution(observed, second, permutations)
    assert not np.allclose(a, b)
    assert np.percentile(a, 95) != pytest.approx(np.percentile(b, 95))


def test_candidate_null_refuses_to_invent_an_exchangeability_scheme(runner):
    with pytest.raises(ValueError, match="explicit.*permutations"):
        runner.candidate_null_distribution(np.arange(4.0), np.ones((4, 1)))


@pytest.mark.parametrize("permutations", [
    np.array([[0, 1, 1, 3]]),
    np.array([[0, 1, 2, 4]]),
    np.array([[0, 1, 2]]),
    np.array([[0.0, 1.0, 2.0, 3.0]]),
    np.empty((0, 4), dtype=int),
])
def test_candidate_null_requires_whole_response_permutations(runner, permutations):
    with pytest.raises(ValueError, match="permutation"):
        runner.candidate_null_distribution(np.arange(4.0), np.ones((4, 1)), permutations)


@pytest.mark.parametrize("design", [
    np.ones((4, 2)),
    np.column_stack([np.ones(4), [0.0, 0.0, 0.0, 1.0]]),
    np.eye(4),
])
def test_rank_deficient_training_design_is_refused_not_pseudoinverted(runner, design):
    with pytest.raises(ValueError, match="rank-deficient"):
        runner._skill(np.arange(4.0), design)


@pytest.mark.parametrize("constant", [1.0, 0.1])
def test_undefined_baseline_is_refused_even_with_mean_roundoff(runner, constant):
    with pytest.raises(ValueError, match="baseline.*zero"):
        runner._skill(np.full(4, constant), np.ones((4, 1)))


def test_held_out_feed_dummies_cannot_forecast_an_unseen_category(runner):
    states = runner.load()
    design = np.column_stack([
        np.ones(len(states)), pd.get_dummies(states.carbon_source, drop_first=True),
    ]).astype(float)
    folds = group_folds(states.carbon_source.to_numpy())
    with pytest.raises(ValueError, match="rank-deficient"):
        runner._skill(np.log(states.q_phb_mmol_per_gdcw_h), design, folds=folds)


def test_incomplete_fold_coverage_cannot_silently_shrink_score_denominator(runner):
    folds = group_folds(np.arange(4))[:-1]
    with pytest.raises(ValueError, match="exactly once"):
        runner._skill(np.arange(4.0), np.ones((4, 1)), folds=folds)


def test_fold_indices_must_address_supplied_observations(runner):
    fold = ValidationFold("outside", np.array([0, 1, 2]), np.array([4]), np.array([], dtype=int))
    with pytest.raises(ValueError, match="out of range"):
        runner._skill(np.arange(4.0), np.ones((4, 1)), folds=(fold,))


def test_failed_candidate_permutations_are_retained_as_undefined(runner):
    permutations = np.array(list(itertools.permutations(range(4))))
    result = runner.candidate_null_distribution(np.arange(4.0), np.ones((4, 2)), permutations)
    assert result.shape == (len(permutations),)
    assert np.isnan(result).all()


def test_random_design_demonstration_preserves_failed_draws(runner, monkeypatch):
    monkeypatch.setattr(runner, "SHUFFLES", 7)
    states = pd.DataFrame({"q_phb_mmol_per_gdcw_h": [1.0, 2.0, 3.0]})
    result = runner.null_distribution(states, 4, np.random.default_rng(4))
    assert result.shape == (7,)
    assert np.isnan(result).all()


def test_random_design_demonstration_is_seeded_and_has_zero_intercept_skill(runner, monkeypatch):
    monkeypatch.setattr(runner, "SHUFFLES", 7)
    states = runner.load()
    first = runner.null_distribution(states, 3, np.random.default_rng(4))
    second = runner.null_distribution(states, 3, np.random.default_rng(4))
    np.testing.assert_array_equal(first, second)
    assert runner.null_distribution(states, 1, np.random.default_rng(4)) == pytest.approx(
        np.zeros(7), abs=1e-14)


@pytest.mark.parametrize("bad_flux", [0.0, -1.0, np.nan, np.inf])
def test_flux_must_be_finite_and_positive_before_log_scoring(runner, bad_flux):
    states = pd.DataFrame({"q_phb_mmol_per_gdcw_h": [1.0, 2.0, bad_flux]})
    with pytest.raises(ValueError, match="finite and positive"):
        runner.null_distribution(states, 2, np.random.default_rng(4))


def test_real_data_metadata_is_conditions_not_biological_replicate_identifiers(runner):
    states = runner.load()
    assert states.carbon_source.value_counts().to_dict() == {
        "glucose": 4, "glucose_ethanol_1to2": 4, "ethanol": 3,
    }
    assert states.dilution_rate_per_h.nunique() == 4
    assert states.state_id.nunique() == len(states) == 11
    assert not {"biological_replicate_id", "cultivation_id", "randomization_block"}.intersection(
        states.columns)


def test_main_separates_demo_from_candidates_and_withholds_inference(runner, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(runner, "SHUFFLES", 9)
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))

    def refuse_unspecified_null(*args, **kwargs):
        pytest.fail("real PHB data has no verified response-permutation scheme")

    monkeypatch.setattr(runner, "candidate_null_distribution", refuse_unspecified_null)
    assert runner.main() == 0
    demo = pd.read_csv(tmp_path / "flux_null_baseline.csv")
    candidates = pd.read_csv(tmp_path / "flux_candidate_scores.csv")
    assert demo.parameters.tolist() == [1, 2, 3, 4, 5]
    assert demo.analysis_kind.eq("random_design_demonstration_not_inference").all()
    assert demo.shuffles.eq(9).all()
    assert demo.seed.eq(runner.SEED).all()
    assert demo.n_failed.eq(0).all()
    assert candidates.model.tolist() == [
        "constant only", "log(mu)", "carbon source", "carbon source + log(mu)",
    ]
    assert candidates.inference_status.eq("pending_exchangeability").all()
    assert candidates.inference_reason.str.contains("cultivation").all()
    assert candidates.score_status.eq("ok").all()
    assert candidates.validation.eq("leave_one_state_out_descriptive").all()
    assert candidates.baseline.eq("fold_training_mean_log_flux").all()
    assert candidates.n_states.eq(11).all()
    assert candidates.n_carbon_sources.eq(3).all()
    assert not {"null_skill_p95", "permutation_p", "verdict"}.intersection(candidates.columns)
    assert candidates.loc[candidates.model == "constant only", "loo_skill"].iloc[0] == pytest.approx(
        0.0, abs=1e-14)
    output = capsys.readouterr().out
    assert "not a candidate-specific threshold" in output
    assert "PENDING" in output
    assert "beats the null" not in output
    assert "INSIDE THE NULL" not in output


def test_main_retains_failed_candidates_and_withholds_partial_demo_quantiles(
        runner, monkeypatch, tmp_path):
    states = runner.load().iloc[:4].copy()
    states["carbon_source"] = ["a", "a", "a", "b"]
    monkeypatch.setattr(runner, "load", lambda: states)
    monkeypatch.setattr(runner, "SHUFFLES", 5)
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    assert runner.main() == 0
    candidates = pd.read_csv(tmp_path / "flux_candidate_scores.csv")
    failed = candidates[candidates.model.str.startswith("carbon source")]
    assert len(failed) == 2
    assert failed.loo_skill.isna().all()
    assert failed.score_status.eq("failed").all()
    assert failed.failure_reason.str.contains("rank-deficient").all()
    demo = pd.read_csv(tmp_path / "flux_null_baseline.csv")
    failed = demo[demo.parameters >= 4]
    assert failed.n_failed.eq(5).all()
    assert failed.null_skill_p95.isna().all()
    assert failed.score_status.eq("failed").all()


def test_a_single_failed_draw_withholds_the_whole_demo_percentile(runner, monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "SHUFFLES", 3)
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    monkeypatch.setattr(runner, "null_distribution", lambda *args: np.array([0.1, np.nan, 0.3]))
    assert runner.main() == 0
    demo = pd.read_csv(tmp_path / "flux_null_baseline.csv")
    assert demo.n_failed.eq(1).all()
    assert demo.null_skill_p95.isna().all()
    assert demo.score_status.eq("failed").all()
