from __future__ import annotations

import ast
import builtins
import copy
import inspect
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.special import expit

from ystwin.analysis.in_silico import InSilicoStudent, fit_student, score_student
from ystwin.reporter import ReporterKinetics, simulate_reporter


LATENT_NAMES = ("upr", "oxidative", "burden")
CONTROL_NAMES = (
    "growth_retention",
    "enzyme_budget_scale",
    "ngam_mmol_per_gdcw_h",
    "allocation_fraction",
)


def _fixture_rhs(state, imposed):
    upr, oxidative, burden = state
    first, second = imposed
    activation = np.array([
        0.035 + 0.85 * first + 0.07 * second + 0.10 * burden,
        0.045 + 0.08 * first + 0.80 * second + 0.08 * upr,
        0.025 + 0.30 * first + 0.35 * second + 0.12 * upr + 0.10 * oxidative,
    ])
    relaxation = np.array([
        0.38 + 0.12 * second,
        0.42 + 0.08 * first,
        0.32 + 0.04 * first + 0.06 * second,
    ])
    return (1.0 - state) * activation - state * relaxation


def _fixture_controls(state):
    upr, oxidative, burden = state.T
    return np.column_stack([
        expit(1.6 - 1.0 * upr - 0.8 * oxidative - 1.1 * burden),
        expit(1.4 - 0.6 * upr - 0.9 * oxidative - 0.8 * burden),
        np.exp(-0.1 + 0.5 * upr + 0.4 * oxidative + 0.6 * burden),
        expit(-2.0 + 0.4 * upr + 0.3 * oxidative + 0.7 * burden),
    ])


def _episodes(n, seed, noise_cv=0.001, degradation=0.05):
    rng = np.random.default_rng(seed)
    times = np.linspace(0.0, 6.0, 61)
    knots = np.linspace(0.0, 6.0, 7)
    loading = np.array([
        [1.30, 0.12, 0.18, 0.65],
        [0.14, 1.20, 0.15, 0.55],
        [0.17, 0.16, 1.25, 0.70],
    ])
    result = []
    for index in range(n):
        doses = rng.uniform(0.02, 1.0, (len(knots), 2))
        imposed = np.column_stack([
            np.interp(times, knots, doses[:, channel]) for channel in range(2)
        ])
        initial = rng.uniform(0.06, 0.88, 3)

        def rhs(time, state):
            known = np.array([
                np.interp(time, times, imposed[:, channel]) for channel in range(2)
            ])
            return _fixture_rhs(state, known)

        solution = solve_ivp(
            rhs, (times[0], times[-1]), initial, t_eval=times,
            rtol=2e-9, atol=1e-11, max_step=0.1,
        )
        assert solution.success
        latent = solution.y.T
        growth = 0.34 * np.exp(-latent @ np.array([0.30, 0.25, 0.40]))
        density = rng.uniform(0.06, 0.20) * np.exp(
            cumulative_trapezoid(growth, times, initial=0.0)
        )
        activity = np.array([0.16, 0.19, 0.14, 0.18]) + latent @ loading
        reporters = np.column_stack([
            simulate_reporter(
                times, activity[:, channel], growth,
                kinetics=ReporterKinetics(k_deg=degradation),
            )
            for channel in range(4)
        ])
        observations = np.column_stack([density, reporters])
        observations *= np.exp(rng.normal(0.0, noise_cv, observations.shape))
        result.append(SimpleNamespace(
            episode_id=f"fixture-{seed}-{index}",
            times_h=times.copy(),
            inputs=imposed,
            observations=observations,
            latent=latent,
            latent_derivative=np.array([
                _fixture_rhs(state, known) for state, known in zip(latent, imposed)
            ]),
            controls=_fixture_controls(latent),
            growth_rate_per_h=growth,
            metadata={
                "family": "fixture_bilinear",
                "protocol": "independent_interpolated_inputs",
                "reporter_degradation_per_h": degradation,
                "reporter_maturation": "instantaneous",
                "latent_names": list(LATENT_NAMES),
                "control_names": list(CONTROL_NAMES),
            },
        ))
    return result


@pytest.fixture(scope="module")
def training():
    return _episodes(24, seed=120)


@pytest.fixture(scope="module")
def holdout():
    return _episodes(6, seed=911)


@pytest.fixture(scope="module")
def student(training):
    return fit_student(training, seed=7)


def _rmse(actual, expected):
    return float(np.sqrt(np.mean((np.asarray(actual) - np.asarray(expected)) ** 2)))


def _corrupt(episodes, field):
    result = copy.deepcopy(episodes)
    combined = np.concatenate([getattr(episode, field) for episode in result])
    permutation = np.random.default_rng(361).permutation(len(combined))
    shuffled = combined[permutation]
    start = 0
    for episode in result:
        stop = start + len(episode.times_h)
        setattr(episode, field, shuffled[start:stop].copy())
        start = stop
    return result


def test_identifiable_whole_episode_holdouts_are_recovered(student, holdout):
    assert isinstance(student, InSilicoStudent)
    scores = score_student(student, holdout, prefix_h=2.0)
    assert isinstance(scores, pd.DataFrame)
    assert set(scores.episode_id) == {episode.episode_id for episode in holdout}
    assert scores.state_rmse.mean() < 0.065
    assert scores.forecast_rmse.mean() < 0.075
    assert scores.control_normalized_rmse.mean() < 0.55
    assert scores.forecast_rmse.mean() < scores.persistence_rmse.mean() * 0.65
    assert (scores.prefix_h == 2.0).all()
    assert (scores.prefix_end_h == 2.0).all()
    assert (scores.forecast_horizon_h == 4.0).all()
    assert (scores.n_forecast_observations == 40).all()
    for name in LATENT_NAMES:
        assert f"state_rmse_{name}" in scores
        assert f"forecast_rmse_{name}" in scores
    for name in CONTROL_NAMES:
        assert f"control_rmse_{name}" in scores


def test_shuffled_latent_labels_destroy_the_encoder(training, student, holdout):
    corrupt = fit_student(_corrupt(training, "latent"), seed=7)
    good = np.mean([
        _rmse(student.infer(episode.times_h, episode.observations), episode.latent)
        for episode in holdout
    ])
    bad = np.mean([
        _rmse(corrupt.infer(episode.times_h, episode.observations), episode.latent)
        for episode in holdout
    ])
    assert bad > good * 2.0
    assert bad - good > 0.04


def test_shuffled_derivative_labels_destroy_open_loop_forecasts(training, student, holdout):
    corrupt = fit_student(_corrupt(training, "latent_derivative"), seed=7)
    good = score_student(student, holdout, prefix_h=2.0).forecast_rmse.mean()
    bad = score_student(corrupt, holdout, prefix_h=2.0).forecast_rmse.mean()
    assert bad > good * 2.0 + 0.02


def test_shuffled_control_labels_destroy_the_control_map(training, student, holdout):
    corrupt = fit_student(_corrupt(training, "controls"), seed=7)
    good = np.mean([
        _rmse(student.predict_controls(episode.latent), episode.controls) for episode in holdout
    ])
    bad = np.mean([
        _rmse(corrupt.predict_controls(episode.latent), episode.controls) for episode in holdout
    ])
    assert good < 0.015
    assert bad > good * 5.0 + 0.03


def test_future_observations_and_time_grid_cannot_change_past_inference(student, holdout):
    episode = holdout[0]
    cutoff = 19
    expected = student.infer(episode.times_h[:cutoff], episode.observations[:cutoff])
    mutated = episode.observations.copy()
    mutated[cutoff:, 0] *= 100.0
    mutated[cutoff:, 1:] *= 0.05
    extended_times = episode.times_h.copy()
    extended_times[cutoff:] += np.arange(len(extended_times) - cutoff) * 0.03
    actual = student.infer(extended_times, mutated)
    np.testing.assert_array_equal(actual[:cutoff], expected)
    for stop in (1, 2, 4, 8, 12):
        np.testing.assert_allclose(
            student.infer(episode.times_h[:stop], episode.observations[:stop]),
            expected[:stop], rtol=0, atol=1e-13,
        )


def test_forecast_has_only_prefix_observations_and_no_episode_argument(student, holdout):
    assert list(inspect.signature(student.infer).parameters) == ["times_h", "observations"]
    assert list(inspect.signature(student.forecast).parameters) == [
        "times_h", "inputs", "prefix_times_h", "prefix_observations",
    ]
    episode = copy.deepcopy(holdout[0])
    prefix = episode.times_h <= 2.0
    prefix_times = episode.times_h[prefix].copy()
    prefix_observations = episode.observations[prefix].copy()
    expected = student.forecast(episode.times_h, episode.inputs, prefix_times, prefix_observations)
    episode.latent[:] = np.nan
    episode.latent_derivative[:] = np.nan
    episode.controls[:] = np.nan
    episode.growth_rate_per_h[:] = np.nan
    episode.observations[~prefix] = np.nan
    actual = student.forecast(episode.times_h, episode.inputs, prefix_times, prefix_observations)
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_allclose(
        actual[prefix], student.infer(prefix_times, prefix_observations), rtol=0, atol=1e-13,
    )
    assert np.isfinite(actual).all()
    assert (actual[~prefix] >= 0).all() and (actual[~prefix] <= 1).all()


def test_known_future_inputs_affect_forecast_but_not_prefix(student, holdout):
    episode = holdout[0]
    prefix = episode.times_h <= 2.0
    changed = episode.inputs.copy()
    changed[~prefix] = 0.0
    first = student.forecast(
        episode.times_h, episode.inputs, episode.times_h[prefix], episode.observations[prefix],
    )
    second = student.forecast(
        episode.times_h, changed, episode.times_h[prefix], episode.observations[prefix],
    )
    np.testing.assert_array_equal(first[prefix], second[prefix])
    assert _rmse(first[~prefix], second[~prefix]) > 0.10


def test_forecast_prefix_queries_are_causal_not_interpolated(student, holdout):
    episode = holdout[0]
    prefix_times = episode.times_h[:15]
    prefix_observations = episode.observations[:15]
    query_times = prefix_times[:-1] + 0.025
    known_inputs = episode.inputs[:len(query_times)]
    actual = student.forecast(query_times, known_inputs, prefix_times, prefix_observations)
    expected = student.infer(prefix_times, prefix_observations)[:-1]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-13)


def test_score_passes_only_the_declared_prefix(student, holdout, monkeypatch):
    original = InSilicoStudent.forecast
    calls = []

    def checked(self, times_h, inputs, prefix_times_h, prefix_observations):
        assert prefix_times_h[-1] <= times_h[0] + 1.75
        assert len(prefix_observations) == len(prefix_times_h)
        assert len(prefix_observations) < len(times_h)
        calls.append(len(prefix_observations))
        return original(self, times_h, inputs, prefix_times_h, prefix_observations)

    monkeypatch.setattr(InSilicoStudent, "forecast", checked)
    scores = score_student(student, holdout[:2], prefix_h=1.75)
    assert calls == [18, 18]
    np.testing.assert_allclose(scores.prefix_end_h, 1.7)
    np.testing.assert_allclose(scores.forecast_horizon_h, 4.3)


def test_duplicate_and_overlapping_episode_ids_are_rejected(training, student, holdout):
    with pytest.raises(ValueError, match="duplicate|unique"):
        fit_student([training[0], training[0]])
    with pytest.raises(ValueError, match="training|overlap|held.out"):
        score_student(student, [training[0]])
    with pytest.raises(ValueError, match="duplicate|unique"):
        score_student(student, [holdout[0], holdout[0]])


def test_mutating_any_episode_after_fit_cannot_change_weights_or_predictions(training, holdout):
    mutable_training = copy.deepcopy(training)
    mutable_holdout = copy.deepcopy(holdout)
    student = fit_student(mutable_training, seed=17)
    metadata = json.dumps(student.metadata, sort_keys=True, allow_nan=False)
    query_times = mutable_holdout[0].times_h.copy()
    query_observations = mutable_holdout[0].observations.copy()
    expected = student.infer(query_times, query_observations)
    for episode in [*mutable_training, *mutable_holdout]:
        for field in ("observations", "latent", "inputs", "latent_derivative", "controls"):
            getattr(episode, field)[:] = 0.123
        episode.metadata["reporter_degradation_per_h"] = 100.0
    np.testing.assert_array_equal(student.infer(query_times, query_observations), expected)
    assert json.dumps(student.metadata, sort_keys=True, allow_nan=False) == metadata


def test_deterministic_fit_and_train_only_scaling(training, student, holdout):
    duplicate = fit_student(training, seed=7)
    assert json.dumps(duplicate.metadata, sort_keys=True) == json.dumps(student.metadata, sort_keys=True)
    episode = holdout[0]
    before = json.dumps(student.metadata, sort_keys=True)
    scores = score_student(student, holdout)
    assert np.isfinite(scores.state_rmse).all()
    assert json.dumps(student.metadata, sort_keys=True) == before
    np.testing.assert_array_equal(
        student.infer(episode.times_h, episode.observations),
        duplicate.infer(episode.times_h, episode.observations),
    )


def test_metadata_names_the_learned_tasks_and_contains_real_matrices(student, training):
    metadata = student.metadata
    json.dumps(metadata, allow_nan=False)
    assert metadata["training_episode_ids"] == [episode.episode_id for episode in training]
    assert metadata["uses_product_targets"] is False
    assert metadata["targets"]["encoder"] == list(LATENT_NAMES)
    assert metadata["targets"]["controller"] == list(CONTROL_NAMES)
    assert set(metadata["learned_models"]) >= {"encoder", "dynamics", "controller"}
    for name in ("encoder", "dynamics", "controller"):
        weights = np.asarray(metadata["learned_models"][name]["coefficients"])
        assert weights.ndim == 2 and weights.size > 3
        assert np.isfinite(weights).all()
        assert np.any(np.abs(weights) > 1e-6)
    assert metadata["input_schema"]
    assert metadata["assumptions"]
    assert metadata["hyperparameters"]["ridge"] == 1e-6
    assert metadata["training_errors"]


def test_control_links_preserve_physical_domains_and_batch_shape(student):
    latent = np.random.default_rng(451).uniform(0.0, 1.0, (2, 9, 3))
    controls = student.predict_controls(latent)
    assert controls.shape == (2, 9, 4)
    assert np.isfinite(controls).all()
    assert (controls > 0).all()
    assert (controls[..., [0, 1, 3]] < 1).all()
    assert student.predict_controls(np.array([0.3, 0.2, 0.1])).shape == (4,)


@pytest.mark.parametrize("times,observations", [
    ([], np.ones((0, 5))),
    ([0, 0], np.ones((2, 5))),
    ([1, 0], np.ones((2, 5))),
    ([0, np.nan], np.ones((2, 5))),
    ([[0, 1]], np.ones((2, 5))),
    ([0, 1], np.ones((2, 4))),
    ([0, 1], np.ones((3, 5))),
    ([0, 1], np.full((2, 5), np.nan)),
    ([0, 1], np.zeros((2, 5))),
    ([0, 1], -np.ones((2, 5))),
    (np.array([0, 1j]), np.ones((2, 5))),
])
def test_infer_validates_inputs(student, times, observations):
    with pytest.raises(ValueError):
        student.infer(times, observations)
    with pytest.raises(ValueError):
        student.forecast([0.0, 1.0], np.zeros((2, 2)), times, observations)


@pytest.mark.parametrize("ridge", [-1.0, 0.0, np.nan, np.inf, "invalid", None, 1j])
def test_fit_rejects_invalid_regularization(training, ridge):
    with pytest.raises(ValueError, match="ridge"):
        fit_student(training, ridge=ridge)


def test_fit_validates_episode_schema(training):
    with pytest.raises(ValueError, match="episode|training|empty"):
        fit_student([])
    for field, value in [
        ("latent", np.ones((61, 2))),
        ("latent", np.full((61, 3), 1.2)),
        ("inputs", np.full((61, 2), -1.0)),
        ("latent_derivative", np.full((61, 3), np.inf)),
        ("controls", np.zeros((61, 4))),
        ("episode_id", ""),
    ]:
        bad = copy.deepcopy(training[0])
        setattr(bad, field, value)
        with pytest.raises(ValueError):
            fit_student([bad])


def test_declared_reporter_kinetics_are_checked(training):
    different = copy.deepcopy(training[:2])
    different[1].metadata["reporter_degradation_per_h"] = 0.2
    with pytest.raises(ValueError, match="degradation|kinetics"):
        fit_student(different)
    bad = copy.deepcopy(training[0])
    bad.metadata["reporter_maturation"] = "finite"
    with pytest.raises(ValueError, match="maturation"):
        fit_student([bad])


def test_forecast_and_controls_validate_shapes_and_domains(student, holdout):
    episode = holdout[0]
    with pytest.raises(ValueError):
        student.predict_controls(np.zeros((4, 4)))
    with pytest.raises(ValueError):
        student.predict_controls([np.nan, 0.0, 0.0])
    for inputs in (episode.inputs[:-1], np.full_like(episode.inputs, -0.1)):
        with pytest.raises(ValueError):
            student.forecast(
                episode.times_h, inputs, episode.times_h[:12], episode.observations[:12],
            )
    with pytest.raises(ValueError, match="prefix|observation"):
        student.forecast(
            episode.times_h, episode.inputs, episode.times_h[1:12], episode.observations[1:12],
        )
    with pytest.raises(ValueError, match="prefix|horizon"):
        score_student(student, holdout, prefix_h=10.0)
    with pytest.raises(ValueError, match="prefix"):
        score_student(student, holdout, prefix_h=-0.1)


def test_no_teacher_or_product_dependencies_or_table_reads(training, holdout, monkeypatch):
    import ystwin.analysis.in_silico as module

    tree = ast.parse(inspect.getsource(module))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert not any(part in ("generator", "reporter", "fba", "bridge", "paths")
                           for part in (node.module or "").split("."))
        if isinstance(node, ast.Import):
            assert not any(alias.name.startswith("ystwin.") for alias in node.names)

    def forbidden(*args, **kwargs):
        raise AssertionError("student attempted an external table/file or reporter/teacher call")

    monkeypatch.setattr(builtins, "open", forbidden)
    monkeypatch.setattr(Path, "open", forbidden)
    monkeypatch.setattr(pd, "read_csv", forbidden)
    monkeypatch.setattr(pd, "read_excel", forbidden)
    monkeypatch.setattr(np, "load", forbidden)
    monkeypatch.setattr(np, "loadtxt", forbidden)
    monkeypatch.setattr("ystwin.reporter.simulate_reporter", forbidden)
    isolated = fit_student(training, seed=19)
    episode = holdout[0]
    isolated.infer(episode.times_h, episode.observations)
    isolated.predict_controls(episode.latent)
    isolated.forecast(
        episode.times_h, episode.inputs, episode.times_h[:21], episode.observations[:21],
    )


@pytest.fixture(scope="module")
def teacher_training():
    from ystwin.generator.in_silico import generate_episodes

    return generate_episodes(24, seed=412)


@pytest.fixture(scope="module")
def teacher_holdout():
    from ystwin.generator.in_silico import generate_episodes

    return generate_episodes(6, seed=813)


@pytest.fixture(scope="module")
def teacher_student(teacher_training):
    return fit_student(teacher_training, seed=11)


def test_actual_teacher_bilinear_system_is_in_the_learned_model_class(teacher_student, teacher_holdout):
    error = np.mean([
        _rmse(teacher_student.predict_derivative(episode.latent, episode.inputs), episode.latent_derivative)
        for episode in teacher_holdout
    ])
    assert error < 0.001
    scores = score_student(teacher_student, teacher_holdout, prefix_h=2.0)
    assert scores.state_rmse.mean() < 0.065
    assert scores.forecast_rmse.mean() < 0.05
    assert scores.forecast_rmse.mean() < scores.persistence_rmse.mean() * 0.7


def test_retention_control_learning_handles_the_healthy_boundary(teacher_student, teacher_holdout):
    error = teacher_student.metadata["training_errors"]["controller_teacher_forced_rmse_per_coordinate"]
    assert max(error[:2]) < 0.002
    assert _rmse(
        teacher_student.predict_controls(teacher_holdout[0].latent), teacher_holdout[0].controls,
    ) < 0.002
    healthy = teacher_student.predict_controls(np.zeros(3))
    np.testing.assert_allclose(healthy[:2], 1.0, atol=1e-4, rtol=0)


def test_controller_links_are_physical_even_for_reported_encoder_excursions(student):
    prediction = student.predict_controls([[-0.1, -0.2, -0.1], [1.1, 1.2, 1.1]])
    assert np.isfinite(prediction).all()
    assert (prediction > 0).all()
    assert (prediction[:, [0, 1, 3]] <= 1).all()


def test_actual_teacher_functions_are_forbidden_after_dataset_generation(
    teacher_training, teacher_holdout, monkeypatch,
):
    import ystwin.generator.in_silico as teacher

    def forbidden(*args, **kwargs):
        raise AssertionError("student called the actual teacher")

    for name, value in vars(teacher).copy().items():
        if inspect.isfunction(value):
            monkeypatch.setattr(teacher, name, forbidden)
    monkeypatch.setattr("ystwin.reporter.simulate_reporter", forbidden)
    student = fit_student(teacher_training)
    episode = teacher_holdout[0]
    student.infer(episode.times_h, episode.observations)
    student.predict_controls(episode.latent)
    result = student.forecast(
        episode.times_h, episode.inputs, episode.times_h[:21], episode.observations[:21],
    )
    assert np.isfinite(result).all()


def test_invariant_dynamics_point_inward_at_every_cube_face(student):
    rng = np.random.default_rng(42)
    weights = np.asarray(student.metadata["learned_models"]["dynamics"]["coefficients"])
    assert (weights >= 0).all()
    for coordinate in range(3):
        states = rng.uniform(0, 1, (100, 3))
        known = rng.uniform(0, 3, (100, 2))
        states[:, coordinate] = 0.0
        assert (student.predict_derivative(states, known)[:, coordinate] >= 0).all()
        states[:, coordinate] = 1.0
        assert (student.predict_derivative(states, known)[:, coordinate] <= 0).all()


def test_bad_observations_are_reported_not_silently_clipped(student, holdout):
    episode = holdout[0]
    prefix_times = episode.times_h[:21]
    prefix = episode.observations[:21].copy()
    prefix[-1, 1:] *= 1000
    inferred = student.infer(prefix_times, prefix)
    assert np.any((inferred[-1] < 0) | (inferred[-1] > 1))
    with pytest.raises(ValueError, match="prefix.*outside.*domain.*violation"):
        student.forecast(episode.times_h, episode.inputs, prefix_times, prefix)


def test_nested_teacher_kinetics_are_authoritative_and_other_truth_metadata_is_ignored(training):
    original = copy.deepcopy(training[:3])
    for episode in original:
        episode.metadata.pop("reporter_degradation_per_h")
        episode.metadata["parameters"] = {
            "reporter_kinetics": {"k_deg": 0.05, "k_mat": None, "k_deg_immature": None},
            "reporter_loadings": "must not be used",
            "growth_cost_coefficients": "must not be used",
        }
    first = fit_student(original)
    changed = copy.deepcopy(original)
    changed[1].metadata["parameters"]["reporter_kinetics"]["k_deg"] = 0.2
    with pytest.raises(ValueError, match="degradation|kinetics"):
        fit_student(changed)
    changed = copy.deepcopy(original)
    changed[0].metadata["parameters"]["reporter_kinetics"]["k_mat"] = 2.0
    with pytest.raises(ValueError, match="maturation"):
        fit_student(changed)
    changed = copy.deepcopy(original)
    for episode in changed:
        episode.metadata["parameters"]["reporter_loadings"] = 1e50
        episode.metadata["parameters"]["growth_cost_coefficients"] = -1e50
        episode.growth_rate_per_h[:] = np.nan
    second = fit_student(changed)
    assert json.dumps(first.metadata, sort_keys=True) == json.dumps(second.metadata, sort_keys=True)


def test_saturating_family_is_a_real_untrained_challenge(teacher_student, teacher_holdout):
    from ystwin.generator.in_silico import generate_episodes

    challenge = generate_episodes(6, seed=814, family="saturating")
    reference = score_student(teacher_student, teacher_holdout).forecast_rmse.mean()
    shifted = score_student(teacher_student, challenge).forecast_rmse.mean()
    assert shifted > reference + 0.015


@pytest.fixture(scope="module")
def reference_training():
    from ystwin.generator.in_silico import generate_episodes

    return generate_episodes(
        64, seed=0, family="bilinear", protocol="random", hours=8.0, dt=0.1, noise_cv=0.001,
    )


@pytest.fixture(scope="module")
def reference_student(reference_training):
    return fit_student(reference_training, seed=0)


@pytest.fixture(scope="module")
def nuisance_validation():
    from ystwin.generator.in_silico import generate_episodes, observe_reporters

    trajectories = tuple(
        episode
        for seed in (52001, 52002)
        for episode in generate_episodes(
            6, seed=seed, family="bilinear", protocol="random", hours=8.0, dt=0.1, noise_cv=0.0,
        )
    )
    result = {}
    for episode in trajectories:
        times = episode.times_h
        phase = (times - times[0]) / (times[-1] - times[0])
        ramp = 0.075 + (0.42 - 0.075) * (1.0 - np.cos(np.pi * phase)) / 2.0
        growth_series = {
            "assay": episode.growth_rate_per_h,
            **{f"constant_{value}": np.full(len(times), value) for value in (0.075, 0.12, 0.24, 0.42)},
            "ramp_up": ramp,
            "ramp_down": ramp[::-1],
        }
        for growth_name, growth in growth_series.items():
            for initial_density in (0.03, 0.1):
                density = initial_density * np.exp(cumulative_trapezoid(growth, times, initial=0.0))
                for noise_cv in (0.0, 0.001):
                    observed = observe_reporters(
                        times, episode.latent, growth, density,
                        seed=episode.metadata["seed"], noise_cv=noise_cv,
                    )
                    key = (noise_cv, growth_name, initial_density)
                    identifier = f"{episode.episode_id}/{growth_name}/{initial_density}/{noise_cv}"
                    result.setdefault(key, []).append(SimpleNamespace(
                        episode_id=identifier,
                        trajectory_id=episode.episode_id,
                        times_h=times,
                        inputs=episode.inputs,
                        observations=observed,
                        latent=episode.latent,
                        latent_derivative=episode.latent_derivative,
                        controls=episode.controls,
                        growth_rate_per_h=growth,
                        metadata={
                            **episode.metadata,
                            "episode_id": identifier,
                            "initial_cell_density": initial_density,
                            "growth_source": "supplied_same_culture",
                            "density_source": "supplied_same_culture",
                            "nuisance_growth": growth_name,
                            "noise": {
                                **episode.metadata["noise"], "cv": noise_cv,
                                "log_variance": float(np.log1p(noise_cv ** 2)),
                            },
                        },
                    ))
    return result


def test_main_encoder_is_refitted_from_only_four_activity_balances(reference_student, reference_training):
    from ystwin.analysis.in_silico import _causal_features

    names = [f"causal_activity_{index}" for index in range(1, 5)]
    encoder = reference_student.metadata["learned_models"]["encoder"]
    assert encoder["feature_names"] == names
    assert reference_student.encoder.coefficients.shape == (4, 3)
    features = np.concatenate([
        _causal_features(episode.times_h, episode.observations, reference_student.reporter_degradation_per_h)[1:]
        for episode in reference_training
    ])
    assert features.shape[1] == 4
    targets = np.concatenate([episode.latent[1:] for episode in reference_training])
    design = (features - features.mean(axis=0)) / features.std(axis=0)
    expected = np.linalg.lstsq(
        np.vstack([design, np.sqrt(len(features) * 1e-6) * np.eye(4)]),
        np.vstack([targets - targets.mean(axis=0), np.zeros((4, 3))]),
        rcond=None,
    )[0]
    np.testing.assert_allclose(reference_student.encoder.coefficients, expected, rtol=0, atol=1e-13)
    residual = reference_student.encoder.predict(features) - targets
    covariance = residual.T @ residual / len(residual) + 1e-6 * np.eye(3)
    np.testing.assert_array_equal(reference_student.encoder_residual_covariance, covariance)
    np.testing.assert_array_equal(
        reference_student.metadata["initialization"]["encoder_residual_covariance"], covariance,
    )
    assert reference_student.metadata["initialization"]["reference_tail_probability"] == 1e-6


@pytest.mark.parametrize("noise_cv", [0.0, 0.001])
def test_same_latent_trajectories_have_no_growth_shortcut_bias(
    reference_student, nuisance_validation, noise_cv,
):
    reference = nuisance_validation[(noise_cv, "assay", 0.1)]
    baseline = np.array([
        reference_student.infer(episode.times_h, episode.observations) for episode in reference
    ])
    for (cv, growth_name, density), episodes in nuisance_validation.items():
        if cv != noise_cv:
            continue
        for original, episode in zip(reference, episodes):
            assert original.trajectory_id == episode.trajectory_id
            np.testing.assert_array_equal(original.latent, episode.latent)
            np.testing.assert_array_equal(original.inputs, episode.inputs)
            np.testing.assert_array_equal(original.controls, episode.controls)
        inferred = np.array([
            reference_student.infer(episode.times_h, episode.observations) for episode in episodes
        ])
        truth = np.array([episode.latent for episode in episodes])
        controls = reference_student.predict_controls(inferred)
        control_truth = np.array([episode.controls for episode in episodes])
        assert _rmse(inferred, truth) < 0.06, (noise_cv, growth_name, density)
        assert _rmse(controls, control_truth) < 0.02, (noise_cv, growth_name, density)
        systematic_shift = np.mean(inferred[:, 8:] - baseline[:, 8:], axis=(0, 1))
        assert np.max(np.abs(systematic_shift)) < 0.01, (noise_cv, growth_name, density, systematic_shift)
        if density == 0.03:
            other = nuisance_validation[(noise_cv, growth_name, 0.1)]
            other_inferred = np.array([
                reference_student.infer(episode.times_h, episode.observations) for episode in other
            ])
            np.testing.assert_allclose(inferred, other_inferred, rtol=0, atol=1e-12)


def test_nuisance_validation_has_no_future_or_teacher_dependency(
    reference_student, nuisance_validation, monkeypatch,
):
    import ystwin.generator.in_silico as teacher

    def forbidden(*args, **kwargs):
        raise AssertionError("prediction called a teacher or reporter function")

    for name, value in vars(teacher).copy().items():
        if inspect.isfunction(value):
            monkeypatch.setattr(teacher, name, forbidden)
    monkeypatch.setattr("ystwin.reporter.simulate_reporter", forbidden)
    before = json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False)
    for (_, _, _), episodes in nuisance_validation.items():
        episode = episodes[0]
        times = episode.times_h
        for stop in (1, 2, 5, 9, 21):
            expected = reference_student.infer(times[:stop], episode.observations[:stop])
            changed = episode.observations.copy()
            changed[stop:, 0] *= 10.0
            changed[stop:, 1:] *= 0.2
            extended_times = times.copy()
            extended_times[stop:] += np.arange(len(times) - stop) * 0.03
            actual = reference_student.infer(extended_times, changed)
            np.testing.assert_array_equal(actual[:stop], expected)
        predicted = reference_student.forecast(
            times, episode.inputs, times[:21], episode.observations[:21],
        )
        assert np.isfinite(predicted).all()
        assert np.isfinite(reference_student.predict_controls(predicted)).all()
    assert json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False) == before


def test_forecast_resolves_each_known_input_knot_including_short_pulses(student, monkeypatch):
    import ystwin.analysis.in_silico as module

    times = np.array([0.0, 0.11321, 0.11322, 0.11323, 0.6])
    inputs = np.zeros((len(times), 2))
    inputs[2] = [1.0, 0.7]
    initial = np.array([0.2, 0.15, 0.1])
    prefix_times = np.array([0.0, 0.05])
    monkeypatch.setattr(
        InSilicoStudent, "infer", lambda self, times_h, observations: np.tile(initial, (len(times_h), 1)),
    )
    expected = [initial]
    start = prefix_times[-1]
    for stop in times[1:]:
        solution = solve_ivp(
            lambda time, state: student.dynamics.predict(
                state, np.array([np.interp(time, times, inputs[:, j]) for j in range(2)]),
            ),
            (start, stop), expected[-1], t_eval=[stop], method="DOP853", rtol=1e-12, atol=1e-14,
        )
        assert solution.success
        expected.append(solution.y[:, -1])
        start = stop
    intervals = []
    original = module.solve_ivp

    def recording(rhs, interval, *args, **kwargs):
        intervals.append(interval)
        return original(rhs, interval, *args, **kwargs)

    monkeypatch.setattr(module, "solve_ivp", recording)
    actual = student.forecast(times, inputs, prefix_times, np.ones((2, 5)))
    assert intervals == list(zip(np.r_[prefix_times[-1], times[1:-1]], times[1:]))
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-9)


@pytest.mark.parametrize("noise_cv", [0.0, 0.001])
def test_every_fresh_prefix_has_explicit_bounded_initialization(
    reference_student, nuisance_validation, noise_cv,
):
    rejected = []
    examined = 0
    projected = 0
    before = json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False)
    for (cv, growth_name, density), episodes in nuisance_validation.items():
        if cv != noise_cv:
            continue
        for episode in episodes:
            raw = reference_student.infer(episode.times_h, episode.observations)
            for stop in range(1, len(episode.times_h)):
                initial = reference_student.initialize(raw[:stop])
                examined += 1
                np.testing.assert_array_equal(initial.raw_state, raw[stop - 1])
                assert np.all((initial.state >= 0.0) & (initial.state <= 1.0))
                diagnostics = initial.diagnostics()
                assert diagnostics["prefix_initialization_accepted"] == initial.accepted
                assert diagnostics["prefix_projection_l2"] == pytest.approx(
                    np.linalg.norm(initial.state - initial.raw_state),
                )
                if not initial.accepted:
                    rejected.append((episode.episode_id, stop, diagnostics))
                if diagnostics["prefix_projection_count"]:
                    projected += 1
                    assert diagnostics["prefix_raw_max_bound_violation"] < 0.1
                    assert diagnostics["prefix_standardized_residual"] <= diagnostics["prefix_residual_limit"]
                    forecast = reference_student.forecast(
                        episode.times_h[:stop + 1], episode.inputs[:stop + 1],
                        episode.times_h[:stop], episode.observations[:stop],
                    )
                    np.testing.assert_allclose(forecast[:stop], raw[:stop], rtol=0, atol=1e-13)
                    assert np.all((forecast[-1] >= 0.0) & (forecast[-1] <= 1.0))
            np.testing.assert_array_equal(
                raw, reference_student.infer(episode.times_h, episode.observations),
            )
    assert examined == 12 * 7 * 2 * 80
    assert rejected == []
    assert (projected > 0) == (noise_cv > 0)
    assert json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False) == before


def test_initialization_is_covariance_weighted_map_not_coordinate_clipping(student):
    covariance = 1e-4 * np.array([[4.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    weighted = replace(student, encoder_residual_covariance=covariance)
    raw = np.array([[0.1, 0.2, 0.3], [-0.01, 0.2, 0.3]])
    initial = weighted.initialize(raw)
    assert initial.accepted
    np.testing.assert_array_equal(initial.raw_state, raw[-1])
    np.testing.assert_allclose(initial.state, [0.0, 0.2025, 0.3], rtol=0, atol=1e-12)
    delta = initial.state - initial.raw_state
    assert initial.standardized_residual == pytest.approx(np.sqrt(delta @ np.linalg.solve(covariance, delta)))
    gross = weighted.initialize(np.array([[0.1, 0.2, 0.3], [-0.5, 0.2, 0.3]]))
    assert not gross.accepted
    assert gross.diagnostics()["prefix_raw_max_bound_violation"] == 0.5
    assert gross.diagnostics()["prefix_projection_l2"] >= 0.5
    for invalid in ([], np.ones(3), np.ones((1, 2)), np.full((2, 3), np.nan)):
        with pytest.raises(ValueError):
            weighted.initialize(invalid)


def test_scoring_keeps_raw_errors_and_initialization_diagnostics(
    reference_student, nuisance_validation, monkeypatch,
):
    episodes = nuisance_validation[(0.001, "constant_0.075", 0.1)][:2]
    original = InSilicoStudent.initialize
    calls = []

    def recording(self, inferred_prefix):
        assert len(inferred_prefix) == 3
        result = original(self, inferred_prefix)
        calls.append(result)
        return result

    monkeypatch.setattr(InSilicoStudent, "initialize", recording)
    before = json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False)
    scores = score_student(reference_student, episodes, prefix_h=0.2)
    assert len(calls) == 2 * len(episodes)
    assert not scores.forecast_rejected.any()
    assert scores.prefix_initialization_accepted.all()
    assert scores.prefix_projection_count.sum() > 0
    for episode, row in zip(episodes, scores.to_dict("records")):
        raw = reference_student.infer(episode.times_h, episode.observations)
        initial = original(reference_student, raw[:3])
        assert row["state_rmse"] == pytest.approx(_rmse(raw, episode.latent))
        for name, value in initial.diagnostics().items():
            assert row[name] == pytest.approx(value)
        assert row["state_bound_excursion_count"] == np.count_nonzero((raw < 0.0) | (raw > 1.0))
    assert json.dumps(reference_student.metadata, sort_keys=True, allow_nan=False) == before


def test_large_incompatible_prefix_is_rejected_and_counted_without_hiding_raw_error(
    student, holdout, monkeypatch,
):
    raw_state = np.array([-10.0, 0.3, 0.2])
    monkeypatch.setattr(
        InSilicoStudent, "infer", lambda self, times_h, observations: np.tile(raw_state, (len(times_h), 1)),
    )
    episode = holdout[0]
    initial = student.initialize(np.tile(raw_state, (21, 1)))
    assert not initial.accepted
    with pytest.raises(ValueError, match="prefix.*outside.*domain.*violation"):
        student.forecast(episode.times_h, episode.inputs, episode.times_h[:21], episode.observations[:21])
    scores = score_student(student, [episode], prefix_h=2.0)
    assert scores.forecast_rejected.sum() == 1
    assert not scores.prefix_initialization_accepted.any()
    assert scores.state_rmse.iloc[0] > 1.0
    assert scores.forecast_rmse.isna().all()
    assert scores.prefix_raw_max_bound_violation.iloc[0] == 10.0


@pytest.mark.parametrize("noise_cv", [0.0, 0.001])
def test_density_scale_does_not_amplify_into_forecast_differences(
    reference_student, nuisance_validation, noise_cv,
):
    for (cv, growth_name, density), episodes in nuisance_validation.items():
        if cv != noise_cv or density != 0.03:
            continue
        other = nuisance_validation[(cv, growth_name, 0.1)]
        for first, second in zip(episodes, other):
            predicted = [
                reference_student.forecast(
                    episode.times_h, episode.inputs, episode.times_h[:21], episode.observations[:21],
                )
                for episode in (first, second)
            ]
            np.testing.assert_allclose(predicted[0], predicted[1], rtol=0, atol=1e-11)
            assert _rmse(predicted[0][21:], first.latent[21:]) < 0.03


def test_observation_refitting_cannot_change_independently_supervised_dynamics_or_controls(
    training, student,
):
    reobserved = fit_student(_corrupt(training, "observations"), seed=7)
    for name in ("dynamics", "controller"):
        assert reobserved.metadata["learned_models"][name] == student.metadata["learned_models"][name]
    np.testing.assert_array_equal(reobserved.control_target_scale, student.control_target_scale)
    assert not np.allclose(reobserved.encoder.coefficients, student.encoder.coefficients)
