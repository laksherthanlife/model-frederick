import builtins
import io
import json
import os
from dataclasses import FrozenInstanceError, fields, replace
from itertools import product
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid
from scipy.special import expit

from ystwin.generator import in_silico
from ystwin.generator.in_silico import (
    CONTROL_FEATURE_NAMES,
    CONTROL_NAMES,
    LATENT_NAMES,
    OBSERVATION_NAMES,
    SyntheticEpisode,
    TeacherParameters,
    generate_episodes,
    observe_reporters,
    simulate_teacher,
    teacher_controls,
    teacher_rhs,
)
from ystwin.reporter import ReporterKinetics, simulate_reporter


def test_frozen_episode_contract_and_json_provenance():
    t = np.linspace(0.0, 1.0, 11)
    inputs = np.full((t.size, 2), 0.4)
    episode = simulate_teacher(t, inputs, episode_id="contract", noise_cv=0.0)
    assert LATENT_NAMES == ("upr", "oxidative", "burden")
    assert CONTROL_NAMES == (
        "growth_retention", "enzyme_budget_scale", "ngam_mmol_per_gdcw_h", "allocation_fraction"
    )
    assert OBSERVATION_NAMES == (
        "cell_density", "reporter_1", "reporter_2", "reporter_3", "reporter_4"
    )
    assert [field.name for field in fields(SyntheticEpisode)] == [
        "episode_id", "times_h", "inputs", "observations", "latent", "latent_derivative",
        "controls", "growth_rate_per_h", "metadata",
    ]
    assert episode.episode_id == "contract"
    for name, shape in (
        ("times_h", (11,)), ("inputs", (11, 2)), ("observations", (11, 5)),
        ("latent", (11, 3)), ("latent_derivative", (11, 3)),
        ("controls", (11, 4)), ("growth_rate_per_h", (11,)),
    ):
        array = getattr(episode, name)
        assert array.shape == shape
        assert np.all(np.isfinite(array))
    with pytest.raises(FrozenInstanceError):
        episode.episode_id = "changed"
    assert not np.shares_memory(episode.observations, episode.latent)
    inputs[:] = 0.9
    t[:] = 9.0
    assert episode.inputs[0, 0] == 0.4
    assert episode.times_h[0] == 0.0
    metadata = json.loads(json.dumps(episode.metadata, allow_nan=False))
    assert metadata["product_targets_used"] is False
    assert metadata["biological_validation"] is False
    assert metadata["truth_available_only_in_simulation"] is True
    assert metadata["input_units"] == "normalized imposed stress levels in [0, 1]; not mM"
    assert metadata["observation_units"] == [
        "normalized synthetic cell-density proxy",
        *["normalized synthetic instrument units per cell"] * 4,
    ]
    provenance = json.loads(TeacherParameters().to_json())
    assert provenance["parameter_source"] == "declared synthetic mathematical teacher and library priors"
    assert provenance["biological_validation"] is False
    assert provenance["product_targets_used"] is False
    assert {"bilinear", "saturating", "reporter", "density", "controls", "noise"} <= set(provenance["equations"])
    assert provenance["parameters"]["reporter_kinetics"]["k_deg"] == 0.05


def test_bilinear_rhs_has_the_declared_equations_and_delayed_burden():
    p = TeacherParameters()
    u, o, b = state = np.array([0.2, 0.3, 0.1])
    s1, s2 = inputs = np.array([0.7, 0.5])
    expected = [
        (s1 * (1.0 - u) - u) / p.tau_u_h + p.cross_uo_per_h * o * (1.0 - u),
        (s2 * (1.0 - o) - o) / p.tau_o_h,
        (p.burden_weights[0] * u + p.burden_weights[1] * o - b) / p.tau_b_h,
    ]
    np.testing.assert_allclose(teacher_rhs(state, inputs, p), expected)
    derivative = teacher_rhs(np.zeros(3), np.ones(2))
    assert np.all(derivative[:2] > 0.0)
    assert derivative[2] == 0.0
    t = np.linspace(0.0, 2.0, 41)
    episode = simulate_teacher(t, np.ones((t.size, 2)), noise_cv=0.0)
    assert 0.0 < episode.latent[1, 2] < 0.1 * episode.latent[1, 0]
    assert episode.latent[-1, 2] > episode.latent[1, 2]
    np.testing.assert_allclose(
        episode.latent_derivative,
        [teacher_rhs(z, s, p) for z, s in zip(episode.latent, episode.inputs)],
    )


@pytest.mark.parametrize("family", ["bilinear", "saturating"])
def test_vector_field_points_inward_on_every_cube_face(family):
    p = TeacherParameters(family=family)
    for state in product((0.0, 1.0), repeat=3):
        z = np.asarray(state)
        for inputs in product((0.0, 1.0), repeat=2):
            derivative = teacher_rhs(z, inputs, p)
            assert np.all(derivative[z == 0.0] >= 0.0)
            assert np.all(derivative[z == 1.0] <= 0.0)


def test_true_ode_solution_uses_elapsed_time_and_declared_nonzero_initial_state():
    t = np.array([2.0, 2.03, 2.15, 2.5, 3.0, 4.0])
    inputs = np.full((t.size, 2), [0.4, 0.6])
    initial = np.array([0.2, 0.3, 0.1])
    p = replace(TeacherParameters(), cross_uo_per_h=0.0)
    episode = simulate_teacher(t, inputs, parameters=p, initial_latent=initial, noise_cv=0.0)
    elapsed = t - t[0]
    for column, tau in ((0, p.tau_u_h), (1, p.tau_o_h)):
        value = inputs[0, column]
        steady = value / (1.0 + value)
        exact = steady + (initial[column] - steady) * np.exp(-(1.0 + value) * elapsed / tau)
        np.testing.assert_allclose(episode.latent[:, column], exact, atol=2e-9, rtol=2e-9)
    np.testing.assert_array_equal(episode.latent[0], initial)
    assert episode.metadata["initial_latent"] == initial.tolist()
    shifted = simulate_teacher(elapsed, inputs, parameters=p, initial_latent=initial, noise_cv=0.0)
    stretched = simulate_teacher(elapsed * 2.0, inputs, parameters=p, initial_latent=initial, noise_cv=0.0)
    np.testing.assert_allclose(episode.latent, shifted.latent, atol=2e-9)
    assert not np.allclose(episode.latent, stretched.latent)


@pytest.mark.parametrize("family", ["bilinear", "saturating"])
def test_piecewise_linear_pulse_cannot_change_truth_before_its_first_rising_interval(family):
    t = np.arange(11) * 0.1
    inputs = np.zeros((t.size, 2))
    inputs[3:6, 0] = 0.9
    inputs[4:7, 1] = 0.7
    episode = simulate_teacher(t, inputs, family=family, noise_cv=0.0)
    np.testing.assert_array_equal(episode.latent[:3], 0.0)
    assert episode.latent[3, 0] > 0.0
    assert episode.latent[4, 1] > 0.0


def test_zero_inputs_leave_zero_truth_but_positive_density_and_basal_reporters():
    t = np.linspace(0.0, 2.0, 21)
    p = TeacherParameters()
    episode = simulate_teacher(t, np.zeros((t.size, 2)), noise_cv=0.0)
    np.testing.assert_array_equal(episode.latent, 0.0)
    np.testing.assert_array_equal(episode.latent_derivative, 0.0)
    np.testing.assert_allclose(episode.growth_rate_per_h, p.assay_growth_rate_per_h)
    np.testing.assert_allclose(
        episode.observations[:, 0], p.initial_cell_density * np.exp(p.assay_growth_rate_per_h * t)
    )
    steady = np.asarray(p.reporter_basal_synthesis_per_h) / (p.assay_growth_rate_per_h + p.reporter_kinetics.k_deg)
    np.testing.assert_allclose(episode.observations[:, 1:], np.broadcast_to(steady, (t.size, 4)))
    assert episode.metadata["growth_source"] == "declared_synthetic_assay"
    assert episode.metadata["density_source"] == "integrated_growth"


def test_common_control_laws_are_frozen_quadratic_state_feature_maps():
    p = TeacherParameters()
    assert CONTROL_FEATURE_NAMES == (
        "1", "upr", "oxidative", "burden", "upr^2", "oxidative^2", "burden^2",
        "upr*oxidative", "upr*burden", "oxidative*burden",
    )
    assert p.growth_cost_coefficients == (0.0, 0.18, 0.14, 0.25, 0.06, 0.05, 0.08, 0.04, 0.0, 0.0)
    assert p.enzyme_cost_coefficients == (0.0, 0.07, 0.06, 0.22, 0.02, 0.02, 0.04, 0.02, 0.0, 0.0)
    assert p.ngam_coefficients == (0.7, 0.10, 0.12, 0.18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert p.allocation_logit_coefficients == (-1.6, -0.25, -0.20, -0.35, -0.08, -0.06, -0.12, -0.08, -0.05, -0.04)
    z = np.random.default_rng(4).uniform(size=(2, 5, 3))
    u, o, b = np.moveaxis(z, -1, 0)
    phi = np.stack((np.ones_like(u), u, o, b, u*u, o*o, b*b, u*o, u*b, o*b), axis=-1)
    controls = teacher_controls(z)
    assert controls.shape == (2, 5, 4)
    np.testing.assert_allclose(controls[..., 0], np.exp(-phi @ p.growth_cost_coefficients))
    np.testing.assert_allclose(controls[..., 1], np.exp(-phi @ p.enzyme_cost_coefficients))
    np.testing.assert_allclose(controls[..., 2], phi @ p.ngam_coefficients)
    np.testing.assert_allclose(controls[..., 3], expit(phi @ p.allocation_logit_coefficients))
    assert np.all((controls[..., :2] > 0.0) & (controls[..., :2] <= 1.0))
    assert np.all(controls[..., 2] >= 0.7)
    assert np.all((controls[..., 3] > 0.0) & (controls[..., 3] < 1.0))
    assert teacher_controls(np.zeros(3)).shape == (4,)
    with pytest.raises((ValueError, TypeError)):
        replace(p, growth_cost_coefficients=(0.0,) * 10)


def test_observations_reuse_reporter_ode_and_same_culture_growth_and_density(monkeypatch):
    t = np.linspace(0.0, 2.0, 21)
    inputs = np.column_stack((t / 2.0, np.full(t.size, 0.4)))
    growth = 0.18 + 0.02 * t
    density = 0.3 * np.exp(cumulative_trapezoid(growth, t, initial=0.0))
    calls = []

    def recording_reporter(*args, **kwargs):
        calls.append((args, kwargs))
        return simulate_reporter(*args, **kwargs)

    monkeypatch.setattr(in_silico, "simulate_reporter", recording_reporter)
    episode = simulate_teacher(
        t, inputs, growth_rate_per_h=growth, cell_density=density, noise_cv=0.0
    )
    assert len(calls) == 4
    np.testing.assert_array_equal(episode.growth_rate_per_h, growth)
    np.testing.assert_array_equal(episode.observations[:, 0], density)
    p = TeacherParameters()
    loadings = np.asarray(p.reporter_loadings)
    assert loadings.shape == (4, 3)
    assert np.linalg.matrix_rank(loadings) == 3
    synthesis = np.asarray(p.reporter_basal_synthesis_per_h) + episode.latent @ loadings.T
    for channel in range(4):
        expected = simulate_reporter(t, synthesis[:, channel], growth, p.reporter_kinetics)
        np.testing.assert_allclose(episode.observations[:, channel + 1], expected)
    assert not np.allclose(episode.observations[:, 1:4], episode.latent)
    assert not np.allclose(episode.observations[:, 1:], synthesis / (growth[:, None] + 0.05))
    assert episode.metadata["growth_source"] == "supplied_same_culture"
    assert episode.metadata["density_source"] == "supplied_same_culture"
    resynthesized = observe_reporters(t, episode.latent, growth, density, noise_cv=0.0)
    np.testing.assert_array_equal(resynthesized, episode.observations)
    synthetic = simulate_teacher(t, inputs, noise_cv=0.0)
    np.testing.assert_array_equal(episode.latent, synthetic.latent)
    assert not np.allclose(episode.observations[:, 1:], synthetic.observations[:, 1:])


def test_supplied_growth_alone_integrates_density_and_zero_growth_is_valid():
    t = np.linspace(0.0, 1.0, 11)
    growth = np.zeros(t.size)
    episode = simulate_teacher(t, np.ones((t.size, 2)), growth_rate_per_h=growth, noise_cv=0.0)
    np.testing.assert_array_equal(episode.growth_rate_per_h, growth)
    np.testing.assert_allclose(episode.observations[:, 0], TeacherParameters().initial_cell_density)
    assert np.all(episode.observations[:, 1:] > 0.0)


def test_reporter_prefix_is_causal_and_maturation_is_exposed():
    t = np.linspace(0.0, 3.0, 31)
    p = replace(TeacherParameters(), reporter_kinetics=ReporterKinetics(k_deg=0.05, k_mat=0.6))
    episode = simulate_teacher(t, np.full((t.size, 2), 0.7), parameters=p, noise_cv=0.0)
    clean = observe_reporters(t, episode.latent, episode.growth_rate_per_h, episode.observations[:, 0], parameters=p, noise_cv=0.0)
    prefix = observe_reporters(t[:16], episode.latent[:16], episode.growth_rate_per_h[:16], episode.observations[:16, 0], parameters=p, noise_cv=0.0)
    np.testing.assert_allclose(prefix, clean[:16], rtol=2e-7, atol=2e-8)
    instantaneous = observe_reporters(t, episode.latent, episode.growth_rate_per_h, episode.observations[:, 0], noise_cv=0.0)
    assert not np.allclose(clean[:, 1:], instantaneous[:, 1:])


def test_noise_is_reproducible_positive_and_does_not_change_truth():
    t = np.linspace(0.0, 1.0, 11)
    inputs = np.full((t.size, 2), 0.5)
    first = simulate_teacher(t, inputs, seed=22, noise_cv=1.0)
    replay = simulate_teacher(t, inputs, seed=22, noise_cv=1.0)
    other = simulate_teacher(t, inputs, seed=23, noise_cv=1.0)
    np.testing.assert_array_equal(first.observations, replay.observations)
    assert first.metadata == replay.metadata
    assert not np.array_equal(first.observations, other.observations)
    for name in ("latent", "latent_derivative", "controls", "growth_rate_per_h"):
        np.testing.assert_array_equal(getattr(first, name), getattr(other, name))
    assert np.all(first.observations > 0.0)
    assert first.metadata["noise"]["model"] == "independent_mean_preserving_lognormal"
    assert first.metadata["noise"]["cv"] == 1.0


@pytest.mark.parametrize("family", ["bilinear", "saturating"])
@pytest.mark.parametrize("protocol", ["random", "pulse", "ramp", "periodic"])
def test_generated_protocols_are_bounded_reproducible_and_independently_seeded(family, protocol):
    first = generate_episodes(2, seed=5, family=family, protocol=protocol, hours=1.0, dt=0.1)
    longer = generate_episodes(3, seed=5, family=family, protocol=protocol, hours=1.0, dt=0.1)
    assert isinstance(first, tuple)
    assert len({episode.episode_id for episode in first}) == 2
    assert len({episode.metadata["seed"] for episode in first}) == 2
    assert not np.array_equal(first[0].inputs, first[1].inputs)
    for index, episode in enumerate(first):
        assert episode.episode_id == f"{family}/{protocol}/5/{index}"
        assert episode.metadata["input_interpolation"] == "piecewise_linear"
        assert episode.metadata["split_unit"] == "episode_id"
        assert episode.metadata["protocol"] == protocol
        assert np.all((episode.inputs >= 0.0) & (episode.inputs <= 1.0))
        assert np.all((episode.latent >= 0.0) & (episode.latent <= 1.0))
        assert np.all(episode.observations > 0.0)
        np.testing.assert_array_equal(episode.observations, longer[index].observations)
        np.testing.assert_array_equal(episode.inputs, longer[index].inputs)
        np.testing.assert_array_equal(episode.controls, teacher_controls(episode.latent))


def test_recorded_episode_and_protocol_seeds_can_be_replayed_directly():
    episode, = generate_episodes(1, seed=37, hours=1.0, dt=0.1)
    metadata = episode.metadata
    protocol_stream, noise_stream = np.random.SeedSequence(metadata["episode_seed"]).spawn(2)
    assert metadata["protocol_seed"] == int(protocol_stream.generate_state(1, dtype=np.uint64)[0])
    assert metadata["seed"] == int(noise_stream.generate_state(1, dtype=np.uint64)[0])
    expected_levels = np.random.default_rng(metadata["protocol_seed"]).uniform(0.0, 1.0, (9, 2))
    np.testing.assert_array_equal(metadata["protocol_parameters"]["knot_values"], expected_levels)


def test_parameter_family_is_not_silently_replaced_by_default_bilinear():
    times = np.linspace(0, 2, 21)
    inputs = np.full((len(times), 2), 0.3)
    parameters = TeacherParameters(family="saturating")
    configured = simulate_teacher(times, inputs, parameters=parameters, noise_cv=0)
    explicit = simulate_teacher(times, inputs, family="saturating", noise_cv=0)
    overridden = simulate_teacher(times, inputs, parameters=parameters, family="bilinear", noise_cv=0)

    np.testing.assert_array_equal(configured.latent, explicit.latent)
    assert configured.metadata["family"] == "saturating"
    assert overridden.metadata["family"] == "bilinear"
    assert not np.allclose(configured.latent, overridden.latent)


def test_episode_splits_have_disjoint_ids_and_saturating_family_changes_dynamics():
    train = generate_episodes(2, seed=1, hours=1.0, dt=0.2, noise_cv=0.0)
    test = generate_episodes(2, seed=2, hours=1.0, dt=0.2, noise_cv=0.0)
    assert not {e.episode_id for e in train} & {e.episode_id for e in test}
    reference = train[0]
    misspecified = simulate_teacher(reference.times_h, reference.inputs, family="saturating", noise_cv=0.0)
    assert not np.allclose(reference.latent, misspecified.latent)
    p = TeacherParameters(family="saturating")
    np.testing.assert_allclose(
        misspecified.latent_derivative,
        [teacher_rhs(z, s, p) for z, s in zip(misspecified.latent, misspecified.inputs)],
    )
    assert misspecified.metadata["family"] == "saturating"
    assert misspecified.metadata["parameters"]["family"] == "saturating"


@pytest.mark.parametrize("changes", [
    {"times_h": [0.0, 0.1, 0.2, 0.3]},
    {"times_h": [[0.0, 0.1, 0.2, 0.3, 0.4]]},
    {"times_h": [0.0, 0.1, 0.1, 0.3, 0.4]},
    {"times_h": [0.4, 0.3, 0.2, 0.1, 0.0]},
    {"times_h": [0.0, 0.1, np.nan, 0.3, 0.4]},
    {"times_h": [0.0, 0.1, 0.2, 0.3, np.inf]},
    {"inputs": np.zeros((5, 3))},
    {"inputs": np.zeros(5)},
    {"inputs": np.full((5, 2), np.nan)},
    {"inputs": np.full((5, 2), 0.5 + 0.1j)},
    {"inputs": np.full((5, 2), -0.01)},
    {"inputs": np.full((5, 2), 1.01)},
    {"initial_latent": [0.0, 0.0]},
    {"initial_latent": [0.0, np.inf, 0.0]},
    {"initial_latent": [0.0, -0.01, 0.0]},
    {"initial_latent": [1.01, 0.0, 0.0]},
    {"growth_rate_per_h": np.ones(4)},
    {"growth_rate_per_h": np.full(5, -0.1)},
    {"growth_rate_per_h": np.full(5, np.nan)},
    {"growth_rate_per_h": np.ones(5), "cell_density": np.ones(4)},
    {"growth_rate_per_h": np.ones(5), "cell_density": np.zeros(5)},
    {"growth_rate_per_h": np.ones(5), "cell_density": np.full(5, np.inf)},
    {"cell_density": np.ones(5)},
    {"noise_cv": -0.01},
    {"noise_cv": np.inf},
    {"family": "unknown"},
    {"seed": -1},
    {"seed": 1.5},
    {"episode_id": ""},
])
def test_simulation_rejects_invalid_shapes_grids_domains_and_nonfinite_values(changes):
    kwargs = {"times_h": np.linspace(0.0, 0.4, 5), "inputs": np.zeros((5, 2))}
    kwargs.update(changes)
    with pytest.raises(ValueError):
        simulate_teacher(**kwargs)


@pytest.mark.parametrize("changes", [
    {"latent": np.zeros((5, 2))},
    {"latent": np.full((5, 3), np.nan)},
    {"latent": np.full((5, 3), 1.01)},
    {"growth_rate_per_h": np.full(5, np.inf)},
    {"cell_density": np.zeros(5)},
])
def test_observer_validates_all_supplied_same_culture_series(changes):
    kwargs = {
        "times_h": np.linspace(0.0, 0.4, 5), "latent": np.zeros((5, 3)),
        "growth_rate_per_h": np.full(5, 0.3), "cell_density": np.ones(5),
    }
    kwargs.update(changes)
    with pytest.raises(ValueError):
        observe_reporters(**kwargs)


@pytest.mark.parametrize("changes", [
    {"tau_u_h": 0.0}, {"tau_o_h": -1.0}, {"tau_b_h": np.inf},
    {"cross_uo_per_h": -0.1}, {"burden_weights": (0.6, 0.6)},
    {"burden_weights": (-0.1, 0.5)}, {"saturation_half": 0.0},
    {"saturation_power": 0}, {"reporter_loadings": ((1.0, 1.0, 1.0),) * 4},
    {"reporter_basal_synthesis_per_h": (0.0, 0.1, 0.1, 0.1)},
    {"assay_growth_rate_per_h": -0.1}, {"initial_cell_density": 0.0},
    {"solver_rtol": 0.0}, {"family": "unknown"},
])
def test_parameters_validate_the_invariant_and_reporter_model(changes):
    with pytest.raises(ValueError):
        TeacherParameters(**changes)


@pytest.mark.parametrize("kinetics", [
    {"k_deg": np.nan}, {"k_deg": -0.1},
    {"k_mat": 0.0}, {"k_mat": -0.1}, {"k_mat": np.inf},
    {"k_mat": 0.6, "k_deg_immature": -0.1},
])
def test_invalid_reporter_kinetics_are_rejected_before_teacher_construction(kinetics):
    # ReporterKinetics now enforces its own domain at construction. Keep that
    # boundary inside the assertion, not in import-time parametrization.
    with pytest.raises(ValueError):
        TeacherParameters(reporter_kinetics=ReporterKinetics(**kinetics))


@pytest.mark.parametrize("bad", [np.zeros(2), [0.0, -0.1, 0.0], [0.0, np.nan, 0.0], [2.0, 0.0, 0.0]])
def test_rhs_and_controls_reject_invalid_latent_states(bad):
    with pytest.raises(ValueError):
        teacher_rhs(bad, [0.0, 0.0])
    with pytest.raises(ValueError):
        teacher_controls(bad)


@pytest.mark.parametrize("changes", [
    {"n": -1}, {"n": 1.5}, {"protocol": "unknown"},
    {"hours": 0.0}, {"hours": np.nan}, {"dt": 0.0}, {"dt": np.inf},
    {"hours": 0.3, "dt": 0.1},
])
def test_generation_rejects_invalid_counts_protocols_and_grids(changes):
    kwargs = {"n": 1}
    kwargs.update(changes)
    with pytest.raises(ValueError):
        generate_episodes(**kwargs)


def test_generation_handles_empty_batches_and_nonintegral_duration():
    assert generate_episodes(0) == ()
    episode, = generate_episodes(1, hours=1.05, dt=0.2)
    np.testing.assert_allclose(episode.times_h, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.05])


def test_integration_errors_are_rejected_instead_of_clipped(monkeypatch):
    def invalid_integration(*args, **kwargs):
        states = np.zeros((3, 5))
        states[0, 1] = -1e-7
        return SimpleNamespace(success=True, y=states, message="invalid test trajectory")

    monkeypatch.setattr(in_silico, "solve_ivp", invalid_integration)
    with pytest.raises(RuntimeError, match="bound|invariant"):
        simulate_teacher(np.linspace(0.0, 0.4, 5), np.zeros((5, 2)))


def test_teacher_does_not_read_any_outcome_or_calibration_files(monkeypatch):
    def forbidden_read(*args, **kwargs):
        raise AssertionError("The synthetic teacher must not read files")

    with monkeypatch.context() as guard:
        guard.setattr(builtins, "open", forbidden_read)
        guard.setattr(io, "open", forbidden_read)
        guard.setattr(os, "open", forbidden_read)
        for name in ("load", "loadtxt", "genfromtxt", "fromfile"):
            guard.setattr(np, name, forbidden_read)
        episode, = generate_episodes(1, hours=1.0, dt=0.2)
        assert episode.metadata["product_targets_used"] is False
        assert episode.metadata["calibration_data_used"] is False
