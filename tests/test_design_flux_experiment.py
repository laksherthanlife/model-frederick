from __future__ import annotations

import importlib.util
import math
import pathlib

import numpy as np
import pandas as pd
import pytest


def _load_script(name):
    filename = pathlib.Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"{name}_regressions", filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def script():
    return _load_script("design_flux_experiment")


@pytest.fixture
def identity_panel(script, monkeypatch):
    def panel(n_strains, alpha, sigma_log, spread, rng):
        expression = np.arange(1, n_strains + 1, dtype=float)
        return expression, expression.copy()

    monkeypatch.setattr(script, "_panel", panel)
    monkeypatch.setattr(script, "_skill",
                        lambda expression, flux, strains: float(np.array_equal(expression, flux)))
    return script


def test_sampled_p_includes_the_observed_assignment(script, monkeypatch):
    monkeypatch.setattr(script, "_EXACT_PERMUTATION_LIMIT", 10)
    calls = []

    def skill(expression, flux, strains):
        calls.append(tuple(expression))
        return 1.0 if len(calls) <= 2 else 0.0

    monkeypatch.setattr(script, "_skill", skill)
    result = script.power_at(4, 1e-3, 0.2, 4.0, 1,
                             np.random.default_rng(0), alpha_level=0.1)

    assert len(calls) == 11
    assert result["power"] == 0.0
    assert result["permutation_method"] == "monte_carlo_with_replacement"
    assert result["permutations_per_trial"] == 10
    assert result["p_value_resolution"] == pytest.approx(1 / 11)


def test_a_zero_sampled_tail_does_not_clear_an_unattainable_level(identity_panel, monkeypatch):
    monkeypatch.setattr(identity_panel, "_EXACT_PERMUTATION_LIMIT", 10)
    result = identity_panel.power_at(8, 1e-3, 0.0, 4.0, 1,
                                     np.random.default_rng(0), alpha_level=0.01)

    assert result["power"] == 0.0
    assert result["permutation_floor"] < 0.01
    assert result["p_value_resolution"] == pytest.approx(1 / 11)
    assert "UNREACHABLE" in result["verdict"]
    assert "resolution" in result["verdict"]


def test_exact_enumeration_counts_the_observed_assignment_once(identity_panel):
    result = identity_panel.power_at(4, 1e-3, 0.0, 4.0, 1,
                                     np.random.default_rng(0), alpha_level=1 / 24)

    assert result["power"] == 1.0
    assert result["permutation_method"] == "exact"
    assert result["permutations_per_trial"] == math.factorial(4)
    assert result["p_value_resolution"] == pytest.approx(1 / 24)
    assert result["permutation_floor"] == pytest.approx(1 / 24)


@pytest.mark.parametrize("n_strains", [4, 7])
def test_all_ties_have_p_one_in_both_modes(script, monkeypatch, n_strains):
    monkeypatch.setattr(script, "_skill", lambda expression, flux, strains: 1.0)
    result = script.power_at(n_strains, 1e-3, 0.2, 4.0, 1,
                             np.random.default_rng(0), alpha_level=0.99)

    assert result["power"] == 0.0


def test_large_panels_sample_without_constructing_factorial_orders(identity_panel, monkeypatch):
    monkeypatch.setattr(identity_panel, "_EXACT_PERMUTATION_LIMIT", 10)

    def forbid_enumeration(*args, **kwargs):
        pytest.fail("large panels must not enumerate permutations before sampling")

    monkeypatch.setattr(identity_panel.itertools, "permutations", forbid_enumeration)
    result = identity_panel.power_at(32, 1e-3, 0.0, 4.0, 1,
                                     np.random.default_rng(0), alpha_level=0.5)

    assert result["permutation_method"] == "monte_carlo_with_replacement"
    assert result["permutations_per_trial"] == 10
    assert result["p_value_resolution"] == pytest.approx(1 / 11)
    assert result["power"] == 1.0


def test_large_theoretical_floor_does_not_overflow(script):
    assert script.permutation_floor(1000) == 0.0


def test_unreachable_exact_design_still_reports_its_resolution(script):
    result = script.power_at(3, 1e-3, 0.2, 4.0, 1, np.random.default_rng(0))

    assert result["power"] == 0.0
    assert result["permutation_method"] == "exact"
    assert result["permutations_per_trial"] == 6
    assert result["p_value_resolution"] == pytest.approx(1 / 6)


def test_monte_carlo_power_is_reproducible_with_the_same_seed(script):
    first = script.power_at(7, 1e-3, 0.2, 4.0, 3, np.random.default_rng(12))
    second = script.power_at(7, 1e-3, 0.2, 4.0, 3, np.random.default_rng(12))

    assert first == second


@pytest.mark.parametrize("n_strains", [0, -1, 2.5, True])
def test_permutation_floor_requires_a_positive_integer(script, n_strains):
    with pytest.raises(ValueError, match="positive integer"):
        script.permutation_floor(n_strains)


@pytest.mark.parametrize("n_trials", [0, -1, 1.5, True])
def test_power_requires_a_positive_integer_trial_count(script, n_trials):
    with pytest.raises(ValueError, match="n_trials must be a positive integer"):
        script.power_at(4, 1e-3, 0.2, 4.0, n_trials, np.random.default_rng(0))


def test_the_output_records_method_draw_count_and_resolution(script, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(script, "_EXACT_PERMUTATION_LIMIT", 24)
    monkeypatch.setattr(script.paths, "outputs_dir", lambda: tmp_path)
    monkeypatch.setattr(script.sys, "argv", ["design_flux_experiment.py", "--trials", "1"])

    assert script.main() == 0
    output = capsys.readouterr().out
    table = pd.read_csv(tmp_path / "flux_experiment_design.csv")
    assert set(table.permutation_method) == {"exact", "monte_carlo_with_replacement"}
    assert table.permutations_per_trial.between(1, 24).all()
    assert (table.p_value_resolution > 0).all()
    assert "permutation_method" in output
    assert "p_value_resolution" in output


def test_error_budget_labels_assumptions_without_changing_residuals(monkeypatch, tmp_path, capsys):
    module = _load_script("error_budget")
    monkeypatch.setattr(module.paths, "outputs_dir", lambda: tmp_path)
    expected = module.residuals(module.states())

    assert module.main() == 0
    output = capsys.readouterr().out
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / "error_budget.csv"), expected)
    assert "fixed-CrtE comparison" in output
    assert "conditional simulation assumption" in output
    assert "estimation error" in output
    assert "leaves prediction error flat" not in output
    assert "MEASUREMENT-LIMITED:" not in output
    assert "not a measurement-limit guarantee" in output
    assert "fixed-CrtE comparison" in module.__doc__


def test_error_budget_remains_a_descriptive_squared_residual_split():
    module = _load_script("error_budget")
    residuals = pd.DataFrame({"strain": ["a", "a", "b", "b"],
                              "log_residual": [0.2, 0.1, -0.3, -0.2]})

    assert module.budget(residuals) == pytest.approx({
        "total": 0.045, "between_strain": 0.0425, "within_strain": 0.0025})
