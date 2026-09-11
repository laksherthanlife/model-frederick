"""Sizing the chemostat run by simulation -- and the line simulation cannot cross.

`scripts/design_flux_experiment.py` answers "how many strains" and refuses to answer "is the
law true". The second refusal is the one worth a test, because the script looks exactly like
a validation and is not one: its generator and its model are the same equation, so it would
report a good fit for a law the organism ignores.
"""

from __future__ import annotations

import importlib.util
import math
import pathlib

import numpy as np
import pytest

_SCRIPT = (pathlib.Path(__file__).resolve().parents[1]
           / "scripts" / "design_flux_experiment.py")


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("design_flux_experiment", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheResolutionBarIsArithmeticNotPower:
    """The exact expression-profile permutation floor is 1/n! before data are collected.
    It bounds this test, not every possible inferential method on the same strains."""

    def test_three_strains_cannot_reach_five_percent(self, script):
        assert script.permutation_floor(3) == pytest.approx(1 / 6)
        assert script.permutation_floor(3) > 0.05

    def test_four_strains_is_the_first_size_that_can(self, script):
        assert script.permutation_floor(4) == pytest.approx(1 / 24)
        assert script.permutation_floor(4) < 0.05

    def test_the_floor_is_one_over_n_factorial(self, script):
        for n in range(2, 8):
            assert script.permutation_floor(n) == pytest.approx(1 / math.factorial(n))


class TestThePowerCurveIsUsableAsAProtocol:
    def test_a_panel_below_the_floor_reports_zero_power_not_a_number(self, script):
        """A design that cannot express the level must say so rather than return a power
        computed against a threshold it can never reach."""
        result = script.power_at(3, alpha=1e-3, sigma_log=0.2, spread=4.0,
                                 n_trials=5, rng=np.random.default_rng(0))

        assert result["power"] == 0.0
        assert "UNREACHABLE" in result["verdict"]

    def test_power_rises_with_panel_size(self, script):
        rng = np.random.default_rng(1)
        small = script.power_at(4, 1e-3, 0.2, 4.0, 40, rng)
        large = script.power_at(6, 1e-3, 0.2, 4.0, 40, rng)

        assert large["power"] > small["power"]

    def test_six_strains_clears_eighty_percent_at_the_measured_scatter(self, script):
        """Conditional simulation at a fixed seed, noise law and expression spread.
        This is not a universal six-strain protocol or measured biological power."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX

        result = script.power_at(6, BETA_CAROTENE_FLUX.alpha,
                                 BETA_CAROTENE_FLUX.loso_rmse_log, 4.0, 60,
                                 np.random.default_rng(2))

        assert result["power"] >= 0.8


class TestTheSimulatorCannotValidateTheLaw:
    """The whole point. Simulation sizes an experiment; it cannot be one."""

    def test_the_generator_and_the_model_are_the_same_equation(self, script):
        """Data drawn from `flux = alpha * expression` and fitted with
        `flux = alpha * expression` recovers alpha whatever the organism does. This asserts
        the circularity exists, so nobody reads a high simulated skill as evidence."""
        rng = np.random.default_rng(3)
        alpha = 2.5e-3

        expression, flux = script._panel(6, alpha, sigma_log=1e-9, spread=4.0, rng=rng)
        recovered = float(np.median(flux / expression))

        assert recovered == pytest.approx(alpha, rel=1e-6)

    def test_noiseless_simulated_data_scores_near_perfectly(self, script):
        """Tier 1 simulator verification, not validation of a biological expression law.
        Even the real-data fixed-CrtE comparison does not validate selecting that gene."""
        rng = np.random.default_rng(4)
        expression, flux = script._panel(6, 1e-3, sigma_log=1e-9, spread=4.0, rng=rng)

        skill = script._skill(expression, flux, [f"s{i}" for i in range(6)])

        assert skill > 0.99


def _comparison_script(name):
    filename = _SCRIPT.with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_target_runner_does_not_call_separate_refits_a_conversion_error(monkeypatch, tmp_path,
                                                                      capsys):
    module = _comparison_script("score_targets")
    monkeypatch.setattr(module.paths, "outputs_dir", lambda: tmp_path)

    assert module.main() == 0
    output = capsys.readouterr().out
    assert "The derivation costs" not in output
    assert "separate target refits" in output
    assert "cancel exactly" in output
    assert "fixed-gene comparison" in output


def test_condition_sweep_identifies_the_fixed_empirical_comparison():
    from ystwin.pathway.spec import load_pathway

    module = _comparison_script("sweep_conditions")
    table = module.sweep(load_pathway("beta_carotene"), feed_g_per_L=20.0)

    assert set(table.model) == {"fixed_CrtE_rate_law"}
    assert table.prediction_scope.str.contains("exploratory").all()
    assert "not a causal" in module.__doc__
    assert "least predictable thing" not in module.__doc__


def test_proteomics_runner_does_not_declare_every_panel_size_equally_accurate(monkeypatch,
                                                                         tmp_path, capsys):
    import sys

    module = _comparison_script("simulate_proteomics_gain")
    monkeypatch.setattr(module.paths, "outputs_dir", lambda: tmp_path)
    monkeypatch.setattr(sys, "argv", ["simulate_proteomics_gain.py", "--trials", "2"])
    monkeypatch.setattr(module, "one_panel",
                        lambda n_strains, rates, sigma_protein, rng:
                        (1 + 1 / n_strains, 1.5 if sigma_protein >= 0.2 else 1.05))

    assert module.main() == 0
    output = capsys.readouterr().out
    assert "Flat." not in output
    assert "not ACCURACY" not in output
    assert "estimation error" in output
    assert "assumed" in output
    assert "first simulated" in output


@pytest.mark.parametrize("pathway", ["beta_carotene", "glycogen"])
def test_pathway_narratives_do_not_claim_all_native_genes_scored_negative(pathway):
    import tomllib

    filename = _SCRIPT.parents[1] / "data" / "pathways" / f"{pathway}.toml"
    text = filename.read_text()
    specification = tomllib.loads(text)

    assert "every native gene NEGATIVE" not in text
    assert "every\n# native mevalonate gene scores NEGATIVE" not in text
    assert "every NATIVE yeast gene scoring -0.46 to -0.86" not in text
    assert "CrtE" in specification["pathway"]["notes"]
    assert "ERG9" in specification["pathway"]["notes"]
    assert "causal" in specification["pathway"]["notes"]


def test_glycogen_scope_is_missing_turnover_calibration_not_a_missing_solver():
    from ystwin.pathway.solve import NodeKinetics, solve_pathway
    from ystwin.pathway.spec import load_pathway

    spec = load_pathway("glycogen")
    with pytest.raises(ValueError, match="degradation_rate_per_h"):
        solve_pathway(spec, 1e-3, 0.025, {})
    solved = solve_pathway(spec, 1e-3, 0.025,
                           {"glycogen": NodeKinetics(degradation_rate_per_h=0.1)})

    assert solved.terminal.content_mmol_per_gdcw == pytest.approx(1e-3 / (0.025 + 0.1))
    assert "degradation_rate_per_h" in spec.notes
    assert "NOT PREDICTABLE by this architecture, in principle" not in spec.notes


def test_feed_scaling_does_not_silently_override_the_product_environment(monkeypatch):
    from ystwin.generator.context import CultureContext
    from ystwin.pathway.spec import load_pathway

    module = _comparison_script("sweep_conditions")
    predict = module.predict_product
    contexts = []

    def record(spec, genotype, environment, *args):
        contexts.append(environment.context)
        return predict(spec, genotype, environment, *args)

    monkeypatch.setattr(module, "predict_product", record)
    table = module.sweep(load_pathway("beta_carotene"), feed_g_per_L=10.0)

    assert table.answered.any()
    assert all(context == CultureContext() for context in contexts)
    assert set(table.biomass_feed_g_per_L) == {10.0}
    assert table.biomass_scope.str.contains("conditional scaling only").all()


@pytest.mark.parametrize("feed", [0.0, -1.0, float("nan"), float("inf")])
def test_conditional_biomass_requires_a_positive_finite_feed(feed):
    from ystwin.pathway.spec import load_pathway

    module = _comparison_script("sweep_conditions")
    with pytest.raises(ValueError, match="positive and finite"):
        module.sweep(load_pathway("beta_carotene"), feed_g_per_L=feed)
