from __future__ import annotations

import math
import os
import sys
from copy import deepcopy
from dataclasses import FrozenInstanceError, fields, is_dataclass
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import RegularGridInterpolator

from ystwin import paths

pytest_plugins = ["pytester"]
_CALLER_ENVIRONMENT_NAMES = ("HOME", "USERPROFILE", "XDG_CACHE_HOME")


def _capture_caller_environment():
    return {name: os.environ.get(name) for name in _CALLER_ENVIRONMENT_NAMES}


_CALLER_ENVIRONMENT = _capture_caller_environment()


def _restore_caller_environment(monkeypatch, environment):
    for name, value in environment.items():
        if value is None:
            monkeypatch.delenv(name, raising=False)
        else:
            monkeypatch.setenv(name, value)


def _owners(model):
    return [model, *model.reactions, *model.metabolites, *model.genes, *model.groups]


def _minimal_model():
    model = cobra.Model("fixture_isolation")
    substrate = cobra.Metabolite("substrate_c", compartment="c")
    supply = cobra.Reaction("supply", lower_bound=0, upper_bound=10)
    supply.add_metabolites({substrate: 1})
    growth = cobra.Reaction("growth", lower_bound=0, upper_bound=float("inf"))
    growth.add_metabolites({substrate: -1})
    growth.gene_reaction_rule = "gene_a"
    model.add_reactions([supply, growth])
    model.objective = growth
    model.add_cons_vars([
        model.problem.Constraint(growth.flux_expression, ub=8, name="capacity"),
        model.problem.Variable("unbounded_probe"),
    ])
    model.add_groups([cobra.core.Group("pathway", members=[growth, substrate, model.genes[0]])])
    model.compartments = {"c": "cytosol"}
    model._sbml = {"fixture_isolation": {"values": ["source"]}}
    for owner in _owners(model):
        owner.annotation = {"fixture_isolation": {"values": ["source"]}}
        owner.notes = {"fixture_isolation": {"values": ["source"]}}
    return model


def _state(model):
    configuration = model.solver.configuration
    return {
        "bounds": {reaction.id: reaction.bounds for reaction in model.reactions},
        "stoichiometry": {
            reaction.id: {met.id: coefficient for met, coefficient in reaction.metabolites.items()}
            for reaction in model.reactions
        },
        "variables": {variable.name: (variable.lb, variable.ub, variable.type)
                      for variable in model.variables},
        "constraint_bounds": {constraint.name: (constraint.lb, constraint.ub)
                              for constraint in model.constraints},
        "constraint_coefficients": str(model.constraints[0].expression),
        "objective": (model.objective.direction, str(model.objective.expression)),
        "configuration": (configuration.presolve, configuration.timeout, configuration.verbosity,
                          configuration.tolerances.to_dict(), model.tolerance),
        "solver_status": model.solver.status,
        "annotations": [deepcopy(owner.annotation) for owner in _owners(model)],
        "notes": [deepcopy(owner.notes) for owner in _owners(model)],
        "compartments": model.compartments,
        "sbml": deepcopy(getattr(model, "_sbml", {})),
    }


def _mutate(model):
    reaction = next(reaction for reaction in model.reactions
                    if reaction.lower_bound < reaction.upper_bound)
    reaction.bounds = (0, 0)
    metabolite = next(iter(reaction.metabolites))
    reaction.add_metabolites({metabolite: 0.5})
    model.constraints[0].lb = -1
    model.constraints[0].set_linear_coefficients({reaction.forward_variable: 2})
    extra = model.problem.Variable("fixture_consumer_only", lb=0, ub=1)
    model.add_cons_vars([extra, model.problem.Constraint(extra, ub=0.5, name="fixture_extra_cap")])
    model.objective = model.problem.Objective(extra, direction="min")
    configuration = model.solver.configuration
    configuration.timeout = 11
    configuration.presolve = True
    model.tolerance *= 10
    for owner in _owners(model):
        for attribute in ("annotation", "notes"):
            getattr(owner, attribute).setdefault("fixture_isolation", {"values": []})[
                "values"].append("consumer")
    compartment = next(iter(model.compartments))
    model.compartments = {compartment: "consumer"}
    model._sbml.setdefault("fixture_isolation", {"values": []})["values"].append("consumer")


@pytest.fixture
def fixture_runner(pytester, monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join(filter(None, [str(repo / "src"), str(repo), os.getenv("PYTHONPATH")])))
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    pytester.makeini("[pytest]\nmarkers =\n    integration: touches the full GSMM\n    slow: long-running\n")
    return pytester


@pytest.mark.parametrize("fixture_name", ["yeast_gem", "ec_yeast_gem"])
@pytest.mark.parametrize("source", ["minimal", pytest.param("native", marks=pytest.mark.integration)])
def test_fixture_consumers_are_isolated_without_reparsing(fixture_runner, fixture_name, source):
    if source == "native" and getattr(paths, fixture_name)() is None:
        pytest.skip(f"{fixture_name} source is absent")
    fixture_runner.makeconftest(f'''
        from contextlib import nullcontext
        from pathlib import Path
        from unittest.mock import patch

        import cobra
        import pytest

        from tests.test_model_fixture_isolation import _minimal_model, _state
        from ystwin import paths

        pytest_plugins = ["tests.conftest"]

        @pytest.fixture(scope="session")
        def observations():
            return {{"reads": [], "consumers": []}}

        @pytest.fixture(scope="session", autouse=True)
        def instrument_source(observations):
            read = cobra.io.read_sbml_model
            def load(path):
                observations["reads"].append(path)
                model = _minimal_model() if {source!r} == "minimal" else read(path)
                observations["source"] = model
                observations["baseline"] = _state(model)
                return model
            resolver = (patch.object(paths, {fixture_name!r}, return_value=Path("minimal.xml"))
                        if {source!r} == "minimal" else nullcontext())
            with resolver, patch.object(cobra.io, "read_sbml_model", side_effect=load):
                yield
            assert len(observations["reads"]) == 1, "SBML must be parsed only once per session"
    ''')
    fixture_runner.makepyfile(test_consumers=f'''
        import pytest
        from tests.test_model_fixture_isolation import _mutate, _state

        @pytest.mark.parametrize("consumer", range(3))
        def test_consumer_gets_the_unmodified_source(request, observations, consumer):
            model = request.getfixturevalue({fixture_name!r})
            assert request.getfixturevalue({fixture_name!r}) is model
            actual = _state(model)
            changed = [key for key, value in observations["baseline"].items() if actual[key] != value]
            assert not changed, f"{fixture_name} inherited mutations in {{changed}}"
            for previous in observations["consumers"]:
                assert model is not previous
                assert model.solver is not previous.solver
                assert model.solver.problem is not previous.solver.problem
                assert model.solver.configuration is not previous.solver.configuration
                assert model.variables[0] is not previous.variables[0]
                assert model.constraints[0] is not previous.constraints[0]
            optimum = model.slim_optimize(error_value=None)
            assert model.solver.status == "optimal"
            assert optimum > 0
            if consumer == 0:
                observations["optimum"] = optimum
            else:
                assert optimum == pytest.approx(observations["optimum"], rel=1e-9)
            observations["consumers"].append(model)
            _mutate(model)
            assert _state(model) != observations["baseline"]
    ''')
    result = fixture_runner.runpytest_subprocess("-q", "--strict-markers", "-W", "error")
    result.assert_outcomes(passed=3)


@pytest.mark.integration
def test_stress_aerobic_fixture_does_not_change_the_next_native_consumer(fixture_runner):
    if paths.yeast_gem() is None:
        pytest.skip("yeast_gem source is absent")
    fixture_runner.makeconftest('''
        import pytest
        pytest_plugins = ["tests.conftest"]

        @pytest.fixture(scope="session")
        def baseline():
            return {}
    ''')
    fixture_runner.makepyfile(test_stress_consumers='''
        from tests.test_stress_energetics import TestTheClaimsAboutTheGEMAreTrue as _Claims

        def bounds(model):
            return {rid: model.reactions.get_by_id(rid).bounds for rid in ("r_1714", "r_1992")}

        class TestStressConsumer:
            aerobic = _Claims.aerobic

            def test_leaves_its_configuration_in_place(self, aerobic, yeast_gem, baseline):
                baseline["bounds"] = bounds(yeast_gem)
                baseline["model"] = yeast_gem
                configured = aerobic(-10.0)
                assert configured is yeast_gem
                assert configured.reactions.get_by_id("r_1714").lower_bound == -10.0
                assert configured.reactions.get_by_id("r_1992").lower_bound == -1000.0
                assert bounds(configured) != baseline["bounds"]

        def test_the_next_consumer_still_gets_the_native_bounds(yeast_gem, baseline):
            assert bounds(yeast_gem) == baseline["bounds"]
            assert yeast_gem is not baseline["model"]
            assert yeast_gem.solver.problem is not baseline["model"].solver.problem
    ''')
    result = fixture_runner.runpytest_subprocess("-q", "--strict-markers", "-W", "error")
    result.assert_outcomes(passed=2)


def test_copies_keep_the_metabolic_graph_internal_to_each_consumer():
    from tests.conftest import _copy_model

    template = _minimal_model()
    baseline = _state(template)
    first, second = _copy_model(template), _copy_model(template)
    for model in (first, second):
        assert all(owner._model is model for owner in _owners(model)[1:])
        assert model.groups.pathway.members == {
            model.reactions.growth, model.metabolites.substrate_c, model.genes.gene_a,
        }
        assert model.reactions.growth.gpr is not template.reactions.growth.gpr
        assert model.reactions.growth.gpr.body is not template.reactions.growth.gpr.body
        assert _state(model) == baseline
    first.genes.gene_a.knock_out()
    first.groups.pathway.remove_members([first.reactions.growth])
    assert first.reactions.growth.bounds == (0, 0)
    assert first.reactions.growth not in first.groups.pathway.members
    for untouched in (template, second):
        assert untouched.genes.gene_a.functional
        assert untouched.reactions.growth in untouched.groups.pathway.members
        assert _state(untouched) == baseline


_MIGRATED_FIXTURES = (
    ("test_carotenoid_ec", ("ec_carotenoid", "installed"), ()),
    ("test_fba_surrogate", ("surrogate",), ()),
    ("test_gem_environment", (), (("TestTheMagnitudeIsNotPredicted", ("kocharin",)),)),
    ("test_insulin_pathway", ("insulin_model",), ()),
    ("test_mech_burden", (), (("TestTheGemPricesTheLoadBesideTheMeasurement", ("sink_model",)),)),
    ("test_mech_chain", (), (("TestOneInputReachesManyReactions", ("acid_model",)),)),
    ("test_mech_oxidative", (), (("TestTheGemPricesIt", ("ros",)),)),
    ("test_stress_biomass", ("deletant_model", "thick_wall_model", "heat_shock_model", "zeroed_models"), ()),
    ("test_stress_genes", ("base_model", "fixed_model"), ()),
    ("test_stress_heterologous", ("phb_model", "phb_default", "carotenoid", "insulin_model"), ()),
    ("test_stress_oxygen", ("gem", "anaerobic"), (("TestAgainstVerduyn1990", ("report",)),)),
    ("test_stress_pathways", ("base_model", "restored_model"), ()),
    ("test_stress_ph", ("gem", "acid_models"), ()),
    ("test_stress_pools", ("pooled",), ()),
    ("test_stress_proteostasis", ("sink_model", "chaperone_model", "composition", "ec_carotenoid", "ec_charged"), ()),
    ("test_stress_ros", ("gem", "ros"), ()),
    ("test_thermo_gate", (), (
        ("TestTheTwoEstimatorsDisagreeFarBeyondTheirQuotedError", ("both",)),
        ("TestARefusedMagnitudeIsCheckedOnDGAndNotOnDG0", ("both",)),
        ("TestTheDisagreementCheckAppliesToEveryReaction", ("both",)),
    )),
)


def _fixture_value_state(value):
    from ystwin.bridge.equilibrator import EquilibratorData, PreferredEnergies

    if isinstance(value, cobra.Model):
        return _state(value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return "NaN" if math.isnan(value) else value
    if isinstance(value, dict):
        return {key: _fixture_value_state(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_fixture_value_state(item) for item in value]
    if isinstance(value, np.ndarray):
        return _fixture_value_state(value.tolist())
    if isinstance(value, np.generic):
        return _fixture_value_state(value.item())
    if isinstance(value, pd.DataFrame):
        return _fixture_value_state(value.to_dict(orient="split"))
    if isinstance(value, RegularGridInterpolator):
        return _fixture_value_state((value.grid, value.values, value.method, value.fill_value))
    if isinstance(value, PreferredEnergies):
        return _fixture_value_state((value.component, value.group))
    if isinstance(value, EquilibratorData):
        return _fixture_value_state({
            "conditions": (value.ph, value.p_mg, value.ionic_strength, value.temperature),
            "contribution_conditions": tuple(str(getattr(value.contribution, name))
                                             for name in ("p_h", "p_mg", "ionic_strength", "temperature")),
            "compounds": {mid: (compound.id, deepcopy(compound.atom_bag))
                          for mid, compound in value.compounds.items()},
            "mismatched": value.mismatched,
        })
    if is_dataclass(value):
        return {field.name: _fixture_value_state(getattr(value, field.name)) for field in fields(value)}
    raise TypeError(type(value))


def _mutate_fixture_value(value):
    from ystwin.bridge.equilibrator import EquilibratorData, PreferredEnergies
    from ystwin.bridge.thermodynamic import ThermodynamicData
    from ystwin.fba.stress_oxygen import AnaerobicReport
    from ystwin.fba.stress_proteostasis import ProteinComposition
    from ystwin.fba.surrogate import FluxSurrogate

    if isinstance(value, cobra.Model):
        _mutate(value)
    elif isinstance(value, dict):
        for item in value.values():
            _mutate_fixture_value(item)
    elif isinstance(value, tuple):
        for item in value:
            _mutate_fixture_value(item)
    elif isinstance(value, pd.DataFrame):
        column = value.select_dtypes(include="number").columns[0]
        value.iloc[0, value.columns.get_loc(column)] += 1
    elif isinstance(value, FluxSurrogate):
        value.axes[0][0] -= 1
        next(iter(value.interpolators.values())).values.flat[0] += 1
        value.inputs["consumer_only"] = (0, 1)
    elif isinstance(value, ProteinComposition):
        value.aa_coefficients[next(iter(value.aa_coefficients))] *= 2
    elif isinstance(value, ThermodynamicData):
        value.by_metabolite[next(iter(value.by_metabolite))] += 1
    elif isinstance(value, EquilibratorData):
        from equilibrator_api import Q_

        value.ph += 1
        value.contribution.p_h = Q_(value.ph)
        compound = next(compound for compound in value.compounds.values()
                        if compound.atom_bag is not None)
        compound.atom_bag["C"] = compound.atom_bag.get("C", 0) + 1
        value.mismatched["consumer_only"] = ("C", "C2")
    elif isinstance(value, PreferredEnergies):
        _mutate_fixture_value(value.component)
        _mutate_fixture_value(value.group)
    elif isinstance(value, AnaerobicReport):
        with pytest.raises(FrozenInstanceError):
            value.glucose_bound = 0
    else:
        raise TypeError(type(value))


@pytest.mark.integration
def test_every_migrated_fixture_caches_builds_but_not_mutable_consumers(fixture_runner, monkeypatch):
    fixture_runner.makeconftest('''
        from collections import Counter
        from unittest.mock import patch
        import pytest
        import tests.test_fba_surrogate as surrogate_tests
        from ystwin.fba import insulin

        pytest_plugins = ["tests.conftest"]
        setups = Counter()

        @pytest.hookimpl(hookwrapper=True)
        def pytest_fixture_setup(fixturedef, request):
            outcome = yield
            if outcome.excinfo is None:
                setups[(request.node.nodeid, fixturedef.argname)] += 1

        @pytest.fixture
        def setup_counts():
            return setups

        @pytest.fixture(scope="session", autouse=True)
        def expensive_build_counts():
            with patch.object(surrogate_tests, "build_surrogate", wraps=surrogate_tests.build_surrogate) as surrogate, \\
                 patch.object(insulin, "add_insulin_precursor_pathway", wraps=insulin.add_insulin_precursor_pathway) as protein:
                yield
                assert surrogate.call_count == 1
                assert protein.call_count == 2
    ''')
    expected = 0
    for index, (module_name, module_fixtures, classes) in enumerate(_MIGRATED_FIXTURES):
        source = [
            "import pytest",
            f"from tests import {module_name} as original",
            "from tests.test_model_fixture_isolation import _fixture_value_state, _mutate_fixture_value",
        ]
        for name in module_fixtures:
            source.extend((f"_{name}_template = original._{name}_template", f"{name} = original.{name}"))
        groups = ([(None, module_fixtures)] if module_fixtures else []) + list(classes)
        for group_index, (class_name, names) in enumerate(groups):
            source.extend(("", f"class TestConsumers{group_index}:"))
            if class_name:
                for name in names:
                    source.extend((
                        f"    _{name}_template = original.{class_name}._{name}_template",
                        f"    {name} = original.{class_name}.{name}",
                    ))
            source.extend((
                "    @pytest.fixture(scope='class')",
                "    def observations(self):",
                "        return {}",
                "",
                "    @pytest.mark.parametrize('consumer', range(2))",
                f"    @pytest.mark.parametrize('fixture_name', {names!r})",
                "    def test_independent_consumers(self, request, observations, fixture_name, consumer):",
                "        value = request.getfixturevalue(fixture_name)",
                "        state = _fixture_value_state(value)",
                "        if consumer == 0:",
                "            observations[fixture_name] = state",
                "        else:",
                "            assert state == observations[fixture_name]",
                "        _mutate_fixture_value(value)",
                "",
                "    def test_templates_remain_pristine_and_build_only_once(self, request, observations, setup_counts):",
                f"        for name in {names!r}:",
                "            value = request.getfixturevalue(name)",
                "            assert _fixture_value_state(value) == observations[name]",
                f"            owner = request.node.getparent(pytest.{'Class' if class_name else 'Module'}).nodeid",
                "            assert setup_counts[(owner, '_' + name + '_template')] == 1",
            ))
            expected += 2 * len(names) + 1
        fixture_runner.makepyfile(**{f"test_cached_{index:02d}": "\n".join(source)})
    monkeypatch.chdir(paths.REPO_ROOT)
    _restore_caller_environment(monkeypatch, _CALLER_ENVIRONMENT)
    result = fixture_runner.runpytest_subprocess(
        "-q", "--strict-markers", "-W", "error", "--tb=short",
        "-c", str(fixture_runner.path / "tox.ini"),
        "--confcutdir", str(fixture_runner.path), str(fixture_runner.path),
    )
    result.assert_outcomes(passed=expected)


def test_repeated_clones_keep_native_unbounded_columns_unbounded():
    from tests.conftest import _copy_model

    template = _minimal_model().copy()
    assert template.reactions.growth.forward_variable.ub == sys.float_info.max
    for _ in range(3):
        model = _copy_model(template)
        assert model.reactions.growth.forward_variable.ub is None
        assert model.variables.unbounded_probe.lb is None
        assert model.variables.unbounded_probe.ub is None
        template = model.copy()


def test_native_clones_preserve_finite_extremes_and_exact_lp_coefficients():
    from tests.conftest import _copy_model

    template = _minimal_model()
    finite = template.problem.Variable("finite_extreme", lb=-sys.float_info.max, ub=sys.float_info.max)
    template.add_cons_vars(finite)
    coefficient = 1.2345678901234567
    template.constraints.substrate_c.set_linear_coefficients({template.variables.supply: coefficient})
    template.objective = template.problem.Objective(
        coefficient * template.reactions.growth.flux_expression + 2.5, direction="max")
    model = _copy_model(template)
    assert model.variables.finite_extreme.lb == -sys.float_info.max
    assert model.variables.finite_extreme.ub == sys.float_info.max
    assert model.constraints.substrate_c.get_linear_coefficients(
        [model.variables.supply])[model.variables.supply] == coefficient
    assert str(model.objective.expression) == str(template.objective.expression)
    assert model.solver.configuration.tolerances.to_dict() == template.solver.configuration.tolerances.to_dict()


def test_model_copies_do_not_inherit_a_solved_basis():
    import swiglpk as glpk

    from tests.conftest import _copy_model

    template = _minimal_model()
    assert template.slim_optimize(error_value=None) == 8
    basis = [glpk.glp_get_col_stat(template.solver.problem, index)
             for index in range(1, len(template.variables) + 1)]
    assert glpk.GLP_BS in basis
    model = _copy_model(template)
    assert model.solver.status is None
    assert all(glpk.glp_get_col_stat(model.solver.problem, index) != glpk.GLP_BS
               for index in range(1, len(model.variables) + 1))
    assert basis == [glpk.glp_get_col_stat(template.solver.problem, index)
                     for index in range(1, len(template.variables) + 1)]
    assert model.slim_optimize(error_value=None) == 8


@pytest.mark.parametrize("provided", [
    pytest.param({}, id="all-absent"),
    pytest.param({"HOME": "declared-home"}, id="home-only"),
    pytest.param({"USERPROFILE": "declared-profile"}, id="profile-only"),
    pytest.param({"XDG_CACHE_HOME": "declared-cache"}, id="cache-only"),
    pytest.param({"HOME": "declared-home", "USERPROFILE": "declared-profile",
                  "XDG_CACHE_HOME": "declared-cache"}, id="distinct-values"),
    pytest.param({"HOME": "", "USERPROFILE": "", "XDG_CACHE_HOME": ""}, id="explicit-empty"),
])
def test_caller_environment_forwarding_preserves_values_and_absence(tmp_path, monkeypatch, provided):
    expected = {name: str(tmp_path / value) if value else "" for name, value in provided.items()}
    for name in _CALLER_ENVIRONMENT_NAMES:
        if name in expected:
            monkeypatch.setenv(name, expected[name])
        else:
            monkeypatch.delenv(name, raising=False)
    captured = _capture_caller_environment()
    for name in _CALLER_ENVIRONMENT_NAMES:
        monkeypatch.setenv(name, str(tmp_path / "pytester" / name))
    unrelated = str(tmp_path / "unrelated")
    monkeypatch.setenv("FIXTURE_ENV_UNRELATED", unrelated)

    _restore_caller_environment(monkeypatch, captured)

    assert captured == {name: expected.get(name) for name in _CALLER_ENVIRONMENT_NAMES}
    assert {name: os.environ[name] for name in _CALLER_ENVIRONMENT_NAMES if name in os.environ} == expected
    assert os.environ["FIXTURE_ENV_UNRELATED"] == unrelated
