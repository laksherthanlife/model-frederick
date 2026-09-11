from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import libsbml
import numpy as np
import pytest

from ystwin.mech import kinetic_sbml
from ystwin.mech.kinetic_sbml import (
    KineticModel,
    KineticSimulationError,
    UnsupportedSBMLError,
)


def _species(model, name, compartment="c", initial=4.0, *, amount=False):
    species = model.createSpecies()
    species.setId(name)
    species.setCompartment(compartment)
    species.setHasOnlySubstanceUnits(amount)
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setInitialConcentration(initial)
    return species


def _parameter(model, name, value=None, *, constant=True):
    parameter = model.createParameter()
    parameter.setId(name)
    parameter.setConstant(constant)
    if value is not None:
        parameter.setValue(value)
    return parameter


def _reaction(model, name, formula, reactants=(), products=()):
    reaction = model.createReaction()
    reaction.setId(name)
    reaction.setReversible(False)
    reaction.setFast(False)
    for names, create in (
        (reactants, reaction.createReactant), (products, reaction.createProduct)
    ):
        for species in names:
            reference = create()
            reference.setSpecies(species)
            reference.setStoichiometry(1.0)
            if model.getLevel() == 3:
                reference.setConstant(True)
    reaction.createKineticLaw().setMath(libsbml.parseL3Formula(formula))
    return reaction


def _document(level=2, version=1):
    document = libsbml.SBMLDocument(level, version)
    model = document.createModel()
    model.setId("analytic_decay")
    if level == 3:
        model.setTimeUnits("second")
        model.setSubstanceUnits("mole")
        model.setExtentUnits("mole")
        model.setVolumeUnits("litre")
    compartment = model.createCompartment()
    compartment.setId("c")
    compartment.setSize(2.0)
    compartment.setSpatialDimensions(3)
    compartment.setConstant(True)
    _species(model, "A")
    _parameter(model, "k", 0.5)
    _reaction(model, "decay", "c * k * A", ("A",))
    return document


def _load(tmp_path, document):
    path = tmp_path / "model.xml"
    assert libsbml.writeSBMLToFile(document, str(path)) == 1
    return KineticModel.from_sbml(path)


def _assignment(model, name, formula):
    if model.getParameter(name) is None:
        _parameter(model, name, constant=False)
    rule = model.createAssignmentRule()
    rule.setVariable(name)
    rule.setMath(libsbml.parseL3Formula(formula))
    return rule


def _function(model, name, formula):
    function = model.createFunctionDefinition()
    function.setId(name)
    function.setMath(libsbml.parseL3Formula(formula))
    return function


def test_analytic_decay_outputs_rates_observables_and_source_provenance(tmp_path):
    model = _load(tmp_path, _document())
    times = np.linspace(0, 8, 41)
    trajectory = model.simulate(times)
    expected = 4 * np.exp(-0.5 * times)

    assert model.model_id == "analytic_decay"
    assert model.species_ids == ("A",)
    assert model.reaction_ids == ("decay",)
    assert model.parameter_values == {"k": 0.5}
    np.testing.assert_allclose(trajectory.states[:, 0], expected, rtol=3e-6)
    np.testing.assert_array_equal(trajectory.times_s, times)
    np.testing.assert_array_equal(trajectory.variables["A"], trajectory.states[:, 0])
    np.testing.assert_allclose(trajectory.reaction_rates[:, 0], expected, rtol=3e-6)
    np.testing.assert_array_equal(trajectory.variables["k"], np.full(times.shape, 0.5))
    assert model.metadata["sha256"] == hashlib.sha256((tmp_path / "model.xml").read_bytes()).hexdigest()
    assert model.metadata["parameter_provenance"] == "imported SBML definitions; not learned here"
    assert model.metadata["learned_parameters"] == {}
    assert model.metadata["initial_state_time_s"] == 0.0
    assert "constant_compartments" in model.metadata["supported_features"]
    assert trajectory.metadata["solver"]["method"] == "BDF"


def test_reversible_transfer_conserves_amount_across_different_volumes(tmp_path):
    document = _document()
    source = document.getModel()
    source.removeReaction(0)
    compartment = source.createCompartment()
    compartment.setId("outside")
    compartment.setSize(5.0)
    compartment.setConstant(True)
    _species(source, "B", "outside", initial=0.0)
    _parameter(source, "reverse", 0.25)
    _reaction(source, "transfer", "c * k * A - outside * reverse * B", ("A",), ("B",))
    model = _load(tmp_path, document)
    times = np.linspace(0, 10, 51)
    trajectory = model.simulate(times)
    amount_a = 8 * (1 / 3 + 2 / 3 * np.exp(-0.75 * times))

    np.testing.assert_allclose(trajectory.states[:, 0], amount_a / 2, rtol=2e-6)
    np.testing.assert_allclose(trajectory.states[:, 1], (8 - amount_a) / 5, atol=2e-7)
    np.testing.assert_allclose(trajectory.states @ [2, 5], 8, atol=1e-12)
    np.testing.assert_allclose(model.rhs(0, model.initial_state()), [-2, 0.8])


def test_amount_and_concentration_initialization_use_native_compartment_size(tmp_path):
    document = _document()
    source = document.getModel()
    unit_definition = source.createUnitDefinition()
    unit_definition.setId("volume")
    unit = unit_definition.createUnit()
    unit.setKind(libsbml.UNIT_KIND_LITRE)
    unit.setExponent(1)
    unit.setScale(-3)
    unit.setMultiplier(1.0)
    source.getSpecies("A").setInitialAmount(6.0)
    _species(source, "amount_pool", initial=4.0, amount=True)
    _reaction(source, "amount_loss", "2 * k * amount_pool", ("amount_pool",))
    model = _load(tmp_path, document)

    np.testing.assert_array_equal(model.initial_state(), [3, 8])
    np.testing.assert_allclose(model.rhs(0, model.initial_state()), [-1.5, -8])
    assert model.evaluate(0, model.initial_state())["c"] == 2
    assert model.metadata["units"]["definitions"]["volume"][0]["scale"] == -3
    assert model.metadata["units"]["numeric_conversion"] == "none; model-native arithmetic"
    assert model.metadata["units"]["species"]["A"]["representation"] == "concentration"
    assert model.metadata["units"]["species"]["amount_pool"]["representation"] == "amount"


def test_boundary_and_constant_species_are_held_without_changing_rates(tmp_path):
    document = _document()
    source = document.getModel()
    source.getSpecies("A").setBoundaryCondition(True)
    species = _species(source, "fixed", initial=7.0)
    species.setConstant(True)
    species.setBoundaryCondition(True)
    _reaction(source, "fixed_loss", "2 * k * fixed", ("fixed",))
    model = _load(tmp_path, document)

    np.testing.assert_array_equal(model.rhs(0, model.initial_state()), [0, 0])
    np.testing.assert_array_equal(model.reaction_rates(0, model.initial_state()), [4, 7])
    np.testing.assert_array_equal(model.simulate([0, 1]).states, [[4, 7], [4, 7]])


@pytest.mark.parametrize("level,version", [(2, 1), (2, 4), (3, 1), (3, 2)])
def test_local_parameters_shadow_globals_only_in_their_reaction(tmp_path, level, version):
    document = _document(level, version)
    source = document.getModel()
    law = source.getReaction(0).getKineticLaw()
    local = law.createLocalParameter() if level == 3 else law.createParameter()
    local.setId("k")
    local.setValue(3.0)
    _reaction(source, "global_rate", "k * A", products=("A",))
    model = _load(tmp_path, document)

    np.testing.assert_allclose(model.reaction_rates(0, [4], {"k": 2}), [24, 8])
    assert model.evaluate(0, [4], {"k": 2})["k"] == 2
    assert model.parameter_values == {"k": 0.5}


def test_functions_and_assignments_are_topologically_ordered_and_lexically_scoped(tmp_path):
    document = _document(2, 4)
    source = document.getModel()
    _parameter(source, "time_scale", 1.0)
    _function(source, "outer", "lambda(k, inner(k) + 1)")
    _function(source, "inner", "lambda(A, A * 2)")
    _assignment(source, "last", "outer(first) + k")
    _assignment(source, "first", "A + time * time_scale")
    source.getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula("last"))
    model = _load(tmp_path, document)
    values = model.evaluate(3, [4], {"k": 0.25})

    assert values["first"] == 7
    assert values["last"] == 15.25
    assert values["A"] == 4
    assert values["k"] == 0.25
    np.testing.assert_allclose(model.reaction_rates(3, [4], {"k": 0.25}), [15.25])
    trajectory = model.simulate([0, 0.01])
    np.testing.assert_allclose(trajectory.variables["first"], trajectory.states[:, 0] + trajectory.times_s)
    assert model.metadata["assignment_order"] == ["first", "last"]
    assert model.metadata["function_order"] == ["inner", "outer"]
    assert set(model.parameter_values) == {"k", "time_scale"}


@pytest.mark.parametrize(
    "formula,expected",
    [
        ("2 + 3 * 4 - 5 / 2", 11.5),
        ("-A + pow(2, 3)", 4),
        ("exp(ln(3)) + log(100) + log(2, 8)", 8),
        ("abs(-3) + floor(1.2) + ceiling(1.1)", 6),
        ("pi + exponentiale", np.pi + np.e),
        ("root(3, -8) + sqrt(9)", 1),
        ("piecewise(7, and(A >= 4, A <= 4), 1 / 0)", 7),
        ("piecewise(7, or(A < 4, A == 4), 1 / 0)", 7),
        ("piecewise(7, not(A != 4), 1 / 0)", 7),
        ("piecewise(7, xor(A == 4, A > 4), 1 / 0)", 7),
        ("piecewise(1 / 0, false, 7, true, ln(-1))", 7),
    ],
)
def test_whitelisted_math_and_lazy_piecewise(tmp_path, formula, expected):
    document = _document()
    document.getModel().getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula(f"k * ({formula})"))
    model = _load(tmp_path, document)

    assert model.reaction_rates(0, [4], {"k": 1})[0] == pytest.approx(expected)


def test_piecewise_knots_resolve_a_narrow_input_and_observables_use_exact_times(tmp_path):
    document = _document()
    source = document.getModel()
    source.removeReaction(0)
    source.getSpecies("A").setInitialConcentration(0)
    _assignment(source, "input", "piecewise(0, time < 10, 1000, time < 10.001, 0)")
    _reaction(source, "pulse", "c * input", products=("A",))
    model = _load(tmp_path, document)
    times = np.array([0, 10, 10.0005, 10.001, 20])
    trajectory = model.simulate(times, breakpoints=(10, 10.001))

    np.testing.assert_allclose(trajectory.states[:, 0], [0, 0, 0.5, 1, 1], atol=2e-7)
    np.testing.assert_array_equal(trajectory.variables["input"], [0, 1000, 1000, 0, 0])
    assert len(trajectory.metadata["solver"]["segments"]) == 3


def test_caller_arrays_mappings_and_exported_metadata_do_not_mutate_the_model(tmp_path):
    model = _load(tmp_path, _document())
    times = np.array([0.0, 1.0])
    state = np.array([2.0])
    parameters = {"k": 0.25}
    overrides = {"A": 3.0}
    model.initial_state(overrides)
    model.rhs(0, state, parameters)
    model.evaluate(0, state, parameters)
    model.reaction_rates(0, state, parameters)
    trajectory = model.simulate(times, parameters=parameters, initial_state=state)
    model.parameter_values["k"] = 99
    model.metadata["units"]["numeric_conversion"] = "corrupted"

    np.testing.assert_array_equal(times, [0, 1])
    np.testing.assert_array_equal(state, [2])
    assert parameters == {"k": 0.25}
    assert overrides == {"A": 3}
    assert model.parameter_values == {"k": 0.5}
    assert model.metadata["units"]["numeric_conversion"] == "none; model-native arithmetic"
    assert trajectory.metadata["parameter_overrides"] == {"k": 0.25}
    assert trajectory.metadata["learned_parameters"] == {}
    times[:] = 7
    state[:] = 7
    parameters["k"] = 7
    np.testing.assert_array_equal(trajectory.times_s, [0, 1])
    assert trajectory.states[0, 0] == 2
    assert trajectory.metadata["parameter_overrides"] == {"k": 0.25}
    np.testing.assert_array_equal(model.initial_state(), [4])


def test_initial_state_overrides_single_time_and_negative_states_are_not_clipped(tmp_path):
    document = _document()
    document.getModel().getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula("2 * k * c"))
    model = _load(tmp_path, document)

    np.testing.assert_array_equal(model.initial_state({"A": 2}), [2])
    np.testing.assert_array_equal(model.simulate([0], initial_state=[2]).states, [[2]])
    assert model.simulate([0, 2], initial_state=[1]).states[-1, 0] == pytest.approx(-1)
    np.testing.assert_array_equal(model.rhs(0, [-2]), [-1])


@pytest.mark.parametrize("parameters", [{"unknown": 1}, {"k": np.nan}, {"k": np.inf}])
@pytest.mark.parametrize("method", ["rhs", "evaluate", "reaction_rates", "simulate"])
def test_parameter_overrides_are_partial_and_validated(tmp_path, parameters, method):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="parameter|finite"):
        if method == "simulate":
            model.simulate([0, 1], parameters=parameters)
        else:
            getattr(model, method)(0, [4], parameters)


def test_assignment_targets_cannot_be_overridden(tmp_path):
    document = _document()
    _assignment(document.getModel(), "derived", "2 * k")
    model = _load(tmp_path, document)

    with pytest.raises(ValueError, match="assignment"):
        model.simulate([0, 1], parameters={"derived": 4})


@pytest.mark.parametrize("overrides", [{"missing": 1}, {"A": np.nan}])
def test_initial_state_overrides_are_validated(tmp_path, overrides):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="species|finite"):
        model.initial_state(overrides)


@pytest.mark.parametrize("state", [[1, 2], [[1]], [np.nan], [np.inf]])
def test_invalid_state_vectors_are_refused(tmp_path, state):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="state"):
        model.rhs(0, state)
    with pytest.raises(ValueError, match="state"):
        model.simulate([0, 1], initial_state=state)


@pytest.mark.parametrize("times", [[], [1, 2], [0, 0], [0, -1], [0, np.nan], [[0, 1]]])
def test_time_grid_requires_initial_time_zero_and_strict_increase(tmp_path, times):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="times_s"):
        model.simulate(times)


@pytest.mark.parametrize("knots", [(np.nan,), (1, 0.5), (1, 1), (-1,)])
def test_invalid_breakpoints_are_refused(tmp_path, knots):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="breakpoints"):
        model.simulate([0, 2], breakpoints=knots)


@pytest.mark.parametrize("options", [{"rtol": 0}, {"atol": np.nan}, {"method": "invented"}])
def test_invalid_solver_settings_are_refused(tmp_path, options):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError):
        model.simulate([0, 1], **options)


def test_solver_failure_is_not_retried_or_hidden(tmp_path, monkeypatch):
    model = _load(tmp_path, _document())
    calls = []

    def failed_solver(*args, **kwargs):
        calls.append(kwargs["method"])
        return SimpleNamespace(success=False, message="forced solver failure")

    monkeypatch.setattr(kinetic_sbml, "solve_ivp", failed_solver)
    with pytest.raises(KineticSimulationError, match="forced solver failure"):
        model.simulate([0, 1])
    assert calls == ["BDF"]


def test_integration_budget_refuses_a_run_instead_of_hanging(tmp_path, monkeypatch):
    model = _load(tmp_path, _document())
    monkeypatch.setattr(kinetic_sbml, "MAX_RHS_EVALUATIONS", 1)
    with pytest.raises(KineticSimulationError, match="evaluation budget"):
        model.simulate([0, 1])


@pytest.mark.parametrize("formula", ["1 / 0", "ln(-1)", "exp(1000)", "piecewise(1, false)"])
def test_undefined_numeric_values_raise_with_reaction_context(tmp_path, formula):
    document = _document()
    document.getModel().getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula(formula))
    model = _load(tmp_path, document)
    with pytest.raises(KineticSimulationError, match="decay"):
        model.reaction_rates(0, [4])


@pytest.mark.parametrize("formula", ["missing * A", "unknown(A)", "sin(A)", "delay(A, 1)"])
def test_unknown_symbols_and_unsupported_math_are_rejected_at_import(tmp_path, formula):
    document = _document()
    document.getModel().getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula(formula))
    with pytest.raises(UnsupportedSBMLError, match="unknown|unsupported|delay"):
        _load(tmp_path, document)


def test_assignment_cycles_are_rejected(tmp_path):
    document = _document()
    source = document.getModel()
    _assignment(source, "first", "last + 1")
    _assignment(source, "last", "first + 1")
    with pytest.raises(UnsupportedSBMLError, match="cycle"):
        _load(tmp_path, document)


@pytest.mark.parametrize("second_body", ["lambda(x, first(x))", "lambda(x, x + k)"])
def test_function_cycles_and_free_global_capture_are_rejected(tmp_path, second_body):
    document = _document()
    source = document.getModel()
    _function(source, "first", "lambda(x, second(x))")
    _function(source, "second", second_body)
    with pytest.raises(UnsupportedSBMLError, match="cycle|unknown|scope"):
        _load(tmp_path, document)


def test_function_arity_is_validated(tmp_path):
    document = _document()
    source = document.getModel()
    _function(source, "f", "lambda(x, x + 1)")
    source.getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula("f(A, k)"))
    with pytest.raises(UnsupportedSBMLError, match="argument|arity"):
        _load(tmp_path, document)


@pytest.mark.parametrize("feature", ["rate", "algebraic", "event", "initial_assignment", "variable_compartment", "stoichiometry_math", "species_assignment", "fast_reaction", "constraint"])
def test_unsupported_model_features_fail_closed(tmp_path, feature):
    document = _document(2, 4)
    source = document.getModel()
    if feature == "rate":
        source.getParameter("k").setConstant(False)
        rule = source.createRateRule()
        rule.setVariable("k")
        rule.setMath(libsbml.parseL3Formula("1"))
    elif feature == "algebraic":
        source.createAlgebraicRule().setMath(libsbml.parseL3Formula("A - 1"))
    elif feature == "event":
        event = source.createEvent()
        event.setId("event")
        event.createTrigger().setMath(libsbml.parseL3Formula("time > 1"))
        assignment = event.createEventAssignment()
        assignment.setVariable("A")
        assignment.setMath(libsbml.parseL3Formula("2"))
    elif feature == "initial_assignment":
        assignment = source.createInitialAssignment()
        assignment.setSymbol("A")
        assignment.setMath(libsbml.parseL3Formula("2"))
    elif feature == "variable_compartment":
        source.getCompartment(0).setConstant(False)
    elif feature == "stoichiometry_math":
        source.getReaction(0).getReactant(0).createStoichiometryMath().setMath(libsbml.parseL3Formula("A"))
    elif feature == "species_assignment":
        rule = source.createAssignmentRule()
        rule.setVariable("A")
        rule.setMath(libsbml.parseL3Formula("2"))
    elif feature == "fast_reaction":
        source.getReaction(0).setFast(True)
    elif feature == "constraint":
        source.createConstraint().setMath(libsbml.parseL3Formula("A > 0"))
    with pytest.raises(UnsupportedSBMLError):
        _load(tmp_path, document)


@pytest.mark.parametrize("feature", ["model_conversion", "species_conversion", "variable_stoichiometry"])
def test_level_three_conversion_factors_and_variable_stoichiometry_are_refused(tmp_path, feature):
    document = _document(3, 2)
    source = document.getModel()
    if feature == "model_conversion":
        source.setConversionFactor("k")
    elif feature == "species_conversion":
        source.getSpecies(0).setConversionFactor("k")
    else:
        source.getReaction(0).getReactant(0).setConstant(False)
    with pytest.raises(UnsupportedSBMLError):
        _load(tmp_path, document)


@pytest.mark.parametrize("feature", ["zero_volume", "nonfinite_parameter", "missing_initial", "missing_kinetics", "nonsecond_time"])
def test_undefined_or_incompatible_source_values_are_refused(tmp_path, feature):
    document = _document()
    source = document.getModel()
    if feature == "zero_volume":
        source.getCompartment(0).setSize(0)
    elif feature == "nonfinite_parameter":
        source.getParameter("k").setValue(np.inf)
    elif feature == "missing_initial":
        source.getSpecies(0).unsetInitialConcentration()
    elif feature == "missing_kinetics":
        source.getReaction(0).unsetKineticLaw()
    elif feature == "nonsecond_time":
        definition = source.createUnitDefinition()
        definition.setId("time")
        unit = definition.createUnit()
        unit.setKind(libsbml.UNIT_KIND_SECOND)
        unit.setExponent(1)
        unit.setMultiplier(60)
    with pytest.raises(UnsupportedSBMLError):
        _load(tmp_path, document)


def test_malformed_xml_is_not_accepted_as_an_empty_model(tmp_path):
    path = tmp_path / "invalid.xml"
    path.write_text("<sbml><model>")
    with pytest.raises(UnsupportedSBMLError):
        KineticModel.from_sbml(path)


def test_small_native_compartment_size_does_not_overflow_a_finite_derivative(tmp_path):
    document = _document()
    source = document.getModel()
    source.getCompartment(0).setSize(1e-300)
    source.getSpecies(0).setInitialConcentration(1)
    source.getParameter("k").setValue(1e-155)
    source.getReaction(0).getReactant(0).setStoichiometry(1e10)
    source.getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula("k * k * A"))
    model = _load(tmp_path, document)

    np.testing.assert_allclose(model.rhs(0, [1]), [-1])


@pytest.mark.parametrize("feature", ["local_value", "compartment_size", "global_value"])
def test_missing_source_values_are_not_replaced_by_libsbml_storage_defaults(tmp_path, feature):
    document = _document()
    source = document.getModel()
    if feature == "local_value":
        parameter = source.getReaction(0).getKineticLaw().createParameter()
        parameter.setId("local")
    elif feature == "compartment_size":
        source.getCompartment(0).unsetSize()
    else:
        source.getParameter("k").unsetValue()
    with pytest.raises(UnsupportedSBMLError, match="value|size"):
        _load(tmp_path, document)


def test_complex_caller_values_are_refused_without_discarding_imaginary_parts(tmp_path):
    model = _load(tmp_path, _document())
    with pytest.raises(ValueError, match="real|state"):
        model.rhs(0, np.array([1 + 2j]))
    with pytest.raises(ValueError, match="real|parameter"):
        model.rhs(0, [4], {"k": np.complex128(1 + 2j)})
    with pytest.raises(ValueError, match="real|times_s"):
        model.simulate(np.array([0j, 1 + 2j]))


def test_wall_time_budget_is_enforced(tmp_path, monkeypatch):
    model = _load(tmp_path, _document())
    monkeypatch.setattr(kinetic_sbml, "MAX_WALL_SECONDS", 0)
    with pytest.raises(KineticSimulationError, match="wall-time budget"):
        model.simulate([0, 1])


def test_real_extension_packages_are_refused(tmp_path):
    document = _document(3, 1)
    document.enablePackage(libsbml.FbcExtension.getXmlnsL3V1V2(), "fbc", True)
    document.setPackageRequired("fbc", False)
    document.getModel().getPlugin("fbc").setStrict(False)
    with pytest.raises(UnsupportedSBMLError, match="package"):
        _load(tmp_path, document)


def test_knots_at_the_final_time_use_the_left_limit_for_integration_only(tmp_path, monkeypatch):
    document = _document()
    source = document.getModel()
    source.getReaction(0).getKineticLaw().setMath(libsbml.parseL3Formula("piecewise(k * A, time < 1, 100 * k * A)"))
    model = _load(tmp_path, document)
    seen = []
    solver = kinetic_sbml.solve_ivp

    def inspect_endpoint(fun, interval, state, **kwargs):
        seen.append(fun(interval[1], state)[0])
        return solver(fun, interval, state, **kwargs)

    monkeypatch.setattr(kinetic_sbml, "solve_ivp", inspect_endpoint)
    trajectory = model.simulate([0, 1], breakpoints=[1, 2])
    assert seen == pytest.approx([-1])
    assert trajectory.reaction_rates[-1, 0] == pytest.approx(50 * trajectory.states[-1, 0])
    assert trajectory.metadata["requested_breakpoints_s"] == [1, 2]
    assert trajectory.metadata["breakpoints_s"] == []


def test_published_hog_source_native_values_simple_rate_and_short_integration():
    path = Path(__file__).resolve().parents[1] / "data" / "hog2013" / "model_wt.xml"
    if not path.exists():
        pytest.skip("parent-managed Petelenz-Kurdziel 2013 WT SBML artifact is absent")
    model = KineticModel.from_sbml(path)
    document = libsbml.readSBMLFromFile(str(path))
    source = document.getModel()
    initial = model.initial_state()
    assert len(model.species_ids) == 29
    assert len(model.reaction_ids) == 58
    np.testing.assert_array_equal(initial, [species.getInitialConcentration() for species in source.getListOfSpecies()])
    state = model.initial_state({"Hog1PP": 1.25})
    values = model.evaluate(0, state)
    rates = model.reaction_rates(0, state)
    assert rates[model.reaction_ids.index("v16r")] == pytest.approx(values["intra"] * values["kv16r_1"] * 1.25)
    isolated_v16r_rhs = np.zeros(len(initial))
    isolated_v16r_rhs[model.species_ids.index("Hog1")] = values["kv16r_1"] * 1.25
    isolated_v16r_rhs[model.species_ids.index("Hog1PP")] = -values["kv16r_1"] * 1.25
    np.testing.assert_allclose(model.rhs(0, state) - model.rhs(0, state, {"kv16r_1": 0}), isolated_v16r_rhs, atol=1e-12)
    expected_rhs = np.zeros(len(initial))
    for reaction, rate in zip(source.getListOfReactions(), rates, strict=True):
        for references, sign in ((reaction.getListOfReactants(), -1), (reaction.getListOfProducts(), 1)):
            for reference in references:
                species = source.getSpecies(reference.getSpecies())
                if species.getBoundaryCondition() or species.getConstant():
                    continue
                volume = source.getCompartment(species.getCompartment()).getSize()
                expected_rhs[model.species_ids.index(species.getId())] += sign * reference.getStoichiometry() * rate / volume
    np.testing.assert_allclose(model.rhs(0, state), expected_rhs, atol=1e-12)
    trajectory = model.simulate([0, 0.01, 0.1, 1])
    assert np.isfinite(trajectory.states).all()
    assert np.isfinite(trajectory.reaction_rates).all()
    assert "cellvol" in trajectory.variables
    assert all(np.isfinite(values).all() for values in trajectory.variables.values())
    for measured, concentration in (("Hog1PP_measured", "Hog1PP"), ("Gpd1_measured", "Gpd1"), ("glycerol_measured", "glycerol_i")):
        np.testing.assert_allclose(trajectory.variables[measured], trajectory.variables[concentration] * trajectory.variables["cellvol"] / trajectory.variables["cellvol_init"])
    assert values["Hog1PP_measured"] > 1
    assert trajectory.metadata["state_clipping"] == "none"


def test_published_hog_trajectory_agrees_between_implicit_solvers():
    path = Path(__file__).resolve().parents[1] / "data" / "hog2013" / "model_wt.xml"
    model = KineticModel.from_sbml(path)
    times = np.unique(np.r_[np.linspace(0, 14400, 145), 3605.0])
    options = {"breakpoints": (3600, 3605), "rtol": 1e-8, "atol": 1e-11}
    bdf = model.simulate(times, method="BDF", **options).states
    radau = model.simulate(times, method="Radau", **options).states
    peak = np.maximum(np.max(np.abs(radau), axis=0), 1e-20)
    assert np.max(np.abs(bdf - radau) / peak) < 1e-5


def test_published_hog_reaction_rates_match_independent_roadrunner():
    roadrunner = pytest.importorskip("roadrunner")
    path = Path(__file__).resolve().parents[1] / "data" / "hog2013" / "model_wt.xml"
    model = KineticModel.from_sbml(path)
    reference = roadrunner.RoadRunner(str(path))
    times = np.array([0.0, 3600.0, 3605.0, 3720.0, 5400.0, 9000.0, 14400.0])
    trajectory = model.simulate(times, breakpoints=(3600, 3605))
    reaction_ids = list(reference.model.getReactionIds())
    indices = [reaction_ids.index(name) for name in model.reaction_ids]
    for time_s, state in zip(times, trajectory.states):
        reference.model.setTime(float(time_s))
        assert reference.getValue("time") == time_s
        for species, concentration in zip(model.species_ids, state):
            reference[f"[{species}]"] = float(concentration)
        expected = np.asarray(reference.getReactionRates())[indices]
        np.testing.assert_allclose(model.reaction_rates(time_s, state), expected, rtol=1e-11, atol=1e-12)
