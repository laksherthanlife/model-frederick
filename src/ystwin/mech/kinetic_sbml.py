from __future__ import annotations

import hashlib
import math
import operator
import time
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import libsbml
import numpy as np
from scipy.integrate import solve_ivp

MAX_RHS_EVALUATIONS = 200_000
MAX_WALL_SECONDS = 30.0


class UnsupportedSBMLError(ValueError):
    pass


class KineticSimulationError(RuntimeError):
    pass


@dataclass(frozen=True)
class Trajectory:
    times_s: np.ndarray
    states: np.ndarray
    variables: dict[str, np.ndarray]
    reaction_rates: np.ndarray
    metadata: dict


def _finite(value, label, error=ValueError):
    if isinstance(value, (complex, np.complexfloating, np.ndarray)):
        raise error(f"{label} must be a finite real scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise error(f"{label} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise error(f"{label} must be finite, got {result}")
    return result


def _real_array(value, label):
    try:
        array = np.asarray(value)
        if array.dtype.kind not in "biuf":
            raise ValueError(f"{label} must contain only real numbers")
        with np.errstate(over="ignore", invalid="ignore"):
            result = np.array(array, dtype=float, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must contain finite real numbers") from exc
    if not np.isfinite(result).all():
        raise ValueError(f"{label} must contain finite real numbers")
    return result


def _named(elements, label):
    result = {}
    for element in elements:
        name = element.getId()
        if not name or name in result:
            raise UnsupportedSBMLError(f"missing or duplicate {label} identifier: {name!r}")
        result[name] = element
    return result


def _references(node, node_type):
    if node is None:
        raise UnsupportedSBMLError("missing mathematical expression")
    names = {node.getName()} if node.getType() == node_type else set()
    for index in range(node.getNumChildren()):
        names.update(_references(node.getChild(index), node_type))
    return names


def _topological_order(graph, label):
    ordered, active, finished = [], [], set()

    def visit(name):
        if name in active:
            raise UnsupportedSBMLError(f"{label} cycle: {' -> '.join(active + [name])}")
        if name in finished:
            return
        if name not in graph:
            raise UnsupportedSBMLError(f"unknown {label}: {name}")
        active.append(name)
        for dependency in sorted(graph[name]):
            visit(dependency)
        active.pop()
        finished.add(name)
        ordered.append(name)

    for name in graph:
        visit(name)
    return ordered


def _real_root(degree, value):
    if value < 0 and degree.is_integer() and int(degree) % 2:
        return -math.pow(-value, 1 / degree)
    return math.pow(value, 1 / degree)


_UNARY_MATH = {
    libsbml.AST_FUNCTION_EXP: math.exp,
    libsbml.AST_FUNCTION_LN: math.log,
    libsbml.AST_FUNCTION_ABS: abs,
    libsbml.AST_FUNCTION_FLOOR: math.floor,
    libsbml.AST_FUNCTION_CEILING: math.ceil,
}
_BINARY_MATH = {
    libsbml.AST_DIVIDE: operator.truediv,
    libsbml.AST_POWER: math.pow,
    libsbml.AST_FUNCTION_POWER: math.pow,
}
_COMPARISONS = {
    libsbml.AST_RELATIONAL_EQ: operator.eq,
    libsbml.AST_RELATIONAL_NEQ: operator.ne,
    libsbml.AST_RELATIONAL_LT: operator.lt,
    libsbml.AST_RELATIONAL_LEQ: operator.le,
    libsbml.AST_RELATIONAL_GT: operator.gt,
    libsbml.AST_RELATIONAL_GEQ: operator.ge,
}
_CONSTANTS = {
    libsbml.AST_CONSTANT_PI: math.pi,
    libsbml.AST_CONSTANT_E: math.e,
    libsbml.AST_CONSTANT_TRUE: True,
    libsbml.AST_CONSTANT_FALSE: False,
}


def _compile(node, symbols, functions, *, constants=None, time_index=None, context="math"):
    if node is None:
        raise UnsupportedSBMLError(f"missing math in {context}")
    kind = node.getType()
    count = node.getNumChildren()
    if kind in (libsbml.AST_INTEGER, libsbml.AST_REAL, libsbml.AST_REAL_E, libsbml.AST_RATIONAL):
        if count:
            raise UnsupportedSBMLError(f"invalid number in {context}")
        value = _finite(node.getValue(), context, UnsupportedSBMLError)
        return lambda values: value
    if kind in _CONSTANTS:
        value = _CONSTANTS[kind]
        return lambda values: value
    if kind == libsbml.AST_NAME:
        name = node.getName()
        if constants is not None and name in constants:
            value = constants[name]
            return lambda values: value
        if name not in symbols:
            raise UnsupportedSBMLError(f"unknown symbol {name!r} in {context}; function scope contains only its arguments")
        index = symbols[name]
        return lambda values: values[index]
    if kind == libsbml.AST_NAME_TIME:
        if time_index is None:
            raise UnsupportedSBMLError(f"time must be passed as a function argument in {context}")
        return lambda values: values[time_index]
    if kind == libsbml.AST_FUNCTION:
        name = node.getName()
        if name not in functions:
            raise UnsupportedSBMLError(f"unknown function {name!r} in {context}")
        arity, body = functions[name]
        if count != arity:
            raise UnsupportedSBMLError(f"function {name!r} requires {arity} arguments, got {count} in {context}")
    elif not node.hasCorrectNumberArguments():
        raise UnsupportedSBMLError(f"invalid math arity in {context}: {libsbml.formulaToL3String(node)}")
    arguments = tuple(
        _compile(node.getChild(index), symbols, functions, constants=constants,
                 time_index=time_index, context=context)
        for index in range(count)
    )
    if kind == libsbml.AST_FUNCTION:
        return lambda values: body(tuple(argument(values) for argument in arguments))
    if kind == libsbml.AST_PLUS and count >= 2:
        return lambda values: sum(argument(values) for argument in arguments)
    if kind == libsbml.AST_TIMES and count >= 2:
        return lambda values: math.prod(argument(values) for argument in arguments)
    if kind == libsbml.AST_MINUS and count in (1, 2):
        if count == 1:
            return lambda values: -arguments[0](values)
        return lambda values: arguments[0](values) - arguments[1](values)
    if kind in _BINARY_MATH and count == 2:
        operation = _BINARY_MATH[kind]
        return lambda values: operation(arguments[0](values), arguments[1](values))
    if kind in _UNARY_MATH and count == 1:
        operation = _UNARY_MATH[kind]
        return lambda values: operation(arguments[0](values))
    if kind == libsbml.AST_FUNCTION_LOG and count in (1, 2):
        if count == 1:
            return lambda values: math.log10(arguments[0](values))
        return lambda values: math.log(arguments[1](values), arguments[0](values))
    if kind == libsbml.AST_FUNCTION_ROOT and count in (1, 2):
        if count == 1:
            return lambda values: math.sqrt(arguments[0](values))
        return lambda values: _real_root(float(arguments[0](values)), arguments[1](values))
    if kind in _COMPARISONS and count >= 2:
        operation = _COMPARISONS[kind]

        def compare(values):
            evaluated = tuple(argument(values) for argument in arguments)
            return all(operation(left, right) for left, right in zip(evaluated, evaluated[1:]))

        return compare
    if kind == libsbml.AST_LOGICAL_AND:
        return lambda values: all(argument(values) for argument in arguments)
    if kind == libsbml.AST_LOGICAL_OR:
        return lambda values: any(argument(values) for argument in arguments)
    if kind == libsbml.AST_LOGICAL_XOR:
        return lambda values: sum(bool(argument(values)) for argument in arguments) % 2 == 1
    if kind == libsbml.AST_LOGICAL_NOT and count == 1:
        return lambda values: not arguments[0](values)
    if kind == libsbml.AST_FUNCTION_PIECEWISE and count >= 1:

        def piecewise(values):
            for index in range(0, count - 1, 2):
                if arguments[index + 1](values):
                    return arguments[index](values)
            if count % 2:
                return arguments[-1](values)
            raise ValueError("piecewise has no matching branch and no otherwise value")

        return piecewise
    raise UnsupportedSBMLError(f"unsupported math node {kind} in {context}: {libsbml.formulaToL3String(node)}")


def _unit_metadata(source):
    definitions = {}
    for name, definition in _named(source.getListOfUnitDefinitions(), "unit").items():
        units = []
        for unit in definition.getListOfUnits():
            values = {
                "kind": libsbml.UnitKind_toString(unit.getKind()),
                "exponent": _finite(unit.getExponentAsDouble(), f"unit {name} exponent", UnsupportedSBMLError),
                "scale": unit.getScale(),
                "multiplier": _finite(unit.getMultiplier(), f"unit {name} multiplier", UnsupportedSBMLError),
                "offset": _finite(unit.getOffset(), f"unit {name} offset", UnsupportedSBMLError),
            }
            if values["offset"] != 0 or values["multiplier"] <= 0:
                raise UnsupportedSBMLError(f"unsupported offset or nonpositive multiplier in unit {name}")
            units.append(values)
        definitions[name] = units
    time_unit = "time" if source.getLevel() == 2 else source.getTimeUnits()
    if time_unit in definitions:
        units = definitions[time_unit]
        try:
            seconds = (len(units) == 1 and units[0]["kind"] == "second"
                       and units[0]["exponent"] == 1
                       and units[0]["multiplier"] * 10.0 ** units[0]["scale"] == 1)
        except OverflowError:
            seconds = False
    else:
        seconds = time_unit == "second" or (source.getLevel() == 2 and time_unit == "time")
    if not seconds:
        raise UnsupportedSBMLError("times_s requires model-native time units equal to one second; time conversion is unsupported")
    return {
        "definitions": definitions,
        "time": time_unit,
        "time_seconds_per_native_unit": 1.0,
        "substance": source.getSubstanceUnits() or ("substance" if source.getLevel() == 2 else "unspecified"),
        "volume": source.getVolumeUnits() or ("volume" if source.getLevel() == 2 else "unspecified"),
        "extent": source.getExtentUnits() or ("substance" if source.getLevel() == 2 else "unspecified"),
        "level_two_defaults": {"time": "second", "substance": "mole", "volume": "litre"},
        "numeric_conversion": "none; model-native arithmetic",
        "interpretation": "declared/default SBML units, not experimentally established mass-specific units; no gDW conversion",
        "compartments": {}, "species": {}, "parameters": {},
    }


def _check_errors(document):
    errors = [document.getError(index) for index in range(document.getNumErrors())
              if document.getError(index).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
    if errors:
        messages = "; ".join(f"{error.getErrorId()}: {error.getMessage().strip()}" for error in errors[:5])
        raise UnsupportedSBMLError(f"invalid SBML: {messages}")


class KineticModel:
    @classmethod
    def from_sbml(cls, path: str | Path) -> KineticModel:
        source_path = Path(path)
        payload = source_path.read_bytes()
        try:
            document = libsbml.readSBMLFromString(payload.decode("utf-8-sig"))
        except UnicodeDecodeError as exc:
            raise UnsupportedSBMLError("only UTF-8 SBML XML is supported") from exc
        _check_errors(document)
        return cls(document, source_path, hashlib.sha256(payload).hexdigest())

    def __init__(self, document, source_path, sha256):
        source = document.getModel()
        if source is None:
            raise UnsupportedSBMLError("SBML document has no model")
        if (document.getLevel(), document.getVersion()) not in {
            (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2)}:
            raise UnsupportedSBMLError("only SBML Level 2 Versions 1-5 and Level 3 Versions 1-2 core are supported")
        plugins = [document.getPlugin(index).getPackageName() for index in range(document.getNumPlugins())
                   if document.getPlugin(index).getURI() != document.getSBMLNamespaces().getURI()]
        empty_legacy_plugins = (document.getLevel() == 2 and set(plugins) <= {"layout", "render"}
                                and source.getListOfAllElementsFromPlugins().getSize() == 0)
        if document.getNumUnknownPackages() or (plugins and not empty_legacy_plugins):
            raise UnsupportedSBMLError("SBML extension packages are unsupported")
        for label, count in (
            ("events", source.getNumEvents()),
            ("initial assignments", source.getNumInitialAssignments()),
            ("constraints", source.getNumConstraints()),
            ("compartment types", source.getNumCompartmentTypes()),
            ("species types", source.getNumSpeciesTypes()),
        ):
            if count:
                raise UnsupportedSBMLError(f"unsupported SBML feature: {label}")
        if source.isSetConversionFactor():
            raise UnsupportedSBMLError("model conversion factors are unsupported")
        compartments = _named(source.getListOfCompartments(), "compartment")
        parameters = _named(source.getListOfParameters(), "parameter")
        species = _named(source.getListOfSpecies(), "species")
        reactions = _named(source.getListOfReactions(), "reaction")
        definitions = _named(source.getListOfFunctionDefinitions(), "function")
        all_names = [name for group in (compartments, parameters, species, reactions, definitions) for name in group]
        if len(set(all_names)) != len(all_names):
            raise UnsupportedSBMLError("duplicate identifiers in the global SBML namespace")
        if not compartments or not species:
            raise UnsupportedSBMLError("at least one compartment and species are required")
        rules = {}
        for rule in source.getListOfRules():
            if not rule.isAssignment():
                raise UnsupportedSBMLError("rate and algebraic rules are unsupported; only parameter assignment rules are supported")
            target = rule.getVariable()
            if target not in parameters or parameters[target].getConstant():
                raise UnsupportedSBMLError(f"assignment target {target!r} must be a nonconstant global parameter")
            if target in rules:
                raise UnsupportedSBMLError(f"duplicate assignment target {target!r}")
            rules[target] = rule.getMath()
        self.model_id = source.getId()
        self.species_ids = tuple(species)
        self.reaction_ids = tuple(reactions)
        self._symbol_ids = tuple(compartments) + tuple(parameters) + self.species_ids
        self._indices = {name: index for index, name in enumerate(self._symbol_ids)}
        self._time_index = len(self._symbol_ids)
        self._species_slice = slice(len(compartments) + len(parameters), self._time_index)
        self._base_values = [0.0] * (self._time_index + 1)
        self._parameter_values = {}
        self._assignment_targets = frozenset(rules)
        units = _unit_metadata(source)
        for name, compartment in compartments.items():
            if not compartment.getConstant() or compartment.getSpatialDimensions() != 3:
                raise UnsupportedSBMLError(f"compartment {name!r} must be constant and three-dimensional")
            if not compartment.isSetSize():
                raise UnsupportedSBMLError(f"compartment {name!r} requires an explicit size")
            size = _finite(compartment.getSize(), f"compartment {name} size", UnsupportedSBMLError)
            if size <= 0:
                raise UnsupportedSBMLError(f"compartment {name!r} size must be positive")
            self._base_values[self._indices[name]] = size
            units["compartments"][name] = {"size": size, "units": compartment.getUnits() or units["volume"]}
        parameter_definitions = {}
        for name, parameter in parameters.items():
            value = (_finite(parameter.getValue(), f"parameter {name}", UnsupportedSBMLError)
                     if parameter.isSetValue() else None)
            if name not in rules:
                if value is None:
                    raise UnsupportedSBMLError(f"parameter {name!r} has no value or assignment")
                self._parameter_values[name] = value
                self._base_values[self._indices[name]] = value
            parameter_definitions[name] = {"value": value, "constant": parameter.getConstant(),
                                           "assignment_target": name in rules}
            units["parameters"][name] = parameter.getUnits() or "unspecified"
        initial, divisors, held = [], [], []
        for name, item in species.items():
            compartment = item.getCompartment()
            if compartment not in compartments:
                raise UnsupportedSBMLError(f"unknown compartment {compartment!r} for species {name!r}")
            if item.isSetConversionFactor() or item.isSetSpatialSizeUnits():
                raise UnsupportedSBMLError(f"species conversion factors and spatialSizeUnits are unsupported: {name}")
            size = self._base_values[self._indices[compartment]]
            amount = item.getHasOnlySubstanceUnits()
            if item.isSetInitialAmount():
                value = item.getInitialAmount() if amount else item.getInitialAmount() / size
            elif item.isSetInitialConcentration():
                value = item.getInitialConcentration() * size if amount else item.getInitialConcentration()
            else:
                raise UnsupportedSBMLError(f"species {name!r} requires an explicit initial amount or concentration")
            initial.append(_finite(value, f"initial state {name}", UnsupportedSBMLError))
            divisors.append(1.0 if amount else size)
            held.append(item.getBoundaryCondition() or item.getConstant())
            units["species"][name] = {
                "representation": "amount" if amount else "concentration",
                "substance_units": item.getSubstanceUnits() or units["substance"],
                "compartment": compartment, "compartment_units": units["compartments"][compartment]["units"],
                "boundary_condition": item.getBoundaryCondition(), "constant": item.getConstant(),
            }
        self._initial = np.array(initial, dtype=float)
        functions = {}
        function_order = _topological_order(
            {name: _references(definition.getBody(), libsbml.AST_FUNCTION)
             for name, definition in definitions.items()}, "function")
        for name in function_order:
            definition = definitions[name]
            if definition.getMath() is None or definition.getMath().getType() != libsbml.AST_LAMBDA:
                raise UnsupportedSBMLError(f"function {name!r} requires a lambda")
            arguments = [definition.getArgument(index) for index in range(definition.getNumArguments())]
            names = [argument.getName() for argument in arguments]
            if any(argument.getType() != libsbml.AST_NAME for argument in arguments) or len(set(names)) != len(names):
                raise UnsupportedSBMLError(f"invalid arguments in function {name!r}")
            body = _compile(definition.getBody(), {argument: index for index, argument in enumerate(names)},
                            functions, context=f"function {name}")
            functions[name] = (len(names), body)
        assignment_order = _topological_order(
            {name: _references(node, libsbml.AST_NAME) & rules.keys() for name, node in rules.items()},
            "assignment")
        self._assignments = [
            (name, self._indices[name], _compile(rules[name], self._indices, functions,
                                               time_index=self._time_index, context=f"assignment {name}"))
            for name in assignment_order
        ]
        self._rate_functions = []
        local_parameters = {}
        stoichiometry = np.zeros((len(species), len(reactions)))
        species_indices = {name: index for index, name in enumerate(self.species_ids)}
        for column, (name, reaction) in enumerate(reactions.items()):
            if reaction.getFast():
                raise UnsupportedSBMLError(f"fast reaction {name!r} is unsupported")
            law = reaction.getKineticLaw()
            if law is None or law.getMath() is None:
                raise UnsupportedSBMLError(f"reaction {name!r} requires a kinetic law")
            if law.isSetTimeUnits() or law.isSetSubstanceUnits():
                raise UnsupportedSBMLError(f"kinetic-law-specific time/substance units are unsupported: {name}")
            locals_list = law.getListOfParameters() if source.getLevel() == 2 else law.getListOfLocalParameters()
            local = {}
            for key, parameter in _named(locals_list, f"local parameter in {name}").items():
                if not parameter.isSetValue():
                    raise UnsupportedSBMLError(f"local parameter {name}.{key} requires an explicit value")
                local[key] = _finite(parameter.getValue(), f"local parameter {name}.{key}", UnsupportedSBMLError)
            local_parameters[name] = local
            self._rate_functions.append(_compile(law.getMath(), self._indices, functions, constants=local,
                                                 time_index=self._time_index, context=f"reaction {name}"))
            for references, sign in ((reaction.getListOfReactants(), -1), (reaction.getListOfProducts(), 1)):
                for reference in references:
                    target = reference.getSpecies()
                    if target not in species_indices:
                        raise UnsupportedSBMLError(f"unknown species {target!r} in reaction {name!r}")
                    if reference.isSetStoichiometryMath() or (source.getLevel() == 3 and not reference.getConstant()):
                        raise UnsupportedSBMLError(f"variable stoichiometry is unsupported in reaction {name!r}")
                    if reference.getDenominator() != 1:
                        raise UnsupportedSBMLError(f"stoichiometry denominators are unsupported in reaction {name!r}")
                    coefficient = _finite(reference.getStoichiometry(), f"stoichiometry in {name}", UnsupportedSBMLError)
                    if coefficient < 0:
                        raise UnsupportedSBMLError(f"negative stoichiometry in reaction {name!r}")
                    row = species_indices[target]
                    stoichiometry[row, column] = _finite(float(stoichiometry[row, column]) + sign * coefficient,
                                                        f"net stoichiometry in {name}", UnsupportedSBMLError)
            for modifier in reaction.getListOfModifiers():
                if modifier.getSpecies() not in species:
                    raise UnsupportedSBMLError(f"unknown modifier species in reaction {name!r}")
        self._stoichiometry = stoichiometry
        self._stoichiometry[np.asarray(held), :] = 0
        self._species_divisors = np.asarray(divisors)
        document.checkConsistency()
        _check_errors(document)
        diagnostics = {}
        for index in range(document.getNumErrors()):
            diagnostic = document.getError(index)
            key = diagnostic.getErrorId()
            if key not in diagnostics:
                diagnostics[key] = {"id": key, "severity": diagnostic.getSeverityAsString(),
                                    "count": 0, "example": diagnostic.getMessage().strip()}
            diagnostics[key]["count"] += 1
        self._metadata = {
            "model_id": self.model_id, "source_path": str(source_path.resolve()), "sha256": sha256,
            "sbml_level": document.getLevel(), "sbml_version": document.getVersion(),
            "libsbml_version": libsbml.getLibSBMLDottedVersion(),
            "species_ids": list(self.species_ids), "reaction_ids": list(self.reaction_ids),
            "supported_features": ["constant_compartments", "constant_stoichiometry", "kinetic_laws",
                                   "local_parameters", "acyclic_parameter_assignments", "lexical_lambda_functions",
                                   "boundary_and_constant_species", "amount_and_concentration_states",
                                   "arithmetic", "real_powers_and_roots", "exp_log_abs_floor_ceiling",
                                   "logical_comparisons", "lazy_piecewise", "pi_e_time"],
            "unsupported_features": ["events", "rate_rules", "algebraic_rules", "delays", "initial_assignments",
                                     "conversion_factors", "variable_compartments", "variable_stoichiometry",
                                     "non_parameter_assignments", "fast_reactions", "constraints", "SBML_packages",
                                     "reaction_rate_symbols", "non_3D_compartments", "non_second_time_units",
                                     "kinetic_law_unit_overrides", "spatialSizeUnits"],
            "units": units, "sbml_diagnostics": list(diagnostics.values()),
            "empty_legacy_plugins": plugins,
            "parameter_provenance": "imported SBML definitions; not learned here",
            "parameter_values_scope": "independent global parameters only; assignment targets are evaluated",
            "parameter_definitions": parameter_definitions, "local_parameters": local_parameters,
            "learned_parameters": {}, "assignment_order": assignment_order, "function_order": function_order,
            "assignment_formulas": {name: libsbml.formulaToL3String(node) for name, node in rules.items()},
            "reaction_formulas": {name: libsbml.formulaToL3String(reaction.getKineticLaw().getMath())
                                  for name, reaction in reactions.items()},
            "initial_state_time_s": 0.0, "time_origin": "times_s must start at zero; initial_state applies at zero",
            "state_equations": "constant-compartment concentration derivatives = stoichiometry @ reaction_rates / native compartment size; amount derivatives omit the divisor; boundary/constant derivatives are zero",
            "dilution": "only explicit SBML reactions; no additional volume or growth dilution",
            "rate_convention": "kinetic laws evaluated exactly as written, including any omitted compartment factors",
            "discontinuities": "caller supplies all integration breakpoints; no automatic knot detection",
            "state_clipping": "none", "solver_fallback": "none",
        }

    @property
    def parameter_values(self) -> dict[str, float]:
        return self._parameter_values.copy()

    @property
    def metadata(self) -> dict:
        return deepcopy(self._metadata)

    def initial_state(self, overrides: Mapping[str, float] | None = None) -> np.ndarray:
        result = self._initial.copy()
        if overrides is not None:
            if not isinstance(overrides, Mapping):
                raise ValueError("initial species overrides must be a mapping")
            for name, value in overrides.items():
                if name not in self.species_ids:
                    raise ValueError(f"unknown initial species {name!r}")
                result[self.species_ids.index(name)] = _finite(value, f"initial species {name}")
        return result

    def _parameters(self, parameters):
        base = self._base_values.copy()
        overrides = {}
        if parameters is not None:
            if not isinstance(parameters, Mapping):
                raise ValueError("parameter overrides must be a mapping")
            for name, value in parameters.items():
                if name in self._assignment_targets:
                    raise ValueError(f"cannot override assignment target parameter {name!r}")
                if name not in self._parameter_values:
                    raise ValueError(f"unknown independent global parameter {name!r}")
                overrides[name] = _finite(value, f"parameter {name}")
                base[self._indices[name]] = overrides[name]
        return base, overrides

    def _state(self, state):
        values = _real_array(state, "state")
        if values.shape != self._initial.shape:
            raise ValueError(f"state must be finite with shape {self._initial.shape}")
        return values

    def _context(self, time_s, state, base):
        values = base.copy()
        values[self._species_slice] = state.tolist()
        values[self._time_index] = time_s
        for name, index, expression in self._assignments:
            values[index] = self._value(expression, values, f"assignment {name}", time_s)
        return values

    @staticmethod
    def _value(expression, values, label, time_s):
        try:
            return _finite(expression(values), label, KineticSimulationError)
        except (ArithmeticError, ValueError, TypeError, KineticSimulationError) as exc:
            raise KineticSimulationError(f"{label} failed at time_s={time_s}: {exc}") from exc

    def _rates(self, time_s, values):
        return np.array([self._value(expression, values, f"reaction {name}", time_s)
                         for name, expression in zip(self.reaction_ids, self._rate_functions, strict=True)])

    def evaluate(self, time_s: float, state, parameters: Mapping[str, float] | None = None) -> dict[str, float]:
        values = self._context(_finite(time_s, "time_s"), self._state(state), self._parameters(parameters)[0])
        return dict(zip(self._symbol_ids, values[:self._time_index], strict=True))

    def reaction_rates(self, time_s: float, state, parameters: Mapping[str, float] | None = None) -> np.ndarray:
        time_s = _finite(time_s, "time_s")
        values = self._context(time_s, self._state(state), self._parameters(parameters)[0])
        return self._rates(time_s, values)

    def _derivative(self, time_s, state, base):
        if not np.isfinite(state).all():
            raise KineticSimulationError(f"nonfinite solver state at time_s={time_s}")
        values = self._context(time_s, state, base)
        with np.errstate(over="ignore", invalid="ignore"):
            derivative = (self._stoichiometry @ self._rates(time_s, values)) / self._species_divisors
        if not np.isfinite(derivative).all():
            raise KineticSimulationError(f"nonfinite derivative at time_s={time_s}")
        return derivative

    def rhs(self, time_s: float, state, parameters: Mapping[str, float] | None = None) -> np.ndarray:
        return self._derivative(_finite(time_s, "time_s"), self._state(state), self._parameters(parameters)[0])

    def simulate(self, times_s: Sequence[float], *, parameters: Mapping[str, float] | None = None,
                 initial_state=None, breakpoints=(), rtol=1e-7, atol=1e-10, method="BDF") -> Trajectory:
        times = _real_array(times_s, "times_s")
        if (times.ndim != 1 or not times.size or times[0] != 0
                or not np.all(times[1:] > times[:-1])):
            raise ValueError("times_s must be a finite, nonempty, strictly increasing vector starting at zero")
        knots = _real_array(breakpoints, "breakpoints")
        if knots.ndim != 1 or np.any(knots < 0) or not np.all(knots[1:] > knots[:-1]):
            raise ValueError("breakpoints must be finite, nonnegative, and strictly increasing")
        used_knots = knots[(knots > 0) & (knots < times[-1])]
        rtol = _finite(rtol, "rtol")
        atol = _real_array(atol, "atol")
        if rtol < 100 * np.finfo(float).eps or np.any(atol <= 0):
            raise ValueError("rtol must be at least 100 machine epsilons and atol must be finite and positive")
        if atol.shape not in ((), self._initial.shape):
            raise ValueError("atol must be scalar or have the species-vector shape")
        if not isinstance(method, str) or method not in {"BDF", "Radau", "LSODA", "RK45", "RK23", "DOP853"}:
            raise ValueError(f"unsupported scipy solver method {method!r}")
        base, overrides = self._parameters(parameters)
        state = self.initial_state() if initial_state is None else self._state(initial_state)
        states = np.empty((times.size, len(state)))
        states[0] = state
        started = time.monotonic()
        calls = 0
        segments = []

        def check_budget():
            if time.monotonic() - started > MAX_WALL_SECONDS:
                raise KineticSimulationError(f"simulation exceeded wall-time budget of {MAX_WALL_SECONDS:g} seconds")

        edges = np.concatenate(([0.0], used_knots, [times[-1]]))
        for left, right in zip(edges, edges[1:]):
            if right == left:
                continue
            left_is_knot, right_is_knot = left in knots, right in knots

            def derivative(time_s, values):
                nonlocal calls
                calls += 1
                if calls > MAX_RHS_EVALUATIONS:
                    raise KineticSimulationError(f"simulation exceeded RHS evaluation budget of {MAX_RHS_EVALUATIONS}")
                check_budget()
                if left_is_knot and time_s == left:
                    time_s = float(np.nextafter(left, right))
                if right_is_knot and time_s == right:
                    time_s = float(np.nextafter(right, left))
                return self._derivative(time_s, values, base)

            requested = np.flatnonzero((times > left) & (times <= right))
            sample_times = times[requested]
            if not sample_times.size or sample_times[-1] != right:
                sample_times = np.append(sample_times, right)
            try:
                solution = solve_ivp(derivative, (left, right), state, t_eval=sample_times,
                                     rtol=rtol, atol=atol, method=method)
            except KineticSimulationError:
                raise
            except (ArithmeticError, ValueError) as exc:
                raise KineticSimulationError(f"{method} failed on [{left}, {right}]: {exc}") from exc
            if not solution.success:
                raise KineticSimulationError(f"{method} failed on [{left}, {right}]: {solution.message}")
            if not np.isfinite(solution.y).all():
                raise KineticSimulationError(f"{method} returned nonfinite states on [{left}, {right}]")
            states[requested] = solution.y[:, :requested.size].T
            state = solution.y[:, -1].copy()
            segments.append({"start_s": float(left), "end_s": float(right), "nfev": solution.nfev,
                             "njev": solution.njev, "nlu": solution.nlu})
        variables = np.empty((times.size, self._time_index))
        rates = np.empty((times.size, len(self.reaction_ids)))
        for index, (time_s, state) in enumerate(zip(times, states, strict=True)):
            check_budget()
            values = self._context(float(time_s), state, base)
            variables[index] = values[:self._time_index]
            rates[index] = self._rates(float(time_s), values)
        metadata = self.metadata
        metadata.update({
            "parameter_overrides": overrides, "initial_state_values": states[0].tolist(),
            "requested_breakpoints_s": knots.tolist(), "breakpoints_s": used_knots.tolist(),
            "breakpoint_evaluation": "one-sided integration at supplied knots, including endpoints; observables evaluated at exact requested times",
            "solver": {"method": method, "rtol": rtol, "atol": atol.tolist(), "segments": segments,
                       "rhs_evaluations": calls, "elapsed_s": time.monotonic() - started,
                       "max_rhs_evaluations": MAX_RHS_EVALUATIONS, "max_wall_seconds": MAX_WALL_SECONDS},
        })
        return Trajectory(times, states, {name: variables[:, index].copy()
                                         for index, name in enumerate(self._symbol_ids)}, rates, metadata)
