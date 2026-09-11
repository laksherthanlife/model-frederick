"""Explicit state allocation, not a learned biological mechanism by construction.

The public COSMIC execution kernel supplies precedent for state-specific uptake,
objective floors, environmental classification, and blending. Normalized CHO tables
have been inspected locally, but the complete fitted execution inputs remain missing;
data/cosmic_sources.json distinguishes those tables from an exact source reproduction.
State kinetics, enzyme capacities, task fractions, and protein budgets always have
caller-declared sources. Only the environmental state distribution is fitted, using
training rows alone. No protein-budget calibration or second stress penalty is applied.

All states solve the same stoichiometric model. Complete flux vectors are blended on
a dry-mass basis; cell-number fractions require explicit state dry masses. A mixture
need not satisfy any one state's capacities, but must satisfy common stoichiometry.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field, fields, is_dataclass
import hashlib
import json
import math
from numbers import Integral, Real
from types import MappingProxyType

import cobra
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

from ..analysis.partial_orders import ForwardEvaluation, Source
from ..bridge.regulation import (
    InfeasibleRegulation, MetabolicTask, _constraint_name, _expression, _set_floor,
    _solve_status, constrain_uptake, exchange_coefficients, specific_rate,
)
from .solver import growth_or_none


__all__ = [
    "Source", "ModelIdentity", "Environment", "Feature", "StateObservation",
    "TrainedStateDistribution", "fit_state_distribution", "FixedDistribution", "FixedTransition",
    "UptakePolicy", "FluxCapacity", "ProteinBudget", "ObjectiveTask", "StatePolicy",
    "MechanisticPolicy", "MixturePolicy", "FixedObjectivePolicy", "StateResult", "AllocationResult",
    "require_matched_models", "validate_host", "solve_allocation", "reactor_derivative",
]


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _number(value, name, minimum=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    value = float(value)
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _source(value):
    if not isinstance(value, Source):
        raise TypeError("every policy or observation needs a declared Source")


def _report(value):
    if is_dataclass(value):
        return {row.name: _report(getattr(value, row.name)) for row in fields(value)}
    if isinstance(value, Mapping):
        return {key: _report(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (tuple, list)):
        return [_report(item) for item in value]
    return value


def _fractions(values, names=None):
    if not isinstance(values, Mapping) or not values:
        raise ValueError("state fractions must be a nonempty mapping")
    if names is not None and set(values) != set(names):
        raise ValueError("state fractions must name every state exactly once")
    output = {_text(key, "state name"): _number(value, "state fraction", 0.0)
              for key, value in values.items()}
    if abs(math.fsum(output.values()) - 1.0) > 1e-8:
        raise ValueError("state fractions must sum to one; no implicit renormalization")
    return output


@dataclass(frozen=True)
class ModelIdentity:
    model_id: str
    organism: str
    strain: str
    network: str
    source: Source
    taxonomy_id: str | None = None

    def __post_init__(self):
        for name in ("model_id", "organism", "strain", "network"):
            _text(getattr(self, name), name)
        _source(self.source)


@dataclass(frozen=True)
class Environment:
    condition_id: str
    organism: str
    strain: str
    values: Mapping[str, float]
    units: Mapping[str, str]
    source: Source

    def __post_init__(self):
        for name in ("condition_id", "organism", "strain"):
            _text(getattr(self, name), name)
        _source(self.source)
        if set(self.values) != set(self.units):
            raise ValueError("environment values and units must name the same features")
        values = {_text(key, "environment feature"): _number(value, key)
                  for key, value in self.values.items()}
        units = {key: _text(value, "environment unit") for key, value in self.units.items()}
        object.__setattr__(self, "values", MappingProxyType(values))
        object.__setattr__(self, "units", MappingProxyType(units))


@dataclass(frozen=True)
class Feature:
    name: str
    unit: str
    role: str = "environment"

    def __post_init__(self):
        _text(self.name, "feature name")
        _text(self.unit, "feature unit")
        if self.role not in {"environment", "growth", "time"}:
            raise ValueError("feature role must be environment, growth, or time")

    def read(self, environment: Environment) -> float:
        if self.name not in environment.values:
            raise ValueError(f"missing environmental feature {self.name!r}")
        if environment.units[self.name] != self.unit:
            raise ValueError(f"feature {self.name!r}: expected unit {self.unit!r}")
        return environment.values[self.name]


def validate_host(model, identity: ModelIdentity, environment: Environment, *, strain_transfer=None):
    if identity.model_id != model.id:
        raise ValueError("host identity belongs to a different model_id")
    annotation = model.annotation.get("taxonomy")
    if identity.taxonomy_id is not None and annotation is not None:
        observed = {str(value) for value in annotation} if isinstance(annotation, list) else {str(annotation)}
        if identity.taxonomy_id not in observed:
            raise ValueError("declared host taxonomy contradicts the model annotation")
    if identity.organism != environment.organism:
        raise ValueError("incompatible host organism: no cross-organism measurement transfer")
    if identity.strain != environment.strain and not strain_transfer:
        raise ValueError("model and measured strain differ; an explicit strain_transfer reason is required")
    if strain_transfer is not None:
        _text(strain_transfer, "strain transfer reason")


def require_matched_models(first: ModelIdentity, second: ModelIdentity):
    if first.network != second.network:
        raise ValueError("plain/ec comparison requires the same base network, not different GEM versions")
    if (first.organism, first.strain) != (second.organism, second.strain):
        raise ValueError("plain/ec comparison requires matched host and strain declarations")


@dataclass(frozen=True)
class StateObservation:
    sample_id: str
    group_id: str
    split: str
    environment: Environment
    fractions: Mapping[str, float] | None
    source: Source

    def __post_init__(self):
        _text(self.sample_id, "sample id")
        _text(self.group_id, "group id")
        _source(self.source)
        if self.split not in {"train", "validation", "test"}:
            raise ValueError("split must be train, validation, or test")
        if self.fractions is not None:
            object.__setattr__(self, "fractions", MappingProxyType(dict(self.fractions)))


@dataclass(frozen=True)
class TrainedStateDistribution:
    states: tuple[str, ...]
    features: tuple[Feature, ...]
    basis: str
    center: np.ndarray
    scale: np.ndarray
    projection: np.ndarray
    coefficients: np.ndarray
    training_minimum: np.ndarray
    training_maximum: np.ndarray
    training_ids: tuple[str, ...]
    training_groups: tuple[str, ...]
    training_sources: tuple[Source, ...]
    training_digest: str
    organism: str
    strain: str
    ridge: float
    projection_dimension: int | None
    training_loss: float
    optimizer_iterations: int

    def __post_init__(self):
        for name in ("center", "scale", "projection", "coefficients", "training_minimum", "training_maximum"):
            value = np.array(getattr(self, name), dtype=float, copy=True)
            if not np.isfinite(value).all():
                raise ValueError("nonfinite fitted distribution parameter")
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    def _vector(self, environment):
        if (environment.organism, environment.strain) != (self.organism, self.strain):
            raise ValueError("state classifier host/strain differs from its training population")
        return np.array([feature.read(environment) for feature in self.features])

    def predict(self, environment: Environment) -> dict[str, float]:
        vector = self._vector(environment)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            design = np.r_[1.0, ((vector - self.center) / self.scale) @ self.projection]
            logits = design @ self.coefficients
            probabilities = np.exp(logits - logsumexp(logits))
        return _fractions(dict(zip(self.states, probabilities)), self.states)

    def extrapolated_features(self, environment: Environment) -> tuple[str, ...]:
        values = self._vector(environment)
        return tuple(feature.name for index, feature in enumerate(self.features)
                     if values[index] < self.training_minimum[index] or values[index] > self.training_maximum[index])

    def report(self):
        result = _report(self)
        result.update(
            method="training-only standardized optional PCA plus regularized multinomial soft-label log loss",
            parameter_selection="projection dimension and ridge explicitly supplied; no held-out optimization",
            biological_validation=False,
            label_source_kinds=sorted({source.kind for source in self.training_sources}),
            claim="statistical fit to declared labels, not learned cellular causes or objectives",
        )
        return result


def fit_state_distribution(
    observations: Sequence[StateObservation], *, features: Sequence[Feature], states: Sequence[str],
    projection_dimension: int | None, ridge: float, basis: str = "biomass",
) -> TrainedStateDistribution:
    features, states, observations = tuple(features), tuple(states), tuple(observations)
    if basis not in {"biomass", "cell"}:
        raise ValueError("fraction basis must be biomass or cell")
    if len(states) < 2 or len(set(states)) != len(states):
        raise ValueError("training needs at least two distinct named states")
    for state in states:
        _text(state, "state name")
    if not features or len({feature.name for feature in features}) != len(features):
        raise ValueError("declare distinct input features")
    if not any(feature.role == "environment" for feature in features):
        raise ValueError("a state distribution needs declared environmental features, not growth alone")
    ridge = _number(ridge, "ridge", 0.0)
    if len({row.sample_id for row in observations}) != len(observations):
        raise ValueError("duplicate observation sample id")
    training = sorted((row for row in observations if row.split == "train"), key=lambda row: row.sample_id)
    if len(training) < 2:
        raise ValueError("at least two training rows are required")
    train_groups = {row.group_id for row in training}
    if train_groups & {row.group_id for row in observations if row.split != "train"}:
        raise ValueError("a group crosses training and held-out partitions")
    hosts = {(row.environment.organism, row.environment.strain) for row in training}
    if len(hosts) != 1:
        raise ValueError("training labels must belong to one explicit host and strain")
    x = np.array([[feature.read(row.environment) for feature in features] for row in training])
    y = np.array([[values[name] for name in states]
                  for row in training for values in [_fractions(row.fractions, states)]])
    if np.any(y.sum(axis=0) <= 0):
        raise ValueError("each trained state needs positive training label support")
    center, scale = x.mean(axis=0), x.std(axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    standardized = (x - center) / scale
    if not np.isfinite(standardized).all():
        raise ValueError("nonfinite training transform")
    if projection_dimension is None:
        projection = np.eye(len(features))
    else:
        if (isinstance(projection_dimension, bool) or not isinstance(projection_dimension, Integral)
                or not 1 <= projection_dimension <= min(x.shape)):
            raise ValueError("projection dimension must be within the training matrix dimensions")
        _, _, vectors = np.linalg.svd(standardized, full_matrices=False)
        projection = vectors[:projection_dimension].T.copy()
        for column in range(projection.shape[1]):
            pivot = np.argmax(np.abs(projection[:, column]))
            if projection[pivot, column] < 0:
                projection[:, column] *= -1
    design = np.column_stack((np.ones(len(training)), standardized @ projection))
    shape = (design.shape[1], len(states) - 1)

    def objective(parameters):
        coefficients = np.column_stack((parameters.reshape(shape), np.zeros(design.shape[1])))
        logits = design @ coefficients
        log_probabilities = logits - logsumexp(logits, axis=1, keepdims=True)
        penalty = 0.5 * ridge * np.sum(coefficients[1:] ** 2)
        loss = -np.sum(y * log_probabilities) / len(training) + penalty
        gradient = design.T @ (np.exp(log_probabilities) - y) / len(training)
        gradient[1:] += ridge * coefficients[1:]
        return float(loss), gradient[:, :-1].ravel()

    fitted = minimize(objective, np.zeros(math.prod(shape)), jac=True, method="L-BFGS-B",
                      options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8})
    if not fitted.success or not np.isfinite(fitted.x).all():
        raise ValueError(f"training optimization failed: {fitted.message}")
    coefficients = np.column_stack((fitted.x.reshape(shape), np.zeros(design.shape[1])))
    digest = hashlib.sha256(json.dumps([_report(row) for row in training], sort_keys=True,
                                       allow_nan=False).encode()).hexdigest()
    organism, strain = next(iter(hosts))
    return TrainedStateDistribution(
        states, features, basis, center, scale, projection, coefficients, x.min(axis=0), x.max(axis=0),
        tuple(row.sample_id for row in training), tuple(sorted(train_groups)),
        tuple(row.source for row in training), digest, organism, strain, ridge, projection_dimension,
        float(fitted.fun), int(fitted.nit),
    )


@dataclass(frozen=True)
class FixedDistribution:
    fractions: Mapping[str, float]
    basis: str
    source: Source

    def __post_init__(self):
        _source(self.source)
        if self.basis not in {"biomass", "cell"}:
            raise ValueError("fraction basis must be biomass or cell")
        object.__setattr__(self, "fractions", MappingProxyType(_fractions(self.fractions)))

    def predict(self, environment):
        return dict(self.fractions)

    def report(self):
        return {"method": "fixed_distribution", **_report(self), "biological_validation": False}


@dataclass(frozen=True)
class FixedTransition:
    states: tuple[str, str]
    feature: Feature
    midpoint: float
    width: float
    source: Source
    basis: str = "biomass"

    def __post_init__(self):
        _source(self.source)
        if len(self.states) != 2 or len(set(self.states)) != 2:
            raise ValueError("fixed transition needs two distinct states")
        if self.basis not in {"biomass", "cell"}:
            raise ValueError("fraction basis must be biomass or cell")
        _number(self.midpoint, "transition midpoint")
        if _number(self.width, "transition width", 0.0) == 0:
            raise ValueError("transition width must be positive")

    def predict(self, environment):
        fraction = float(expit((self.feature.read(environment) - self.midpoint) / self.width))
        return {self.states[0]: 1.0 - fraction, self.states[1]: fraction}

    def report(self):
        return {"method": "fixed_transition", **_report(self), "biological_validation": False}


@dataclass(frozen=True)
class UptakePolicy:
    exchange: str
    maximum: float
    source: Source
    unit: str = "mmol/gDW/h"
    feature: Feature | None = None
    half_saturation: float | None = None
    mode: str = "cap"

    def __post_init__(self):
        _text(self.exchange, "uptake exchange")
        _source(self.source)
        specific_rate(self.maximum, self.unit)
        if self.unit == "1/h":
            raise ValueError("uptake is not a growth-rate unit")
        if self.mode not in {"cap", "fixed"}:
            raise ValueError("uptake mode must be cap or fixed")
        if (self.feature is None) != (self.half_saturation is None):
            raise ValueError("kinetic uptake needs both a concentration feature and half_saturation")
        if self.half_saturation is not None and _number(self.half_saturation, "half saturation", 0.0) == 0:
            raise ValueError("half saturation must be positive in the feature's declared concentration unit")

    def magnitude(self, environment):
        if self.feature is None:
            return float(self.maximum)
        concentration = _number(self.feature.read(environment), "uptake concentration", 0.0)
        return self.maximum * concentration / (self.half_saturation + concentration)


@dataclass(frozen=True)
class FluxCapacity:
    reaction: str
    lower: float | None
    upper: float | None
    source: Source
    unit: str = "model_flux"

    def __post_init__(self):
        _text(self.reaction, "capacity reaction")
        _text(self.unit, "capacity unit")
        _source(self.source)
        if self.lower is None and self.upper is None:
            raise ValueError("a capacity must supply a bound")
        for name in ("lower", "upper"):
            if getattr(self, name) is not None:
                _number(getattr(self, name), name)
        if self.lower is not None and self.upper is not None and self.lower > self.upper:
            raise ValueError("capacity lower bound exceeds upper bound")


@dataclass(frozen=True)
class ProteinBudget:
    exchange: str
    amount: float
    unit: str
    model_unit: str
    source: Source

    def __post_init__(self):
        _text(self.exchange, "protein exchange")
        _source(self.source)
        _number(self.amount, "protein budget", 0.0)
        if self.unit not in {"g/gDW", "mg/gDW"} or self.model_unit not in {"g/gDW", "mg/gDW"}:
            raise ValueError("protein amount and model pool units must be explicit g/gDW or mg/gDW")

    def in_model_units(self):
        grams = self.amount * (1e-3 if self.unit == "mg/gDW" else 1.0)
        return grams * (1000.0 if self.model_unit == "mg/gDW" else 1.0)


@dataclass(frozen=True)
class ObjectiveTask:
    reaction: str
    fraction: float
    source: Source
    direction: str = "forward"

    def __post_init__(self):
        _text(self.reaction, "objective reaction")
        _source(self.source)
        if not 0 <= _number(self.fraction, "objective fraction") <= 1:
            raise ValueError("objective fraction must be in [0, 1]")
        if self.direction not in {"forward", "uptake", "secretion"}:
            raise ValueError("objective direction must be forward, uptake, or secretion")

    def coefficients(self, model):
        return ({self.reaction: 1.0} if self.direction == "forward" else
                exchange_coefficients(model, self.reaction, direction=self.direction))


@dataclass(frozen=True)
class StatePolicy:
    name: str
    source: Source
    uptakes: tuple[UptakePolicy, ...] = ()
    enzyme_capacities: tuple[FluxCapacity, ...] = ()
    tasks: tuple[ObjectiveTask, ...] = ()
    measured_tasks: tuple[MetabolicTask, ...] = ()
    protein_budget: ProteinBudget | None = None
    dry_mass_g_per_cell: float | None = None
    dry_mass_source: Source | None = None
    biomass_reaction: str = "r_2111"
    selection: str = "pfba"

    def __post_init__(self):
        _text(self.name, "state name")
        _source(self.source)
        _text(self.biomass_reaction, "biomass reaction")
        if self.selection not in {"fba", "pfba"}:
            raise ValueError("flux selection must be fba or pfba")
        for name, key in (("uptakes", "exchange"), ("enzyme_capacities", "reaction"),
                          ("tasks", "reaction"), ("measured_tasks", "reaction")):
            rows = tuple(getattr(self, name))
            if len({getattr(row, key) for row in rows}) != len(rows):
                raise ValueError(f"duplicate {name} would obscure a repeated constraint or penalty")
            object.__setattr__(self, name, rows)
        if self.dry_mass_g_per_cell is not None:
            if _number(self.dry_mass_g_per_cell, "cell dry mass", 0.0) == 0:
                raise ValueError("cell dry mass must be positive")
            _source(self.dry_mass_source)


@dataclass(frozen=True)
class MechanisticPolicy:
    state: StatePolicy


@dataclass(frozen=True)
class FixedObjectivePolicy:
    state: StatePolicy


@dataclass(frozen=True)
class MixturePolicy:
    states: tuple[StatePolicy, ...]
    distribution: TrainedStateDistribution | FixedDistribution | FixedTransition

    def __post_init__(self):
        states = tuple(self.states)
        if not states or len({state.name for state in states}) != len(states):
            raise ValueError("mixture states must be nonempty and distinct")
        if not isinstance(self.distribution, (TrainedStateDistribution, FixedDistribution, FixedTransition)):
            raise TypeError("mixture needs an explicit fitted or fixed distribution")
        object.__setattr__(self, "states", states)


@dataclass(frozen=True)
class StateResult:
    name: str
    status: str
    fluxes: dict[str, float] | None
    growth_rate: float | None
    mass_balance_residual: float | None
    bounds_residual: float | None
    uptake_constraints: tuple[dict, ...] = ()
    objective_floors: tuple[dict, ...] = ()
    message: str = ""


def _mass_balance_residual(model, fluxes):
    balances = {metabolite.id: [] for metabolite in model.metabolites}
    for reaction in model.reactions:
        for metabolite, coefficient in reaction.metabolites.items():
            balances[metabolite.id].append(coefficient * fluxes[reaction.id])
    return max((abs(math.fsum(values)) for values in balances.values()), default=0.0)


def _apply_protein_budget(model, budget):
    amount = budget.in_model_units()
    coefficients = exchange_coefficients(model, budget.exchange)
    for rid, coefficient in coefficients.items():
        reaction = model.reactions.get_by_id(rid)
        lower, upper = reaction.bounds
        if coefficient > 0 and upper > 0:
            upper = amount / coefficient
        elif coefficient < 0 and lower < 0:
            lower = amount / coefficient
        if lower > upper:
            raise InfeasibleRegulation("explicit protein budget contradicts an existing required draw")
        reaction.bounds = lower, upper
    model.add_cons_vars([model.problem.Constraint(_expression(model, coefficients), ub=amount,
                                                  name=_constraint_name("allocation_protein_", coefficients))])


def _solve_state(model, environment, state, fixed_objective=False):
    uptakes, floors = [], []

    def failed(status, message, mass_balance_residual=None, bounds_residual=None):
        return StateResult(state.name, status, None, None, mass_balance_residual, bounds_residual,
                           tuple(uptakes), tuple(floors), message)

    with model as scoped:
        try:
            scoped.reactions.get_by_id(state.biomass_reaction)
            for uptake in state.uptakes:
                report = constrain_uptake(scoped, uptake.exchange, uptake.magnitude(environment),
                                          unit=uptake.unit, mode=uptake.mode, source=uptake.source.reference)
                uptakes.append(report.report())
            if state.protein_budget is not None:
                _apply_protein_budget(scoped, state.protein_budget)
            for capacity in state.enzyme_capacities:
                reaction = scoped.reactions.get_by_id(capacity.reaction)
                lower = reaction.lower_bound if capacity.lower is None else max(reaction.lower_bound, capacity.lower)
                upper = reaction.upper_bound if capacity.upper is None else min(reaction.upper_bound, capacity.upper)
                if lower > upper:
                    return failed("infeasible", f"{capacity.reaction}: declared capacity contradicts native bounds")
                reaction.bounds = lower, upper
            for measured in state.measured_tasks:
                if measured.condition_id != environment.condition_id or measured.host != environment.organism:
                    raise ValueError("measured state task has an unmatched condition or host")
                floor = measured.require_measurement()
                _set_floor(scoped, measured.coefficients(scoped), floor)
                floors.append({"reaction": measured.reaction, "floor": floor,
                               "kind": "measured_lower_bound", "source": measured.source})
            tasks = ((ObjectiveTask(state.biomass_reaction, 1.0, state.source),)
                     if fixed_objective or not state.tasks else state.tasks)
            for task in tasks:
                coefficients = task.coefficients(scoped)
                scoped.objective = scoped.problem.Objective(_expression(scoped, coefficients), direction="max")
                optimum = growth_or_none(scoped)
                status = _solve_status(scoped, optimum)
                if status != "optimal":
                    return failed(status, f"{task.reaction}: optimization did not solve")
                if optimum < -scoped.tolerance:
                    raise ValueError("state objective is negative; declare its physical uptake/secretion direction")
                floor = task.fraction * max(0.0, optimum)
                _set_floor(scoped, coefficients, floor)
                floors.append({"reaction": task.reaction, "fraction": task.fraction, "optimum": optimum,
                               "floor": floor, "kind": "fraction_of_optimum_lower_bound",
                               "source": asdict(task.source)})
            optimum = growth_or_none(scoped)
            status = _solve_status(scoped, optimum)
            if status != "optimal":
                return failed(status, "state objective floors are jointly unsolved")
            if state.selection == "pfba":
                solution = cobra.flux_analysis.pfba(scoped)
            else:
                solution = cobra.core.solution.get_solution(scoped)
            if solution.status != "optimal":
                return failed(solution.status, "final flux selection did not solve")
            fluxes = {rid: float(value) for rid, value in solution.fluxes.items()}
            if not all(math.isfinite(value) for value in fluxes.values()):
                return failed("numerical_failure", "nonfinite full flux vector")
            residual = _mass_balance_residual(scoped, fluxes)
            bound_error = max((max(0.0, reaction.lower_bound - fluxes[reaction.id],
                                   fluxes[reaction.id] - reaction.upper_bound) for reaction in scoped.reactions),
                              default=0.0)
            tolerance = max(1e-8, 10 * scoped.tolerance)
            if residual > tolerance or bound_error > tolerance:
                return failed("numerical_failure", f"flux residuals exceed {tolerance:g}: S={residual:g}, bounds={bound_error:g}",
                              residual, bound_error)
            return StateResult(state.name, "optimal", fluxes, fluxes[state.biomass_reaction], residual,
                               bound_error, tuple(uptakes), tuple(floors), "complete flux vector; no extra stress penalty")
        except InfeasibleRegulation as exc:
            return failed(exc.status, str(exc))
        except cobra.exceptions.OptimizationError as exc:
            status = str(scoped.solver.status)
            return failed(status if status in {"infeasible", "unbounded"} else "numerical_failure", str(exc))


@dataclass(frozen=True)
class AllocationResult:
    policy: str
    status: str
    condition_id: str
    fluxes: dict[str, float] | None
    growth_rate: float | None
    biomass_fractions: dict[str, float]
    cell_fractions: dict[str, float] | None
    cell_population_growth_rate: float | None
    mean_cell_dry_mass_g: float | None
    mass_balance_residual: float | None
    states: tuple[StateResult, ...]
    strain_transfer: str | None
    provenance: dict = field(default_factory=dict)

    def report(self, *, include_fluxes=True):
        report = _report(self)
        report["biological_validation"] = False
        report["flux_convention"] = "complete per-gDW model flux vector; biomass is 1/h; enzyme pseudofluxes retain declared model units"
        report["blending"] = "biomass fractions multiply every flux once; no post-blend penalty or infeasible-state renormalization"
        if not include_fluxes:
            report.pop("fluxes")
            for state in report["states"]:
                state.pop("fluxes")
        return report

    def forward(self, reaction=None):
        value = None if self.fluxes is None else (self.growth_rate if reaction is None else self.fluxes[reaction])
        status = self.status if self.status in {"optimal", "infeasible", "unbounded", "numerical_failure"} else "unsupported"
        return ForwardEvaluation(status, value, {"policy": self.policy, "condition_id": self.condition_id})


def solve_allocation(
    model: cobra.Model, identity: ModelIdentity, environment: Environment,
    policy: MechanisticPolicy | MixturePolicy | FixedObjectivePolicy, *, strain_transfer: str | None = None,
) -> AllocationResult:
    """One interface, mutually exclusive mechanisms. The caller's model is left unchanged."""
    validate_host(model, identity, environment, strain_transfer=strain_transfer)
    if isinstance(policy, MixturePolicy):
        states = policy.states
        fractions = _fractions(policy.distribution.predict(environment), [state.name for state in states])
        basis = policy.distribution.basis
        kind = "fixed_transition" if isinstance(policy.distribution, FixedTransition) else "mixture"
        distribution = policy.distribution.report()
    elif isinstance(policy, (MechanisticPolicy, FixedObjectivePolicy)):
        states, fractions, basis = (policy.state,), {policy.state.name: 1.0}, "biomass"
        kind = "mechanistic" if isinstance(policy, MechanisticPolicy) else "fixed_objective"
        distribution = {"method": "single_state", "biological_validation": False}
    else:
        raise TypeError("choose one mechanistic, mixture, or fixed-objective allocation policy")
    if len({state.biomass_reaction for state in states}) != 1:
        raise ValueError("full-flux blending requires the same biomass definition in every state")
    masses = {state.name: state.dry_mass_g_per_cell for state in states}
    have_masses = all(mass is not None for mass in masses.values())
    cell_fractions, mean_mass = None, None
    if basis == "cell":
        if not have_masses:
            raise ValueError("cell fractions require a sourced dry-mass-per-cell for every state")
        cell_fractions = fractions
        mean_mass = math.fsum(fractions[name] * masses[name] for name in fractions)
        biomass_fractions = {name: fraction * masses[name] / mean_mass for name, fraction in fractions.items()}
    else:
        biomass_fractions = fractions
        if have_masses:
            cells = {name: fraction / masses[name] for name, fraction in fractions.items()}
            total_cells = math.fsum(cells.values())
            cell_fractions = {name: value / total_cells for name, value in cells.items()}
            mean_mass = 1.0 / total_cells
    results = []
    for state in states:
        if biomass_fractions[state.name] == 0:
            results.append(StateResult(state.name, "inactive", None, None, None, None,
                                       message="exact zero mixture weight; no LP solved"))
        else:
            results.append(_solve_state(model, environment, state, isinstance(policy, FixedObjectivePolicy)))
    active = [result for result in results if biomass_fractions[result.name] > 0]
    failures = [result.status for result in active if result.status != "optimal"]
    status = failures[0] if failures else "optimal"
    fluxes, growth, cell_growth, residual = None, None, None, None
    if not failures:
        fluxes = {reaction.id: math.fsum(biomass_fractions[result.name] * result.fluxes[reaction.id]
                                         for result in active) for reaction in model.reactions}
        growth = fluxes[states[0].biomass_reaction]
        residual = _mass_balance_residual(model, fluxes)
        if cell_fractions is not None:
            cell_growth = math.fsum(cell_fractions[result.name] * result.growth_rate for result in active)
        if residual > max(1e-8, 10 * model.tolerance):
            status, fluxes, growth, cell_growth = "numerical_failure", None, None, None
    provenance = {
        "model": _report(identity), "model_taxonomy_annotation": model.annotation.get("taxonomy"),
        "environment": _report(environment), "state_policies": [_report(state) for state in states],
        "distribution": distribution, "protein_budget_calibration": "never automatic; unchanged unless explicitly declared",
        "solver": model.solver.interface.__name__, "solver_tolerance": model.tolerance,
        "selected_state_optima": "FBA/pFBA selections, not identified unique fluxes or uncertainty bounds",
    }
    if isinstance(policy, MixturePolicy) and isinstance(policy.distribution, TrainedStateDistribution):
        provenance["extrapolated_features"] = list(policy.distribution.extrapolated_features(environment))
    return AllocationResult(kind, status, environment.condition_id, fluxes, growth, biomass_fractions,
                            cell_fractions, cell_growth, mean_mass, residual, tuple(results), strain_transfer, provenance)


def reactor_derivative(
    model: cobra.Model, result: AllocationResult, *, biomass_gdw_per_l: float,
    concentrations_mmol_per_l: Mapping[str, float], exchanges: Mapping[str, str],
    dilution_per_h: float = 0.0, feed_mmol_per_l: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Instantaneous dFBA mass balance, not an integration or fitted reactor model.

    Medium keys map to boundary reaction ids. Biomass is gDW/L, all medium amounts
    are mmol/L, and time is hours. Cells are not substituted for dry biomass.
    """
    if result.status != "optimal" or result.fluxes is None:
        raise ValueError("an unsolved allocation has no reactor derivative")
    biomass = _number(biomass_gdw_per_l, "biomass concentration", 0.0)
    dilution = _number(dilution_per_h, "dilution rate", 0.0)
    feed = {} if feed_mmol_per_l is None else dict(feed_mmol_per_l)
    if set(concentrations_mmol_per_l) != set(exchanges) or not set(feed) <= set(exchanges):
        raise ValueError("medium concentrations, exchanges, and feed must use matching component names")
    rates = {"biomass_gDW_per_L_per_h": (result.growth_rate - dilution) * biomass}
    for name, exchange in exchanges.items():
        concentration = _number(concentrations_mmol_per_l[name], "medium concentration", 0.0)
        inlet = _number(feed.get(name, 0.0), "feed concentration", 0.0)
        coefficients = exchange_coefficients(model, exchange, direction="secretion")
        specific = math.fsum(coefficient * result.fluxes[rid] for rid, coefficient in coefficients.items())
        rates[name] = biomass * specific + dilution * (inlet - concentration)
    return rates
