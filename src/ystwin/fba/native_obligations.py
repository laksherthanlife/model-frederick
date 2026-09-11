from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, dataclass, field
import math

import cobra
from cobra.core.formula import elements_and_molecular_weights
from cobra.util.array import create_stoichiometric_matrix
import numpy as np

from .dynamic_rates import _validated_optimize
from .solver import configure

__all__ = [
    "NativeCarbonObligation", "NativeObligationRoles", "NativeObligationReport",
    "impose_native_obligation",
]


def _identifier(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be an explicit nonempty identifier")
    return value


def _identifiers(values, label, *, nonempty=False):
    if isinstance(values, (str, bytes, Mapping)):
        raise ValueError(f"{label} must be a sequence of explicit identifiers")
    try:
        result = tuple(_identifier(value, label) for value in values)
    except TypeError as exc:
        raise ValueError(f"{label} must be a sequence of explicit identifiers") from exc
    if len(set(result)) != len(result) or (nonempty and not result):
        raise ValueError(f"{label} must contain unique identifiers and cannot omit required roles")
    return result


@dataclass(frozen=True)
class NativeCarbonObligation:
    carbon_fraction: float
    source: str
    source_medium: str
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        try:
            fraction = float(self.carbon_fraction)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("carbon_fraction must be a finite scalar in [0, 1]") from exc
        if (isinstance(self.carbon_fraction, (bool, np.bool_, np.ndarray))
                or not math.isfinite(fraction) or not 0 <= fraction <= 1):
            raise ValueError("carbon_fraction must be in [0, 1]; carbon credit/remobilization is unsupported")
        _identifier(self.source, "source provenance")
        _identifier(self.source_medium, "source_medium")
        if not isinstance(self.metadata, Mapping):
            raise ValueError("native coefficient metadata must be a mapping")
        metadata = deepcopy(dict(self.metadata))
        if metadata.get("native_training_scope", "glucose_only") != "glucose_only":
            raise ValueError("this native carbon bridge only supports glucose_only training scope")
        object.__setattr__(self, "carbon_fraction", fraction)
        object.__setattr__(self, "metadata", metadata)


@dataclass(frozen=True)
class NativeObligationRoles:
    native_metabolite_id: str
    glucose_metabolite_id: str
    glucose_exchange_ids: tuple[str, ...]
    native_export_ids: tuple[str, ...] = ()
    resource_metabolite_ids: tuple[str, ...] = ()

    def __post_init__(self):
        _identifier(self.native_metabolite_id, "native_metabolite_id")
        _identifier(self.glucose_metabolite_id, "glucose_metabolite_id")
        if self.native_metabolite_id == self.glucose_metabolite_id:
            raise ValueError("native glycerol and glucose must have distinct species roles")
        for name in ("glucose_exchange_ids", "native_export_ids", "resource_metabolite_ids"):
            object.__setattr__(self, name, _identifiers(getattr(self, name), name,
                                                       nonempty=name == "glucose_exchange_ids"))
        if set(self.glucose_exchange_ids).intersection(self.native_export_ids):
            raise ValueError("native exports and glucose boundaries must have distinct roles")


@dataclass(frozen=True)
class NativeObligationReport:
    demand_reaction_id: str
    coupling_constraint_id: str
    uptake_constraint_id: str
    native_moles_per_glucose_mole: float
    glucose_flux_coefficients: dict[str, float]
    terminal_flux_coefficients: dict[str, float]
    metadata: dict


def _boundary(model, rid):
    reaction = model.reactions.get_by_id(rid)
    if len(reaction.metabolites) != 1:
        raise ValueError(f"{rid}: native boundary roles require a one-metabolite reaction")
    metabolite, coefficient = next(iter(reaction.metabolites.items()))
    if not math.isfinite(coefficient) or coefficient == 0:
        raise ValueError(f"{rid}: boundary stoichiometry must be finite and nonzero")
    return reaction, metabolite, float(coefficient)


def _can_supply(reaction, coefficient):
    return (coefficient > 0 and reaction.upper_bound > 0
            or coefficient < 0 and reaction.lower_bound < 0)


def _chemical_role(metabolite, elements, label):
    if not metabolite.formula or metabolite.elements != elements or metabolite.charge not in (None, 0):
        raise ValueError(f"{metabolite.id}: {label} role requires its neutral molecular formula {elements}")


def _validate_roles_and_medium(model, roles):
    native = model.metabolites.get_by_id(roles.native_metabolite_id)
    glucose = model.metabolites.get_by_id(roles.glucose_metabolite_id)
    _chemical_role(native, {"C": 3, "H": 8, "O": 3}, "native glycerol")
    _chemical_role(glucose, {"C": 6, "H": 12, "O": 6}, "glucose")
    glucose_coefficients, native_coefficients = {}, {}
    native_species = {native.id}
    for rid in roles.glucose_exchange_ids:
        reaction, metabolite, coefficient = _boundary(model, rid)
        if metabolite != glucose:
            raise ValueError(f"{rid}: glucose boundary does not contain the declared glucose species")
        glucose_coefficients[rid] = coefficient
    for rid in roles.native_export_ids:
        reaction, metabolite, coefficient = _boundary(model, rid)
        _chemical_role(metabolite, native.elements, "native glycerol export")
        if _can_supply(reaction, coefficient):
            raise ValueError(f"{rid}: native uptake/carbon credit requires an explicit inventory model and is unsupported")
        native_coefficients[rid] = -coefficient
        native_species.add(metabolite.id)
    for mid in roles.resource_metabolite_ids:
        resource = model.metabolites.get_by_id(mid)
        if (resource.elements or {}).get("C", 0):
            raise ValueError(f"{mid}: a carbon-containing species cannot be a noncarbon resource role")
    boundary_inputs = {}
    for reaction in model.boundary:
        if len(reaction.metabolites) == 1 and reaction.bounds != (0., 0.):
            metabolite = next(iter(reaction.metabolites))
            if metabolite == glucose and reaction.id not in glucose_coefficients:
                raise ValueError(f"{reaction.id}: every active glucose boundary needs an explicit role")
            if metabolite.id in native_species and reaction.id not in native_coefficients:
                raise ValueError(f"{reaction.id}: every active native terminal boundary needs an explicit role")
        for metabolite, coefficient in reaction.metabolites.items():
            if not math.isfinite(coefficient) or coefficient == 0:
                raise ValueError(f"{reaction.id}: invalid boundary stoichiometry")
            if not _can_supply(reaction, coefficient):
                continue
            if reaction.id in glucose_coefficients:
                continue
            elements = metabolite.elements or {}
            if elements.get("C", 0):
                raise ValueError(f"{reaction.id}: other carbon uptake or inventory remobilization is outside the glucose-only medium scope")
            if not metabolite.formula and metabolite.id not in roles.resource_metabolite_ids:
                raise ValueError(f"{reaction.id}: an unformulated supply requires an explicit noncarbon resource role")
            if metabolite.formula and (not elements or not elements.keys() <= elements_and_molecular_weights.keys()
                                       or any(not math.isfinite(n) or n <= 0 for n in elements.values())):
                raise ValueError(f"{reaction.id}: unknown supply formula cannot establish a carbon-free medium")
        if any(_can_supply(reaction, c) for c in reaction.metabolites.values()):
            boundary_inputs[reaction.id] = tuple(reaction.bounds)
    return glucose_coefficients, native_coefficients, boundary_inputs


def impose_native_obligation(
    model: cobra.Model, *, roles: NativeObligationRoles, obligation: NativeCarbonObligation,
    demand_reaction_id: str,
) -> tuple[cobra.Model, NativeObligationReport]:
    if not isinstance(model, cobra.Model):
        raise ValueError("native obligations require a cobra model")
    if not isinstance(roles, NativeObligationRoles) or not isinstance(obligation, NativeCarbonObligation):
        raise ValueError("explicit native roles and a native carbon obligation are required")
    _identifier(demand_reaction_id, "demand_reaction_id")
    coupling_id = f"{demand_reaction_id}__carbon_commitment"
    uptake_id = f"{demand_reaction_id}__net_glucose"
    if (demand_reaction_id in model.reactions or demand_reaction_id in model.variables
            or coupling_id in model.constraints or uptake_id in model.constraints):
        raise ValueError("native obligation identifiers are already installed or have a solver collision")
    glucose_coefficients, native_coefficients, boundary_inputs = _validate_roles_and_medium(model, roles)
    coefficient = obligation.carbon_fraction * 6. / 3.
    out = model.copy()
    out.notes, out.annotation = deepcopy(model.notes), deepcopy(model.annotation)
    for collection in (out.metabolites, out.reactions, out.genes):
        for entity in collection:
            entity.notes, entity.annotation = deepcopy(entity.notes), deepcopy(entity.annotation)
    demand = cobra.Reaction(demand_reaction_id, lower_bound=0.,
                            upper_bound=math.inf if coefficient > 0 else 0.)
    demand.add_metabolites({out.metabolites.get_by_id(roles.native_metabolite_id): -1.})
    out.add_reactions([demand])
    terminal_coefficients = {demand_reaction_id: 1., **native_coefficients}
    uptake = sum(c * out.reactions.get_by_id(rid).flux_expression for rid, c in glucose_coefficients.items())
    terminal = sum(c * out.reactions.get_by_id(rid).flux_expression for rid, c in terminal_coefficients.items())
    out.add_cons_vars([
        out.problem.Constraint(uptake, lb=0., name=uptake_id),
        out.problem.Constraint(terminal - coefficient * uptake, lb=0., name=coupling_id),
    ])
    configure(out)
    with out:
        out.objective = out.problem.Objective(0., direction="min")
        feasible = _validated_optimize(out, create_stoichiometric_matrix(out, array_type="lil"),
                                       "native obligation feasibility", allow_infeasible=True)
        if feasible is None:
            raise ValueError("native obligation is infeasible under the existing hard constraints")
    metadata = {
        "source": obligation.source,
        "source_medium": obligation.source_medium,
        "source_metadata": deepcopy(obligation.metadata),
        "native_training_scope": "glucose_only",
        "medium_dependence": "the native coefficient may depend on medium, stress and time; an unchanged coefficient is not validated for a changed medium. Only glucose carbon input is supported here.",
        "applied_boundary_input_bounds": boundary_inputs,
        "roles": asdict(roles),
        "carbon_fraction": obligation.carbon_fraction,
        "coefficient_unit": "mol native glycerol/mol utilized glucose (equivalently mmol/mmol)",
        "constraint_equation": "native_terminal >= (6/3)*carbon_fraction*net_glucose_input; net_glucose_input >= 0",
        "uptake_basis": "sum of stoichiometric molecular inputs across all declared glucose boundaries; net solved flux, not uptake capacity or gross exchange cycling",
        "demand_interpretation": "total native terminal glycerol commitment; the new sink and declared existing native exports count together, with no fate-specific rate assignment",
        "zero_coefficient": "the new sink is closed at zero commitment so it cannot open an unforced disposal route",
        "fate_specific": False,
        "absolute_source_flux_identified": False,
        "absolute_conversion": "no source mol/L to mmol/gDW conversion is performed; demand and uptake share the host GEM specific-flux basis",
        "uncertainty": "source-model/learned coefficient conditional on the recorded native condition; uncertainty must be supplied by native evidence, not a downstream objective",
        "remobilization_supported": False,
        "external_native_carbon_credit_supported": False,
        "cost_basis": "one nonnegative native-species demand in the existing stoichiometric solve; existing synthesis/cofactor chemistry, enzyme competition, shared resource budgets and hard constraints are unchanged",
        "resource_roles": "caller-identified non-molecular resources permit existing unformulated supplies; no chemical carbon source is exempted and no resource bound is changed",
        "feasibility_scope": "checked under installation-time hard constraints using the existing pinned physical solver policy; subsequent constraint changes still require a feasible solve",
        "biological_validation": False,
    }
    return out, NativeObligationReport(demand_reaction_id, coupling_id, uptake_id, coefficient,
                                       glucose_coefficients, terminal_coefficients, metadata)
