from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

import cobra
from cobra.util.array import create_stoichiometric_matrix

from ..fba.chemical_task import ChemicalTask, install_chemical_task
from ..fba.dynamic_rates import _nonnegative_flux, _validated_optimize
from ..fba.native_obligations import NativeCarbonObligation, NativeObligationRoles, impose_native_obligation
from ..fba.solver import configure
from .native_physiology import NativeExchangeModel

__all__ = ["BalancedCondition", "HostFluxRoles", "content_from_flux", "predict_balanced_envelope"]


@dataclass(frozen=True)
class BalancedCondition:
    experiment_id: str
    growth_rate_per_h: float
    carbon_input: str
    readout_basis: str
    native_context_assumption: str


@dataclass(frozen=True)
class HostFluxRoles:
    biomass_reaction_id: str
    native: NativeObligationRoles
    oxygen_exchange_ids: tuple[str, ...]
    intracellular_compartment_ids: tuple[str, ...]
    extracellular_compartment_ids: tuple[str, ...]

    def __post_init__(self):
        if not isinstance(self.biomass_reaction_id, str) or not self.biomass_reaction_id.strip():
            raise ValueError("an explicit biomass_reaction_id is required")
        if not isinstance(self.native, NativeObligationRoles):
            raise ValueError("explicit NativeObligationRoles are required")
        for name in ("oxygen_exchange_ids", "intracellular_compartment_ids", "extracellular_compartment_ids"):
            values = getattr(self, name)
            if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
                raise ValueError(f"{name} requires an explicit sequence of identifiers")
            values = tuple(values)
            if any(not isinstance(value, str) or not value.strip() for value in values):
                raise ValueError(f"{name} requires nonempty identifiers")
            if len(set(values)) != len(values) or (name != "extracellular_compartment_ids" and not values):
                raise ValueError(f"{name} requires unique identifiers and cannot omit required roles")
            object.__setattr__(self, name, values)
        if set(self.intracellular_compartment_ids).intersection(self.extracellular_compartment_ids):
            raise ValueError("intracellular and extracellular compartment roles must be disjoint")
        if set(self.oxygen_exchange_ids).intersection(
                (*self.native.glucose_exchange_ids, *self.native.native_export_ids)):
            raise ValueError("oxygen, glucose and native terminal boundary roles must be distinct")


def _molecular_uptake_cap(model, identifiers, capacity, *, label, metabolite_id=None):
    if not math.isfinite(capacity) or capacity < 0:
        raise ValueError(f"{label} molecular uptake capacity must be finite and nonnegative")

    def matches(metabolite):
        if metabolite_id is not None:
            return metabolite.id == metabolite_id
        return metabolite.elements == {"O": 2} and metabolite.charge in (None, 0)

    coefficients = {}
    for identifier in identifiers:
        reaction = model.reactions.get_by_id(identifier)
        if len(reaction.metabolites) != 1:
            raise ValueError(f"{label} boundary {identifier} must contain exactly one molecular species")
        metabolite, coefficient = next(iter(reaction.metabolites.items()))
        if not matches(metabolite):
            raise ValueError(f"{label} boundary {identifier} does not contain the declared molecular species")
        if not math.isfinite(coefficient) or coefficient == 0:
            raise ValueError(f"{label} boundary {identifier} requires finite nonzero stoichiometry")
        coefficients[identifier] = float(coefficient)
    for reaction in model.boundary:
        if (reaction.bounds != (0, 0) and reaction.id not in coefficients
                and any(matches(metabolite) for metabolite in reaction.metabolites)):
            raise ValueError(f"every active {label} boundary needs an explicit role: {reaction.id}")
    constraint_id = f"__native_{label}_molecular_uptake_cap"
    if constraint_id in model.constraints:
        raise ValueError(f"{label} molecular uptake constraint is already installed")
    uptake = sum(coefficient * model.reactions.get_by_id(identifier).flux_expression
                 for identifier, coefficient in coefficients.items())
    model.add_cons_vars([model.problem.Constraint(uptake, lb=0.0, ub=capacity, name=constraint_id)])
    return {"constraint_id": constraint_id, "flux_coefficients": coefficients,
            "lower_mmol_per_base_gdw_h": 0.0, "upper_mmol_per_base_gdw_h": capacity,
            "basis": "net stoichiometric molecular input across explicitly declared boundaries"}


def content_from_flux(flux_mmol_per_base_gdw_h, growth_rate_per_h, molar_mass_g_per_mol):
    values = (flux_mmol_per_base_gdw_h, growth_rate_per_h, molar_mass_g_per_mol)
    if any(isinstance(value, bool) or not math.isfinite(float(value)) for value in values):
        raise ValueError("balanced content inputs must be finite real numbers")
    flux, growth, mass = map(float, values)
    if flux < 0 or growth <= 0 or mass <= 0:
        raise ValueError("balanced content requires nonnegative flux and positive growth/molar mass")
    ratio = flux * mass / (1000.0 * growth)
    if not math.isfinite(ratio):
        raise ValueError("balanced product/base-host mass ratio is nonfinite")
    return 1000.0 * ratio / (1.0 + ratio)


def _failure(condition, status, reason, details):
    return ({"experiment_id": condition.experiment_id, "status": status,
             "lower": None, "upper": None, "reason": reason}, details)


def predict_balanced_envelope(
    host: cobra.Model, *, task: ChemicalTask, condition: BalancedCondition,
    physiology: NativeExchangeModel, obligations: tuple[NativeCarbonObligation, ...], roles: HostFluxRoles,
):
    if not isinstance(condition.experiment_id, str) or not condition.experiment_id.strip():
        raise ValueError("an explicit experiment identifier is required")
    if not isinstance(roles, HostFluxRoles):
        raise ValueError("explicit host flux and compartment roles are required")
    details = {"prediction_scope": "feasible_envelope_not_realized_titre", "scenarios": []}
    if condition.carbon_input != "glucose":
        return _failure(condition, "unsupported", "native exchange/obligation training supports glucose carbon only", details)
    if condition.readout_basis != "balanced_growth_intracellular_content":
        return _failure(condition, "unsupported", "this frozen readout requires balanced-growth intracellular content", details)
    if not isinstance(condition.native_context_assumption, str) or not condition.native_context_assumption.strip():
        return _failure(condition, "unsupported", "native strain/medium transfer assumptions must be explicit", details)
    try:
        growth = float(condition.growth_rate_per_h)
    except (TypeError, ValueError, OverflowError):
        return _failure(condition, "unsupported", "growth is not a finite positive imposed condition", details)
    if not math.isfinite(growth) or growth <= 0:
        return _failure(condition, "unsupported", "growth is not a finite positive imposed condition", details)
    if not obligations or not all(isinstance(item, NativeCarbonObligation) for item in obligations):
        raise ValueError("explicit native-obligation scenarios are required")
    physiological_predictions = {row["observable_id"]: row for row in physiology.predict([growth])}
    required = ("GlucoseUptake", "O2uptake")
    if any(physiological_predictions[name]["prediction_status"] != "quantified" for name in required):
        return _failure(condition, "unsupported", "growth lies outside quantified native uptake support", details)
    glucose, oxygen = (float(physiological_predictions[name]["prediction"]) for name in required)
    details.update(
        imposed_growth_per_h=growth, native_glucose_cap_mmol_per_base_gdw_h=glucose,
        native_oxygen_cap_mmol_per_base_gdw_h=oxygen,
        native_context_assumption=condition.native_context_assumption,
        physiological_inputs={name: physiological_predictions[name] for name in required},
        mass_basis="native rates per base-host dry mass; new intracellular product contributes to the total-dry-mass denominator",
        unmodeled=["realized heterologous allocation", "condition-matched active enzyme abundances",
                   "additional retained-native-pool mass", "full source strain/medium transfer"],
    )
    bounds = []
    for index, obligation in enumerate(obligations):
        scenario = {"carbon_fraction": obligation.carbon_fraction, "native_source": obligation.source}
        try:
            model = host.copy()
            configure(model)
            compartment_ids = {*roles.intracellular_compartment_ids, *roles.extracellular_compartment_ids}
            if not compartment_ids <= model.compartments.keys():
                raise ValueError("compartment roles contain identifiers absent from the host model")
            scenario["molecular_uptake_caps"] = {
                "glucose": _molecular_uptake_cap(
                    model, roles.native.glucose_exchange_ids, glucose, label="glucose",
                    metabolite_id=roles.native.glucose_metabolite_id),
                "oxygen": _molecular_uptake_cap(model, roles.oxygen_exchange_ids, oxygen, label="oxygen"),
            }
            biomass = model.reactions.get_by_id(roles.biomass_reaction_id)
            if not biomass.lower_bound <= growth <= biomass.upper_bound:
                raise ValueError("imposed growth conflicts with an existing biomass bound")
            biomass.bounds = (growth, growth)
            model, native_report = impose_native_obligation(
                model, roles=roles.native, obligation=obligation,
                demand_reaction_id=f"__native_commitment_{index}")
            model, installed = install_chemical_task(model, task)
            compartment = model.metabolites.get_by_id(installed.metabolite_id).compartment
            if compartment in roles.extracellular_compartment_ids:
                raise ValueError("intracellular content cannot be inferred from an extracellular output task")
            if compartment not in roles.intracellular_compartment_ids:
                raise ValueError("output compartment lacks an explicit intracellular role")
            stoichiometry = create_stoichiometric_matrix(model, array_type="lil")
            product = model.reactions.get_by_id(installed.reaction_id)
            fluxes = []
            for direction in ("min", "max"):
                model.objective = model.problem.Objective(product.flux_expression, direction=direction)
                solution = _validated_optimize(model, stoichiometry, f"unseen task {direction}", allow_infeasible=True)
                if solution is None:
                    scenario.update(status="infeasible", reason="chemical task is infeasible with native obligations")
                    details["scenarios"].append(scenario)
                    return _failure(condition, "infeasible", scenario["reason"], details)
                fluxes.append(_nonnegative_flux(solution[0], "task flux"))
            if fluxes[1] < fluxes[0]:
                raise RuntimeError("task flux envelope is inverted")
            content = [content_from_flux(value, growth, installed.molar_mass_g_per_mol) for value in fluxes]
            bounds.append(content)
            scenario.update(status="ok", flux_lower_mmol_per_base_gdw_h=fluxes[0],
                            flux_upper_mmol_per_base_gdw_h=fluxes[1], lower=content[0], upper=content[1],
                            native_obligation=native_report.metadata, chemical_task=installed.metadata)
            details["scenarios"].append(scenario)
        except (ValueError, RuntimeError, KeyError) as exc:
            scenario.update(status="failed", reason=str(exc))
            details["scenarios"].append(scenario)
            return _failure(condition, "failed", str(exc), details)
    return ({"experiment_id": condition.experiment_id, "status": "ok",
             "lower": min(pair[0] for pair in bounds), "upper": max(pair[1] for pair in bounds),
             "reason": None}, details)
