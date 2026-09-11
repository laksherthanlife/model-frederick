from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, fields

import cobra

__all__ = ["ChemicalTask", "EnzymeDefinition", "MetaboliteDefinition", "ReactionDefinition",
           "InstalledChemicalTask", "install_chemical_task"]


@dataclass(frozen=True)
class MetaboliteDefinition:
    id: str
    formula: str
    charge: float
    compartment: str


@dataclass(frozen=True)
class EnzymeDefinition:
    id: str
    molecular_weight_g_per_mol: float
    abundance_mmol_per_gdw: float | None
    source: str


@dataclass(frozen=True)
class ReactionDefinition:
    id: str
    stoichiometry: tuple[tuple[str, float], ...]
    source: str
    enzyme_id: str | None = None
    kcat_per_s: float | None = None


@dataclass(frozen=True)
class ChemicalTask:
    metabolites: tuple[MetaboliteDefinition, ...]
    reactions: tuple[ReactionDefinition, ...]
    output_metabolite_id: str
    enzymes: tuple[EnzymeDefinition, ...] = ()
    protein_pool_metabolite_id: str | None = None

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, payload):
        if not isinstance(payload, dict) or set(payload) - {item.name for item in fields(cls)}:
            raise ValueError("unknown chemical task fields")
        try:
            metabolites = tuple(MetaboliteDefinition(**item) for item in payload["metabolites"])
            enzymes = tuple(EnzymeDefinition(**item) for item in payload.get("enzymes", ()))
            reactions = tuple(ReactionDefinition(
                **{**item, "stoichiometry": tuple(tuple(pair) for pair in item["stoichiometry"])}
            ) for item in payload["reactions"])
            return cls(metabolites, reactions, payload["output_metabolite_id"], enzymes,
                       payload.get("protein_pool_metabolite_id"))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("invalid or unknown chemical task definition fields") from exc


@dataclass(frozen=True)
class InstalledChemicalTask:
    reaction_id: str
    metabolite_id: str
    molar_mass_g_per_mol: float
    carbon_atoms: float
    metadata: dict


def _number(value, label, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a real finite number")
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{label} must be a real finite number") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{label} must be {'positive and ' if positive else ''}finite")
    return result


def _named(items, label):
    result = {}
    for item in items:
        if not isinstance(item.id, str) or not item.id.strip():
            raise ValueError(f"{label} requires a nonempty identifier")
        if item.id in result:
            raise ValueError(f"duplicate {label} identifier: {item.id}")
        result[item.id] = item
    return result


def _chemical(metabolite):
    if not isinstance(metabolite.formula, str) or not metabolite.formula.strip():
        raise ValueError(f"chemical participant {metabolite.id} requires a formula")
    elements = metabolite.elements
    if not elements or any(not math.isfinite(float(value)) or value <= 0 for value in elements.values()):
        raise ValueError(f"invalid chemical formula for {metabolite.id}")
    _number(metabolite.charge, f"charge of {metabolite.id}")
    _number(metabolite.formula_weight, f"formula mass of {metabolite.id}", positive=True)


def install_chemical_task(host: cobra.Model, task: ChemicalTask):
    if not isinstance(task, ChemicalTask):
        raise ValueError("an explicit ChemicalTask is required")
    definitions = _named(task.metabolites, "metabolite")
    reactions = _named(task.reactions, "reaction")
    enzymes = _named(task.enzymes, "enzyme")
    out = host.copy()
    for definition in definitions.values():
        if definition.compartment not in host.compartments:
            raise ValueError(f"unknown host compartment {definition.compartment}")
        proposed = cobra.Metabolite(definition.id, formula=definition.formula,
                                    charge=_number(definition.charge, "metabolite charge"),
                                    compartment=definition.compartment)
        _chemical(proposed)
        if definition.id in out.metabolites:
            current = out.metabolites.get_by_id(definition.id)
            if (current.elements != proposed.elements or current.charge != proposed.charge
                    or current.compartment != proposed.compartment):
                raise ValueError(f"existing metabolite identity differs: {definition.id}")
        else:
            out.add_metabolites([proposed])
    if task.output_metabolite_id not in out.metabolites:
        raise ValueError("task output metabolite is missing")
    output = out.metabolites.get_by_id(task.output_metabolite_id)
    _chemical(output)
    if output.elements.get("C", 0) <= 0:
        raise ValueError("this carbon-budget task requires a carbon-containing output")
    resources = {}
    if enzymes:
        if task.protein_pool_metabolite_id not in out.metabolites:
            raise ValueError("enzyme constraints require an explicitly identified host protein pool")
        pool = out.metabolites.get_by_id(task.protein_pool_metabolite_id)
        for enzyme in enzymes.values():
            if not isinstance(enzyme.source, str) or not enzyme.source.strip():
                raise ValueError("enzyme parameters require an independent source")
            mass = _number(enzyme.molecular_weight_g_per_mol, "enzyme molecular mass", positive=True)
            abundance = (float("inf") if enzyme.abundance_mmol_per_gdw is None
                         else _number(enzyme.abundance_mmol_per_gdw, "enzyme abundance"))
            if abundance < 0:
                raise ValueError("enzyme abundance cannot be negative")
            mid, rid = f"__task_enzyme_{enzyme.id}", f"__task_draw_{enzyme.id}"
            if mid in out.metabolites or rid in out.reactions:
                raise ValueError("task enzyme bookkeeping collides with the host")
            protein = cobra.Metabolite(mid, compartment=pool.compartment)
            draw = cobra.Reaction(rid, lower_bound=0.0, upper_bound=abundance)
            draw.add_metabolites({pool: -mass / 1000.0, protein: 1.0})
            out.add_reactions([draw])
            resources[enzyme.id] = protein
    used_enzymes, unconstrained = set(), []
    for definition in reactions.values():
        if definition.id in out.reactions:
            raise ValueError(f"reaction identity collides with host: {definition.id}")
        if not isinstance(definition.source, str) or not definition.source.strip():
            raise ValueError("chemical reactions require a source")
        coefficients = {}
        for identifier, value in definition.stoichiometry:
            if identifier in coefficients:
                raise ValueError("duplicate chemical stoichiometry participant")
            if identifier not in out.metabolites:
                raise ValueError(f"unknown chemical participant {identifier}")
            coefficient = _number(value, "stoichiometric coefficient")
            if coefficient == 0:
                raise ValueError("chemical stoichiometry cannot contain zero coefficients")
            metabolite = out.metabolites.get_by_id(identifier)
            _chemical(metabolite)
            coefficients[identifier] = coefficient
        if not any(value < 0 for value in coefficients.values()) or not any(value > 0 for value in coefficients.values()):
            raise ValueError("chemical reactions must consume and produce material")
        reaction = cobra.Reaction(definition.id, lower_bound=0.0, upper_bound=float("inf"))
        reaction.add_metabolites({out.metabolites.get_by_id(key): value for key, value in coefficients.items()})
        imbalance = reaction.check_mass_balance()
        if any(abs(value) > 1e-9 for value in imbalance.values()):
            raise ValueError(f"chemical mass/charge balance failed for {definition.id}: {imbalance}")
        if definition.enzyme_id is not None:
            if definition.enzyme_id not in resources:
                raise ValueError("catalysis refers to an undeclared enzyme")
            kcat = _number(definition.kcat_per_s, "catalytic turnover", positive=True)
            reaction.add_metabolites({resources[definition.enzyme_id]: -1.0 / (3600.0 * kcat)})
            used_enzymes.add(definition.enzyme_id)
        elif definition.kcat_per_s is not None:
            raise ValueError("turnover requires an identified enzyme")
        else:
            unconstrained.append(definition.id)
        reaction.notes["chemical_task_source"] = definition.source
        out.add_reactions([reaction])
    if set(enzymes) != used_enzymes:
        raise ValueError("declared task enzymes must all participate in catalysis")
    if not any(reaction.metabolites[output] > 0 or
               (reaction.metabolites[output] < 0 and reaction.lower_bound < 0)
               for reaction in output.reactions):
        raise ValueError("no host or supplied reaction produces the task output")
    sink_id = "__chemical_task_output"
    if sink_id in out.reactions:
        raise ValueError("a chemical task output is already installed")
    out.add_boundary(output, type="demand", reaction_id=sink_id, lb=0.0, ub=float("inf"))
    digest = hashlib.sha256(json.dumps(task.to_dict(), sort_keys=True, allow_nan=False).encode()).hexdigest()
    metadata = {
        "task_sha256": digest, "new_chemical_reactions": list(reactions),
        "reactions_without_supplied_enzyme_constraints": unconstrained,
        "prediction_scope": "feasible_capacity_not_realized_allocation",
        "mass_charge_balance_checked": True,
        "source_scope": "Supplied chemistry and independent molecular inputs, not outcome-derived flux or allocation constants",
        "unverified": ["reaction thermodynamics", "in vivo enzyme activity", "realized heterologous allocation"],
        "protein_pool_metabolite_id": task.protein_pool_metabolite_id,
    }
    return out, InstalledChemicalTask(sink_id, output.id, float(output.formula_weight),
                                      float(output.elements["C"]), metadata)
