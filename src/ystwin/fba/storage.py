from __future__ import annotations

import math
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from sys import float_info

import cobra
from cobra.util.solver import linear_reaction_coefficients

from .product_panel import biomass_pool_quota


__all__ = ["StorageBasis", "separate_storage_biomass"]

_STORAGE_SPECIES = {"glycogen": "s_0773[c]", "trehalose": "s_1520[c]"}
_CARBOHYDRATE = "s_3718[c]"
_BIOMASS = "s_0450[c]"
_PROTEIN_POOL = "prot_pool[c]"
_UNCHANGED_RATES = ("r_4048", "r_4041", "r_2111", "r_2111_REV")


class _FrozenDict(dict):
    def _refuse_mutation(self, *args, **kwargs):
        raise TypeError("StorageBasis mappings are immutable")

    __setitem__ = __delitem__ = clear = pop = popitem = setdefault = update = __ior__ = _refuse_mutation

    def __reduce__(self):
        return type(self), (dict(self),)

    def __deepcopy__(self, memo):
        return self


def _freeze(value):
    if isinstance(value, Mapping):
        return _FrozenDict({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    return deepcopy(value)


def _finite_number(value, label: str, *, nonnegative: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite number") from exc
    if isinstance(value, bool) or not math.isfinite(number):
        raise ValueError(f"{label} must be a finite number, not a boolean")
    if nonnegative and number < 0.0:
        raise ValueError(f"{label} must be non-negative")
    return number


def _finite_sum(values, label: str) -> float:
    try:
        total = math.fsum(values)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{label} must remain finite; arithmetic overflow") from exc
    if not math.isfinite(total):
        raise ValueError(f"{label} must remain finite; arithmetic overflow")
    return total


def _selected_products(products) -> tuple[str, ...]:
    if isinstance(products, (str, bytes, Mapping)):
        raise ValueError("storage products must be a nonempty sequence of supported product names")
    try:
        names = tuple(products)
    except TypeError as exc:
        raise ValueError("storage products must be a nonempty sequence of supported product names") from exc
    if not names or any(not isinstance(p, str) or p not in _STORAGE_SPECIES for p in names):
        raise ValueError(f"supported storage products are {tuple(_STORAGE_SPECIES)}; got {names!r}")
    if len(set(names)) != len(names):
        raise ValueError("duplicate storage products are not allowed")
    return tuple(p for p in _STORAGE_SPECIES if p in names)


@dataclass(frozen=True)
class StorageBasis:
    original_model_id: str
    structural_mass_fraction: float
    removed_quotas_mmol_per_reference_g: dict[str, float]
    molecular_masses: dict[str, float]
    metadata: dict

    def __post_init__(self):
        if not isinstance(self.original_model_id, str) or not self.original_model_id:
            raise ValueError("original_model_id must identify the reference reconstruction")
        fraction = _finite_number(self.structural_mass_fraction, "structural_mass_fraction")
        if not 0.0 < fraction < 1.0:
            raise ValueError("structural_mass_fraction must be strictly between zero and one")
        if not isinstance(self.removed_quotas_mmol_per_reference_g, Mapping):
            raise ValueError("removed storage quotas must be a product mapping")
        products = _selected_products(tuple(self.removed_quotas_mmol_per_reference_g))
        if not isinstance(self.molecular_masses, Mapping) or set(self.molecular_masses) != set(products):
            raise ValueError("molecular masses must cover exactly the separated storage products")
        quotas = {
            p: _finite_number(self.removed_quotas_mmol_per_reference_g[p], f"{p} quota", nonnegative=True)
            for p in products
        }
        masses = {
            p: _finite_number(self.molecular_masses[p], f"{p} molecular mass", nonnegative=True)
            for p in products
        }
        if any(quotas[p] == 0.0 or masses[p] == 0.0 for p in products):
            raise ValueError("removed storage quotas and molecular masses must be positive")
        removed_mass = _finite_sum(
            (quotas[p] * (masses[p] / 1000.0) for p in products), "removed reference storage mass"
        )
        if not math.isclose(fraction, 1.0 - removed_mass, rel_tol=1e-12, abs_tol=0.0):
            raise ValueError("structural_mass_fraction is inconsistent with the removed storage mass")
        if not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")
        object.__setattr__(self, "structural_mass_fraction", fraction)
        object.__setattr__(self, "removed_quotas_mmol_per_reference_g", _freeze(quotas))
        object.__setattr__(self, "molecular_masses", _freeze(masses))
        object.__setattr__(self, "metadata", _freeze(self.metadata))

    def to_structural_specific(self, value_per_reference_g: float) -> float:
        value = _finite_number(value_per_reference_g, "reference-specific value")
        converted = value / self.structural_mass_fraction
        if not math.isfinite(converted):
            raise ValueError("structural-specific value must remain finite; arithmetic overflow")
        return converted

    def _mass_components(
        self, structural_biomass_g_l: float, pools_mmol_l: Mapping[str, float]
    ) -> tuple[float, dict[str, float]]:
        structural = _finite_number(structural_biomass_g_l, "structural biomass mass", nonnegative=True)
        if not isinstance(pools_mmol_l, Mapping):
            raise ValueError("pool inventories must be a product-to-mmol/L mapping")
        unknown = set(pools_mmol_l) - self.molecular_masses.keys()
        if unknown:
            raise ValueError(f"unknown or unseparated pool products: {sorted(map(repr, unknown))}")
        masses = {
            p: _finite_number(value, f"{p} pool inventory", nonnegative=True)
            * (self.molecular_masses[p] / 1000.0)
            for p, value in pools_mmol_l.items()
        }
        total = _finite_sum((structural, *masses.values()), "total dry mass")
        return total, masses

    def total_dry_mass(
        self, structural_biomass_g_l: float, pools_mmol_l: Mapping[str, float]
    ) -> float:
        return self._mass_components(structural_biomass_g_l, pools_mmol_l)[0]

    def content_mg_per_g_total(
        self, product: str, structural_biomass: float, pools: Mapping[str, float]
    ) -> float:
        if not isinstance(product, str) or product not in self.molecular_masses:
            raise ValueError(f"unknown or unseparated storage product {product!r}")
        total, masses = self._mass_components(structural_biomass, pools)
        if total == 0.0:
            raise ValueError("content is undefined at zero total dry mass")
        return (masses.get(product, 0.0) / total) * 1000.0


def _validate_standard_program(model: cobra.Model) -> dict[cobra.Reaction, float]:
    if set(model.constraints.keys()) != set(model.metabolites.list_attr("id")):
        raise ValueError("custom solver constraints must be installed after the storage basis conversion")
    expected_variables = {
        variable.name for r in model.reactions for variable in (r.forward_variable, r.reverse_variable)
    }
    if set(model.variables.keys()) != expected_variables or any(v.type != "continuous" for v in model.variables):
        raise ValueError("only standard continuous reaction variables support the storage basis conversion")
    for metabolite in model.metabolites:
        constraint = metabolite.constraint
        if constraint.lb != 0.0 or constraint.ub != 0.0:
            raise ValueError("storage conversion requires homogeneous steady-state metabolite balances")
        expected = {
            variable: sign * r.metabolites[metabolite]
            for r in metabolite.reactions
            for variable, sign in ((r.forward_variable, 1.0), (r.reverse_variable, -1.0))
        }
        if any(not math.isfinite(c) for c in expected.values()):
            raise ValueError(f"nonfinite stoichiometry in balance {metabolite.id!r}")
        actual = {v: float(c) for v, c in constraint.expression.as_coefficients_dict().items() if c}
        if actual != expected:
            raise ValueError(f"custom metabolite balance constraint {metabolite.id!r} cannot be converted")
    for r in model.reactions:
        lb, ub = r.bounds
        if math.isnan(lb) or math.isnan(ub) or lb == math.inf or ub == -math.inf or lb > ub:
            raise ValueError(f"reaction {r.id!r} has invalid flux bounds")
        expected = (
            (max(lb, 0.0), max(ub, 0.0)),
            (max(-ub, 0.0), max(-lb, 0.0)),
        )
        for variable, bounds in zip((r.forward_variable, r.reverse_variable), expected):
            actual = (
                -math.inf if variable.lb is None else variable.lb,
                math.inf if variable.ub is None else variable.ub,
            )
            if bounds[1] == math.inf and actual[1] == float_info.max:
                actual = (actual[0], math.inf)
            if not all(math.isclose(a, b, rel_tol=8.0 * float_info.epsilon, abs_tol=0.0) for a, b in zip(actual, bounds)):
                raise ValueError(f"custom solver variable bounds on {r.id!r} cannot be converted")
    coefficients = linear_reaction_coefficients(model)
    expression = sum(c * r.flux_expression for r, c in coefficients.items())
    if (model.objective.expression - expression).expand() != 0:
        raise ValueError("storage conversion requires a linear objective in net reaction fluxes")
    return {r: _finite_number(c, f"objective coefficient for {r.id}") for r, c in coefficients.items()}


def _accumulation_reactions(model: cobra.Model, products: tuple[str, ...]) -> dict[str, tuple[str, ...]]:
    sinks = {}
    for product in products:
        metabolite = model.metabolites.get_by_id(_STORAGE_SPECIES[product])
        ids = []
        for reaction in metabolite.reactions:
            if not reaction.boundary:
                continue
            if reaction.metabolites == {metabolite: -1.0} and reaction.lower_bound >= 0.0:
                ids.append(reaction.id)
            elif reaction.bounds != (0.0, 0.0):
                raise ValueError(f"storage boundary {reaction.id!r} must be a nonnegative unit demand, not an inventory supply")
        sinks[product] = tuple(sorted(ids))
    return sinks


def separate_storage_biomass(
    model: cobra.Model, products: tuple[str, ...] = ("glycogen", "trehalose")
) -> tuple[cobra.Model, StorageBasis]:
    products = _selected_products(products)
    if not isinstance(model, cobra.Model):
        raise ValueError("storage conversion requires a frozen cobra model")
    quotas = {p: biomass_pool_quota(model, p) for p in products}
    species = {p: model.metabolites.get_by_id(_STORAGE_SPECIES[p]) for p in products}
    masses = {p: float(metabolite.formula_weight) for p, metabolite in species.items()}
    fraction = 1.0 - _finite_sum(
        (quotas[p] * (masses[p] / 1000.0) for p in products), "removed reference storage mass"
    )
    if not 0.0 < fraction < 1.0:
        raise ValueError("removed storage mass must leave a positive structural_mass_fraction below one")
    coefficients = _validate_standard_program(model)
    accumulation = _accumulation_reactions(model, products)
    pool = model.metabolites.get_by_id(_PROTEIN_POOL)
    enzyme_draws = tuple(sorted(
        r.id for r in pool.reactions
        if r.metabolites[pool] < 0.0
        and len(r.metabolites) == 2
        and all(c == 1.0 for m, c in r.metabolites.items() if m != pool)
    ))
    pool_supply = tuple(sorted(r.id for r in pool.reactions if r.metabolites[pool] > 0.0))
    converted_ids = tuple(r.id for r in model.reactions if r.id not in _UNCHANGED_RATES)
    transformed_id = f"{model.id}__structural_storage__{'__'.join(products)}"
    basis = StorageBasis(
        original_model_id=model.id,
        structural_mass_fraction=fraction,
        removed_quotas_mmol_per_reference_g=quotas,
        molecular_masses=masses,
        metadata={
            "mode": "structural_biomass_variable_storage",
            "transformed_model_id": transformed_id,
            "species": {
                p: {
                    "metabolite_id": metabolite.id,
                    "formula": metabolite.formula,
                    "molecular_mass_units": "g/mol (numerically mg/mmol)",
                    "mass_convention": (
                        "frozen-model glucose equivalents, not actual anhydroglucose polymer residue mass"
                        if p == "glycogen" else "frozen-model trehalose species formula mass"
                    ),
                }
                for p, metabolite in species.items()
            },
            "reference_mass_convention": "One original biomass drain unit is one reference gram; the remaining mass is 1 minus the selected model-species reserve masses, not a refit of all biomass chemistry.",
            "structural_mass_fraction_equation": "f = 1 - sum(quota_mmol_per_reference_g[p] * molecular_mass_g_per_mol[p] / 1000)",
            "total_dry_mass_equation": "X_total_g_l = X_structural_g_l + sum(pool_mmol_l[p] * molecular_mass_g_per_mol[p] / 1000)",
            "unchanged_rate_reaction_ids": _UNCHANGED_RATES,
            "converted_bound_reaction_ids": converted_ids,
            "model_bounds_already_converted": True,
            "specific_amount_conversion_factor": 1.0 / fraction,
            "specific_to_reference_conversion_factor": fraction,
            "enzyme_amount_reaction_ids": enzyme_draws,
            "enzyme_amount_units": "mmol enzyme/g structural biomass; not mmol/g/h",
            "protein_pool_supply_reaction_ids": pool_supply,
            "protein_pool_units": "g protein/g structural biomass",
            "metabolic_flux_units": "mmol/g structural biomass/h",
            "growth_rate_units": "1/h",
            "biomass_reaction_id": "r_2111",
            "stoichiometry_convention": "Remove selected reserve consumption from r_4048 and divide its other inputs by f; leave its carbohydrate output unchanged. Divide every r_4041 coefficient except carbohydrate input and biomass output by f. All other reaction stoichiometry, including enzyme kcat costs, is unchanged.",
            "energy_convention": "All original growth-associated ATP, water, ADP, proton and phosphate terms are retained per reference-equivalent structural growth and divided by f. No energy refund is attributed to reserve removal. Native reserve biosynthetic costs remain; non-growth maintenance bounds are divided by f.",
            "objective_convention": "Existing linear net-flux objectives are reexpressed with v_reference = f*v_structural for converted reactions; unchanged-rate coefficients, including growth, retain their original values.",
            "external_values_requiring_conversion": (
                "reference-specific uptake, secretion, product and maintenance rates: divide by f before setting new bounds",
                "reference-specific enzyme mmol/g and protein g/g budgets, including separate proteomics provenance tables: divide by f; copied model bounds are already converted",
                "reference-composition inoculum: X_structural = f*X_reference and pools[p] = quota[p]*X_reference only when that initial composition is explicitly assumed",
            ),
            "values_not_requiring_conversion": (
                "specific growth rates in 1/h",
                "already structural-specific rates and budgets",
                "volumetric concentrations and inventories, time, kcat, molecular masses, dimensionless multipliers",
            ),
            "accumulation_reaction_ids": accumulation,
            "inventory_convention": "Inventories are separate nonnegative mmol/L states; omitted pool entries mean zero, never an implicit biomass quota. Only selected pools are outside structural biomass. Existing forward unit demands measure net inventory addition, not gross synthesis.",
            "initial_mass_accounting": "Given total initial dry mass and explicit inventories, structural mass is total mass minus the sum of inventory species masses, and must be nonnegative. Zero inventories therefore make structural and total mass identical. Multiplication by f applies only to an explicitly assumed original reference composition.",
            "inventory_balance_convention": "At constant volume dX_structural/dt = mu*X_structural and dP[p]/dt = v_net[p]*X_structural. No growth-dilution term is subtracted from volumetric inventory; P/X_structural instead has dilution -mu*(P/X_structural).",
            "reference_equivalence": "For balanced original composition impose each inventory demand = quota[p]*mu/f, with both pools when both quotas are removed. The remaining specific fluxes and budgets map by 1/f, while the four unchanged rates retain their values.",
            "inventory_remobilization_supported": False,
            "native_turnover_convention": "Native turnover reactions remain with unchanged stoichiometry and converted specific bounds. No inventory withdrawal, signed demand, free-carbon source, loss kinetics or time integrator is installed.",
            "solver_extension_convention": "Only standard continuous reaction variables and homogeneous stoichiometric balances are supported; add custom coupling constraints after the conversion in compatible units.",
        },
    )
    out = model.copy()
    out.notes = deepcopy(model.notes)
    out.annotation = deepcopy(model.annotation)
    for collection in (out.metabolites, out.reactions, out.genes):
        for entity in collection:
            entity.notes = deepcopy(entity.notes)
            entity.annotation = deepcopy(entity.annotation)
    carbohydrate = out.reactions.get_by_id("r_4048")
    removed_ids = {species[p].id for p in products}
    carbohydrate.add_metabolites({
        metabolite: 0.0 if metabolite.id in removed_ids else basis.to_structural_specific(coefficient)
        for metabolite, coefficient in carbohydrate.metabolites.items()
        if metabolite.id != _CARBOHYDRATE
    }, combine=False)
    assembly = out.reactions.get_by_id("r_4041")
    assembly.add_metabolites({
        metabolite: basis.to_structural_specific(coefficient)
        for metabolite, coefficient in assembly.metabolites.items()
        if metabolite.id not in {_CARBOHYDRATE, _BIOMASS}
    }, combine=False)
    for rid in converted_ids:
        reaction = out.reactions.get_by_id(rid)
        reaction.bounds = tuple(
            bound if math.isinf(bound) else basis.to_structural_specific(bound)
            for bound in reaction.bounds
        )
    out.objective = {
        out.reactions.get_by_id(r.id): c if r.id in _UNCHANGED_RATES else c * fraction
        for r, c in coefficients.items()
    }
    out.objective.direction = model.objective.direction
    out.id = transformed_id
    out.notes.update({
        "ystwin_storage_basis": "structural_biomass_variable_storage",
        "ystwin_reference_model_id": model.id,
        "ystwin_structural_mass_fraction": repr(fraction),
        "ystwin_separated_storage_products": ", ".join(products),
    })
    return out, basis
