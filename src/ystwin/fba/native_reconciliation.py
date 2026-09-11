from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from decimal import Decimal
import math
from numbers import Real

import cobra
from cobra.util.array import create_stoichiometric_matrix
import numpy as np

from ..analysis.native_physiology import load_native_chemostat_data
from .dynamic_rates import _validated_optimize
from .solver import PINNED_TOLERANCE, configure


_PRIMARY_URL = "https://repository.tudelft.nl/file/File_0981d011-0515-40cd-bc5f-b5004cd9ca42"
_PRIMARY_TOKENS = {
    "GlucoseUptake": "0.3 0.6 1.1 1.7 2.3 2.8 3.4 4.5 8.6 11.1".split(),
    "O2uptake": "0.8 1.3 2.5 3.9 5.3 7.0 7.4 6.1 5.1 3.7".split(),
    "CO2production": "0.8 1.4 2.7 4.2 5.7 7.5 8.0 8.8 14.9 18.9".split(),
    "Ethanol": "0 0 0 0 0 0 0.11 2.3 9.5 13.9".split(),
    "Acetate": "0 0 0 0 0 0 0.08 0.41 0.62 0.60".split(),
    "Pyruvate": "0 0 0 0 0 0 0.01 0.01 0.03 0.05".split(),
    "Glycerol": "0 0 0 0 0 0 0 0 0.05 0.15".split(),
    "BiomassYield_g_per_g": "0.45 0.47 0.48 0.49 0.48 0.48 0.46 0.37 0.23 0.20".split(),
    "CarbonRecovery_pct": "98.9 95.0 96.0 102.4 100.9 102.6 97.0 99.1 99.4 97.9".split(),
}


def _number(value, label):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite real number")
    return float(value)


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be explicit and nonempty")
    return value


def bound_metadata(value):
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or math.isnan(value):
        raise ValueError("metadata bound must be real and cannot be NaN")
    if math.isinf(value):
        return {"kind": "unbounded", "direction": "positive" if value > 0 else "negative"}
    return float(value)


@dataclass(frozen=True)
class FluxObservable:
    key: str
    unit: str
    coefficients: Mapping[str, float]
    source: str

    def __post_init__(self):
        for field in ("key", "unit", "source"):
            _text(getattr(self, field), field)
        if not isinstance(self.coefficients, Mapping) or not self.coefficients:
            raise ValueError("observable requires explicit nonempty reaction coefficients")
        coefficients = {_text(k, "reaction id"): _number(v, "flux coefficient")
                        for k, v in self.coefficients.items()}
        if any(v == 0 for v in coefficients.values()):
            raise ValueError("observable coefficients must be nonzero")
        object.__setattr__(self, "coefficients", coefficients)


@dataclass(frozen=True)
class ObservationConstraint:
    observable: FluxObservable
    lower: float | None
    upper: float | None
    source: str

    def __post_init__(self):
        if not isinstance(self.observable, FluxObservable):
            raise ValueError("constraint requires an explicit physical observable")
        _text(self.source, "constraint provenance")
        for name in ("lower", "upper"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _number(value, name))
        if self.lower is None and self.upper is None:
            raise ValueError("an observation constraint needs a finite bound")
        if self.lower is not None and self.upper is not None and self.lower > self.upper:
            raise ValueError("observation lower bound exceeds upper bound")


def boundary_observable(model, *, key, metabolite_id, reaction_ids, direction, formula, source):
    if direction not in ("uptake", "production"):
        raise ValueError("boundary direction must be uptake or production")
    if (not isinstance(reaction_ids, Sequence) or isinstance(reaction_ids, (str, bytes))
            or not reaction_ids or len(set(reaction_ids)) != len(reaction_ids)):
        raise ValueError("boundary roles require unique explicit reaction ids")
    species = model.metabolites.get_by_id(metabolite_id)
    if not formula or species.formula != formula:
        raise ValueError(f"{metabolite_id}: molecular formula disagrees with the declared role")
    coefficients = {}
    for rid in reaction_ids:
        reaction = model.reactions.get_by_id(rid)
        if len(reaction.metabolites) != 1 or species not in reaction.metabolites:
            raise ValueError(f"{rid}: declared role is not a one-species boundary for {metabolite_id}")
        coefficient = _number(reaction.metabolites[species], "boundary stoichiometry")
        coefficients[rid] = coefficient if direction == "uptake" else -coefficient
    for reaction in species.reactions:
        if (len(reaction.metabolites) == 1 and reaction.bounds != (0.0, 0.0)
                and reaction.id not in coefficients):
            raise ValueError(f"every active boundary for {metabolite_id} needs a role: {reaction.id}")
    return FluxObservable(key, "mmol/gDW/h", coefficients, source)


def biomass_flux_equivalence(model, *, growth_reaction_id, assembly_reaction_id, biomass_metabolite_id):
    if growth_reaction_id == assembly_reaction_id:
        raise ValueError("growth and assembly require distinct roles")
    species = model.metabolites.get_by_id(biomass_metabolite_id)
    growth = model.reactions.get_by_id(growth_reaction_id)
    assembly = model.reactions.get_by_id(assembly_reaction_id)
    consumed = growth.get_coefficient(species.id)
    produced = assembly.get_coefficient(species.id)
    if consumed >= 0 or produced <= 0:
        raise ValueError("biomass role stoichiometry must identify production and consumption")
    closed = []
    for reaction in species.reactions:
        if reaction.id in (growth_reaction_id, assembly_reaction_id):
            continue
        if reaction.bounds != (0.0, 0.0):
            raise ValueError(f"unaccounted active biomass source or sink: {reaction.id}")
        closed.append(reaction.id)
    return {"growth_reaction_id": growth_reaction_id, "assembly_reaction_id": assembly_reaction_id,
            "biomass_metabolite_id": biomass_metabolite_id,
            "assembly_flux_per_growth_flux": -consumed / produced,
            "closed_other_reaction_ids": sorted(closed)}


def ec_native_observables(model):
    source = "Explicit species/reaction roles checked against native EC SBML stoichiometry and formulas"
    roles = {
        "GlucoseUptake": ("s_0565[e]", ("r_1714", "r_1714_REV"), "uptake", "C6H12O6"),
        "O2uptake": ("s_1277[e]", ("r_1992", "r_1992_REV"), "uptake", "O2"),
        "CO2production": ("s_0458[e]", ("r_1672",), "production", "CO2"),
        "Ethanol": ("s_0681[e]", ("r_1761",), "production", "C2H6O"),
        "Glycerol": ("s_0766[e]", ("r_1808",), "production", "C3H8O3"),
        "Acetate": ("s_0364[e]", ("r_1634",), "production", "C2H3O2"),
        "Pyruvate": ("s_1400[e]", ("r_2033",), "production", "C3H3O3"),
    }
    result = {
        key: boundary_observable(model, key=key, metabolite_id=mid, reaction_ids=rids,
                                 direction=direction, formula=formula, source=source)
        for key, (mid, rids, direction, formula) in roles.items()
    }
    biomass = model.reactions.get_by_id("r_2111")
    if {m.id: c for m, c in biomass.metabolites.items()} != {"s_0450[c]": -1.0}:
        raise ValueError("native biomass drain role disagrees with the SBML")
    result["growth"] = FluxObservable("growth", "1/h", {"r_2111": 1.0},
                                      "Native SBML biomass drain; model dry-mass normalization")
    return result


def native_conditions():
    data = load_native_chemostat_data()
    result = []
    for index, (growth, rows) in enumerate(data.observations.groupby("growth_rate_per_h", sort=True)):
        observations = {}
        for record in rows.to_dict("records"):
            key = record["observable_id"]
            token = _PRIMARY_TOKENS[key][index]
            if float(token) != record["reported_value"]:
                raise ValueError("primary printed token disagrees with the verified native TSV")
            quantified = record["observation_status"] == "quantified"
            half_width = float(Decimal(10) ** Decimal(token).as_tuple().exponent / 2) if quantified else None
            observations[key] = {
                "reported_value": float(record["reported_value"]),
                "value": float(record["value"]) if quantified else None,
                "observation_status": record["observation_status"], "unit": record["unit"],
                "detection_limit": None, "source_line": int(record["source_line"]),
                "primary_printed_token": token, "rounding_half_width": half_width,
                "rounding_is_confidence_interval": False,
            }
        result.append({
            "condition_id": f"vanHoek1998:D={growth:g}", "growth_rate_per_h": float(growth),
            "evidence_scope": ("native_development_training" if rows.iloc[0]["split"] == "train"
                               else "native_development_interpolation"),
            "observations": observations, "source": data.metadata["source"],
            "precision_provenance": {
                "url": _PRIMARY_URL, "page": 4227, "table": 1,
                "verification": "Printed Table 1 tokens rechecked directly in the primary PDF; not TSV or pandas decimal formatting",
                "interpretation": "Half a last printed decimal under a nearest-rounding assumption; not confidence bounds or a measured SD",
                "growth": "Dilution rate remains the imposed reported setpoint, not an inferred uncertainty range",
            },
        })
    return result


def condition_constraints(condition, observables, *, mode, impose_growth=False):
    if mode not in ("uptake_caps", "quantified_point", "printed_rounding", "rounded_uptake_caps"):
        raise ValueError("unknown native constraint mode")
    constraints = []
    for key, record in condition["observations"].items():
        if key not in observables or record["value"] is None:
            continue
        value = record["value"]
        if mode in ("uptake_caps", "rounded_uptake_caps"):
            if key not in ("GlucoseUptake", "O2uptake"):
                continue
            lower, upper = 0.0, value
            if mode == "rounded_uptake_caps":
                upper += record["rounding_half_width"]
        elif mode == "printed_rounding":
            half = record["rounding_half_width"]
            lower, upper = value - half, value + half
        else:
            lower = upper = value
        constraints.append(ObservationConstraint(
            observables[key], lower, upper, f"{condition['condition_id']}; {mode}; primary Table 1",
        ))
    if impose_growth:
        growth = condition["growth_rate_per_h"]
        constraints.append(ObservationConstraint(observables["growth"], growth, growth,
                                                  "Native reported chemostat dilution rate"))
    return tuple(constraints)


def _expression(model, observable):
    return sum(c * model.reactions.get_by_id(rid).flux_expression
               for rid, c in observable.coefficients.items())


def _observed(expression, primal):
    return math.fsum(float(coefficient) * (1.0 if term == 1 else primal[term.name])
                     for term, coefficient in expression.as_coefficients_dict().items())


def _bound_residual(value, lower, upper):
    value = _number(value, "solved physical value")
    lower_residual = value - lower if lower is not None and math.isfinite(lower) else None
    upper_residual = upper - value if upper is not None and math.isfinite(upper) else None
    violation = max(0.0, -lower_residual if lower_residual is not None else 0.0,
                    -upper_residual if upper_residual is not None else 0.0)
    return {"value": value, "lower": bound_metadata(lower), "upper": bound_metadata(upper),
            "lower_residual": lower_residual, "upper_residual": upper_residual,
            "violation": violation}


def _audit(model, stoichiometry, fluxes):
    primal = model.solver.primal_values
    mass = np.asarray(stoichiometry @ fluxes).ravel()
    reactions = {r.id: _bound_residual(float(v), *r.bounds) for r, v in zip(model.reactions, fluxes)}
    variables = {v.name: _bound_residual(primal[v.name], v.lb, v.ub) for v in model.variables}
    constraints = {c.name: _bound_residual(_observed(c.expression, primal), c.lb, c.ub)
                   for c in model.constraints}
    summary = {
        "max_abs_mass_balance": float(np.max(np.abs(mass), initial=0.0)),
        "max_reaction_bound_violation": max((v["violation"] for v in reactions.values()), default=0.0),
        "max_variable_bound_violation": max((v["violation"] for v in variables.values()), default=0.0),
        "max_constraint_violation": max((v["violation"] for v in constraints.values()), default=0.0),
        "n_mass_balance_rows": len(mass), "n_reactions": len(reactions),
        "n_variables": len(variables), "n_constraints": len(constraints),
        "physical_acceptance_tolerance": PINNED_TOLERANCE,
        "row_units": "Original unscaled native SBML rows; chemical rows in mmol/gDW/h, pseudo-resource rows in their native resource units",
    }
    if not np.isfinite(mass).all() or any(summary[k] > PINNED_TOLERANCE for k in (
        "max_abs_mass_balance", "max_reaction_bound_violation",
        "max_variable_bound_violation", "max_constraint_violation",
    )):
        raise RuntimeError("full original-unit native residual audit failed")
    return {"summary": summary, "mass_balance": {m.id: float(v) for m, v in zip(model.metabolites, mass)},
            "reaction_bounds": reactions, "variables": variables, "constraints": constraints}


def _prepare(model, constraints, timeout_s):
    if not isinstance(model, cobra.Model):
        raise ValueError("native reconciliation requires a cobra model")
    if isinstance(timeout_s, bool) or not isinstance(timeout_s, int) or not 1 <= timeout_s <= 30:
        raise ValueError("per-solve timeout must be an integer from 1 to 30 seconds")
    if any(not isinstance(c, ObservationConstraint) for c in constraints):
        raise ValueError("native constraints must be explicit ObservationConstraint records")
    keys = [c.observable.key for c in constraints]
    if len(set(keys)) != len(keys):
        raise ValueError("native observable constraints must have unique keys")
    out = model.copy()
    configure(out)
    out.solver.configuration.timeout = timeout_s
    return out


def _install(out, constraints, relaxations=None, slack=None):
    for index, constraint in enumerate(constraints):
        expression = _expression(out, constraint.observable)
        scales = (relaxations or {}).get(constraint.observable.key, (0.0, 0.0))
        for side, bound, scale in (("lower", constraint.lower, scales[0]),
                                    ("upper", constraint.upper, scales[1])):
            if bound is None:
                continue
            name = f"native_reconciliation_{index}_{side}"
            if name in out.constraints:
                raise ValueError("native reconciliation constraint identifier collision")
            adjusted = expression
            if slack is not None:
                adjusted += (scale if side == "lower" else -scale) * slack
            out.add_cons_vars([out.problem.Constraint(
                adjusted, lb=bound if side == "lower" else None,
                ub=bound if side == "upper" else None, name=name,
            )])


def _solve(out, constraints, objective, direction):
    if direction not in ("min", "max"):
        raise ValueError("objective direction must be min or max")
    out.objective = out.problem.Objective(objective, direction=direction)
    stoichiometry = create_stoichiometric_matrix(out, array_type="lil")
    result = {"status": None, "reason": None, "objective_value": None, "fluxes": None,
              "observables": None, "audit": None,
              "observation_constraints": [asdict(c) for c in constraints],
              "solver": {"interface": out.solver.interface.__name__,
                         "physical_tolerance": PINNED_TOLERANCE,
                         "timeout_seconds": out.solver.configuration.timeout,
                         "policy": "existing _validated_optimize; no tolerance or global policy changes"}}
    try:
        answer = _validated_optimize(out, stoichiometry, "native reconciliation", allow_infeasible=True)
        if answer is None:
            result.update(status="infeasible", reason="infeasible under original hard bounds and the explicitly recorded scenario")
            return result
        value, fluxes = answer
        value = _number(value, "physical objective")
        if not np.isfinite(fluxes).all():
            raise RuntimeError("nonfinite physical reaction fluxes")
        audit = _audit(out, stoichiometry, fluxes)
        result.update(status="optimal", objective_value=value,
                      fluxes={r.id: float(v) for r, v in zip(out.reactions, fluxes)},
                      observables={c.observable.key: _observed(_expression(out, c.observable), out.solver.primal_values)
                                   for c in constraints}, audit=audit)
    except (RuntimeError, ValueError) as exc:
        result.update(status="failed", reason=str(exc))
    return result


def solve_native(model, *, constraints, objective=None, direction="max", timeout_s=20):
    constraints = tuple(constraints)
    out = _prepare(model, constraints, timeout_s)
    _install(out, constraints)
    expression = _expression(out, objective) if objective is not None else 0.0
    result = _solve(out, constraints, expression, direction)
    if result["status"] == "optimal" and objective is not None:
        result["observables"][objective.key] = math.fsum(
            c * result["fluxes"][rid] for rid, c in objective.coefficients.items())
    return result


def minimum_relaxation(model, *, constraints, relaxations, max_relaxation, relaxation_unit, timeout_s=20):
    constraints = tuple(constraints)
    maximum = _number(max_relaxation, "maximum diagnostic relaxation")
    if maximum <= 0:
        raise ValueError("maximum diagnostic relaxation must be positive and finite")
    _text(relaxation_unit, "diagnostic objective unit")
    if not isinstance(relaxations, Mapping) or not relaxations:
        raise ValueError("explicit nonempty observation relaxation scales are required")
    keyed = {c.observable.key: c for c in constraints}
    if not set(relaxations) <= set(keyed):
        raise ValueError("only declared observation constraints can be relaxed")
    scales = {}
    for key, pair in relaxations.items():
        if not isinstance(pair, Sequence) or len(pair) != 2:
            raise ValueError("relaxation requires lower and upper physical-unit scales")
        lower, upper = (_number(v, "relaxation scale") for v in pair)
        if min(lower, upper) < 0 or max(lower, upper) <= 0:
            raise ValueError("relaxation scales must be nonnegative with at least one positive side")
        if (lower and keyed[key].lower is None) or (upper and keyed[key].upper is None):
            raise ValueError("cannot relax an absent observation bound")
        scales[key] = (lower, upper)
    out = _prepare(model, constraints, timeout_s)
    slack_id = "native_reconciliation_discrepancy"
    if slack_id in out.variables:
        raise ValueError("native reconciliation diagnostic variable collision")
    slack = out.problem.Variable(slack_id, lb=0.0, ub=maximum)
    out.add_cons_vars([slack])
    _install(out, constraints, scales, slack)
    result = _solve(out, constraints, slack, "min")
    result.update(minimum_relaxation=None, relaxation_unit=relaxation_unit,
                  max_relaxation=maximum, relaxation_scales=scales, adjustments=None,
                  interpretation="Bounded minimum discrepancy diagnostic; not a replacement of measurements, a parameter fit, or a confidence interval")
    if result["status"] == "optimal":
        amount = _number(slack.primal, "physical diagnostic slack")
        adjustments = {}
        for key, pair in scales.items():
            constraint, value = keyed[key], result["observables"][key]
            adjustments[key] = {
                "value": value, "unit": constraint.observable.unit,
                "lower_excess": max(0.0, constraint.lower - value) if constraint.lower is not None else None,
                "upper_excess": max(0.0, value - constraint.upper) if constraint.upper is not None else None,
                "allowed_lower_adjustment": pair[0] * amount,
                "allowed_upper_adjustment": pair[1] * amount,
            }
        result.update(minimum_relaxation=amount, adjustments=adjustments)
    return result


class NativeConsistencyError(ValueError):
    pass


def require_native_consistency(result, *, growth_key, target_growth_per_h):
    target = _number(target_growth_per_h, "native target growth")
    if target < 0:
        raise ValueError("native growth cannot be negative")
    if result["status"] != "optimal":
        raise NativeConsistencyError(f"native solve did not establish consistency: {result['status']}")
    value = _number(result["observables"][growth_key], "native achievable growth")
    if target - value > PINNED_TOLERANCE:
        raise NativeConsistencyError(f"native growth shortfall {target - value:.12g} 1/h; target {target:g}, achievable {value:.12g}")
    return value


@dataclass(frozen=True)
class EnergyRoles:
    ngam_reaction_id: str
    biomass_assembly_reaction_id: str
    atp_metabolite_id: str
    water_metabolite_id: str
    adp_metabolite_id: str
    proton_metabolite_id: str
    phosphate_metabolite_id: str

    def __post_init__(self):
        values = asdict(self)
        for key, value in values.items():
            _text(value, key)
        if len(set(values.values())) != len(values):
            raise ValueError("energy roles must have distinct explicit identifiers")


def energy_diagnostic(model, *, roles, source, gam_mmol_atp_per_gdw=None, ngam_mmol_atp_per_gdw_h=None):
    if not isinstance(roles, EnergyRoles):
        raise ValueError("explicit stoichiometric energy roles are required")
    _text(source, "diagnostic parameter provenance")
    roles_and_formulas = ((roles.atp_metabolite_id, "C10H12N5O13P3", -1.0),
                         (roles.water_metabolite_id, "H2O", -1.0),
                         (roles.adp_metabolite_id, "C10H12N5O10P2", 1.0),
                         (roles.proton_metabolite_id, "H", 1.0),
                         (roles.phosphate_metabolite_id, "HO4P", 1.0))
    hydrolysis = {}
    for mid, formula, coefficient in roles_and_formulas:
        if model.metabolites.get_by_id(mid).formula != formula:
            raise ValueError("ATP hydrolysis species formula does not match the declared role")
        hydrolysis[mid] = coefficient
    maintenance = model.reactions.get_by_id(roles.ngam_reaction_id)
    if {m.id: c for m, c in maintenance.metabolites.items()} != hydrolysis:
        raise ValueError("NGAM stoichiometry is not the complete declared ATP hydrolysis vector")
    assembly = model.reactions.get_by_id(roles.biomass_assembly_reaction_id)
    original_gam = -assembly.get_coefficient(roles.atp_metabolite_id)
    if original_gam < 0 or any(not math.isclose(assembly.get_coefficient(mid), original_gam * coefficient,
                                                rel_tol=0.0, abs_tol=1e-12)
                                for mid, coefficient in hydrolysis.items()):
        raise ValueError("biomass energy stoichiometry is not a consistent hydrolysis multiple")
    out = model.copy()
    if gam_mmol_atp_per_gdw is not None:
        gam = _number(gam_mmol_atp_per_gdw, "GAM mmol ATP/gDW")
        if gam < 0:
            raise ValueError("GAM must be nonnegative")
        out.reactions.get_by_id(assembly.id).add_metabolites({
            out.metabolites.get_by_id(mid): coefficient * (gam - original_gam)
            for mid, coefficient in hydrolysis.items()
        })
    else:
        gam = original_gam
    if ngam_mmol_atp_per_gdw_h is not None:
        ngam = _number(ngam_mmol_atp_per_gdw_h, "NGAM mmol ATP/gDW/h")
        if ngam < 0:
            raise ValueError("NGAM must be nonnegative")
        out.reactions.get_by_id(maintenance.id).bounds = (ngam, ngam)
    return out, {
        "kind": "explicit_native_energy_diagnostic", "calibrated": False, "source": source,
        "roles": asdict(roles), "original_gam_mmol_atp_per_gdw": original_gam,
        "gam_mmol_atp_per_gdw": gam,
        "original_ngam_bounds_mmol_atp_per_gdw_h": [bound_metadata(v) for v in maintenance.bounds],
        "ngam_bounds_mmol_atp_per_gdw_h": [bound_metadata(v) for v in out.reactions.get_by_id(maintenance.id).bounds],
        "hydrolysis_vector": hydrolysis,
        "identifiability": "At a fixed biomass assembly flux, GAM and NGAM contribute through GAM*v_assembly + v_NGAM; one condition cannot separate them",
        "scope": "Only the declared energy vector or NGAM bound pair is changed, on a model copy. This is not an adopted physiological calibration.",
    }
