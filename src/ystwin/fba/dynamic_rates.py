"""The closure `simulate_batch` integrates: GEM sets growth, enzymes set product.

`fba/dynamic.py` deliberately decides no biology -- it integrates a
``f(substrate) -> (mu, q_S, q_P)`` it is handed. This module builds that function, and in
doing so it is where four layers stop auditing and start setting the number:

    substrate  -> uptake kinetics                    -> q_S
    q_S        -> GEM, solved under a stress-constrained NGAM  -> mu
    stress     -> regulation -> enzyme abundance     -> [E]
    [E], kcat  -> enzyme capacity                    -> q_P

In the steady-state chain none of these could reach the answer: `mu` came from the pump and
the product flux came from a fitted scalar. Once the vessel is a batch, `mu` is computed at
every step and the titre is its integral, so the GEM and everything constraining it are in
the chain by construction rather than by permission.

WHAT IS DIFFERENT ABOUT q_P HERE, and it is the whole point. The steady-state chain gets its
entry flux from ``v_in = alpha * expression``, where alpha is FITTED on one product's own
measured states -- which is why the chain cannot answer for a product nobody has calibrated.
Here ``q_P = kcat * [E] * saturation``, and neither kcat nor [E] is fitted to the product:
kcat is a property of a protein and [E] comes from the proteome. That is the difference
between a model that needs a calibration per product and one that needs an enzyme.

WHAT IT COSTS, measured rather than asserted. The capacity route was checked against native
enzymes where kcat, abundance and flux are each known independently (see
`data/proteome/SOURCE.md`): it lands within 6% for one enzyme and short by up to twentyfold
for others, and the error is ONE-SIDED LOW because in vitro turnover understates the in vivo
rate. **So a titre from this module carries about an order of magnitude, biased low.** It is
not a replacement for a fitted alpha where one exists; it is an answer where none does.

THE REGULATION ROUTE IS ABUNDANCE, NOT CEILINGS. `bridge/regulation.py`'s E-Flux formulation
scales reaction upper bounds, and this repository has measured that it never binds: every
stress-touched module is INDUCED, and raising a ceiling above a flux already below it does
nothing. The route that works is the other one -- expression scales how much enzyme there
is, and enzyme is what ``kcat * [E]`` is made of. That is why `regulation_scale` multiplies
an abundance here and not a bound.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from swiglpk import GLP_DUAL, GLP_PRIMAL

from .physiology import _uptake_boundary, cap_uptake
from .solver import PINNED_TOLERANCE, configure, growth_or_none

__all__ = ["CoupledRates", "RateClosure", "NetworkRateClosure", "build_rates", "build_network_rates"]

#: Resting non-growth-associated maintenance, mmol ATP/gDW/h, and the stress-driven demand on
#: top of it. Both are `bridge/latent_bridge.py`'s, not new constants -- imported by value
#: here so that a reader of this module can see what a stress of 1.0 actually does, and
#: changed in one place if that module's bound moves.
RESTING_NGAM = 0.7
MAX_STRESS_NGAM = 6.5

_GLPK_WORKING_FEASIBILITY = 1e-10
_GLPK_SIMPLEX_ITERATIONS = 10000
_COUPLED_SOLVER_METHOD = (
    "pinned GLPK scaled warm-start primal simplex with dual recovery on nonterminal status; "
    f"presolve off; at most {_GLPK_SIMPLEX_ITERATIONS} iterations per attempt; full physical validation"
)


@dataclass(frozen=True)
class RateClosure:
    """What was used to build the closure, kept so a result can state its own inputs.

    Carried rather than discarded because a titre computed from a kcat and an abundance is
    only interpretable beside them, and a number whose provenance has to be reconstructed
    from a call site is a number nobody can check.
    """

    kcat_per_s: float
    enzyme_mmol_per_gdcw: float
    regulation_scale: float
    stress: float
    ngam: float
    uptake_vmax_mmol_per_gdcw_h: float
    uptake_km_mmol_per_l: float
    product_reaction_id: str | None = None
    biomass_reaction_id: str | None = None
    growth_fraction: float = 0.9
    growth_retention: float = 1.0
    maximum_growth_per_h: float | None = None
    allocation_assumption: str = "fixed fraction of attainable growth; maximize product; minimize uptake"
    solver_method: str = _COUPLED_SOLVER_METHOD
    solver_working_feasibility_tolerance: float = _GLPK_WORKING_FEASIBILITY
    physical_validation_tolerance: float = PINNED_TOLERANCE

    @property
    def q_product_max_mmol_per_gdcw_h(self) -> float:
        """The ceiling on product formation: kcat times the enzyme actually present."""
        return self.kcat_per_s * 3600.0 * self.enzyme_mmol_per_gdcw * self.regulation_scale

    def summary(self) -> str:
        return (f"kcat {self.kcat_per_s:.3g}/s x [E] {self.enzyme_mmol_per_gdcw:.3g} mmol/gDCW "
                f"x regulation {self.regulation_scale:.3g} -> q_P_max "
                f"{self.q_product_max_mmol_per_gdcw_h:.4g} mmol/gDCW/h; stress {self.stress:.2f} "
                f"-> NGAM {self.ngam:.2f}")


@dataclass(frozen=True)
class NetworkRateClosure:
    substrate_exchange_id: str
    oxygen_exchange_id: str
    ngam_reaction_id: str
    product_reaction_id: str
    biomass_reaction_id: str
    ngam: float
    uptake_vmax_mmol_per_gdcw_h: float
    uptake_km_mmol_per_l: float
    oxygen_lower_bound: float
    growth_fraction: float
    growth_retention: float
    maximum_growth_per_h: float | None
    allocation_fraction: float
    fixed_growth_per_h: float | None
    attainable_growth_per_h: float
    growth_limit_per_h: float
    minimum_feasible_uptake_mmol_per_gdcw_h: float
    allocation_assumption: str = (
        "growth floor gamma * attainable; q_P = q_min + alpha * (q_max - q_min); "
        "fix product, maximize growth, fix growth, minimize actual carbon uptake")
    product_rate_basis: str = "allocation-conditioned, not measured"
    enzyme_capacity_basis: str = "supplied network enzyme constraints and shared protein budgets"
    interpolation_assumption: str = "piecewise-linear mixtures of feasible coupled flux states"
    solver_method: str = _COUPLED_SOLVER_METHOD
    solver_working_feasibility_tolerance: float = _GLPK_WORKING_FEASIBILITY
    physical_validation_tolerance: float = PINNED_TOLERANCE

    @property
    def physiology_conditioned(self) -> bool:
        return self.fixed_growth_per_h is not None

    def summary(self) -> str:
        growth = (f"fixed measured growth {self.fixed_growth_per_h:g}/h"
                  if self.physiology_conditioned else f"growth floor fraction {self.growth_fraction:g}")
        return (f"network enzyme budgets; {growth}; product allocation {self.allocation_fraction:g} "
                f"of [q_min, q_max]; q_P {self.product_rate_basis}; NGAM {self.ngam:g}")


@dataclass(frozen=True)
class CoupledRates:
    uptake_grid: np.ndarray
    grid_fluxes: np.ndarray
    feasible: np.ndarray
    substrate_km_mmol_per_l: float
    product_reaction_id: str
    statuses: tuple[str, ...]
    minimum_feasible_uptake_mmol_per_gdcw_h: float = 0.0
    _solve_at_uptake: Callable[[float], tuple[float, float, float] | None] | None = field(
        default=None, repr=False, compare=False)
    _fixed_growth_per_h: float | None = field(default=None, repr=False)
    _exact_cache: dict[float, tuple[float, float, float] | None] = field(
        default_factory=dict, repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in ("uptake_grid", "grid_fluxes", "feasible"):
            values = np.array(getattr(self, name), copy=True)
            values.setflags(write=False)
            object.__setattr__(self, name, values)

    def _infeasible_rates(self, allowed: float) -> tuple[float, float, float]:
        if self._fixed_growth_per_h is not None:
            raise ValueError(
                f"fixed_growth_per_h={self._fixed_growth_per_h:g} is infeasible at "
                f"allowed molecular uptake {allowed:g}")
        return 0.0, 0.0, 0.0

    def at_uptake(self, allowed_molecular_uptake: float) -> tuple[float, float, float]:
        if not np.isfinite(allowed_molecular_uptake) or allowed_molecular_uptake < 0:
            raise ValueError("allowed_molecular_uptake must be finite and non-negative")
        allowed = float(allowed_molecular_uptake)
        if allowed > self.uptake_grid[-1]:
            raise ValueError("allowed_molecular_uptake exceeds the configured uptake_vmax")
        right = int(np.searchsorted(self.uptake_grid, allowed))
        nearest = min((right, max(0, right - 1)),
                      key=lambda index: abs(float(self.uptake_grid[index]) - allowed))
        node = float(self.uptake_grid[nearest])
        if abs(node - allowed) <= 4.0 * np.spacing(max(node, allowed)):
            if self.feasible[nearest]:
                return tuple(float(value) for value in self.grid_fluxes[nearest])
            return self._infeasible_rates(allowed)
        left = right - 1
        if self.feasible[left] and self.feasible[right]:
            fraction = (allowed - self.uptake_grid[left]) / (
                self.uptake_grid[right] - self.uptake_grid[left])
            return tuple(float(value) for value in (
                (1.0 - fraction) * self.grid_fluxes[left] + fraction * self.grid_fluxes[right]))
        if allowed < self.minimum_feasible_uptake_mmol_per_gdcw_h:
            return self._infeasible_rates(allowed)
        if self._solve_at_uptake is not None:
            if allowed not in self._exact_cache:
                self._exact_cache[allowed] = self._solve_at_uptake(allowed)
            result = self._exact_cache[allowed]
            if result is not None:
                return result
        return self._infeasible_rates(allowed)

    def __call__(self, substrate_mmol_per_l: float) -> tuple[float, float, float]:
        if not np.isfinite(substrate_mmol_per_l) or substrate_mmol_per_l < 0:
            raise ValueError("substrate_mmol_per_l must be finite and non-negative")
        substrate = float(substrate_mmol_per_l)
        if substrate == 0:
            return self.at_uptake(0.0)
        if substrate >= self.substrate_km_mmol_per_l:
            saturation = 1.0 / (1.0 + self.substrate_km_mmol_per_l / substrate)
        else:
            ratio = substrate / self.substrate_km_mmol_per_l
            saturation = ratio / (1.0 + ratio)
        return self.at_uptake(float(self.uptake_grid[-1]) * saturation)


def _close_other_carbon(model, retained: set[str]) -> None:
    for reaction in model.boundary:
        if reaction.id in retained:
            continue
        for metabolite, coefficient in reaction.metabolites.items():
            if not metabolite.elements.get("C", 0):
                continue
            if coefficient < 0:
                reaction.lower_bound = max(0.0, reaction.lower_bound)
            else:
                reaction.upper_bound = min(0.0, reaction.upper_bound)


def _validate_network_options(uptake_vmax, uptake_km, oxygen_lower_bound, growth_fraction,
                              growth_retention, maximum_growth, allocation_fraction,
                              fixed_growth, grid_points) -> int:
    for name, value in (("uptake_vmax_mmol_per_gdcw_h", uptake_vmax),
                        ("uptake_km_mmol_per_l", uptake_km)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    for name, value in (("growth_fraction", growth_fraction),
                        ("growth_retention", growth_retention),
                        ("allocation_fraction", allocation_fraction)):
        if not np.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be finite and in [0, 1]")
    if not np.isfinite(oxygen_lower_bound) or oxygen_lower_bound > 0:
        raise ValueError("oxygen_lower_bound must be finite and non-positive")
    if maximum_growth is not None and (
            not np.isfinite(maximum_growth) or maximum_growth <= 0):
        raise ValueError("maximum_growth_per_h must be positive and finite")
    if fixed_growth is not None and (not np.isfinite(fixed_growth) or fixed_growth < 0):
        raise ValueError("fixed_growth_per_h must be finite and non-negative")
    if (isinstance(grid_points, (bool, np.bool_)) or not np.isfinite(grid_points)
            or int(grid_points) != grid_points or grid_points < 3):
        raise ValueError("grid_points must be a finite integer of at least 3")
    return int(grid_points)


def _solve_with_physical_precision(model):
    """Use a tighter working tolerance, not a looser physical acceptance criterion.

    GLPK checks feasibility in internally scaled coordinates; passing that test
    does not establish original-unit balance in an enzyme-constrained network.
    Use the declared working tolerance while solving, then validate the unscaled
    fluxes and bounds against PINNED_TOLERANCE below. Keep the warm basis across
    allocation stages, with a bounded primal attempt and dual recovery for a
    nonterminal status. No LP coefficients or bounds are modified, and the
    caller's configuration is restored even when solving fails.
    """
    configuration = model.solver.configuration
    feasibility = configuration.tolerances.feasibility
    presolve = configuration.presolve
    simplex_parameters = configuration._smcp
    method, iteration_limit = simplex_parameters.meth, simplex_parameters.it_lim
    try:
        configuration.tolerances.feasibility = _GLPK_WORKING_FEASIBILITY
        configuration.presolve = False
        simplex_parameters.it_lim = min(iteration_limit, _GLPK_SIMPLEX_ITERATIONS)
        simplex_parameters.meth = GLP_PRIMAL
        value = growth_or_none(model)
        if model.solver.status not in ("optimal", "infeasible", "unbounded"):
            simplex_parameters.meth = GLP_DUAL
            value = growth_or_none(model)
        return value
    finally:
        configuration.tolerances.feasibility = feasibility
        configuration.presolve = presolve
        simplex_parameters.meth, simplex_parameters.it_lim = method, iteration_limit


def _validated_optimize(model, stoichiometry, stage: str, *, allow_infeasible=False):
    value = _solve_with_physical_precision(model)
    status = str(model.solver.status)
    if status == "infeasible" and allow_infeasible:
        return None
    if status != "optimal" or value is None or not np.isfinite(value):
        raise RuntimeError(f"{stage}: metabolic solve failed ({status}, objective={value})")
    primal = model.solver.primal_values
    fluxes = np.array([
        primal[reaction.forward_variable.name] - primal[reaction.reverse_variable.name]
        for reaction in model.reactions
    ], dtype=float)
    if not np.all(np.isfinite(fluxes)):
        raise RuntimeError(f"{stage}: solver returned nonfinite reaction fluxes")
    bounds = np.array([reaction.bounds for reaction in model.reactions], dtype=float)
    violations = np.maximum(bounds[:, 0] - fluxes, fluxes - bounds[:, 1])
    if np.any(violations > PINNED_TOLERANCE):
        index = int(np.argmax(violations))
        raise RuntimeError(
            f"{stage}: flux bound violation for {model.reactions[index].id}: "
            f"{fluxes[index]} outside {tuple(bounds[index])}")
    residual = np.asarray(stoichiometry @ fluxes).ravel()
    if not np.all(np.isfinite(residual)) or np.any(np.abs(residual) > PINNED_TOLERANCE):
        index = int(np.argmax(np.abs(residual)))
        raise RuntimeError(
            f"{stage}: full flux solution violates S*v=0 at {model.metabolites[index].id}: "
            f"residual={residual[index]:.12g} (physical tolerance {PINNED_TOLERANCE:g})")
    for variable in model.solver.variables:
        observed = primal[variable.name]
        if (not np.isfinite(observed)
                or (variable.lb is not None and observed < variable.lb - PINNED_TOLERANCE)
                or (variable.ub is not None and observed > variable.ub + PINNED_TOLERANCE)):
            raise RuntimeError(f"{stage}: solver variable bounds violated for {variable.name}")
    for constraint in model.solver.constraints:
        observed = constraint.primal
        if constraint.name not in model.metabolites:
            observed = 0.0
            for term, coefficient in constraint.expression.as_coefficients_dict().items():
                observed += float(coefficient) * (1.0 if term == 1 else primal[term.name])
        if (observed is None or not np.isfinite(observed)
                or (constraint.lb is not None and observed < constraint.lb - PINNED_TOLERANCE)
                or (constraint.ub is not None and observed > constraint.ub + PINNED_TOLERANCE)):
            raise RuntimeError(f"{stage}: solver constraint violated for {constraint.name}")
    return float(value), fluxes


def _nonnegative_flux(value: float, name: str) -> float:
    if not np.isfinite(value) or value < -PINNED_TOLERANCE:
        raise RuntimeError(f"{name} must be a finite non-negative solved flux, got {value}")
    return max(0.0, float(value))


def _build_coupled_rates(model, substrate_exchange_id, oxygen_exchange_id, ngam_reaction_id,
                         *, product_reaction_id, biomass_reaction_id, uptake_vmax,
                         uptake_km, oxygen_lower_bound, growth_fraction, growth_retention,
                         maximum_growth, allocation_fraction, fixed_growth, grid_points,
                         requested_ngam=None, product_capacity=None, fixed_fraction_growth=False):
    if not product_reaction_id:
        raise ValueError("product_reaction_id is required for coupled metabolic production")
    scoped = model.copy()
    configure(scoped)
    substrate = _uptake_boundary(scoped, substrate_exchange_id)
    coefficient = float(next(iter(substrate.metabolites.values())))
    uptake_id = cap_uptake(scoped, substrate_exchange_id, uptake_vmax)
    oxygen_id = cap_uptake(scoped, oxygen_exchange_id, abs(oxygen_lower_bound))
    biomass = scoped.reactions.get_by_id(biomass_reaction_id)
    product = scoped.reactions.get_by_id(product_reaction_id)
    maintenance = scoped.reactions.get_by_id(ngam_reaction_id)
    roles = (uptake_id, oxygen_id, biomass_reaction_id, product_reaction_id, ngam_reaction_id)
    if len(set(roles)) != len(roles) or product_reaction_id == substrate_exchange_id:
        raise ValueError("product, substrate, oxygen, biomass and maintenance must be different reactions")
    if len(product.metabolites) != 1:
        raise ValueError("product_reaction_id must name the actual one-metabolite product demand")
    product_coefficient = -float(next(iter(product.metabolites.values())))
    if not np.isfinite(product_coefficient) or product_coefficient <= 0:
        raise ValueError("the declared product demand cannot carry forward production")
    _close_other_carbon(scoped, {substrate_exchange_id, uptake_id})
    ngam = max(0.0, maintenance.lower_bound,
               requested_ngam if requested_ngam is not None else 0.0)
    if not np.isfinite(ngam) or ngam > maintenance.upper_bound:
        raise ValueError(f"maintenance {ngam_reaction_id} cannot meet NGAM {ngam} within its capacity")
    maintenance.lower_bound = ngam
    product_lower = max(0.0, product.lower_bound)
    product_upper = product.upper_bound
    if product_capacity is not None:
        product_upper = min(product_upper, product_capacity / product_coefficient)
    if product_lower > product_upper:
        raise ValueError("product demand bounds are infeasible under its kinetic capacity")
    product.bounds = (product_lower, product_upper)
    growth_lower = max(0.0, biomass.lower_bound)
    if growth_lower > biomass.upper_bound:
        raise ValueError("biomass reaction cannot carry non-negative growth")
    biomass.bounds = (growth_lower, biomass.upper_bound)
    stoichiometry = create_stoichiometric_matrix(scoped, array_type="lil").tocsr()
    uptake_expression = coefficient * substrate.flux_expression
    product_expression = product_coefficient * product.flux_expression
    scoped.objective = biomass
    attainable, _ = _validated_optimize(scoped, stoichiometry, "reference attainable growth")
    attainable = _nonnegative_flux(attainable, "attainable growth")
    limit = min(attainable, maximum_growth) if maximum_growth is not None else attainable
    limit *= growth_retention
    if limit < growth_lower:
        raise ValueError("growth_retention or maximum_growth_per_h conflicts with required growth")
    if fixed_growth is not None and not growth_lower <= fixed_growth <= limit:
        raise ValueError(
            f"fixed_growth_per_h={fixed_growth:g} is infeasible within growth bounds "
            f"[{growth_lower:g}, {limit:g}]")
    growth_bounds = ((fixed_growth, fixed_growth) if fixed_growth is not None
                     else (growth_lower, limit))
    biomass.bounds = growth_bounds
    uptake_index = scoped.reactions.index(substrate)
    growth_index = scoped.reactions.index(biomass)
    product_index = scoped.reactions.index(product)

    with scoped:
        name = "__ystwin_allowed_uptake"
        while name in scoped.solver.variables:
            name += "_"
        allowed_variable = scoped.problem.Variable(name, lb=0.0, ub=uptake_vmax)
        constraints = [scoped.problem.Constraint(
            uptake_expression - allowed_variable, ub=0.0, name=f"{name}_carbon")]
        if product_capacity is not None:
            capacity_per_uptake = product_capacity / uptake_vmax
            if not np.isfinite(capacity_per_uptake):
                raise ValueError("product capacity per molecular uptake must be finite")
            constraints.append(scoped.problem.Constraint(
                product_expression - capacity_per_uptake * allowed_variable,
                ub=0.0, name=f"{name}_product"))
        scoped.add_cons_vars([allowed_variable, *constraints])
        scoped.objective = scoped.problem.Objective(allowed_variable, direction="min")
        minimum, _ = _validated_optimize(scoped, stoichiometry, "minimum feasible uptake")
        minimum = _nonnegative_flux(minimum, "minimum feasible uptake")

    def solve_at_uptake(allowed):
        if allowed < minimum:
            return None
        with scoped:
            cap_uptake(scoped, substrate_exchange_id, allowed)
            biomass.bounds = growth_bounds
            upper = product_upper
            if product_capacity is not None:
                upper = min(upper, product_capacity * (allowed / uptake_vmax) / product_coefficient)
            if upper < product_lower:
                return None
            product.bounds = (product_lower, upper)
            scoped.objective = biomass
            maximum = _validated_optimize(scoped, stoichiometry, "attainable growth at uptake",
                                          allow_infeasible=True)
            if maximum is None:
                return None
            maximum_growth_at_uptake = _nonnegative_flux(maximum[0], "attainable growth")
            if fixed_growth is None:
                floor = max(growth_lower, maximum_growth_at_uptake * growth_fraction)
                biomass.bounds = ((floor, floor) if fixed_fraction_growth else (floor, limit))
            scoped.objective = scoped.problem.Objective(product_expression, direction="min")
            minimum_product, _ = _validated_optimize(scoped, stoichiometry, "minimum product")
            scoped.objective = scoped.problem.Objective(product_expression, direction="max")
            maximum_product, _ = _validated_optimize(scoped, stoichiometry, "maximum product")
            minimum_product = _nonnegative_flux(minimum_product, "minimum product")
            maximum_product = _nonnegative_flux(maximum_product, "maximum product")
            if maximum_product < minimum_product - PINNED_TOLERANCE:
                raise RuntimeError("solved product interval has q_max below q_min")
            target = minimum_product + allocation_fraction * (maximum_product - minimum_product)
            product.bounds = (target / product_coefficient, target / product_coefficient)
            scoped.objective = biomass
            growth, _ = _validated_optimize(scoped, stoichiometry, "growth at allocated product")
            growth = _nonnegative_flux(growth, "allocated growth")
            biomass.bounds = (growth, growth)
            scoped.objective = scoped.problem.Objective(uptake_expression, direction="min")
            _, fluxes = _validated_optimize(scoped, stoichiometry, "minimum actual carbon uptake")
            return (
                _nonnegative_flux(fluxes[growth_index], "growth"),
                _nonnegative_flux(coefficient * fluxes[uptake_index], "actual carbon uptake"),
                _nonnegative_flux(product_coefficient * fluxes[product_index], "product"),
            )

    grid = np.unique(np.append(np.linspace(0.0, uptake_vmax, grid_points), minimum))
    fluxes = np.zeros((len(grid), 3))
    feasible = np.zeros(len(grid), dtype=bool)
    statuses = []
    for index, allowed in enumerate(grid):
        observed = solve_at_uptake(float(allowed))
        if observed is None:
            statuses.append("infeasible")
        else:
            fluxes[index] = observed
            feasible[index] = True
            statuses.append("optimal")
    if not np.any(feasible):
        raise ValueError("no feasible coupled metabolic state exists within the uptake range")
    rates = CoupledRates(grid, fluxes, feasible, uptake_km, product_reaction_id,
                         tuple(statuses), minimum, solve_at_uptake, fixed_growth)
    return rates, ngam, attainable, limit, minimum


def build_rates(model,
                glucose_exchange_id: str,
                oxygen_exchange_id: str,
                ngam_reaction_id: str,
                kcat_per_s: float,
                enzyme_mmol_per_gdcw: float,
                *,
                stress: float = 0.0,
                regulation_scale: float = 1.0,
                uptake_vmax_mmol_per_gdcw_h: float = 10.0,
                uptake_km_mmol_per_l: float = 0.5,
                oxygen_lower_bound: float = -1000.0,
                grid_points: int = 21,
                product_reaction_id: str | None = None,
                biomass_reaction_id: str = "r_2111",
                growth_fraction: float = 0.9,
                growth_retention: float = 1.0,
                maximum_growth_per_h: float | None = None):
    """Build the ``(mu, q_S, q_P)`` closure `simulate_batch` needs.

    The GEM is solved once per grid point rather than once per integration step, and the
    growth rate is interpolated between them. That is not an approximation of convenience:
    an integration with a few hundred steps across several conditions is thousands of LP
    solves on a four-thousand-reaction model, which does not finish, and the map from one
    uptake bound to one growth rate is smooth and one-dimensional. `fba/surrogate.py` exists
    for exactly this reason and says so in its own docstring.

    Args:
        model: A cobra model with the pathway already installed.
        glucose_exchange_id, oxygen_exchange_id, ngam_reaction_id: Reaction ids, passed
            rather than looked up, because a name lookup by substring is how this project
            once matched 2-phenylethanol for ethanol.
        kcat_per_s: Turnover of the enzyme that sets the product rate.
        enzyme_mmol_per_gdcw: How much of it there is. From `pathway/proteome.py` for a
            native enzyme; for a heterologous one there is no measurement and a promoter
            anchor is a CEILING, which the caller must label.
        stress: 0-1. Constrains NGAM before every solve, which is how the latent stress
            layer reaches growth.
        regulation_scale: Multiplier on the enzyme abundance. This is the route by which
            regulation acts -- on how much enzyme exists, not on a reaction ceiling.

    Returns:
        ``(closure, RateClosure)`` -- the function, and what built it.
    """
    if (not np.isfinite(kcat_per_s) or not np.isfinite(enzyme_mmol_per_gdcw)
            or kcat_per_s <= 0 or enzyme_mmol_per_gdcw <= 0):
        raise ValueError(
            f"kcat_per_s and enzyme_mmol_per_gdcw must both be positive, got "
            f"{kcat_per_s} and {enzyme_mmol_per_gdcw}. There is no default for either: a "
            f"turnover comes from a measurement or a sequence-based predictor, and an "
            f"abundance from proteomics. Defaulting them would make the titre a property "
            f"of this module rather than of the strain")
    if not 0.0 <= stress <= 1.0:
        raise ValueError(f"stress is a fraction of full induction, 0-1, got {stress}")
    if not product_reaction_id:
        raise ValueError(
            "product_reaction_id is required: product must consume precursors inside "
            "the same metabolic solve as growth, not be added after optimization")
    if not np.isfinite(regulation_scale) or regulation_scale <= 0:
        raise ValueError("regulation_scale must be positive and finite")
    if not 0 < growth_fraction <= 1:
        raise ValueError("growth_fraction must be in (0, 1]")
    count = _validate_network_options(
        uptake_vmax_mmol_per_gdcw_h, uptake_km_mmol_per_l, oxygen_lower_bound,
        growth_fraction, growth_retention, maximum_growth_per_h, 1.0, None, grid_points)
    q_p_max = float(kcat_per_s) * 3600.0 * float(enzyme_mmol_per_gdcw) * float(regulation_scale)
    if not np.isfinite(q_p_max) or q_p_max <= 0:
        raise ValueError("the combined kcat/enzyme/regulation product capacity must be positive and finite")
    rates, ngam, _, _, _ = _build_coupled_rates(
        model, glucose_exchange_id, oxygen_exchange_id, ngam_reaction_id,
        product_reaction_id=product_reaction_id, biomass_reaction_id=biomass_reaction_id,
        uptake_vmax=uptake_vmax_mmol_per_gdcw_h, uptake_km=uptake_km_mmol_per_l,
        oxygen_lower_bound=oxygen_lower_bound, growth_fraction=growth_fraction,
        growth_retention=growth_retention, maximum_growth=maximum_growth_per_h,
        allocation_fraction=1.0, fixed_growth=None, grid_points=count,
        requested_ngam=RESTING_NGAM + stress * MAX_STRESS_NGAM,
        product_capacity=q_p_max, fixed_fraction_growth=True)

    # Product formation saturates with the same carbon availability that limits growth.
    # A product rate held constant while the culture starves is the step function that
    # made an earlier version of this simulation depend on its own step count.
    return rates, RateClosure(
        kcat_per_s, enzyme_mmol_per_gdcw, regulation_scale,
        stress, ngam, uptake_vmax_mmol_per_gdcw_h, uptake_km_mmol_per_l,
        product_reaction_id, biomass_reaction_id, growth_fraction, growth_retention,
        maximum_growth_per_h)


def build_network_rates(model, substrate_exchange_id, oxygen_exchange_id, ngam_reaction_id,
                        *, product_reaction_id, biomass_reaction_id="r_2111",
                        uptake_vmax_mmol_per_gdcw_h=10.0, uptake_km_mmol_per_l=0.5,
                        oxygen_lower_bound=-1000.0, growth_fraction=0.9,
                        growth_retention=1.0, maximum_growth_per_h=None,
                        allocation_fraction=0.25, fixed_growth_per_h=None,
                        grid_points=13) -> tuple[CoupledRates, NetworkRateClosure]:
    count = _validate_network_options(
        uptake_vmax_mmol_per_gdcw_h, uptake_km_mmol_per_l, oxygen_lower_bound,
        growth_fraction, growth_retention, maximum_growth_per_h,
        allocation_fraction, fixed_growth_per_h, grid_points)
    rates, ngam, attainable, limit, minimum = _build_coupled_rates(
        model, substrate_exchange_id, oxygen_exchange_id, ngam_reaction_id,
        product_reaction_id=product_reaction_id, biomass_reaction_id=biomass_reaction_id,
        uptake_vmax=uptake_vmax_mmol_per_gdcw_h, uptake_km=uptake_km_mmol_per_l,
        oxygen_lower_bound=oxygen_lower_bound, growth_fraction=growth_fraction,
        growth_retention=growth_retention, maximum_growth=maximum_growth_per_h,
        allocation_fraction=allocation_fraction, fixed_growth=fixed_growth_per_h,
        grid_points=count)
    return rates, NetworkRateClosure(
        substrate_exchange_id, oxygen_exchange_id, ngam_reaction_id,
        product_reaction_id, biomass_reaction_id, ngam, uptake_vmax_mmol_per_gdcw_h,
        uptake_km_mmol_per_l, oxygen_lower_bound, growth_fraction, growth_retention,
        maximum_growth_per_h, allocation_fraction, fixed_growth_per_h,
        attainable, limit, minimum)
