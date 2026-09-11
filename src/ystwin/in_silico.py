"""Explicit synthetic controls driving the network/fixed-volume comparison.

``simulate_controlled_comparison`` -- retained under its older name
``simulate_controlled_product`` -- is a comparison route, not the shared physical engine.
It does not map reporter stress to biological enzyme abundances. A caller supplies the
existing four-control schedule; same-teacher recovery tests only its approximation within
that synthetic world, and the trajectory's own summary reports ``biological_validation``
as ``False``. It integrates through ``fba.dynamic.simulate_fixed_volume_comparison``,
named in full below so that no call site here reads as if it were the shared engine.

The shared engine (:mod:`ystwin.mech.engine`, reached through
``mech.adapters.simulate_protocol``) is deliberately NOT re-exported here, unlike in
``generator/culture.py``, ``fba/dynamic.py`` and ``hybrid.py``. It consumes ``Protocol``
and ``PhysicalState`` and refuses to convert anything else, so a four-control GSMM
schedule has no route into it and a re-export would only imply one.
``docs/MECHANISTIC_LAYER.md`` section 0 carries the repository-wide entry map.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import cobra
import numpy as np

from .fba.dynamic import simulate_fixed_volume_comparison
from .fba.dynamic_rates import build_network_rates
from .fba.physiology import _uptake_boundary, cap_uptake
from .fba.product_panel import ProductTask

__all__ = [
    "InSilicoEnvironment", "MetabolicControl", "ControlledTrajectory",
    "apply_synthetic_control", "simulate_controlled_comparison",
    "simulate_controlled_product",
]

_CARBON = {"glucose": ("r_1714", "C6H12O6"), "ethanol": ("r_1761", "C2H6O")}


@dataclass(frozen=True)
class MetabolicControl:
    """The teacher/student control order; NGAM is per reference gDCW, others are fractions."""

    growth_retention: float
    enzyme_budget_scale: float
    ngam_mmol_per_gdcw_h: float
    allocation_fraction: float

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if not 0 <= self.growth_retention <= 1:
            raise ValueError("growth_retention must be in [0, 1]")
        if not 0 < self.enzyme_budget_scale <= 1:
            raise ValueError("enzyme_budget_scale must be in (0, 1]")
        if self.ngam_mmol_per_gdcw_h < 0:
            raise ValueError("ngam_mmol_per_gdcw_h must be nonnegative")
        if not 0 <= self.allocation_fraction <= 1:
            raise ValueError("allocation_fraction must be in [0, 1]")


@dataclass(frozen=True)
class InSilicoEnvironment:
    """Carbon uptake is mmol C/reference gDCW/h; oxygen is mmol/reference gDCW/h.

    The inoculum is structural biomass. The caller must supply a network on the
    same structural basis; ``structural_mass_fraction`` converts new reference
    uptake/maintenance bounds, not the network or the inoculum itself.
    """

    name: str
    carbon_source: str = "glucose"
    substrate_g_per_l: float = 20.0
    carbon_uptake: float = 20.0
    oxygen_uptake: float = 1000.0
    initial_structural_biomass_g_per_l: float = 0.05
    hours: float = 8.0

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("environment name is required")
        if self.carbon_source not in _CARBON:
            raise ValueError("carbon_source must be glucose or ethanol")
        for name in ("substrate_g_per_l", "carbon_uptake", "initial_structural_biomass_g_per_l", "hours"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not np.isfinite(self.oxygen_uptake) or self.oxygen_uptake < 0:
            raise ValueError("oxygen_uptake must be finite and nonnegative")


def apply_synthetic_control(
    model: cobra.Model, control: MetabolicControl, *, structural_mass_fraction: float = 1.0,
) -> tuple[cobra.Model, dict]:
    if not np.isfinite(structural_mass_fraction) or not 0 < structural_mass_fraction <= 1:
        raise ValueError("structural_mass_fraction must be in (0, 1]")
    out = model.copy()
    maintenance = out.reactions.get_by_id("r_4046")
    previous_maintenance = maintenance.bounds
    demand = control.ngam_mmol_per_gdcw_h / structural_mass_fraction
    maintenance.bounds = (demand, demand)
    pool = _uptake_boundary(out, "prot_pool_exchange")
    coefficient = float(next(iter(pool.metabolites.values())))
    old_supply = coefficient * (pool.upper_bound if coefficient > 0 else pool.lower_bound)
    if not np.isfinite(old_supply) or old_supply <= 0:
        raise ValueError("the reference protein supply must have a positive finite capacity")
    new_supply = old_supply * control.enzyme_budget_scale
    cap_uptake(out, "prot_pool_exchange", new_supply)
    return out, {
        "control_source": "declared or learned in-silico controller",
        "control": asdict(control),
        "structural_mass_fraction": structural_mass_fraction,
        "ngam_reference_units": "mmol ATP per reference gDCW per hour",
        "ngam_applied_units": "mmol ATP per structural g per hour",
        "ngam_before": list(previous_maintenance), "ngam_applied": demand,
        "protein_supply_before": old_supply, "protein_supply_applied": new_supply,
        "intervention": "explicit simulated maintenance setpoint and protein budget; not inferred wet-lab calibration",
    }


@dataclass(frozen=True)
class ControlledTrajectory:
    times_h: np.ndarray
    structural_biomass_g_per_l: np.ndarray
    substrate_mmol_per_l: np.ndarray
    product_mmol_per_l: np.ndarray
    total_dry_biomass_g_per_l: np.ndarray
    growth_rate_per_h: np.ndarray
    initial_substrate_mmol_per_l: float
    total_substrate_consumed_mmol_per_l: float
    structural_biomass_exposure_g_h_per_l: float
    product: ProductTask
    substrate_molar_mass_g_per_mol: float
    intracellular: bool
    stage_provenance: tuple[dict, ...]

    def summary(self) -> dict:
        mass = self.product.molar_mass_g_per_mol
        net_product = float(self.product_mmol_per_l[-1] - self.product_mmol_per_l[0])
        titre = float(self.product_mmol_per_l[-1] * mass)
        duration = float(self.times_h[-1] - self.times_h[0])
        consumed_mass = self.total_substrate_consumed_mmol_per_l * self.substrate_molar_mass_g_per_mol
        total_exposure = self.structural_biomass_exposure_g_h_per_l
        if self.intracellular:
            total_exposure += float(np.trapezoid(self.product_mmol_per_l, self.times_h)) * mass / 1000
        return {
            "product": self.product.name,
            "simulation_route": "in_silico_gsmm_control_comparison",
            "biological_validation": False,
            "prediction_basis": self.product.output_basis,
            "product_molar_mass_g_per_mol": mass,
            "substrate_molar_mass_g_per_mol": self.substrate_molar_mass_g_per_mol,
            "rate_basis": "net product accumulation divided by biomass-time exposure",
            "yield_basis": "net product mass / consumed substrate mass; initial product inventory excluded",
            "total_dry_mass_basis": "structural biomass plus tracked intracellular product only; extracellular product excluded",
            "total_dry_mass_exposure_basis": "integrated structural biomass plus trapezoidal tracked intracellular product mass",
            "growth_rate_sampling": "preceding interval at control knots; sampled reporter driver, not a resolved discontinuity",
            "titre_mg_per_l": titre,
            "productivity_mg_per_l_h": net_product * mass / duration,
            "yield_g_per_g": net_product * mass / consumed_mass if consumed_mass > 0 else float("nan"),
            "mean_specific_rate_mmol_per_g_struct_h": net_product / self.structural_biomass_exposure_g_h_per_l,
            "mean_specific_rate_mmol_per_g_total_dw_h": net_product / total_exposure,
            "content_mg_per_g_total_dw": titre / self.total_dry_biomass_g_per_l[-1],
            "final_structural_biomass_g_per_l": float(self.structural_biomass_g_per_l[-1]),
            "final_total_dry_biomass_g_per_l": float(self.total_dry_biomass_g_per_l[-1]),
            "total_substrate_consumed_mmol_per_l": self.total_substrate_consumed_mmol_per_l,
            "pool_location": "intracellular" if self.intracellular else "extracellular",
            "truth_scope": "explicit in-silico controller and network, not measured biological production",
            "storage_assumption": "tracked pool starts at its declared inventory; no remobilization in this synthesis-only protocol",
        }


def simulate_controlled_comparison(
    model: cobra.Model,
    product: ProductTask,
    control_times_h,
    controls,
    environment: InSilicoEnvironment,
    *,
    structural_mass_fraction: float = 1.0,
    initial_pool_mmol_per_l: float = 0.0,
    growth_fraction: float = 0.9,
    grid_points: int = 5,
    output_dt: float = 0.1,
) -> ControlledTrajectory:
    times = np.asarray(control_times_h, dtype=float)
    values = np.asarray(controls, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all()
            or times[0] != 0 or np.any(np.diff(times) <= 0)
            or not np.isclose(times[-1], environment.hours, rtol=0, atol=1e-10)
            or values.shape != (len(times) - 1, 4)):
        raise ValueError("control schedule requires ascending [0, hours] times and one four-control row per interval")
    control_objects = tuple(MetabolicControl(*row) for row in values)
    if not np.isfinite(output_dt) or output_dt <= 0:
        raise ValueError("output_dt must be positive and finite")
    if not np.isfinite(initial_pool_mmol_per_l) or initial_pool_mmol_per_l < 0:
        raise ValueError("initial_pool_mmol_per_l must be finite and nonnegative")
    if not np.isfinite(structural_mass_fraction) or not 0 < structural_mass_fraction <= 1:
        raise ValueError("structural_mass_fraction must be in (0, 1]")
    exchange_id, formula = _CARBON[environment.carbon_source]
    exchange = model.reactions.get_by_id(exchange_id)
    if len(exchange.metabolites) != 1:
        raise ValueError("the carbon exchange must identify one substrate")
    substrate = next(iter(exchange.metabolites))
    if substrate.formula != formula:
        raise ValueError("carbon substrate identity disagrees with the declared environment")
    substrate_mass = float(substrate.formula_weight)
    carbon_atoms = float(substrate.elements["C"])
    compartment = model.metabolites.get_by_id(product.metabolite_id).compartment
    if compartment not in ("c", "e"):
        raise ValueError("product pool location must be explicitly cytosolic or extracellular")
    intracellular = compartment == "c"
    initial_substrate = environment.substrate_g_per_l * 1000 / substrate_mass
    biomass, remaining, pool = environment.initial_structural_biomass_g_per_l, initial_substrate, initial_pool_mmol_per_l
    parts = {name: [] for name in ("time", "biomass", "substrate", "product", "growth")}
    consumed = exposure = 0.0
    provenance = []
    for index, (start, end, control) in enumerate(zip(times[:-1], times[1:], control_objects)):
        controlled, record = apply_synthetic_control(
            model, control, structural_mass_fraction=structural_mass_fraction)
        rates, info = build_network_rates(
            controlled, exchange_id, "r_1992", "r_4046", product_reaction_id=product.reaction_id,
            uptake_vmax_mmol_per_gdcw_h=environment.carbon_uptake / carbon_atoms / structural_mass_fraction,
            oxygen_lower_bound=-environment.oxygen_uptake / structural_mass_fraction,
            growth_retention=control.growth_retention, allocation_fraction=control.allocation_fraction,
            growth_fraction=growth_fraction, grid_points=grid_points)
        result = simulate_fixed_volume_comparison(
            rates, biomass, remaining, float(end - start), product.molar_mass_g_per_mol,
            steps=max(2, int(np.ceil((end - start) / output_dt))),
            substrate_molar_mass_g_per_mol=substrate_mass)
        keep = slice(None) if index == 0 else slice(1, None)
        parts["time"].append(start + result.hours[keep])
        parts["biomass"].append(result.biomass_g_per_l[keep])
        parts["substrate"].append(result.substrate_mmol_per_l[keep])
        parts["product"].append(pool + result.product_mmol_per_l[keep])
        # Keep the preceding interval's endpoint growth sample. The reporter API
        # linearly interpolates this driver; a right-hand value at the knot would
        # leak the next control backwards into the observation prefix. This is a
        # sampled comparison, not an exactly resolved discontinuous growth driver.
        parts["growth"].append(result.growth_rate_per_h[keep])
        biomass = float(result.biomass_g_per_l[-1])
        remaining = float(result.substrate_mmol_per_l[-1])
        pool += float(result.product_mmol_per_l[-1])
        consumed += result.substrate_consumed_mmol_per_l
        exposure += result.biomass_exposure_g_h_per_l
        provenance.append({"start_h": float(start), "end_h": float(end),
                           **record, "rate_model": asdict(info),
                           "integrator_route": result.simulation_route,
                           "specific_rate_units": "mmol per structural g per hour; growth per hour"})
    joined = {name: np.concatenate(arrays) for name, arrays in parts.items()}
    total_mass = joined["biomass"].copy()
    if intracellular:
        total_mass += joined["product"] * product.molar_mass_g_per_mol / 1000
    return ControlledTrajectory(
        joined["time"], joined["biomass"], joined["substrate"], joined["product"], total_mass,
        joined["growth"], initial_substrate, float(consumed), float(exposure), product,
        substrate_mass, intracellular, tuple(provenance))


# The retained legacy name for the comparison route above. It is not, and must not become,
# an alias for the shared engine: the two run different physics on different contracts.
simulate_controlled_product = simulate_controlled_comparison
