"""Two entry points, two different physics. Neither of them runs silently.

``simulate_protocol`` is the shared physical engine (:mod:`ystwin.mech.engine`), re-exported
here unchanged from :mod:`ystwin.mech.adapters`. It takes the shared ``Protocol`` and
``PhysicalState`` contracts, converts nothing, and is scored against no dataset in this
repository.

``simulate_hybrid_comparison`` -- retained under its older name ``simulate_hybrid_case`` --
is the allocation-conditioned comparison route. Only the existing empirical
growth-retention factor can couple the sensor estimate to the comparison. It supplies
neither gene-specific enzyme regulation nor a learned ATP-maintenance law, its own summary
reports ``biological_validation`` as ``False``, and it is what produced every committed
``outputs/hybrid_benchmark/`` table. It integrates through
``fba.dynamic.simulate_fixed_volume_comparison``, which is named in full below so that no
call site here reads as if it were the shared engine.

``docs/MECHANISTIC_LAYER.md`` section 0 carries the repository-wide entry map.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass

import cobra
import numpy as np

from .analysis.hybrid_stress import CalibratedSensorModel, StressEstimate
from .fba.dynamic import BatchResult, simulate_fixed_volume_comparison
from .fba.dynamic_rates import NetworkRateClosure, build_network_rates
from .fba.product_panel import ProductTask
from .generator.context import CultureContext, context_growth_rate
from .mech.adapters import simulate_protocol

__all__ = ["Environment", "HybridCase", "ScenarioUnsupported", "simulate_hybrid_case",
           "simulate_hybrid_comparison", "simulate_protocol"]

_CARBON_EXCHANGES = {"glucose": "r_1714", "ethanol": "r_1761"}
_CARBON_FORMULAS = {"glucose": "C6H12O6", "ethanol": "C2H6O"}


class ScenarioUnsupported(ValueError):
    pass


@dataclass(frozen=True)
class Environment:
    name: str
    carbon_source: str = "glucose"
    substrate_g_per_l: float = 20.0
    carbon_uptake: float = 20.0
    oxygen_uptake: float = 1000.0
    temperature_c: float = 30.0
    initial_biomass_g_per_l: float = 0.05
    hours: float = 24.0
    stressor: str | None = None
    dose_mM: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("environment name must be a nonempty string")
        if self.carbon_source not in _CARBON_EXCHANGES:
            raise ValueError(f"unsupported carbon_source {self.carbon_source!r}")
        for name in ("substrate_g_per_l", "carbon_uptake", "initial_biomass_g_per_l", "hours"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not np.isfinite(self.oxygen_uptake) or self.oxygen_uptake < 0:
            raise ValueError("oxygen_uptake must be finite and nonnegative")
        if not np.isfinite(self.temperature_c) or not 0 < self.temperature_c < 50:
            raise ValueError("temperature_c must be finite and between 0 and 50")
        if self.stressor not in (None, "DTT", "H2O2"):
            raise ValueError("stressor must be DTT, H2O2 or None")
        if (not np.isfinite(self.dose_mM) or self.dose_mM < 0
                or (self.stressor is None and self.dose_mM != 0)):
            raise ValueError("dose_mM must be finite, nonnegative and paired with a stressor")


@dataclass(frozen=True)
class HybridCase:
    product: ProductTask
    environment: Environment
    stress: StressEstimate
    trajectory: BatchResult
    rate_info: NetworkRateClosure
    substrate_carbon_atoms: float
    couple_stress: bool
    grid_points_requested: int
    grid_points_feasible: int
    grid_points_total: int

    def summary(self) -> dict:
        trajectory = self.trajectory
        consumed_carbon = (trajectory.substrate_consumed_mmol_per_l
                           * self.substrate_carbon_atoms)
        product_carbon = trajectory.product_mmol_per_l[-1] * self.product.carbon_atoms
        return {
            "product": self.product.name,
            "environment": self.environment.name,
            "status": "predicted",
            "prediction_kind": "allocation-conditioned scenario; not validated realized production",
            "simulation_route": "gsmm_allocation_comparison",
            "biological_validation": False,
            "prediction_basis": self.product.output_basis,
            "biomass_basis": "supplied network dry-mass basis; extra product is not added to biomass",
            "product_targets_used_for_fitting": False,
            "carbon_source": self.environment.carbon_source,
            "initial_substrate_g_per_l": self.environment.substrate_g_per_l,
            "carbon_uptake_bound_mmol_c_per_gdcw_h": self.environment.carbon_uptake,
            "oxygen_uptake_bound_mmol_per_gdcw_h": self.environment.oxygen_uptake,
            "temperature_c": self.environment.temperature_c,
            "initial_biomass_g_per_l": self.environment.initial_biomass_g_per_l,
            "hours": self.environment.hours,
            "stressor": self.environment.stressor or "none",
            "dose_mM": self.environment.dose_mM,
            "latent_upr": self.stress.upr,
            "latent_oxidative": self.stress.oxidative,
            "learned_growth_retention": self.stress.growth_retention,
            "applied_growth_retention": self.rate_info.growth_retention,
            "state_source": self.stress.source,
            "state_extrapolated": self.stress.extrapolated,
            "couple_stress": self.couple_stress,
            "allocation_fraction": self.rate_info.allocation_fraction,
            "growth_fraction": self.rate_info.growth_fraction,
            "titre_mg_per_l": trajectory.titre_mg_per_l,
            "productivity_mg_per_l_h": trajectory.average_productivity_mg_per_l_h,
            "mean_specific_rate_mmol_per_gdcw_h": trajectory.mean_specific_rate_mmol_per_gdcw_h,
            "yield_g_per_g": trajectory.yield_g_per_g,
            "carbon_yield": float(product_carbon / consumed_carbon)
            if consumed_carbon > 0 else float("nan"),
            "extra_content_mg_per_gdcw": trajectory.content_mg_per_gdcw,
            "final_biomass_g_per_l": trajectory.final_biomass_g_per_l,
            "substrate_consumed_mmol_per_l": trajectory.substrate_consumed_mmol_per_l,
            "grid_points_requested": self.grid_points_requested,
            "grid_points_feasible": self.grid_points_feasible,
            "grid_points_total": self.grid_points_total,
            "assumptions": "; ".join((*self.product.assumptions,
                "fixed-volume well-mixed batch; zero initial extra product pool",
                ("calibrated sensor growth retention transfers as an empirical growth-cap factor"
                 if self.couple_stress else "sensor state is reported only; growth retention is not applied"),
                "no learned ATP-maintenance or gene-specific protein-regulation coefficient",
                "temperature/carbon growth ceiling comes from the existing culture-context prior",
                "integrated product is net extra accumulation, not all preexisting cellular product")),
        }

    def provenance(self) -> dict:
        return {"environment": asdict(self.environment), "product": asdict(self.product),
                "stress": asdict(self.stress), "rate_model": asdict(self.rate_info)}


def simulate_hybrid_comparison(
    sensor_model: CalibratedSensorModel,
    metabolic_model: cobra.Model,
    product: ProductTask,
    environment: Environment,
    *,
    couple_stress: bool = True,
    observed_features: Mapping[str, float] | None = None,
    allocation_fraction: float = 0.25,
    growth_fraction: float = 0.9,
    grid_points: int = 9,
    steps: int = 200,
) -> HybridCase:
    stress = (sensor_model.predict_condition(environment.stressor, environment.dose_mM)
              if observed_features is None else sensor_model.infer(observed_features))
    if stress.extrapolated:
        raise ScenarioUnsupported("sensor condition is outside the calibrated support")
    exchange_id = _CARBON_EXCHANGES[environment.carbon_source]
    exchange = metabolic_model.reactions.get_by_id(exchange_id)
    if len(exchange.metabolites) != 1:
        raise ScenarioUnsupported("carbon exchange must identify exactly one substrate")
    substrate = next(iter(exchange.metabolites))
    if substrate.formula != _CARBON_FORMULAS[environment.carbon_source]:
        raise ScenarioUnsupported("carbon exchange formula disagrees with the declared carbon source")
    substrate_mass = float(substrate.formula_weight)
    carbon_atoms = float(substrate.elements["C"])
    context = CultureContext(carbon_source=environment.carbon_source,
                             temperature_c=environment.temperature_c)
    growth_limit = context_growth_rate(context)
    if growth_limit <= 0:
        raise ScenarioUnsupported("culture-context prior supports no growth in this environment")
    rates, info = build_network_rates(
        metabolic_model, exchange_id, "r_1992", "r_4046",
        product_reaction_id=product.reaction_id,
        uptake_vmax_mmol_per_gdcw_h=environment.carbon_uptake / carbon_atoms,
        oxygen_lower_bound=-environment.oxygen_uptake,
        growth_fraction=growth_fraction,
        growth_retention=stress.growth_retention if couple_stress else 1.0,
        maximum_growth_per_h=growth_limit,
        allocation_fraction=allocation_fraction,
        grid_points=grid_points)
    trajectory = simulate_fixed_volume_comparison(
        rates, environment.initial_biomass_g_per_l,
        environment.substrate_g_per_l * 1000.0 / substrate_mass,
        environment.hours, product.molar_mass_g_per_mol,
        steps=steps, substrate_molar_mass_g_per_mol=substrate_mass)
    return HybridCase(product, environment, stress, trajectory, info, carbon_atoms,
                      couple_stress, grid_points, int(rates.feasible.sum()), len(rates.feasible))


# The retained legacy name for the comparison route above. It is not, and must not become,
# an alias for ``simulate_protocol``: the two run different physics.
simulate_hybrid_case = simulate_hybrid_comparison
