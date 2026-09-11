"""Branch A of g_psi: measured physiology -> constraints -> flux. Load-bearing.

Inputs are quantities the plate measures. Uptake is bounded by supply, never by what the measured
rate demands, so the model keeps the ability to call a measurement impossible. Fluorescence is
not an input, so G4's verdict does not gate this branch.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cobra
import pandas as pd
from cobra.exceptions import OptimizationError

from ..fba.solver import growth_or_none

__all__ = [
    "MeasuredState",
    "PhysiologicalConstraints",
    "FluxPrediction",
    "physiology_constraints",
    "predict_flux",
    "predict_flux_trajectory",
]

BIOMASS_REACTION = "r_2111"
GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"

# Uptake must be bounded by supply, never by what the measured rate demands.
#
# 22.0 was attributed to "van Hoek et al. 1998 aerobic batch". That number IS in PMID
# 9797269 and it is a different quantity: "the fermentative capacity showed only a small
# further increase, up to 22.0 mmol of ETHANOL per g dry yeast biomass per h at D = 0.40",
# measured offline by moving chemostat-grown cells to ANAEROBIC conditions under CO2 with
# 2% glucose. So it was an anaerobic ethanol production capacity read as an aerobic
# glucose uptake ceiling -- wrong metabolite, wrong gas regime.
#
# The highest in-situ glucose uptake that paper reports is 11.1 mmol/gDW/h, at D = 0.40.
# A Vmax must sit above the highest observed rate rather than at it, and no measured yeast
# glucose Vmax was found in either van Hoek 1998 paper, so this is ASSERTED at twice the
# highest observed uptake and labelled as such. It is a supply ceiling whose only job is
# to stop the bound being set by demand; it is not a fitted kinetic parameter.
_GLUCOSE_QMAX = 22.2   # mmol/gDW/h, ASSERTED: 2 x van Hoek's highest observed 11.1
_GLUCOSE_KM_MM = 0.5


@dataclass(frozen=True)
class MeasuredState:
    """What the plate actually measured at one time.

    Deliberately contains no latent quantity. Anything inferred from fluorescence
    belongs to Branch B, and keeping it out of this type is what makes the
    separation checkable rather than merely intended.

    Args:
        biomass_gl: Biomass, g/L.
        growth_rate: Specific growth rate, 1/h. May be negative for a dying culture.
        glucose_mM: Residual glucose, if assayed.
        time_h: Time of the measurement.
    """

    biomass_gl: float
    growth_rate: float
    glucose_mM: float
    time_h: float


@dataclass(frozen=True)
class PhysiologicalConstraints:
    """Pre-solve bounds derived from measurement alone."""

    growth_lower_bound: float
    glucose_uptake: float
    oxygen_uptake: float | None
    time_h: float
    branch: str = "physiology"

    def summary(self) -> str:
        return (
            f"[{self.branch}] t={self.time_h:.2f} h  growth >= {self.growth_lower_bound:.4f}/h, "
            f"glucose <= {self.glucose_uptake:.3f} mmol/gDW/h"
            + (f", oxygen <= {self.oxygen_uptake:.3f}" if self.oxygen_uptake else "")
        )


@dataclass(frozen=True)
class FluxPrediction:
    """Fluxes solved under one set of constraints, with their provenance."""

    fluxes: dict[str, float]
    constraints: PhysiologicalConstraints
    branch: str
    max_growth: float
    ranges: dict[str, tuple[float, float]] = field(default_factory=dict)


def physiology_constraints(
    state: MeasuredState,
    oxygen_uptake: float | None = None,
    glucose_qmax: float = _GLUCOSE_QMAX,
    glucose_km_mM: float = _GLUCOSE_KM_MM,
) -> PhysiologicalConstraints:
    """Turn a measured state into pre-solve bounds.

    Substrate uptake is bounded by supply, not by demand:

        q_glucose <= q_max * S / (Km + S)

    so the model retains the ability to declare a measured growth rate impossible.
    That check is the point of the branch: it is what makes a flux prediction a
    claim rather than a restatement of the measurement.

    Args:
        state: Measured biomass, growth rate and residual glucose.
        oxygen_uptake: Measured oxygen uptake, if the platform provides one. A plate
            reader does not, so this stays ``None`` for H1 data and the model's own
            limit applies -- which is why no oxygen-linked claim rests on it.
        glucose_qmax: Maximum specific uptake, mmol/gDW/h.
        glucose_km_mM: Half-saturation for uptake.
    """
    availability = (
        state.glucose_mM / (glucose_km_mM + state.glucose_mM) if state.glucose_mM > 0 else 0.0
    )
    return PhysiologicalConstraints(
        growth_lower_bound=float(state.growth_rate),
        glucose_uptake=float(glucose_qmax * availability),
        oxygen_uptake=oxygen_uptake,
        time_h=float(state.time_h),
    )


def predict_flux(
    model: cobra.Model,
    constraints: PhysiologicalConstraints,
    reactions: list[str],
    with_ranges: bool = False,
    biomass_reaction: str = BIOMASS_REACTION,
) -> FluxPrediction:
    """Solve the model under measured-physiology constraints. Does not modify ``model``.

    Args:
        model: GSMM to solve.
        constraints: Bounds from :func:`physiology_constraints`.
        reactions: Reactions whose flux is wanted.
        with_ranges: Also run FVA, so the report can show what the solve leaves free
            rather than a point value that looks more determined than it is.
        biomass_reaction: Growth reaction id.

    Raises:
        ValueError: if the model cannot support the measured growth rate.
    """
    with model as m:
        m.reactions.get_by_id(GLUCOSE_EXCHANGE).lower_bound = -abs(constraints.glucose_uptake)
        if constraints.oxygen_uptake is not None:
            m.reactions.get_by_id(OXYGEN_EXCHANGE).lower_bound = -abs(constraints.oxygen_uptake)

        m.objective = biomass_reaction
        # `growth_or_none`, not `slim_optimize`: an infeasible LP returns nan and every
        # comparison against nan is False, so the refusal below could never fire.
        feasible_growth = growth_or_none(m)
        if feasible_growth is None or feasible_growth < constraints.growth_lower_bound - 1e-6:
            raise ValueError(
                f"the model cannot support the measured growth rate "
                f"{constraints.growth_lower_bound:.4f}/h under these constraints "
                f"(its maximum is {feasible_growth if feasible_growth else 0.0:.4f}/h); "
                "either the uptake bound is too tight or the measured rate is wrong"
            )

        # Fix growth at the measured rate; pFBA picks among the rest.
        growth = m.reactions.get_by_id(biomass_reaction)
        growth.bounds = (constraints.growth_lower_bound, constraints.growth_lower_bound)
        try:
            solution = cobra.flux_analysis.pfba(m)
        except OptimizationError as exc:
            raise ValueError(f"no feasible flux distribution: {exc}") from exc

        fluxes = {r: float(solution.fluxes[r]) for r in reactions}
        ranges: dict[str, tuple[float, float]] = {}
        if with_ranges:
            fva = cobra.flux_analysis.flux_variability_analysis(
                m, reaction_list=reactions, fraction_of_optimum=0.0
            )
            ranges = {
                r: (float(fva.loc[r, "minimum"]), float(fva.loc[r, "maximum"]))
                for r in reactions
            }

    return FluxPrediction(
        fluxes=fluxes, constraints=constraints, branch=constraints.branch,
        max_growth=float(feasible_growth), ranges=ranges,
    )


def predict_flux_trajectory(
    model: cobra.Model,
    states: list[MeasuredState],
    reactions: list[str],
    **kwargs,
) -> pd.DataFrame:
    """Solve once per measured timepoint. One row per time, labelled by branch."""
    rows = []
    for state in states:
        constraints = physiology_constraints(state)
        prediction = predict_flux(model, constraints, reactions, **kwargs)
        rows.append({
            "time_h": state.time_h, "branch": prediction.branch,
            "biomass_gl": state.biomass_gl, "growth_rate": state.growth_rate,
            "glucose_uptake": constraints.glucose_uptake,
            **{f"flux_{r}": v for r, v in prediction.fluxes.items()},
        })
    return pd.DataFrame(rows)
