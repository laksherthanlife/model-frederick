"""Experiment design: break the collinearity between stress and growth rate.

Every perturbation run so far slows growth, so stress and 1/mu are one axis and no latent
stress state is separable from growth at any sample size. The fix is a design: nutrient
limitation slows growth without stressing, sub-inhibitory doses stress without slowing
growth, and crossing them fills the off-diagonal cells.

The environment is a third set of axes
--------------------------------------
Dose and nutrient limitation are two knobs on one plate in one incubator. Carbon source,
temperature, aeration and pH move the same growth rate for reasons that have nothing to do
with stress, and a model that has only ever seen 30-degree aerobic glucose reads any of them
as stress. :func:`environment_grid` crosses the environment with the dose ladder and records
every axis on every row, so a downstream fit can be conditioned on the conditions instead of
having to infer them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

import numpy as np

from .context import CultureContext, baseline_activity, context_growth_rate
from .culture import ChemostatPhysiology, CultureParameters, chemostat_physiology

__all__ = [
    "Condition",
    "EnvironmentCondition",
    "collinearity",
    "context_grid",
    "decoupling_grid",
    "environment_grid",
    "stressor_only_series",
]


@dataclass(frozen=True)
class Condition:
    """One planned culture, with the two axes it is meant to separate."""

    dose_mM: float
    nutrient_factor: float
    growth_rate: float
    promoter_activity: float


def _axis(name: str, values) -> tuple:
    """Materialise one design axis, refusing an empty one.

    Every builder below is a cross product, and a cross product with one empty factor is
    empty. That is the worst possible failure here: a sweep handed a mistyped or filtered-out
    axis runs over zero conditions, writes a valid empty table and exits 0, and the run reads
    as "swept, nothing happened" rather than "never swept". An axis with no values is a design
    error and is refused where it enters, not where the emptiness eventually shows up.
    """
    materialised = tuple(values)
    if not materialised:
        raise ValueError(
            f"design axis {name!r} is empty, so the cross product would be empty and the "
            "sweep would report zero conditions as a successful run; give it at least one "
            "value")
    return materialised


def _condition(params: CultureParameters, dose: float, nutrient: float) -> Condition:
    return Condition(
        dose_mM=float(dose), nutrient_factor=float(nutrient),
        growth_rate=float(params.growth_rate_at(dose) * nutrient),
        promoter_activity=float(params.promoter_activity_at(dose)),
    )


def stressor_only_series(params: CultureParameters, doses) -> list[Condition]:
    """A dose ladder at full nutrients: the design already run, and a collinear one."""
    return [_condition(params, d, 1.0) for d in _axis("doses", doses)]


def decoupling_grid(params: CultureParameters, doses, nutrient_factors) -> list[Condition]:
    """Cross the dose ladder with nutrient limitation to fill all four quadrants."""
    return [_condition(params, d, n)
            for n in _axis("nutrient_factors", nutrient_factors)
            for d in _axis("doses", doses)]


@dataclass(frozen=True)
class EnvironmentCondition:
    """One planned culture, with its whole environment recorded beside its two axes.

    Everything a downstream fit could want to condition on, and nothing inferred. The point of
    carrying the context rather than a summary of it is that two contexts can produce the same
    growth rate for different reasons -- a cool aerobic glucose culture and a warm
    poorly-aerated one -- and a model handed only the rate cannot tell them apart.

    Args:
        context: The physiological baseline, verbatim.
        dose_mM: Stressor dose.
        nutrient_factor: Growth-only limitation, the axis that is not stress.
        growth_rate: Realised rate, 1/h, after context, dose and nutrient limitation.
        promoter_activity: Reporter synthesis rate at this dose.
        baseline_modules: Module activities the context alone produces, before dosing.
        physiology: Measured carbon and oxygen fluxes at this growth rate, or ``None``.
        physiology_refusal: Why ``physiology`` is ``None``. Empty when it is not.
    """

    context: CultureContext
    dose_mM: float
    nutrient_factor: float
    growth_rate: float
    promoter_activity: float
    baseline_modules: dict[str, float] = field(repr=False, default_factory=dict)
    physiology: ChemostatPhysiology | None = field(repr=False, default=None)
    physiology_refusal: str = ""

    def to_row(self) -> dict[str, float | str | bool | None]:
        """One flat record: every condition, then every outcome.

        Flat because the consumer is a dataframe and a fit, and a nested column is a column
        somebody has to unpack before they can condition on it. Refusals travel as data --
        ``physiology_refusal`` is a column, not an exception -- so that a sweep over a wide
        environment reports which corners it could not ground instead of stopping at the first.
        """
        physiology = self.physiology
        row: dict[str, float | str | bool | None] = {
            "carbon_source": self.context.carbon_source,
            "growth_phase": self.context.growth_phase,
            "temperature_c": float(self.context.temperature_c),
            "oxygen": float(self.context.oxygen),
            "anaerobic_supplements": bool(self.context.anaerobic_supplements),
            "glucose_g_per_L": float(self.context.glucose_g_per_L),
            "ph_medium": float(self.context.ph_medium),
            "ph_cytosolic": float(self.context.ph_cytosolic),
            "strain": self.context.strain,
            "dose_mM": float(self.dose_mM),
            "nutrient_factor": float(self.nutrient_factor),
            "growth_rate": float(self.growth_rate),
            "promoter_activity": float(self.promoter_activity),
            "physiology_refusal": self.physiology_refusal,
            "fermentative": None if physiology is None else physiology.fermentative,
            "biomass_yield_g_per_g": None if physiology is None else physiology.biomass_yield,
            "q_glucose": None if physiology is None else physiology.glucose_uptake,
            "q_oxygen": None if physiology is None else physiology.oxygen_uptake,
            "q_ethanol": None if physiology is None else physiology.ethanol,
            "q_co2": None if physiology is None else physiology.co2_production,
        }
        row.update({f"baseline_{name}": value
                    for name, value in self.baseline_modules.items()})
        return row


def _environment_condition(
    params: CultureParameters, context: CultureContext, dose: float, nutrient: float
) -> EnvironmentCondition:
    """Compose the three growth axes, and ground the physiology or say why not.

    The composition is multiplicative: the context sets the unstressed rate, the dose removes a
    *fraction* of it, and nutrient limitation removes another. Taking the dose's effect as a
    fraction rather than as an absolute rate is what makes it transferable between contexts, and
    it is the same position `context.py` already takes -- context shifts the baseline and the
    growth rate, not the response. **Asserted.** Nothing here measures whether a stressor's
    fractional inhibition is the same at 37 degrees as at 30, and it is the assumption a sweep
    across environments is most likely to break.
    """
    unstressed = context_growth_rate(context)
    inhibition = params.growth_rate_at(dose) / params.mu_max
    growth = unstressed * inhibition * float(nutrient)
    try:
        physiology = chemostat_physiology(growth)
        refusal = ""
    except ValueError as reason:
        physiology, refusal = None, str(reason)
    return EnvironmentCondition(
        context=context,
        dose_mM=float(dose),
        nutrient_factor=float(nutrient),
        growth_rate=float(growth),
        promoter_activity=float(params.promoter_activity_at(dose)),
        baseline_modules=baseline_activity(context),
        physiology=physiology,
        physiology_refusal=refusal,
    )


def context_grid(
    carbon_sources=("glucose",),
    growth_phases=("exponential",),
    temperatures_c=(30.0,),
    oxygen_fractions=(0.21,),
    ph_medium=(5.5,),
    ph_cytosolic=(7.0,),
    glucose_g_per_L=(20.0,),
    anaerobic_supplements=(False,),
    strain: str = "BY4741",
) -> list[CultureContext]:
    """Full factorial over the environment axes the generator supports.

    Every default is the corner every simulated plate in this repository was generated in, so
    calling this with no arguments returns exactly that one context. Widening one argument
    widens the design along one axis and leaves the rest pinned, which is what makes a sweep
    attributable.

    Raises:
        ValueError: from ``CultureContext`` for any combination it has no physiology for. Not
            caught here: an unsupported axis value is a design error, and returning a smaller
            grid than was asked for would hide it.
    """
    axes = product(
        _axis("carbon_sources", carbon_sources), _axis("growth_phases", growth_phases),
        _axis("temperatures_c", temperatures_c), _axis("oxygen_fractions", oxygen_fractions),
        _axis("anaerobic_supplements", anaerobic_supplements),
        _axis("glucose_g_per_L", glucose_g_per_L), _axis("ph_medium", ph_medium),
        _axis("ph_cytosolic", ph_cytosolic))
    return [
        CultureContext(
            carbon_source=carbon, growth_phase=phase, temperature_c=float(temperature),
            oxygen=float(oxygen), anaerobic_supplements=bool(supplements),
            glucose_g_per_L=float(glucose), ph_medium=float(medium_ph),
            ph_cytosolic=float(cytosolic_ph), strain=strain,
        )
        for carbon, phase, temperature, oxygen, supplements, glucose, medium_ph, cytosolic_ph
        in axes
    ]


def environment_grid(
    params: CultureParameters, contexts, doses, nutrient_factors=(1.0,)
) -> list[EnvironmentCondition]:
    """Cross the environment with the dose ladder and the nutrient axis.

    The order of the loops is the order a wet lab would run them -- one environment is one
    incubator setting, and the doses inside it are one plate -- so a truncated sweep is a
    coherent subset of plates rather than a scatter of half-finished ones.
    """
    return [
        _environment_condition(params, context, dose, nutrient)
        for context in _axis("contexts", contexts)
        for nutrient in _axis("nutrient_factors", nutrient_factors)
        for dose in _axis("doses", doses)
    ]


def collinearity(conditions: list[Condition]) -> float:
    """|correlation| between stress and growth across a design; 1 means unidentifiable.

    Takes anything carrying ``promoter_activity`` and ``growth_rate``, which is both
    :class:`Condition` and :class:`EnvironmentCondition` -- the question it answers is about the
    two numbers, and an environment axis is worth measuring on exactly the same scale as a
    nutrient axis.
    """
    if len(conditions) < 3:
        raise ValueError("need at least 3 conditions to measure collinearity")
    stress = np.array([c.promoter_activity for c in conditions], dtype=float)
    growth = np.array([c.growth_rate for c in conditions], dtype=float)
    if np.std(stress) < 1e-12 or np.std(growth) < 1e-12:
        return 0.0
    return float(abs(np.corrcoef(stress, growth)[0, 1]))
