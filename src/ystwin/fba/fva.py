"""Flux-variability analysis focused on one product.

PARKED: no strain carries the pathway yet. See docs/PARKED.md.

What matters is the *width* of the feasible product flux at a realistic growth rate: a range
spanning its own maximum means the model bounds the product rather than predicting it.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra
import pandas as pd

__all__ = ["FluxRange", "capacity_sweep", "product_flux_range"]

from .physiology import cap_uptake
from .solver import FVA_PROCESSES, configure as pin_solver, growth_or_none

GLUCOSE_EXCHANGE = "r_1714"
BIOMASS_REACTION = "r_2111"


@dataclass(frozen=True)
class FluxRange:
    """Feasible interval for one reaction under a fixed growth requirement."""

    reaction: str
    minimum: float
    maximum: float
    max_growth: float
    enforced_growth: float
    glucose_uptake: float
    pathway_capacity: float | None

    @property
    def width(self) -> float:
        return self.maximum - self.minimum

    @property
    def relative_width(self) -> float:
        """Width as a fraction of the maximum. 1.0 means fully undetermined.

        **NaN where ``maximum <= 0``, and that NaN is a trap for any caller that COMPARES
        it.** Every comparison against NaN is False, so a NaN row is silently counted as
        "did not widen", "did not narrow" and "is not zero" all at once -- it never demotes
        anything and never equals anything. `scripts/run_regulation.py` does exactly such a
        comparison when it builds its `widened` mask, so rows with a non-positive ceiling
        fall out of that count without being reported as missing.

        The NaN itself is defensible: with a ceiling at or below zero there is no positive
        scale to divide by, and returning 0.0 would say "fully determined" about a range
        nobody can normalise. What is not defensible is comparing it without checking.
        Callers must test ``math.isnan`` explicitly.

        This was NOT changed to a different return value on 2026-09-02, when it was found,
        and the reason is worth recording: the branch is live in committed artefacts --
        `outputs/regulation_depth_sweep.csv` has 1 NaN row of 8, and
        `outputs/regulation_product_range.csv` has 2 and 3 of 28 in its two width columns --
        so any new value would move numbers that are cited elsewhere. Changing it is a
        deliberate call with a re-run behind it, not a tidy-up.
        """
        if self.maximum <= 0:
            return float("nan")
        return self.width / self.maximum

    def summary(self) -> str:
        return (
            f"{self.reaction}: [{self.minimum:.4g}, {self.maximum:.4g}] "
            f"mmol/gDW/h at mu={self.enforced_growth:.4f} "
            f"({self.enforced_growth / self.max_growth:.0%} of max); "
            f"relative width {self.relative_width:.3f}"
        )


def _configure(model, glucose_uptake, pathway_capacity, product_reaction):
    # An ec model sets its own uptake; forcing a bound overrides it.
    if glucose_uptake is not None:
        cap_uptake(model, GLUCOSE_EXCHANGE, glucose_uptake)
    if pathway_capacity is not None:
        model.reactions.get_by_id(product_reaction).upper_bound = pathway_capacity


def product_flux_range(
    model: cobra.Model,
    product_reaction: str,
    glucose_uptake: float | None = 10.0,
    growth_fraction: float = 0.9,
    pathway_capacity: float | None = None,
    biomass_reaction: str = BIOMASS_REACTION,
) -> FluxRange:
    """Feasible range of the product flux with growth held at a fraction of its max.

    Args:
        model: Model carrying the product reaction. Its constraints are not
            modified; its solver and tolerance are pinned (see :mod:`ystwin.fba.solver`),
            because the interval this returns is exactly the case where which optimum a
            solver happens to pick is arbitrary.
        product_reaction: Reaction whose range is wanted.
        glucose_uptake: Magnitude of the glucose uptake bound, mmol/gDW/h.
            ``None`` leaves the model's own uptake constraint in place.
        growth_fraction: Fraction of maximum growth to enforce, in (0, 1].
        pathway_capacity: Optional ceiling on the product reaction, standing in
            for an effective in-vivo enzyme capacity.
        biomass_reaction: Growth reaction id.
    """
    if not 0.0 < growth_fraction <= 1.0:
        raise ValueError(f"growth_fraction must be in (0, 1], got {growth_fraction}")

    with model as m:
        pin_solver(m)
        _configure(m, glucose_uptake, pathway_capacity, product_reaction)
        m.objective = biomass_reaction
        max_growth = growth_or_none(m)
        if max_growth is None or max_growth <= 0:
            raise ValueError("model cannot grow under these constraints")

        enforced = max_growth * growth_fraction
        m.reactions.get_by_id(biomass_reaction).lower_bound = enforced
        fva = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=[product_reaction], fraction_of_optimum=0.0,
            processes=FVA_PROCESSES,
        )

    return FluxRange(
        reaction=product_reaction,
        minimum=float(fva.loc[product_reaction, "minimum"]),
        maximum=float(fva.loc[product_reaction, "maximum"]),
        max_growth=float(max_growth),
        enforced_growth=float(enforced),
        glucose_uptake=float("nan") if glucose_uptake is None else float(glucose_uptake),
        pathway_capacity=pathway_capacity,
    )


def capacity_sweep(
    model: cobra.Model,
    product_reaction: str,
    capacities,
    glucose_uptake: float | None = 10.0,
    growth_fraction: float = 0.9,
) -> pd.DataFrame:
    """Product range across candidate effective pathway capacities.

    This is the quantitative form of the tiered capacity policy: it shows how much
    of the prediction is coming from stoichiometry and how much from the capacity
    bound that has to be estimated from data.
    """
    rows = []
    for cap in capacities:
        rng = product_flux_range(
            model, product_reaction, glucose_uptake=glucose_uptake,
            growth_fraction=growth_fraction, pathway_capacity=cap,
        )
        rows.append({
            "capacity": cap,
            "minimum": rng.minimum,
            "maximum": rng.maximum,
            "width": rng.width,
            "relative_width": rng.relative_width,
            "max_growth": rng.max_growth,
        })
    return pd.DataFrame(rows)
