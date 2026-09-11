"""Audit a predicted product flux against the genome-scale model. Never supply one.

This is the layer `predict.py`'s docstring has always described and never called. It exists
because the GEM has a real job here and it is not the job it was being credited with.

**What it cannot do, established rather than assumed.** Capping every reaction of the
installed beta-carotene pathway at exactly 1e-3, then at 1e-4 -- the measured flux
magnitude -- gives FVA ranges of [0, 1e-3] and [0, 1e-4] on Yeast9 v9.0.2. The ceiling
becomes exactly the cap; the FLOOR NEVER LEAVES ZERO. Same on ecYeastGEM_batch. That is not
a parameterisation failure that a better kcat would fix: a capacity constraint is an upper
bound, an upper bound cannot make a flux mandatory, and GECKO's ``v <= kcat*[E]`` under a
shared protein pool has the zero vector feasible too (PMID 28779005). No FBA formulation
supplies a heterologous flux from environment alone.

**What it can do**, and all three are checks on a number that arrives from elsewhere:

1. FEASIBILITY -- is the predicted flux inside what the network could carry at this growth
   rate and uptake? A prediction above the ceiling is refuted by stoichiometry alone.
2. PRECURSOR BUDGET -- the native drain on the precursor node is a floor the pathway has to
   out-compete. Reported, not enforced, because the floor itself is conditional.
3. GROWTH COST -- what carrying this flux costs in growth rate.

**The growth cost is regime-dependent and this module makes that unavoidable.** The same
pathway at the same flux costs 0.366% of growth under a glucose bound of 10 and 3.500%
under a measured uptake of 1.22 -- a ten-fold swing driven entirely by an input. So
``glucose_uptake`` is a required argument here rather than a defaulted one, and it is
carried out in the result. A burden number without its uptake regime is not a number.

**Where regulation belongs, and why it is here rather than in the prediction.** E-Flux
(Colijn 2009, PMID 19714220) turns a stress module's activity into scaled reaction upper
bounds. Wired into the PREDICTION it would be inert -- scaling a ceiling cannot set a flux,
which is the same argument that puts the whole GEM in this module rather than upstream of
it. Wired into the AUDIT it is not inert: stress-scaled capacities tighten the envelope the
prediction is checked against, so a flux that clears the unstressed envelope by 100x may not
clear the stressed one. That is a real question and it is the one this layer can answer.

Pass ``module_activity`` to ask it. Leaving it out asks the unstressed question, and the
result records which was asked.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra

from .fva import GLUCOSE_EXCHANGE, product_flux_range
from .physiology import cap_uptake
from .solver import FVA_PROCESSES, configure as pin_solver, growth_or_none

__all__ = ["FluxAudit", "audit_predicted_flux", "precursor_floor"]


@dataclass(frozen=True)
class FluxAudit:
    """What the GEM says about a flux it did not produce.

    Args:
        predicted_flux: The number being audited, mmol/gDCW/h.
        feasible_min: FVA minimum for the product reaction. Expected to be 0.
        feasible_max: FVA maximum -- the stoichiometric ceiling.
        within_envelope: Whether the prediction sits inside [min, max].
        headroom: ``feasible_max / predicted_flux``, how far below the ceiling it sits.
        precursor_floor: Minimum flux the native competition forces through the precursor
            node, mmol/gDCW/h, or ``None`` when it was not asked for.
        growth_cost_fraction: Fractional drop in maximum growth from carrying the flux.
        glucose_uptake: The uptake bound this was computed under. Load-bearing: see the
            module docstring.
        max_growth: Maximum growth with the product flux at zero.
        regulation: How the envelope was constrained -- ``"unstressed"`` when no module
            activity was supplied, otherwise a summary of the E-Flux layer applied. Carried
            because two audits of the same flux are not comparable unless this matches.
        notes: Anything a reader needs.
    """

    predicted_flux: float
    feasible_min: float
    feasible_max: float
    within_envelope: bool
    headroom: float
    precursor_floor: float | None
    growth_cost_fraction: float
    glucose_uptake: float
    max_growth: float
    regulation: str = "unstressed"
    notes: tuple[str, ...] = ()

    def summary(self) -> str:
        verdict = "inside" if self.within_envelope else "OUTSIDE"
        return (f"{self.predicted_flux:.4g} is {verdict} the feasible envelope "
                f"[{self.feasible_min:.4g}, {self.feasible_max:.4g}] "
                f"({self.headroom:.0f}x headroom) at glucose <= {self.glucose_uptake:g}; "
                f"carrying it costs {self.growth_cost_fraction:.2%} of growth")


def precursor_floor(
    model: cobra.Model,
    precursor_metabolite: str,
    growth_fraction: float = 0.9,
    biomass_reaction: str = "r_2111",
) -> float:
    """Minimum total flux the network is forced to push through a precursor node.

    Every irreversible consumer of the metabolite is FVA-minimised with growth held at a
    fraction of maximum, and the minima are summed. What comes back is the demand the host
    itself makes on the node -- ergosterol and dolichol off FPP, for instance -- which the
    heterologous pathway competes with rather than adds to.

    It is REPORTED and never enforced, and the reason is worth stating: the floor is a
    property of the medium as much as of the network. Open a sterol uptake reaction and the
    FPP floor collapses to exactly zero, because the cell stops having to make its own. A
    number that moves to zero on a medium change is not a constraint to build on.

    Raises:
        KeyError: if the metabolite is not in the model, naming it rather than returning a
            zero that would read as "no native demand".
    """
    if precursor_metabolite not in {m.id for m in model.metabolites}:
        # GECKO exports suffix the compartment; try the bare id both ways before failing.
        alternatives = [m.id for m in model.metabolites
                        if m.id.split("[")[0] == precursor_metabolite]
        if not alternatives:
            raise KeyError(
                f"no metabolite {precursor_metabolite!r} in this model; the precursor "
                "node named by the pathway spec has to exist before its budget can be read")
        precursor_metabolite = alternatives[0]

    metabolite = model.metabolites.get_by_id(precursor_metabolite)
    consumers = [r.id for r in metabolite.reactions
                 if r.metabolites[metabolite] < 0 and r.lower_bound >= 0]
    if not consumers:
        return 0.0

    with model as m:
        pin_solver(m)
        m.objective = biomass_reaction
        max_growth = growth_or_none(m)
        if max_growth is None or max_growth <= 0:
            raise ValueError("model cannot grow, so no floor can be read from it")
        m.reactions.get_by_id(biomass_reaction).lower_bound = max_growth * growth_fraction
        fva = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=consumers, fraction_of_optimum=0.0, processes=FVA_PROCESSES)
    return float(fva["minimum"].clip(lower=0.0).sum())


def audit_predicted_flux(
    model: cobra.Model,
    product_reaction: str,
    predicted_flux: float,
    glucose_uptake: float,
    precursor_metabolite: str | None = None,
    module_activity: dict[str, float] | None = None,
    growth_fraction: float = 0.9,
    biomass_reaction: str = "r_2111",
) -> FluxAudit:
    """Check a flux the GEM did not produce.

    Args:
        model: Model carrying the installed pathway. Not modified.
        product_reaction: The reaction whose flux is being audited.
        predicted_flux: From the pathway flux layer, mmol/gDCW/h.
        glucose_uptake: Uptake bound magnitude. **Required**, because the growth cost
            swings ten-fold with it and a burden number without its regime is meaningless.
        precursor_metabolite: Optional, for the budget check.
        module_activity: Optional stress module activities, from
            ``generator.stress_panel.module_response``. When given, E-Flux scales the
            capacity of every reaction those modules reach and the envelope is computed
            under those bounds -- so the question becomes "does this flux fit in a STRESSED
            cell". Omit it to ask the unstressed question.
        growth_fraction: Growth held at this fraction of maximum during FVA.
        biomass_reaction: Growth reaction id.

    Raises:
        ValueError: on a negative flux or a non-positive uptake.
    """
    if predicted_flux < 0:
        raise ValueError(f"predicted_flux must be non-negative, got {predicted_flux}")
    if glucose_uptake <= 0:
        raise ValueError(
            f"glucose_uptake must be positive, got {glucose_uptake}. It is required rather "
            "than defaulted because the growth cost below is a property of it")

    notes: list[str] = []
    regulation = "unstressed"

    if module_activity:
        # Regulation narrows the envelope; it cannot move the prediction. Applied inside a
        # context so the caller's model is unchanged, and the audit is then run against the
        # constrained copy rather than against the original.
        from ..bridge.regulation import apply_eflux, eflux_layer
        from ..generator.stress_panel import METABOLITE_POOLS

        # The metabolite pools -- redox, peroxide, atp, ph, nadh -- carry activity and are
        # not transcriptional, so they have no regulon and E-Flux refuses them rather than
        # dropping them silently. That refusal is right and it is the caller's job to
        # answer it, which is what this does: they are excluded, deliberately, and the
        # exclusion is reported on the result rather than buried.
        with model as constrained:
            pin_solver(constrained)
            layer = eflux_layer(constrained, module_activity,
                                allow_unmapped=METABOLITE_POOLS)
            if layer.binding:
                apply_eflux(constrained, layer, biomass_reaction=biomass_reaction,
                            acknowledge_infeasible=True)
                regulation = layer.summary()
                if layer.modules_excluded:
                    notes.append(
                        f"excluded from the E-Flux layer as non-transcriptional: "
                        f"{list(layer.modules_excluded)}. These are metabolite pools with "
                        "no regulon; they can still be read, they just cannot scale a "
                        "reaction bound")
                notes.append(
                    f"envelope computed under E-Flux with {len(layer.binding)} reaction "
                    f"bound(s) tightened by the supplied module activities; the ceiling "
                    "moves and the floor cannot, because E-Flux scales upper bounds only")
            else:
                regulation = "no reaction reached"
                notes.append(
                    "module activities were supplied and E-Flux tightened NOTHING -- every "
                    "gene those modules reach is either absent from this model or already "
                    "at its bound. The envelope below is the unstressed one")
            envelope = product_flux_range(
                constrained, product_reaction, glucose_uptake=glucose_uptake,
                growth_fraction=growth_fraction, biomass_reaction=biomass_reaction)
    else:
        envelope = product_flux_range(
            model, product_reaction, glucose_uptake=glucose_uptake,
            growth_fraction=growth_fraction, biomass_reaction=biomass_reaction)

    if envelope.minimum > 1e-9:
        notes.append(
            f"the feasible minimum is {envelope.minimum:.3g}, not zero. Something in this "
            "model makes the product mandatory, which is worth understanding before "
            "trusting the ceiling")

    with model as m:
        pin_solver(m)
        # THE SAME UPTAKE AS THE ENVELOPE. Without this line the cost is computed against
        # the model's own default uptake while `unconstrained` came from `glucose_uptake`,
        # and the subtraction compares two different cultures: it reported 90% growth cost
        # for a flux of 9.5e-4 mmol/gDCW/h, which is four orders of magnitude below the
        # ceiling it had just cleared by 106x.
        cap_uptake(m, GLUCOSE_EXCHANGE, glucose_uptake)
        m.reactions.get_by_id(biomass_reaction).lower_bound = 0.0
        product = m.reactions.get_by_id(product_reaction)
        product.bounds = (predicted_flux, predicted_flux)
        m.objective = biomass_reaction
        constrained = growth_or_none(m)
    unconstrained = envelope.max_growth
    if constrained is None:
        cost = float("nan")
        notes.append(
            f"the model is infeasible with the product pinned at {predicted_flux:.4g}, so "
            "the prediction is refuted by stoichiometry rather than merely expensive")
    else:
        cost = (unconstrained - float(constrained)) / unconstrained

    within = envelope.minimum - 1e-9 <= predicted_flux <= envelope.maximum + 1e-9
    if not within:
        notes.append(
            f"predicted {predicted_flux:.4g} lies outside [{envelope.minimum:.4g}, "
            f"{envelope.maximum:.4g}]; the kinetic layer is asking for a flux the network "
            "cannot carry at this growth rate and uptake")

    floor = None
    if precursor_metabolite is not None:
        floor = precursor_floor(model, precursor_metabolite, growth_fraction,
                                biomass_reaction)
        notes.append(
            f"native demand on {precursor_metabolite} is at least {floor:.4g} mmol/gDCW/h "
            "under this medium; open a sterol uptake and it collapses to zero, so this is "
            "context and not a constraint")

    headroom = (envelope.maximum / predicted_flux) if predicted_flux > 0 else float("inf")
    return FluxAudit(
        predicted_flux=float(predicted_flux),
        feasible_min=envelope.minimum,
        feasible_max=envelope.maximum,
        within_envelope=within,
        headroom=headroom,
        precursor_floor=floor,
        growth_cost_fraction=cost,
        glucose_uptake=float(glucose_uptake),
        max_growth=unconstrained,
        regulation=regulation,
        notes=tuple(notes),
    )
