"""Retained intracellular osmolyte pools -- what the cell HOLDS, not what it exports.

Yeast9 can only let glycerol go: the channel ``r_1172`` (``s_0765 -> s_0766``) and the exchange
``r_1808``, and nothing that lets glycerol sit in the cytosol. So the central event of osmotic
adaptation -- Fps1 shuts and glycerol ACCUMULATES -- is inexpressible, and a solve that prices
"the osmotic response" is pricing an export instead.

A retained pool is a demand on the cytosolic species carrying the pool's growth dilution::

    d[X]/dt = v_retain - mu*[X] = 0    ->    v_retain = mu*[X]

the DILUTED terminal balance of :mod:`ystwin.pathway.solve` (``X = v_in/mu``) and the same
``q = mu * content`` :mod:`ystwin.bridge.stress_diversion` uses on Hakkaart's contents. It is
installed as an equality COUPLING the demand to biomass rather than as a pinned flux, so the
pool holds its content at whatever growth rate pays for it, with no outer iteration.

Measured on Yeast9 at glucose <= 10, oxygen free: holding a pool and exporting the identical
flux cost the same growth to 5e-15, because ``r_1172`` is free and the GEM prices the
SYNTHESIS, not the destination. What Fps1 changes is WHICH flux -- open, the cell pays
``mu*[X] + leak`` for the same ``[X]``; shut, ``mu*[X]``. The leak term is the audit's export
ladder (-3.10% / -6.20% / -31.02% at 0.5 / 1.0 / 5.0 mmol/gDCW/h), so the two compare directly.
A retained pool is also a fixed FRACTION of growth at any uptake where a forced export is not.

WHERE THE DEMAND ITSELF EARNS ITS PLACE, because with Fps1 open it does not. The same
coupling put on the PRE-EXISTING ``r_1808`` reaches 0.799871 /h to 4e-15 on a model carrying
no sinks at all, so installing them changes no answer this repository already had. Shut
``r_1172`` -- Fps1 closed, the event this module is about -- and glycerol can no longer reach
the extracellular compartment: the ``r_1808`` route forces growth to ZERO, while the demand
still holds the pool at 0.799871. That case, and only that case, is what the sink adds.

Every cost here is quoted at ONE cytosolic volume and is not one number. The measured osmotic
pool costs 5.20% / 9.89% / 12.91% across `thermo_gate`'s (1.0, 2.0, 2.7) mL/gDCW, a 2.48x
band, and the volume-free John 2012 anchor lands at 2.62 -- near the TOP of it.

The aerobic caveat is this module's reason to exist. `stress_diversion` retracted its osmotic
glycerol row because Jouhten 2008 (PMID 18613954) measures q_glycerol = 0.00 +/- 0.00 at every
oxygenated level: aerobically there is no glycerol EXPORT to divert. A retained pool is not
that number and is not refuted by it: on one set of samples Petelenz-Kurdziel 2013 measured
the internal pool rising 18.2x while the supernatant rose 1.47x and an unstressed control fell
0.70x. Those four numbers are :data:`FPS1_CLOSURE_MOL_PER_L`, and they are the whole argument.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra

from .solver import growth_or_none

__all__ = [
    "BIOMASS_REACTION",
    "FPS1_CLOSURE_MOL_PER_L",
    "GLYCEROL_C",
    "GLYCEROL_EXCHANGE",
    "MEASURED_POOLS",
    "OsmolytePool",
    "RetentionCost",
    "TREHALOSE_C",
    "add_retained_pool_sinks",
    "cost_of_exporting",
    "fps1_closure_folds",
    "hold_pool",
    "retained_demand_id",
    "retention_flux",
    "retention_versus_export",
]

GLYCEROL_C = "s_0765"
TREHALOSE_C = "s_1520"
GLYCEROL_EXCHANGE = "r_1808"
TREHALOSE_EXCHANGE = "r_1650"
BIOMASS_REACTION = "r_2111"
GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"

#: Cytosolic species a stress response holds against dilution, by Yeast9 metabolite id.
RETAINED_POOL_METABOLITES = (GLYCEROL_C, TREHALOSE_C)


@dataclass(frozen=True)
class OsmolytePool:
    """One MEASURED intracellular pool, in the unit its measurement was made in.

    Two units and never a silent conversion between them. ``concentration_m`` is what the HOG
    literature reports and needs a cytosolic volume to become a content;
    ``content_mmol_per_gdcw`` is what a chemostat carbon balance reports and needs nothing.

    Args:
        name: Short label, used in reports.
        metabolite_id: Cytosolic Yeast9 metabolite the pool is made of.
        molar_mass_g_per_mol: For the mg/gDCW cross-check the literature also quotes.
        condition: The state the number was measured in. Load-bearing: a pool without its
            stress condition is not a number.
        source: Citation. Required -- this module refuses an uncited pool.
        concentration_m: Measured intracellular concentration, mol/L.
        content_mmol_per_gdcw: Measured content per dry weight, mmol/gDCW.
    """

    name: str
    metabolite_id: str
    molar_mass_g_per_mol: float
    condition: str
    source: str
    concentration_m: float | None = None
    content_mmol_per_gdcw: float | None = None

    def content(self, cytosolic_volume_ml_per_gdcw: float | None = None) -> float:
        """The pool in the GEM's unit, mmol/gDCW.

        A molar measurement REQUIRES the volume and does not default one, for the reason
        :data:`ystwin.pathway.thermo_gate.CYTOSOLIC_VOLUMES_ML_PER_GDCW` gives: nothing in
        this repository measures that volume, and its plausible range 1.0-2.7 mL/gDCW is a
        2.7-fold spread that lands directly on the answer. A per-dry-weight measurement
        refuses one, because dividing by a volume it never used would invent a claim.

        Raises:
            ValueError: on a missing volume, a volume supplied where none is meaningful, or
                a pool carrying neither unit or both.
        """
        if (self.concentration_m is None) == (self.content_mmol_per_gdcw is None):
            raise ValueError(
                f"pool {self.name!r} must carry exactly one of concentration_m and "
                "content_mmol_per_gdcw; it is a measurement in one unit, not a conversion")
        if self.content_mmol_per_gdcw is not None:
            if cytosolic_volume_ml_per_gdcw is not None:
                raise ValueError(
                    f"pool {self.name!r} is measured per gram dry weight "
                    f"({self.content_mmol_per_gdcw:g} mmol/gDCW) and needs no cytosolic "
                    "volume; supplying one would divide a measurement by an assumption")
            return float(self.content_mmol_per_gdcw)
        if cytosolic_volume_ml_per_gdcw is None:
            raise ValueError(
                f"pool {self.name!r} is measured as a concentration ({self.concentration_m:g} "
                "M) and becomes a content only in a stated cytosolic volume, mL/gDCW. "
                "ystwin.pathway.thermo_gate puts the plausible range at 1.0-2.7 and this "
                "repository measures none of it, so the volume is an argument, never a default")
        if cytosolic_volume_ml_per_gdcw <= 0:
            raise ValueError(
                f"cytosolic_volume_ml_per_gdcw must be positive, got "
                f"{cytosolic_volume_ml_per_gdcw}")
        # mmol/mL is mol/L, so the conversion is one multiplication -- the same arithmetic
        # ystwin.pathway.thermo_gate.molar_from_content inverts.
        return float(self.concentration_m * cytosolic_volume_ml_per_gdcw)

    def content_mg_per_gdcw(self, cytosolic_volume_ml_per_gdcw: float | None = None) -> float:
        """The pool as a fraction of dry weight, the unit an osmolyte paper quotes."""
        return self.content(cytosolic_volume_ml_per_gdcw) * self.molar_mass_g_per_mol

    def retention_flux(self, growth_rate_per_h: float,
                       cytosolic_volume_ml_per_gdcw: float | None = None) -> float:
        """``mu * [X]``, mmol/gDCW/h: what holding this pool costs in flux."""
        if growth_rate_per_h < 0:
            raise ValueError(f"growth rate must be non-negative, got {growth_rate_per_h}")
        return growth_rate_per_h * self.content(cytosolic_volume_ml_per_gdcw)


#: Pools with a measurement behind them. Nothing enters this table without a citation.
MEASURED_POOLS: dict[str, OsmolytePool] = {
    # Petelenz-Kurdziel 2013 Dataset S2, sheet `glycerol_i`, column "Wild Type4": W303, HPLC,
    # 0.4 M NaCl at t = 3600 s; its Readme puts these in a 50 fL cell, half solute-accessible.
    "glycerol_unstressed": OsmolytePool(
        name="glycerol_unstressed",
        metabolite_id=GLYCEROL_C,
        molar_mass_g_per_mol=92.09,
        condition="W303 in YPD, immediately before the shock (t = 3600 s)",
        source="Petelenz-Kurdziel et al. 2013, PLoS Comput Biol 9:e1003084 (PMID 23762021), "
               "Dataset S2 sheet glycerol_i, Wild Type4",
        concentration_m=0.05483825479084679,
    ),
    "glycerol_osmotic": OsmolytePool(
        name="glycerol_osmotic",
        metabolite_id=GLYCEROL_C,
        molar_mass_g_per_mol=92.09,
        condition="W303 in YPD, 90 min after 0.4 M NaCl (t = 9000 s), the time-course peak",
        source="Petelenz-Kurdziel et al. 2013, PLoS Comput Biol 9:e1003084 (PMID 23762021), "
               "Dataset S2 sheet glycerol_i, Wild Type4",
        concentration_m=0.9968547595808357,
    ),
    # Measured per DRY WEIGHT under salt, so it needs no volume: the check on the molar pool
    # above. Different strain, 1 M NaCl and 72 h, so it brackets rather than calibrates.
    "glycerol_osmotic_dry_weight": OsmolytePool(
        name="glycerol_osmotic_dry_weight",
        metabolite_id=GLYCEROL_C,
        molar_mass_g_per_mol=92.09,
        condition="S. cerevisiae MTCC 2918, 1 M NaCl, 72 h; 240.87 +/- 0.38 mg/gDW",
        source="John, Gayathiri, Rose & Mandal 2012, Curr Microbiol 64:100 (PMID 22038037)",
        content_mmol_per_gdcw=240.87 / 92.09,
    ),
    # Same experiment, sheet `trehalose_i`: 72x smaller than the glycerol pool beside it.
    "trehalose_osmotic": OsmolytePool(
        name="trehalose_osmotic",
        metabolite_id=TREHALOSE_C,
        molar_mass_g_per_mol=342.30,
        condition="W303 in YPD, 3 h after 0.4 M NaCl (t = 14400 s), the time-course peak",
        source="Petelenz-Kurdziel et al. 2013, PLoS Comput Biol 9:e1003084 (PMID 23762021), "
               "Dataset S2 sheet trehalose_i, Wild Type4",
        concentration_m=0.01375812596331335,
    ),
    # Measured per dry weight, so it needs no volume at all -- the anchor the molar pools are
    # checked against. Same source `bridge/stress_diversion` already uses for this content.
    "trehalose_chemostat": OsmolytePool(
        name="trehalose_chemostat",
        metabolite_id=TREHALOSE_C,
        molar_mass_g_per_mol=342.30,
        condition="CEN.PK113-7D, aerobic glucose-limited chemostat, D = 0.025 /h, pH 5.0, "
                  "0.04% CO2; 19.4 +/- 3.7 mg/gDCW",
        source="Hakkaart et al. 2020, Biotechnol Bioeng 117:721 (PMID 31654410, PMC7028085), "
               "Tables 1-2",
        content_mmol_per_gdcw=19.4 / 342.30,
    ),
}


#: One set of Wild Type4 samples inside and outside the cell, mol/L, at the shock
#: (t = 3600 s) and the internal peak (t = 9000 s). Same source as MEASURED_POOLS.
FPS1_CLOSURE_MOL_PER_L: dict[str, float] = {
    "internal_before": 0.05483825479084679,
    "internal_peak": 0.9968547595808357,
    "external_before": 0.0026083179498316866,
    "external_peak": 0.0038337713106743406,
    # Column "No Stress2" at the same two times: the pool does not rise on the clock, so
    # the 18-fold is the shock and not the culture getting older.
    "internal_unstressed_before": 0.05237585671452344,
    "internal_unstressed_peak": 0.03688540108548501,
}
"""Why this module exists, in numbers somebody measured on one set of samples.

Over the 90 minutes after 0.4 M NaCl the internal glycerol pool rises 18.2-fold while the
supernatant rises 1.47-fold, and an unstressed culture sampled at the same two times FALLS
0.70-fold. Glycerol is being HELD, not made and released -- which is the one thing Yeast9
cannot express, because ``r_1172`` and ``r_1808`` only let it go.
"""


def fps1_closure_folds() -> dict[str, float]:
    """Internal, external and unstressed-internal fold change over the same 90 minutes."""
    m = FPS1_CLOSURE_MOL_PER_L
    return {
        "internal": m["internal_peak"] / m["internal_before"],
        "external": m["external_peak"] / m["external_before"],
        "unstressed": m["internal_unstressed_peak"] / m["internal_unstressed_before"],
    }


def retention_flux(content_mmol_per_gdcw: float, growth_rate_per_h: float) -> float:
    """``mu * [X]``, mmol/gDCW/h -- the whole model of a retained pool.

    Free function so a caller with a content from somewhere else does not have to build an
    :class:`OsmolytePool` to use the convention.
    """
    if content_mmol_per_gdcw < 0 or growth_rate_per_h < 0:
        raise ValueError(
            f"content and growth rate must be non-negative, got {content_mmol_per_gdcw} "
            f"and {growth_rate_per_h}")
    return float(content_mmol_per_gdcw * growth_rate_per_h)


def _resolve(model: cobra.Model, base_id: str, compartment: str = "c") -> str:
    """Map a bare Yeast9 metabolite id onto this model's convention, as `carotenoid.py` does."""
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(
        f"osmolyte {base_id!r} not found in {model.id!r} "
        f"(tried {base_id}, {base_id}[{compartment}], {base_id}_{compartment})")


def retained_demand_id(metabolite_id: str) -> str:
    """The demand reaction id for a retained pool of one metabolite."""
    return f"DM_{metabolite_id}_retained"


def add_retained_pool_sinks(
    model: cobra.Model,
    metabolite_ids: tuple[str, ...] = RETAINED_POOL_METABOLITES,
) -> cobra.Model:
    """Return a copy of ``model`` carrying one retained-pool demand per osmolyte.

    The demands are added OPEN, ``(0, 1000)``, and carry no flux on their own: a drain the
    objective does not want sits at zero, which is why installing them cannot change any
    existing answer. :func:`hold_pool` is what makes one carry ``mu*[X]``.

    Args:
        model: Host model. Not modified.
        metabolite_ids: Cytosolic osmolytes to give a pool to.

    Raises:
        ValueError: if a pool demand is already installed.
        KeyError: naming an osmolyte the model does not carry.
    """
    out = model.copy()
    for base_id in metabolite_ids:
        resolved = _resolve(out, base_id)
        rid = retained_demand_id(resolved)
        if out.reactions.has_id(rid):
            raise ValueError(f"retained pool {rid!r} is already installed on this model")
        metabolite = out.metabolites.get_by_id(resolved)
        demand = out.add_boundary(metabolite, type="demand", reaction_id=rid)
        demand.name = f"retained intracellular pool of {metabolite.name}"
        demand.subsystem = "retained osmolyte pool"
    return out


def hold_pool(
    model: cobra.Model,
    metabolite_id: str,
    content_mmol_per_gdcw: float,
    biomass_reaction: str = BIOMASS_REACTION,
) -> None:
    """Couple a pool's demand to growth so its flux is ``mu * content``, in place.

    Call inside ``with model:`` -- the constraint is registered with cobrapy's context
    manager and is undone on exit. Coupling rather than pinning the flux is the point: the
    pool stays at its stated content whatever growth rate the solver ends up at, so the
    answer is the steady state rather than a guess at one.

    Raises:
        KeyError: if the pool demand is absent; run :func:`add_retained_pool_sinks` first.
        ValueError: on a negative content.
    """
    if content_mmol_per_gdcw < 0:
        raise ValueError(f"content must be non-negative, got {content_mmol_per_gdcw}")
    resolved = _resolve(model, metabolite_id)
    demand = model.reactions.get_by_id(retained_demand_id(resolved))
    biomass = model.reactions.get_by_id(biomass_reaction)
    coupling = model.problem.Constraint(
        demand.flux_expression - content_mmol_per_gdcw * biomass.flux_expression,
        lb=0.0, ub=0.0, name=f"retain_{resolved}",
    )
    model.add_cons_vars([coupling])


def _apply_regime(model: cobra.Model, glucose_uptake: float,
                  oxygen_uptake: float | None) -> None:
    """Pin the uptake regime. Required, never defaulted -- see `fba/audit.py`."""
    model.reactions.get_by_id(GLUCOSE_EXCHANGE).lower_bound = -abs(glucose_uptake)
    if oxygen_uptake is not None:
        model.reactions.get_by_id(OXYGEN_EXCHANGE).lower_bound = -abs(oxygen_uptake)


def _growth(model: cobra.Model) -> float:
    # Uses the package guard: slim_optimize signals failure with nan, not None.
    value = growth_or_none(model)
    return 0.0 if value is None else float(value)


@dataclass(frozen=True)
class RetentionCost:
    """Holding a pool beside exporting the identical flux, in one regime.

    Args:
        metabolite: The cytosolic osmolyte held.
        content_mmol_per_gdcw: The pool held against dilution.
        growth_free: Maximum growth with no pool and no export.
        growth_held: Maximum growth with the pool coupled to growth.
        retention_flux: ``growth_held * content``, the flux that holds it.
        growth_exported: Maximum growth with that same flux forced out of the cell.
        leak_flux: Extra export forced ON TOP of the retained pool, an open Fps1.
        growth_held_leaking: Maximum growth holding the pool while leaking.
        glucose_uptake: Uptake bound the whole comparison was run at, mmol/gDCW/h.
        oxygen_uptake: Oxygen bound, or ``None`` for unrestricted.
    """

    metabolite: str
    content_mmol_per_gdcw: float
    growth_free: float
    growth_held: float
    retention_flux: float
    growth_exported: float
    leak_flux: float
    growth_held_leaking: float
    glucose_uptake: float
    oxygen_uptake: float | None

    @property
    def cost_of_holding(self) -> float:
        return 1.0 - self.growth_held / self.growth_free

    @property
    def cost_of_exporting(self) -> float:
        return 1.0 - self.growth_exported / self.growth_free

    @property
    def cost_of_leaking(self) -> float:
        return 1.0 - self.growth_held_leaking / self.growth_free

    def summary(self) -> str:
        return (
            f"{self.content_mmol_per_gdcw:.4g} mmol/gDCW of {self.metabolite} at glucose <= "
            f"{self.glucose_uptake:g}: holding costs {self.cost_of_holding:.2%} "
            f"({self.retention_flux:.4g} mmol/gDCW/h), exporting the same flux costs "
            f"{self.cost_of_exporting:.2%}, holding it with a {self.leak_flux:g} leak costs "
            f"{self.cost_of_leaking:.2%}")


def cost_of_exporting(
    model: cobra.Model,
    flux_mmol_per_gdcw_h: float,
    glucose_uptake: float,
    oxygen_uptake: float | None = None,
    exchange_reaction: str = GLYCEROL_EXCHANGE,
) -> float:
    """Maximum growth with ``exchange_reaction`` forced to carry a flux out of the cell.

    Raises:
        ValueError: on a negative flux. A secretion exchange pinned below zero is an
            UPTAKE: at -5.0 on ``r_1808`` growth reached 1.1398 /h against a free maximum
            of 0.8877 at glucose <= 10, oxygen free -- the model fed on an osmolyte the
            medium never contained.
    """
    if flux_mmol_per_gdcw_h < 0:
        raise ValueError(
            f"export flux must be non-negative, got {flux_mmol_per_gdcw_h}; pinning "
            f"{exchange_reaction!r} negative makes it an uptake and the cost comes back "
            "negative because the model is being fed, not taxed")
    with model as scoped:
        _apply_regime(scoped, glucose_uptake, oxygen_uptake)
        exchange = scoped.reactions.get_by_id(exchange_reaction)
        exchange.bounds = (flux_mmol_per_gdcw_h, flux_mmol_per_gdcw_h)
        return _growth(scoped)


def retention_versus_export(
    model: cobra.Model,
    metabolite_id: str,
    content_mmol_per_gdcw: float,
    glucose_uptake: float,
    oxygen_uptake: float | None = None,
    leak_flux_mmol_per_gdcw_h: float = 0.0,
    exchange_reaction: str = GLYCEROL_EXCHANGE,
    biomass_reaction: str = BIOMASS_REACTION,
) -> RetentionCost:
    """Price a retained pool, the same flux exported, and the pool held with a leak.

    ``model`` must already carry the pool demand (:func:`add_retained_pool_sinks`). Nothing
    here is modified permanently: every solve runs inside ``with model:``.

    Raises:
        ValueError: on a negative leak, for the reason :func:`cost_of_exporting` gives -- a
            leak of -5.0 reported ``cost_of_leaking`` of -16.35%, an open Fps1 scored as a
            growth benefit.
    """
    if leak_flux_mmol_per_gdcw_h < 0:
        raise ValueError(
            f"leak flux must be non-negative, got {leak_flux_mmol_per_gdcw_h}; a negative "
            f"lower bound on {exchange_reaction!r} opens an uptake rather than a leak")
    resolved = _resolve(model, metabolite_id)
    with model as scoped:
        _apply_regime(scoped, glucose_uptake, oxygen_uptake)
        growth_free = _growth(scoped)
    if growth_free <= 0.0:
        raise ValueError(
            f"the model does not grow at glucose <= {glucose_uptake:g}, oxygen "
            f"{oxygen_uptake}; every cost here is a fraction of that growth and there is none")
    with model as scoped:
        _apply_regime(scoped, glucose_uptake, oxygen_uptake)
        hold_pool(scoped, resolved, content_mmol_per_gdcw, biomass_reaction)
        growth_held = _growth(scoped)
    flux = retention_flux(content_mmol_per_gdcw, growth_held)
    growth_exported = cost_of_exporting(
        model, flux, glucose_uptake, oxygen_uptake, exchange_reaction)
    with model as scoped:
        _apply_regime(scoped, glucose_uptake, oxygen_uptake)
        hold_pool(scoped, resolved, content_mmol_per_gdcw, biomass_reaction)
        leak = scoped.reactions.get_by_id(exchange_reaction)
        leak.bounds = (leak_flux_mmol_per_gdcw_h, leak.upper_bound)
        growth_held_leaking = _growth(scoped)
    return RetentionCost(
        metabolite=resolved,
        content_mmol_per_gdcw=content_mmol_per_gdcw,
        growth_free=growth_free,
        growth_held=growth_held,
        retention_flux=flux,
        growth_exported=growth_exported,
        leak_flux=leak_flux_mmol_per_gdcw_h,
        growth_held_leaking=growth_held_leaking,
        glucose_uptake=glucose_uptake,
        oxygen_uptake=oxygen_uptake,
    )
