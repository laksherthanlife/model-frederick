"""Make the proton budget bind, then price a weak-acid load against it.

Shipped Yeast9 disposes of a cytosolic proton for 1/3 ATP by excreting pyrophosphate through
the r_4460 symport (0.61 ATP in the shipped default medium, 1/3 in every glucose-limited regime
measured here), so its H+-ATPase r_0227 carries **zero** flux at the optimum and pH costs
nothing. Closing that route makes a cytosolic proton cost one maintenance-ATP equivalent. A
weak acid then enters undissociated, dissociates in the cytosol, and the released proton is
priced against that. Read :func:`proton_price` before quoting the "1 ATP" figure: r_0227 is
DISPENSABLE here, so the price is not evidence of its 1:1 stoichiometry.

Measured on Yeast9 v9.0.2 with GLPK: one propionate load of 3 mmol/gDCW/h moves 441 reactions
across all eight compartments -- one input reaching many reactions at different depths. See
:func:`proton_load_depth`, which names every one of them rather than asserting the claim.

Aerobically this is NOT the whole story and the module does not pretend otherwise: acetate and
lactate are carbon sources, so at glucose 10 with O2 free they BUY growth (-4.43 and -6.07 ATP
per mol) while propionate, which Yeast9 cannot assimilate, costs exactly 1.00. The uncoupling
cost is only visible on its own anaerobically. Numbers in PROVENANCE.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra

from .solver import growth_or_none

__all__ = [
    "CYTOSOLIC_PH",
    "DEFAULT_PERMEABILITY",
    "PMA1_REACTION",
    "PROTON_ESCAPES",
    "ProtonBudgetReport",
    "PROTON_EXCHANGE",
    "PYROPHOSPHATE_EXCHANGE",
    "WEAK_ACIDS",
    "WeakAcid",
    "acid_reaction_ids",
    "add_weak_acid_uncoupling",
    "anaerobic_constraints",
    "constrain_proton_budget",
    "acid_load_response",
    "atp_per_acid",
    "marginal_atp_price",
    "ph_response_curve",
    "proton_budget_cost",
    "proton_escape_routes",
    "proton_load_depth",
    "proton_price",
    "set_acid_load",
    "undissociated_fraction",
    "verduyn_anaerobic_yield",
    "weak_acid_influx",
    "VERDUYN_YATP",
    "VERDUYN_YSX",
    "VERDUYN_DILUTION",
]

PROTON_EXCHANGE = "r_1832"
PYROPHOSPHATE_EXCHANGE = "r_4527"
PYROPHOSPHATE_SYMPORT = "r_4460"
PMA1_REACTION = "r_0227"
BIOMASS_REACTION = "r_2111"
GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"

_PROTON_C, _PROTON_E = "s_0794", "s_0796"

# PMA1/PMA2 with the PMP1/PMP2 proteolipids; YND1 is a Golgi apyrase, not a plasma-membrane pump.
PMA1_GPR = "(YCR024C-A and YEL017C-A and YGL008C) or (YCR024C-A and YEL017C-A and YPL036W)"
APYRASE_GENE = "YER005W"

# Measured cytosolic pH of exponentially growing S. cerevisiae at medium pH 5.0.
CYTOSOLIC_PH = 7.08

# PLACEHOLDER. mmol/gDCW/h of undissociated acid taken up per mM of it in the medium. No
# permeability coefficient for S. cerevisiae was found, so callers must supply a measured one.
DEFAULT_PERMEABILITY = 1.0

# Anaerobic medium: sterols, unsaturated fatty acid, heme a, and the only two vitamins that
# turn out to be required -- pantothenate for CoA and nicotinate for NAD; see PROVENANCE.
ANAEROBIC_LIPIDS = ("r_1757", "r_1915", "r_1994", "r_2106", "r_2134", "r_2189")
ANAEROBIC_HEME = "r_4780"
ANAEROBIC_VITAMINS = ("r_1548", "r_1967")
VITAMIN_TRACE_UPTAKE = 0.01

# Verduyn et al. 1990's measured anaerobic chemostat, the check this module is scored against.
VERDUYN_DILUTION = 0.10
VERDUYN_YSX = 0.10
VERDUYN_YATP = 16.0
GLUCOSE_G_PER_MMOL = 0.18016

# Free proton sinks found with proton_escape_routes, not guessed: each lets a cytosolic proton
# leave attached to an excreted metabolite. Value is the direction to shut.
PROTON_ESCAPES = {
    "r_4527": ("upper", "diphosphate exchange -- 2 H+ leave per PPi through the r_4460 symport, "
                        "and yeast does not excrete pyrophosphate"),
}


@dataclass(frozen=True)
class WeakAcid:
    """A medium acid that crosses the membrane undissociated and acidifies the cytosol."""

    key: str
    name: str
    pka: float
    formula: str
    source: str
    anion_c: str | None = None
    anion_e: str | None = None
    acid_c: str | None = None
    anion_exchange: str | None = None
    add_anion_efflux: bool = True
    escapes: tuple[tuple[str, str, str], ...] = ()

    @property
    def acid_e_id(self) -> str:
        return f"{self.key}_ha_e"

    @property
    def acid_c_id(self) -> str:
        return self.acid_c or f"{self.key}_ha_c"

    @property
    def anion_c_id(self) -> str:
        return self.anion_c or f"{self.key}_an_c"

    @property
    def anion_e_id(self) -> str:
        return self.anion_e or f"{self.key}_an_e"

    @property
    def anion_exchange_id(self) -> str:
        return self.anion_exchange or f"EX_{self.anion_e_id}"


_PKA_SOURCE = "Martell & Smith, Critical Stability Constants (Plenum, 1976), 25 C, zero ionic strength"

WEAK_ACIDS: dict[str, WeakAcid] = {
    # Acetate exists in Yeast9 in both compartments with a reversible carrier (r_1106, ADY2).
    "acetic": WeakAcid(
        key="acetic", name="acetic acid", pka=4.757, formula="C2H4O2", source=_PKA_SOURCE,
        anion_c="s_0362", anion_e="s_0364", anion_exchange="r_1634", add_anion_efflux=False,
    ),
    # Yeast9's "propionate" s_4250 is C3H6O2 charge 0 -- the undissociated acid, mislabelled.
    # It is NOT reused: r_4678 would turn a propionate load into an anaerobic redox sink.
    "propionic": WeakAcid(
        key="propionic", name="propionic acid", pka=4.874, formula="C3H6O2", source=_PKA_SOURCE,
    ),
    # (S)-lactate exists in both compartments but r_1207 is inward-only, so efflux is added.
    "lactic": WeakAcid(
        key="lactic", name="L-lactic acid", pka=3.86, formula="C3H6O3", source=_PKA_SOURCE,
        anion_c="s_0063", anion_e="s_0064", anion_exchange="r_1551",
        escapes=(("r_4711", "upper",
                  "ethyl-(2S)-lactate esterase consumes the cytosolic proton and the ester is "
                  "excreted, so a lactate load costs nothing until this direction is shut"),),
    ),
}


def undissociated_fraction(medium_ph: float, pka: float) -> float:
    """Henderson-Hasselbalch fraction of the acid that is undissociated at ``medium_ph``."""
    return 1.0 / (1.0 + 10.0 ** (medium_ph - pka))


def weak_acid_influx(
    total_acid_mM: float,
    medium_ph: float,
    pka: float,
    permeability: float = DEFAULT_PERMEABILITY,
) -> float:
    """Undissociated-acid uptake, mmol/gDCW/h. ``permeability`` is a PLACEHOLDER by default."""
    return permeability * total_acid_mM * undissociated_fraction(medium_ph, pka)


def _resolve(model: cobra.Model, base_id: str, compartment: str) -> str:
    """Map a bare Yeast9 metabolite id onto this model's naming (GECKO uses ``s_0362[c]``)."""
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(f"metabolite {base_id!r} not found in {model.id!r}")


def _add_reaction(model, rid, name, stoichiometry, gene_rule="", bounds=(0.0, 1000.0)):
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.subsystem = "weak acid uncoupling"
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gene_rule
    return rxn


def acid_reaction_ids(acid: WeakAcid) -> tuple[str, ...]:
    """Reaction ids ``add_weak_acid_uncoupling`` installs for ``acid``."""
    ids = [f"EX_{acid.acid_e_id}", f"{acid.key.upper()}_DIFF", f"{acid.key.upper()}_DISS"]
    if acid.add_anion_efflux:
        ids.append(f"{acid.key.upper()}_ANEFF")
    if acid.anion_exchange is None:
        ids.append(f"EX_{acid.anion_e_id}")
    return tuple(ids)


def add_weak_acid_uncoupling(
    model: cobra.Model,
    acid: str | WeakAcid = "acetic",
    block_escapes: bool = True,
) -> cobra.Model:
    """Return a copy of ``model`` that can take up ``acid`` undissociated and dissociate it.

    Four steps, none of them charging ATP: the undissociated acid enters by passive diffusion
    (Casal et al. 1996 found it is the only route in glucose-grown cells), dissociates in the
    cytosol at pH 7.08 well above every pKa here, the anion leaves, and the proton is left for
    the host model to price. Pair with :func:`constrain_proton_budget` or the proton is disposed
    of at 1/3 ATP by pyrophosphate excretion and the load is under-priced.

    Args:
        model: Host model to extend. Not modified.
        acid: Key into :data:`WEAK_ACIDS`, or a :class:`WeakAcid`.
        block_escapes: Shut this acid's verified free proton sinks (``acid.escapes``).

    Raises:
        ValueError: if this acid is already installed.
        KeyError: naming a metabolite the host model does not carry.
    """
    spec = WEAK_ACIDS[acid] if isinstance(acid, str) else acid
    if model.metabolites.has_id(spec.acid_e_id):
        raise ValueError(f"{spec.name} uncoupling is already installed on this model")

    out = model.copy()
    key = spec.key.upper()
    proton_c = _resolve(out, _PROTON_C, "c")

    new_mets = [cobra.Metabolite(spec.acid_e_id, name=f"{spec.name} (undissociated)",
                                 formula=spec.formula, charge=0, compartment="e")]
    if spec.acid_c is None:
        new_mets.append(cobra.Metabolite(spec.acid_c_id, name=f"{spec.name} (undissociated)",
                                         formula=spec.formula, charge=0, compartment="c"))
        acid_c = spec.acid_c_id
    else:
        acid_c = _resolve(out, spec.acid_c, "c")
    anion_formula = _deprotonate(spec.formula)
    if spec.anion_c is None:
        new_mets.append(cobra.Metabolite(spec.anion_c_id, name=f"{spec.name} anion",
                                         formula=anion_formula, charge=-1, compartment="c"))
        anion_c = spec.anion_c_id
    else:
        anion_c = _resolve(out, spec.anion_c, "c")
    if spec.anion_e is None:
        new_mets.append(cobra.Metabolite(spec.anion_e_id, name=f"{spec.name} anion",
                                         formula=anion_formula, charge=-1, compartment="e"))
        anion_e = spec.anion_e_id
    else:
        anion_e = _resolve(out, spec.anion_e, "e")
    out.add_metabolites(new_mets)

    # Uptake-only: the medium is held at a fixed pH, so the model does not set it.
    out.add_boundary(out.metabolites.get_by_id(spec.acid_e_id), type="exchange",
                     reaction_id=f"EX_{spec.acid_e_id}", lb=0.0, ub=0.0)
    _add_reaction(out, f"{key}_DIFF", f"{spec.name} passive diffusion",
                  {spec.acid_e_id: -1, acid_c: 1}, bounds=(-1000.0, 1000.0))
    _add_reaction(out, f"{key}_DISS", f"{spec.name} dissociation, cytosol",
                  {acid_c: -1, anion_c: 1, proton_c: 1})
    if spec.add_anion_efflux:
        _add_reaction(out, f"{key}_ANEFF", f"{spec.name} anion efflux",
                      {anion_c: -1, anion_e: 1}, bounds=(-1000.0, 1000.0))
    if spec.anion_exchange is None:
        out.add_boundary(out.metabolites.get_by_id(anion_e), type="exchange",
                         reaction_id=f"EX_{anion_e}", lb=0.0, ub=1000.0)
    if block_escapes:
        for rid, side, _ in spec.escapes:
            _shut(out, rid, side)
    return out


def _shut(model: cobra.Model, rid: str, side: str) -> bool:
    """Close one direction of ``rid``. Returns False if the model does not carry it."""
    if not model.reactions.has_id(rid):
        return False
    r = model.reactions.get_by_id(rid)
    if side == "upper":
        r.upper_bound = 0.0
    else:
        r.lower_bound = 0.0
    return True


def _deprotonate(formula: str) -> str:
    """Strip one hydrogen from a neutral acid formula, e.g. C2H4O2 -> C2H3O2."""
    import re

    def sub(match):
        n = int(match.group(1) or 1)
        return "H" if n == 2 else f"H{n - 1}"

    out, count = re.subn(r"H(\d*)", sub, formula, count=1)
    if not count:
        raise ValueError(f"formula {formula!r} has no hydrogen to lose")
    return out


def set_acid_load(model: cobra.Model, acid: str | WeakAcid, influx: float) -> cobra.Reaction:
    """Force ``influx`` mmol/gDCW/h of undissociated acid into the cell. Modifies ``model``."""
    spec = WEAK_ACIDS[acid] if isinstance(acid, str) else acid
    rxn = model.reactions.get_by_id(f"EX_{spec.acid_e_id}")
    rxn.bounds = (-abs(influx), -abs(influx))
    return rxn


@dataclass(frozen=True)
class ProtonBudgetReport:
    """What :func:`constrain_proton_budget` actually changed, and what it could not find."""

    changed: dict[str, tuple[float, float]]
    missing: tuple[str, ...]

    def summary(self) -> str:
        lines = [f"{rid}: bounds -> {b}" for rid, b in self.changed.items()]
        if self.missing:
            lines.append(f"not present in this model: {', '.join(self.missing)}")
        return "\n".join(lines)


def constrain_proton_budget(
    model: cobra.Model,
    block_pyrophosphate_efflux: bool = True,
    block_proton_uptake: bool = True,
    restrict_pma1_gpr: bool = True,
) -> ProtonBudgetReport:
    """Close the free routes that let Yeast9 dispose of a cytosolic proton below 1 ATP.

    Modifies ``model`` in place; use inside a ``with model:`` block to scope it.

    Args:
        model: Model to constrain.
        block_pyrophosphate_efflux: Shut r_4527. Yeast does not excrete pyrophosphate, and with
            r_4527 open the model exports 2 H+ per PPi through the symporter r_4460 at 1/3 ATP.
            This is the constraint that makes the proton cost what Pma1 charges.
        block_proton_uptake: Shut the uptake direction of r_1832. The medium is a proton sink at
            a fixed pH, not a reservoir the cell may draw on to power its symporters.
        restrict_pma1_gpr: Drop YND1 (a Golgi apyrase) from r_0227's gene rule so that deleting
            PMA1/PMA2 actually removes the pump.
    """
    changed: dict[str, tuple[float, float]] = {}
    missing: list[str] = []
    if block_pyrophosphate_efflux:
        for rid, (side, _reason) in PROTON_ESCAPES.items():
            if _shut(model, rid, side):
                changed[rid] = model.reactions.get_by_id(rid).bounds
            else:
                missing.append(rid)
    if block_proton_uptake:
        if model.reactions.has_id(PROTON_EXCHANGE):
            r = model.reactions.get_by_id(PROTON_EXCHANGE)
            r.lower_bound = 0.0
            changed[PROTON_EXCHANGE] = r.bounds
        else:
            missing.append(PROTON_EXCHANGE)
    if restrict_pma1_gpr:
        if model.reactions.has_id(PMA1_REACTION):
            model.reactions.get_by_id(PMA1_REACTION).gene_reaction_rule = PMA1_GPR
        else:
            missing.append(PMA1_REACTION)
    return ProtonBudgetReport(changed=changed, missing=tuple(missing))


def proton_escape_routes(
    model: cobra.Model,
    load: float = 1.0,
    tolerance: float = 1e-6,
) -> list[tuple[str, str, float]]:
    """Which boundary reactions move when a proton is dumped in the cytosol.

    This is how the free sinks in :data:`PROTON_ESCAPES` were found rather than guessed: a
    healthy model answers with H+ exchange, a leaky one answers with something else leaving.
    Does not modify ``model``.
    """
    from cobra.flux_analysis import pfba

    with model as m:
        base = pfba(m).fluxes
    with model as m:
        src = cobra.Reaction("_proton_escape_probe", lower_bound=load, upper_bound=load)
        m.add_reactions([src])
        src.add_metabolites({m.metabolites.get_by_id(_resolve(m, _PROTON_C, "c")): 1})
        loaded = pfba(m).fluxes
    delta = (loaded - base).fillna(0.0)
    rows = [(rid, model.reactions.get_by_id(rid).name, float(v)) for rid, v in delta.items()
            if abs(v) > tolerance and model.reactions.has_id(rid)
            and model.reactions.get_by_id(rid).boundary]
    rows.sort(key=lambda row: abs(row[2]), reverse=True)
    return rows


def anaerobic_constraints(
    model: cobra.Model,
    glucose_uptake: float = 10.0,
    vitamin_uptake: float = VITAMIN_TRACE_UPTAKE,
) -> cobra.Model:
    """Apply an anaerobic glucose-limited regime in place, as Verduyn et al. (1990) ran it.

    Yeast9 cannot grow with oxygen shut: heme a, NAD and coenzyme A all need O2 somewhere in
    their biosynthesis. Rather than edit the fixed biomass composition, this supplies them from
    the medium the way an anaerobic mineral medium does.
    """
    model.reactions.get_by_id(GLUCOSE_EXCHANGE).lower_bound = -abs(glucose_uptake)
    model.reactions.get_by_id(OXYGEN_EXCHANGE).bounds = (0.0, 0.0)
    for rid in ANAEROBIC_LIPIDS + (ANAEROBIC_HEME,):
        model.reactions.get_by_id(rid).lower_bound = -1000.0
    for rid in ANAEROBIC_VITAMINS:
        model.reactions.get_by_id(rid).lower_bound = -abs(vitamin_uptake)
    return model



def _growth(model: cobra.Model, doing: str) -> float:
    # slim_optimize signals failure with nan, and both obvious guards are False against nan.
    value = growth_or_none(model)
    if value is None:
        raise ValueError(
            f"the LP is infeasible while {doing}; refusing to read a price off it rather than "
            "letting nan flow into a reported number")
    return float(value)


def marginal_atp_price(model: cobra.Model) -> float:
    """Growth lost per mmol/gDCW/h of extra cytosolic ATP hydrolysis. Does not modify ``model``."""
    with model as m:
        base = _growth(m, "pricing cytosolic ATP hydrolysis")
        r = cobra.Reaction("_atp_price_probe", lower_bound=1.0, upper_bound=1.0)
        m.add_reactions([r])
        r.add_metabolites({m.metabolites.get_by_id(_resolve(m, "s_0434", "c")): -1,
                           m.metabolites.get_by_id(_resolve(m, "s_0803", "c")): -1,
                           m.metabolites.get_by_id(_resolve(m, "s_0394", "c")): 1,
                           m.metabolites.get_by_id(_resolve(m, "s_1322", "c")): 1,
                           m.metabolites.get_by_id(_resolve(m, _PROTON_C, "c")): 1})
        loaded = _growth(m, "pricing cytosolic ATP hydrolysis under load")
    return base - loaded


def proton_price(model: cobra.Model) -> float:
    """Growth cost of one cytosolic proton, in maintenance-ATP equivalents. Does not modify.

    NOT a coupling ratio. The denominator is an ATP-hydrolysis probe that itself releases one
    cytosolic proton, so once that proton dominates the probe's cost this returns 1.0 whatever
    r_0227 does: MEASURED unchanged at 1.0000 with r_0227 doctored to 2:1 and 3:1 and with it
    deleted outright. Shipped Yeast9 answers 1/3 because pyrophosphate carries the proton out.
    """
    with model as m:
        base = _growth(m, "pricing a cytosolic proton")
        r = cobra.Reaction("_proton_price_probe", lower_bound=1.0, upper_bound=1.0)
        m.add_reactions([r])
        r.add_metabolites({m.metabolites.get_by_id(_resolve(m, _PROTON_C, "c")): 1})
        loaded = _growth(m, "pricing a cytosolic proton under load")
    return (base - loaded) / marginal_atp_price(model)


def proton_budget_cost(model: cobra.Model) -> dict[str, float]:
    """What closing the budget costs, measured rather than asserted. Does not modify ``model``.

    Apply the regime to ``model`` first; the answer is regime-dependent and this returns the
    regime's numbers, not a constant.
    """
    with model as m:
        before = _growth(m, "reading the proton budget before constraining it")
        price_before = proton_price(m)
        pump_before = float(m.optimize().fluxes[PMA1_REACTION])
        constrain_proton_budget(m)
        after = _growth(m, "reading the proton budget after constraining it")
        price_after = proton_price(m)
        pump_after = float(m.optimize().fluxes[PMA1_REACTION])
    return {"growth_before": float(before), "growth_after": float(after),
            "relative_growth_cost": float((before - after) / before),
            "proton_price_before": price_before, "proton_price_after": price_after,
            "pma1_before": pump_before, "pma1_after": pump_after}


def verduyn_anaerobic_yield(
    model: cobra.Model,
    dilution: float = VERDUYN_DILUTION,
) -> dict[str, float]:
    """Anaerobic biomass yield at a fixed growth rate, against Verduyn et al. 1990's 0.10 g/g.

    Growth is PINNED and glucose minimised, which is what a chemostat does; reading a yield off
    a growth-maximising solve at an arbitrary uptake would compare a different quantity.
    Modifies nothing. Apply :func:`anaerobic_constraints` with a non-binding glucose bound and
    :func:`constrain_proton_budget` before calling.
    """
    with model as m:
        m.reactions.get_by_id(BIOMASS_REACTION).bounds = (dilution, dilution)
        m.objective = m.problem.Objective(
            m.reactions.get_by_id(GLUCOSE_EXCHANGE).flux_expression, direction="max")
        sol = m.optimize()
        qs = -float(sol.fluxes[GLUCOSE_EXCHANGE])
        ethanol = float(sol.fluxes["r_1761"])
    ysx = dilution / (qs * GLUCOSE_G_PER_MMOL)
    return {"dilution": dilution, "glucose_uptake": qs, "ysx": ysx, "ethanol": ethanol,
            "measured_ysx": VERDUYN_YSX, "relative_error": (ysx - VERDUYN_YSX) / VERDUYN_YSX}


def acid_load_response(
    model: cobra.Model,
    acid: str | WeakAcid,
    loads,
) -> list[dict[str, float]]:
    """Growth, Pma1 flux and anion export against undissociated-acid influx.

    Does not modify ``model``. Below the load at which the anion's own metabolic value
    saturates the response is not linear; above it every further mmol is a futile cycle.
    """
    spec = WEAK_ACIDS[acid] if isinstance(acid, str) else acid
    anion_ex = spec.anion_exchange_id
    rows = []
    for load in loads:
        with model as m:
            set_acid_load(m, spec, load)
            sol = m.optimize()
            rows.append({"load": float(load), "growth": float(sol.objective_value),
                         "pma1_flux": float(sol.fluxes[PMA1_REACTION]),
                         "anion_export": float(sol.fluxes[anion_ex]),
                         "proton_export": float(sol.fluxes[PROTON_EXCHANGE])})
    return rows


def atp_per_acid(model: cobra.Model, acid: str | WeakAcid, low: float, high: float) -> float:
    """ATP-equivalents charged per mol of acid over ``[low, high]``, from the growth slope.

    In the same units as :func:`proton_price`, so for an inert acid it counts PROTONS delivered
    per mol, not a coupling ratio: MEASURED 1/2/3 when the dissociation step is doctored to
    release 1/2/3 protons, and unchanged at 1.0000 with r_0227 deleted.
    Verduyn et al. (1990) measured the growth penalty as linear in the acid added, and it is
    1.00 here for an acid the cell cannot eat. NEGATIVE means the acid pays: aerobically acetate and lactate are
    carbon sources and return -4.43 and -6.07, so a bare "acid costs ATP" claim is wrong in the
    regime a bioprocess actually runs in. Uncoupling alone is the anaerobic number.
    """
    rows = acid_load_response(model, acid, [low, high])
    slope = (rows[1]["growth"] - rows[0]["growth"]) / (high - low)
    return -slope / marginal_atp_price(model)


def ph_response_curve(
    model: cobra.Model,
    acid: str | WeakAcid,
    total_acid_mM: float,
    ph_values,
    permeability: float = DEFAULT_PERMEABILITY,
) -> list[dict[str, float]]:
    """Growth and ATP bill against medium pH at fixed total acid. Does not modify ``model``.

    This is the curve the ODE pH layer prices against. Its SHAPE is Henderson-Hasselbalch and
    carries no free parameter; only its SCALE does, through the placeholder ``permeability``.

    Prefer ``atp_burden`` -- the load's cost in ATP equivalents read off the growth lost against
    the unloaded reference, which also catches an acid that pays for itself. ``atp_cost`` is
    only r_0227's flux in ONE optimal solution: it carries a growth-dependent baseline on top of
    the load (0.478 against an influx of 0.037 at pH 7), and it reads 0.000 with r_0227 deleted
    while the growth penalty is unchanged, so it is not the quantity the load costs.
    """
    spec = WEAK_ACIDS[acid] if isinstance(acid, str) else acid
    per_atp = marginal_atp_price(model)
    with model as m:
        set_acid_load(m, spec, 0.0)
        reference = _growth(m, "reading the zero-acid reference growth")
    rows = []
    for ph in ph_values:
        frac = undissociated_fraction(ph, spec.pka)
        influx = weak_acid_influx(total_acid_mM, ph, spec.pka, permeability)
        with model as m:
            set_acid_load(m, spec, influx)
            sol = m.optimize()
            ok = sol.status == "optimal"
            growth = float(sol.objective_value) if ok else float("nan")
            pma1 = float(sol.fluxes[PMA1_REACTION]) if ok else float("nan")
            proton = float(sol.fluxes[PROTON_EXCHANGE]) if ok else float("nan")
        rows.append({"ph": float(ph), "undissociated_fraction": frac, "acid_influx": influx,
                     "growth": growth, "pma1_flux": pma1, "atp_cost": pma1,
                     "proton_export": proton,
                     "relative_growth": growth / reference,
                     "atp_burden": (reference - growth) / per_atp})
    return rows


def proton_load_depth(
    model: cobra.Model,
    acid: str | WeakAcid,
    influx: float,
    tolerance: float = 1e-6,
) -> list[tuple[str, str, float]]:
    """Every reaction whose parsimonious flux moves when the acid load is applied.

    Returns ``(reaction id, name, delta)`` sorted by magnitude. Does not modify ``model``.
    """
    from cobra.flux_analysis import pfba

    spec = WEAK_ACIDS[acid] if isinstance(acid, str) else acid
    with model as m:
        base = pfba(m).fluxes
    with model as m:
        set_acid_load(m, spec, influx)
        loaded = pfba(m).fluxes
    delta = (loaded - base).reindex(loaded.index).fillna(0.0)
    moved = [(rid, model.reactions.get_by_id(rid).name, float(v))
             for rid, v in delta.items() if abs(v) > tolerance]
    moved.sort(key=lambda row: abs(row[2]), reverse=True)
    return moved


PROVENANCE = {
    "pKa acetic 4.757, propionic 4.874": _PKA_SOURCE,
    "pKa lactic 3.86": "widely tabulated value, same reference family as above",
    "cytosolic pH 7.08 +- 0.05": (
        "Orij, Urbanus, Vink, Brul, Smits et al. 2012, Genome Biol 13:R80 -- pHluorin, "
        "2088 biological replicates at medium pH 5.0"
    ),
    "Pma1 translocates protons at the expense of ATP": (
        "Malpartida & Serrano 1981, FEBS Lett 131:351-354 -- purified yeast plasma membrane "
        "ATPase reconstituted in liposomes"
    ),
    "Pma1 stoichiometry 1 H+ per ATP": (
        "Perlin, San Francisco, Slayman & Rosen 1986, Arch Biochem Biophys 248:53 -- coupling "
        "ratio of 1 H+ per ATP for the fungal plasma-membrane P-type pump. Venema & Palmgren "
        "1995, J Biol Chem 270:19659 show the yeast ratio is metabolically modulated and that "
        "the pump from glucose-STARVED cells is essentially uncoupled, so 1:1 is the "
        "glucose-fed value and every regime in this module is glucose-fed. Yeast9's r_0227 "
        "(s_0434 + s_0803 --> s_0394 + s_0796 + s_1322) already encodes it: the hydrolysis "
        "proton is deposited outside, so no new pump reaction is added. VERIFIED by reading "
        "the reaction, and check_mass_balance() on it is empty."
    ),
    "undissociated acid is the only route in": (
        "Casal, Cardoso & Leao 1996, Microbiology 142:1385 -- in glucose-grown S. cerevisiae "
        "acetic acid enters only by simple diffusion of the undissociated form"
    ),
    "weak-acid uncoupling costs ATP, linearly in acid added": (
        "Verduyn, Postma, Scheffers & van Dijken 1990, 'Energetics of Saccharomyces cerevisiae "
        "in anaerobic glucose-limited chemostat cultures', J Gen Microbiol 136:405-412 -- YATP "
        "about 16 g biomass/mol ATP, lowered by acetic or propionic acid, with a linear "
        "correlation between the energy required to compensate proton import and the acid "
        "added. This module reproduces the linearity exactly (Pma1 rises 1.000 mmol/gDCW/h per "
        "mmol of propionate) and under-states the magnitude: its marginal ATP yield anaerobic "
        "at glucose 10 is 10.32 g biomass/mol ATP against the measured 16, so a growth penalty "
        "computed here is 0.65x the one Verduyn's YATP implies."
    ),
    "anaerobic biomass yield 0.10 g/g glucose at D = 0.10 /h": (
        "Verduyn, Postma, Scheffers & van Dijken 1990, J Gen Microbiol 136:395-403 "
        "(PMID 1975265). Reproduced by verduyn_anaerobic_yield at 0.1068 g/g, +6.8%."
    ),
    "DEFAULT_PERMEABILITY = 1.0": (
        "PLACEHOLDER. No permeability coefficient for S. cerevisiae was found, and a bilayer "
        "coefficient in cm/s would still need a membrane-area-per-gDCW conversion that is a "
        "second unmeasured parameter here. Only the SCALE of the pH curve depends on it; its "
        "shape is Henderson-Hasselbalch and has no free parameter."
    ),
    "growth cost of closing the budget, by regime": (
        "MEASURED on Yeast9 v9.0.2, GLPK: shipped default -0.726% (0.085844 -> 0.085221 /h), "
        "glucose<=10 with O2 free -0.726% (0.887685 -> 0.881241), REFERENCE_AEROBIC_BATCH "
        "glucose 11.1 / O2 3.7 -1.669% (0.345984 -> 0.340208), anaerobic glucose 10 -1.527% "
        "(0.202264 -> 0.199175). r_0227 goes from zero flux to 0.20/2.08/0.80/0.44 in the same "
        "four, and r_4527 from 0.20/2.10/0.82/0.45 to zero."
    ),
    "one pH input moves 441 reactions": (
        "MEASURED by proton_load_depth on the anaerobic regime with a 3 mmol/gDCW/h propionate "
        "load: 441 reactions move, touching all eight compartments (c, e, m, mm, er, erm, gm, "
        "ce). Depths, in order: the four installed reactions, then r_1832 and r_0227, then "
        "ethanol and CO2 export, pyruvate decarboxylase and alcohol dehydrogenase, then "
        "pyrophosphatase r_0569, polyphosphate hydrolysis r_4333, glutamine synthetase r_0476 "
        "and ammonium exchange r_1654, then 36 glycerophospholipid and 274 ER-membrane "
        "reactions rescaled with biomass."
    ),
    "anaerobic supplements are the ones the model actually needs": (
        "MEASURED by dropping each in turn: pantothenate (r_1548) and nicotinate (r_1967) each kill growth outright, and every other candidate costs under 0.3%. Yeast9 cannot grow with O2 shut because heme a, NAD and CoA all need oxygen somewhere, so an anaerobic mineral medium supplies them."
    ),
    "VITAMIN_TRACE_UPTAKE = 0.01": (
        "Not binding: pantothenate and nicotinate draw 4e-5 and 1.2e-3 mmol/gDCW/h, so anything "
        "from 0.005 upwards gives the same answer and the regime stays glucose-limited"
    ),
    "ATP charged per mol of acid, by acid and regime": (
        "MEASURED by atp_per_acid over loads 2 -> 5 mmol/gDCW/h. Budget closed: propionic "
        "1.0000 aerobic and 1.0000 anaerobic, lactic 1.0000 anaerobic, acetic 1.0348 "
        "anaerobic. Budget shipped: propionic 0.3333 aerobic -- the threefold under-charge "
        "this module exists to remove. Aerobically acetate and lactate are SUBSTRATES and "
        "return -4.4347 and -6.0668, i.e. they buy growth; only propionate, which Yeast9 "
        "cannot assimilate, is a pure uncoupler in both regimes."
    ),
    "the price is NOT evidence of Pma1's stoichiometry": (
        "MEASURED, and it refutes the obvious reading of the 1.0000. With the budget closed at "
        "glucose 10 / O2 free, deleting r_0227 outright changes growth by 3e-15, leaves "
        "proton_price at 1.0000 and atp_per_acid(propionic) at 1.0000, and makes ph_response_"
        "curve's atp_cost read 0.000 for an unchanged growth penalty. proton_price also returns "
        "1.0000 with r_0227 doctored to 2 and 3 H+ per ATP, and atp_per_acid returns 1/2/3 when "
        "the DISSOCIATION step is doctored to release 1/2/3 protons. So these numbers count "
        "protons in units of one maintenance-ATP equivalent; the pump is a correlate at the "
        "optimum, not the payer. The 3x change on closing the budget IS real: the growth cost "
        "per mol of propionate goes 0.001547 -> 0.004609 /h per mmol/gDCW/h, a factor 2.979."
    ),
    "anion efflux carries no ATP cost": (
        "ASSERTED. Yeast9 already models acetate anion transport as a free carrier (r_1106, "
        "ADY2); the real propionate exporter Pdr12 (Piper et al. 1998, EMBO J 17:4257) is "
        "ATP-driven with no measured stoichiometry, so the cost reported here is a lower bound."
    ),
}
