"""Make biomass composition a parameter, so a stress state can move it.

Yeast9 pins the six carbohydrate components of r_4048 as constants. Nothing can accumulate
trehalose, deplete glycogen or thicken the wall, and deleting any one of them reads as an exact
-100% lethal -- a pseudoreaction artefact that makes an ablation study meaningless. These helpers
rewrite the pseudoreaction, and the pool pseudoreactions above it, from measured mass fractions
while conserving the gram the biomass reaction is defined to mean.

RENORMALISATION, exactly. Mass is read from each pseudoreaction's own elemental imbalance:
``sum(coefficient * formula_weight)``, which is the mg/gDCW it absorbs. Yeast9's eight pools sum
to 956.9167 mg, not the 1000 the reaction is named for, so what is conserved here is the SHIPPED
total rather than an assumed gram -- normalising to 1000 would move every coefficient in the model
for no measured reason. Every component and pool the caller does not name is then multiplied by one
common factor ``(conserved - demanded) / (conserved - shipped_named)``, so unnamed quantities keep
their shipped proportions to each other. ``renormalise="pool"`` conserves the carbohydrate pool
mass (nothing outside r_4048 moves); ``renormalise="whole_cell"`` conserves the total, so the
carbohydrate fraction itself may rise or fall. :func:`set_pool_masses` does the same one level up.
The proportional closure is a modelling convention, not a measurement -- see MODELLING_CLOSURES.

No default stressed composition is provided. MEASURED_STATES carries the states that have a
citation; anything else is the caller's to supply and sweep.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cobra
from cobra.core.formula import elements_and_molecular_weights as _ATOMIC

__all__ = [
    "CARBOHYDRATE_COMPONENTS",
    "CARBOHYDRATE_REACTION",
    "COMPONENT_SYNTHASES",
    "MEASURED_STATES",
    "MODELLING_CLOSURES",
    "POOL_REACTIONS",
    "SLIME_COUPLED_POOLS",
    "STORAGE_COMPONENTS",
    "TREHALASE_REACTION",
    "TREHALOSE_CYCLE_RATE",
    "WALL_COMPONENTS",
    "MeasuredComposition",
    "apply_measured",
    "component_masses",
    "masses_for_protein_ratios",
    "pin_trehalose_cycle",
    "pool_masses",
    "reaction_mass",
    "scale_components",
    "set_component_masses",
    "set_pool_masses",
    "total_biomass_mass",
]

CARBOHYDRATE_REACTION = "r_4048"
BIOMASS_REACTION = "r_4041"
GROWTH_REACTION = "r_2111"

POOL_REACTIONS = {
    "protein": "r_4047",
    "carbohydrate": "r_4048",
    "RNA": "r_4049",
    "DNA": "r_4050",
    "lipid_backbone": "r_4063",
    "lipid_chain": "r_4065",
    "cofactor": "r_4598",
    "ion": "r_4599",
}

# The formula-less pool metabolite each pseudoreaction produces; it carries no mass, so it is the
# one species a rescaling must leave alone.
POOL_PRODUCTS = {
    "r_4047": "s_3717",
    "r_4048": "s_3718",
    "r_4049": "s_3719",
    "r_4050": "s_3720",
    "r_4063": "s_3746",
    "r_4065": "s_3747",
    "r_4598": "s_4205",
    "r_4599": "s_4206",
}

# One quantity, not two: the SLIME reactions co-produce backbone and chain in a fixed ratio,
# r_2108 merges them 1:1, and both overflow exchanges (r_4062, r_4064) ship clamped to (0, 0).
SLIME_COUPLED_POOLS = ("lipid_backbone", "lipid_chain")

CARBOHYDRATE_COMPONENTS = {
    "beta_1_3_glucan": "s_0001",
    "beta_1_6_glucan": "s_0004",
    "chitin": "s_0509",
    "glycogen": "s_0773",
    "mannan": "s_1107",
    "trehalose": "s_1520",
}
WALL_COMPONENTS = ("beta_1_3_glucan", "beta_1_6_glucan", "chitin", "mannan")
STORAGE_COMPONENTS = ("glycogen", "trehalose")

COMPONENT_SYNTHASES = {
    "beta_1_3_glucan": "r_0005",
    "beta_1_6_glucan": "r_0006",
    "chitin": "r_0272",
    "glycogen": "r_0510",
    "mannan": "r_0362",
    "trehalose": "r_1051",
}
"""The gene-bearing reaction that makes each component, read from Yeast9 rather than assumed.

Mannan is the odd one: Yeast9 has NO mannan synthase. ``s_1107`` is filled only by transport
``r_1932`` from the ER species that ``r_0362`` (protein O-mannosyltransferase, PMT1-PMT4) makes,
so an ablation of "mannan synthesis" is an ablation of protein mannosylation.
"""

TREHALASE_REACTION = "r_0194"
TREHALOSE_PHOSPHATASE = "r_1051"
TREHALOSE_EXCHANGE = "r_1650"

TREHALOSE_CYCLE_RATE = None
"""REFUSED. Hottiger 1987 shows the cycle runs -- trehalose "turned over rapidly", synthase up 6x
and neutral trehalase up 3x at once -- but reports enzyme activities, not an in vivo flux, and no
measured cycling rate in mmol/gDCW/h was found. Sweep it; do not default it.
"""

MODELLING_CLOSURES = (
    "ASSERTED: components and pools the caller does not name are scaled by one common factor, "
    "i.e. they keep their shipped proportions to each other. No measurement fixes this closure.",
    "ASSERTED: the growth-associated maintenance term (55.3 ATP in r_4041) is left untouched when "
    "pool masses are rescaled. Yeast9 fits it to chemostat data and no measurement resolves it "
    "per macromolecule, so refusing to rescale it is the smaller error than inventing a split.",
    "ASSERTED: what is conserved is Yeast9's own 956.9167 mg total, not 1000 mg. Renormalising to "
    "the nominal gram would move every coefficient in the model without a measurement asking it "
    "to; the shipped total is at least the number the vendored GAM was fitted against.",
)

# --- raw measurements -------------------------------------------------------------------------

_GLUCOSYL_MW = 162.1406  # anhydroglucose residue C6H10O5, the unit Yeast9 uses for the polymers
_TREHALOSE_MW = 342.29648  # C12H22O11; two glucose released per molecule assayed

# Silljé et al. 1999, J Bacteriol 181:396-400 (PMID 9882651, PMC93391). Galactose-limited
# chemostat, CEN-PK113-7D; contents from Fig. 1, dry weight per cell from Fig. 2C.
_SILLJE1999 = {
    "fast": {"dilution_rate": 0.20, "pg_per_cell": 10.6, "trehalose_fmol_glc": 0.0,
             "glycogen_fmol_glc": 2.0},
    "slow": {"dilution_rate": 0.033, "pg_per_cell": 17.4, "trehalose_fmol_glc": 6.2,
             "glycogen_fmol_glc": 11.0},
}

# Ram et al. 1998, J Bacteriol 180:1418-1424 (PMID 9515908, PMC107039) Table 2: mg sugar per g dry
# weight of WHOLE CELLS. GlcN reports chitin, Man reports mannan, Glc lumps both beta-glucans.
_RAM1998 = {
    "FY834": {"GlcN": 2.8, "Glc": 76.6, "Man": 120.7},
    "gas1": {"GlcN": 4.7, "Glc": 65.3, "Man": 167.5},
    "fks1": {"GlcN": 5.3, "Glc": 38.5, "Man": 184.9},
}

# Hottiger, Schmutz & Wiemken 1987, J Bacteriol 169:5518-5522 (PMID 2960663, PMC213980), abstract:
# "the trehalose content of the cells increased from 0.01 to 1 g/g of protein within 1 h".
_HOTTIGER1987 = {"exponential_27C": 0.01, "heat_shock_40C": 1.0}


def _fmol_glucose_to_mass(fmol_glc: float, pg_per_cell: float, mw_per_glucose: float) -> float:
    """mg per gDCW from a per-cell glucose-equivalent assay. Units cancel to mg/g exactly."""
    return fmol_glc * mw_per_glucose / pg_per_cell


def _sillje_masses(phase: str) -> dict[str, float]:
    d = _SILLJE1999[phase]
    return {
        "trehalose": _fmol_glucose_to_mass(
            d["trehalose_fmol_glc"], d["pg_per_cell"], _TREHALOSE_MW / 2.0),
        "glycogen": _fmol_glucose_to_mass(
            d["glycogen_fmol_glc"], d["pg_per_cell"], _GLUCOSYL_MW),
    }


def _ram_factors(mutant: str) -> dict[str, float]:
    wt, mut = _RAM1998["FY834"], _RAM1998[mutant]
    glucan = mut["Glc"] / wt["Glc"]
    return {
        "chitin": mut["GlcN"] / wt["GlcN"],
        "mannan": mut["Man"] / wt["Man"],
        "beta_1_3_glucan": glucan,
        "beta_1_6_glucan": glucan,
    }


@dataclass(frozen=True)
class MeasuredComposition:
    """One cited composition state. Exactly one of ``masses``, ``factors`` or ``protein_ratios``."""

    name: str
    regime: str
    source: str
    derivation: str
    masses: dict[str, float] | None = None
    factors: dict[str, float] | None = None
    protein_ratios: dict[str, float] | None = None
    asserted: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        given = [f for f in (self.masses, self.factors, self.protein_ratios) if f is not None]
        if len(given) != 1:
            raise ValueError(
                f"{self.name}: give exactly one of masses, factors or protein_ratios")


MEASURED_STATES: dict[str, MeasuredComposition] = {
    "sillje1999_fast": MeasuredComposition(
        name="sillje1999_fast",
        regime="galactose-limited chemostat, D = 0.20 /h, CEN-PK113-7D",
        source="Silljé et al. 1999, J Bacteriol 181:396-400, Fig. 1 and Fig. 2C "
               "(PMID 9882651, PMC93391)",
        derivation="trehalose 0 and glycogen 2 fmol glucose/cell at 10.6 pg dry weight/cell",
        masses=_sillje_masses("fast"),
    ),
    "sillje1999_slow": MeasuredComposition(
        name="sillje1999_slow",
        regime="galactose-limited chemostat, D = 0.033 /h, CEN-PK113-7D",
        source="Silljé et al. 1999, J Bacteriol 181:396-400, Fig. 1 and Fig. 2C "
               "(PMID 9882651, PMC93391)",
        derivation="trehalose 6.2 and glycogen 11 fmol glucose/cell at 17.4 pg dry weight/cell",
        masses=_sillje_masses("slow"),
    ),
    "sillje1999_scu10": MeasuredComposition(
        name="sillje1999_scu10",
        regime="SCU10, tps1 gsy1 gsy2 triple deletant, galactose batch and chemostat "
               "D 0.033-0.20 /h",
        source="Silljé et al. 1999, J Bacteriol 181:396-400 (PMID 9882651, PMC93391)",
        derivation="trehalose and glycogen both undetectable at every dilution rate; SCU10 grew "
                   "with doubling times similar to the wild type under all conditions tested",
        masses={"trehalose": 0.0, "glycogen": 0.0},
    ),
    "cwi_gas1_ram1998": MeasuredComposition(
        name="cwi_gas1_ram1998",
        regime="gas1 deletant vs FY834, early log, whole-cell hydrolysate",
        source="Ram et al. 1998, J Bacteriol 180:1418-1424, Table 2 (PMID 9515908, PMC107039)",
        derivation="GlcN 4.7/2.8, Man 167.5/120.7, Glc 65.3/76.6 mg per g dry weight",
        factors=_ram_factors("gas1"),
        asserted=("the assay reports one Glc figure for both beta-glucans, so the same factor is "
                  "applied to (1->3)- and (1->6)-beta-glucan",),
    ),
    "cwi_fks1_ram1998": MeasuredComposition(
        name="cwi_fks1_ram1998",
        regime="fks1 deletant vs FY834, early log, whole-cell hydrolysate",
        source="Ram et al. 1998, J Bacteriol 180:1418-1424, Table 2 (PMID 9515908, PMC107039)",
        derivation="GlcN 5.3/2.8, Man 184.9/120.7, Glc 38.5/76.6 mg per g dry weight",
        factors=_ram_factors("fks1"),
        asserted=("the assay reports one Glc figure for both beta-glucans, so the same factor is "
                  "applied to (1->3)- and (1->6)-beta-glucan",),
    ),
    "hottiger1987_27C": MeasuredComposition(
        name="hottiger1987_27C",
        regime="exponential glucose culture at 27 C, immediately before the shift",
        source="Hottiger, Schmutz & Wiemken 1987, J Bacteriol 169:5518-5522, abstract "
               "(PMID 2960663, PMC213980)",
        derivation="trehalose 0.01 g per g of protein",
        protein_ratios={"trehalose": _HOTTIGER1987["exponential_27C"]},
    ),
    "hottiger1987_heat_shock": MeasuredComposition(
        name="hottiger1987_heat_shock",
        regime="same culture 1 h after a 27 -> 40 C shift",
        source="Hottiger, Schmutz & Wiemken 1987, J Bacteriol 169:5518-5522, abstract "
               "(PMID 2960663, PMC213980)",
        derivation="trehalose 1 g per g of protein",
        protein_ratios={"trehalose": _HOTTIGER1987["heat_shock_40C"]},
        asserted=("the ratio is to protein, so it fixes trehalose only once protein is fixed; "
                  "under renormalise='whole_cell' both move and the state is the fixed point",),
    ),
}

# --- reading the shipped composition ------------------------------------------------------------


def _resolve(model: cobra.Model, base_id: str) -> str:
    """Map a bare Yeast9 metabolite id onto whatever suffix convention this model uses."""
    if model.metabolites.has_id(base_id):
        return base_id
    for met in model.metabolites:
        if met.id.startswith(f"{base_id}["):
            return met.id
    raise KeyError(f"metabolite {base_id!r} not found in {model.id!r}")


def reaction_mass(model: cobra.Model, reaction_id: str) -> float:
    """mg per gDCW absorbed by a pool pseudoreaction, from its elemental imbalance."""
    imbalance = model.reactions.get_by_id(reaction_id).check_mass_balance()
    mass = 0.0
    for element, amount in imbalance.items():
        if element == "charge":
            continue
        if element not in _ATOMIC:
            raise ValueError(
                f"{reaction_id} is unbalanced in {element!r}, which has no atomic weight; "
                "its mass cannot be accounted for"
            )
        mass -= _ATOMIC[element] * amount
    return mass


def pool_masses(model: cobra.Model) -> dict[str, float]:
    """mg per gDCW held by each biomass pool present in ``model``."""
    return {
        name: reaction_mass(model, rid)
        for name, rid in POOL_REACTIONS.items()
        if model.reactions.has_id(rid)
    }


def total_biomass_mass(model: cobra.Model) -> float:
    """mg per unit biomass flux, summed over the pools. Yeast9 ships 956.917, not 1000."""
    return sum(pool_masses(model).values())


def component_masses(model: cobra.Model) -> dict[str, float]:
    """mg per gDCW of each carbohydrate component of r_4048; zero once it has been driven out."""
    rxn = model.reactions.get_by_id(CARBOHYDRATE_REACTION)
    stoichiometry = rxn.metabolites
    out = {}
    for name, base_id in CARBOHYDRATE_COMPONENTS.items():
        met = model.metabolites.get_by_id(_resolve(model, base_id))
        coefficient = stoichiometry.get(met, 0.0)
        out[name] = -coefficient * met.formula_weight if coefficient else 0.0
    return out


# --- writing a new composition ------------------------------------------------------------------


def _set_component_mass(model: cobra.Model, met_id: str, mass_mg: float) -> None:
    """cobrapy drops a metabolite once its coefficient reaches zero, so read it back defensively."""
    rxn = model.reactions.get_by_id(CARBOHYDRATE_REACTION)
    met = model.metabolites.get_by_id(met_id)
    target = -mass_mg / met.formula_weight
    rxn.add_metabolites({met: target - rxn.metabolites.get(met, 0.0)}, combine=True)


def _scale_pool(model: cobra.Model, reaction_id: str, factor: float) -> None:
    rxn = model.reactions.get_by_id(reaction_id)
    product = POOL_PRODUCTS[reaction_id]
    deltas = {
        met: coef * (factor - 1.0)
        for met, coef in rxn.metabolites.items()
        if met.id.split("[")[0] != product
    }
    rxn.add_metabolites(deltas, combine=True)


def _renormalisation_factor(conserved: float, shipped_named: float, demanded: float,
                            what: str, mode: str) -> float:
    """The one factor every unnamed quantity is multiplied by. Refuses rather than clamps."""
    free_now = conserved - shipped_named
    if free_now <= 0:
        raise ValueError(
            f"the named {what} already account for the whole conserved mass "
            f"({conserved:.4f} mg/gDCW); nothing is left to renormalise against"
        )
    if demanded >= conserved:
        raise ValueError(
            f"requested {demanded:.4f} mg/gDCW of named {what} but only {conserved:.4f} "
            f"mg/gDCW is conserved under renormalise={mode!r}"
        )
    return (conserved - demanded) / free_now


def set_component_masses(
    model: cobra.Model,
    masses: dict[str, float],
    renormalise: str = "whole_cell",
) -> cobra.Model:
    """Return a copy of ``model`` whose carbohydrate components hold ``masses`` mg/gDCW.

    Args:
        model: Host model. Not modified.
        masses: mg per gDCW keyed by :data:`CARBOHYDRATE_COMPONENTS` name. Components left out
            are rescaled, not held.
        renormalise: ``"pool"`` conserves the carbohydrate pool mass, so only r_4048 changes;
            ``"whole_cell"`` conserves the total biomass mass, so the other pools are rescaled
            and the carbohydrate fraction is free to move.

    Raises:
        KeyError: naming a component this module does not know.
        ValueError: for a negative mass, an unknown mode, or a request that leaves no mass for
            the components that were not named -- which is a refusal, not a clamp.
    """
    if renormalise not in ("pool", "whole_cell"):
        raise ValueError(f"renormalise must be 'pool' or 'whole_cell', got {renormalise!r}")
    unknown = set(masses) - set(CARBOHYDRATE_COMPONENTS)
    if unknown:
        raise KeyError(f"unknown carbohydrate components {sorted(unknown)}")
    negative = sorted(k for k, v in masses.items() if v < 0)
    if negative:
        raise ValueError(f"negative mass requested for {negative}")

    out = model.copy()
    shipped = component_masses(out)
    if renormalise == "pool":
        conserved = sum(shipped.values())
        scaled_pools: tuple[str, ...] = ()
    else:
        conserved = total_biomass_mass(out)
        scaled_pools = tuple(
            rid for name, rid in POOL_REACTIONS.items()
            if name != "carbohydrate" and out.reactions.has_id(rid)
        )
    factor = _renormalisation_factor(
        conserved, sum(shipped[k] for k in masses), sum(masses.values()),
        "components", renormalise)

    for name, base_id in CARBOHYDRATE_COMPONENTS.items():
        target = masses[name] if name in masses else shipped[name] * factor
        _set_component_mass(out, _resolve(out, base_id), target)
    for rid in scaled_pools:
        _scale_pool(out, rid, factor)
    return out


def scale_components(
    model: cobra.Model,
    factors: dict[str, float],
    renormalise: str = "whole_cell",
) -> cobra.Model:
    """Return a copy of ``model`` with named components multiplied by ``factors``."""
    shipped = component_masses(model)
    unknown = set(factors) - set(shipped)
    if unknown:
        raise KeyError(f"unknown carbohydrate components {sorted(unknown)}")
    return set_component_masses(
        model, {k: shipped[k] * f for k, f in factors.items()}, renormalise=renormalise
    )


def _refuse_broken_slime_coupling(shipped: dict[str, float],
                                  targets: dict[str, float]) -> None:
    """Refuse to move the two SLIME lipid pools apart; unequal factors kill growth silently."""
    present = [name for name in SLIME_COUPLED_POOLS if name in shipped]
    if len(present) < 2:
        return
    first, second = (targets[name] / shipped[name] for name in present)
    if abs(first - second) > 1e-9 * max(abs(first), abs(second), 1.0):
        raise ValueError(
            f"{present[0]} would scale by {first:.6g} and {present[1]} by {second:.6g}. They are "
            "one quantity in yeast-GEM's SLIME formalism -- r_2108 merges them 1:1 and the "
            "overflow exchanges r_4062/r_4064 are clamped to (0, 0) -- so an unequal move makes "
            "the model infeasible and it reports growth 0. Name both in their shipped ratio, "
            "or neither."
        )


def set_pool_masses(model: cobra.Model, masses: dict[str, float]) -> cobra.Model:
    """Return a copy of ``model`` whose named biomass POOLS hold ``masses`` mg/gDCW.

    The same closure one level up: protein, lipid or RNA becomes a parameter, every pool the
    caller does not name keeps its shipped proportion, and the total is conserved. Within a
    rescaled pool the components keep their proportions too.

    There is deliberately no measured protein or lipid state in :data:`MEASURED_STATES`. This is
    the swept parameter the brief asks for and it has no default.

    Raises:
        KeyError: naming a pool this module does not know, or one this model does not carry.
        ValueError: for a negative mass, a pool that ships at zero mass (nothing to scale), a
            request that leaves no mass for the pools that were not named, or one that would
            move the two :data:`SLIME_COUPLED_POOLS` by different factors.
    """
    unknown = set(masses) - set(POOL_REACTIONS)
    if unknown:
        raise KeyError(f"unknown biomass pools {sorted(unknown)}")
    negative = sorted(k for k, v in masses.items() if v < 0)
    if negative:
        raise ValueError(f"negative mass requested for {negative}")

    out = model.copy()
    shipped = pool_masses(out)
    absent = set(masses) - set(shipped)
    if absent:
        raise KeyError(f"{out.id!r} does not carry the pools {sorted(absent)}")
    empty = sorted(k for k in masses if shipped[k] == 0.0)
    if empty:
        raise ValueError(f"pools {empty} ship at zero mass; there are no proportions to scale")

    total = sum(shipped.values())
    factor = _renormalisation_factor(
        total, sum(shipped[k] for k in masses), sum(masses.values()), "pools", "whole_cell")
    targets = {
        name: masses[name] if name in masses else shipped[name] * factor
        for name in POOL_REACTIONS if name in shipped
    }
    _refuse_broken_slime_coupling(shipped, targets)
    for name, rid in POOL_REACTIONS.items():
        if name in shipped:
            _scale_pool(out, rid, targets[name] / shipped[name])
    return out


def masses_for_protein_ratios(
    model: cobra.Model,
    ratios: dict[str, float],
    renormalise: str = "whole_cell",
) -> dict[str, float]:
    """mg/gDCW of each component that satisfies component:protein = ``ratio`` AFTER renormalising.

    Hottiger 1987 reports trehalose per gram of protein, and under ``"whole_cell"`` protein itself
    moves, so the answer is a fixed point rather than ratio times shipped protein. With
    ``sum_x = R*P*T / (T - shipped_named + R*P)`` for R the summed ratio, P protein and T the
    conserved total. Under ``"pool"`` protein does not move and the answer is ``ratio * protein``.
    """
    if renormalise not in ("pool", "whole_cell"):
        raise ValueError(f"renormalise must be 'pool' or 'whole_cell', got {renormalise!r}")
    unknown = set(ratios) - set(CARBOHYDRATE_COMPONENTS)
    if unknown:
        raise KeyError(f"unknown carbohydrate components {sorted(unknown)}")
    negative = sorted(k for k, v in ratios.items() if v < 0)
    if negative:
        raise ValueError(f"ratios must be non-negative, got negatives for {negative}")

    protein = pool_masses(model)["protein"]
    if renormalise == "pool":
        return {k: r * protein for k, r in ratios.items()}
    shipped = component_masses(model)
    total = total_biomass_mass(model)
    summed = sum(ratios.values())
    free = total - sum(shipped[k] for k in ratios)
    scaled_protein = protein * total / (free + summed * protein)
    return {k: r * scaled_protein for k, r in ratios.items()}


def apply_measured(
    model: cobra.Model,
    state: str,
    renormalise: str = "whole_cell",
) -> cobra.Model:
    """Return a copy of ``model`` holding the cited composition state ``state``."""
    if state not in MEASURED_STATES:
        raise KeyError(f"no measured state {state!r}; have {sorted(MEASURED_STATES)}")
    record = MEASURED_STATES[state]
    if record.masses is not None:
        return set_component_masses(model, record.masses, renormalise=renormalise)
    if record.factors is not None:
        return scale_components(model, record.factors, renormalise=renormalise)
    return set_component_masses(
        model,
        masses_for_protein_ratios(model, record.protein_ratios, renormalise=renormalise),
        renormalise=renormalise,
    )


# --- the cost composition alone cannot express ---------------------------------------------------


def pin_trehalose_cycle(model: cobra.Model, rate: float) -> cobra.Model:
    """Return a copy of ``model`` forced to hydrolyse ``rate`` mmol/gDCW/h of trehalose.

    Composition alone prices trehalose accumulation as a small GAIN -- it is the cheapest
    carbohydrate per gram -- so the cost Hottiger 1987 actually observed is the TURNOVER, synthase
    and trehalase running at once. Forcing trehalase makes synthesis follow by mass balance, and
    the GEM prices the ATP. It adds no chemistry: both reactions are Yeast9's own.

    ``rate`` is required and has no default; see :data:`TREHALOSE_CYCLE_RATE`.

    Raises:
        KeyError: if this model does not carry both reactions under their Yeast9 ids -- GECKO
            splits ``r_0194`` into ``r_0194No1`` and this refuses rather than guessing.
        ValueError: for a negative rate, or if trehalose can be taken up, which would let the
            forced flux be an import instead of a cycle.
    """
    if rate < 0:
        raise ValueError(f"rate must be non-negative, got {rate}")
    for rid in (TREHALASE_REACTION, TREHALOSE_PHOSPHATASE):
        if not model.reactions.has_id(rid):
            raise KeyError(f"{model.id!r} has no reaction {rid!r}; the cycle cannot be pinned")
    if model.reactions.has_id(TREHALOSE_EXCHANGE):
        lower = model.reactions.get_by_id(TREHALOSE_EXCHANGE).lower_bound
        if lower < 0:
            raise ValueError(
                f"{TREHALOSE_EXCHANGE} allows uptake (lower bound {lower}); a forced trehalase "
                "flux would be fed by the medium and would not be a cycle"
            )
    out = model.copy()
    out.reactions.get_by_id(TREHALASE_REACTION).lower_bound = rate
    return out
