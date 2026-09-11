"""Anaerobic and hypoxic mode for Yeast9, and the diagnosis that says what anoxia costs.

Vendored yeast-GEM v9.0.2 cannot grow without oxygen: at ``r_1714 >= -10`` with
``r_1992 = 0`` the LP optimum is exactly 0, and opening ergosterol uptake alone leaves it
at 0. That is not a medium bug, it is five separate requirements, found by zeroing each of
the 74 reactions that can consume O2 and keeping the 19 whose loss is lethal:

    sterol ring         ERG1, ERG11, ERG25    6 reactions   -> ergosterol in the medium
    C16:1 / C18:1       OLE1                  2 reactions   -> Tween 80 in the medium
    NAD de novo         BNA1, BNA2, BNA4      3 reactions   -> nicotinic acid in the medium
    CoA via beta-Ala    FMS1                  1 reaction    -> pantothenate in the medium
    heme and heme a     HEM13, HEM14, COX15   3 reactions   -> NOT in any anaerobic medium
    O2 carriers                               4 reactions   -> not requirements of their own

The first four are supplied by the medium the measurement was made in: Verduyn 1990's
mineral medium carries nicotinic acid and Ca-pantothenate as vitamins and was supplemented
with ergosterol and Tween 80. The fifth cannot be. Anaerobic yeast does not make heme --
both porphyrin steps and heme A synthesis consume O2 in this model, and Snoek & Steensma
2007 (PMID 17192845) review the same requirement in the cell -- and it contains no heme a,
because it assembles no cytochrome oxidase. So heme a leaves the anaerobic biomass rather
than arriving in the feed. Importing heme a instead would fabricate a medium component;
that route is available here only through an argument that names itself.

WHAT IS ASSERTED HERE, so it is never mistaken for a measurement: taking heme a out of the
cofactor pseudoreaction is a BIOMASS COMPOSITION CHANGE, not a medium component, and it is
the one thing in this module that is a modelling judgement rather than a bound. The
coefficient removed is 1e-06 mmol/gDCW; the two routes agree to 13 decimal places, which is
what :func:`anaerobic_constraints`'s ``import_heme_a`` exists to show.

**Every number below was measured on this model, not asserted.** Anaerobic growth at
glucose <= 10 is 0.199918 /h -- with the IAH1 esterases held to hydrolysis, a restriction
that is not free and costs 1.16%; released, the optimum is 0.202264 /h and 10.5% of the
glucose carbon leaves as diethyl succinate. The maximal-yield solution at mu = 0.10 /h needs 5.177
mmol/gDCW/h of glucose, a biomass yield of 0.1072 g/g against Verduyn 1990's measured
0.10 g/g at D = 0.10 (+7.2%), and reaching that paper's mu_max of 0.31 /h needs 15.314
mmol/gDCW/h at an ethanol yield of 0.4444 g/g, 86.9% of the Gay-Lussac maximum.

TWO MEASURED DEPARTURES, both recorded rather than patched:

1. The anaerobic optimum secretes NO GLYCEROL -- FVA pins ``r_1808`` to [0, 7e-13] at the
   optimum, though it can reach 10.53 mmol/gDCW/h when growth is not maximised. Every
   anaerobic yeast culture makes glycerol, because it is the redox sink for the NADH that
   biosynthesis leaves over; Verduyn 1990 measures it in the very culture this module is
   checked against. This module does NOT force one, because no glycerol yield from that
   paper is encoded here and inventing one would be worse. :func:`anaerobic_glycerol_price`
   prices a flux the ODE layer supplies instead: 0.5 / 1.0 / 2.0 / 5.0 mmol/gDCW/h costs
   -1.87% / -4.38% / -10.73% / -42.03% of growth, cheaper at low flux than the same export
   costs aerobically (-3.10% / -6.20% / -31.02% at 0.5 / 1.0 / 5.0).
2. The NADH goes somewhere, and it is nameable. Blocking both NADH-linked ammonium routes
   (``r_0470``, ``r_0472``) changes growth and yield by nothing at all, so those two are not
   it. The route that IS carrying it is a cytosolic transhydrogenase in disguise: GPD makes
   glycerol 3-phosphate, ``r_0489`` dephosphorylates it, and GCY1 ``r_0487`` oxidises the
   glycerol back with NADP+ -- so the glycerol is made and re-eaten instead of leaving. Close
   all three and the model secretes 0.418 mmol/gDCW/h of glycerol at a cost of 6.99% growth.
   Nothing here closes them: ``r_0487`` is vendored, irreversible as written, and 0.418 is
   not Verduyn's measured yield either.

HYPOXIA IS NOT A MILD VERSION OF THIS, and an earlier draft of this module had it exactly
backwards. Below qO2 ~ 0.105 mmol/gDCW/h the unsupplemented model's growth is EXACTLY
linear in oxygen, mu = 1.8246 * qO2, because oxygen there is a biosynthetic reagent for the
five requirements above and not an energy source: 0.548 mmol O2 per gDCW is the model's
minimum. Supplementing lifts qO2 = 0.001 from 0.001825 to 0.202429 /h, a factor of 110.9.
So which regime a hypoxic number belongs to is the whole question, and
:func:`hypoxic_supplement_effect` reports both.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra

from .physiology import BIOMASS_REACTION, GLUCOSE_EXCHANGE, OXYGEN_EXCHANGE, cap_uptake
from .solver import growth_or_none

__all__ = [
    "ANAEROBIC_COFACTOR_PSEUDOREACTION",
    "ANAEROBIC_SUPPLEMENTS",
    "COFACTOR_PSEUDOREACTION",
    "ESSENTIAL_OXYGEN_REACTIONS",
    "GLYCEROL_EXCHANGE",
    "HEME_A_CYTOSOL",
    "MOLECULAR_OXYGEN_METABOLITES",
    "NGAM_REACTION",
    "REVERSIBLE_ESTERASES",
    "VENDORED_NGAM",
    "VERDUYN_1990_ANAEROBIC",
    "AnaerobicReference",
    "AnaerobicReport",
    "AnaerobicSupplement",
    "OxygenRequirement",
    "add_anaerobic_cofactor_pseudoreaction",
    "anaerobic_constraints",
    "anaerobic_glycerol_price",
    "anaerobic_model",
    "biomass_yield",
    "constrain_esterases_to_hydrolysis",
    "diagnose_oxygen_essentiality",
    "hypoxic_model",
    "hypoxic_supplement_effect",
    "maintenance_sensitivity",
    "oxygen_consuming_reactions",
    "oxygen_metabolite_ids",
    "validate_anaerobic",
]

COFACTOR_PSEUDOREACTION = "r_4598"
ANAEROBIC_COFACTOR_PSEUDOREACTION = "r_4598_anaerobic"
HEME_A_CYTOSOL = "s_3714"
HEME_A_EXCHANGE = "r_4780"
ETHANOL_EXCHANGE = "r_1761"
GLYCEROL_EXCHANGE = "r_1808"
NGAM_REACTION = "r_4046"

# Vendored yeast-GEM v9.0.2 pins r_4046 at (0.7, 0.7) mmol ATP/gDCW/h, aerobic or not.
VENDORED_NGAM = 0.7

_GLUCOSE_G_PER_MMOL = 0.18016
_ETHANOL_G_PER_MMOL = 0.04607

# What oxygen_metabolite_ids() resolves to on yeast-GEM v9.0.2: c, er, e, m, p, erm.
MOLECULAR_OXYGEN_METABOLITES: tuple[str, ...] = (
    "s_1275", "s_1276", "s_1277", "s_1278", "s_1279", "s_2817",
)


def oxygen_metabolite_ids(model: cobra.Model) -> tuple[str, ...]:
    """Every molecular-oxygen species in ``model``, read from formula and charge.

    Not a hardcoded id list, because superoxide carries the SAME formula ``O2`` -- s_3813
    and s_3931 in v9.0.2 -- and charge is what separates it. An id list would also go
    stale silently on a model with a compartment this one lacks.
    """
    return tuple(sorted(
        met.id for met in model.metabolites
        if met.formula == "O2" and (met.charge or 0) == 0
    ))


@dataclass(frozen=True)
class OxygenRequirement:
    """One O2-consuming reaction whose loss stops growth, and what it is required for."""

    reaction_id: str
    name: str
    genes: str
    requirement: str
    resolved_by: str


# Measured on yeast-GEM v9.0.2 by zeroing each O2 consumer in turn; every entry took growth
# from 0.887685 /h to below 1e-12, and the set is identical at glucose <= 10, 11.1 and 1000.
ESSENTIAL_OXYGEN_REACTIONS: tuple[OxygenRequirement, ...] = (
    OxygenRequirement("r_1011", "squalene epoxidase (NADP)", "ERG1 (YGR175C)",
                      "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_0317", "lanosterol 14-alpha-demethylase",
                      "ERG11 (YHR007C) and NCP1 (YHR042W)",
                      "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_0238", "C-4 methyl sterol oxidase", "ERG25 (YGR060W)",
                      "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_0239", "C-4 methyl sterol oxidase", "ERG25 (YGR060W)",
                      "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_0240", "C-4 methyl sterol oxidase", "ERG25 (YGR060W)",
                      "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_0241", "C-4 sterol methyl oxidase (4,4-dimethylzymosterol)",
                      "ERG25 (YGR060W)", "ergosterol", "ergosterol uptake r_1757"),
    OxygenRequirement("r_2182", "palmitoyl-CoA desaturase (C16:0 -> C16:1)", "OLE1 (YGL055W)",
                      "unsaturated fatty acid", "palmitoleate uptake r_1994"),
    OxygenRequirement("r_2183", "stearoyl-CoA desaturase (C18:0 -> C18:1)", "OLE1 (YGL055W)",
                      "unsaturated fatty acid", "oleate uptake r_2189"),
    OxygenRequirement("r_0694", "tryptophan 2,3-dioxygenase", "BNA2 (YJR078W)",
                      "NAD de novo", "nicotinate uptake r_1967"),
    OxygenRequirement("r_0671", "kynurenine 3-monooxygenase", "BNA4 (YBL098W)",
                      "NAD de novo", "nicotinate uptake r_1967"),
    OxygenRequirement("r_0058", "3-hydroxyanthranilate 3,4-dioxygenase", "BNA1 (YJR025C)",
                      "NAD de novo", "nicotinate uptake r_1967"),
    OxygenRequirement("r_0937", "polyamine oxidase", "FMS1 (YMR020W)",
                      "beta-alanine -> pantothenate -> CoA", "pantothenate uptake r_1548"),
    OxygenRequirement("r_0304", "coproporphyrinogen oxidase (O2 required)", "HEM13 (YDR044W)",
                      "heme", "heme a leaves the anaerobic biomass"),
    OxygenRequirement("r_0942", "protoporphyrinogen oxidase", "HEM14 (YER014W)",
                      "heme", "heme a leaves the anaerobic biomass"),
    OxygenRequirement("r_0530", "heme O monooxygenase",
                      "COX15 (YER141W) or (ARH1 (YDR376W) and YAH1 (YPL252C))",
                      "heme a", "heme a leaves the anaerobic biomass"),
    OxygenRequirement("r_1977", "O2 transport, cytoplasm-ER", "",
                      "carrier for the sterol pathway", "no requirement of its own"),
    OxygenRequirement("r_1978", "O2 transport, cytoplasm-mitochondrion", "",
                      "carrier for the heme pathway", "no requirement of its own"),
    OxygenRequirement("r_1979", "O2 transport, extracellular-cytoplasm", "",
                      "carrier for every cytosolic consumer", "no requirement of its own"),
    OxygenRequirement("r_3531", "O2 transport, cytoplasm-ER membrane", "",
                      "carrier for OLE1", "no requirement of its own"),
)


@dataclass(frozen=True)
class AnaerobicSupplement:
    """A medium component an anaerobic culture is actually given, and why."""

    exchange_id: str
    component: str
    requirement: str
    source: str


ANAEROBIC_SUPPLEMENTS: tuple[AnaerobicSupplement, ...] = (
    AnaerobicSupplement(
        "r_1757", "ergosterol", "sterol ring; ERG1/ERG11/ERG25 all consume O2",
        "Andreasen & Stier 1953, PMID 13034889; supplemented in Verduyn 1990, PMID 1975265",
    ),
    AnaerobicSupplement(
        "r_2189", "oleate (C18:1), from Tween 80", "OLE1 desaturation consumes O2",
        "Andreasen & Stier 1954, J Cell Comp Physiol 43:271-281, "
        "doi 10.1002/jcp.1030430303; "
        "Tween 80 supplement in Verduyn 1990, PMID 1975265. Polysorbate 80 is >=58% oleate "
        "by the USP-NF monograph. CAVEAT: Dekker et al. 2019, PMID 31425603, grew "
        "CEN.PK113-7D anaerobically with neither UFA synthesis nor supplementation",
    ),
    AnaerobicSupplement(
        "r_1994", "palmitoleate (C16:1), from Tween 80",
        "yeast-GEM's lipid chain pseudoreaction r_4065 demands a C16:1 chain that OLE1 "
        "cannot make without O2",
        "Tween 80 supplement in Verduyn 1990, PMID 1975265; polysorbate 80 carries up to "
        "8% palmitoleate by the USP-NF monograph. This one is a property of the model's "
        "fixed lipid composition as much as of the cell -- see the oleate caveat",
    ),
    AnaerobicSupplement(
        "r_1967", "nicotinic acid", "NAD de novo from tryptophan needs O2 three times",
        "vitamin of the Verduyn mineral medium, Verduyn et al. 1992, PMID 1523884",
    ),
    AnaerobicSupplement(
        "r_1548", "(R)-pantothenate", "CoA needs beta-alanine, and FMS1 consumes O2",
        "vitamin of the Verduyn mineral medium, Verduyn et al. 1992, PMID 1523884",
    ),
)

# IAH1 (YOR126C) is irreversible-hydrolytic in r_0369/r_0656/r_0657; these nine are written
# in the synthesis direction and left reversible, so the optimum makes esters through them.
REVERSIBLE_ESTERASES: tuple[str, ...] = (
    "r_4711", "r_4712", "r_4713", "r_4714", "r_4715", "r_4716", "r_4717", "r_4720", "r_4721",
)


@dataclass(frozen=True)
class AnaerobicReference:
    """A measured anaerobic phenotype."""

    name: str
    mu_max: float
    biomass_yield: float
    yield_dilution_rate: float
    source: str

    def implied_glucose_uptake(self, substrate_g_per_mmol: float = _GLUCOSE_G_PER_MMOL) -> float:
        """The uptake this paper's own pair implies at mu_max, mmol/gDCW/h.

        ``mu_max / (Y * M)``. Arithmetic on two reported numbers, not a third measurement,
        and it assumes the yield measured at D = 0.10 still holds at mu_max -- which the
        model itself says is optimistic, since a fixed NGAM makes yield rise with mu.
        """
        return self.mu_max / (self.biomass_yield * substrate_g_per_mmol)


VERDUYN_1990_ANAEROBIC = AnaerobicReference(
    name="anaerobic glucose-limited chemostat, supplemented mineral medium",
    mu_max=0.31,
    biomass_yield=0.10,
    yield_dilution_rate=0.10,
    source="Verduyn, Postma, Scheffers & van Dijken 1990, PMID 1975265; S. cerevisiae "
           "CBS 8066, mineral medium with ergosterol and Tween 80",
)


@dataclass(frozen=True)
class AnaerobicReport:
    """What the anaerobic model gives beside what was measured."""

    reference: AnaerobicReference
    growth_at_glucose_bound: float
    glucose_bound: float
    yield_glucose_uptake: float
    model_biomass_yield: float
    yield_relative_error: float
    glucose_for_reference_mu_max: float
    ethanol_yield_g_per_g: float
    implied_glucose_uptake: float
    mu_at_implied_uptake: float
    mu_relative_error: float
    glycerol_at_optimum: float

    def summary(self) -> str:
        return (
            f"{self.reference.name}\n"
            f"  growth at glucose <= {self.glucose_bound:g}: "
            f"{self.growth_at_glucose_bound:.6f} /h\n"
            f"  biomass yield at mu = {self.reference.yield_dilution_rate:g}: model "
            f"{self.model_biomass_yield:.4f} g/g on {self.yield_glucose_uptake:.3f} "
            f"mmol/gDCW/h, measured {self.reference.biomass_yield:.2f} g/g, "
            f"rel.err {self.yield_relative_error:+.1%}\n"
            f"  glucose needed for the measured mu_max {self.reference.mu_max:g} /h: "
            f"{self.glucose_for_reference_mu_max:.3f} mmol/gDCW/h\n"
            f"  ethanol yield there: {self.ethanol_yield_g_per_g:.4f} g/g\n"
            f"  mu at the uptake the measured pair implies "
            f"({self.implied_glucose_uptake:.3f} mmol/gDCW/h): "
            f"{self.mu_at_implied_uptake:.6f} /h against {self.reference.mu_max:g}, "
            f"rel.err {self.mu_relative_error:+.1%}\n"
            f"  glycerol secreted at the optimum: {self.glycerol_at_optimum:.6f} "
            f"mmol/gDCW/h -- an anaerobic culture makes some, this model makes none\n"
            f"  source: {self.reference.source}"
        )


def oxygen_consuming_reactions(model: cobra.Model) -> tuple[str, ...]:
    """Every non-boundary reaction whose bounds let it consume O2 in some compartment."""
    oxygen = set(oxygen_metabolite_ids(model))
    out = []
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        for met, coefficient in rxn.metabolites.items():
            if met.id not in oxygen:
                continue
            if (coefficient < 0 and rxn.upper_bound > 0) or (
                    coefficient > 0 and rxn.lower_bound < 0):
                out.append(rxn.id)
                break
    return tuple(out)


def diagnose_oxygen_essentiality(
    model: cobra.Model,
    glucose_uptake: float = 10.0,
    lethal_fraction: float = 1e-6,
) -> tuple[tuple[str, str, float], ...]:
    """Zero each O2 consumer in turn and return the ones whose loss is lethal.

    This is the measurement :data:`ESSENTIAL_OXYGEN_REACTIONS` records. Does not modify
    ``model``.

    Args:
        model: Aerobic model.
        glucose_uptake: Uptake bound, mmol/gDCW/h. Regime-dependent, so it is carried out.
        lethal_fraction: Fraction of base growth below which a deletion counts as lethal.

    Returns:
        ``(reaction_id, name, growth_after_deletion)`` for each lethal deletion.
    """
    with model as scoped:
        cap_uptake(scoped, GLUCOSE_EXCHANGE, glucose_uptake)
        scoped.objective = BIOMASS_REACTION
        base = growth_or_none(scoped)
        if not base:
            raise RuntimeError(
                f"the aerobic reference does not grow (growth={base!r}) at glucose "
                f"<= {glucose_uptake:g}, so nothing can be called essential against it"
            )
        lethal = []
        for rid in oxygen_consuming_reactions(scoped):
            reaction = scoped.reactions.get_by_id(rid)
            with scoped:
                reaction.bounds = (0.0, 0.0)
                growth = growth_or_none(scoped) or 0.0
            if growth / base < lethal_fraction:
                lethal.append((rid, reaction.name, float(growth)))
    return tuple(lethal)


def add_anaerobic_cofactor_pseudoreaction(model: cobra.Model) -> cobra.Reaction:
    """Add a copy of ``r_4598`` without heme a, and close the vendored one. In place.

    Anaerobic yeast synthesises no heme -- HEM13 (r_0304) and HEM14 (r_0942) both consume
    O2 here, and COX15 (r_0530) consumes another to reach heme a -- and it contains no
    heme a, because it assembles no cytochrome oxidase. Snoek & Steensma 2007
    (PMID 17192845) review the O2-requiring steps anaerobic yeast has to do without.

    Added as a separate reaction rather than by editing ``r_4598`` so the change is one
    reaction id a reader can find, and the vendored stoichiometry stays intact. It is the
    only composition change here, and the ONE ASSERTED thing in this module: the
    coefficient removed is 1e-06 mmol/gDCW of heme a.

    Neither this reaction nor the vendored one mass-balances, and that is not a defect
    introduced here: ``s_4205`` ("cofactor") carries no formula, so a pseudoreaction
    pointing at it cannot balance. What is checkable is the DIFFERENCE, which is exactly
    1e-06 of C49H55FeN4O6 at charge -3, and the tests check that.

    Raises:
        ValueError: if it is already installed.
    """
    if model.reactions.has_id(ANAEROBIC_COFACTOR_PSEUDOREACTION):
        raise ValueError(
            f"{ANAEROBIC_COFACTOR_PSEUDOREACTION} is already on this model"
        )
    aerobic = model.reactions.get_by_id(COFACTOR_PSEUDOREACTION)
    stoichiometry = {
        met: coefficient for met, coefficient in aerobic.metabolites.items()
        if met.id != HEME_A_CYTOSOL
    }
    if len(stoichiometry) == len(aerobic.metabolites):
        raise KeyError(
            f"{COFACTOR_PSEUDOREACTION} does not contain {HEME_A_CYTOSOL} (heme a), so this "
            "model is not the one the anaerobic mode was measured on"
        )
    anaerobic = cobra.Reaction(
        ANAEROBIC_COFACTOR_PSEUDOREACTION,
        name="cofactor pseudoreaction, anaerobic (no heme a)",
        lower_bound=aerobic.lower_bound,
        upper_bound=aerobic.upper_bound,
    )
    anaerobic.subsystem = aerobic.subsystem
    model.add_reactions([anaerobic])
    anaerobic.add_metabolites(stoichiometry)
    aerobic.bounds = (0.0, 0.0)
    return anaerobic


def constrain_esterases_to_hydrolysis(model: cobra.Model) -> tuple[str, ...]:
    """Restrict the nine reversible IAH1 esterases to hydrolysis. In place.

    These nine are written ``acid + ethanol + H+ <=> H2O + ester``, so SYNTHESIS is the
    forward direction and ``upper_bound = 0`` is what blocks it; hydrolysis is the negative
    flux that survives. Left free, the anaerobic optimum at glucose <= 10 runs ``r_4713``
    forward at +0.7848 mmol/gDCW/h and dumps that much diethyl succinate -- 10.5% of the
    glucose carbon -- and the ethanol yield reads 0.3816 g/g instead of 0.4453.

    THE CONSTRAINT IS NOT FREE: it costs 1.16% of anaerobic growth, 0.202264 -> 0.199918 /h,
    so the headline anaerobic number is measured with it applied. The evidence that the ester
    dump is the artefact rather than the constraint is inside the model: the same gene's other
    three reactions (r_0369, r_0656, r_0657) are irreversible in the hydrolysis direction.
    Ester synthesis stays available by the route that pays for it, the ATF1/ATF2 acyl-CoA
    transferases.

    Returns:
        The reaction ids constrained.
    """
    constrained = []
    for rid in REVERSIBLE_ESTERASES:
        if not model.reactions.has_id(rid):
            continue
        reaction = model.reactions.get_by_id(rid)
        reaction.upper_bound = 0.0
        constrained.append(rid)
    return tuple(constrained)


def anaerobic_constraints(
    model: cobra.Model,
    glucose_uptake: float = 10.0,
    supplements: tuple[AnaerobicSupplement, ...] = ANAEROBIC_SUPPLEMENTS,
    supplement_bound: float = 1000.0,
    hydrolysis_only_esterases: bool = True,
    import_heme_a: bool = False,
) -> cobra.Model:
    """Put ``model`` into anaerobic mode in place. Use inside a ``with model:`` block.

    Args:
        model: Model to constrain.
        glucose_uptake: Uptake bound, mmol/gDCW/h.
        supplements: Medium components to open. Each is justified in
            :data:`ANAEROBIC_SUPPLEMENTS`; passing an empty tuple reproduces the
            unsupplemented medium, which does not grow.
        supplement_bound: Uptake bound for each supplement. Unlimited by default, because
            no supplement concentration was measured here. Harmless in every regime this
            module uses -- with O2 off or capped the optimum draws under 0.05 of each, and
            growth is identical at a bound of 0.05 -- but three of the five are lipids, so
            reopening oxygen on a supplemented model turns them into fuel: at glucose <= 10
            with O2 free, growth goes 0.887685 -> 8.332145 /h. Do not reuse this medium
            aerobically without a measured bound.
        hydrolysis_only_esterases: Apply :func:`constrain_esterases_to_hydrolysis`.
        import_heme_a: FABRICATES A MEDIUM COMPONENT. Opens the heme a exchange r_4780
            instead of taking heme a out of the biomass. No anaerobic medium contains
            heme a; this exists only so the two routes can be compared, and it gives the
            same growth to 13 decimal places.
    """
    cap_uptake(model, GLUCOSE_EXCHANGE, glucose_uptake)
    model.reactions.get_by_id(OXYGEN_EXCHANGE).bounds = (0.0, 0.0)
    for supplement in supplements:
        cap_uptake(model, supplement.exchange_id, supplement_bound)
    if import_heme_a:
        cap_uptake(model, HEME_A_EXCHANGE, supplement_bound)
    else:
        add_anaerobic_cofactor_pseudoreaction(model)
    if hydrolysis_only_esterases:
        constrain_esterases_to_hydrolysis(model)
    model.objective = BIOMASS_REACTION
    return model


def anaerobic_model(model: cobra.Model, **kwargs) -> cobra.Model:
    """A copy of ``model`` in anaerobic mode. The input is not modified.

    Keyword arguments are :func:`anaerobic_constraints`'.
    """
    out = model.copy()
    return anaerobic_constraints(out, **kwargs)


def hypoxic_model(
    model: cobra.Model,
    oxygen_uptake: float,
    glucose_uptake: float = 10.0,
    supplements: tuple[AnaerobicSupplement, ...] = (),
) -> cobra.Model:
    """A copy of ``model`` at a capped oxygen uptake. The input is not modified.

    The default is the bare mineral medium every other regime in this package uses, and in
    it low oxygen is a BIOSYNTHETIC limitation rather than an energetic one: measured at
    glucose <= 10, qO2 of 0.001 / 0.01 / 0.05 / 0.1 gives 0.001825 / 0.018246 / 0.091231 /
    0.182462 /h, exactly linear at 1.8246 /h per mmol O2, and the line ends at qO2 ~ 0.105
    where growth reaches 0.1913 /h. Pass :data:`ANAEROBIC_SUPPLEMENTS` for the supplemented
    medium, where that limitation is gone -- see :func:`hypoxic_supplement_effect`.

    THE SUPPLEMENTS ARE OPENED UNBOUNDED, and three of the five are lipids, so at high oxygen
    they stop being a requirement and become fuel: measured supplement carbon is 1.1% of
    glucose carbon at qO2 = 0.001 and 2.5% at qO2 = 5, but with oxygen free the optimum eats
    40.8 mmol/gDCW/h of oleate and grows at 8.33 /h. Do not read a supplemented number above
    the oxygen range this module measures without checking the uptakes.

    Raises:
        ValueError: on a non-positive uptake, which is anaerobiosis and needs
            :func:`anaerobic_model` instead of a bound that silently returns zero growth.
    """
    if oxygen_uptake <= 0.0:
        raise ValueError(
            f"oxygen_uptake={oxygen_uptake!r} is anaerobic, not hypoxic. Use "
            "anaerobic_model(), which supplies what the cell can no longer make; a zero "
            "bound here would return growth 0 and read as a result"
        )
    out = model.copy()
    cap_uptake(out, GLUCOSE_EXCHANGE, glucose_uptake)
    cap_uptake(out, OXYGEN_EXCHANGE, oxygen_uptake)
    for supplement in supplements:
        cap_uptake(out, supplement.exchange_id, 1000.0)
    out.objective = BIOMASS_REACTION
    return out


def hypoxic_supplement_effect(
    model: cobra.Model,
    oxygen_uptake: float,
    glucose_uptake: float = 10.0,
) -> tuple[float, float, float]:
    """Growth at one oxygen uptake, bare and supplemented, and the ratio.

    The number that says which of the two a hypoxic result belongs to. At qO2 = 0.001 the
    ratio is 110.9; by qO2 = 1.0 it is 1.09.

    Returns:
        ``(bare, supplemented, supplemented / bare)``.
    """
    bare = growth_or_none(hypoxic_model(model, oxygen_uptake, glucose_uptake)) or 0.0
    supplemented = growth_or_none(hypoxic_model(
        model, oxygen_uptake, glucose_uptake, supplements=ANAEROBIC_SUPPLEMENTS)) or 0.0
    if bare <= 0.0:
        raise RuntimeError(
            f"the bare model does not grow at qO2={oxygen_uptake:g}, so the ratio would "
            "be a division by zero rather than a number"
        )
    return bare, supplemented, supplemented / bare


def anaerobic_glycerol_price(
    model: cobra.Model,
    fluxes: tuple[float, ...] = (0.5, 1.0, 2.0, 5.0),
    glucose_uptake: float = 10.0,
    **kwargs,
) -> tuple[tuple[float, float, float], ...]:
    """What an imposed anaerobic glycerol export costs in growth.

    The anaerobic optimum makes no glycerol and a real anaerobic culture does, so the flux
    has to come from outside the GEM; this prices it rather than inventing a yield. Measured
    at glucose <= 10: -1.87% / -4.38% / -10.73% / -42.03% at 0.5 / 1.0 / 2.0 / 5.0.

    Keyword arguments are :func:`anaerobic_constraints`'.

    Returns:
        ``(flux, growth, relative_cost)`` per requested flux, cost negative for a loss.
    """
    out = anaerobic_model(model, glucose_uptake=glucose_uptake, **kwargs)
    base = growth_or_none(out)
    if not base:
        raise RuntimeError(
            f"the anaerobic model does not grow (growth={base!r}), so nothing can be "
            "priced against it"
        )
    priced = []
    for flux in fluxes:
        with out as scoped:
            scoped.reactions.get_by_id(GLYCEROL_EXCHANGE).bounds = (flux, flux)
            growth = growth_or_none(scoped) or 0.0
        priced.append((float(flux), float(growth), float(growth) / base - 1.0))
    return tuple(priced)


def maintenance_sensitivity(
    model: cobra.Model,
    ngam_values: tuple[float, ...] = (0.7, 1.0, 1.5, 2.0),
    growth_rate: float = 0.10,
    **kwargs,
) -> tuple[tuple[float, float, float], ...]:
    """Biomass yield at a fixed growth rate as the vendored NGAM is varied. DIAGNOSTIC ONLY.

    Nothing here adopts a value: r_4046 stays at the vendored 0.7 everywhere else in this
    module. It exists because the anaerobic yield comes out 7.2% high and this says where
    that could live -- measured 0.1072 / 0.1042 / 0.0995 / 0.0953 g/g at NGAM 0.7 / 1.0 /
    1.5 / 2.0. Reading a maintenance coefficient off that agreement would be fitting, not
    measuring, so it is not done.

    Keyword arguments are :func:`anaerobic_constraints`'.

    Returns:
        ``(ngam, yield_g_per_g, glucose_uptake)`` per value.
    """
    out = anaerobic_model(model, glucose_uptake=1000.0, **kwargs)
    rows = []
    for ngam in ngam_values:
        with out as scoped:
            scoped.reactions.get_by_id(NGAM_REACTION).bounds = (ngam, ngam)
            yield_g_per_g, uptake, _ = biomass_yield(scoped, growth_rate)
        rows.append((float(ngam), yield_g_per_g, uptake))
    return tuple(rows)


def biomass_yield(
    model: cobra.Model,
    growth_rate: float,
    substrate_exchange: str = GLUCOSE_EXCHANGE,
    substrate_g_per_mmol: float = _GLUCOSE_G_PER_MMOL,
) -> tuple[float, float, cobra.Solution]:
    """Maximal biomass yield at a fixed growth rate: the least substrate that reaches it.

    Fixing growth and minimising uptake is the chemostat question -- at a dilution rate the
    rate is set and the yield is what is being measured -- and it is not the same LP as
    maximising growth at a fixed uptake.

    Args:
        model: Model already in the regime of interest.
        growth_rate: Growth to hold, /h.
        substrate_exchange: Exchange to minimise uptake through.
        substrate_g_per_mmol: Molar mass, g/mmol, for the g/g yield.

    Returns:
        ``(yield_g_per_g, substrate_uptake, solution)``.

    Raises:
        RuntimeError: if the growth rate is infeasible in this regime.
    """
    with model as scoped:
        biomass = scoped.reactions.get_by_id(BIOMASS_REACTION)
        biomass.bounds = (growth_rate, growth_rate)
        exchange = scoped.reactions.get_by_id(substrate_exchange)
        scoped.objective = exchange
        # Uptake is negative on an unsplit exchange, so least uptake is the maximum.
        scoped.objective.direction = "max"
        if growth_or_none(scoped) is None:
            raise RuntimeError(
                f"growth {growth_rate:g} /h is infeasible in this regime, so no yield "
                "can be reported for it"
            )
        solution = scoped.optimize()
        uptake = abs(float(solution.fluxes[substrate_exchange]))
    if uptake <= 0.0:
        raise RuntimeError(
            f"{substrate_exchange} carries no uptake at growth {growth_rate:g} /h; the "
            "yield would be a division by zero rather than a number"
        )
    return growth_rate / (uptake * substrate_g_per_mmol), uptake, solution


def validate_anaerobic(
    model: cobra.Model,
    reference: AnaerobicReference = VERDUYN_1990_ANAEROBIC,
    glucose_bound: float = 10.0,
    **kwargs,
) -> AnaerobicReport:
    """Put ``model`` in anaerobic mode and compare it with a measured phenotype.

    Reports rather than asserts, and the reason is that only one of the two reference
    numbers is a like-for-like comparison. The yield at a dilution rate is: fix mu, minimise
    uptake. mu_max is not -- this model has no maximum growth rate of its own, only one set
    by the uptake bound it is given -- so mu_max is answered twice, once as the uptake it
    would take to reach the measured value and once as the growth reached at the uptake the
    paper's own two numbers imply.

    Keyword arguments are :func:`anaerobic_constraints`'.
    """
    out = anaerobic_model(model, glucose_uptake=glucose_bound, **kwargs)
    growth = growth_or_none(out)
    if growth is None:
        raise RuntimeError("the anaerobic model is infeasible, so there is nothing to check")
    glycerol = float(out.optimize().fluxes[GLYCEROL_EXCHANGE])

    yield_g_per_g, uptake, _ = biomass_yield(out, reference.yield_dilution_rate)
    # mu_max is above what glucose_bound allows, so that bound comes off for both questions.
    implied = reference.implied_glucose_uptake()
    with out as unlimited:
        cap_uptake(unlimited, GLUCOSE_EXCHANGE, 1000.0)
        _, uptake_at_mu_max, solution = biomass_yield(unlimited, reference.mu_max)
    with out as at_implied:
        cap_uptake(at_implied, GLUCOSE_EXCHANGE, implied)
        mu_at_implied = growth_or_none(at_implied) or 0.0
    ethanol_yield = (
        float(solution.fluxes[ETHANOL_EXCHANGE]) * _ETHANOL_G_PER_MMOL
        / (uptake_at_mu_max * _GLUCOSE_G_PER_MMOL)
    )
    return AnaerobicReport(
        reference=reference,
        growth_at_glucose_bound=float(growth),
        glucose_bound=glucose_bound,
        yield_glucose_uptake=uptake,
        model_biomass_yield=yield_g_per_g,
        yield_relative_error=yield_g_per_g / reference.biomass_yield - 1.0,
        glucose_for_reference_mu_max=uptake_at_mu_max,
        ethanol_yield_g_per_g=ethanol_yield,
        implied_glucose_uptake=implied,
        mu_at_implied_uptake=float(mu_at_implied),
        mu_relative_error=float(mu_at_implied) / reference.mu_max - 1.0,
        glycerol_at_optimum=glycerol,
    )
