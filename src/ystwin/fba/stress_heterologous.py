"""PHB installation, an export cost, and a precursor drain that is not vacuous.

Three gaps a four-auditor GEM audit found on the product side.

1. `data/pathways/phb.toml` has declared the phaA/phaB/phaC chemistry since 2026-09-04 and
   nothing ever installed it. :func:`add_phb_pathway` does, following `fba/carotenoid.py`.
2. NOTHING CHARGED FOR EXPORT -- not a protein, not a small molecule. The audit's wider
   claim that there is no secretion cost at all is half wrong and the half that survives is
   the load-bearing one: `fba/secretion.py` already charges the ER oxidative FOLDING term
   (one O2 in, one H2O2 out per disulfide, PMID 16407158), and `fba/insulin.py` says in its
   own docstring that leaving the product on a demand "is the largest omission here"
   because "adding an extracellular species and a costless transporter ... would name a
   transport step and give it a price of zero". :func:`add_atp_coupled_export` gives that
   step a price, at yeast-GEM's own ABC-system rate, and refuses rather than inventing.
3. `fba/audit.py` advertises a precursor budget and measures exactly zero for the declared
   beta-carotene node. :func:`native_drain` says why and reports it as vacuous rather than
   as a zero that reads like "no native demand".

Every id below was read out of the vendored `data/gem/yeast-GEM.xml.gz`.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra

from .fva import GLUCOSE_EXCHANGE
from .physiology import cap_uptake
from .solver import FVA_PROCESSES, configure as pin_solver, growth_or_none

__all__ = [
    "ATP_EQUIVALENTS_PER_RESIDUE",
    "ATP_PER_ABC_EXPORT",
    "HBCOA_ID",
    "MFALPHA1_PREPRO_RESIDUES",
    "MFALPHA1_PREPRO_WITH_SPACER_RESIDUES",
    "PHB_DEMAND_ID",
    "PHB_REACTION_IDS",
    "PHB_REPEAT_UNIT_ID",
    "SECRETION_NOT_CHARGED",
    "PrecursorDrain",
    "SecretionCostUnknown",
    "add_atp_coupled_export",
    "add_phb_pathway",
    "leader_peptide_atp_equivalents",
    "native_drain",
]

# ---------------------------------------------------------------------------
# PHB
# ---------------------------------------------------------------------------

PHB_REPEAT_UNIT_ID = "phb_c"
HBCOA_ID = "hbcoa_c"
PHB_DEMAND_ID = "DM_phb_c"
PHB_REACTION_IDS = ("PHAA", "PHAB", "PHAC", PHB_DEMAND_ID)

# Cytosolic ids for the precursor and cofactors phb.toml declares.
_ACCOA_C, _AACCOA_C, _COA_C = "s_0373", "s_0367", "s_0529"
_NADPH_C, _NADP_C, _PROTON_C = "s_1212", "s_1207", "s_0794"

_PHB_FORMULA, _PHB_CHARGE = "C4H6O2", 0
"""One 3-hydroxybutyrate repeat unit, KEGG C06143 ``(C4H6O2)n``, per data/pathways/phb.toml.

Checked against, never used to build the species: the formula is derived below from the
host's own cofactors. The polymer has no molar mass, so everything is per repeat unit.
"""

_HILL_FIRST = ("C", "H")

# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

ATP_PER_ABC_EXPORT = 1.0
"""ATP per molecule exported by an ABC transporter -- yeast-GEM v9.0.2's OWN convention.

Not a measured turnover. All twelve reactions the vendored SBML names "via ABC system"
(r_0964, r_1028, r_2219-r_2228) are written `ATP + H2O + X_in -> ADP + H+ + phosphate +
X_out` with an ATP coefficient of exactly -1, and the test file pins that. NONE OF THE
TWELVE IS A PLASMA-MEMBRANE EXPORT -- they are vacuolar and peroxisomal transports -- so
this is the host model's own convention carried across, not a measurement of secretion.
PMID 33062054 shows beta-carotene export in S. cerevisiae runs through an endogenous ABC
transporter and rises with ATP supply; it reports no per-molecule stoichiometry.
"""

ATP_EQUIVALENTS_PER_RESIDUE = 4.0
"""ATP-equivalents to polymerise one residue that is later cleaved off and discarded.

Martin WF 2025, PMID 40562331: "The pure cost of synthesizing peptide bonds is 4 ATP each:
2 from PPi formation at aminoacyl tRNA synthesis ... and 1 GTP each for the two elongation
factors." All four are net here, unlike `fba/insulin.py`'s GTP_PER_RESIDUE = 2.0, because a
leader charged this way never draws on the host's aminoacyl-tRNA pool and so never pays the
synthetase half through the model's own reactions.
"""

MFALPHA1_PREPRO_RESIDUES = 85
"""Residues of the S. cerevisiae MFalpha1 pre-pro leader, to the Kex2 site.

FETCHED, NOT REMEMBERED: UniProt P01149 (MFAL1_YEAST, 165 aa) gives `SIGNAL 1..19` and
`PROPEP 20..89`, and its sequence carries the dibasic Kex2 site at Lys84-Arg85 with the
Ste13 spacer EAEA at 86-89. So 85 to the Kex2 cut and
:data:`MFALPHA1_PREPRO_WITH_SPACER_RESIDUES` with the spacer. WHICH ONE A CONSTRUCT CARRIES
IS A PROPERTY OF THE CONSTRUCT and no measurement here settles it, which is why nothing
defaults to either.
"""

MFALPHA1_PREPRO_WITH_SPACER_RESIDUES = 89

SECRETION_NOT_CHARGED = (
    ("COPII vesicle budding", "Sar1 GTP is per vesicle; a vesicle carries an unmeasured "
                              "number of cargo molecules, so there is no per-molecule rate"),
    ("Kar2/BiP translocation", "post-translational translocation costs one Hsp70 ATP cycle "
                               "per BiP binding site and the sites per chain are not measured"),
    ("N-glycosylation", "the dolichol cycle is in no vendored model here and the glycan "
                        "occupancy of a given construct is not measured"),
    ("signal peptidase and Kex2", "proteolysis, so no cofactor; the cost is the leader's "
                                  "translation, which leader_peptide_atp_equivalents gives"),
    ("the leader's amino acids", "only the peptide-bond energy is charged, not the carbon "
                                 "and nitrogen of residues whose sequence the caller has "
                                 "not supplied; this biases a ceiling UPWARD"),
)
"""What an export cost SHOULD include and this module refuses to invent, with the reason.

The one secretory term that IS charged elsewhere is oxidative folding -- see
`fba/secretion.py`, which this module deliberately does not duplicate.
"""

_ATP_C, _H2O_C, _ADP_C, _PI_C = "s_0434", "s_0803", "s_0394", "s_1322"


class SecretionCostUnknown(ValueError):
    """Raised when export is asked for with no stoichiometry anybody has measured."""


# ---------------------------------------------------------------------------


def _resolve(model: cobra.Model, base_id: str, compartment: str = "c") -> str:
    """Map a bare Yeast9 id onto whatever convention this model uses, or fail by name."""
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(
        f"required metabolite {base_id!r} not found in {model.id!r} "
        f"(tried {base_id}, {base_id}[{compartment}], {base_id}_{compartment})"
    )


def _hill(elements: dict[str, float]) -> str:
    """Render an element count in Hill order, dropping zeros."""
    def render(symbol: str) -> str:
        count = elements[symbol]
        count = int(count) if float(count).is_integer() else count
        return symbol if count == 1 else f"{symbol}{count}"

    ordered = [s for s in _HILL_FIRST if elements.get(s)]
    ordered += sorted(s for s in elements if s not in _HILL_FIRST and elements[s])
    return "".join(render(s) for s in ordered)


def _combine(*terms) -> tuple[dict[str, float], float]:
    """Sum ``(metabolite | (elements, charge), sign)`` pairs into an element count and charge."""
    elements: dict[str, float] = {}
    charge = 0.0
    for term, sign in terms:
        counts, term_charge = term if isinstance(term, tuple) else (term.elements, term.charge)
        for symbol, count in counts.items():
            elements[symbol] = elements.get(symbol, 0.0) + sign * count
        charge += sign * (term_charge or 0)
    return {s: v for s, v in elements.items() if v}, charge


def _add_reaction(model, rid, name, stoichiometry, gene_rule, bounds=(0.0, 1000.0),
                  subsystem=""):
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.subsystem = subsystem
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gene_rule
    return rxn


def add_phb_pathway(
    model: cobra.Model,
    accoa_id: str = _ACCOA_C,
    aaccoa_id: str = _AACCOA_C,
    coa_id: str = _COA_C,
    nadph_id: str = _NADPH_C,
    nadp_id: str = _NADP_C,
    proton_id: str = _PROTON_C,
    heterologous_thiolase: bool = False,
) -> cobra.Model:
    """Return a COPY of ``model`` carrying phaB/phaC, a PHB demand, and optionally phaA.

    Chemistry from data/pathways/phb.toml, which sources it to PMID 23009357.

    **phaA is OFF by default and that is a correctness fix, not a simplification.** Its
    chemistry is identical to the host's own thiolase r_0103, and a plain GEM has no enzyme
    identity, so a duplicate can only ever add a cycle. Measured at glucose <= 10 / O2 free:
    installing it moves r_0103's FVA range from [-14.358, 0] to [-1000, 10.961], so anyone
    reading ERG10 off a PHB-carrying model reads the loop rather than the biology. It buys
    nothing for that price -- the PHB ceiling is 10.960766 mmol/gDCW/h with phaA installed,
    zeroed, or absent, identical to twelve figures.

    The analogy to CrtE in `fba/carotenoid.py` does NOT hold and used to be claimed here:
    CrtE duplicates r_0461, which is irreversible, so it closes no cycle. r_0103 is
    reversible, so phaA does.

    Set ``heterologous_thiolase=True`` for an enzyme-constrained host, where phaA carries its
    own kcat and protein cost and is therefore distinguishable from the native enzyme.

    The two new species' formulas are DERIVED by solving each reaction's own mass balance
    against the host's cofactors, never typed and never copied off another species, so the
    pathway is balanced by construction on a host with any protonation convention.

    Args:
        model: Host model to extend. Not modified.
        accoa_id, aaccoa_id, coa_id, nadph_id, nadp_id, proton_id: Cytosolic base ids.

    Raises:
        ValueError: if the pathway is already installed, if the host's CoA bookkeeping does
            not yield the C4H6O2 repeat unit, or if any added reaction comes out unbalanced.
        KeyError: naming any precursor or cofactor the host model does not carry.
    """
    if model.metabolites.has_id(PHB_REPEAT_UNIT_ID):
        raise ValueError("the PHB pathway is already installed on this model")

    out = model.copy()
    accoa = _resolve(out, accoa_id)
    aaccoa = _resolve(out, aaccoa_id)
    coa = _resolve(out, coa_id)
    nadph = _resolve(out, nadph_id)
    nadp = _resolve(out, nadp_id)
    proton = _resolve(out, proton_id)
    get = out.metabolites.get_by_id

    # phaB's balance solved for its unknown product, and phaC's for the repeat unit. The
    # hydride plus the proton is the host's own statement of what a 2-electron reduction adds.
    hbcoa_elements, hbcoa_charge = _combine(
        (get(aaccoa), 1), (get(nadph), 1), (get(proton), 1), (get(nadp), -1))
    phb_elements, phb_charge = _combine(
        ((hbcoa_elements, hbcoa_charge), 1), (get(coa), -1))
    if (_hill(phb_elements), phb_charge) != (_PHB_FORMULA, _PHB_CHARGE):
        raise ValueError(
            f"this host's thioester bookkeeping yields {_hill(phb_elements)!r} at charge "
            f"{phb_charge:g} for the 3-hydroxybutyrate repeat unit, not {_PHB_FORMULA} at "
            f"{_PHB_CHARGE} (KEGG C06143). Check {aaccoa!r} and {coa!r} before installing")

    out.add_metabolites([
        cobra.Metabolite(HBCOA_ID, name="(R)-3-hydroxybutanoyl-CoA",
                         formula=_hill(hbcoa_elements), charge=hbcoa_charge, compartment="c"),
        cobra.Metabolite(PHB_REPEAT_UNIT_ID, name="poly-3-hydroxybutyrate repeat unit",
                         formula=_PHB_FORMULA, charge=_PHB_CHARGE, compartment="c"),
    ])

    # phaA duplicates the reversible native r_0103, so in a plain GEM it only adds a cycle.
    if heterologous_thiolase:
        _add_reaction(
            out, "PHAA", "acetyl-CoA C-acetyltransferase (phaA)",
            {accoa: -2, aaccoa: 1, coa: 1},
            "phaA", subsystem="PHB biosynthesis",
        )
    # phaB: NADPH-dependent acetoacetyl-CoA reductase, EC 1.1.1.36.
    _add_reaction(
        out, "PHAB", "acetoacetyl-CoA reductase (phaB)",
        {aaccoa: -1, nadph: -1, proton: -1, HBCOA_ID: 1, nadp: 1},
        "phaB", subsystem="PHB biosynthesis",
    )
    # phaC: PHA synthase; the polymer grows by one repeat unit and CoA leaves.
    _add_reaction(
        out, "PHAC", "PHA synthase (phaC)",
        {HBCOA_ID: -1, PHB_REPEAT_UNIT_ID: 1, coa: 1},
        "phaC", subsystem="PHB biosynthesis",
    )
    # Assayed from the biomass and not secreted (PMID 23009357): demand, not exchange.
    out.add_boundary(out.metabolites.get_by_id(PHB_REPEAT_UNIT_ID),
                     type="demand", reaction_id=PHB_DEMAND_ID)

    unbalanced = {rid: out.reactions.get_by_id(rid).check_mass_balance()
                  for rid in PHB_REACTION_IDS
                  if rid in {r.id for r in out.reactions}
                  and not out.reactions.get_by_id(rid).boundary
                  and out.reactions.get_by_id(rid).check_mass_balance()}
    if unbalanced:
        raise ValueError(
            f"refusing to leave unbalanced reactions in {model.id!r}: {unbalanced}. The "
            "host's cofactor formulas and this chemistry do not agree")
    return out


# ---------------------------------------------------------------------------


def leader_peptide_atp_equivalents(
    residues: int, per_residue: float = ATP_EQUIVALENTS_PER_RESIDUE
) -> float:
    """ATP-equivalents to translate a leader of ``residues`` residues that is then discarded.

    There is deliberately NO default leader. Which leader a construct carries is a property
    of the construct, so the caller supplies the length -- :data:`MFALPHA1_PREPRO_RESIDUES`
    is the one this repository can cite -- and this supplies the rate.

    Raises:
        ValueError: for a negative length or rate.
    """
    if residues < 0:
        raise ValueError(f"residues must be >= 0, got {residues}")
    if per_residue < 0:
        raise ValueError(f"per_residue must be >= 0, got {per_residue}")
    return float(residues) * float(per_residue)


def add_atp_coupled_export(
    model: cobra.Model,
    metabolite_id: str,
    *,
    atp_per_molecule: float | None,
    provenance: str,
    leader_residues: int = 0,
    reaction_id: str | None = None,
    exported_id: str | None = None,
    exchange: bool = True,
) -> cobra.Reaction:
    """Charge ATP for moving one metabolite out of the cell. Modifies ``model`` in place.

    Installs ``X + n ATP + n H2O -> X_e + n ADP + n phosphate + n H+`` and, by default, an
    exchange on the extracellular species so the export is not a dead end. ``n`` is
    ``atp_per_molecule`` plus the leader's translation, if any. One LUMPED step: it prices
    the export, it does not resolve the route, and what it leaves out is
    :data:`SECRETION_NOT_CHARGED`.

    Args:
        model: Model carrying the cargo. Modified in place; the caller is mid-construction.
        metabolite_id: The cargo. Its extracellular twin is created if absent.
        atp_per_molecule: ATP per molecule exported. ``None`` REFUSES -- pass
            :data:`ATP_PER_ABC_EXPORT` to adopt the host model's own convention, or a
            measured number with its citation.
        provenance: Where the number came from. Required, and must say ASSERTED or
            PLACEHOLDER plainly when it is neither cited nor read off the host model. It is
            written onto the reaction's notes so it stays greppable.
        leader_residues: Residues of signal/pro peptide translated and then cleaved off.
        reaction_id, exported_id: Override the generated ids. The default strips the
            cargo's own compartment suffix rather than appending to it.
        exchange: Add an exchange on the extracellular species.

    Raises:
        SecretionCostUnknown: if ``atp_per_molecule`` is ``None``.
        ValueError: on a negative cost, an empty provenance, or a second install.
    """
    if atp_per_molecule is None:
        raise SecretionCostUnknown(
            f"no ATP-per-molecule cost was supplied for exporting {metabolite_id!r}, and "
            "this module will not invent one. yeast-GEM v9.0.2 charges 1 ATP for every "
            "reaction it names 'via ABC system' -- pass ATP_PER_ABC_EXPORT to adopt that "
            "convention -- or supply a measured turnover with its citation. What is NOT "
            f"charged either way, and why: {SECRETION_NOT_CHARGED}")
    if atp_per_molecule < 0:
        raise ValueError(f"atp_per_molecule must be >= 0, got {atp_per_molecule}")
    if not provenance.strip():
        raise ValueError(
            "provenance is required: a cost without one is worse here than a refusal. "
            "Say ASSERTED or PLACEHOLDER plainly if that is what it is")

    cargo = model.metabolites.get_by_id(metabolite_id)
    stem = metabolite_id
    for suffix in (f"_{cargo.compartment}", f"[{cargo.compartment}]"):
        stem = stem[: -len(suffix)] if stem.endswith(suffix) else stem
    exported_id = exported_id or f"{stem}_e"
    if model.metabolites.has_id(exported_id):
        raise ValueError(f"{exported_id!r} already exists; export is already installed")
    exported = cobra.Metabolite(exported_id, name=f"{cargo.name} [extracellular]",
                                formula=cargo.formula, charge=cargo.charge, compartment="e")
    model.add_metabolites([exported])

    leader = leader_peptide_atp_equivalents(leader_residues)
    cost = float(atp_per_molecule) + leader
    # Not "EX_": that prefix means a boundary exchange, and this is a transport reaction
    # standing beside the real exchange the block below adds.
    reaction = cobra.Reaction(reaction_id or f"SEC_{stem}",
                              name=f"ATP-coupled export of {cargo.name} ({cost:g} ATP)")
    reaction.bounds = (0.0, 1000.0)
    reaction.subsystem = "heterologous export"
    stoichiometry = {cargo: -1.0, exported: 1.0}
    if cost:
        for base, coefficient in ((_ATP_C, -cost), (_H2O_C, -cost),
                                  (_ADP_C, cost), (_PI_C, cost), (_PROTON_C, cost)):
            stoichiometry[model.metabolites.get_by_id(_resolve(model, base))] = coefficient
    reaction.add_metabolites(stoichiometry)
    model.add_reactions([reaction])
    reaction.notes["secretion_cost_provenance"] = provenance
    reaction.notes["secretion_cost_atp"] = {"transport": float(atp_per_molecule),
                                            "leader_translation": leader,
                                            "leader_residues": int(leader_residues)}
    reaction.notes["secretion_cost_not_charged"] = [what for what, _ in SECRETION_NOT_CHARGED]
    if exchange:
        # Secretion only. cobra defaults an exchange to (-1000, 1000), which would let a
        # heterologous product be taken up from a medium that never contained it.
        model.add_boundary(exported, type="exchange", lb=0.0)
    return reaction


# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PrecursorDrain:
    """What the host itself pulls through a precursor node, and whether that can be asked.

    Args:
        metabolite: The node, as resolved in this model.
        floor: Minimum total consumption the network is forced through it, mmol/gDCW/h.
        glucose_uptake: The uptake bound it was measured under. Load-bearing -- the FPP
            floor moves 10.3x between the shipped default and glucose 10.
        consumers: Reactions that can consume the node, reversible ones included.
        vacuous: True when the node has fewer than two consumers, so no competition exists
            to measure and the zero is topology rather than biology.
        why: One sentence a reader can act on.
    """

    metabolite: str
    floor: float
    glucose_uptake: float
    consumers: tuple[str, ...]
    vacuous: bool
    why: str

    def summary(self) -> str:
        verdict = "VACUOUS" if self.vacuous else f"{self.floor:.4g} mmol/gDCW/h"
        return f"{self.metabolite}: {verdict} at glucose <= {self.glucose_uptake:g} -- {self.why}"


def native_drain(
    model: cobra.Model,
    precursor_metabolite: str,
    glucose_uptake: float,
    growth_fraction: float = 0.9,
    biomass_reaction: str = "r_2111",
) -> PrecursorDrain:
    """The precursor budget `fba/audit.py` advertises, in a form that is not vacuous.

    Four changes against :func:`ystwin.fba.audit.precursor_floor`, each measured rather
    than argued. ``glucose_uptake`` is REQUIRED, because the floor moves an order of
    magnitude with it; consumption is weighted by the stoichiometric coefficient instead
    of read off the raw flux, which ALONE doubles the FPP floor because squalene synthase
    takes two; reversible consumers are counted rather than dropped by a
    ``lower_bound >= 0`` filter, which on FPP adds 8.0e-7 and on acetyl-CoA adds exactly
    nothing; and a node with fewer than two consumers is reported as vacuous instead of as
    a zero that reads like "no native demand".

    Raises:
        KeyError: if the metabolite is not in the model, naming it.
        ValueError: on a non-positive uptake, or if the model cannot grow.
    """
    if glucose_uptake <= 0:
        raise ValueError(
            f"glucose_uptake must be positive, got {glucose_uptake}. It is required rather "
            "than defaulted because this floor is a property of it")

    resolved = precursor_metabolite
    if not model.metabolites.has_id(resolved):
        alternatives = [m.id for m in model.metabolites
                        if m.id.split("[")[0] == precursor_metabolite]
        if not alternatives:
            raise KeyError(
                f"no metabolite {precursor_metabolite!r} in this model; the precursor node "
                "named by the pathway spec has to exist before its budget can be read")
        resolved = alternatives[0]

    metabolite = model.metabolites.get_by_id(resolved)
    # A reversible reaction consumes the node in one of its two directions whichever way it
    # is written, so both signs qualify and the filter is on reversibility, not on the sign.
    consumers = {r.id: r.metabolites[metabolite] for r in metabolite.reactions
                 if r.metabolites[metabolite] < 0 or r.lower_bound < 0}
    if len(consumers) < 2:
        return PrecursorDrain(
            metabolite=resolved, floor=0.0, glucose_uptake=float(glucose_uptake),
            consumers=tuple(sorted(consumers)), vacuous=True,
            why=(f"{len(consumers)} reaction(s) can consume this node, so nothing competes "
                 "with the heterologous route for it and the floor is zero by topology. "
                 "Ask at the nearest upstream branch point instead"))

    with model as m:
        pin_solver(m)
        cap_uptake(m, GLUCOSE_EXCHANGE, glucose_uptake)
        m.objective = biomass_reaction
        max_growth = growth_or_none(m)
        if max_growth is None or max_growth <= 0:
            raise ValueError("model cannot grow under this uptake, so no floor can be read")
        m.reactions.get_by_id(biomass_reaction).lower_bound = max_growth * growth_fraction
        fva = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=sorted(consumers), fraction_of_optimum=0.0,
            processes=FVA_PROCESSES)

    # Net consumption is -coef * v: increasing in v when the node is a reactant, decreasing
    # when it is a product, so the guaranteed drain sits at opposite ends of the interval.
    floor = 0.0
    for rid, coefficient in consumers.items():
        extreme = fva.loc[rid, "minimum" if coefficient < 0 else "maximum"]
        floor += max(0.0, -coefficient * float(extreme))

    return PrecursorDrain(
        metabolite=resolved, floor=floor, glucose_uptake=float(glucose_uptake),
        consumers=tuple(sorted(consumers)), vacuous=False,
        why=(f"{len(consumers)} native consumers compete for this node; open a medium that "
             "supplies what they make and the floor falls, so this is context, not a bound"))
