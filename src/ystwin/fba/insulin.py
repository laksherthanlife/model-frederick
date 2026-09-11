"""Install a human insulin precursor demand on Yeast9, as a PROTEIN and not a metabolite.

Every product this repository has modelled so far -- beta-carotene, PHB, glycogen -- is a
small molecule made by a linear chain of enzymes, and `fba/carotenoid.py` is shaped for
exactly that: name the precursor, write the enzymes, drain the end. A protein is a
different kind of thing. It is not made by an enzyme chain at all; it is polymerised on a
ribosome from charged tRNAs at a cost paid in ATP and GTP, and the sequence -- not a
stoichiometric coefficient -- is what fixes the reaction.

So this module is deliberately NOT a copy of the carotenoid one below the surface. It keeps
that file's contract (return a copy, resolve ids rather than assume them, refuse a second
install, one demand at the end) and changes the thing underneath: the stoichiometry is
COMPUTED from a sequence, against the aminoacyl-tRNA pool that Yeast9's own protein
pseudoreaction `r_4047` uses.

WHAT THIS PRODUCES IS A CEILING, AND A CEILING IS NOT A PREDICTION. That is this
repository's standing finding about FBA and it is not softened here. On beta-carotene the
ceiling is non-binding; on glycogen it stands 2,967x above the measured content. There is
NO insulin expression measurement, NO insulin titre and NO insulin plate in this
repository, so nothing here is calibrated and nothing here is predicted. What an FBA answer
can honestly say is "not more than this", and the number below is that and nothing else.

THERE IS DELIBERATELY NO `data/pathways/insulin.toml`, and the reason is the spec schema's
own, not a shortcut. `pathway/spec.py` declares its scope in one line -- "INTRACELLULAR,
NON-SECRETED products, whose only outlet besides the next step is dilution by growth" --
and refuses `fate = "secreted"` at load time. A recombinant insulin precursor in yeast is
secreted BY DESIGN: it is expressed behind the MFalpha1 leader precisely so it leaves the
cell, which is what makes downstream purification possible. Its balance is therefore
`d[X]/dt = v_in - v_secretion - mu*[X]`, and `solve_pathway`'s terminal line `X = v_in/mu`
omits the term that does most of the work. Three further reasons, each on its own
sufficient:

  * a `PathwaySpec` is a chain of METABOLITE nodes with a `molar_mass_g_per_mol` and a
    `rate_law` each. Translation has one node and no intermediates to declare;
  * `pathway/flux.py`'s law is `flux = alpha * expression` on an ENTRY ENZYME. The entry
    "enzyme" here is the ribosome, which is not a construct gene and whose dosage nobody
    varied;
  * `measurable` would be `false` on the only node, so `Calibratability` would report the
    spec unfittable anyway -- the same verdict `phb.toml` already carries.

A file that cannot load, describing a fit nobody can make, would be a claim rather than a
record. The refusal is the deliverable.
"""

from __future__ import annotations

import cobra

from .secretion import ER_PEROXIDE_ID, add_oxidative_folding

__all__ = [
    "DISULFIDE_BONDS",
    "GTP_PER_RESIDUE",
    "PATHWAY_REACTION_IDS",
    "PRODUCT_DEMAND_ID",
    "PROINSULIN_SEQUENCE",
    "TRANSLOCATION_ID",
    "add_insulin_precursor_pathway",
    "residue_composition",
]

TRANSLATION_ID = "INSPRE_TRANSLATION"
TRANSLOCATION_ID = "INSPRE_TRANSLOCATION"
FOLDING_ID = "INSPRE_FOLDING"
PRODUCT_DEMAND_ID = "DM_insulin_precursor_er"
PATHWAY_REACTION_IDS = (TRANSLATION_ID, TRANSLOCATION_ID, FOLDING_ID, PRODUCT_DEMAND_ID)

REDUCED_METABOLITE_ID = "insulin_precursor_red_c"
TRANSLOCATED_METABOLITE_ID = "insulin_precursor_red_er"
FOLDED_METABOLITE_ID = "insulin_precursor_er"

PROTEIN_PSEUDOREACTION_ID = "r_4047"
"""Yeast9's protein pseudoreaction, and the authority this module reads its ids out of.

Read from the vendored SBML this session rather than remembered. It consumes twenty
AMINOACYL-tRNAs and releases twenty free tRNAs plus one unit of the pseudo-metabolite
`s_3717` "protein", e.g.

    0.527012401964609 s_0404 (Ala-tRNA(Ala)) + ... --> 0.527012401964609 s_1582 (tRNA(Ala))
    + ... + s_3717

Two things about it decide the shape of everything below.

**It is written against charged tRNAs, not free amino acids.** That is why this module is
too. Writing insulin off the free amino acid pool would silently skip the aminoacyl-tRNA
synthetases and with them the 2 ATP-equivalents per residue they cost -- `r_0157` is
`ATP + L-alanine + tRNA(Ala) --> Ala-tRNA(Ala) + AMP + diphosphate`, i.e. ATP to AMP, two
high-energy phosphate bonds, and it is mass- and charge-balanced in the model as written.
Going through the tRNA pool means the host GEM charges that cost itself and this module
does not have to assert it.

**It carries no GTP, and that is not because translation is free.** yeast-GEM puts the
elongation energy in the growth-associated maintenance of the BIOMASS pseudoreaction
`r_4041` -- `55.3 ATP + 55.3 H2O --> 55.3 ADP + 55.3 H+ + 55.3 phosphate` -- aggregated
over every macromolecule at once. A heterologous protein does not pass through `r_4041`,
so it would pay no elongation cost at all unless one is written here. See
:data:`GTP_PER_RESIDUE`.
"""

PROINSULIN_SEQUENCE = (
    "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"   # B chain, P01308 25-54
    "RR"                               # dibasic Kex2/PC1 site, 55-56
    "EAEDLQVGQVELGGGPGAGSLQPLALEGSLQ"  # C peptide, 57-87
    "KR"                               # dibasic Kex2/PC2 site, 88-89
    "GIVEQCCTSICSLYQLENYCN"            # A chain, 90-110
)
"""Human proinsulin: UniProt P01308 (INS_HUMAN), residues 25-110. 86 residues.

FETCHED, NOT REMEMBERED. `https://rest.uniprot.org/uniprotkb/P01308.fasta` and
`.../P01308.txt` were retrieved this session; the entry is 110 aa and its feature table
reads

    FT   SIGNAL          1..24
    FT   PEPTIDE         25..54     /note="Insulin B chain"
    FT   PROPEP          57..87     /note="C peptide"
    FT   PEPTIDE         90..110    /note="Insulin A chain"
    FT   DISULFID        31..96     /note="Interchain (between B and A chains)"
    FT   DISULFID        43..109    /note="Interchain (between B and A chains)"
    FT   DISULFID        95..100

so residues 25-110 are the signal-cleaved PROINSULIN and the two residues at 55-56 and
88-89 that the feature table leaves unassigned are the dibasic processing sites, RR and KR,
which is what this string's line breaks record. The four segments concatenate back to
`P01308[24:110]` character for character, and
`tests/test_insulin_pathway.py` asserts that against the fetched FASTA's own residues where
the file is present.

WHICH PRECURSOR THIS IS, AND WHICH IT IS NOT. Industrial yeast processes generally express
a SHORTER single-chain precursor -- a mini-proinsulin in which the 31-residue C peptide is
replaced by a short linker -- rather than full proinsulin. That construct is NOT modelled
here, and the reason is that no sequence for one could be established from an accessible
source this session: the searches that returned it returned patents and closed-access
papers. Modelling proinsulin is the version that is fully sourced, and it is the
CONSERVATIVE direction for a ceiling in the only sense that matters -- a shorter chain
costs less per molecule, so the ceiling on a mini-proinsulin is HIGHER than the number this
module reports, not lower. Anyone wanting that number needs a cited sequence, and then this
module computes it unchanged: the stoichiometry is a function of the string.

THE INITIATOR METHIONINE AND THE SECRETORY LEADER ARE NOT HERE, and their absence is a
real omission rather than a modelling convention. What a ribosome actually translates is
`MFalpha1-leader + spacer + precursor`, tens of residues longer, all of them paid for and
then proteolytically discarded. Their cost is omitted, which again biases the ceiling
UPWARD. Worth noticing on its own account: proinsulin contains no methionine and no
tryptophan at all, so two of `r_4047`'s twenty aminoacyl-tRNA nodes carry a coefficient of
exactly zero in the reaction this module writes -- and a translated construct would have at
least the initiator Met.
"""

DISULFIDE_BONDS = 3
"""Three, from the same P01308 feature table quoted above: 31-96, 43-109 and 95-100.

Six of the 86 residues are cysteine, and all six are paired, which the composition confirms
independently: `residue_composition(PROINSULIN_SEQUENCE)["C"] == 6 == 2 * DISULFIDE_BONDS`.
"""

GTP_PER_RESIDUE = 2.0
"""GTP hydrolysed per residue incorporated. THE ONE LOAD-BEARING ASSUMPTION IN THIS FILE.

Martin WF (2025), "ATP requirements for growth reveal the bioenergetic impact of
mitochondrial symbiosis", Biochim Biophys Acta Bioenerg 149564, doi
10.1016/j.bbabio.2025.149564, PMID 40562331, PMC7617979, CC-BY -- fetched as Europe PMC
full-text XML this session, not recalled. Verbatim:

    "The pure cost of synthesizing peptide bonds is 4 ATP each: 2 from PPi formation at
    aminoacyl tRNA synthesis, which renders the reaction irreversible, and 1 GTP each for
    the two elongation factors."

and, of the same figure, "The cost assumed for protein synthesis at the ribosome is
uncontested, 4 ATP per peptide bond", tracing it to Stouthamer 1973 (PMID 4148026), which
is closed access and was therefore NOT used as the source.

The split is what makes the number usable here. Of the four, the first two are the
synthetase step -- ATP to AMP plus PPi -- which the host GEM already charges through its own
twenty aminoacyl-tRNA synthetase reactions the moment this module draws on the charged
pool. Only the remaining two, the elongation factor GTPs (eEF1A delivery and eEF2
translocation in yeast), are missing, and they are what this constant adds.

**Why adding them is a correction and not a surcharge.** Native yeast protein pays the
elongation cost through the 55.3 ATP of growth-associated maintenance in `r_4041`.
Heterologous insulin never touches `r_4041`. Charging nothing would make a heterologous
protein strictly cheaper than the host's own, which is wrong in a direction that inflates
the ceiling. It is still an assumption, so it is a parameter: pass `gtp_per_residue=0.0`
to `add_insulin_precursor_pathway` and the model reproduces yeast-GEM's own convention
exactly, with translation energy accounted nowhere. Both numbers are reported rather than
one, because how much of the answer rests on this is the reader's business.

PER RESIDUE, NOT PER PEPTIDE BOND, and the two differ. An 86-residue chain has 85 peptide
bonds; the first residue is placed by initiation, which costs its own GTP (eIF2) and ATP,
neither of which is modelled. Martin 2025 itself moves between the two -- "4.78 mmol of
peptide bonds or 4.78 mmol of amino acids" -- because at proteome scale the difference is
invisible. Here it is 1.2%, smaller than any other uncertainty in this file, and per
residue is the form that does not need an initiation term bolted on to be self-consistent.
"""

_HILL_ORDER_FIRST = ("C", "H")

#: IUPAC three-letter to one-letter codes. Not a judgement and not derived from anything in
#: this repository; it is the naming convention `r_4047`'s metabolites are written in
#: ("Ala-tRNA(Ala)") and the one a FASTA sequence is written in, and this dictionary is the
#: only thing joining them. Twenty entries, checked for completeness at install time.
_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Cytosolic cofactor ids, Yeast9 v9.0.2, read out of data/gem/yeast-GEM.xml.gz this session.
# Resolved rather than trusted at install time -- see `_resolve`.
_GTP_C, _GDP_C, _PI_C = "s_0785", "s_0739", "s_1322"
_H2O_C, _PROTON_C = "s_0803", "s_0794"
#: ER oxygen. The folding step consumes it THERE and puts its peroxide into the ER pool
#: that `fba/secretion.py` installs, so no cytosolic O2 or H2O2 id belongs here. Both were
#: defined and resolved until 2026-08-31, left over from the version of this step that ran
#: in the cytosol -- an unused requirement that would have refused a host model missing a
#: species this pathway never touches.
_O2_ER = "s_1276"


def _resolve(model: cobra.Model, base_id: str, compartment: str = "c") -> str:
    """Map a bare Yeast9 metabolite id onto whatever convention this model uses.

    Same contract as `carotenoid.py::_resolve`, and deliberately a second copy rather than
    an import of a sibling module's private name: which ids a product needs is that
    product's business, and coupling two unrelated pathways through an underscore-prefixed
    helper neither of them owns buys nothing. Yeast9 SBML uses ``s_0785``; GECKO exports
    use ``s_0785[c]``.
    """
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(
        f"required metabolite {base_id!r} not found in {model.id!r} "
        f"(tried {base_id}, {base_id}[{compartment}], {base_id}_{compartment})")


def residue_composition(sequence: str) -> dict[str, int]:
    """Count residues by one-letter code, refusing anything that is not one of the twenty.

    Args:
        sequence: One-letter amino acid sequence, upper case.

    Returns:
        Mapping from one-letter code to count. Only residues that occur appear.

    Raises:
        ValueError: naming every character that is not one of the twenty canonical
            residues. Selenocysteine (U), pyrrolysine (O) and the ambiguity codes (B, Z,
            X) are refused rather than mapped: `r_4047` carries no tRNA for any of them,
            so there is nothing to write the reaction against, and a guess would be an
            invented residue in a product formula.
    """
    known = set(_THREE_TO_ONE.values())
    unknown = sorted({c for c in sequence if c not in known})
    if unknown:
        raise ValueError(
            f"sequence contains {len(unknown)} residue code(s) with no tRNA in the model: "
            f"{unknown}. Only the twenty canonical residues can be written against "
            f"{PROTEIN_PSEUDOREACTION_ID}.")
    return {code: sequence.count(code) for code in sorted(known) if code in sequence}


def _translation_species(model: cobra.Model, pseudoreaction_id: str):
    """Read the aminoacyl-tRNA / free-tRNA pairs out of the host model's own translation.

    The alternative -- a hardcoded table of twenty `s_` ids -- is exactly the kind of
    remembered constant this repository does not accept: it would be right for one release
    of one GEM and silently wrong for the next, and nothing would say so. `r_4047` IS the
    model's statement of which species translation runs on, so it is what gets read.

    Returns:
        ``(charged, free)``, both keyed by one-letter residue code.

    Raises:
        KeyError: if the pseudoreaction is absent, naming it.
        ValueError: if it does not resolve to exactly the twenty canonical residues, with
            what was found. A GEM whose protein pseudoreaction is shaped differently is a
            GEM this module has not been checked against, and proceeding on a partial parse
            would put an unbalanced reaction into a model that still solves.
    """
    if not model.reactions.has_id(pseudoreaction_id):
        raise KeyError(
            f"{pseudoreaction_id!r} (the protein pseudoreaction) is not in {model.id!r}; "
            f"this module reads the aminoacyl-tRNA pool out of it and cannot proceed "
            f"without it")
    reaction = model.reactions.get_by_id(pseudoreaction_id)
    charged: dict[str, cobra.Metabolite] = {}
    free: dict[str, cobra.Metabolite] = {}
    for metabolite, coefficient in reaction.metabolites.items():
        name = metabolite.name
        if coefficient < 0 and "-tRNA(" in name:
            code = _THREE_TO_ONE.get(name.split("-tRNA(")[0])
            if code:
                charged[code] = metabolite
        elif coefficient > 0 and name.startswith("tRNA(") and name.endswith(")"):
            code = _THREE_TO_ONE.get(name[len("tRNA("):-1])
            if code:
                free[code] = metabolite
    expected = set(_THREE_TO_ONE.values())
    if set(charged) != expected or set(free) != expected:
        raise ValueError(
            f"{pseudoreaction_id!r} in {model.id!r} did not yield all twenty "
            f"aminoacyl-tRNA/tRNA pairs: charged={sorted(charged)}, free={sorted(free)}")
    return charged, free


def _hill_formula(elements: dict[str, float]) -> str:
    """Render an element count as a Hill-order formula string, dropping zeros.

    Counts are integers by construction here -- a residue count times an integer element
    count -- so they are rendered as integers; a non-integral count would be a bug and is
    left to render as a float rather than being rounded into invisibility.
    """
    def render(symbol: str) -> str:
        count = elements[symbol]
        count = int(count) if float(count).is_integer() else count
        return symbol if count == 1 else f"{symbol}{count}"

    ordered = [s for s in _HILL_ORDER_FIRST if elements.get(s)]
    ordered += sorted(s for s in elements if s not in _HILL_ORDER_FIRST and elements[s])
    return "".join(render(s) for s in ordered)


def _chain_composition(counts: dict[str, int], charged, free):
    """Elemental formula and charge of the reduced polypeptide, in the host's own convention.

    Derived, never asserted. Yeast9 writes a charged tRNA as the free tRNA plus the
    aminoacyl group -- ``Ala-tRNA(Ala)`` is ``C3H7NOR`` at charge +1 against ``tRNA(Ala)``
    at ``RH`` and charge 0 -- so the difference between the two IS the residue as this
    model accounts it, and summing those differences is what makes the reaction below
    balanced by construction rather than by a formula somebody typed.

    The chemistry that fixes the two remaining terms:

      * ONE water is consumed, at termination. Peptidyl-tRNA hydrolysis is what releases
        the finished chain from the last tRNA;
      * ``L`` protons are released, one per residue. Each of the ``L - 1`` peptide bonds
        needs the incoming alpha-amino group deprotonated to attack, and the C-terminal
        carboxylate gives up the last one at pH 7. The N-terminal ammonium keeps its.

    Both are checkable rather than asserted, and the test file checks them the hard way:
    the formula this returns must equal the sum of the model's own FREE amino acid formulas
    minus ``L - 1`` waters -- the textbook polypeptide formula, reached by a completely
    different route through the same GEM.

    Returns:
        ``(elements, charge)`` for the REDUCED chain, cysteines as free thiols.
    """
    length = sum(counts.values())
    elements: dict[str, float] = {}
    charge = 0.0
    for code, n in counts.items():
        for symbol, count in charged[code].elements.items():
            elements[symbol] = elements.get(symbol, 0.0) + n * count
        for symbol, count in free[code].elements.items():
            elements[symbol] = elements.get(symbol, 0.0) - n * count
        charge += n * charged[code].charge
    elements["H"] = elements.get("H", 0.0) + 2.0 - length   # + H2O, - L protons
    elements["O"] = elements.get("O", 0.0) + 1.0
    charge -= length
    return {s: v for s, v in elements.items() if v}, charge


def _add_reaction(model, rid, name, stoichiometry, gene_rule, bounds=(0.0, 1000.0),
                  subsystem=""):
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.subsystem = subsystem
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gene_rule
    return rxn


def add_insulin_precursor_pathway(
    model: cobra.Model,
    sequence: str = PROINSULIN_SEQUENCE,
    disulfide_bonds: int = DISULFIDE_BONDS,
    gtp_per_residue: float = GTP_PER_RESIDUE,
    protein_pseudoreaction_id: str = PROTEIN_PSEUDOREACTION_ID,
) -> cobra.Model:
    """Return a COPY of ``model`` able to make and drain a human insulin precursor.

    The input model is not modified. Three reactions are added:

    ``INSPRE_TRANSLATION``
        the chain, off the aminoacyl-tRNA pool, plus ``gtp_per_residue`` GTP per residue.
    ``INSPRE_FOLDING``
        oxidative folding: ``disulfide_bonds`` disulfides formed, electrons to O2, one
        hydrogen peroxide out per bond. See below -- this is a mechanism choice and it is
        stated as one.
    ``DM_insulin_precursor_c``
        a demand, not an exchange, for the reason given below.

    Args:
        model: Host model to extend.
        sequence: One-letter residue sequence of the precursor. Defaults to
            :data:`PROINSULIN_SEQUENCE`; anything else needs its own source.
        disulfide_bonds: How many disulfides the folded product carries. Must not exceed
            half the cysteine count.
        gtp_per_residue: See :data:`GTP_PER_RESIDUE`. Pass ``0.0`` for yeast-GEM's own
            convention, in which translation energy is accounted nowhere outside biomass.
        protein_pseudoreaction_id: Where the aminoacyl-tRNA ids are read from.

    Returns:
        A new model carrying the three reactions and two metabolites.

    Raises:
        ValueError: if the pathway is already installed, if ``disulfide_bonds`` exceeds the
            cysteines available, or if the host's protein pseudoreaction is not shaped as
            expected.
        KeyError: naming any metabolite or reaction the host model does not carry.

    OXIDATIVE FOLDING IS WRITTEN, AND WHAT IT COSTS IS NOT FREE. Insulin's three disulfides
    are not a structural footnote; forming each one removes two hydrogens from the chain
    and those electrons have to go somewhere. In yeast that somewhere is Ero1p, and Gross
    E, Sevier CS, Heldman N, Vitu E, Bentzur M, Kaiser CA, Thorpe C, Fass D (2006),
    "Generating disulfides enzymatically: reaction products and electron acceptors of the
    endoplasmic reticulum thiol oxidase Ero1p", PNAS 103:299-304, PMID 16407158,
    PMC1326156, doi 10.1073/pnas.0506448103, measured the product on RECOMBINANT YEAST
    Ero1p and reported that "under aerobic conditions, reduction of molecular oxygen by
    Ero1p yielded stoichiometric hydrogen peroxide". One O2 in, one H2O2 out, per
    disulfide, and the two hydrogens balance exactly.

    That is the same shape as the correction `carotenoid.py` records for CrtI, and it lands
    in the same place: **three hydrogen peroxide per insulin precursor**, made by the
    oxidative stress this repository's HyPer7 and Yap1 sensors exist to read. That is a
    testable consequence of expressing a disulfide-rich protein and it is not free.

    Two caveats on it, both structural.

    The same paper's headline is that Ero1p "can transfer electrons to a variety of small
    and macromolecular electron acceptors in addition to molecular oxygen", and it exists
    because the anaerobic fate was unknown. O2 is the aerobic route and it is the one
    written; **this pathway therefore cannot run anaerobically as modelled**, which is a
    property of the reaction as written and not a measurement of the cell.

    THE COMPARTMENT WAS WRONG AND IS NOW RIGHT. Oxidative folding happens in the ER. Yeast9
    v9.0.2 carries ER oxygen (`s_1276`) but had **no ER hydrogen peroxide species at all** --
    the only H2O2 pools are cytosol `s_0837`, mitochondrion, nucleus and peroxisome -- and no
    ER-to-cytosol peroxide transport. So this step was first written in the cytosol with the
    compartment declared wrong, on the argument that inventing a species and a transporter
    would be two fabrications to gain nothing.

    `fba/secretion.py` then supplied both properly, with the export to the compartment
    carrying catalase and the two-hydrogens-per-bond formula correction, and this module now
    installs it and folds in the ER. **The flux ceiling is unchanged at 0.009898
    mmol/gDCW/h**, which is what the cytosolic version predicted would happen -- neither O2
    nor H2O2 is compartment-limited in the solution. What changed is that the reaction is now
    where the chemistry is, and `INSPRE_TRANSLOCATION` moves the chain there.

    That translocation carries CARBON ONLY. SRP targeting, the Sec61 translocon and signal
    peptidase are machinery Yeast9 does not model, so their ATP cost cannot be charged and is
    not pretended to be.

    THE PRODUCT LEAVES BY A DEMAND, NOT AN EXCHANGE, AND THAT IS THE LARGEST OMISSION HERE.
    A real insulin precursor is SECRETED, behind the MFalpha1 leader. Everything that makes
    that happen -- SRP targeting, translocation, signal peptidase, PDI-assisted folding,
    Kex2 and Ste13 processing, COPII vesicles, the Golgi -- is machinery Yeast9 does not
    model, at an ATP cost it therefore cannot charge. A demand reaction says "carbon and
    energy leave the balance here" and claims nothing about how. Adding an extracellular
    species and a costless transporter would look more like secretion and assert more:
    it would name a transport step and give it a price of zero. The ceiling is the same
    number either way, and the omitted secretory cost biases it UPWARD.
    """
    if model.metabolites.has_id(FOLDED_METABOLITE_ID):
        raise ValueError("insulin precursor pathway is already installed on this model")

    counts = residue_composition(sequence)
    length = sum(counts.values())
    cysteines = counts.get("C", 0)
    if 2 * disulfide_bonds > cysteines:
        raise ValueError(
            f"{disulfide_bonds} disulfide bond(s) need {2 * disulfide_bonds} cysteines and "
            f"the sequence has {cysteines}")

    out = model.copy()
    charged, free = _translation_species(out, protein_pseudoreaction_id)
    gtp = _resolve(out, _GTP_C)
    gdp = _resolve(out, _GDP_C)
    phosphate = _resolve(out, _PI_C)
    water = _resolve(out, _H2O_C)
    proton = _resolve(out, _PROTON_C)

    elements, charge = _chain_composition(counts, charged, free)
    reduced_formula = _hill_formula(elements)
    folded_elements = dict(elements)
    folded_elements["H"] = folded_elements["H"] - 2.0 * disulfide_bonds
    # The ER peroxide pool and its export to the cytosol, from `fba/secretion.py`. Installed
    # here rather than assumed, and idempotent-by-refusal: a model that already carries it
    # says so, which is how a second secreted product reuses one pool instead of adding two.
    if not out.metabolites.has_id(ER_PEROXIDE_ID):
        out = add_oxidative_folding(out)

    out.add_metabolites([
        cobra.Metabolite(
            REDUCED_METABOLITE_ID,
            name=f"human insulin precursor, reduced ({length} aa)",
            formula=reduced_formula, charge=charge, compartment="c"),
        cobra.Metabolite(
            FOLDED_METABOLITE_ID,
            name=f"human insulin precursor ({length} aa, {disulfide_bonds} disulfides)",
            formula=_hill_formula(folded_elements), charge=charge, compartment="er"),
        cobra.Metabolite(
            TRANSLOCATED_METABOLITE_ID,
            name=f"human insulin precursor, reduced, in the ER ({length} aa)",
            formula=reduced_formula, charge=charge, compartment="er"),
    ])
    # UniProt is the identity of the PRODUCT, the way `kegg.compound` is the identity of a
    # carotenoid in `carotenoid.py`. Nothing in this repository resolves a thermodynamic
    # energy through it -- there is no formation energy for a polypeptide in any vendored
    # table, and `bridge/thermodynamic.py` will return None for both of these, which is the
    # correct answer and not a gap to paper over. The annotation is here so that a reader
    # holding this model can check what the formula is a formula OF.
    for mid in (REDUCED_METABOLITE_ID, FOLDED_METABOLITE_ID):
        out.metabolites.get_by_id(mid).annotation.update(
            {"uniprot": "P01308", "sbo": "SBO:0000252"})

    # Translation. The aminoacyl-tRNA terms make the host's own synthetases pay the 2
    # ATP-equivalents per residue; the GTP terms are the elongation factors -- see
    # GTP_PER_RESIDUE for why they are here at all and why they are a parameter.
    #   1 H2O for the termination hydrolysis, L protons from the peptide bonds,
    #   plus one H2O and one proton per GTP hydrolysed.
    gtp_total = gtp_per_residue * length
    stoichiometry: dict[str, float] = {}
    for code, n in counts.items():
        stoichiometry[charged[code].id] = -float(n)
        stoichiometry[free[code].id] = float(n)
    stoichiometry[water] = -(1.0 + gtp_total)
    stoichiometry[proton] = float(length) + gtp_total
    if gtp_total:
        stoichiometry[gtp] = -gtp_total
        stoichiometry[gdp] = gtp_total
        stoichiometry[phosphate] = gtp_total
    stoichiometry[REDUCED_METABOLITE_ID] = 1.0
    _add_reaction(
        out, TRANSLATION_ID,
        f"insulin precursor translation ({length} aa)",
        stoichiometry,
        # No gene rule. `carotenoid.py` names crtE, crtYB and crtI because those are the
        # construct's own genes and a reaction there IS an enzyme. This reaction is the
        # ribosome plus two elongation factors plus twenty synthetases -- roughly two
        # hundred gene products, none of them heterologous. Naming one would be a fiction;
        # naming the construct's INS coding sequence would be a fiction of a different kind,
        # since the coding sequence is the template and not the catalyst. Left empty, which
        # is what `spec.py` already uses to mean "no single enzyme owns this step".
        "",
        subsystem="insulin precursor production",
    )
    # Translocation into the ER. Carbon only: SRP targeting, the Sec61 translocon and signal
    # peptidase are machinery Yeast9 does not model, so their ATP cost cannot be charged and
    # is not pretended to be. What this reaction does is put the chain in the compartment
    # where folding actually happens, which is the difference between a declared-wrong
    # compartment and a right one.
    _add_reaction(
        out, TRANSLOCATION_ID,
        f"insulin precursor translocation into the ER ({length} aa)",
        {REDUCED_METABOLITE_ID: -1.0, TRANSLOCATED_METABOLITE_ID: 1.0},
        "",
        subsystem="insulin precursor production",
    )
    # Oxidative folding, Ero1p route: 2 RSH + O2 -> RSSR + H2O2 per disulfide, IN THE ER.
    #
    # This was written in the cytosol and declared wrong, because Yeast9 v9.0.2 has no ER
    # peroxide species and no ER-to-cytosol peroxide transport, and inventing both to gain
    # nothing would have been two fabrications. `fba/secretion.py` supplies them properly --
    # the pool, the export to the compartment carrying catalase, and the two-hydrogens-per-
    # bond formula correction -- so the compartment is now right rather than declared wrong.
    _add_reaction(
        out, FOLDING_ID,
        f"insulin precursor oxidative folding ({disulfide_bonds} disulfides, Ero1p/O2)",
        {TRANSLOCATED_METABOLITE_ID: -1.0, _O2_ER: -float(disulfide_bonds),
         FOLDED_METABOLITE_ID: 1.0, ER_PEROXIDE_ID: float(disulfide_bonds)},
        "",
        subsystem="insulin precursor production",
    )
    out.add_boundary(
        out.metabolites.get_by_id(FOLDED_METABOLITE_ID),
        type="demand", reaction_id=PRODUCT_DEMAND_ID,
    )
    return out
