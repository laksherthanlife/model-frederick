"""The oxidative cost of secretory folding, which yeast-GEM v9.0.2 does not carry.

**Why this module exists.** Every product this repository had modelled was a small molecule.
A SECRETED protein is different, and the difference is a cost the vendored model cannot
express: forming a disulfide bond in the endoplasmic reticulum consumes molecular oxygen and
produces hydrogen peroxide, and neither the reaction nor the peroxide is in the GEM.

Checked against yeast-GEM v9.0.2 rather than assumed, before any of this was written:

- Searching reaction names for *disulfide*, *thiol oxidation*, *Ero1* or *oxidoreductin*
  returns four hits and all four are cofactor recycling -- glutathione oxidoreductase twice
  (`r_0481`, `r_0482`), a glutathione disulfide exchange (`r_1806`) and a
  methionine/thioredoxin oxidoreductase (`r_4186`). **None forms a disulfide in a substrate
  protein.**
- There is **no `ERO1` and no `PDI1` gene** in the model.
- Hydrogen peroxide exists as `s_0837` (cytosol), `s_0838` (mitochondrion), `s_0839`
  (nucleus) and `s_0840` (peroxisome). The ER carries 135 metabolites and **peroxide is not
  one of them**.

So the burden that makes a disulfide-rich secreted protein expensive was absent, not merely
unparameterised. This adds it.

**The chemistry.** Protein disulfide isomerase oxidises a pair of cysteine thiols and is
itself re-oxidised by Ero1, a flavoprotein that passes the electrons to molecular oxygen. The
net per disulfide formed is

    2 R-SH + O2  ->  R-S-S-R + H2O2

which is the same two-electron flavin chemistry as CrtI in `fba/carotenoid.py`, and it is
written the same way here for the same reason: a flavoprotein reduces O2 at one site, so the
product is peroxide and not water.

**What is deliberately NOT modelled.** The protein thiols themselves. A cobra model has no
species for "reduced cysteine pair in a folding polypeptide", and inventing one would mean
inventing its formation energy, its pool size and its turnover. What is installed instead is
the NET cost, carried on the folding reaction of a named protein: *n* disulfides consume *n*
O2 and make *n* H2O2. That is an accounting of the oxygen and the peroxide, which is what the
flux balance can use and what this lab's Yap1 sensors can read. It is not a mechanism, and a
caller reading it as one would be over-reading it.

**Why the peroxide is exported to the cytosol.** It has to go somewhere or the ER species is a
dead end and the folding reaction cannot carry flux. Peroxide crosses membranes readily, and
the cytosol is where the model puts its catalase (`r_0255`) and where this project's
`NativeYap1` and `AlteredYap1` sensors report from. `docs/research/SECRETED_PROTEIN.md`
carries the argument.
"""
from __future__ import annotations

import cobra

__all__ = [
    "DISULFIDE_O2_PER_BOND",
    "ER_PEROXIDE_ID",
    "PEROXIDE_EXPORT_ID",
    "add_oxidative_folding",
    "folding_reaction",
]

#: ER-compartment hydrogen peroxide. Absent from yeast-GEM v9.0.2, which carries peroxide in
#: the cytosol, mitochondrion, nucleus and peroxisome only.
ER_PEROXIDE_ID = "h2o2_er"

#: ER -> cytosol peroxide, so the folding reaction is not a dead end and the peroxide reaches
#: the compartment carrying catalase and the sensors.
PEROXIDE_EXPORT_ID = "H2O2ter"

#: Molecular oxygen consumed, and peroxide produced, per disulfide bond formed. One, because
#: Ero1 is a two-electron flavoprotein oxidase and a disulfide is a two-electron oxidation.
DISULFIDE_O2_PER_BOND = 1.0

_O2_ER, _H2O2_C = "s_1276", "s_0837"


def add_oxidative_folding(model: cobra.Model) -> cobra.Model:
    """Return a COPY of ``model`` carrying ER peroxide and its export to the cytosol.

    Installed once per model and independent of any particular protein, so a second secreted
    product reuses it rather than adding a second peroxide pool. :func:`folding_reaction` is
    what attaches a specific protein's disulfide count to it.

    Args:
        model: A cobra model carrying an ``er`` compartment and cytosolic peroxide.

    Returns:
        A copy. The input is not modified, matching `fba/carotenoid.py`.

    Raises:
        ValueError: if it is already installed, or if the model lacks ER oxygen or cytosolic
            peroxide to connect to -- a silent no-op here would look like a closed gap.
    """
    if model.metabolites.has_id(ER_PEROXIDE_ID):
        raise ValueError("oxidative folding is already installed on this model")
    for required in (_O2_ER, _H2O2_C):
        if not model.metabolites.has_id(required):
            raise ValueError(
                f"this model has no {required!r}, so the folding cost cannot be connected "
                "to it. ER oxygen and cytosolic peroxide are both required; check the model "
                "is yeast-GEM v9.0.2 or supply the ids explicitly")

    out = model.copy()
    peroxide = cobra.Metabolite(ER_PEROXIDE_ID, name="hydrogen peroxide",
                                formula="H2O2", charge=0, compartment="er")
    out.add_metabolites([peroxide])

    export = cobra.Reaction(PEROXIDE_EXPORT_ID, name="hydrogen peroxide transport, ER to cytosol")
    export.bounds = (0.0, 1000.0)
    export.add_metabolites({peroxide: -1.0, out.metabolites.get_by_id(_H2O2_C): 1.0})
    export.subsystem = "secretory folding"
    out.add_reactions([export])
    return out


def _less_thiol_hydrogens(formula: str | None, disulfides: int) -> str | None:
    """The folded protein's formula: two hydrogens lighter per disulfide bond.

    ``2 R-SH + O2 -> R-S-S-R + H2O2`` moves two hydrogens off the protein and into the
    peroxide for every bond formed. Giving the folded species the SAME formula as the
    unfolded one leaves the reaction unbalanced by exactly ``2 * disulfides`` hydrogen, which
    is what `Reaction.check_mass_balance` reported the first time this was written and what
    made the bookkeeping honest.

    Returns ``None`` unchanged when the protein carries no formula, because a formula that
    cannot be read cannot be adjusted and inventing one would be worse than leaving the
    balance uncheckable.
    """
    import re

    if not formula or not disulfides:
        return formula
    # ACCUMULATED, not assigned. A comprehension keyed by element silently drops every
    # repeat: "C2H4O2H2" parses to C2O2 and all six hydrogens vanish, and "CH3COOH" reports
    # one hydrogen and raises a refusal whose message is a lie. Non-Hill spellings are
    # ordinary in a hand-written formula, and this function takes whatever a caller put on
    # the metabolite. Found by the 2026-08-31 code review.
    counts: dict[str, int] = {}
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element:
            counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    hydrogens = counts.get("H")
    if hydrogens is None or hydrogens < 2 * disulfides:
        raise ValueError(
            f"formula {formula!r} carries {hydrogens} hydrogen, fewer than the "
            f"{2 * disulfides} that {disulfides} disulfide bond(s) remove. Either the "
            "formula or the disulfide count is wrong")
    counts["H"] = hydrogens - 2 * disulfides
    return "".join(f"{element}{count if count != 1 else ''}"
                   for element, count in counts.items() if count)


def folding_reaction(model: cobra.Model, protein_id: str, disulfides: int, *,
                     reaction_id: str | None = None) -> cobra.Reaction:
    """Attach the oxidative cost of *disulfides* bonds to a protein already in the model.

    The reaction consumes the unfolded species and produces the folded one, paying
    :data:`DISULFIDE_O2_PER_BOND` oxygen and one peroxide per bond. It is added to the model
    in place, because the caller is mid-construction of a pathway.

    Args:
        model: A model :func:`add_oxidative_folding` has already been run on.
        protein_id: The metabolite whose folding is being paid for. It is consumed and
            re-produced as ``<protein_id>_folded``.
        disulfides: How many disulfide bonds the protein carries. **Zero is legal and means
            zero cost** -- not every secreted protein has any, and refusing zero would make
            the caller special-case it.
        reaction_id: Override the generated id.

    Raises:
        ValueError: for a negative disulfide count, or if :func:`add_oxidative_folding` has
            not been run.
    """
    if disulfides < 0:
        raise ValueError(f"disulfides must be >= 0, got {disulfides}")
    if not model.metabolites.has_id(ER_PEROXIDE_ID):
        raise ValueError(
            "run add_oxidative_folding(model) first; without ER peroxide the cost has "
            "nowhere to go and the reaction would silently be free")

    unfolded = model.metabolites.get_by_id(protein_id)
    folded = cobra.Metabolite(f"{protein_id}_folded", name=f"{unfolded.name}, folded",
                              formula=_less_thiol_hydrogens(unfolded.formula, disulfides),
                              charge=unfolded.charge, compartment=unfolded.compartment)
    model.add_metabolites([folded])

    reaction = cobra.Reaction(reaction_id or f"FOLD_{protein_id}",
                              name=f"oxidative folding of {unfolded.name} "
                                   f"({disulfides} disulfide bonds)")
    reaction.bounds = (0.0, 1000.0)
    reaction.subsystem = "secretory folding"
    cost = float(disulfides) * DISULFIDE_O2_PER_BOND
    stoichiometry = {unfolded: -1.0, folded: 1.0}
    if cost:
        stoichiometry[model.metabolites.get_by_id(_O2_ER)] = -cost
        stoichiometry[model.metabolites.get_by_id(ER_PEROXIDE_ID)] = cost
    reaction.add_metabolites(stoichiometry)
    model.add_reactions([reaction])
    return reaction
