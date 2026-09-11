"""Install the strain's crtE / crtYB(PSY) / crtI / crtYB(LCY) pathway into Yeast9 or GECKO.

Reached from the live chain by `fba/audit.py`, which asks this model whether a predicted
flux fits and what it costs. It does NOT supply the flux, and cannot: cap every reaction
here at the measured magnitude and FVA returns [0, cap] -- the ceiling moves and the floor
stays at exactly zero, because an upper bound cannot make a flux mandatory.

Still specific to one product. The chain around it is not: `pathway/spec.py` declares a
pathway as data and `pathway/solve.py` solves any of them. What has not been generalised is
INSTALLING arbitrary chemistry into a GEM, which is why this file still names crt genes.

crtYB is bifunctional, so two reactions carry the same gene. tHMG1 adds no chemistry and is a
bound change, not a reaction. Beta-carotene accumulates in membranes, so the model carries a
demand rather than an extracellular exchange.
"""

from __future__ import annotations

import cobra

__all__ = [
    "PATHWAY_REACTION_IDS",
    "PRODUCT_DEMAND_ID",
    "add_beta_carotene_pathway",
]

PRODUCT_DEMAND_ID = "DM_betacarotene_c"
PATHWAY_REACTION_IDS = ("CRTE", "CRTYB_PSY", "CRTI", "CRTYB_LCY", PRODUCT_DEMAND_ID)

# Yeast9 metabolite ids for the native precursors this pathway taps.
_FPP_C, _IPP_C, _GGPP_C, _PPI_C = "s_0190", "s_0943", "s_0189", "s_0633"
#: Cytosolic O2 and H2O2. CrtI is a flavin oxidase: its bound FAD is reoxidised by molecular
#: oxygen, two electrons at a time, so hydrogen peroxide is the product. See CRTI below.
#:
#: There are deliberately no FAD, FADH2 or proton ids here. Free flavin was the terminal
#: acceptor until 2026-08-30, and when CRTI was corrected the ids, the parameters that
#: defaulted to them and the `_resolve` calls that consumed them were all left behind --
#: so this function went on REFUSING any host model that did not carry FAD, FADH2 and a
#: cytosolic proton, for chemistry it no longer performs. Removed 2026-08-31.
_O2_C, _H2O2_C = "s_1275", "s_0837"

_NEW_METABOLITES = (
    # id, name, formula, charge, KEGG compound, ModelSEED compound
    # -- all neutral C40 hydrocarbons in the cytosol
    ("phytoene_c", "phytoene", "C40H64", 0, "C05413", "cpd03205"),
    ("lycopene_c", "lycopene", "C40H56", 0, "C05432", "cpd03217"),
    ("betacarotene_c", "beta-carotene", "C40H56", 0, "C02094", "cpd01420"),
)
"""The three species this module adds to a host GEM, with what they are.

The last two columns are the load-bearing addition and they are cross-references, not
inventions. Every published metabolite in Yeast9 carries a `kegg.compound` annotation and
`bridge/thermodynamic.py` resolves formation energies through it; these three carried none,
so `ThermodynamicData.load().dgf()` returned None for all three and
`data/pathways/beta_carotene.toml` recorded that as the vendored tables having no energy
for a C40 carotenoid. They have all three. The gap was that the ids in this tuple are
invented HERE, an hour before they were first asked about, and no alias table maps an id
this repository made up. `docs/research/CAROTENOID_ENERGIES.md` is the corrected record.

Provenance, and it is a lookup rather than a judgement: ModelSEED's own
`data/thermo/compounds.tsv` gives the KEGG cross-reference for each compound in its
`aliases` field (`KEGG: C05413` on cpd03205, `C05432` on cpd03217, `C02094` on cpd01420),
each row carrying the InChIKey that identifies the stereoisomer, and the reverse rows
already sit in the vendored `data/thermo/aliases.tsv`. The same three pairs are declared in
`data/thermo/heterologous_metabolites.tsv`, which is what lets the gate resolve them for a
model that does not have this pathway installed;
`tests/test_carotenoid_ec.py::TestTheDeclaredIdentityIsTheSameInBothPlaces` fails if the
two ever disagree. Declaring the identity is not endorsing the number behind it -- see
`data/thermo/refuted_energies.tsv`, which refuses all three energies with a citation.
"""


def _resolve(model: cobra.Model, base_id: str, compartment: str = "c") -> str:
    """Map a bare Yeast9 metabolite id onto whatever convention this model uses.

    Yeast9 SBML uses ``s_0189``; GECKO exports use ``s_0189[c]``. Resolving rather
    than assuming lets one pathway definition serve both, and fails by name when a
    precursor genuinely is not there.
    """
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(
        f"required precursor {base_id!r} not found in {model.id!r} "
        f"(tried {base_id}, {base_id}[{compartment}], {base_id}_{compartment})"
    )


def _add_reaction(model, rid, name, stoichiometry, gene_rule, bounds=(0.0, 1000.0), subsystem=""):
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.subsystem = subsystem
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gene_rule
    return rxn


def add_beta_carotene_pathway(
    model: cobra.Model,
    fpp_id: str = _FPP_C,
    ipp_id: str = _IPP_C,
    ggpp_id: str = _GGPP_C,
    ppi_id: str = _PPI_C,
    oxygen_id: str = _O2_C,
    peroxide_id: str = _H2O2_C,
) -> cobra.Model:
    """Return a copy of ``model`` carrying the crtE/crtYB/crtI pathway.

    The input model is not modified. Works on both plain Yeast9 and GECKO
    enzyme-constrained exports; precursor ids are resolved, not assumed.

    Args:
        model: Host model to extend.
        fpp_id, ipp_id, ggpp_id, ppi_id: Base metabolite ids for the native precursors.
        oxygen_id, peroxide_id: Base metabolite ids for CrtI's electron acceptor and its
            two-electron product. These are the cofactors the pathway actually touches;
            see the note beside the id table for the flavin ids that used to be here.

    Raises:
        ValueError: if the pathway is already installed.
        KeyError: naming any precursor the host model does not carry.
    """
    if model.metabolites.has_id("betacarotene_c"):
        raise ValueError("beta-carotene pathway is already installed on this model")

    out = model.copy()
    fpp = _resolve(out, fpp_id)
    ipp = _resolve(out, ipp_id)
    ggpp = _resolve(out, ggpp_id)
    ppi = _resolve(out, ppi_id)
    oxygen = _resolve(out, oxygen_id)
    peroxide = _resolve(out, peroxide_id)
    out.add_metabolites(
        [
            cobra.Metabolite(mid, name=name, formula=formula, charge=charge,
                             compartment="c")
            for mid, name, formula, charge, _kegg, _seed in _NEW_METABOLITES
        ]
    )
    for mid, _name, _formula, _charge, kegg, seed in _NEW_METABOLITES:
        # An annotation, because that is where a cross-reference belongs and where every
        # other consumer of this model already looks: `bridge/thermodynamic.py` and
        # `bridge/equilibrator.py` both resolve energies through `kegg.compound` and
        # neither knows this module exists. `seed.compound` is not read by anything here
        # and is carried so a reader can check the KEGG id against
        # data/thermo/compounds.tsv without holding the alias table open.
        out.metabolites.get_by_id(mid).annotation.update(
            {"kegg.compound": kegg, "seed.compound": seed})

    # crtE duplicates native BTS1 chemistry; present for construct-level capacity.
    _add_reaction(
        out, "CRTE", "GGPP synthase (crtE)",
        {fpp: -1, ipp: -1, ggpp: 1, ppi: 1},
        "crtE", subsystem="beta-carotene biosynthesis",
    )
    # crtYB, phytoene synthase domain: 2 GGPP -> phytoene + 2 PPi
    _add_reaction(
        out, "CRTYB_PSY", "phytoene synthase (crtYB)",
        {ggpp: -2, "phytoene_c": 1, ppi: 2},
        "crtYB", subsystem="beta-carotene biosynthesis",
    )
    # crtI: four desaturations, electrons to O2, four H2O2 out.
    #
    # **This was written with free FAD as the terminal acceptor until 2026-08-30, and that
    # form is thermodynamically impossible.** Scored by component contribution the reaction
    # as written came to **+166.2 +- 12.8 kJ/mol** -- uphill by thirteen standard errors, on
    # a step a strain measured by Verwaal 2007 (PMID 17496128) runs to completion. The
    # thermodynamic gate reported it `cannot_say`, and the reason was read for a day as
    # "two estimators disagree about carotenoids". They do not. Splitting the reaction shows
    # the disagreement sits on the FLAVIN pair, not the C40 pair.
    #
    # The defect is that FAD in CrtI is a **prosthetic group, not a substrate**. A bacterial
    # -type phytoene desaturase holds one FAD and reoxidises it; writing four free FADH2 as
    # products makes the cell pay for four flavin reductions it never performs, and free
    # FAD/FADH2 is a poor acceptor -- E'0 = -0.219 V -- which is where the +166 comes from.
    #
    # That eQuilibrator is right about flavins was checked rather than assumed, against the
    # one cytosolic FAD reaction yeast-GEM carries. Soluble fumarate reductase `r_0455`:
    # from E'0(fumarate/succinate) = +0.031 V and E'0(FAD/FADH2) = -0.219 V the textbook
    # value is -2 * 96.485 * 0.250 = **-48.2 kJ/mol**; component contribution returns -43.5
    # and the ModelSEED group-contribution table returns -12.6.
    #
    # Written as the flavin oxidase it is, with O2 as terminal acceptor and hydrogen
    # peroxide as the two-electron product, the same chemistry is **-292.7 +- 17.4 kJ/mol**.
    #
    # THE TWO ALTERNATIVES AND WHY NEITHER WAS TAKEN, so the choice is checkable:
    #
    #   4-electron, 2 O2 -> 4 H2O, dGr'0 = -676.5. Rejected on chemistry, not on the number.
    #   Lycopene is C40H56 and phytoene C40H64: no oxygen appears in the product, so O2 is
    #   REDUCED and not INSERTED. That makes this an oxidase and not a monooxygenase, and a
    #   flavoprotein oxidase reduces O2 by two electrons at a single site. H2O2 is what it
    #   makes. Picking the 4-electron form would have been picking the more favourable
    #   number for a mechanism the product formula rules out.
    #
    #   A quinone acceptor, which is the other physiologically argued route for CrtI.
    #   Not reachable here: yeast-GEM v9.0.2 carries ubiquinone (s_1537) in the MITOCHONDRION
    #   only, and CrtI is cytosolic. The three cytosolic species whose names contain "quinone"
    #   are s_3762 benzoquinone (C6H4O2), s_3763 benzosemiquinone and s_4195 hydroquinone --
    #   a xenobiotic detoxification substrate with one reaction between them
    #   (r_4158, NADPH:quinone oxidoreductase), not a respiratory carrier pool. Writing CrtI
    #   onto that couple would be worse biology than writing it onto O2.
    #
    # The choice of acceptor is a modelling decision and is stated as one. What is not a
    # decision is that the shipped form was impossible. Two consequences worth knowing: the
    # pathway now consumes O2, so it cannot run anaerobically, and it produces **four H2O2
    # per beta-carotene** -- the same oxidative stress this repository's Yap1 sensors are
    # built to read. The FBA ceiling barely moves, 0.0873 -> 0.0863 mmol/gDCW/h.
    _add_reaction(
        out, "CRTI", "phytoene desaturase (crtI)",
        {"phytoene_c": -1, oxygen: -4, "lycopene_c": 1, peroxide: 4},
        "crtI", subsystem="beta-carotene biosynthesis",
    )
    # crtYB, lycopene cyclase domain: bicyclisation, no net atom change.
    _add_reaction(
        out, "CRTYB_LCY", "lycopene beta-cyclase (crtYB)",
        {"lycopene_c": -1, "betacarotene_c": 1},
        "crtYB", subsystem="beta-carotene biosynthesis",
    )
    # Accumulates in membranes rather than being exported: demand, not exchange.
    out.add_boundary(
        out.metabolites.get_by_id("betacarotene_c"),
        type="demand", reaction_id=PRODUCT_DEMAND_ID,
    )
    return out
