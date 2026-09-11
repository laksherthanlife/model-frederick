"""Price protein synthesis -- chaperones and a heterologous product -- as a proteome cost.

No chaperone or HSP gene exists in Yeast9, ecYeastGEM or Yeast8, so the Hsf1 -> chaperone ->
burden edge cannot be written as chemistry at all. What a GEM can price is the polymer: the
amino acids, drawn in the model's own ratios, plus the translation ATP. The ODE layer supplies
the flux; this module only says what it costs, which is the direction the GEM works in.

Every stoichiometric coefficient is READ OUT of the host model -- ``r_4047`` for composition,
``r_4041`` for the ATP-hydrolysis convention. The one number typed in is the translation
energy: four ATP equivalents per peptide bond (Stouthamer 1973, PMID 4148026). Two are
already in the model -- every aminoacyl-tRNA synthetase runs ATP -> AMP + PPi and ``r_0568``
hydrolyses the PPi -- so only the two elongation GTPs are added. ``fba/insulin.py`` charges
those two as GTP for one named sequence; they are charged here as ATP, per gram of any
protein, because ``r_0800`` interconverts them one for one and the ec model cannot balance a
GTP hydrolysis.

The price this yields is roughly half the measured one, and that gap is reported rather than
tuned away: see :func:`burden_against_literature`. :func:`scale_protein_pool` prices the same
protein the other way, against the GECKO enzyme pool, and lands near the measured magnitude --
but it does NOT independently predict it. Growth in that model is exactly affine in the pool
cap, so removing a fraction f of the pool costs 1.0051*f of growth at every f: the agreement
with a proportional bound is arithmetic, not a result. What the two prices jointly show is
only that the chemistry alone is ~2.2x too cheap to account for the measured burden.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import cobra
import pandas as pd

__all__ = [
    "ATP_EQUIVALENTS_PER_PEPTIDE_BOND",
    "BIOMASS_REACTION",
    "CHARGING_ATP_EQUIVALENTS",
    "CRT_ENZYME_MASSES_KDA",
    "CRT_ENZYME_OF_REACTION",
    "CRT_KCAT_PROVENANCE",
    "EGUCHI_BURDEN_LIMIT",
    "ELONGATION_ATP_EQUIVALENTS",
    "GAM_REACTION",
    "HAC1_CHAPERONE_ORFS",
    "HSF1_CHAPERONE_ORFS",
    "KAFRI_COST_OF_PROTEIN",
    "METZL_RAZ_RIBOSOME_RESERVE",
    "METZL_RAZ_RIBOSOME_SLOPE",
    "METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION",
    "PROTEIN_PSEUDOREACTION",
    "PROT_POOL_EXCHANGE",
    "PROT_POOL_ID",
    "ProteinComposition",
    "add_chaperone_sink",
    "add_heterologous_protein_sink",
    "add_protein_sink",
    "attach_crt_enzyme_demands",
    "attach_enzyme_demand",
    "burden_against_literature",
    "burden_curve",
    "charge_balance",
    "couple_sink_to_growth",
    "formula_weight",
    "grams_per_gdcw_at_protein_fraction",
    "mass_balance",
    "protein_pool_burden_curve",
    "read_protein_composition",
    "ribosome_allocation_penalty",
    "scale_protein_pool",
    "sink_reaction_id",
]

PROTEIN_PSEUDOREACTION = "r_4047"
GAM_REACTION = "r_4041"
BIOMASS_REACTION = "r_2111"
PROT_POOL_ID = "prot_pool[c]"
PROT_POOL_EXCHANGE = "prot_pool_exchange"

# Stouthamer 1973, PMID 4148026. The split between the two halves is in the module docstring.
ATP_EQUIVALENTS_PER_PEPTIDE_BOND = 4.0
CHARGING_ATP_EQUIVALENTS = 2.0
ELONGATION_ATP_EQUIVALENTS = ATP_EQUIVALENTS_PER_PEPTIDE_BOND - CHARGING_ATP_EQUIVALENTS

# Hsf1 regulon chaperones (SSA1, SSA2, HSP104, HSP82). None is in any of the three models.
HSF1_CHAPERONE_ORFS = ("YAL005C", "YLL024C", "YLL026W", "YPL240C")
# Hac1/UPR targets: KAR2 and PDI1. Also absent from all three.
HAC1_CHAPERONE_ORFS = ("YJL034W", "YCL043C")

#: Eguchi et al. 2018, eLife 7:e34595: the overexpression limits of *some* glycolytic
#: proteins "were up to 15% of the total cellular protein". A threshold, not a slope.
EGUCHI_BURDEN_LIMIT = 0.15

#: Metzl-Raz et al. 2017, eLife 6:e28034 (PMID 28857745), Fig. 2A: phi_R = 0.35*X + 0.08 with X
#: in GENERATIONS/h, not mu. The paper converts it itself -- mu = X*ln(2), "Delta_r/Delta_mu =
#: 21/ln(2) [min]" -- so the slope against mu in 1/h is 0.35/ln(2), not 0.35.
METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION = 0.35
METZL_RAZ_RIBOSOME_SLOPE = METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION / math.log(2.0)
#: The excess ribosome reserve, "a constant ~8% of the proteome", carried at every growth rate.
METZL_RAZ_RIBOSOME_RESERVE = 0.08

KAFRI_COST_OF_PROTEIN = (
    "Kafri, Metzl-Raz, Jona & Barkai 2016, Cell Reports 14:22-31 (PMID 26725116): the cost "
    "of unneeded protein is transcription-limited in low phosphate, translation-limited in "
    "low nitrogen and both in rich media, so no single ribosome-allocation coefficient "
    "applies across media and none is used here."
)

#: No kcat is published for any of the four heterologous carotenoid steps, so none is
#: shipped. The repo's fitted cyclase vmax is per gDCW, not per enzyme, so it is not one.
CRT_KCAT_PROVENANCE = {
    "CRTE": None,
    "CRTYB_PSY": None,
    "CRTI": None,
    "CRTYB_LCY": None,
}

#: UniProt sequence masses of the *X. dendrorhous* genes Verwaal 2007 expressed -- crtE
#: Q1L6K3, crtYB Q7Z859, crtI O13506. A lookup, unlike the kcats, which are refused above.
CRT_ENZYME_MASSES_KDA = {"crtE": 42.153, "crtYB": 74.736, "crtI": 65.091}

#: crtYB is bifunctional, so both domains draw on one enzyme.
CRT_ENZYME_OF_REACTION = {
    "CRTE": "crtE",
    "CRTYB_PSY": "crtYB",
    "CRTI": "crtI",
    "CRTYB_LCY": "crtYB",
}

_ATP_FORMULA_ELEMENTS = {"C": 10, "H": 12, "N": 5, "O": 13, "P": 3}
_ADP_FORMULA_ELEMENTS = {"C": 10, "H": 12, "N": 5, "O": 10, "P": 2}
_WATER_ELEMENTS = {"H": 2, "O": 1}
_PHOSPHATE_ELEMENTS = {"H": 1, "O": 4, "P": 1}
_PROTON_ELEMENTS = {"H": 1}


def _elements(metabolite: cobra.Metabolite) -> dict[str, float]:
    """cobra warns on the fractional counts a lumped polymer needs; that warning is expected."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return dict(metabolite.elements or {})


def mass_balance(reaction: cobra.Reaction) -> dict[str, float]:
    """``reaction.check_mass_balance()`` with the fractional-formula warning silenced.

    Formula-less metabolites -- GECKO's ``prot_*`` species -- are invisible to this, so on an
    ec model it balances the metabolic half only.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return {k: v for k, v in reaction.check_mass_balance().items() if k != "charge"}


def formula_weight(metabolite: cobra.Metabolite) -> float:
    """``metabolite.formula_weight`` with the fractional-formula warning silenced."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return float(metabolite.formula_weight)


def charge_balance(reaction: cobra.Reaction) -> float:
    """Net charge of a reaction, products minus reactants."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return float(reaction.check_mass_balance().get("charge", 0.0))


def _formula(elements: dict[str, float]) -> str:
    return "".join(f"{el}{elements[el]}" for el in sorted(elements) if elements[el] != 0)


@dataclass(frozen=True)
class ProteinComposition:
    """Per-gram protein composition read out of a model's own protein pseudoreaction."""

    source_reaction: str
    grams_per_unit: float
    residues_per_unit: float
    aa_coefficients: dict[str, float] = field(default_factory=dict)
    trna_coefficients: dict[str, float] = field(default_factory=dict)
    polymer_elements: dict[str, float] = field(default_factory=dict)
    polymer_charge: float = 0.0

    @property
    def residues_per_gram(self) -> float:
        """mmol of peptide-bonded residue per gram of protein."""
        return self.residues_per_unit / self.grams_per_unit

    @property
    def mean_residue_mass(self) -> float:
        """g/mol of the average peptide-bonded residue. ~110 for any real proteome."""
        return 1000.0 / self.residues_per_gram

    @property
    def polymer_elements_per_gram(self) -> dict[str, float]:
        """Elemental composition of one gram of the polymer, in mmol."""
        return {el: n / self.grams_per_unit for el, n in self.polymer_elements.items()}

    @property
    def polymer_charge_per_gram(self) -> float:
        return self.polymer_charge / self.grams_per_unit

    def summary(self) -> str:
        return (
            f"{self.source_reaction}: {self.grams_per_unit:.6f} g protein per unit, "
            f"{self.residues_per_unit:.6f} mmol residues per unit "
            f"({self.residues_per_gram:.5f} mmol/g, mean residue {self.mean_residue_mass:.3f} "
            f"g/mol), polymer {_formula(self.polymer_elements)} charge {self.polymer_charge:+.4f}"
        )


def _atomic_mass_table() -> dict[str, float]:
    from cobra.core.formula import elements_and_molecular_weights

    return elements_and_molecular_weights


def read_protein_composition(
    model: cobra.Model, reaction_id: str = PROTEIN_PSEUDOREACTION
) -> ProteinComposition:
    """Read the amino-acid composition out of the model's protein pseudoreaction.

    Every reactant must be a charged tRNA (formula carries the ``R`` placeholder) and exactly
    one product must be the formula-less protein pseudo-metabolite.

    Args:
        model: Yeast9 or an ecYeastGEM export.
        reaction_id: Protein pseudoreaction to read.

    Raises:
        ValueError: if the reaction is not shaped like a protein pseudoreaction.
    """
    rxn = model.reactions.get_by_id(reaction_id)
    reactants = {m: -c for m, c in rxn.metabolites.items() if c < 0}
    products = {m: c for m, c in rxn.metabolites.items() if c > 0}
    bad = [m.id for m in reactants if "R" not in _elements(m)]
    if bad:
        raise ValueError(f"{reaction_id} reactants are not all charged tRNAs: {bad}")
    protein = [m for m in products if not m.formula]
    if len(protein) != 1:
        raise ValueError(
            f"{reaction_id} must have exactly one formula-less protein product, got "
            f"{[m.id for m in protein]}"
        )
    protein_met = protein[0]
    per_unit = products[protein_met]
    residues = sum(reactants.values()) / per_unit

    # Everything not the protein must be its elements; the R placeholders cancel out.
    raw: dict[str, float] = {}
    raw_charge = 0.0
    for met, coef in rxn.metabolites.items():
        if met is protein_met:
            continue
        for el, n in _elements(met).items():
            raw[el] = raw.get(el, 0.0) - coef * n / per_unit
        raw_charge -= coef * (met.charge or 0.0) / per_unit
    raw.pop("R", None)

    # A charged tRNA still carries its alpha-amino proton; the peptide bond releases it.
    proton = _proton_metabolite(model)
    polymer = dict(raw)
    for el, n in _PROTON_ELEMENTS.items():
        polymer[el] = polymer.get(el, 0.0) - residues * n
    polymer = {k: v for k, v in polymer.items() if abs(v) > 0}
    polymer_charge = raw_charge - residues * (proton.charge or 0.0)

    masses = _atomic_mass_table()
    grams = sum(polymer[el] * masses[el] for el in polymer) / 1000.0
    return ProteinComposition(
        source_reaction=reaction_id,
        grams_per_unit=grams,
        residues_per_unit=residues,
        aa_coefficients={m.id: c / per_unit for m, c in reactants.items()},
        trna_coefficients={m.id: c / per_unit for m, c in products.items() if m is not protein_met},
        polymer_elements=polymer,
        polymer_charge=polymer_charge,
    )


def _match(metabolite: cobra.Metabolite, elements: dict[str, float]) -> bool:
    return _elements(metabolite) == elements


def _atp_hydrolysis_terms(
    model: cobra.Model, gam_reaction: str = GAM_REACTION
) -> dict[cobra.Metabolite, float]:
    """The model's own ATP -> ADP + Pi + H+ convention, normalised to one ATP.

    Taken from the growth-associated maintenance term rather than written out, so the ec
    model's uncharged formula convention and Yeast9's charged one both come out balanced.
    """
    rxn = model.reactions.get_by_id(gam_reaction)
    wanted = {
        "atp": _ATP_FORMULA_ELEMENTS,
        "h2o": _WATER_ELEMENTS,
        "adp": _ADP_FORMULA_ELEMENTS,
        "pi": _PHOSPHATE_ELEMENTS,
        "h": _PROTON_ELEMENTS,
    }
    found: dict[str, tuple[cobra.Metabolite, float]] = {}
    for met, coef in rxn.metabolites.items():
        for key, elements in wanted.items():
            if _match(met, elements):
                found[key] = (met, coef)
    missing = sorted(set(wanted) - set(found))
    if missing:
        raise ValueError(f"{gam_reaction} carries no {missing} term; cannot read the ATP cost")
    gam = abs(found["atp"][1])
    if gam <= 0:
        raise ValueError(f"{gam_reaction} has a zero ATP coefficient")
    terms = {met: coef / gam for met, coef in found.values()}
    residual: dict[str, float] = {"charge": 0.0}
    for met, coef in terms.items():
        residual["charge"] += coef * (met.charge or 0.0)
        for el, n in _elements(met).items():
            residual[el] = residual.get(el, 0.0) + coef * n
    unbalanced = {k: v for k, v in residual.items() if abs(v) > 1e-9}
    if unbalanced:
        raise ValueError(f"{gam_reaction}'s ATP terms are not balanced: {unbalanced}")
    return terms


def _proton_metabolite(model: cobra.Model, gam_reaction: str = GAM_REACTION) -> cobra.Metabolite:
    for met in model.reactions.get_by_id(gam_reaction).metabolites:
        if _match(met, _PROTON_ELEMENTS):
            return met
    raise ValueError(f"{gam_reaction} carries no proton term")


def sink_reaction_id(protein_id: str) -> str:
    return f"PROTSYN_{protein_id}"


def add_protein_sink(
    model: cobra.Model,
    protein_id: str,
    name: str = "",
    gene_rule: str = "",
    atp_per_residue: float = ELONGATION_ATP_EQUIVALENTS,
    composition_reaction: str = PROTEIN_PSEUDOREACTION,
    gam_reaction: str = GAM_REACTION,
    add_demand: bool = True,
) -> cobra.Model:
    """Return a copy of ``model`` carrying a protein-synthesis sink priced per gram.

    One unit of flux through ``PROTSYN_<protein_id>`` is one gram of that protein per gDCW per
    hour: it draws the same charged tRNAs the biomass protein draws, in the model's own ratios,
    releases the free tRNAs and one proton per peptide bond, and hydrolyses
    ``atp_per_residue`` ATP equivalents per residue.

    Args:
        model: Host model. Not modified.
        protein_id: Short id; becomes ``<id>_c``, ``PROTSYN_<id>``, ``DM_<id>_c``.
        name: Human-readable metabolite name.
        gene_rule: GPR for the sink, e.g. the Hsf1 chaperone regulon.
        atp_per_residue: ATP equivalents per peptide bond charged here. Default is the two
            elongation GTPs; the two charging equivalents are already in the model.
        composition_reaction: Protein pseudoreaction supplying the composition.
        gam_reaction: Reaction supplying the model's ATP-hydrolysis convention.
        add_demand: Add ``DM_<id>_c`` so the protein can accumulate.

    Raises:
        ValueError: if the sink is already installed or ``atp_per_residue`` is negative.
    """
    if atp_per_residue < 0:
        raise ValueError(f"atp_per_residue must be >= 0, got {atp_per_residue}")
    met_id = f"{protein_id}_c"
    if model.metabolites.has_id(met_id):
        raise ValueError(f"protein sink {protein_id!r} is already installed on this model")

    out = model.copy()
    comp = read_protein_composition(out, composition_reaction)
    atp_terms = _atp_hydrolysis_terms(out, gam_reaction)
    proton = _proton_metabolite(out, gam_reaction)
    per_gram = 1.0 / comp.grams_per_unit
    residues_per_gram = comp.residues_per_gram

    protein = cobra.Metabolite(
        met_id,
        name=name or f"{protein_id} (heterologous protein)",
        formula=_formula(comp.polymer_elements_per_gram),
        charge=comp.polymer_charge_per_gram,
        compartment=proton.compartment,
    )
    out.add_metabolites([protein])

    stoich: dict[cobra.Metabolite, float] = {}
    for mid, coef in comp.aa_coefficients.items():
        stoich[out.metabolites.get_by_id(mid)] = -coef * per_gram
    for mid, coef in comp.trna_coefficients.items():
        stoich[out.metabolites.get_by_id(mid)] = coef * per_gram
    for met, coef in atp_terms.items():
        stoich[met] = stoich.get(met, 0.0) + coef * atp_per_residue * residues_per_gram
    stoich[proton] = stoich.get(proton, 0.0) + residues_per_gram
    stoich[protein] = 1.0

    rxn = cobra.Reaction(
        sink_reaction_id(protein_id),
        name=f"{protein_id} synthesis (g/gDCW/h)",
        lower_bound=0.0,
        upper_bound=1000.0,
    )
    rxn.subsystem = "heterologous protein synthesis"
    out.add_reactions([rxn])
    rxn.add_metabolites(stoich)
    if gene_rule:
        rxn.gene_reaction_rule = gene_rule
    if add_demand:
        out.add_boundary(protein, type="demand", reaction_id=f"DM_{met_id}")
    return out


def add_chaperone_sink(
    model: cobra.Model,
    protein_id: str = "chaperone",
    orfs: tuple[str, ...] = HSF1_CHAPERONE_ORFS + HAC1_CHAPERONE_ORFS,
    **kwargs,
) -> cobra.Model:
    """Install a chaperone-synthesis sink carrying the Hsf1/Hac1 regulon as its GPR.

    The ORFs are absent from every yeast model, so the GPR exists only so the ODE layer can
    address the reaction by the transcription factor driving it.
    """
    return add_protein_sink(
        model,
        protein_id,
        name=kwargs.pop("name", "chaperone protein (Hsf1/Hac1 regulon)"),
        gene_rule=kwargs.pop("gene_rule", " or ".join(orfs)),
        **kwargs,
    )


def add_heterologous_protein_sink(
    model: cobra.Model, protein_id: str = "hetprot", **kwargs
) -> cobra.Model:
    """Install a sink for a recombinant product protein of average yeast composition."""
    return add_protein_sink(
        model,
        protein_id,
        name=kwargs.pop("name", "heterologous product protein"),
        **kwargs,
    )


def grams_per_gdcw_at_protein_fraction(
    composition: ProteinComposition, fraction_of_total_protein: float
) -> float:
    """Convert 'x% of total cellular protein' into grams per gDCW.

    Biomass composition is pinned, so a heterologous protein is additive rather than
    displacing: ``het / (het + native) = f``.
    """
    if not 0.0 <= fraction_of_total_protein < 1.0:
        raise ValueError(f"fraction must be in [0, 1), got {fraction_of_total_protein}")
    f = fraction_of_total_protein
    return composition.grams_per_unit * f / (1.0 - f)


def couple_sink_to_growth(
    model: cobra.Model,
    sink_reaction: str,
    grams_per_gdcw: float,
    biomass_reaction: str = BIOMASS_REACTION,
):
    """Pin sink flux to ``grams_per_gdcw * mu``, in place. Scope it with ``with model:``.

    In balanced exponential growth a cell holding g grams of the protein per gDCW must
    synthesise it at g*mu to keep up with dilution; no chemostat is implied, only that mu is
    the same mu the biomass reaction carries.
    """
    sink = model.reactions.get_by_id(sink_reaction)
    biomass = model.reactions.get_by_id(biomass_reaction)
    constraint = model.problem.Constraint(
        sink.flux_expression - grams_per_gdcw * biomass.flux_expression,
        lb=0.0,
        ub=0.0,
        name=f"burden_{sink_reaction}",
    )
    model.add_cons_vars([constraint])
    return constraint


def burden_curve(
    model: cobra.Model,
    grams_per_gdcw: Sequence[float],
    sink_reaction: str = "PROTSYN_hetprot",
    biomass_reaction: str = BIOMASS_REACTION,
    constrain: Callable[[cobra.Model], object] | None = None,
) -> pd.DataFrame:
    """Growth as a function of grams of protein carried per gDCW. Does not modify ``model``.

    Args:
        model: Model already carrying the sink.
        grams_per_gdcw: Protein loads to score.
        sink_reaction: Sink to load.
        biomass_reaction: Objective.
        constrain: Optional callable applied to the scoped model to set the uptake regime.

    Raises:
        ValueError: if no zero load is included; ``relative_growth`` needs that baseline.
    """
    if not any(load == 0.0 for load in grams_per_gdcw):
        raise ValueError(
            "grams_per_gdcw must include 0.0: relative_growth is measured against the "
            f"unloaded model, and {list(grams_per_gdcw)} carries no such reference"
        )
    rows = []
    for load in grams_per_gdcw:
        with model as m:
            if constrain is not None:
                constrain(m)
            m.objective = biomass_reaction
            couple_sink_to_growth(m, sink_reaction, load, biomass_reaction)
            solution = m.optimize()
            growth = (
                float(solution.objective_value) if solution.status == "optimal" else float("nan")
            )
            rows.append(
                {
                    "grams_per_gdcw": load,
                    "growth": growth,
                    "sink_flux": growth * load,
                    "status": solution.status,
                }
            )
    frame = pd.DataFrame(rows)
    base = float(frame.loc[frame["grams_per_gdcw"] == 0.0, "growth"].iloc[0])
    frame["relative_growth"] = frame["growth"] / base
    return frame


def ribosome_allocation_penalty(
    fraction_of_total_protein: float,
    growth_rate: float,
    slope: float = METZL_RAZ_RIBOSOME_SLOPE,
    reserve: float = METZL_RAZ_RIBOSOME_RESERVE,
) -> float:
    """Relative growth loss the measured growth law predicts for a diverted proteome fraction.

    phi_R = reserve + slope*mu is measured (Metzl-Raz 2017), with ``slope`` already converted
    out of the paper's generations/h axis into mu in 1/h. What is ASSUMED here, and is the
    only assumption, is that a heterologous fraction f displaces every native sector including
    the ribosomes, so phi_R scales by (1-f) and the loss is f*(reserve + slope*mu)/(slope*mu).
    With reserve = 0 it collapses to the naive bound, loss = f.

    Raises:
        ValueError: on a fraction outside [0, 1) or a non-positive growth rate.
    """
    if not 0.0 <= fraction_of_total_protein < 1.0:
        raise ValueError(f"fraction must be in [0, 1), got {fraction_of_total_protein}")
    if growth_rate <= 0.0:
        raise ValueError(f"growth rate must be > 0, got {growth_rate}")
    return fraction_of_total_protein * (reserve + slope * growth_rate) / (slope * growth_rate)


def burden_against_literature(
    model: cobra.Model,
    fraction_of_total_protein: float = EGUCHI_BURDEN_LIMIT,
    sink_reaction: str = "PROTSYN_hetprot",
    composition_reaction: str = PROTEIN_PSEUDOREACTION,
    biomass_reaction: str = BIOMASS_REACTION,
    constrain: Callable[[cobra.Model], object] | None = None,
) -> pd.DataFrame:
    """What this GEM charges for a protein load, beside what the burden literature measures.

    Four rows, and the point is that they disagree: the GEM prices amino acids and translation
    ATP only, which is roughly half of the measured cost. ``gem_over_reference`` below 1 means
    the GEM is too cheap. Kafri's row carries no number on purpose -- see
    :data:`KAFRI_COST_OF_PROTEIN`.

    Args:
        model: Model already carrying the sink.
        fraction_of_total_protein: Load, as a fraction of total cellular protein.
        sink_reaction: Sink to load.
        composition_reaction: Protein pseudoreaction supplying grams per gDCW.
        biomass_reaction: Objective.
        constrain: Optional callable applied to the scoped model to set the uptake regime.
    """
    composition = read_protein_composition(model, composition_reaction)
    load = grams_per_gdcw_at_protein_fraction(composition, fraction_of_total_protein)
    curve = burden_curve(
        model, [0.0, load], sink_reaction=sink_reaction,
        biomass_reaction=biomass_reaction, constrain=constrain,
    )
    base_growth = float(curve.iloc[0]["growth"])
    gem_loss = 1.0 - float(curve.iloc[1]["relative_growth"])
    proportional = fraction_of_total_protein
    growth_law = ribosome_allocation_penalty(fraction_of_total_protein, base_growth)
    rows = [
        {
            "source": "this GEM: amino acids + translation ATP",
            "basis": "measured here",
            "relative_growth_loss": gem_loss,
            "note": f"{load:.6f} g/gDCW at mu = {base_growth:.6f} /h",
        },
        {
            "source": "proportional allocation, Eguchi 2018 eLife 7:e34595 limit",
            "basis": "measured limit + proportional-displacement assumption",
            "relative_growth_loss": proportional,
            "note": "expression limit of some glycolytic proteins, up to 15% of total protein",
        },
        {
            "source": "growth law, Metzl-Raz 2017 (PMID 28857745)",
            "basis": "measured phi_R = 0.08 + (0.35/ln2) mu + displacement assumption",
            "relative_growth_loss": growth_law,
            "note": "the reserve term makes it steeper than proportional at every finite mu",
        },
        {
            "source": "Kafri 2016 (PMID 26725116)",
            "basis": "refused: no single coefficient",
            "relative_growth_loss": float("nan"),
            "note": KAFRI_COST_OF_PROTEIN,
        },
    ]
    frame = pd.DataFrame(rows)
    frame.insert(0, "fraction_of_total_protein", fraction_of_total_protein)
    frame["gem_over_reference"] = gem_loss / frame["relative_growth_loss"]
    return frame


def scale_protein_pool(
    model: cobra.Model,
    fraction_of_total_protein: float,
    pool_exchange: str = PROT_POOL_EXCHANGE,
) -> cobra.Model:
    """Take a heterologous protein fraction out of a GECKO enzyme pool. Returns a copy.

    The other half of the burden, and the half a plain GEM cannot express: a protein that is
    not an enzyme still occupies proteome, so the enzyme budget available to everything else
    scales by (1-f). Same displacement assumption as :func:`ribosome_allocation_penalty`.

    Raises:
        ValueError: on a fraction outside [0, 1), or a model with no protein pool.
    """
    if not 0.0 <= fraction_of_total_protein < 1.0:
        raise ValueError(f"fraction must be in [0, 1), got {fraction_of_total_protein}")
    if not model.reactions.has_id(pool_exchange):
        raise ValueError(
            f"{model.id!r} has no {pool_exchange}: this is not an enzyme-constrained model, so "
            "there is no proteome budget to take a share of"
        )
    out = model.copy()
    pool = out.reactions.get_by_id(pool_exchange)
    pool.upper_bound = pool.upper_bound * (1.0 - fraction_of_total_protein)
    return out


def protein_pool_burden_curve(
    model: cobra.Model,
    fractions: Sequence[float],
    pool_exchange: str = PROT_POOL_EXCHANGE,
    biomass_reaction: str = BIOMASS_REACTION,
) -> pd.DataFrame:
    """Growth against heterologous proteome fraction, priced by the GECKO pool alone.

    Reported beside :func:`burden_against_literature` rather than added to it: the two price
    different things, and adding them would double-charge the same protein. Read the result as
    proportional displacement restated, not as a prediction -- growth is affine in the pool cap,
    so ``relative_growth_loss`` is 1.0051*fraction identically, whatever the fraction.
    """
    if not model.reactions.has_id(pool_exchange):
        raise ValueError(f"{model.id!r} has no {pool_exchange}")
    if not any(fraction == 0.0 for fraction in fractions):
        raise ValueError(
            "fractions must include 0.0: relative_growth is measured against the unloaded "
            f"pool, and {list(fractions)} carries no such reference"
        )
    cap = model.reactions.get_by_id(pool_exchange).upper_bound
    rows = []
    for fraction in fractions:
        if not 0.0 <= fraction < 1.0:
            raise ValueError(f"fraction must be in [0, 1), got {fraction}")
        with model as m:
            m.objective = biomass_reaction
            m.reactions.get_by_id(pool_exchange).upper_bound = cap * (1.0 - fraction)
            solution = m.optimize()
            growth = (
                float(solution.objective_value) if solution.status == "optimal" else float("nan")
            )
        rows.append(
            {
                "fraction_of_total_protein": fraction,
                "pool_cap": cap * (1.0 - fraction),
                "growth": growth,
                "status": solution.status,
            }
        )
    frame = pd.DataFrame(rows)
    base = float(frame.loc[frame["fraction_of_total_protein"] == 0.0, "growth"].iloc[0])
    frame["relative_growth"] = frame["growth"] / base
    frame["relative_growth_loss"] = 1.0 - frame["relative_growth"]
    return frame


def _declared_without_citation(provenance: str) -> bool:
    return any(word in provenance.upper() for word in ("PLACEHOLDER", "ASSERTED", "REFUSED"))


def attach_enzyme_demand(
    model: cobra.Model,
    reaction_id: str,
    enzyme_id: str,
    kcat_per_s: float,
    molecular_weight_kda: float,
    provenance: str,
    pool_id: str = PROT_POOL_ID,
) -> cobra.Model:
    """Charge a heterologous reaction against a GECKO protein pool. Returns a copy.

    The installed carotenoid reactions carry no protein cost at all, so in the ec model they
    are free. This is the missing edge: a draw reaction ``MW * prot_pool -> prot_<enzyme>``
    plus a ``1/kcat`` enzyme coefficient on the reaction, which is GECKO's own convention. An
    enzyme id already present is REUSED, so a bifunctional enzyme carries one pool draw.

    Args:
        model: A GECKO export. Not modified.
        reaction_id: Reaction to charge.
        enzyme_id: Enzyme tag; becomes ``prot_<id>[c]`` and ``draw_prot_<id>``.
        kcat_per_s: Turnover number, 1/s.
        molecular_weight_kda: Enzyme mass, kDa (= g/mmol, GECKO's draw unit).
        provenance: Citation, or a declaration containing PLACEHOLDER/ASSERTED. Required, and
            for the four crt reactions a citation is REFUSED because no kcat is published.

    Raises:
        ValueError: if the model has no protein pool, the reaction is already charged, the
            reused enzyme's mass disagrees, or kcat/MW/provenance are not usable.
    """
    if not provenance:
        raise ValueError("provenance is required: cite the kcat or declare it PLACEHOLDER")
    if reaction_id in CRT_KCAT_PROVENANCE and CRT_KCAT_PROVENANCE[reaction_id] is None:
        if not _declared_without_citation(provenance):
            raise ValueError(
                f"no kcat is published for {reaction_id}; provenance {provenance!r} reads as a "
                "citation. Declare it PLACEHOLDER or ASSERTED instead"
            )
    if kcat_per_s <= 0 or molecular_weight_kda <= 0:
        raise ValueError(
            f"kcat ({kcat_per_s} 1/s) and MW ({molecular_weight_kda} kDa) must both be positive"
        )
    if not model.metabolites.has_id(pool_id):
        raise ValueError(
            f"{model.id!r} has no {pool_id}: this is not an enzyme-constrained model, so there "
            "is no proteome to charge"
        )
    target = model.reactions.get_by_id(reaction_id)
    if any(m.id.startswith("prot_") for m in target.metabolites):
        raise ValueError(f"{reaction_id} already carries an enzyme cost")

    out = model.copy()
    pool = out.metabolites.get_by_id(pool_id)
    enzyme_met_id, draw_id = f"prot_{enzyme_id}[c]", f"draw_prot_{enzyme_id}"
    if out.metabolites.has_id(enzyme_met_id):
        enzyme = out.metabolites.get_by_id(enzyme_met_id)
        draw = out.reactions.get_by_id(draw_id)
        existing = -draw.metabolites[pool]
        if not math.isclose(existing, molecular_weight_kda, rel_tol=1e-9):
            raise ValueError(
                f"{draw_id} already draws {existing} kDa from the pool, not {molecular_weight_kda}"
            )
    else:
        enzyme = cobra.Metabolite(
            enzyme_met_id, name=f"prot_{enzyme_id} [cytoplasm]", compartment=pool.compartment
        )
        out.add_metabolites([enzyme])
        draw = cobra.Reaction(draw_id, name=draw_id, lower_bound=0.0, upper_bound=1000.0)
        out.add_reactions([draw])
        # GECKO's own convention: the draw coefficient is the mass in kDa, i.e. g per mmol.
        draw.add_metabolites({pool: -molecular_weight_kda, enzyme: 1.0})
    draw.notes["kcat_provenance"] = provenance
    charged = out.reactions.get_by_id(reaction_id)
    charged.add_metabolites({enzyme: -1.0 / (kcat_per_s * 3600.0)})
    charged.notes["kcat_provenance"] = provenance
    return out


def attach_crt_enzyme_demands(
    model: cobra.Model,
    kcats_per_s: dict[str, float],
    provenance: str,
    masses_kda: dict[str, float] = CRT_ENZYME_MASSES_KDA,
) -> cobra.Model:
    """Charge the installed crt pathway against the GECKO pool. Returns a copy.

    The masses are UniProt sequence masses; the kcats are not published and must be supplied
    by the caller with a provenance declaring them PLACEHOLDER or ASSERTED. crtYB's two
    domains share one enzyme, so the pathway costs three proteins and not four.

    Args:
        model: A GECKO export carrying the pathway. Not modified.
        kcats_per_s: Reaction id -> kcat, one entry per reaction to charge.
        provenance: Declaration for all of them; a citation is refused.
        masses_kda: Enzyme tag -> mass in kDa.

    Raises:
        ValueError: naming any reaction that is not one of the four crt reactions.
    """
    unknown = sorted(set(kcats_per_s) - set(CRT_ENZYME_OF_REACTION))
    if unknown:
        raise ValueError(f"not carotenoid reactions: {unknown}")
    # An empty request still owes the caller a copy; attach_enzyme_demand makes one per entry.
    out = model.copy() if not kcats_per_s else model
    for reaction_id, kcat in kcats_per_s.items():
        enzyme = CRT_ENZYME_OF_REACTION[reaction_id]
        out = attach_enzyme_demand(
            out, reaction_id, enzyme, kcat, masses_kda[enzyme], provenance
        )
    return out
