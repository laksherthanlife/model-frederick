"""The pathway-level audit: the owner's stress diagram, item by item, against what exists.

The reference table is ``data/gem/stress_pathway_coverage.tsv``: 97 stable item keys across
SENSOR, CASCADE, TF and ADAPTATION, preserving the ORIGINAL DIAGRAM's placement rather than
asserting a corrected molecular ontology. Its verdicts, reaction ids and ablation results
refer to yeast-GEM v9.0.2 and the older phenomenological panel. They must not be transferred
to the EC model without a fresh query. ``data/gem/stress_response_map.json`` is a separate
conceptual atlas: Hog1 is a MAPK, Pma1/V-ATPase are effectors, Rim101 primarily responds to
alkaline pH, and SPS senses extracellular amino acids rather than quorum. Neither catalogue
membership nor that expanded map supplies quantitative signaling kinetics.

THE VERDICT THAT CARRIES THE HISTORICAL AUDIT IS ``NOT_METABOLIC``. It identifies processes
outside this metabolic reconstruction's represented chemistry, not proteins that cannot
catalyse reactions: kinases phosphorylate proteins, Ire1 processes RNA, and recombinases
act on DNA. A signaling process can be legitimately absent from the GEM while still being
an essential missing mechanism in a whole-cell stress model. Catalogue membership and a
GPR association do not establish dynamic signaling coverage. ``model_gene_coverage`` queries
the supplied model directly and reports only the associations it actually observes.
Three catalogued sensor/cascade rows are absent from the model AND unrepresented at the panel layer,
which is a hole of a different kind: :func:`unrepresented_pathways` is what returns them, and
it is disjoint from :func:`gaps` by construction.

WHAT THE AUDIT FOUND, in descending order of how much it should change someone's plans.

1. FOUR PANEL MODULES REACH ZERO REACTIONS. Measured by running
   ``bridge.regulation.load_regulons`` against the shipped GEM: UPR, calcium, zinc and
   xenobiotic reach no reaction at all, and cell_wall and alkaline_ph reach exactly one each.
   calcium is the sharpest -- CRZ1 has no target whatsoever at the MacIsaac evidence level, so
   its regulon is empty and ``gene_scales`` can scale nothing. In every one of those four cases
   the effector reactions ARE in the model, and :func:`panel_module_reactions` reaches 4, 3, 6
   and 3 of them: Ena1/Pmc1/Pmr1 with Mid1-Cch1 after the restore for calcium, Zrt1/Zrt2 and
   Zrc1/Cot1/Yke4 for zinc, Tpo1's six polyamine transports for xenobiotic, and the three lipid
   pseudoreactions that price ER expansion for UPR. The regulon route cannot see any of them.
   :data:`UNREACHED_MODULES` names them and :func:`panel_module_reactions` is the declared map
   that connects them, since fixing the evidence level is a change to ``bridge/regulation.py``.

2. TWO ITEMS A SYMBOL SWEEP CALLS ABSENT ARE PRESENT. GPX3 is in yeast-GEM under its systematic
   name HYR1 (YIR037W) on r_0483, and RHO1 -- the CWI cascade's terminal GTPase -- is in it as a
   required conjunct of every clause of the 1,3-beta-glucan synthase rule r_0005. Both are found
   by resolving ORFs rather than symbols. :func:`verify_coverage` checks the historical
   reaction declarations when explicitly called; :func:`model_gene_coverage` queries the
   associations in the supplied model independently of those declarations.

3. THE cAMP SECOND MESSENGER IS FULLY IN THE GEM AND NO MODULE READS IT. Adenylate cyclase
   (r_0147, CYR1) and five phosphodiesterase reactions (r_0033-r_0037, PDE1/PDE2) sit in Purine
   metabolism. No ``MODULES`` entry carries cAMP-PKA, so Gpr1's whole pathway is unrepresented
   at the panel layer while its chemistry is represented at the model layer. Blocking all six
   moves growth by 1.7e-8 /h shipped and 0.000000 /h at reference batch, so this is reported
   as coverage and not proposed as a block; :func:`unrepresented_pathways` returns the two
   regulator rows, and deliberately not this one.

4. THE V-ATPase IS A SECOND CYTOSOLIC PROTON SINK, AND ``fba/stress_ph.py`` IS RIGHT TO IGNORE
   IT. r_1085 and r_1086 each consume one cytosolic H+ (s_0794) per ATP; r_0227 consumes none,
   exporting the hydrolysis proton to s_0796 instead, so the two are not the same reaction
   written twice. The suspicion this raises -- a second escape that undercuts the proton budget
   -- is REFUTED by measurement: with ``stress_ph.constrain_proton_budget`` applied under the
   reference batch regime, ``proton_price`` is 1.0000 ATP with the V-ATPases open and 1.0000
   with both blocked, both carry zero flux at the optimum, and blocking them moves growth by
   0.000000 /h. The omission is correct, not a defect, and this row is the audit's own
   false-positive caught before it became a bug report.

5. THE BASELINE PHENOMENOLOGICAL PANEL HAS NO NEGATIVE TRANSCRIPTIONAL TARGET WEIGHTS.
   Across its 25 stressors there are exactly five negative target weights and all five are
   on metabolite pools, which have no transcription factor or regulon. This describes
   ``generator/stress_panel.py``, not every model in the repository: the family registry
   has a negative feedback edge, and the synthetic teacher has inhibitory global controls.
   Neither supplies measured gene-specific repression. Dot6 and Tod6 are absent from
   ``MODULES`` and the capacity coupling is REFUSED here -- see
   :data:`REFUSED_REPRESSION_COUPLING`, which raises on ``float()``.

6. THREE ADAPTATIONS THAT THE GENE LAYER RESTORES HAD NO PATHWAY ROW, AND NOW DO. The seam
   between this table and ``fba/stress_genes.py`` is checkable: every gene that module actually
   puts into the model should belong to some item here, or a gene was restored for a function
   nothing represents. It did not hold. ``adaptation.manganese_homeostasis`` (r_4590, SMF1 shipped
   plus SMF2 restored) closes SMF2; ``adaptation.n_glycosylation`` (ALG3_ER) closes ALG3, the step
   DTT and tunicamycin actually hit; and ``adaptation.cell_wall_mannan`` named six of MANPOL_G's
   ten transferases and now names all ten. Manganese is the sharper of the two new rows: blocking
   r_4590 is LETHAL in both regimes, the same biomass-requirement signature as Nha1 and Trk1, and
   SOD2 on r_4190 is the MANGANESE superoxide dismutase whose cofactor the GEM's stoichiometry
   does not carry. A seventh item, ``adaptation.cell_wall_crosslinking`` (Gas1-5, Crh1, Utr2), is
   REFUSED: a transglycosylase is an isomerisation of the polymer, so against yeast-GEM's two
   lumped glucan components its reaction is A -> A, a free cycle carrying no flux.

WHAT THE ABLATION MEASURED, which is criterion (a) and is the reason this is not a spreadsheet.
:func:`pathway_essentiality` blocks each row's reactions and reports the growth change against
``generator/panel_experiment.MEASURED_GROWTH_RATE_SE`` = 0.0117 /h. Thirty-one rows carry
reactions the shipped model has; EIGHTEEN of them are ADAPTATION rows that are lethal when
blocked, plus cascade.CWI_RHO1 on the glucan synthase, and the same nineteen are lethal in BOTH
regimes and by the whole of the growth. That regime-independence is the tell: they are lethal as BIOMASS
REQUIREMENTS, not as stress responses. Trehalose, glycogen, chitin and the two glucans are fixed
components of the carbohydrate pseudoreaction r_4048 at 0.13655, 0.35689, 0.02361, 0.73914 and
0.24696 mmol/gDCW; Nha1, Trk1 and Smf1/Smf2 are the only routes by which sodium, potassium and
manganese enter a cell whose ion pseudoreaction demands all three. Of the twelve non-lethal rows, NINE move growth by less
than 1e-7 /h in both regimes -- Pma1, the V-ATPase, cAMP turnover, SOD, glutathione peroxidase,
Tsa1, glycophagy, Ena1 and Tpo1 -- and the largest of the other three is glycerol at
0.000086 /h, 0.0074x the floor. Undosed, this whole table is either lethal or invisible, and
neither of those is a stress response.

ONE ROW CLEARS THE FLOOR, and only under a dose. With ``stress_ros.dose_peroxide`` at 1.0
mmol/gDCW/h under the reference aerobic batch regime, blocking catalase moves growth by
**0.017238 /h = 1.47x** the 0.0117 /h growth SE; at 5.0 it is 0.086187 /h = 7.37x. That is the
audit's own ablation result and it is a conditional one: the pathway is invisible until the
stressor it answers is applied, which is the single most useful thing this table says about how
to test any of the others.

WHAT :func:`apply_pathway_coverage` DOES, AND WHAT IT DELIBERATELY DOES NOT. It writes one
cobra ``Group`` per represented item, so the pathway-to-reaction map lives in the model object
and a caller can ask for "the reactions that implement trehalose accumulation" without going
through an empty regulon. It changes no bound, no gene rule and no stoichiometry, and growth is
identical before and after to the solver's tolerance -- verified, not assumed. It does NOT write
a second E-Flux layer: relating a module's activity to its enzymes' capacity needs a slope, and
``bridge/regulation.py`` already records that no such slope is measured anywhere in this
repository. Inventing one here to make the map look load-bearing is the failure this table exists
to document.

FREE SCALARS, criterion (e). The catalogue loader, validator, group installer, association
query and ablation fit no biological parameter. Supplied models and catalogue records are
inputs, not learned regulatory couplings. The one graded constant this module carries,
:data:`REFUSED_REPRESSION_COUPLING`, counts as one free scalar under ``mech/params.FREE_TAGS``
against zero independent targets -- 1 > 0, the gate fails, and that failure is the reason the
repression arm is declared and not built. ``float()`` on it raises.
"""

from __future__ import annotations

import csv
import importlib
import pathlib
from dataclasses import dataclass, field

import cobra

from ..mech.params import Param
from ..paths import data_dir
from .solver import growth_or_none

__all__ = [
    "COVERAGE_PATH",
    "LAYERS",
    "REFUSED_REPRESSION_COUPLING",
    "UNREACHED_MODULES",
    "VERDICTS",
    "PathwayItem",
    "ModelGeneCoverage",
    "model_gene_coverage",
    "apply_pathway_coverage",
    "coverage_by_layer",
    "gaps",
    "group_name",
    "load_coverage",
    "panel_module_reactions",
    "pathway_essentiality",
    "refusals",
    "represented",
    "unrepresented_pathways",
    "verify_coverage",
]

COVERAGE_PATH = data_dir() / "gem" / "stress_pathway_coverage.tsv"

LAYERS = ("SENSOR", "CASCADE", "TF", "ADAPTATION")

VERDICTS = (
    "IN_GEM", "IN_GEM_RESTORED", "IN_MODULE", "NOT_METABOLIC", "GAP", "REFUSED",
)

#: Verdicts asserting the item has reactions the caller can address.
_WITH_REACTIONS = ("IN_GEM", "IN_GEM_RESTORED")
#: Verdicts asserting the item has no reaction of its own.
_WITHOUT_REACTIONS = ("IN_MODULE", "NOT_METABOLIC", "GAP", "REFUSED")

UNREACHED_MODULES = ("UPR", "calcium", "zinc", "xenobiotic")
"""Panel modules whose regulon reaches zero GEM reactions, measured this session.

Counted by resolving ``bridge.regulation.load_regulons`` against yeast-GEM v9.0.2's gene
list: UPR 7 regulon genes and 0 in the model, calcium 0 genes at all, zinc 6 and 0,
xenobiotic 6 and 0. cell_wall and alkaline_ph reach exactly one reaction each and are not
listed here because one is not zero. ``tests/test_stress_pathways.py`` re-measures it.
"""

REFUSED_REPRESSION_COUPLING = Param.refused(
    name="c_T",
    units="dimensionless (relative enzyme capacity per unit repressor activity)",
    source=(
        "Searched for a Dot6/Tod6 PAC-motif repression coupling with a slope: Lippman & "
        "Broach 2009 (PMID 19901341) and Huber 2011 (PMID 21730963) both establish the "
        "regulon and neither publishes a coupling to metabolic capacity. "
        "REVISED_BUILD_LIST.md section 3 D2 records the same search returning nothing in "
        "any organism-appropriate source."
    ),
    reason=(
        "No published value relates RiBi/RP repression to a reaction bound in any "
        "organism-appropriate source, and this repository's own additive operator cannot "
        "carry one: gene scales add across regulons, so a two-regulon 50% repression "
        "clips an essential reaction to zero and the model dies. It is a knockout "
        "operator, not a dial."
    ),
    missing=(
        "A PAC-motif regulon at a stated evidence level, plus one RiBi transcript "
        "timecourse at a dose that moves it -- the committed RNA is a single post-read "
        "timepoint at 0, 0.2 and 0.5 mM, the three doses where nothing moves."
    ),
)
"""The Dot6/Tod6 repression coupling, REFUSED. ``float()`` raises :class:`RefusedValue`."""

_COLUMNS = (
    "layer", "item", "label", "orfs", "verdict", "gem_reactions",
    "panel_module", "code_module", "basis", "source",
)


@dataclass(frozen=True)
class PathwayItem:
    """One row: a diagram item, a verdict, and the evidence for it."""

    layer: str
    item: str
    label: str
    orfs: str
    verdict: str
    gem_reactions: str
    panel_module: str
    code_module: str
    basis: str
    source: str

    @property
    def orf_ids(self) -> tuple[str, ...]:
        return () if self.orfs == "-" else tuple(self.orfs.split("|"))

    @property
    def reaction_ids(self) -> tuple[str, ...]:
        return () if self.gem_reactions == "-" else tuple(self.gem_reactions.split("|"))

    @property
    def is_gap(self) -> bool:
        """A GAP or REFUSED verdict in the historical catalogue, not a whole-cell gap test."""
        return self.verdict in ("GAP", "REFUSED")


@dataclass(frozen=True)
class ModelGeneCoverage:
    item: str
    model_id: str
    declared_orfs: tuple[str, ...]
    present_orfs: tuple[str, ...]
    absent_orfs: tuple[str, ...]
    reaction_ids: tuple[str, ...]
    associations: tuple[tuple[str, tuple[str, ...]], ...]
    evidence_scope: str = field(default="gene_reaction_associations_only", init=False)


def model_gene_coverage(
    model: cobra.Model, items: tuple[PathwayItem, ...] | None = None,
) -> tuple[ModelGeneCoverage, ...]:
    records = load_coverage() if items is None else items
    if len({record.item for record in records}) != len(records):
        raise ValueError("duplicate pathway item in model coverage audit")
    coverage = []
    for record in records:
        identifiers = tuple(dict.fromkeys(record.orf_ids))
        present = tuple(orf for orf in identifiers if orf in model.genes)
        absent = tuple(orf for orf in identifiers if orf not in model.genes)
        associations = tuple(
            (orf, tuple(sorted(reaction.id for reaction in model.genes.get_by_id(orf).reactions)))
            for orf in present
        )
        reactions = tuple(sorted({reaction for _, attached in associations for reaction in attached}))
        coverage.append(ModelGeneCoverage(
            record.item, model.id, identifiers, present, absent, reactions, associations))
    return tuple(coverage)


def load_coverage(path: pathlib.Path | None = None) -> tuple[PathwayItem, ...]:
    """Read the table, refusing anything self-inconsistent rather than half-applying it.

    Raises:
        ValueError: on an unknown layer or verdict, a duplicated item key, an empty basis,
            or a verdict that disagrees with its own columns -- a row claiming reactions
            while its verdict says it has none, or the reverse.
    """
    source = pathlib.Path(path) if path is not None else COVERAGE_PATH
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(
            (ln for ln in handle if not ln.startswith("#")), delimiter="\t"
        ))
    items: list[PathwayItem] = []
    seen: set[str] = set()
    for row in rows:
        missing = set(_COLUMNS) - set(row)
        if missing:
            raise ValueError(f"{source} is missing columns {sorted(missing)}")
        record = PathwayItem(**{k: (row[k] or "").strip() for k in _COLUMNS})
        _check_row(record, source, seen)
        seen.add(record.item)
        items.append(record)
    return tuple(items)


def _check_row(record: PathwayItem, source: pathlib.Path, seen: set[str]) -> None:
    if record.layer not in LAYERS:
        raise ValueError(f"{record.item}: unknown layer {record.layer!r}")
    if record.verdict not in VERDICTS:
        raise ValueError(f"{record.item}: unknown verdict {record.verdict!r}")
    if record.item in seen:
        raise ValueError(f"{record.item} appears twice in {source}")
    if not record.basis:
        raise ValueError(
            f"{record.item}: basis is required for every verdict. An item with no stated "
            "ground is the padding this table exists to avoid"
        )
    if not record.source:
        raise ValueError(f"{record.item}: source is required, including for a refusal")
    if record.verdict in _WITH_REACTIONS and not record.reaction_ids:
        raise ValueError(
            f"{record.item} is {record.verdict} and names no reaction. That verdict IS the "
            "claim that reactions exist"
        )
    if record.verdict in _WITHOUT_REACTIONS and record.reaction_ids:
        raise ValueError(
            f"{record.item} is {record.verdict} and names reactions {record.reaction_ids}. "
            "Use IN_GEM if the reactions are real"
        )
    named = [v for v in (record.code_module, record.panel_module) if v and v != "-"]
    if record.verdict == "IN_MODULE" and not named:
        raise ValueError(
            f"{record.item} is IN_MODULE and names no module. The verdict is the claim "
            "that some named module carries it"
        )


def coverage_by_layer(
    items: tuple[PathwayItem, ...] | None = None
) -> dict[str, tuple[PathwayItem, ...]]:
    """The table split into the diagram's four layers, in table order."""
    records = load_coverage() if items is None else items
    out: dict[str, tuple[PathwayItem, ...]] = {}
    for layer in LAYERS:
        out[layer] = tuple(r for r in records if r.layer == layer)
    return out


def represented(items: tuple[PathwayItem, ...] | None = None) -> tuple[PathwayItem, ...]:
    """Rows the model or the code actually carries."""
    records = load_coverage() if items is None else items
    return tuple(r for r in records if r.verdict in (*_WITH_REACTIONS, "IN_MODULE"))


def gaps(items: tuple[PathwayItem, ...] | None = None) -> tuple[PathwayItem, ...]:
    """Rows marked GAP or REFUSED in the historical catalogue.

    NOT_METABOLIC is excluded by this classification, not declared biologically dispensable.
    Missing signaling can still be a gap in a whole-cell stress model.
    """
    records = load_coverage() if items is None else items
    return tuple(r for r in records if r.is_gap)


def refusals(items: tuple[PathwayItem, ...] | None = None) -> tuple[PathwayItem, ...]:
    """Rows that deliberately add nothing. Read the ``basis`` column."""
    records = load_coverage() if items is None else items
    return tuple(r for r in records if r.verdict == "REFUSED")


def unrepresented_pathways(
    items: tuple[PathwayItem, ...] | None = None
) -> tuple[PathwayItem, ...]:
    """NOT_METABOLIC rows with neither a panel nor a code-module reference in the table.

    This is a catalogue classification, not a search for every implementation in the
    repository. The historical panel lacks cAMP-PKA (``sensor.GPR1`` and
    ``cascade.PKA_kinases``) and the SPS amino-acid sensor (``sensor.SPS``). Adding a panel
    key could name them, but implementing their signaling requires more than a name.

    IN_GEM rows without a panel reference, including ``sensor.VATPASE`` and
    ``cascade.cAMP_second_messenger``, are excluded by definition. Their historical
    yeast-GEM growth ablations were inert under the tested conditions; that does not prove
    their biological pathways dispensable or their regulation implemented. REFUSED rows
    remain available through :func:`refusals` with the missing evidence named.
    """
    records = load_coverage() if items is None else items
    return tuple(
        r for r in records
        if r.verdict == "NOT_METABOLIC"
        and not [v for v in (r.panel_module, r.code_module) if v and v != "-"]
    )


def panel_module_reactions(
    model: cobra.Model, items: tuple[PathwayItem, ...] | None = None
) -> dict[str, frozenset[str]]:
    """Panel module -> the GEM reactions this table declares for it, present in ``model``.

    This is the map that gives ``UNREACHED_MODULES`` reactions. It is curated and cited row
    by row rather than derived from a motif map, and it is deliberately NOT a regulon: it
    carries no direction, no weight and no slope, because none of those is measured.
    """
    records = load_coverage() if items is None else items
    out: dict[str, set[str]] = {}
    for record in records:
        if not record.panel_module or record.panel_module == "-":
            continue
        present = {r for r in record.reaction_ids if model.reactions.has_id(r)}
        if present:
            out.setdefault(record.panel_module, set()).update(present)
    return {k: frozenset(v) for k, v in out.items()}


def verify_coverage(
    model: cobra.Model,
    items: tuple[PathwayItem, ...] | None = None,
    restored: cobra.Model | None = None,
) -> tuple[str, ...]:
    """Re-check every claim in the table against a real model. Empty is the pass.

    Absence is checked by ORF against the model's own gene list, never by symbol: cobrapy
    leaves ``Gene.name`` empty, so a symbol lookup returns nothing for every gene and the
    whole diagram would read as absent. That trap is why two rows here exist at all -- GPX3
    is in the model as HYR1 and RHO1 is in it on the glucan synthase.

    Args:
        model: The shipped yeast-GEM, unmodified.
        items: Table rows. Defaults to the file.
        restored: The same model after ``stress_genes.restore_stress_genes``. Required to
            check IN_GEM_RESTORED rows; without it they are skipped rather than passed, and
            the skip is reported.

    Returns:
        One string per failed claim, naming the row and what was expected.
    """
    records = load_coverage() if items is None else items
    problems: list[str] = []
    for record in records:
        if record.verdict == "IN_GEM":
            problems += _check_in_gem(model, record)
        elif record.verdict == "IN_GEM_RESTORED":
            problems += _check_restored(model, restored, record)
        else:
            problems += _check_absent(model, record)
        problems += _check_modules(record)
    return tuple(problems)


def _check_in_gem(model: cobra.Model, record: PathwayItem) -> list[str]:
    out = [
        f"{record.item}: IN_GEM names reaction {rid}, which {model.id!r} does not have"
        for rid in record.reaction_ids if not model.reactions.has_id(rid)
    ]
    out += [
        f"{record.item}: IN_GEM names ORF {orf}, which {model.id!r} does not have. Either "
        "the row is wrong or the gene belongs in absent_stress_genes.tsv"
        for orf in record.orf_ids if not model.genes.has_id(orf)
    ]
    return out


def _check_restored(
    model: cobra.Model, restored: cobra.Model | None, record: PathwayItem
) -> list[str]:
    if restored is None:
        return [f"{record.item}: IN_GEM_RESTORED not checked, no restored model was given"]
    out = [
        f"{record.item}: IN_GEM_RESTORED names reaction {rid}, absent even after restore"
        for rid in record.reaction_ids if not restored.reactions.has_id(rid)
    ]
    out += [
        f"{record.item}: IN_GEM_RESTORED names ORF {orf}, absent even after restore"
        for orf in record.orf_ids if not restored.genes.has_id(orf)
    ]
    before = [o for o in record.orf_ids if not model.genes.has_id(o)]
    before += [r for r in record.reaction_ids if not model.reactions.has_id(r)]
    if not before:
        out.append(
            f"{record.item}: IN_GEM_RESTORED, but every ORF and reaction it names is "
            "already in the shipped model. It is IN_GEM"
        )
    return out


def _check_absent(model: cobra.Model, record: PathwayItem) -> list[str]:
    """A NOT_METABOLIC or REFUSED row whose gene turns out to be in the model is the one
    error this table cannot survive: it would be padding the gap list with a real gene."""
    return [
        f"{record.item}: {record.verdict} names ORF {orf}, which {model.id!r} DOES have "
        f"(label {model.genes.get_by_id(orf).name or '?'}). Absence claimed by symbol is "
        "the trap this table warns about"
        for orf in record.orf_ids if model.genes.has_id(orf)
    ]


def _check_modules(record: PathwayItem) -> list[str]:
    out: list[str] = []
    if record.panel_module and record.panel_module != "-":
        from ..generator.stress_panel import MODULES
        if record.panel_module not in MODULES:
            out.append(
                f"{record.item}: panel_module {record.panel_module!r} is not a key of "
                "generator.stress_panel.MODULES"
            )
    if record.code_module and record.code_module != "-":
        try:
            importlib.import_module(f"ystwin.{record.code_module}")
        except ImportError as exc:
            out.append(f"{record.item}: code_module {record.code_module!r} does not import ({exc})")
    return out


def group_name(record: PathwayItem) -> str:
    """The cobra ``Group`` id a row installs. Derived, so the guard cannot drift from it."""
    return f"stress_pathway__{record.item}"


def apply_pathway_coverage(
    model: cobra.Model, items: tuple[PathwayItem, ...] | None = None
) -> cobra.Model:
    """Return a copy of ``model`` carrying one ``Group`` per row whose reactions it has.

    This is the declaration made operative: after it, ``model.groups`` answers "which
    reactions implement trehalose accumulation" directly, which the regulon route cannot do
    for the four modules in :data:`UNREACHED_MODULES`.

    Nothing else changes. No bound, no gene rule, no stoichiometry and no objective is
    touched, so the optimum is identical -- ``tests/test_stress_pathways.py`` asserts that
    rather than trusting it. Rows naming a reaction the model lacks are skipped, so this is
    safe to call on the shipped model and on a restored one alike.

    Raises:
        ValueError: if the groups are already installed.
    """
    records = load_coverage() if items is None else items
    installed = [
        group_name(r) for r in records
        if r.reaction_ids and any(g.id == group_name(r) for g in model.groups)
    ]
    if installed:
        raise ValueError(f"stress-pathway groups already installed: {installed}")

    out = model.copy()
    for record in records:
        present = [
            out.reactions.get_by_id(rid) for rid in record.reaction_ids
            if out.reactions.has_id(rid)
        ]
        if not present:
            continue
        group = cobra.core.Group(group_name(record), name=record.label, members=present)
        group.kind = "collection"
        out.add_groups([group])
    return out


def pathway_essentiality(
    model: cobra.Model, items: tuple[PathwayItem, ...] | None = None
) -> dict[str, float]:
    """Item -> the growth lost when its reactions are blocked, /h. Does not modify ``model``.

    Criterion (a) for this table, and the floor to score it against is
    ``generator.panel_experiment.MEASURED_GROWTH_RATE_SE`` = 0.0117 /h, which is the
    standard error of a growth rate estimated from a real optical-density trace. Do not
    borrow the plate reader's activity CV for this: it measures a different quantity.

    An infeasible or unbounded solve is reported as the whole of the growth, which is what
    lethal means here. Rows with no reaction in ``model`` are absent from the result rather
    than reported as zero, so a missing row cannot read as a measured null.
    """
    records = load_coverage() if items is None else items
    base = growth_or_none(model) or 0.0
    out: dict[str, float] = {}
    for record in records:
        present = [rid for rid in record.reaction_ids if model.reactions.has_id(rid)]
        if not present:
            continue
        with model as scoped:
            for rid in present:
                scoped.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            value = growth_or_none(scoped)
        out[record.item] = base - (value or 0.0)
    return out
