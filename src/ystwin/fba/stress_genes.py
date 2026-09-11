"""Restore the stress genes yeast-GEM v9.0.2 does not have, from a declared table.

The table is ``data/gem/absent_stress_genes.tsv`` and it, not this file, is the work. Each of
its 41 rows carries a GO molecular-function term, the SGD evidence line behind it, and -- for
the 20 that were REFUSED -- why adding the gene would have been wrong.

Four verdicts. GPR_EXTEND rewrites the rule of a reaction that already exists and is the cheap
and correct fix (8 rows). NEW_REACTION adds chemistry the model genuinely lacks (11 rows over 2
reactions). REFUSED and ALREADY_PRESENT change nothing and are carried so the next reader does
not re-litigate them.

What this repairs is a FALSE NEGATIVE, not a flux: a regulon naming TSA1 scaled nothing, because
YML028W was not a gene in the model. Two repairs also move a prediction -- r_4589/r_4587 carried
only SMF1 and r_4590 only ALR1, so deleting either read as exactly lethal.

That reach is narrow, and saying so is part of the claim. Of the 19 genes added, FOUR are named
by any module regulon -- TSA1, CTR1, FTR1, MNN1 -- and three more are load-bearing for a deletion
call rather than a regulon: SMF2 carries the alr1 correction and MID1 with CCH1 carry half of the
smf1 one. Measured gene by gene, the remaining TWELVE move no single-gene observable at all --
no regulon scale, no deletion call, no flux -- and ALG3 is the extreme case: ALG3_ER measured FVA
[0, 0] in every regime probed and no regulon names YBL082C. They are carried because the gaps are
real, not because they move a number, and that is the honest reading of criterion (a).

The refusals have a cost, and it is four regulon genes that stay out of reach: GAS5 in UPR,
CRH1 in cell_wall, ERG28 and CYB5 in hypoxia are all named by a module regulon and all correctly
refused, so those modules are still short a gene each. CCS1 is named by the oxidative regulon and
costs it nothing, because SOD1 (YJR104C) is already in the model on r_4270 and in that regulon.

The prediction ledger is two corrected and one broken, and it is the whole genome, not a spot
check: a full single-gene deletion sweep moves exactly three calls and takes the lethal count
from 188 to 185. smf1 and alr1 stop reading lethal, which is right, and dpm1 starts reading
viable, which is wrong. The two corrections come from the GPR_EXTEND rows alone and the one
regression from MANPOL_G, so the halves are separable.

Measured through the repository's own gene-level layer, ``bridge/regulation.eflux_layer``: at
oxidative activity +0.5 the YAP1/SKN7 regulon scales r_0550 by 1.000 before this module and 1.500
after, and the CUP2/MAC1 copper regulon reaches r_4589 only after. The layer's own count of
unreachable regulon genes falls with it -- 101 absent to 99 on the oxidative module.

Absence is established by ORF against the SBML's own ``fbc:geneProduct`` labels, never by symbol:
SGD resolves the bare string "CTR1" to YGL077C (HNM1), which IS in the model.

Measured under the model's shipped default bounds: growth 0.085844 /h before and after. ALG3_ER
carries zero flux and MANPOL_G flips DPM1 from lethal to viable; both are declared in the table.

And these additions are not a free lunch, which is measured rather than argued. With every
exchange opened maximum growth is 77.556090 /h before and after; on a closed medium neither model
can drain ATP, NADPH or NADH, and both new reactions measure FVA [0, 0] there. MANPOL_G cannot be
a shortcut because its net is the route it parallels: r_1803 + MANPOL_G + r_1801 and
r_0361 + r_0362 + r_1932, with r_1748 returning the dolichol carrier, both spend one
GDP-mannose[c] per mannan, differing only in which compartment the proton lands in, and r_1826
moves that for free. Measured rather than argued again: a whole-model FVA before and after moves
four ranges, and only r_1801 and r_1803 by more than 1e-7 -- each to the new route's own
0.060266 mmol/gDCW/h -- and no reaction becomes unbounded that was not unbounded already.
"""

from __future__ import annotations

import csv
import gzip
import pathlib
import re
from dataclasses import dataclass

import cobra

from ..paths import data_dir

__all__ = [
    "DECLARATION_PATH",
    "NEW_REACTION_IDS",
    "SGD_GAF_PATH",
    "VERDICTS",
    "Declaration",
    "additions",
    "gene_label_map",
    "load_declarations",
    "new_reaction_targets",
    "refusals",
    "restore_stress_genes",
    "sgd_annotations",
    "unsupported_claims",
    "verify_absence",
]

DECLARATION_PATH = data_dir() / "gem" / "absent_stress_genes.tsv"
SGD_GAF_PATH = data_dir() / "external" / "sgd" / "gene_association.sgd.gaf.gz"
NEW_REACTION_IDS = ("MANPOL_G", "ALG3_ER")
VERDICTS = ("GPR_EXTEND", "NEW_REACTION", "REFUSED", "ALREADY_PRESENT")

_GENE_PRODUCT = re.compile(r'fbc:id="([^"]+)"\s+fbc:label="([^"]+)"')
_COLUMNS = (
    "orf", "symbol", "verdict", "target", "reaction_name", "stoichiometry",
    "bounds", "old_rule", "new_rule", "go_term", "evidence", "basis",
)


@dataclass(frozen=True)
class Declaration:
    """One row of the table: a gene, a verdict, and the evidence for it."""

    orf: str
    symbol: str
    verdict: str
    target: str
    reaction_name: str
    stoichiometry: str
    bounds: str
    old_rule: str
    new_rule: str
    go_term: str
    evidence: str
    basis: str

    @property
    def parsed_bounds(self) -> tuple[float, float]:
        low, high = self.bounds.split(",")
        return float(low), float(high)

    @property
    def parsed_stoichiometry(self) -> dict[str, float]:
        pairs = (item.split(":") for item in self.stoichiometry.split("|"))
        return {mid: float(coef) for mid, coef in pairs}


def gene_label_map(sbml_path: pathlib.Path) -> dict[str, str]:
    """ORF -> gene symbol, read straight out of the SBML. Reads ``.xml`` or ``.xml.gz``.

    cobrapy leaves ``Gene.name`` empty, so a model object cannot answer "is TSA1 here?".
    """
    opener = gzip.open if str(sbml_path).endswith(".gz") else open
    with opener(sbml_path, "rt") as handle:
        return dict(_GENE_PRODUCT.findall(handle.read()))


def load_declarations(path: pathlib.Path | None = None) -> tuple[Declaration, ...]:
    """Read the table, refusing anything self-inconsistent rather than half-applying it.

    Raises:
        ValueError: on an unknown verdict, a duplicated ORF, or two rows that name the same
            new reaction and disagree about its chemistry.
    """
    source = pathlib.Path(path) if path is not None else DECLARATION_PATH
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(
            (ln for ln in handle if not ln.startswith("#")), delimiter="\t"
        ))
    declarations = []
    seen_orfs: set[str] = set()
    for row in rows:
        missing = set(_COLUMNS) - set(row)
        if missing:
            raise ValueError(f"{source} is missing columns {sorted(missing)}")
        record = Declaration(**{k: (row[k] or "").strip() for k in _COLUMNS})
        if record.verdict not in VERDICTS:
            raise ValueError(f"{record.orf}: unknown verdict {record.verdict!r}")
        if record.orf in seen_orfs:
            raise ValueError(f"{record.orf} appears twice in {source}")
        seen_orfs.add(record.orf)
        declarations.append(record)
    _check_reaction_rows_agree(declarations, source)
    return tuple(declarations)


def _check_reaction_rows_agree(declarations, source) -> None:
    """Rows sharing a new reaction repeat its definition; a disagreement is an error, not a vote."""
    by_target: dict[str, list[Declaration]] = {}
    for record in declarations:
        if record.verdict == "NEW_REACTION":
            by_target.setdefault(record.target, []).append(record)
    for target, group in by_target.items():
        shapes = {(r.reaction_name, r.stoichiometry, r.bounds, r.new_rule) for r in group}
        if len(shapes) != 1:
            raise ValueError(f"{source}: rows for {target} disagree about its definition")


def additions(declarations: tuple[Declaration, ...] | None = None) -> tuple[Declaration, ...]:
    """The rows that change the model."""
    records = load_declarations() if declarations is None else declarations
    return tuple(r for r in records if r.verdict in ("GPR_EXTEND", "NEW_REACTION"))


def refusals(declarations: tuple[Declaration, ...] | None = None) -> tuple[Declaration, ...]:
    """The rows that deliberately change nothing. Read the ``basis`` column."""
    records = load_declarations() if declarations is None else declarations
    return tuple(r for r in records if r.verdict == "REFUSED")


def new_reaction_targets(
    declarations: tuple[Declaration, ...] | None = None
) -> tuple[str, ...]:
    """Reaction ids the table adds, in table order. Derived, so the guard cannot drift from it."""
    records = load_declarations() if declarations is None else declarations
    return tuple(dict.fromkeys(r.target for r in records if r.verdict == "NEW_REACTION"))


def sgd_annotations(
    path: pathlib.Path | None = None
) -> dict[str, set[tuple[str, str, str]]] | None:
    """ORF -> {(GO term, ``code:reference``, aspect)} from the SGD GAF, or ``None`` if absent.

    The GAF is gitignored, so callers skip on ``None`` rather than failing. SGD writes the ORF
    in the synonym column and the symbol in the object column, so both are indexed.

    A NOT-qualified line is dropped rather than indexed. SGD uses it to record that a gene does
    NOT have a function, and indexing it would let a row cite, as its evidence, the annotation
    that refutes it.
    """
    source = pathlib.Path(path) if path is not None else SGD_GAF_PATH
    if not source.exists():
        return None
    annotations: dict[str, set[tuple[str, str, str]]] = {}
    with gzip.open(source, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15 or "NOT" in cols[3].split("|"):
                continue
            entry = (cols[4], f"{cols[6]}:{cols[5]}", cols[8])
            for key in {cols[1], cols[2], *(s for s in cols[10].split("|") if s)}:
                annotations.setdefault(key, set()).add(entry)
    return annotations


def unsupported_claims(
    declarations: tuple[Declaration, ...] | None = None,
    annotations: dict[str, set[tuple[str, str, str]]] | None = None,
) -> tuple[tuple[str, str, str], ...]:
    """Rows whose GO term and evidence line are not an aspect-F annotation in the GAF.

    Empty is the pass, and it is how the one invented citation in this table was found: SMF2's
    manganese term is IBA in SGD, not the IDA this file first claimed.
    """
    records = load_declarations() if declarations is None else declarations
    table = sgd_annotations() if annotations is None else annotations
    if table is None:
        raise FileNotFoundError(f"SGD GAF not found at {SGD_GAF_PATH}")
    return tuple(
        (r.orf, r.go_term, r.evidence)
        for r in records
        if r.verdict != "ALREADY_PRESENT"
        and (r.go_term, r.evidence, "F") not in table.get(r.orf, set())
    )


def verify_absence(
    model: cobra.Model, declarations: tuple[Declaration, ...] | None = None
) -> tuple[str, ...]:
    """ORFs the table calls absent that ``model`` actually carries. Empty is the pass.

    ALREADY_PRESENT rows are excluded, since those are the ones asserted to be there.
    """
    records = load_declarations() if declarations is None else declarations
    return tuple(
        r.orf for r in records
        if r.verdict != "ALREADY_PRESENT" and model.genes.has_id(r.orf)
    )


def restore_stress_genes(
    model: cobra.Model, declarations: tuple[Declaration, ...] | None = None
) -> cobra.Model:
    """Return a copy of ``model`` with the declared GPRs extended and reactions added.

    The input model is not modified. Every GPR_EXTEND row asserts the reaction's current rule
    against the table's ``old_rule`` first, so a yeast-GEM whose annotation has moved on
    refuses the edit instead of silently overwriting it.

    Raises:
        ValueError: if the additions are already installed, or a rule has drifted.
        KeyError: naming a reaction or metabolite the host model does not carry.
    """
    records = load_declarations() if declarations is None else declarations
    installed = [rid for rid in new_reaction_targets(records) if model.reactions.has_id(rid)]
    if installed:
        raise ValueError(f"stress-gene additions already installed: {installed}")

    out = model.copy()
    for record in records:
        if record.verdict == "GPR_EXTEND":
            _extend_gpr(out, record)
    for target, group in _group_new_reactions(records).items():
        _add_declared_reaction(out, target, group[0])
    return out


def _extend_gpr(model: cobra.Model, record: Declaration) -> None:
    if not model.reactions.has_id(record.target):
        raise KeyError(f"{record.symbol}: {model.id!r} has no reaction {record.target}")
    reaction = model.reactions.get_by_id(record.target)
    # Two genes can share a target (CTR1/CTR3, MID1/CCH1), so the second row sees new_rule already.
    if reaction.gene_reaction_rule not in (record.old_rule, record.new_rule):
        raise ValueError(
            f"{record.target}: rule is {reaction.gene_reaction_rule!r}, table expected "
            f"{record.old_rule!r}. Refusing to overwrite an annotation that has moved."
        )
    reaction.gene_reaction_rule = record.new_rule


def _group_new_reactions(records) -> dict[str, list[Declaration]]:
    grouped: dict[str, list[Declaration]] = {}
    for record in records:
        if record.verdict == "NEW_REACTION":
            grouped.setdefault(record.target, []).append(record)
    return grouped


def _add_declared_reaction(model: cobra.Model, target: str, record: Declaration) -> None:
    stoichiometry = record.parsed_stoichiometry
    unknown = [mid for mid in stoichiometry if not model.metabolites.has_id(mid)]
    if unknown:
        raise KeyError(f"{target}: {model.id!r} has no metabolites {unknown}")
    low, high = record.parsed_bounds
    reaction = cobra.Reaction(target, name=record.reaction_name, lower_bound=low, upper_bound=high)
    reaction.subsystem = "restored stress genes"
    model.add_reactions([reaction])
    reaction.add_metabolites(
        {model.metabolites.get_by_id(mid): coef for mid, coef in stoichiometry.items()}
    )
    reaction.gene_reaction_rule = record.new_rule
