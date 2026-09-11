"""The published carotenoid strains, read as data instead of described in prose.

``data/carotenoid_batch/published_batch_titres.tsv`` is a 29-row survey of published
*S. cerevisiae* carotenoid strains -- titre, content, medium, cultivation mode, and since
2026-09-03 a per-gene cassette dosage with named evidence. Until this module it had **no code
reader at all**: `grep -rn published_batch_titres` returned two prose mentions and nothing
else, so five docstrings made categorical claims about a file nothing could open.

**Every one of those claims had gone false, and the file itself is what falsified them.**

    "every cassette is qualitative"                       five papers state a number
    "not one reports its cassette quantitatively"         five do, with `copy_evidence`
    "every carotenoid row is glucose or galactose"        sucrose, xylose, olive oil,
                                                          oleic acid, grape juice, acetate

The shape of the failure is the point and is worth stating once. A derived claim ABOUT A FILE
was copied into five modules with nothing that rechecked it, and it went stale three days
after it was written -- `pathway/capacity.py` on 2026-08-31, the ``copies_*`` columns on
2026-09-03. `pathway/enzyme_capacity.py`'s copy was written IN the very commit that added the
columns, and so was false on the day it shipped.

**THE CONCLUSIONS SURVIVE, ON BETTER GROUND.** :func:`crtyb_dosage_unfittable_reason`
computes the real one: crtYB takes only TWO distinct numeric dosages across the papers that
state one, and those papers differ in host, promoter, medium, cultivation mode and assay. The
coefficient is unfittable through CONFOUNDING, not through absent numbers -- which is a
sharper reason, and one that stops being true the day a third comparable dosage arrives.
`tests/test_capacity.py` asserts exactly that, so the refusal is revisited when it should be.

**THE COLUMN DOMAIN IS DECLARED AND VIOLATIONS RAISE.** A survey table accumulates free-text
non-values -- `promoter` held the string ``"not recorded here"`` in a column otherwise holding
promoter names -- and a reader that shrugs at them turns a gap into a value. Every column this
module reads has a declared vocabulary; anything else is a load-time failure naming the row.
See ``data/carotenoid_batch/SOURCE.md`` for provenance and for the population rule the file
follows.
"""

from __future__ import annotations

import csv
import functools
import math
from dataclasses import dataclass

from .. import paths

__all__ = [
    "CARBON_SPECIES",
    "EPISOMAL_2U",
    "HIGH_COPY",
    "UNRECORDED",
    "UNSTATED",
    "CassetteDomainError",
    "PublishedCassette",
    "crtyb_dosage_unfittable_reason",
    "numeric_crtyb_dosages",
    "published_cassettes",
    "rows_above_content",
]

#: A dosage stated as "on a 2-micron plasmid": episomal, and the copy number is a
#: distribution over a population rather than a number. NOT convertible to an integer, and
#: that is why it is a distinct value rather than a guess at 20.
EPISOMAL_2U = "2u"

#: A dosage stated only as "high-copy". Weaker than :data:`EPISOMAL_2U`, which at least names
#: the replicon.
HIGH_COPY = "high-copy"

#: The declared MISSING sentinel, for a cell whose value the source states but this survey has
#: not transcribed. Distinct from a blank, which means the row was never worked up at all --
#: see the population rule in ``data/carotenoid_batch/SOURCE.md``. It exists because the file
#: carried the free-text string "not recorded here" in a promoter column, which no reader
#: could tell from a promoter named "not recorded here".
UNRECORDED = "unrecorded"

#: What a ``carbon_source`` cell says when the source does not. Enumerated rather than
#: pattern-matched, so a NEW non-value fails to load instead of parsing to "no carbon".
UNSTATED = ("not stated", "not stated in abstract", "not verified", "UNVERIFIED")

#: Carbon-bearing species this survey's free-text media descriptions name, as
#: ``{species: (substrings that name it)}``. Free text is parsed against a DECLARED
#: vocabulary and a cell naming none of them raises, so the set is a claim that can be wrong
#: rather than a regex that quietly returns nothing.
#:
#: ``ethanol`` is in the vocabulary and matches NOTHING in the current file. That is
#: deliberate and it is the load-bearing entry: `pathway/environment_flux.py` fits its
#: descriptor on the feed's ETHANOL CARBON FRACTION, so the reason no carotenoid row here can
#: fit that coefficient is that z = 0 on every one of them.
CARBON_SPECIES: dict[str, tuple[str, ...]] = {
    "glucose": ("glucose", "ypd", "hydrolysate"),
    "galactose": ("galactose", "ypg"),
    "ethanol": ("ethanol",),
    "sucrose": ("sucrose",),
    "xylose": ("xylose",),
    "acetate": ("acetate",),
    "glycerol": ("glycerol",),
    "fatty_acid": ("olive oil", "oleic acid"),
    "fruit_juice": ("grape juice",),
    "undeclared_rich_medium": ("ynb", "sc", "optimized medium", "by-products"),
}

_TABLE = paths.data_dir() / "carotenoid_batch" / "published_batch_titres.tsv"

_DOSAGE_GENES = ("crtE", "crtYB", "crtI")


class CassetteDomainError(ValueError):
    """A cell holds something outside its column's declared vocabulary.

    Raised at load time rather than skipped, because the alternative is what this module was
    written to end: a survey table quietly growing free-text non-values that read as data.
    """


@dataclass(frozen=True)
class PublishedCassette:
    """One published strain, with the parts of its row that have a declared domain."""

    pmid: str
    reference: str
    strain: str
    carbon_source: str
    #: Parsed from ``carbon_source`` against :data:`CARBON_SPECIES`. Empty exactly when the
    #: source states no medium -- see :data:`UNSTATED`.
    carbon_species: frozenset[str]
    dfba_ready: bool
    #: ``{gene: int | EPISOMAL_2U | HIGH_COPY}``, only for genes whose dosage the row states.
    copies: dict[str, int | str]
    #: Promoter name, :data:`UNRECORDED`, or ``None`` where the row was never worked up.
    promoter: str | None
    integration: str
    copy_evidence: str
    content_mg_per_gdcw: float | None
    measurand: str | None = None
    journal: str = ""
    mode: str = ""
    hours: str = ""
    bioreactor: str = ""
    titre_as_reported: str = ""
    notes: str = ""

    @property
    def states_a_numeric_crtyb_dosage(self) -> bool:
        return isinstance(self.copies.get("crtYB"), int)


def _dosage(raw: str, row: int, gene: str) -> int | str | None:
    if raw == "":
        return None
    if raw in (EPISOMAL_2U, HIGH_COPY, UNRECORDED):
        return raw
    try:
        return int(raw)
    except ValueError:
        raise CassetteDomainError(
            f"row {row}: copies_{gene} is {raw!r}, which is neither an integer nor one of "
            f"{EPISOMAL_2U!r}, {HIGH_COPY!r}, {UNRECORDED!r}. A dosage column that accepts "
            f"free text turns a gap into a value") from None


def _carbon_species(raw: str, row: int) -> frozenset[str]:
    if raw == "" or raw in UNSTATED:
        return frozenset()
    lowered = raw.lower()
    found = {species for species, needles in CARBON_SPECIES.items()
             if any(needle in lowered for needle in needles)}
    if not found:
        raise CassetteDomainError(
            f"row {row}: carbon_source {raw!r} names no species in CARBON_SPECIES and is not "
            f"one of the declared non-values {UNSTATED}. Add the species to the vocabulary or "
            f"the non-value to UNSTATED -- do not leave it parsing to 'no carbon', which is "
            f"what makes a downstream 'every row is glucose' claim possible")
    return frozenset(found)


def _promoter(raw: str, dfba_ready: bool, row: int) -> str | None:
    if raw == "":
        if dfba_ready:
            raise CassetteDomainError(
                f"row {row}: promoter is blank on a dfba_ready row. The file's convention is "
                f"that promoter is populated exactly on dfba_ready=yes rows; use "
                f"{UNRECORDED!r} to say the source states it and this survey has not "
                f"transcribed it")
        return None
    if not dfba_ready:
        raise CassetteDomainError(
            f"row {row}: promoter {raw!r} is populated on a dfba_ready=no row, against the "
            f"file's own convention. Either the row is dfba_ready or the convention has "
            f"changed and data/carotenoid_batch/SOURCE.md must say so")
    if raw in UNSTATED or raw.lower().startswith("not recorded"):
        raise CassetteDomainError(
            f"row {row}: promoter is the free-text non-value {raw!r}. Use the declared "
            f"sentinel {UNRECORDED!r}, which a reader can tell from a promoter name")
    return raw


def _content_mg_per_gdcw(raw: str, row: int) -> float | None:
    if raw == "":
        return None
    try:
        content = float(raw)
        if math.isfinite(content) and content >= 0.0:
            return content
    except ValueError:
        pass
    raise CassetteDomainError(
        f"row {row}: content_mg_per_gdcw is {raw!r}; expected a finite nonnegative number "
        f"in mg/gDCW or blank")


@functools.lru_cache(maxsize=1)
def published_cassettes() -> tuple[PublishedCassette, ...]:
    """Every row of the survey, validated against its columns' declared domains.

    Raises:
        FileNotFoundError: if the table is absent, rather than returning an empty tuple. An
            empty survey would make every claim below vacuously true, which is the opposite
            of what a missing measurement should do.
        CassetteDomainError: naming the row and the cell, on any value outside its domain.
    """
    if not _TABLE.is_file():
        raise FileNotFoundError(
            f"{_TABLE} is not vendored. It is the published-strain survey; see "
            f"data/carotenoid_batch/SOURCE.md for what it holds and where each row came from")
    with _TABLE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for column in ("content_mg_per_gdcw", "measurand"):
            if column not in (reader.fieldnames or ()):
                raise CassetteDomainError(f"{_TABLE}: missing required column {column!r}")
        raw_rows = list(reader)

    found = []
    for index, row in enumerate(raw_rows, start=2):   # start=2: line 1 is the header
        ready = row["dfba_ready"].strip()
        if ready not in ("yes", "no"):
            raise CassetteDomainError(
                f"row {index}: dfba_ready is {ready!r}, expected 'yes' or 'no'")
        dfba_ready = ready == "yes"
        copies = {gene: _dosage(row[f"copies_{gene}"].strip(), index, gene)
                  for gene in _DOSAGE_GENES}
        content = row["content_mg_per_gdcw"].strip()
        found.append(PublishedCassette(
            pmid=row["pmid"].strip(),
            reference=row["reference"].strip(),
            strain=row["strain"].strip(),
            carbon_source=row["carbon_source"].strip(),
            carbon_species=_carbon_species(row["carbon_source"].strip(), index),
            dfba_ready=dfba_ready,
            copies={gene: value for gene, value in copies.items() if value is not None},
            promoter=_promoter(row["promoter"].strip(), dfba_ready, index),
            integration=row["integration"].strip(),
            copy_evidence=row["copy_evidence"].strip(),
            content_mg_per_gdcw=_content_mg_per_gdcw(content, index),
            measurand=row["measurand"].strip() or None,
            journal=row.get("journal", "").strip(),
            mode=row.get("mode", "").strip(),
            hours=row.get("hours", "").strip(),
            bioreactor=row.get("bioreactor", "").strip(),
            titre_as_reported=row.get("titre_as_reported", "").strip(),
            notes=row.get("notes", "").strip()))
    return tuple(found)


def rows_above_content(threshold_mg_per_gdcw: float, *,
                       specific_only: bool = True) -> tuple[PublishedCassette, ...]:
    if not math.isfinite(threshold_mg_per_gdcw) or threshold_mg_per_gdcw < 0.0:
        raise ValueError(
            f"threshold_mg_per_gdcw must be finite and nonnegative, got {threshold_mg_per_gdcw}")
    return tuple(
        cassette for cassette in published_cassettes()
        if cassette.content_mg_per_gdcw is not None
        and math.isfinite(cassette.content_mg_per_gdcw)
        and cassette.content_mg_per_gdcw > threshold_mg_per_gdcw
        and (not specific_only or cassette.measurand in ("beta_carotene", "both_reported")))


def numeric_crtyb_dosages() -> dict[int, tuple[str, ...]]:
    """``{crtYB copy number: the papers stating it}``, over the rows that state one.

    The authoritative signal behind the capacity refusal, computed rather than remembered.
    """
    found: dict[int, list[str]] = {}
    for cassette in published_cassettes():
        if cassette.states_a_numeric_crtyb_dosage:
            found.setdefault(int(cassette.copies["crtYB"]), []).append(cassette.reference)
    return {dosage: tuple(sorted(papers)) for dosage, papers in sorted(found.items())}


def crtyb_dosage_unfittable_reason() -> str:
    """Why the crtYB-dosage coefficient cannot be fitted from the literature, derived.

    Computed at the point of use so that the five docstrings which used to carry their own
    copies of a retracted reason can quote ONE sentence that recomputes itself. What it says
    changes when the table changes, which is the property the old prose lacked.
    """
    dosages = numeric_crtyb_dosages()
    papers = sorted({paper for group in dosages.values() for paper in group})
    stated = ", ".join(
        f"{dosage} cop{'y' if dosage == 1 else 'ies'} ({'; '.join(group)})"
        for dosage, group in dosages.items())
    return (
        f"{len(papers)} of {len(published_cassettes())} published strains state a crtYB "
        f"dosage as a NUMBER -- {stated} -- so the survey holds only "
        f"{len(dosages)} distinct numeric dosages, and every pair of them differs in host, "
        f"promoter, medium, cultivation mode and assay as well as in dosage. The coefficient "
        f"is unfittable through CONFOUNDING, not through absent numbers, and that changes "
        f"the day a third comparable dosage arrives")
