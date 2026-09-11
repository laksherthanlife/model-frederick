"""Find every literature identifier this repository cites, in code and in prose.

One scanner, imported by both ``scripts/refresh_citations.py`` and
``tests/test_citations.py``. It used to be two copies of a regex and a
``_logical_text`` helper, one in each file, and that is the wrong shape for the thing
the guard is protecting: a scanner that misses a citation makes the citation invisible
to the table, and two scanners can drift apart so that the refresh vendors an
identifier the test never asks about, or the reverse. There is one definition here.

What "cites" means is deliberately narrow: an identifier is a *claim of provenance*, and
this module's job is to find every one of them so that something else can check where it
points. It does not judge whether the paper supports the sentence. Nothing here can.

Four identifier kinds, because documentation cites things code does not:

``pmid``
    PubMed. Needs a label -- a bare seven-digit number is indistinguishable from a page
    range or an accession -- so it is found by scanning forward from ``PMID``/``PMIDs``,
    and from the cells of a markdown table whose header names a PMID column.
``pmcid``
    PubMed Central. Self-identifying (``PMC`` prefix), so no label scan is needed. Worth
    resolving separately from PMIDs rather than folded into them: the one PMCID error
    this repository has had was a PMCID whose *PMID mapping* was the tell.
``doi``
    Self-identifying by shape (``10.<registrant>/<suffix>``). Most documentation
    citations here carry a DOI and no PMID at all -- conference proceedings, statistics
    journals, book chapters -- so a PMID-only guard covers none of them.
``arxiv``
    Preprints. Two traps this repository has already met in the wild: an arXiv title
    that differs from the published title, so citing "the paper" by its published name
    and its arXiv identifier reads as consistent when the two are different documents;
    and a result that changed between preprint versions, so an unversioned identifier
    does not pin the formula being cited.

Markdown hides citations in ways Python does not
------------------------------------------------
The Python surface had one such trick and it is guarded: implicit string concatenation
splits a citation across a line continuation, so a scanner reading the file rather than
the value sees a quote where the digits should be. :func:`logical_text` rejoins adjacent
string literals before matching.

Markdown's version is the **soft line break**. A reference list wraps, and the label ends
one line while the digits begin the next::

    2021;87(14):e0030121. doi:[10.1128/AEM.00301-21](...). PMID
    34047633. PMC8276805. -- A third `HAC1` primer set, ...

The old pattern was ``PMID[: ]+(\\d{6,8})``, and ``[: ]`` does not match a newline, so
that citation -- and one more like it in the same file -- was invisible to the guard
while being perfectly legible to a reader. The fix is not to unwrap paragraphs before
matching, which is the tempting analogue of the Python fix and is the wrong tool: markdown
renders a soft break *as whitespace*, so a citation split mid-token is broken in the
rendered document too, not merely hidden, and unwrapping cannot repair it. The fix is that
the separator between a label and its identifier is whitespace **including newlines**,
which is exactly what markdown says it is.

Two further markdown-only forms, both real in this repository:

* emphasis between the label and the digits -- ``PMID **12676948**``;
* a plural label introducing a list -- ``PMIDs 11102521, 11179418, 9712873`` -- where
  ``PMIDs`` does not even match ``PMID[: ]`` and the second and third identifiers have no
  label of their own.

Hence :data:`_RUN_SEP`: the run of characters allowed between a PMID label and an
identifier, and between successive identifiers in one labelled list. It is a closed set
of things that carry no meaning -- whitespace, markdown emphasis and code marks, list
punctuation -- rather than "any short gap", so it cannot wander into the next sentence.

What this scanner still cannot see is recorded in ``docs/CITATION_AUDIT.md``.
"""
from __future__ import annotations

import pathlib
import re
from dataclasses import dataclass

__all__ = ["Citation", "KINDS", "logical_text", "scan_text", "scan_file", "scan_paths"]

KINDS = ("pmid", "pmcid", "doi", "arxiv")


@dataclass(frozen=True, order=True)
class Citation:
    """One identifier, and the line of the file that claims it."""

    kind: str
    ident: str
    path: str
    line: int

    @property
    def where(self) -> str:
        return f"{self.path}:{self.line}"


# --------------------------------------------------------------------------- #
# logical text
# --------------------------------------------------------------------------- #

# Python implicit string concatenation splits a citation across lines -- `"... PMID "`
# followed by `"33654827 ..."` -- and a scanner reading the file rather than the value
# sees a quote where the digits should be, so the citation is invisible to it. Rejoining
# adjacent string literals before matching reconstructs the logical string, and is exact
# rather than a guess: the pattern being removed is precisely a closing quote, whitespace,
# and a reopening quote of the same kind.
_CONCAT = re.compile(r'"\s*"|\'\s*\'', re.S)


def logical_text(path: pathlib.Path) -> tuple[str, list[int]]:
    """The text to match against, plus the source line of each character.

    The line map is what lets a finding be reported as ``file:line`` after the text has
    been transformed. Without it a scanner that rejoins concatenated strings can say
    *that* a citation exists but not *where*, which is most of the value.
    """
    raw = path.read_text(encoding="utf-8", errors="replace")
    drop = bytearray(len(raw))
    if path.suffix == ".py":
        for match in _CONCAT.finditer(raw):
            for index in range(match.start(), match.end()):
                drop[index] = 1
    chars: list[str] = []
    lines: list[int] = []
    line = 1
    for index, char in enumerate(raw):
        if not drop[index]:
            chars.append(char)
            lines.append(line)
        if char == "\n":
            line += 1
    return "".join(chars), lines


# --------------------------------------------------------------------------- #
# identifier patterns
# --------------------------------------------------------------------------- #

_PMID_LABEL = re.compile(r"PMIDs?\b")

# Characters permitted between a PMID label and its digits, and between successive
# identifiers under one label. Whitespace covers the markdown soft line break, which is
# the whole reason this is not `[: ]`. The rest is markdown emphasis and code marks
# (`*`, `_`, backtick), bracketing, and list punctuation. Letters are excluded, so the
# run stops at the first word -- a label followed by prose yields nothing rather than
# reaching forward into an unrelated number.
_RUN_SEP = re.compile(r"(?:[\s:;,.()\[\]{}*_`=&#|/-]|\band\b)*")
_ID_DIGITS = re.compile(r"(?<!\d)(\d{6,8})(?!\d)")

_PMCID = re.compile(r"\bPMC(\d{6,8})\b")

# A DOI is self-identifying by shape. The suffix is deliberately permissive -- Elsevier
# suffixes contain parentheses (`10.1016/S0092-8674(00)81360-4`) and journals use dots,
# slashes and hyphens freely -- so the closing delimiter is decided afterwards by
# `_trim_doi` rather than by excluding characters DOIs legitimately contain.
_DOI = re.compile(r"\b10\.\d{4,9}/[^\s\"'<>\[\],;|`]+")

_ARXIV = re.compile(
    r"(?:arxiv\s*[:.]?\s*(?:abs/)?|arxiv\.org/(?:abs|pdf)/)"
    r"((?:\d{4}\.\d{4,5}|[a-z-]+(?:\.[A-Z]{2})?/\d{7})(?:v\d+)?)",
    re.I,
)

# A header cell naming PMID means "this column holds a PubMed identifier". A markdown
# table's header is an authoritative statement about its column, which makes it a better
# signal than any guess about which bare numbers in a table are identifiers -- and the
# audit documents in docs/ put wrong-and-right PMIDs in exactly such columns ("PMID",
# "PMID cited", "Correct PMID", "PMID / ID"), with no inline label anywhere on the row.
_PMID_COLUMN = re.compile(r"\bPMID", re.I)
# The rule under a table header: cells of dashes and optional alignment colons. A `|` row
# is a header if and only if this is the next row, which is what makes header detection
# structural rather than a guess about which row looks like a heading.
_TABLE_RULE = re.compile(r"^\|?(?:\s*:?-{1,}:?\s*\|)+\s*:?-*:?\s*\|?$")


def _trim_doi(text: str) -> str:
    """Drop the delimiters a DOI collected from the prose around it.

    Trailing sentence punctuation is never part of a DOI. A trailing ``)`` may be: it is
    kept while the parentheses balance, which is what distinguishes the DOI inside
    ``(doi:10.1016/S0092-8674(00)81360-4)`` from its closing bracket.
    """
    text = text.rstrip(".,;:*_")
    while text.endswith(")") and text.count("(") < text.count(")"):
        text = text[:-1].rstrip(".,;:*_")
    return text


def _labelled_pmids(text: str, lines: list[int]) -> list[tuple[str, int]]:
    """Every PMID reachable from a ``PMID``/``PMIDs`` label, with its line."""
    out: list[tuple[str, int]] = []
    for label in _PMID_LABEL.finditer(text):
        position = label.end()
        while True:
            gap = _RUN_SEP.match(text, position)
            digits = _ID_DIGITS.match(text, gap.end())
            if digits is None:
                break
            out.append((digits.group(1), lines[digits.start(1)]))
            position = digits.end()
    return out


def _table_pmids(raw: str) -> list[tuple[str, int]]:
    """PMIDs in the cells of any markdown table whose header names a PMID column.

    ``docs/SOURCE_AUDIT.md`` tabulates identifiers against what they resolve to, with the
    identifier alone in its own column and no inline label on the row. Those are the most
    load-bearing citations in the repository -- they are the record of which identifiers
    were checked -- and a label scan cannot see one of them.
    """
    out: list[tuple[str, int]] = []
    lines = raw.split("\n")
    columns: set[int] = set()
    for number, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped.startswith("|"):
            columns = set()  # the table ended
            continue
        if _TABLE_RULE.match(stripped):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        following = lines[number].strip() if number < len(lines) else ""
        if _TABLE_RULE.match(following):
            # A header, and the only row whose cells name their columns.
            columns = {i for i, cell in enumerate(cells) if _PMID_COLUMN.search(cell)}
            continue
        for index in sorted(columns):
            if index < len(cells):
                # A PMC-prefixed number is a PMCID by definition and never a PMID, so its
                # digits are masked before the column is read. Without this, a row that
                # tabulates a *PMCID* error under a "PMID" heading -- which is exactly what
                # an audit of a PMCID error looks like -- mints a PMID that does not exist.
                cell = _PMCID.sub(lambda m: "PMC" + "#" * len(m.group(1)), cells[index])
                for match in _ID_DIGITS.finditer(cell):
                    out.append((match.group(1), number))
    return out


def scan_text(text: str, lines: list[int], raw: str, path: str) -> list[Citation]:
    """Every identifier in one already-loaded file."""
    found: set[Citation] = set()

    def add(kind: str, ident: str, line: int) -> None:
        found.add(Citation(kind, ident, path, line))

    for pmid, line in _labelled_pmids(text, lines):
        add("pmid", pmid, line)
    if path.endswith(".md"):
        for pmid, line in _table_pmids(raw):
            add("pmid", pmid, line)
    for match in _PMCID.finditer(text):
        add("pmcid", "PMC" + match.group(1), lines[match.start()])
    for match in _DOI.finditer(text):
        add("doi", _trim_doi(match.group(0)).lower(), lines[match.start()])
    for match in _ARXIV.finditer(text):
        add("arxiv", match.group(1).lower(), lines[match.start()])
    return sorted(found)


def scan_file(path: pathlib.Path, root: pathlib.Path) -> list[Citation]:
    text, lines = logical_text(path)
    raw = path.read_text(encoding="utf-8", errors="replace")
    return scan_text(text, lines, raw, str(path.relative_to(root)))


def scan_paths(root: pathlib.Path, bases: tuple[str, ...],
               suffixes: tuple[str, ...] = (".py", ".md")) -> list[Citation]:
    """Every identifier cited anywhere under ``bases``, sorted and de-duplicated."""
    out: list[Citation] = []
    for base in bases:
        directory = root / base
        if not directory.exists():
            continue
        for path in sorted(directory.rglob("*")):
            if path.suffix not in suffixes or not path.is_file():
                continue
            if "__pycache__" in path.parts or ".pytest_cache" in path.parts:
                continue
            out.extend(scan_file(path, root))
    return sorted(set(out))


def by_identifier(citations: list[Citation]) -> dict[str, dict[str, set[str]]]:
    """``{kind: {identifier: {"file:line", ...}}}``."""
    out: dict[str, dict[str, set[str]]] = {kind: {} for kind in KINDS}
    for citation in citations:
        out[citation.kind].setdefault(citation.ident, set()).add(citation.where)
    return out
