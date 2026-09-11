"""Every identifier in the source *and in the documentation* must resolve to a paper
somebody has looked at.

The panel's citation culture had a hole in exactly the shape of its own test. Reporters
and modules carry a ``source`` string naming a PMID, and the test that guarded it
asserted ``bool(module.source)`` -- that a citation *exists*. A PMID pointing at a
butterfly wing-pattern paper passes that assertion, and several did: the DNA-damage
regulon cited a butterfly, the calcium regulon cited a retinoblastoma study, the
retrograde regulon cited cauliflower mosaic virus, and the ATP sensor cited a
dermoscopy paper while asserting on its authority that the sensor had been shown in
yeast. Named authors and years were right throughout; only the digits drifted, which is
why reading the prose never caught it.

Checking existence catches nothing. Checking the digits against the internet would put
a network call in the test suite, which fails offline and makes the build depend on
NCBI being up. So the identifiers are pinned to vendored tables instead -- one per
resolver, under ``data/citations/`` -- refreshed deliberately by
``scripts/refresh_citations.py``, and the tests assert that every identifier cited
anywhere appears in the table for its kind. Adding a citation then requires refreshing,
which puts the resolved title in the diff where a reviewer sees it beside the claim it is
attached to.

``docs/`` is in scope, and is the larger surface: about six thousand lines of research
notes carrying the literature the project's conclusions rest on -- the identifiability
bound, the simulation-optimism theorem, the qPCR correctability ceiling, the UPR anchor
recommendation. Those documents also cite DOIs and arXiv preprints, which code does not,
and which a PubMed-only guard could not see at all.

That is the most this can do automatically. Whether a paper *supports* the sentence next
to it is a judgement, and no test makes it; ``docs/CITATION_AUDIT.md`` records the ones
that were checked by hand and what was wrong. What these tests remove is the failure mode
where nobody could have noticed.
"""

from __future__ import annotations

import pathlib

import pandas as pd
import pytest

from ystwin import citations as scan

REPO = pathlib.Path(__file__).resolve().parents[1]
TABLES = REPO / "data" / "citations"
BASES = ("src", "scripts", "tests", "docs")

# The column holding the identifier in each table, and the label used in messages.
KEY = {"pmid": "pmid", "pmcid": "pmcid", "doi": "doi", "arxiv": "arxiv_id"}


def _table(kind: str) -> pd.DataFrame:
    path = TABLES / f"{kind}_titles.csv"
    if not path.exists():
        pytest.skip(f"citation table not vendored at {path}")
    return pd.read_csv(path, dtype=str).fillna("")


@pytest.fixture(scope="module")
def cited() -> dict[str, dict[str, set[str]]]:
    return scan.by_identifier(scan.scan_paths(REPO, BASES))


@pytest.fixture(scope="module")
def tables() -> dict[str, pd.DataFrame]:
    return {kind: _table(kind) for kind in scan.KINDS}


# --------------------------------------------------------------------------- #
# the tables themselves
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", scan.KINDS)
def test_the_table_exists_and_is_populated(tables, kind):
    frame = tables[kind]
    assert len(frame) > 0
    assert KEY[kind] in frame.columns
    assert "title" in frame.columns
    assert "cited_in" in frame.columns


@pytest.mark.parametrize("kind", scan.KINDS)
def test_every_cited_identifier_is_in_the_table(tables, cited, kind):
    """The check the old one should have been.

    A new identifier that nobody has resolved fails here, which is the point: the fix is
    to refresh the table, and the refresh puts the title in the diff.
    """
    known = set(tables[kind][KEY[kind]])
    missing = {i: sorted(w) for i, w in cited[kind].items() if i not in known}
    assert not missing, (
        f"{kind} identifiers cited but absent from {kind}_titles.csv: {missing}. "
        "Run scripts/refresh_citations.py and check each resolved title against the "
        "claim it sits beside."
    )


@pytest.mark.parametrize("kind", scan.KINDS)
def test_no_identifier_resolves_to_an_empty_title(tables, kind):
    """An unresolvable identifier is worse than none: it looks like provenance.

    The one admissible exception is a DOI that exists but has no bibliographic record at
    any registry -- a supplementary-file DOI, for instance. That is recorded as
    ``registry == "handle-only"`` rather than left blank, so the distinction between
    "cannot be described" and "does not exist" stays visible instead of being smoothed
    over by whichever answer the test wanted.
    """
    frame = tables[kind]
    blank = frame[frame.title.str.strip() == ""]
    if kind == "doi":
        blank = blank[blank.registry != "handle-only"]
    assert blank.empty, f"unresolvable {kind}s: {sorted(blank[KEY[kind]])}"


def test_every_doi_records_where_it_resolved(tables):
    """A DOI with no registry is a DOI nobody can re-check.

    Only the DOI table needs this: it is the only one with more than one possible
    registry, because a DOI's registration agency is not knowable from the string.
    """
    blank = tables["doi"][tables["doi"].registry.str.strip() == ""]
    assert blank.empty, (
        f"DOIs that resolved at no registry and have no handle: {sorted(blank.doi)}. "
        "An identifier that resolves nowhere is not provenance."
    )


@pytest.mark.parametrize("kind", scan.KINDS)
def test_the_table_has_no_duplicate_entries(tables, kind):
    frame = tables[kind]
    duplicated = frame[frame[KEY[kind]].duplicated()][KEY[kind]].tolist()
    assert not duplicated, f"duplicate rows for {duplicated}"


@pytest.mark.parametrize("kind", scan.KINDS)
def test_the_table_carries_no_stale_entries(tables, cited, kind):
    """An identifier in the table but nowhere in the repository means a citation was
    removed.

    Harmless, and worth surfacing: the tables are meant to be a record of what this
    repository cites, and an entry with no citer is a leftover that quietly grows.
    """
    stale = sorted(set(tables[kind][KEY[kind]]) - set(cited[kind]))
    assert not stale, (
        f"{kind}_titles.csv carries {len(stale)} entries no longer cited: {stale}. "
        "Re-run scripts/refresh_citations.py to prune them."
    )


@pytest.mark.parametrize("kind", scan.KINDS)
def test_documentation_is_actually_covered(cited, kind):
    """The guard used to stop at the ``tests/`` boundary, and prose was the larger surface.

    This asserts the scan reaches ``docs/`` at all. Without it, a refactor that narrowed
    the walked directories back to code would leave every table passing vacuously --
    which is how the previous hole survived: nothing failed when nothing was looked at.
    """
    from_docs = [i for i, where in cited[kind].items()
                 if any(w.startswith("docs/") for w in where)]
    assert from_docs, f"no {kind} identifiers found in docs/ -- is the scan still walking it?"


# --------------------------------------------------------------------------- #
# regression pins on the scanner
# --------------------------------------------------------------------------- #


def test_a_citation_wrapped_across_a_markdown_line_is_still_found():
    """The markdown analogue of the string-concatenation hole, pinned.

    ``docs/research/UPR_ANCHOR.md`` wraps two reference entries so that the word ``PMID``
    ends one line and the digits begin the next. The original pattern was
    ``PMID[: ]+(\\d{6,8})``, and ``[: ]`` does not match a newline, so both citations were
    invisible to the guard while being perfectly legible to a reader. Narrowing the
    separator back to a literal space and colon fails here.
    """
    found = scan.by_identifier(scan.scan_paths(REPO, ("docs",)))["pmid"]
    for pmid in ("34047633", "12184808"):
        where = found.get(pmid, set())
        assert any("UPR_ANCHOR" in w for w in where), (
            f"PMID {pmid} wraps across a line break in UPR_ANCHOR.md and the scanner no "
            f"longer sees it; found only at {sorted(where)}"
        )


def test_a_plural_label_introducing_a_list_yields_every_identifier():
    """``PMIDs 11102521, 11179418, 9712873`` is three citations, not zero.

    ``PMIDs`` does not match ``PMID[: ]``, and the second and third identifiers carry no
    label of their own. All three describe the same claim and all three must be checkable.
    """
    text = "conventional dose | PMIDs 11102521, 11179418, 9712873 | High"
    lines = [1] * len(text)
    found = {c.ident for c in scan.scan_text(text, lines, text, "x.md")}
    assert found == {"11102521", "11179418", "9712873"}


def test_emphasis_between_the_label_and_the_digits_does_not_hide_a_citation():
    text = "Young 2003, PMID **12676948** (note: the file still cites something else)"
    lines = [1] * len(text)
    found = {c.ident for c in scan.scan_text(text, lines, text, "x.md")}
    assert "12676948" in found


def test_a_pmid_column_in_a_table_is_read_even_with_no_inline_label():
    """``docs/SOURCE_AUDIT.md`` tabulates identifiers against what they resolve to.

    The identifier sits alone in its own column with no ``PMID`` anywhere on the row, and
    those rows are the repository's record of which identifiers were checked. A label
    scan cannot see one of them; the table header can.
    """
    text = ("| PMID | Resolves to |\n"
            "| --- | --- |\n"
            "| 25955212 | Zhao 2015, SoNar sensor |\n")
    lines = [1] * len(text)
    found = {c.ident for c in scan.scan_text(text, lines, text, "x.md")}
    assert "25955212" in found


def test_a_table_row_is_not_mistaken_for_a_header():
    """Only the row above the ``---`` rule names its columns.

    Treating any ``|`` row mentioning PMID as a header would make page ranges and sample
    sizes in an unrelated column look like citations, and every false identifier here
    costs somebody a lookup.
    """
    text = ("| Location | Note |\n"
            "| --- | --- |\n"
            "| stress_panel.py | the PMID here is fine |\n"
            "| pages 1234567 | not an identifier |\n")
    lines = [1] * len(text)
    found = {c.ident for c in scan.scan_text(text, lines, text, "x.md")}
    assert "1234567" not in found


def test_a_doi_containing_parentheses_survives_its_delimiters():
    """Elsewhere-legal DOIs contain brackets, and markdown wraps them in more brackets."""
    text = "doi:[10.1016/S0092-8674(00)81360-4](https://doi.org/10.1016/S0092-8674(00)81360-4)."
    lines = [1] * len(text)
    found = {c.ident for c in scan.scan_text(text, lines, text, "x.md") if c.kind == "doi"}
    assert found == {"10.1016/s0092-8674(00)81360-4"}


def test_python_string_concatenation_still_does_not_hide_a_citation(tmp_path):
    """The original hole, kept pinned now that the helper lives in one place."""
    path = tmp_path / "m.py"
    path.write_text('SOURCE = ("Takaine 2019, PMID "\n          "30858198 -- yeast ATP")\n')
    found = {c.ident for c in scan.scan_file(path, tmp_path)}
    assert "30858198" in found


# --------------------------------------------------------------------------- #
# regression pins on specific errors that have happened here
# --------------------------------------------------------------------------- #


def test_validprime_is_cited_by_the_identifier_that_is_actually_validprime(tables):
    """``PMC3326325`` and ``PMC3326333`` are adjacent identifiers in one journal issue.

    ``docs/research/G4_STATISTICS.md`` cited the first for the ValidPrime paper -- the
    primary reference for the 60% correctability ceiling that the G4 gate threshold is
    derived from -- and it is an unrelated ribosomal-protein paper. Both identifiers are
    still cited, because the document now explains the error; this pins which one is
    which so the explanation cannot quietly invert.
    """
    frame = tables["pmcid"].set_index("pmcid")
    assert "ValidPrime" in frame.loc["PMC3326333", "title"]
    assert "ValidPrime" not in frame.loc["PMC3326325", "title"]
    assert "ribosomal protein S26" in frame.loc["PMC3326325", "title"]


def test_the_queen_sensor_no_longer_claims_a_yeast_demonstration_it_lacks():
    """A regression pin on the specific failure that motivated all of this.

    QUEEN-2m's source asserted the sensor was 'shown in yeast' on the authority of
    Yaginuma 2014 -- which is an E. coli paper, cited by an identifier belonging to a
    dermoscopy study. The sensor genuinely does work in yeast; the demonstration is
    Takaine, and the file did not carry it. This fails if either error returns.
    """
    from ystwin.generator.stress_panel import REPORTERS

    source = REPORTERS["QUEEN-2m"].source
    assert "24815987" not in source, "the dermoscopy PMID is back"
    assert "25283467" in source, "Yaginuma's real PMID should be cited"
    assert "30858198" in source or "Takaine" in source, (
        "the yeast demonstration must be cited, not implied"
    )


def test_every_reporter_and_module_still_carries_a_source():
    """The original assertion, kept.

    Existence is necessary and was never sufficient; it is retained so that removing a
    citation altogether still fails, rather than being quietly permitted by the newer
    and stricter checks passing vacuously on an empty string.
    """
    from ystwin.generator.stress_panel import MODULES, REPORTERS

    for name, item in list(MODULES.items()) + list(REPORTERS.items()):
        assert getattr(item, "source", "").strip(), f"{name} has no source"
