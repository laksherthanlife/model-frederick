"""Re-resolve every identifier this repository cites, and vendor what it resolves to.

The point is not the files; it is the diff. An identifier drifts by a digit and still
looks like provenance, and prose review never catches it because the author and year
beside it are right. Resolving the identifiers into tracked tables means the next person
to add or change a citation produces a diff containing the paper's actual title, sitting
next to the claim it is attached to, where a reviewer can see whether the two match.

An audit of one source file found 14 of its 31 PMIDs pointing at unrelated papers -- a
butterfly wing-pattern study cited for a DNA-damage regulon, a dermoscopy paper cited for
an ATP sensor -- and that is what this exists to stop happening again unnoticed.

What this cannot do is decide whether the paper supports the sentence. That stays a
judgement. What it removes is the state where nobody could have noticed.

Four tables, one per resolver, because they are four different questions:

``pmid_titles.csv``
    PubMed, via NCBI E-utilities. PMIDs cited from code *and* from prose live in one
    table: same resolver, same schema, and a PMID cited from both surfaces is one paper
    that should appear once, with both citers listed.
``pmcid_titles.csv``
    PubMed Central, via E-utilities against ``db=pmc``. Separate because it carries the
    PMC-to-PMID mapping, and that mapping is the check: the one PMCID error this
    repository has had (``PMC3326325`` for ValidPrime, which is a ribosomal-protein
    paper) is caught by looking at what PMID the PMCID belongs to.
``doi_titles.csv``
    Crossref, then DataCite, then the DOI handle system. Most of the literature the
    documentation rests on -- statistics journals, proceedings, book chapters -- carries a
    DOI and no PMID at all, so a PubMed-only guard covers none of it. Three tiers because
    a DOI's registration agency is not knowable from the string: the dataset DOIs this
    repository cites (Zenodo, Edinburgh DataShare) are registered with DataCite and are
    simply absent from Crossref, and calling them unresolvable would be wrong. The
    ``registry`` column records which one answered, so a DOI that exists but carries no
    bibliographic record anywhere -- a supplementary-file DOI, for instance -- is written
    down as exactly that rather than as a resolution failure or a fake title.
``arxiv_titles.csv``
    The arXiv API, plus Crossref for the published version where the preprint declares
    one. Two extra columns, and they are the entire reason this is not folded into the
    DOI table: ``version``, because a result can change between preprint versions so an
    unversioned identifier does not pin the formula being cited; and ``published_title``,
    because an arXiv title that differs from the published title lets a citation name one
    document and point at another while reading as perfectly consistent.

Run it after adding or changing any identifier; ``tests/test_citations.py`` fails until
you do.

Usage: python scripts/refresh_citations.py [--check] [--kind pmid|pmcid|doi|arxiv]
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import citations as scan

REPO = pathlib.Path(__file__).resolve().parents[1]
TABLES = REPO / "data" / "citations"
BASES = ("src", "scripts", "tests", "docs")

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CROSSREF = "https://api.crossref.org/works/"
DATACITE = "https://api.datacite.org/dois/"
HANDLE = "https://doi.org/api/handles/"
ARXIV = "http://export.arxiv.org/api/query"

# Crossref asks callers to identify themselves. It gets a project name and nothing else:
# a contact address here would put a personal identifier into every outbound request.
AGENT = "ystwin-citation-guard/1.0 (https://github.com/; academic reproducibility check)"

SCHEMA = {
    "pmid": ["pmid", "title", "journal", "year", "cited_in"],
    "pmcid": ["pmcid", "pmid", "title", "journal", "year", "cited_in"],
    "doi": ["doi", "registry", "title", "journal", "year", "cited_in"],
    "arxiv": ["arxiv_id", "version", "title", "published_title", "published_doi",
              "cited_in"],
}


def table_path(kind: str) -> pathlib.Path:
    return TABLES / f"{kind}_titles.csv"


def _get(url: str, timeout: int = 45) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


# --------------------------------------------------------------------------- #
# resolvers
# --------------------------------------------------------------------------- #


def _esummary(db: str, ids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for start in range(0, len(ids), 50):
        chunk = ids[start:start + 50]
        url = f"{ESUMMARY}?db={db}&id={','.join(chunk)}&retmode=json"
        payload = json.loads(_get(url))["result"]
        out.update({i: payload.get(i, {}) for i in chunk})
        time.sleep(0.4)  # NCBI asks for no more than three requests a second
    return out


def resolve_pmid(ids: list[str]) -> list[dict]:
    records = _esummary("pubmed", ids)
    return [{
        "pmid": i,
        "title": (records[i].get("title", "") or "").rstrip("."),
        "journal": records[i].get("source", ""),
        "year": (records[i].get("pubdate", "") or "")[:4],
    } for i in ids]


def resolve_pmcid(ids: list[str]) -> list[dict]:
    """PMCIDs, and the PMID each one belongs to.

    The mapping is the point. A PMCID that resolves to a title is only half the check;
    the failure this repository met was a PMCID one digit off, and the fastest way to see
    that is that its PMID is not the PMID the prose names beside it.
    """
    records = _esummary("pmc", [i.removeprefix("PMC") for i in ids])
    rows = []
    for cited in ids:
        i = cited.removeprefix("PMC")
        record = records[i]
        pmid = ""
        for article in record.get("articleids", []):
            if article.get("idtype") == "pmid":
                pmid = str(article.get("value", ""))
        rows.append({
            "pmcid": cited,
            "pmid": pmid,
            "title": (record.get("title", "") or "").rstrip("."),
            "journal": record.get("source", ""),
            "year": (record.get("pubdate", "") or "")[:4],
        })
    return rows


def _crossref(doi: str) -> dict | None:
    try:
        payload = json.loads(_get(CROSSREF + urllib.parse.quote(doi, safe="")))
    except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError):
        return None
    return payload.get("message")


def _crossref_fields(message: dict | None) -> dict:
    if not message:
        return {"title": "", "journal": "", "year": ""}
    titles = message.get("title") or []
    containers = message.get("container-title") or []
    issued = (message.get("issued") or {}).get("date-parts") or [[]]
    year = str(issued[0][0]) if issued and issued[0] else ""
    return {
        "title": (titles[0] if titles else "").strip().rstrip("."),
        "journal": (containers[0] if containers else "").strip(),
        "year": year,
    }


def _datacite_fields(doi: str) -> dict | None:
    """Bibliographic fields for a DOI registered with DataCite rather than Crossref."""
    try:
        payload = json.loads(_get(DATACITE + urllib.parse.quote(doi, safe="/")))
    except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError):
        return None
    attributes = (payload.get("data") or {}).get("attributes") or {}
    titles = attributes.get("titles") or []
    return {
        "title": (titles[0].get("title", "") if titles else "").strip().rstrip("."),
        "journal": attributes.get("publisher") or "",
        "year": str(attributes.get("publicationYear") or ""),
    }


def _handle_exists(doi: str) -> bool:
    """Whether the DOI resolves at all, when no registry will describe it.

    A supplementary-file DOI has a handle and no metadata record. Distinguishing that
    from a fabricated identifier matters: one is a citation this guard cannot describe,
    the other is a citation that does not exist.
    """
    try:
        payload = json.loads(_get(HANDLE + urllib.parse.quote(doi, safe="/")))
    except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError):
        return False
    return payload.get("responseCode") == 1


def resolve_doi(ids: list[str]) -> list[dict]:
    rows = []
    for index, doi in enumerate(ids):
        fields = _crossref_fields(_crossref(doi))
        registry = "crossref" if fields["title"] else ""
        if not registry:
            from_datacite = _datacite_fields(doi)
            if from_datacite and from_datacite["title"]:
                fields, registry = from_datacite, "datacite"
        if not registry and _handle_exists(doi):
            registry = "handle-only"
        rows.append({"doi": doi, "registry": registry, **fields})
        if index % 10 == 9:
            print(f"    ... {index + 1}/{len(ids)} DOIs", file=sys.stderr)
        time.sleep(0.15)
    return rows


_ATOM = {"a": "http://www.w3.org/2005/Atom", "arxiv": "http://arxiv.org/schemas/atom"}


def resolve_arxiv(ids: list[str]) -> list[dict]:
    """Preprint titles, the version cited, and the published title where one exists.

    Entries are matched to the identifier that asked for them by the ``<id>`` each entry
    carries, never by position. The arXiv API does not promise to answer an ``id_list`` in
    the order it was given, and it does not: keying on position silently pairs every title
    with the wrong identifier, which is precisely the failure this whole guard exists to
    catch, committed by the guard itself.
    """
    resolved: dict[str, dict] = {}
    for start in range(0, len(ids), 25):
        chunk = ids[start:start + 25]
        query = urllib.parse.urlencode({"id_list": ",".join(chunk),
                                        "max_results": len(chunk)})
        tree = ET.fromstring(_get(f"{ARXIV}?{query}"))
        for entry in tree.findall("a:entry", _ATOM):
            served = (entry.findtext("a:id", "", _ATOM) or "").rsplit("/", 1)[-1]
            doi = (entry.findtext("arxiv:doi", "", _ATOM) or "").strip()
            resolved[served.rsplit("v", 1)[0]] = {
                "served": served,
                "title": " ".join((entry.findtext("a:title", "", _ATOM) or "").split()),
                "doi": doi,
            }
        time.sleep(3.0)  # arXiv asks for one request every three seconds

    rows: list[dict] = []
    for cited in ids:
        record = resolved.get(cited.rsplit("v", 1)[0], {})
        doi = record.get("doi", "")
        published = _crossref_fields(_crossref(doi))["title"] if doi else ""
        if doi:
            time.sleep(0.15)
        served = record.get("served", "")
        rows.append({
            "arxiv_id": cited,
            # What the citation pins versus what arXiv currently serves. A citation with
            # no version pins nothing, and that is written down rather than implied: a
            # result that moved between versions is cited by a string that still looks
            # right.
            "version": served if "v" in cited else f"unpinned:{served}",
            "title": record.get("title", ""),
            "published_title": published,
            "published_doi": doi,
        })
    return rows


RESOLVERS = {"pmid": resolve_pmid, "pmcid": resolve_pmcid,
             "doi": resolve_doi, "arxiv": resolve_arxiv}
KEY = {"pmid": "pmid", "pmcid": "pmcid", "doi": "doi", "arxiv": "arxiv_id"}


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #


def _cited() -> dict[str, dict[str, set[str]]]:
    return scan.by_identifier(scan.scan_paths(REPO, BASES))


def _rows_for(kind: str, cited: dict[str, set[str]]) -> list[dict]:
    ids = sorted(cited)
    rows = RESOLVERS[kind](ids)
    for row in rows:
        row["cited_in"] = ";".join(sorted(cited[row[KEY[kind]]]))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report drift and exit non-zero; do not write")
    parser.add_argument("--kind", choices=scan.KINDS, action="append",
                        help="refresh only these tables (default: all)")
    parser.add_argument("--missing-only", action="store_true",
                        help="resolve only new identifiers; preserve existing bibliographic records")
    args = parser.parse_args()

    cited = _cited()
    kinds = args.kind or list(scan.KINDS)
    status = 0

    for kind in kinds:
        if not cited[kind]:
            print(f"{kind}: none cited")
            continue
        path = table_path(kind)
        retained = {}
        if args.missing_only and path.exists():
            with path.open(newline="", encoding="utf-8") as handle:
                retained = {row[KEY[kind]]: row for row in csv.DictReader(handle)}
        requested = {identifier: locations for identifier, locations in cited[kind].items()
                     if identifier not in retained}
        print(f"{kind}: resolving {len(requested)} identifiers")
        try:
            resolved = _rows_for(kind, requested) if requested else []
        except OSError as exc:
            print(f"could not resolve {kind}: {exc}", file=sys.stderr)
            return 2
        for row in resolved:
            if row.get("year") in (None, "None"):
                row["year"] = ""
        combined = {**retained, **{row[KEY[kind]]: row for row in resolved}}
        rows = [{key: value.replace("\r\n", "\n") if isinstance(value, str) else value
                 for key, value in combined[identifier].items()} for identifier in sorted(combined)]

        # A DOI that exists but has no bibliographic record anywhere is recorded, not
        # counted as a failure; a DOI that resolves nowhere is the failure.
        unresolved = [r[KEY[kind]] for r in rows
                      if not r["title"] and r.get("registry") != "handle-only"]
        if unresolved:
            print(f"  UNRESOLVED, and an unresolvable identifier looks like provenance "
                  f"while carrying none: {unresolved}", file=sys.stderr)
            status = 1

        path = table_path(kind)
        if args.check:
            if not path.exists():
                print(f"  no vendored table at {path}; run without --check",
                      file=sys.stderr)
                status = 1
                continue
            existing = {r[KEY[kind]]: r for r in csv.DictReader(path.open())}
            drift = [r for r in rows
                     if r[KEY[kind]] not in existing
                     or existing[r[KEY[kind]]]["title"] != r["title"]]
            for row in drift:
                was = existing.get(row[KEY[kind]], {}).get("title", "<absent>")
                print(f"  {row[KEY[kind]]}: {was!r} -> {row['title']!r}")
            print(f"  {len(drift)} drifted or new of {len(rows)}")
            if drift:
                status = 1
            continue

        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SCHEMA[kind], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"  wrote {path.relative_to(REPO)} with {len(rows)} entries")

    if not args.check:
        print("Read the diff: each title should support the claim it sits beside.")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
