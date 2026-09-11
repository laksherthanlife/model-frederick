"""Check repository references and the declared current scientific claim inventory.

References, dependency pins and marked test counts are checked against the repository.
Scientific records and explicit value markers are evaluated once by
``ystwin.analysis.claims.audit_claims`` using ``data/current_claims.json``. Unmarked
numbers in prose are not a scientific inventory. Matching a saved cell does not establish
support without the declared method, configuration, coverage and provenance checks.

Historical records are nonblocking SKIPs; a verified scientific refusal is a PASS, not a
positive scientific result. Unknown provenance blocks current support.

Verification is read-only by default. An explicit ``--output-dir`` writes an audit report.
A blocking record always produces a nonzero exit.

Usage: python scripts/audit_claims.py [--quiet] [--output-dir DIRECTORY]
"""
from __future__ import annotations

import argparse
import importlib
import pathlib
import re
import subprocess
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

pd = importlib.import_module("pandas")
audit_claims = importlib.import_module("ystwin.analysis.claims").audit_claims

REPO = pathlib.Path(__file__).resolve().parents[1]
README = REPO / "README.md"
SRC = REPO / "src" / "ystwin"


def _row(check: str, expected: object, actual: object, passed: bool, detail: str = ""):
    return {
        "check": check,
        "record_type": "repository",
        "expected": str(expected),
        "actual": str(actual),
        "status": "PASS" if passed else "FAIL",
        "blocking": not passed,
        "detail": detail,
    }


# ---------------------------------------------------------------------------
# individual checks
# ---------------------------------------------------------------------------


_SEARCH_ROOTS = ("src/ystwin", "scripts", "tests")
_OWNED_PREFIXES = ("src/", "scripts/", "tests/")


def _module_exists(name: str) -> bool:
    """Whether a documented ``module.py`` reference resolves anywhere it could live.

    Documentation names modules three ways -- package-relative (``gates/g1_optical.py``),
    repo-relative (``scripts/run_g4.py``) and bare (``reporter.py``) -- and the
    architecture map legitimately names test modules too. Searching only the package
    would report those as missing, which is a false alarm in the check whose entire
    value is that its alarms are trustworthy.
    """
    direct = [SRC / name, REPO / name] + [REPO / root / name for root in _SEARCH_ROOTS]
    if any(path.exists() for path in direct):
        return True
    if "/" in name:
        return False
    # A bare basename may sit at any depth under one of the roots.
    return any(
        any((REPO / root).rglob(name))
        for root in _SEARCH_ROOTS
        if (REPO / root).exists()
    )


def _collected_test_count() -> int | None:
    """How many tests pytest collects here, or ``None`` if collection could not run."""
    result = subprocess.run(
        [sys.executable, "-B", "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider"],
        cwd=REPO, capture_output=True, text=True,
    )
    found = re.search(r"(\d+)\s+tests? collected", result.stdout)
    return int(found.group(1)) if result.returncode == 0 and found else None


def check_test_count(text: str) -> list[dict]:
    """Any marked test count in the README must match what pytest collects.

    Collection rather than a full run: the question is how many tests the document
    thinks exist, and collecting is seconds where running is minutes.
    """
    # Read a marked claim, not free prose. Scanning for "N tests" anywhere also matched
    # narrative mentions -- "the README claiming 149 tests when there were 1,204" -- and
    # flagged the history of this very check as a live false claim. A marker is the
    # authoritative signal; a regex over English is a guess about intent.
    markers = len(re.findall(r"<!--\s*audit:test_count\s*-->", text))
    claimed = [int(m.replace(",", ""))
               for m in re.findall(r"<!--\s*audit:test_count\s*-->[^\n]*?(\d[\d,]*)", text)]
    if markers != len(claimed):
        # A MARKER WITH NO NUMBER AFTER IT IS NOT THE SAME AS NO MARKER, and this returned
        # PASS for it. A one-character slip in a sed left `The suite is ** tests**.` and the
        # gate reported "none found" and went green -- the claim had not been withdrawn, it
        # had been emptied, which is the state a reader is least able to notice.
        return [_row("test count claimed in README", "a number after each marker",
                     "marker present, no number", False,
                     f"{markers - len(claimed)} <!-- audit:test_count --> marker(s) with no digits after "
                     "them on the same line. Restore the count or remove the marker")]
    if not claimed:
        return [_row("test count claimed in README", "a marked number", "none found", True,
                     "README carries no <!-- audit:test_count --> marker; nothing to drift")]

    result = subprocess.run(
        [sys.executable, "-B", "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider"],
        cwd=REPO, capture_output=True, text=True,
    )
    # The summary line, not a per-file "path: N" tally. The per-file format needs -qq;
    # under -q pytest prints node ids and the tally never matched, so this check could
    # only ever fail -- and it read as passing because the README carried no claim.
    summary = re.search(r"(\d+)\s+tests? collected", result.stdout)
    actual = int(summary.group(1)) if summary else 0
    if result.returncode != 0 or actual == 0:
        return [_row("test count", claimed, "could not collect", False,
                     f"pytest exited {result.returncode}: "
                     f"{(result.stdout + result.stderr)[-500:]!r}")]

    rows = []
    for value in claimed:
        # EXACT, and it was not. The tolerance here read
        # `abs(value - actual) <= max(5, 0.02 * actual)` with the comment "allow a little
        # slack: tests get added between a README edit and a commit" -- which is precisely
        # the drift the check exists to catch, and 2% of a suite this size is 55 tests. It
        # let the README sit at 2721 against a collected 2762 and report PASS with both
        # numbers printed side by side in `outputs/audit_claims.csv`, which is worse than
        # no check: a gate that prints a mismatch and calls it a pass teaches a reader to
        # stop reading the column. The README says of this check that it "fails if this
        # line drifts", so either the tolerance goes or that sentence does.
        exact = value == actual
        rows.append(_row("test count claimed in README", value, actual, exact,
                         "" if exact else
                         f"README says {value}, pytest collects {actual} -- update the "
                         "number beside the <!-- audit:test_count --> marker"))
    return rows


def check_module_references(text: str) -> list[dict]:
    """Every ``some/module.py`` the README names in backticks must exist.

    This is the check that catches a deleted module still documented in the contents
    table -- the failure that makes a reader distrust the whole document.
    """
    rows = []
    referenced = sorted(set(re.findall(r"`([a-z_][a-z_0-9/]*\.py)`", text)))
    for name in referenced:
        found = _module_exists(name)
        rows.append(_row(f"module referenced: {name}", "exists",
                         "exists" if found else "MISSING", found,
                         "" if found else "named in README, absent from the tree"))
    return rows


def check_output_tables(text: str) -> list[dict]:
    """Every ``outputs/...`` file the README names must be present and non-empty.

    A result nobody can open is not a result. This is the check that makes the
    difference between citing a table and having one.
    """
    rows = []
    referenced = sorted({_trim_link_paren(n)
                         for n in re.findall(r"`?(outputs/[\w./\-()]+)`?", text)})
    if not referenced:
        return [_row("outputs referenced in README", "any", "none", True,
                     "README cites no output tables; consider citing them")]
    for name in referenced:
        path = REPO / name
        if not path.exists():
            rows.append(_row(f"output table: {name}", "exists", "MISSING", False))
        elif path.stat().st_size == 0:
            rows.append(_row(f"output table: {name}", "non-empty", "empty", False))
        else:
            rows.append(_row(f"output table: {name}", "exists", "exists", True))
    return rows


def _trim_link_paren(name: str) -> str:
    """Drop a closing paren belonging to a Markdown link rather than the filename.

    Balanced parentheses in plate-reader export filenames are retained.
    """
    # SENTENCE PUNCTUATION FIRST, because a link at the end of a sentence hands back
    # ``outputs/x.csv).`` and the balance rule below only ever looks at the LAST character.
    # A trailing full stop therefore hid the paren from it, the name never matched a file,
    # and the check reported a committed table as MISSING -- an alarm that is wrong for a
    # reason having nothing to do with the table, which is the one gate whose whole job is
    # to be believed. No output filename here ends in a full stop or a comma.
    name = name.rstrip(".,;:")
    while name.endswith(")") and name.count(")") > name.count("("):
        name = name[:-1]
        name = name.rstrip(".,;:")
    return name


def check_documented_commands(text: str) -> list[dict]:
    """Every ``python scripts/x.py`` in the README must name a script that exists."""
    rows = []
    referenced = sorted(set(re.findall(r"python3?\s+(scripts/[\w/]+\.py)", text)))
    for name in referenced:
        path = REPO / name
        rows.append(_row(f"documented command: {name}", "exists",
                         "exists" if path.exists() else "MISSING", path.exists()))
    return rows


# Knowledge of where a particular user keeps their files belongs in exactly one place.
# ``paths.py`` is that place: it resolves each asset by environment variable first and
# falls back to a home-relative default only as a documented last resort, so a machine
# that stores its data elsewhere overrides rather than edits. This audit names the
# pattern it looks for, so it necessarily matches itself.
#: Files allowed to contain a LITERAL absolute home path. Each needs a reason.
_LITERAL_PATH_EXEMPT = {
    "src/ystwin/paths.py",       # the resolver's own documented fallbacks
    "scripts/audit_claims.py",   # this file names the pattern it searches for
    # Quotes ANOTHER repository's hardcoded path, as evidence about that codebase rather
    # than as knowledge this one uses. The invariant is "only the resolver knows where a
    # file lives on a machine"; a comparison document reporting that a different project
    # got this wrong does not violate it.
    "docs/COMPARISON.md",
}

#: Python files allowed to CALL ``Path.home()``. Prose that merely names the call is not
#: checked -- see :func:`check_no_absolute_home_paths` for why the two are separate.
_HOME_CALL_EXEMPT = {
    "src/ystwin/paths.py",       # the resolver, which is the whole point
    "scripts/audit_claims.py",   # names the pattern it searches for
    # Calls it to REDACT: `_portable` turns a resolved location into `~/...` so the tracked
    # CSV does not carry a username. That is the opposite of locating data by convention,
    # and exempting it is narrower than weakening the pattern.
    "scripts/audit_reproducibility.py",
}


def check_no_absolute_home_paths() -> list[dict]:
    """Literal home paths cannot be committed; only the resolver may call Path.home().

    Every tracked file is inspected for machine-specific data. Only Python sources are
    inspected for resolver calls: prose naming a call is not executing it.
    """
    # EVERY TRACKED FILE, from git, rather than `*.py` under three directories. The old
    # form could not see the twelve home paths that `outputs/audit_reproducibility.csv`
    # was committing on every audit run -- a CSV, written by this gate's sibling, holding
    # exactly what this gate forbids. A glob over the directories somebody remembered is a
    # guess about where the repository is; `git ls-files` is the repository.
    offenders = []
    listing = subprocess.run(["git", "ls-files", "-z"], cwd=REPO,
                             capture_output=True, text=True)
    if listing.returncode != 0:
        return [_row("home paths confined to ystwin.paths", "none outside the resolver",
                     "REFUSED", False,
                     "git ls-files failed, so the set of tracked files is unknown; a gate "
                     "that cannot enumerate what it checks must not report a pass")]
    calls = []
    for relative in sorted(f for f in listing.stdout.split("\0") if f):
        path = REPO / relative
        if not path.is_file():
            continue
        try:
            body = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if (relative not in _LITERAL_PATH_EXEMPT
                and re.search(r"/Users/|/home/[a-z]", body)):
            offenders.append(relative)
        if (relative.endswith(".py") and relative not in _HOME_CALL_EXEMPT
                and re.search(r"Path\.home\(\)", body)):
            calls.append(relative)
    return [
        _row("no literal home path outside the resolver", "none",
             f"{len(offenders)} found", not offenders, ", ".join(sorted(offenders))),
        _row("no Path.home() call outside the resolver", "none",
             f"{len(calls)} found", not calls, ", ".join(sorted(calls))),
    ]


def check_dependencies_pinned() -> list[dict]:
    """Runtime dependencies must be pinned, not ranged."""
    text = (REPO / "pyproject.toml").read_text()
    block = re.search(r"^dependencies\s*=\s*\[(.*?)\]", text, re.S | re.M)
    if not block:
        return [_row("dependencies pinned", "a dependencies list", "none found", False)]
    entries = re.findall(r'"([^"]+)"', block.group(1))
    unpinned = [e for e in entries if "==" not in e]
    return [_row("runtime dependencies pinned with ==", "all",
                 f"{len(entries) - len(unpinned)}/{len(entries)}", not unpinned,
                 ", ".join(unpinned))]


def check_git_history() -> list[dict]:
    """The repository must have a commit."""
    result = subprocess.run(["git", "rev-list", "--count", "HEAD"],
                            cwd=REPO, capture_output=True, text=True)
    count = int(result.stdout.strip()) if result.returncode == 0 else 0
    return [_row("git history exists", ">= 1 commit", f"{count} commits", count >= 1)]


def check_stale_bytecode() -> list[dict]:
    """Report orphaned bytecode as a local hygiene hint, not a clone-independent check.

    The load-bearing invariant is checked by the module-reference audits.
    """
    stale = []
    for cache in SRC.rglob("__pycache__"):
        for compiled in cache.glob("*.pyc"):
            source = cache.parent / (compiled.name.split(".")[0] + ".py")
            if not source.exists():
                stale.append(source.stem)
    return [_row("no orphaned bytecode in src (local hint only)", "none",
                 f"{len(stale)} found", not stale, ", ".join(sorted(set(stale))))]


def check_docs_module_references() -> list[dict]:
    """Repository modules named in documents under ``docs/`` must exist."""
    docs = REPO / "docs"
    if not docs.exists():
        return []
    rows = []
    for document in sorted(docs.rglob("*.md")):
        text = document.read_text(encoding="utf-8", errors="replace")
        missing = []
        for name in sorted(set(re.findall(r"`([a-z_][a-z_0-9/]*\.py)`", text))):
            # Only paths rooted in this repository are existence claims about it. A
            # research report naming pytfa/thermo/metabolite.py or
            # flapjack_api/analysis/analysis.py is discussing another codebase, and
            # failing on those would make the check unusable in exactly the documents
            # that survey alternatives.
            if "/" in name and not name.startswith(_OWNED_PREFIXES):
                continue
            if not _module_exists(name):
                missing.append(name)
        relative = document.relative_to(REPO)
        rows.append(_row(f"modules referenced in {relative}", "all exist",
                         "all exist" if not missing else f"{len(missing)} missing",
                         not missing, ", ".join(missing)))
    return rows


_RETRACTION = re.compile(r"<!--\s*audit:retracted\s+(.+?)\s*-->")


def check_retracted_claims() -> list[dict]:
    """An explicitly retracted phrase must not stand in another current document.

    The author supplies the exact phrase in an audit:retracted marker; no claim is
    inferred from an unmarked number. Superseded documents retain historical context.
    """
    documents = sorted((REPO / "docs").rglob("*.md")) + [REPO / "README.md"]
    retractions: list[tuple[pathlib.Path, str]] = []
    for document in documents:
        for phrase in _RETRACTION.findall(document.read_text()):
            retractions.append((document, phrase))

    if not retractions:
        return [_row("retracted claims do not stand elsewhere", "any retraction marker",
                     "none declared", True,
                     "no <!-- audit:retracted ... --> marker anywhere; nothing to enforce")]

    rows = []
    for source, phrase in retractions:
        offenders = []
        for document in documents:
            if document == source or "superseded" in document.parts:
                continue
            if phrase.lower() in document.read_text().lower():
                offenders.append(document.relative_to(REPO).as_posix())
        rows.append(_row(
            f"retracted claim absent elsewhere: {phrase!r}",
            f"only in {source.relative_to(REPO).as_posix()} or docs/superseded/",
            "nowhere else" if not offenders else f"also in {offenders}",
            not offenders,
            "" if not offenders else
            f"{source.relative_to(REPO).as_posix()} retracts this and {offenders} still "
            "asserts it. Delete or rewrite it there, or move the passage into "
            "docs/superseded/ if it is being kept as history"))
    return rows


#: How far the architecture map's LOC figure may sit from the truth. It is written with a
#: leading "~" and moves with every edit, so an exact match would fail on a comment; the
#: module and script counts carry no tilde and are checked exactly.
_SCOPE_LOC_TOLERANCE = 0.05


def check_architecture_scope() -> list[dict]:
    """The five explicit counts in the architecture scope line must describe this tree."""
    doc = REPO / "docs" / "ARCHITECTURE.md"
    if not doc.exists():
        return [_row("architecture scope counts", "a scope line", "no ARCHITECTURE.md",
                     True, "nothing to check")]
    text = doc.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"Scope: every module under `src/ystwin/` \((\d[\d,]*) code modules \+ (\d[\d,]*) "
        r"package\s+`__init__\.py`, ~(\d[\d,]*) LOC\), every script under `scripts/` "
        r"\((\d[\d,]*)\), and the contracts the\s+(\d[\d,]*)-test suite",
        text)
    if not match:
        return [_row("architecture scope counts", "a parseable scope line", "REFUSED",
                     False,
                     "the scope line in docs/ARCHITECTURE.md no longer matches the shape "
                     "this check reads. Restore it or update the pattern -- a scope line "
                     "nobody can parse is one nobody is checking")]
    claimed = [int(g.replace(",", "")) for g in match.groups()]

    src = REPO / "src" / "ystwin"
    actual = [
        len([q for q in src.rglob("*.py") if q.name != "__init__.py"]),
        len(list(src.rglob("__init__.py"))),
        sum(len(q.read_text(encoding="utf-8", errors="replace").splitlines())
            for q in (REPO / "src").rglob("*.py")),
        len(list((REPO / "scripts").glob("*.py"))),
        _collected_test_count(),
    ]
    names = ["code modules", "package __init__.py", "src LOC", "scripts/*.py",
             "tests collected"]
    rows = []
    for name, said, is_ in zip(names, claimed, actual):
        if is_ is None:                       # the collector could not run
            rows.append(_row(f"architecture scope: {name}", said, "could not count", False,
                             "pytest collection unavailable here"))
            continue
        if name == "src LOC":
            ok = abs(is_ - said) <= _SCOPE_LOC_TOLERANCE * max(is_, 1)
            detail = (f"within {_SCOPE_LOC_TOLERANCE:.0%}" if ok else
                      f"the map says ~{said:,} and the tree holds {is_:,}")
        else:
            ok = is_ == said
            detail = "" if ok else f"the map says {said:,} and the tree holds {is_:,}"
        rows.append(_row(f"architecture scope: {name}", f"{said:,}", f"{is_:,}", ok, detail))
    return rows


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true", help="only print blocking records")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="write diagnostics here; default verification is read-only")
    args = parser.parse_args(argv)

    if not README.exists():
        print("no README.md to audit", file=sys.stderr)
        return 2
    text = README.read_text(encoding="utf-8")

    rows: list[dict] = []
    rows += check_git_history()
    rows += check_no_absolute_home_paths()
    rows += check_dependencies_pinned()
    rows += check_stale_bytecode()
    rows += check_module_references(text)
    rows += check_docs_module_references()
    rows += check_documented_commands(text)
    rows += check_output_tables(text)
    rows += check_test_count(text)
    rows += check_retracted_claims()
    rows += check_architecture_scope()
    claim_rows = audit_claims(REPO)
    rows += claim_rows

    report = None
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        report = args.output_dir / "audit_claims.csv"
        pd.DataFrame(rows).to_csv(report, index=False)

    failures = [row for row in rows if row["blocking"]]
    width = 72
    print("=" * width)
    print("Claim audit: repository references and declared scientific evidence")
    print("=" * width)
    for row in rows:
        if args.quiet and not row["blocking"]:
            continue
        mark = "FAIL" if row["blocking"] else {"PASS": "ok  ", "SKIP": "skip"}.get(row["status"], row["status"])
        print(f"  [{mark}] {row['check']}")
        if row["blocking"]:
            print(f"         expected {row['expected']}, got {row['actual']}")
            if row["detail"]:
                print(f"         {row['detail']}")

    # Say how many checks compared a number, because "N/N checks pass" has been quoted as
    # "N claims traceable" and those are different statements: most repository checks confirm
    # a reference resolves, which says nothing about the figure quoted from it.
    markers = [row for row in claim_rows if row["record_type"] == "marked_value"]
    compared = sum(row.get("comparison_status") in {"PASS", "FAIL"} for row in markers)
    passed = sum(row["status"] == "PASS" for row in rows)
    skipped = sum(row["status"] == "SKIP" for row in rows)
    refusals = sum(row.get("verdict") == "REFUSAL_CONFIRMED" for row in claim_rows)
    print()
    destination = f" -> {report}" if report else " (read-only; no report written)"
    print(f"  {passed} PASS, {skipped} SKIP, {len(failures)} blocking of {len(rows)} records{destination}")
    print(f"  {compared} marked cell comparisons; {refusals} verified scientific refusals.")
    print("  Cell agreement and reference resolution are not scientific support; see each claim verdict.")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
