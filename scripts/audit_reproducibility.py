"""Ask whether a stranger could reproduce this repository's results.

``audit_claims.py`` asks whether the README describes the code that exists. This asks the
next question, which is harder and less often asked: if someone who has never met this
project clones it, installs it, and runs the scripts, do they get these numbers?

The failure modes are all silent, and none of them show up in a test run.

* A path that only resolves on the author's machine. The results still appear; they are
  just not reproducible, and nothing says so.
* A dependency pinned to one version and installed as another. The pin exists precisely
  to stop that, and a pin nobody verifies is decoration.
* A different LP solver. ``outputs/d1_capacity_sweep_ec.csv`` holds FVA widths, and a
  degenerate LP has many optima -- GLPK reports one vertex and Gurobi another. Both are
  correct. Only one matches the committed table.
* A result file nobody else can see, because it was never committed.
* A table older than the code that writes it, which then documents behaviour the library
  no longer has.

Absent data is **not** a failure. This repository is meant to work on a machine that has
never seen a plate reader, and every real-data path resolves to ``None`` and skips. What
*is* a failure is an asset that resolves to somewhere outside the checkout with no
environment variable pointing there, because that is a path that works here and nowhere
else.

Read-only by default; an explicit temporary ``--output-dir`` receives CSV and JSON reports.
Exits non-zero if any check fails. Changed identities require investigation, not a claim
that the scientific values are wrong. Blind restamping is refused.

Usage: python scripts/audit_reproducibility.py [--quiet] [--output-dir TEMP] [--offline-suite]

``--offline-suite`` runs the test suite twice in subprocesses -- once as configured, once
with every data-locating variable pointed at a directory that does not exist -- and is off
by default because it costs several minutes. It is the check that actually proves the
"works without the wet-lab data" claim, so run it before a release.
"""
from __future__ import annotations

import argparse
import ast
import fnmatch
import importlib.metadata
import os
import pathlib
import json
import re
import subprocess
import sys
import tomllib
import uuid

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

pd = importlib.import_module("pandas")
_artifacts = importlib.import_module("ystwin.artifacts")
ArtifactError = _artifacts.ArtifactError
REGISTRY_PATH = _artifacts.REGISTRY_PATH
audit_registry = _artifacts.audit_registry
artifact_producer_digest = _artifacts.producer_digest
read_json = _artifacts.read_json

__all__ = [
    "REPO",
    "check_asset_resolvers",
    "check_dependency_pins",
    "check_env_vars_documented",
    "check_offline_suite",
    "check_output_freshness",
    "check_output_tables",
    "check_outputs_tracked",
    "library_dependencies",
    "check_solver",
    "check_writer_map_matches_doc",
    "describe_resolvers",
    "writers_from_doc",
    "writers_from_scripts",
]

REPO = pathlib.Path(__file__).resolve().parents[1]

#: Filesystem timestamps have a resolution and a checkout writes many files at once. A
#: fresh clone gives every file effectively the same mtime, and the freshness check must
#: not read that as a stale table.
FRESHNESS_TOLERANCE_S = 2.0


def _row(check: str, expected: object, actual: object, status: object, detail: str = "") -> dict:
    """One CSV row, in the shape ``audit_claims.py`` established.

    ``status`` may be a bool, for the usual pass/fail, or the literal string ``"SKIP"``
    for a check that could not run here -- an absent optional dependency, a repository
    with no git. A skip is reported and does not fail the audit; it is distinct from a
    pass because a check that never ran has proved nothing.
    """
    if isinstance(status, str):
        text = status
    else:
        text = "PASS" if status else "FAIL"
    return {
        "check": check,
        "expected": str(expected),
        "actual": str(actual),
        "status": text,
        "detail": detail,
    }


# ---------------------------------------------------------------------------
# where the data lives
# ---------------------------------------------------------------------------


def describe_resolvers(paths_module) -> list[dict]:
    """Every asset resolver in :mod:`ystwin.paths`, its variable, and where it lands.

    The environment variable is recovered by watching ``paths._resolve`` actually get
    called rather than by pattern-matching the source. The call is the authority: a
    resolver that consults a variable does so by passing its name to ``_resolve``, and a
    resolver that consults none -- ``data_dir``, ``outputs_dir`` -- reports none, which is
    the correct answer for a vendored directory and for a destination.

    Args:
        paths_module: The imported :mod:`ystwin.paths`.

    Returns:
        One record per resolver with ``name``, ``env_var`` (``None`` if it consults no
        variable), ``resolved`` (a path or ``None``) and ``defaults``.
    """
    import inspect

    seen: list[tuple[str, tuple]] = []

    def spy(env_var, *defaults):
        seen.append((env_var, defaults))
        return None

    records = []
    original = paths_module._resolve
    for name in getattr(paths_module, "__all__", []):
        function = getattr(paths_module, name, None)
        if not inspect.isfunction(function):
            continue
        signature = inspect.signature(function)
        required = [
            p for p in signature.parameters.values()
            if p.default is inspect.Parameter.empty
            and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
        ]
        if required:
            # require() and resolve_or_exit() take a resolved path; they are not
            # resolvers themselves.
            continue
        if function is getattr(paths_module, "outputs_dir", None):
            default = paths_module.REPO_ROOT / "outputs"
            records.append({"name": name, "env_var": "YSTWIN_OUTPUTS",
                            "resolved": pathlib.Path(os.environ.get("YSTWIN_OUTPUTS") or default),
                            "defaults": [str(default)]})
            continue
        seen.clear()
        paths_module._resolve = spy
        try:
            function()
        except Exception:  # a resolver that cannot even be called is reported below
            pass
        finally:
            paths_module._resolve = original
        env_var = seen[0][0] if seen else None
        defaults = seen[0][1] if seen else ()
        try:
            resolved = function()
        except Exception as exc:
            resolved = f"error: {type(exc).__name__}: {exc}"
        records.append({
            "name": name,
            "env_var": env_var,
            "resolved": resolved,
            "defaults": [str(d) for d in defaults],
        })
    return records


def _portable(resolved: pathlib.Path, repo: pathlib.Path) -> str:
    """How to name a resolved location in a file that gets committed.

    **The rows this function feeds are written to `outputs/audit_reproducibility.csv`,
    which is tracked**, and until 2026-08-31 they carried the absolute path -- twelve of
    them, naming one person's home directory, their Desktop, their Downloads folder and
    the private workbook in it. `.gitignore` refuses to commit the plate workbooks
    precisely because publishing a path into git history costs a rewrite to remove; this
    audit was committing the paths themselves, once per run.

    What the check is actually about is whether an asset resolves *portably*, and the
    portable name for each case is not the absolute path:

    * inside the checkout -> the repo-relative path, which is the same on every machine
    * found by an override -> the variable, which is the reproducible handle; the path it
      points at is machine-specific by definition and says nothing a reader can act on
    * found by convention outside the checkout -> the home-relative form, because that row
      is a FAILURE whose whole content is "this location exists here and nowhere else",
      and the username is not part of that statement

    Found by the 2026-08-31 audit, which also found that `audit_claims.py`'s home-path gate
    could not see this: it scans `*.py` under src/, scripts/ and tests/, and this is a CSV.
    """
    try:
        return str(resolved.relative_to(repo)) + " (in-repo)"
    except ValueError:
        pass
    try:
        return "~/" + str(resolved.relative_to(pathlib.Path.home()))
    except ValueError:
        return str(resolved)


def check_asset_resolvers(records: list[dict], repo: pathlib.Path,
                          environ: dict | None = None) -> list[dict]:
    """Report where each asset resolves, and fail only on the unportable case.

    A resolver returning ``None`` is a documented state, not a failure: it means the
    machine does not have that dataset and every consumer will skip. A resolver landing
    *outside* the checkout while its variable is unset is the failure, because the path
    was found by a convention -- a sibling directory, a folder on someone's Desktop --
    that no one else's machine satisfies. The result computed from it will differ
    elsewhere with nothing in the repository to explain why.

    Args:
        records: Output of :func:`describe_resolvers`.
        repo: The checkout root; anything under it is portable by definition.
        environ: Environment to read overrides from. Defaults to the real one.
    """
    environ = os.environ if environ is None else environ
    rows = []
    for record in records:
        name, env_var, resolved = record["name"], record["env_var"], record["resolved"]
        label = f"asset resolver: {name}()"
        variable = env_var or "none"
        if isinstance(resolved, str):
            rows.append(_row(label, "resolves or None", resolved, False,
                             "the resolver itself raised"))
            continue
        if resolved is None:
            rows.append(_row(label, f"a path or None (override: {variable})", "absent",
                             True, f"not on this machine; set {variable} to supply it"))
            continue
        inside = pathlib.Path(resolved).is_relative_to(repo)
        overridden = bool(env_var and environ.get(env_var))
        if inside:
            rows.append(_row(label, f"a path or None (override: {variable})",
                             _portable(pathlib.Path(resolved), repo), True))
        elif overridden:
            # The path itself is deliberately NOT recorded: it is machine-specific by
            # definition, this row is committed, and `_portable` says why.
            rows.append(_row(label, f"a path or None (override: {variable})",
                             f"resolved via {env_var}", True,
                             "location supplied by the variable; not recorded here "
                             "because it is specific to this machine"))
        else:
            rows.append(_row(
                label, f"in-repo, or located by {variable}",
                _portable(pathlib.Path(resolved), repo), False,
                f"resolved outside the checkout by convention with {variable} unset; "
                f"this path exists here and on no other machine"))
    return rows


def check_env_vars_documented(records: list[dict], doc_text: str,
                              doc_name: str = "docs/REPRODUCING.md") -> list[dict]:
    """Every variable the resolvers consult must be named in the reproduction guide.

    A variable nobody documents is a variable nobody sets, which leaves the sibling-
    directory convention as the only way in -- and that convention is what makes a
    checkout personal.
    """
    rows = []
    for record in sorted(records, key=lambda r: r["name"]):
        env_var = record["env_var"]
        if not env_var:
            continue
        present = env_var in doc_text
        rows.append(_row(f"override documented: {env_var}", f"named in {doc_name}",
                         "named" if present else "MISSING", present))
    return rows


# ---------------------------------------------------------------------------
# the environment the numbers were produced in
# ---------------------------------------------------------------------------

_REQUIREMENT = re.compile(r"^\s*([A-Za-z0-9][A-Za-z0-9._-]*)\s*(\[[^\]]*\])?\s*(.*)$")


def _requirement_groups(pyproject_text: str) -> tuple[list[tuple[str, list[str], bool]], str]:
    """The requirement lists in a ``pyproject.toml``, read as TOML rather than as text.

    Returns ``[(group, entries, required)]`` and an error string, empty when the file
    parsed. PEP 621 says exactly which keys hold requirements -- ``project.dependencies``
    and each list under ``project.optional-dependencies`` -- so there is no guessing to do
    and no other key to accidentally include.

    **This was a regex over the whole file with a blocklist of key names**, and it broke
    the day a ``[tool.ruff]`` section was added: it matched every ``name = [...]`` anywhere
    in the document, so ``select = ["E", "F"]`` was read as two dependencies named ``E``
    and ``F``, ``ignore = ["E501"]`` as a third, and ``extend-exclude = ["archive"]`` as a
    fourth -- four FAIL rows saying a linter rule was not pinned to a version. The
    blocklist named seven keys it had to skip and would have needed an eighth, a ninth and
    a tenth for every future tool section, each one added only after it had already fired.

    The failure direction is the bad one for a gate: adding an unrelated config section
    turns the audit red, so the pressure is to widen the blocklist until it stops
    complaining, and a blocklist wide enough to stay quiet is one that has stopped reading
    the file. The format already knows which keys are requirements. Ask it.
    """
    try:
        data = tomllib.loads(pyproject_text)
    except tomllib.TOMLDecodeError as exc:
        return [], f"pyproject.toml is not valid TOML: {exc}"
    project = data.get("project", {})
    groups: list[tuple[str, list[str], bool]] = [
        ("runtime", list(project.get("dependencies", [])), True)]
    optional = project.get("optional-dependencies", {})
    for name in sorted(optional):
        groups.append((name, list(optional[name]), False))
    return groups, ""


def check_dependency_pins(pyproject_text: str, installed) -> list[dict]:
    """Runtime dependencies must be pinned with ``==`` *and* be what is installed.

    ``audit_claims.py`` checks the first half. The second half is the one that matters at
    run time: a pin says which version produced the committed tables, and an environment
    holding a different version produces different tables while the pin still reads
    correct. That divergence is exactly what pinning exists to prevent and it is invisible
    from the file alone.

    Optional groups are checked for the pin and, when the package is absent, reported as
    skipped rather than failed -- ``[thermo]`` is documented as optional and its tests
    skip without it.

    Args:
        pyproject_text: Contents of ``pyproject.toml``.
        installed: Callable taking a distribution name and returning its installed
            version, or ``None``/raising :class:`importlib.metadata.PackageNotFoundError`
            when it is not installed.
    """
    rows = []
    groups, error = _requirement_groups(pyproject_text)
    if error:
        return [_row("dependency pins readable", "parseable pyproject.toml", "REFUSED",
                     False, error)]

    for group, entries, required in groups:
        for entry in entries:
            parsed = _REQUIREMENT.match(entry)
            if not parsed:
                rows.append(_row(f"dependency pinned: {entry}", "name==version",
                                 "unparseable", False))
                continue
            name, _extras, rest = parsed.groups()
            if not rest.startswith("=="):
                rows.append(_row(f"dependency pinned: {name} [{group}]", "name==version",
                                 entry, False,
                                 "a range lets a fresh install disagree with outputs/"))
                continue
            pinned = rest[2:].split(";")[0].strip()
            try:
                actual = installed(name)
            except Exception:
                actual = None
            if actual is None:
                rows.append(_row(
                    f"pin matches installed: {name} [{group}]", pinned, "not installed",
                    False if required else "SKIP",
                    "a runtime dependency the environment does not have" if required
                    else f"optional group [{group}]; its consumers skip"))
                continue
            rows.append(_row(f"pin matches installed: {name} [{group}]", pinned, actual,
                             actual == pinned,
                             "" if actual == pinned else
                             "the installed version is not the one the tables were made with"))
    return rows


def check_python_version(pyproject_text: str, version_info=None) -> list[dict]:
    """The running interpreter must satisfy ``requires-python``."""
    version_info = sys.version_info if version_info is None else version_info
    match = re.search(r'requires-python\s*=\s*"([^"]+)"', pyproject_text)
    if not match:
        return [_row("requires-python declared", "a constraint", "none", False)]
    constraint = match.group(1).strip()
    floor = re.search(r">=\s*(\d+)\.(\d+)", constraint)
    running = f"{version_info[0]}.{version_info[1]}.{version_info[2]}"
    if not floor:
        return [_row("interpreter satisfies requires-python", constraint, running, "SKIP",
                     "constraint is not a simple lower bound")]
    wanted = (int(floor.group(1)), int(floor.group(2)))
    ok = tuple(version_info[:2]) >= wanted
    return [_row("interpreter satisfies requires-python", constraint, running, ok)]


def check_solver(solver_name: str | None, available: list[str] | None = None,
                 expected: str = "glpk") -> list[dict]:
    """Which LP solver ``cobra`` resolves to, and whether it is the documented one.

    Not a correctness question -- every solver here returns *an* optimum. It is a
    reproducibility question: an LP with a degenerate optimum has a face of solutions and
    each solver picks a different vertex of it, so the FVA widths committed under GLPK are
    not the widths a machine with Gurobi installed will print. The committed tables name
    no solver, so the only way anyone notices is this check.

    Args:
        solver_name: ``cobra.Configuration().solver.__name__``, or ``None`` if cobra is
            not importable.
        available: Solver names cobra can see, for the detail column.
        expected: The interface the committed tables were produced with.
    """
    if solver_name is None:
        return [_row("cobra solver", expected, "cobra not importable", "SKIP",
                     "install the runtime dependencies to check this")]
    ok = expected in solver_name.lower()
    detail = f"available: {', '.join(sorted(available))}" if available else ""
    if not ok:
        detail = (f"the committed FVA widths are GLPK's; a different solver returns a "
                  f"different vertex of the same degenerate optimum. {detail}")
    return [_row("cobra solver", f"an interface containing {expected!r}", solver_name,
                 ok, detail)]


# ---------------------------------------------------------------------------
# what is in outputs/, and whether anyone else can see it
# ---------------------------------------------------------------------------


def _git(repo: pathlib.Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)


def check_outputs_tracked(repo: pathlib.Path) -> list[dict]:
    """Every result in ``outputs/`` must be tracked, or deliberately ignored.

    A result on one person's disk is not a result. The repository's ``.gitignore`` states
    the policy -- CSV tables in, binaries and figures out -- so the invariant is that every
    file is on one side of it or the other. A file that is neither tracked nor ignored is
    one that was produced, was not committed, and nobody has noticed.
    """
    outputs = repo / "outputs"
    if not outputs.is_dir():
        return [_row("outputs/ exists", "a directory", "missing", False)]
    if _git(repo, "rev-parse", "--git-dir").returncode != 0:
        return [_row("outputs tracked by git", "all", "not a git repository", "SKIP")]

    rows = []
    tracking = _artifacts.git_tracking(repo)
    current = None
    if (repo / REGISTRY_PATH).is_file():
        try:
            registry = _artifacts.load_registry(repo)
            current = {name for name, item in registry["artifacts"].items() if item["lifecycle"] == "current"}
        except (OSError, ValueError, KeyError, TypeError):
            pass
    for relative in sorted(tracking["present"]):
        if relative in tracking["tracked"]:
            rows.append(_row(f"output tracked: {relative}", "tracked", "tracked", True))
            continue
        ignored = relative in tracking["ignored"]
        if current is not None and relative not in current:
            rows.append(_row(f"local artifact visibility: {relative}", "explicit local decision", "ignored local" if ignored else "untracked local", "SKIP",
                             "Not retained by Git and not declared current; the registry reports its local decision separately."))
            continue
        rows.append(_row(f"output tracked: {relative}", "tracked or ignored", "ignored" if ignored else "UNTRACKED", ignored,
                         "" if ignored else "declared current but not committed; retained ownership remains blocking"))
    return rows


def check_output_tables(repo: pathlib.Path) -> list[dict]:
    """Every ``outputs/*.csv`` must be non-empty and parse with at least one data row.

    A zero-byte or header-only table reads, in a diff and in a citation, exactly like a
    result. It is the shape a crashed script leaves behind.
    """
    outputs = repo / "outputs"
    if not outputs.is_dir():
        return [_row("outputs/ exists", "a directory", "missing", False)]
    rows = []
    try:
        selected, failed_exports = _artifacts.validation_scope(repo)
    except (OSError, ValueError, KeyError, TypeError) as error:
        selected, failed_exports = None, {}
        rows.append(_row("artifact validation scope", "readable registry", "invalid", False, str(error)))
    for path in sorted(outputs.rglob("*.csv")):
        relative = path.relative_to(repo).as_posix()
        if selected is not None and relative not in selected:
            continue
        if path.stat().st_size == 0:
            rows.append(_row(f"table parses: {relative}", ">= 1 data row", "empty file", False))
            continue
        try:
            frame = pd.read_csv(path)
        except Exception as exc:
            rows.append(_row(f"table parses: {relative}", ">= 1 data row",
                             f"{type(exc).__name__}", False, str(exc)[:200]))
            continue
        count = len(frame)
        rows.append(_row(f"table parses: {relative}", ">= 1 data row", f"{count} rows",
                         count >= 1,
                         "" if count >= 1 else "headers only; the script wrote nothing"))
    for path in sorted(outputs.rglob("*.json")):
        relative = path.relative_to(repo).as_posix()
        if relative in failed_exports:
            rows.append(_row(f"historical failed export: {relative}", "original error preserved", failed_exports[relative], "SKIP",
                             "Not parsed as a current prediction; origin bytes, failure evidence and original runtime are validated by the registry."))
            continue
        if selected is not None and relative not in selected:
            continue
        try:
            read_json(path)
        except (OSError, ValueError) as error:
            rows.append(_row(f"JSON parses: {relative}", "valid unambiguous JSON", "invalid", False, str(error)))
        else:
            rows.append(_row(f"JSON parses: {relative}", "valid unambiguous JSON", "valid", True))
    return rows


# ---------------------------------------------------------------------------
# which script writes which file, and whether the table predates it
# ---------------------------------------------------------------------------

_WRITE_CALLS = {"to_csv", "to_excel", "savez", "savez_compressed", "save", "savefig"}
"""Library methods that put a file on disk, naming it in their FIRST ARGUMENT.

``savefig`` was missing, and its absence hid all twelve figures: every ``fig*.svg``,
``fig*.png`` and ``atp_sensor_*.png`` in ``outputs/`` was reported as having no writer at
all -- "nobody can regenerate it" -- while ``make_figures.py`` regenerates them on demand.
A false alarm from an integrity gate is worse than a quiet one, because the honest
response to it is to stop believing the gate.

It happened again, at thirty-seven times the scale, and adding a name to this set could
not have fixed it: the native-evidence family writes through ``open("x")``, which names
its file in the RECEIVER, and through local ``_write(path, payload)`` helpers, which name
it in the caller. Those two shapes are :data:`_PATH_WRITE_METHODS` and
:func:`_write_helpers`; this set stays what it says it is."""

_PATH_WRITE_METHODS = {"write_text", "write_bytes"}
"""``pathlib`` writes that name their destination in the RECEIVER rather than an argument.

``(output / "ledger.json").open("xb")`` is the same shape: the filename is to the left of
the dot, where ``_patterns_in`` used to look only to the right of the open bracket."""

_WRITE_MODES = frozenset("wxa")
"""``open`` mode letters that create or replace a file. ``x`` -- exclusive create -- is how
every write-once evidence directory in this repository is written, and a mode string is
the authoritative statement of intent, so it is read rather than guessed at."""

_READ_CALLS = {"read_csv", "read_excel", "load", "read_parquet", "imread"}
"""Calls that make a script a CONSUMER of a file rather than a producer of one.

Symmetric with :data:`_WRITE_CALLS`, and the reason is a defect this gate could not see.
``run_calibration.py`` reads ``outputs/sensor_characterisation.csv`` and writes
``outputs/panel_calibration.csv``. When the fourth plate rewrote the first, the second
described a superseded input -- and the freshness check passed it, because the check
modelled script -> output and never output -> output. Found by the 2026-08-30 audit, which
verified the stale table was byte-identical to the one built from the 251-row predecessor.
"""


def _pattern_from(node: ast.AST) -> str | None:
    """The filename a write call's destination expression names, as an fnmatch pattern.

    ``OUT / "power_analysis.csv"`` gives ``power_analysis.csv``; ``OUT / f"g1_{tag}.csv"``
    gives ``g1_*.csv``, because the tag is a plate name not known until the script runs.
    """
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left, right = _pattern_from(node.left), _pattern_from(node.right)
        return f"{left}/{right}" if left and right else right
    if isinstance(node, ast.Call):
        name = getattr(node.func, "attr", getattr(node.func, "id", ""))
        if name in {"Path", "PurePosixPath"} and node.args:
            return _pattern_from(node.args[0])
        if name == "outputs_dir":
            return "outputs"
        if name == "data_dir":
            return "data"
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                parts.append(value.value)
            else:
                parts.append("*")
        return "".join(parts)
    return None


def _writing_mode(node: ast.Call, index: int) -> bool:
    """Whether the argument at ``index`` is a literal file mode that creates or replaces.

    ``path.open("xb")`` puts the mode first; ``open(path, "w")`` and
    ``zipfile.ZipFile(out / "model_source.zip", "x")`` put it second; ``mode=`` says the
    same thing by name. No mode at all means read, and a computed one is not evidence.

    A mode is a short string over ``rwxab+t``, and that shape is what stops an ordinary
    second argument from being read as one -- the alternative, a list of the constructors
    that take a mode, is the kind of list this scan exists to avoid keeping.
    """
    value = None
    if len(node.args) > index:
        argument = node.args[index]
        if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
            value = argument.value
    for keyword in node.keywords:
        if keyword.arg == "mode" and isinstance(keyword.value, ast.Constant):
            if isinstance(keyword.value.value, str):
                value = keyword.value.value
    if not value or len(value) > 3 or not set(value) <= set("rwxab+t"):
        return False
    return bool(set(value) & _WRITE_MODES)


def _write_destinations(tree: ast.AST, helpers: dict[str, int],
                        restricted: bool = False) -> list[ast.AST]:
    """Every expression a write in this tree names its destination with.

    Four shapes, and only the first was ever read. ``frame.to_csv(OUT / "power.csv")``
    names the file in its first argument. ``(output / "ledger.json").open("xb")`` names it
    in the RECEIVER -- to the left of the dot -- and so does ``path.write_bytes(data)``.
    ``_write(output / "ledger.json", payload)`` names it in whichever argument the local
    helper opens, which is what :func:`_write_helpers` works out. And
    ``zipfile.ZipFile(out / "model_source.zip", "x")`` names it first with the mode second.

    Args:
        tree: Module or function to scan.
        helpers: Write helpers, as ``{name: index of the destination argument}``.
        restricted: Only count writes that can be attributed to a parameter of the
            enclosing function -- direct filesystem writes, and helper calls made as a
            bare name or through ``self``. Used while DETECTING helpers, where
            ``stream.write(content)`` on a file handle would otherwise look like a call to
            the ``ArtifactStore.write`` method that happens to share its name.
    """
    found: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute):
            name = node.func.attr
            through_self = isinstance(node.func.value, ast.Name) and node.func.value.id == "self"
            if name in helpers and (through_self or not restricted):
                if helpers[name] < len(node.args):
                    found.append(node.args[helpers[name]])
                continue
            if name in _WRITE_CALLS and node.args:
                found.append(node.args[0])
                continue
            if name in _PATH_WRITE_METHODS:
                found.append(node.func.value)
                continue
            if name == "open" and _writing_mode(node, 0):
                found.append(node.func.value)
                continue
        elif isinstance(node.func, ast.Name):
            if node.func.id == "open" and node.args and _writing_mode(node, 1):
                found.append(node.args[0])
                continue
            if node.func.id in helpers:
                if helpers[node.func.id] < len(node.args):
                    found.append(node.args[helpers[node.func.id]])
                continue
        # The generic (destination, mode) constructor. Only where the first argument spells
        # the filename out, so a second argument that merely looks like a mode cannot claim
        # a file: at that point the mode is confirming a destination, not proposing one.
        if restricted or not node.args or not _writing_mode(node, 1):
            continue
        if _pattern_from(node.args[0]) is not None:
            found.append(node.args[0])
    return found


def _parameters_reached(node: ast.AST, assignments: dict[str, list[ast.AST]],
                        parameters: set[str]) -> set[str]:
    """Which of ``parameters`` a destination expression is built from, through locals.

    ``ArtifactStore.write_bytes`` opens ``safe_path(self.root, self.directory / name)``
    three assignments away from its ``name`` argument, so asking only what the receiver is
    called answers ``path`` and learns nothing.
    """
    reached: set[str] = set()
    seen: set[str] = set()
    queue: list[ast.AST] = [node]
    while queue:
        for sub in ast.walk(queue.pop()):
            if not isinstance(sub, ast.Name) or sub.id in seen:
                continue
            seen.add(sub.id)
            if sub.id in parameters:
                reached.add(sub.id)
            queue.extend(assignments.get(sub.id, ()))
    return reached


def _destination_parameter(function: ast.AST, helpers: dict[str, int]) -> int | None:
    """Which argument a caller of ``function`` names the file with, if exactly one.

    Answered only when one parameter reaches every write in the body. ``_write_once(path,
    root, payload)`` builds its destination from BOTH ``root`` and ``path``, and a helper
    map that guessed between them would attribute files to a script by coin flip -- worse
    than the gap, because it would be believed. Returns the index AS CALLED, so a method's
    ``self`` is already discounted.
    """
    arguments = getattr(function, "args", None)
    if arguments is None:
        return None
    parameters = [a.arg for a in (*arguments.posonlyargs, *arguments.args)]
    if not parameters:
        return None
    offset = 1 if parameters[0] in ("self", "cls") else 0
    named = set(parameters[offset:])
    assignments: dict[str, list[ast.AST]] = {}
    for node in ast.walk(function):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                assignments.setdefault(target.id, []).append(node.value)
    candidates: set[str] = set()
    for destination in _write_destinations(function, helpers, restricted=True):
        # A write whose filename is written out here is the function's own, not the
        # caller's; `_patterns_in` reads it directly and there is nothing to forward.
        if _pattern_from(destination) is not None:
            continue
        candidates |= _parameters_reached(destination, assignments, named)
    if len(candidates) != 1:
        return None
    return parameters.index(candidates.pop()) - offset


def _write_helpers(tree: ast.AST) -> dict[str, int]:
    """Locally defined functions that write whatever file their caller names.

    ``_write(path, payload)``, ``write_freeze(result, path)`` and
    ``ArtifactStore.write(self, name, value, role)`` are each one ``open("x")`` behind a
    name the caller supplies, and the whole native-evidence family writes through one of
    them. Read as a plain call the destination is a bare parameter, so the scan gave up and
    453 committed files were reported as having no writer at all -- the same false alarm
    ``savefig`` once caused for the twelve figures, at thirty-seven times the size.

    Names are collected per module, not globally: ``write`` means ``ArtifactStore.write``
    inside ``native_population_development.py`` and a file handle everywhere else.

    Returns:
        ``{"_write": 0, "write": 0}`` -- helper name to the destination's argument index.
    """
    functions = [n for n in ast.walk(tree)
                 if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    helpers: dict[str, int] = {}
    # A helper may be written in terms of another (`write` calls `write_bytes`), so this
    # runs to a fixpoint rather than once.
    for _ in range(len(functions) + 1):
        found = False
        for function in functions:
            if function.name in helpers:
                continue
            index = _destination_parameter(function, helpers)
            if index is not None and index >= 0:
                helpers[function.name] = index
                found = True
        if not found:
            break
    return helpers


def _destination_bindings(tree: ast.AST) -> dict[str, list[str]]:
    """Filenames each variable can hold, for a write whose destination is a bare name.

    ``destination = OUT / "result.csv"`` is one, and ``viz/figures.py`` writes that way.
    Both order-robustness scripts write another way entirely --

        REPORT_NAME = "partial_order_synthetic_demo.json"
        paths = (args.output_dir / REPORT_NAME, args.output_dir / SUMMARY_NAME)
        for path, payload in zip(paths, payloads):
            with path.open("x", ...) as stream:

    -- where the destination is a loop variable over a tuple assembled three statements
    earlier out of two module constants. Nothing about that is computed at run time; it is
    simply spread out, and reading only the statement the write is in loses it.

    Assignments are read in source order so a name can be built from an earlier one, and a
    name assigned twice keeps both patterns: widening attribution costs an extra row, and
    narrowing it costs a false "nobody can regenerate this".
    """
    bindings: dict[str, list[str]] = {}

    def patterns_of(node: ast.AST) -> list[str]:
        if isinstance(node, ast.Name):
            return list(bindings.get(node.id, []))
        if isinstance(node, (ast.Tuple, ast.List)):
            return [p for element in node.elts for p in patterns_of(element)]
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            left, right = patterns_of(node.left), patterns_of(node.right)
            return [f"{a}/{b}" for a in left for b in right] if left else right
        if isinstance(node, (ast.ListComp, ast.SetComp, ast.GeneratorExp)):
            # `[output / name for name in (REPORT_NAME, SUMMARY_NAME, OBSERVATIONS_NAME)]`
            # -- the same three names, one syntax further from the write.
            for generator in node.generators:
                if isinstance(generator.target, ast.Name):
                    record(generator.target.id, patterns_of(generator.iter))
            return patterns_of(node.elt)
        pattern = _pattern_from(node)
        return [pattern] if pattern else []

    def record(name: str, values: list[str]) -> None:
        held = bindings.setdefault(name, [])
        held.extend(v for v in values if v and v not in held)

    ordered = sorted((n for n in ast.walk(tree) if isinstance(n, (ast.Assign, ast.For))),
                     key=lambda n: (n.lineno, n.col_offset))
    for node in ordered:
        if isinstance(node, ast.Assign):
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                record(node.targets[0].id, patterns_of(node.value))
            continue
        source = node.iter
        zipped = (isinstance(source, ast.Call)
                  and getattr(source.func, "id", getattr(source.func, "attr", None)) == "zip")
        if isinstance(node.target, ast.Name):
            record(node.target.id, patterns_of(source.args[0] if zipped else source))
        elif isinstance(node.target, (ast.Tuple, ast.List)) and zipped:
            for element, iterated in zip(node.target.elts, source.args):
                if isinstance(element, ast.Name):
                    record(element.id, patterns_of(iterated))
    return bindings


def _shadowed_names(tree: ast.AST) -> set[int]:
    """Name nodes a function parameter shadows, by node identity.

    Bindings are read across the whole module, which widens attribution on purpose -- but
    a parameter is not the module's variable of the same name.
    ``run_native_reconciliation.py`` binds ``path = root / FROZEN_CHECKPOINT`` at module
    level to READ the frozen checkpoint, and its ``_json(path, payload)`` helper opens a
    parameter that happens to be called ``path`` too. Resolving the second through the
    first made the reconciliation script the writer of the checkpoint it reads -- a false
    "you can regenerate this" over a frozen input, which is the worst direction to be
    wrong in.
    """
    shadowed: set[int] = set()
    for function in ast.walk(tree):
        if not isinstance(function, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        arguments = function.args
        names = {a.arg for a in (*arguments.posonlyargs, *arguments.args,
                                 *arguments.kwonlyargs)}
        for node in ast.walk(function):
            if isinstance(node, ast.Name) and node.id in names:
                shadowed.add(id(node))
    return shadowed


def _bound_patterns(node: ast.AST, bindings: dict[str, list[str]],
                    shadowed: set[int] | None = None) -> list[str]:
    """The filenames a destination expression can name, resolving module-level bindings."""
    if isinstance(node, ast.Name):
        if shadowed is not None and id(node) in shadowed:
            return []
        return bindings.get(node.id, [])
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left = _bound_patterns(node.left, bindings, shadowed)
        right = _bound_patterns(node.right, bindings, shadowed)
        return [f"{a}/{b}" for a in left for b in right] if left else right
    pattern = _pattern_from(node)
    return [pattern] if pattern else []


def _patterns_in(tree: ast.AST, helpers: dict[str, int] | None = None) -> list[str]:
    """Every output filename one module's write calls name.

    Three forms have to be resolved, not one. ``frame.to_csv(OUT / "power_analysis.csv")``
    names its destination inline, and ``_pattern_from`` reads it straight off the call.
    ``(output / "ledger.json").open("xb")`` names it in the receiver, and
    ``_write(output / "ledger.json", payload)`` in a helper's argument -- see
    :func:`_write_destinations`. But ``viz/figures.py`` binds the path first --

        svg = out_dir / f"{stem}.svg"
        fig.savefig(svg, ...)

    -- and a walk that only looks at the call sees a bare ``Name`` and gives up. That is
    how all twelve committed figures came to be reported as having no writer: the pattern
    was three lines above the call rather than inside it.

    Bindings are collected across the whole module rather than per scope. The cost of a
    name reused in two functions is an extra pattern, which only ever widens attribution;
    the cost of missing one is a false "nobody can regenerate this".

    Args:
        tree: The parsed module.
        helpers: Write helpers visible here, defaulting to this module's own. The caller
            passes the union over everything a script reaches, because
            ``native_population_scoring.py`` writes through the ``ArtifactStore`` defined
            in ``native_population_development.py`` and neither file alone shows that.
    """
    bindings = _destination_bindings(tree)
    loop_values: dict[str, list[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.For):
            # ``for frame, name in ((a, "training_recovery"), ...)`` -- the destination is
            # ``OUT / f"{name}.csv"``, so without the loop the pattern is ``*.csv`` and
            # claims every table in the repository. The literal names are right there.
            for variable, literals in _loop_literals(node).items():
                loop_values.setdefault(variable, []).extend(literals)

    patterns: list[str] = []
    if helpers is None:
        helpers = _write_helpers(tree)
    shadowed = _shadowed_names(tree)
    for argument in _write_destinations(tree, helpers):
        candidates = _bound_patterns(argument, bindings, shadowed)
        for candidate in candidates:
            # `store.write(f"checkpoints/{model_id}_start_{n}.json", ...)` names a
            # directory too; `_writers_of` matches on the filename alone, so only that
            # half of the pattern can be used.
            candidate = pathlib.PurePosixPath(candidate).as_posix()
            if "." not in candidate:
                continue
            for resolved in _expand_loop_names(argument, candidate, loop_values):
                if resolved not in patterns:
                    patterns.append(resolved)
    return patterns


def _loop_literals(node: ast.For) -> dict[str, list[str]]:
    """String constants a ``for`` binds to each loop variable, where they are literal.

    Handles the two shapes this repository writes: ``for name in ("a", "b")`` and
    ``for frame, name in ((x, "a"), (y, "b"))``. Anything computed yields nothing, and the
    caller then falls back to the wildcard.
    """
    items = node.iter
    if not isinstance(items, (ast.Tuple, ast.List)):
        return {}
    found: dict[str, list[str]] = {}
    if isinstance(node.target, ast.Name):
        values = [e.value for e in items.elts
                  if isinstance(e, ast.Constant) and isinstance(e.value, str)]
        if values:
            found[node.target.id] = values
        return found
    if not isinstance(node.target, (ast.Tuple, ast.List)):
        return {}
    names = [e.id if isinstance(e, ast.Name) else None for e in node.target.elts]
    for element in items.elts:
        if not isinstance(element, (ast.Tuple, ast.List)):
            continue
        for name, value in zip(names, element.elts):
            if name and isinstance(value, ast.Constant) and isinstance(value.value, str):
                found.setdefault(name, []).append(value.value)
    return found


def _expand_loop_names(argument: ast.AST, pattern: str,
                       loop_values: dict[str, list[str]]) -> list[str]:
    """Substitute a loop variable's literal values into ``pattern``, or keep the wildcard.

    ``OUT / f"{name}.csv"`` with ``name`` running over four literals is four filenames, not
    one glob over every CSV in the repository -- and the difference decides whether an
    unrelated table is reported as having a writer it does not have.
    """
    # ``OUT / f"{name}.csv"`` is a BinOp whose right side is the f-string, and
    # ``_pattern_from`` descends it the same way. Reading ``values`` off the BinOp would
    # find nothing and silently fall back to the wildcard.
    node = argument
    while isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        node = node.right
    fields = [v.value.id for v in getattr(node, "values", [])
              if isinstance(v, ast.FormattedValue) and isinstance(v.value, ast.Name)]
    if len(fields) != 1 or fields[0] not in loop_values:
        return [pattern]
    head, _, tail = pattern.partition("*")
    return [f"{head}{value}{tail}" for value in loop_values[fields[0]]]


def _too_weak_to_attribute(pattern: str) -> bool:
    """Whether a pattern claims a whole file type rather than a particular file.

    ``viz/figures.py`` writes ``out_dir / f"{stem}.svg"``, and the stem is the caller's.
    Resolving that binding yields ``*.svg``, which is true and useless: it matches every
    SVG in the repository, so every script reaching that module would claim every figure
    and the "nobody can regenerate this" signal would never fire again for a figure.

    Which script wrote which figure is not recoverable from the library at all -- but it is
    recorded, per file, in ``figures_manifest.csv``, and :func:`_writers_from_manifest`
    reads it. A bare extension is dropped here so that the recorded answer is the only one.
    """
    stem, _, extension = pathlib.PurePosixPath(pattern).name.rpartition(".")
    return bool(extension) and stem.strip("*") == ""


def _writers_from_manifest(repo: pathlib.Path) -> dict[str, list[str]]:
    """Figure -> producing script, from the manifest the figure run writes.

    ``make_figures.py`` records ``figure_produced_by`` beside every ``svg`` and ``png`` it
    emits. That is the authoritative statement of what produced what, written by the thing
    that produced it, and it is exactly the attribution a glob over the source cannot
    recover. Absent or malformed, this contributes nothing rather than raising: the audit
    runs on repositories that have never generated a figure.
    """
    manifest = repo / "outputs" / "figures_manifest.csv"
    if not manifest.exists():
        return {}
    try:
        table = pd.read_csv(manifest)
    except (pd.errors.EmptyDataError, ValueError):
        return {}
    if not {"figure_produced_by", "svg", "png"} <= set(table.columns):
        return {}
    found: dict[str, list[str]] = {}
    for row in table.itertuples():
        command = str(row.figure_produced_by)
        match = re.search(r"(scripts/[\w/]+\.py)", command)
        if not match:
            continue
        script = match.group(1)[len("scripts/"):]
        for column in ("svg", "png"):
            value = getattr(row, column, None)
            if isinstance(value, str) and value:
                name = pathlib.PurePosixPath(value).as_posix()
                found.setdefault(script, [])
                if name not in found[script]:
                    found[script].append(name)
    return {script: sorted(names) for script, names in found.items()}


def writers_from_scripts(repo: pathlib.Path) -> dict[str, list[str]]:
    """Map each script under ``scripts/`` to the output filenames it actually writes.

    Read from the syntax tree rather than from documentation, and rather than from a list
    kept here, so a script added tomorrow is covered without anyone remembering to add it.
    Only genuine write calls count: ``run_calibration.py`` names
    ``sensor_characterisation.csv`` too, but it reads it.

    Returns:
        ``{"run_power.py": ["power_analysis.csv"], ...}``, patterns possibly containing
        ``*`` where the name is computed.
    """
    found: dict[str, list[str]] = {}
    scripts = repo / "scripts"
    src = repo / "src" / "ystwin"
    if not scripts.is_dir():
        return found
    for path in sorted(scripts.rglob("*.py")):
        # A script's own body is not the whole story. `make_figures.py` writes nothing
        # itself: every figure leaves through `fig.savefig` inside `ystwin.viz.figures`,
        # so scanning only `scripts/` attributed twelve committed files to nobody. The
        # reachability walk that the freshness check already runs is the right authority
        # here too -- if a script reaches a module, that module's writes are the script's.
        reachable = [path, *sorted(library_dependencies(path, src))]
        trees = []
        for source in reachable:
            try:
                trees.append(ast.parse(source.read_text(encoding="utf-8")))
            except (SyntaxError, OSError):
                continue
        # Helpers are pooled over everything the script reaches, for the same reason the
        # write scan is: `native_population_scoring.py` writes every one of its artefacts
        # through the `ArtifactStore` defined in `native_population_development.py`.
        helpers: dict[str, int] = {}
        for tree in trees:
            helpers.update(_write_helpers(tree))
        patterns = []
        for tree in trees:
            for pattern in _patterns_in(tree, helpers):
                if pattern not in patterns and not _too_weak_to_attribute(pattern):
                    patterns.append(pattern)
        if patterns:
            found[path.relative_to(scripts).as_posix()] = sorted(patterns)
    for script, names in _writers_from_manifest(repo).items():
        merged = sorted(set(found.get(script, [])) | set(names))
        found[script] = merged
    # The provenance stamp is written by THIS file, through `write_text` rather than any
    # of the calls the scan above recognises, so it is declared here. Added AFTER the
    # loop and not before it: seeded first, the loop's own pass over
    # `audit_reproducibility.py` overwrote it, and the stamp reported as a result
    # nobody can regenerate -- a false instance of exactly what this map exists to find.
    # Only where that script actually exists. Added unconditionally it appeared in
    # every synthetic tree a test builds, and a map that claims a writer for a file
    # in a repository holding neither is worse than the gap it was closing.
    return found

def writers_from_doc(text: str) -> dict[str, list[str]]:
    """The same map as the reproduction guide states it, for cross-checking.

    The guide's table is prose and drifts like prose. Reading it here means the audit can
    say *which* of the two is wrong instead of only that something is.
    """
    found: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2:
            continue
        script = re.fullmatch(r"`([\w/]+\.py)`", cells[0])
        if not script:
            continue
        # The alternation is the artifact extension set the table's scripts actually write
        # (measured); without `json` no writer row could ever name `gates_manifest.json`.
        names = re.findall(r"`([\w<>.\-]+\.(?:csv|json|npz|png|xlsx))`", cells[1])
        if names:
            found[script.group(1)] = sorted(re.sub(r"<[^>]+>", "*", n) for n in names)
    return found


def _patterns_agree(one: str, other: str) -> bool:
    """Whether two patterns could name the same file.

    The guide is free to be more specific than the source can be: ``run_d2.py`` builds its
    filename from a tag that happens to contain the reporter, so the syntax tree can only
    say ``d2_*.csv`` where the guide says ``d2_<plate>__<reporter>.csv``. Those are the
    same set of files described at two resolutions and flagging it would be a false alarm.
    Matching each pattern against the other as a glob accepts exactly that case and still
    rejects a name neither side mentions.
    """
    one = one if "/" in one else f"outputs/{one}"
    other = other if "/" in other else f"outputs/{other}"
    return _path_pattern_matches(one, other) or _path_pattern_matches(other, one)


def check_writer_map_matches_doc(from_scripts: dict[str, list[str]],
                                 from_doc: dict[str, list[str]],
                                 doc_name: str = "docs/REPRODUCING.md") -> list[dict]:
    """The documented script-to-output map must match the one in the code.

    Only scripts the guide lists are compared. A script the guide omits is a documentation
    gap that ``audit_claims.py`` is the right place to catch; here the question is narrower
    -- does the table that tells a stranger what to run tell them the truth?
    """
    rows = []
    for script in sorted(from_doc):
        documented = set(from_doc[script])
        actual = set(from_scripts.get(script, []))
        if not actual:
            rows.append(_row(f"writer map: {script}", sorted(documented), "no write calls found",
                             False, f"{doc_name} says it writes these; the source does not"))
            continue
        missing = {d for d in documented if not any(_patterns_agree(d, a) for a in actual)}
        extra = {a for a in actual if not any(_patterns_agree(a, d) for d in documented)}
        ok = not missing and not extra
        detail = ""
        if missing:
            detail += f"documented but not written: {sorted(missing)}. "
        if extra:
            detail += f"written but not documented: {sorted(extra)}."
        rows.append(_row(f"writer map: {script}", sorted(documented), sorted(actual), ok,
                         detail.strip()))
    return rows


def _module_path(dotted: str, src: pathlib.Path) -> pathlib.Path | None:
    """``ystwin.analysis.power`` -> the file that defines it, if it is in this package."""
    if not dotted.startswith("ystwin"):
        return None
    parts = dotted.split(".")[1:]
    base = src.joinpath(*parts)
    for candidate in (base.with_suffix(".py"), base / "__init__.py"):
        if candidate.exists():
            return candidate
    return None


def _dotted_of(path: pathlib.Path, src: pathlib.Path) -> str:
    relative = path.relative_to(src).with_suffix("")
    parts = [p for p in relative.parts if p != "__init__"]
    return ".".join(["ystwin", *parts])


def _imports_in(path: pathlib.Path, src: pathlib.Path) -> set[str]:
    """Dotted ``ystwin`` module names one file imports, relative imports resolved.

    ``importlib.import_module("ystwin.analysis.native_population_development")`` counts.
    Three scripts load their library that way and an import-statement-only walk called
    them dependency-free, so every file they write -- the forty-one training checkpoints
    among them -- was attributed to no script at all. A string literal handed to
    ``import_module`` states the dependency exactly as an ``import`` line does; a computed
    one (``f"ystwin.{record.code_module}"``) states nothing this can read, and is skipped.
    """
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return set()
    dotted = _dotted_of(path, src) if path.is_relative_to(src) else ""
    package = dotted.rsplit(".", 1)[0] if (dotted and path.name != "__init__.py") else dotted
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package
                for _ in range(node.level - 1):
                    base = base.rsplit(".", 1)[0] if "." in base else base
                root = f"{base}.{node.module}" if node.module else base
            else:
                root = node.module or ""
            if not root.startswith("ystwin"):
                continue
            names.add(root)
            names.update(f"{root}.{a.name}" for a in node.names)
        elif isinstance(node, ast.Call) and node.args:
            called = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(
                node.func, "id", None)
            argument = node.args[0]
            if called == "import_module" and isinstance(argument, ast.Constant) and isinstance(
                    argument.value, str):
                names.add(argument.value)
    return {n for n in names if n.startswith("ystwin")}


def library_dependencies(script: pathlib.Path, src: pathlib.Path) -> set[pathlib.Path]:
    """Every file under ``src/ystwin`` a script reaches, transitively.

    A script is usually thin -- it parses arguments and calls the library -- so its own
    modification time says almost nothing about whether its output is current. What
    matters is the code the numbers came out of. Import statements are the authority on
    what that is, including the ones inside functions, which this package uses to defer
    the expensive imports.
    """
    if not src.is_dir():
        return set()
    seen: set[pathlib.Path] = set()
    frontier = [script]
    while frontier:
        current = frontier.pop()
        for dotted in _imports_in(current, src):
            module = _module_path(dotted, src)
            if module is not None and module not in seen:
                seen.add(module)
                frontier.append(module)
    return seen


def outputs_read_by(source: pathlib.Path) -> set[str]:
    """Filenames under ``outputs/`` that this script READS, by the same AST walk as writes.

    A name is taken as read when it is the first argument of a :data:`_READ_CALLS` call.
    Names the script also writes are excluded by the caller: a script that rewrites a table
    it read is not stale against itself.
    """
    try:
        tree = ast.parse(source.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return set()
    bindings = _destination_bindings(tree)
    names: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        attr = node.func.attr if isinstance(node.func, ast.Attribute) else None
        argument = None
        if attr in _READ_CALLS and node.args:
            argument = node.args[0]
        elif attr in {"read_text", "read_bytes"}:
            argument = node.func.value
        elif attr == "open" and not _writing_mode(node, 0):
            argument = node.func.value
        elif isinstance(node.func, ast.Name) and node.func.id == "open" and node.args:
            if not _writing_mode(node, 1):
                argument = node.args[0]
        if argument is not None:
            names.update(pattern for pattern in _bound_patterns(argument, bindings)
                         if pattern.startswith("outputs/") and "*" not in pattern)
    return names


def _path_pattern_matches(identity: str, pattern: str) -> bool:
    """Whether ``pattern`` names ``identity``, anchored at the RIGHT.

    A pattern carries as much of the path as the source spelled out and no more:
    ``store.write(f"checkpoints/{model_id}_start_{n}.json", ...)`` says two components and
    ``frame.to_csv(OUT / "power_analysis.csv")`` says one, because the run directory is a
    command-line argument in both cases. Requiring equal component counts instead made
    every one-component pattern unmatchable against a file inside a run directory, which
    is most of ``outputs/`` -- 575 files reported as having no writer at all, the same
    false alarm in a new place. Matching from the right keeps the extra precision a
    directory-qualified pattern buys and still lets a bare filename find its file.
    """
    path_parts, pattern_parts = identity.split("/"), pattern.split("/")
    if len(pattern_parts) > len(path_parts):
        return False
    return all(fnmatch.fnmatchcase(part, rule)
               for part, rule in zip(path_parts[-len(pattern_parts):], pattern_parts))


def _writers_of(name: str, from_scripts: dict[str, list[str]]) -> list[str]:
    """The script that writes ``name``, resolving overlapping patterns by specificity.

    ``run_gates.py`` writes ``d2_g1passed_<plate>.csv`` and ``run_d2.py`` writes
    ``d2_<plate>.csv``; the second glob also matches the first's files. Attributing a table
    to both would report the same staleness twice and blame the wrong script once, so the
    longest literal match wins, which is the usual rule for overlapping globs.
    """
    matches = [(s, p) for s, patterns in from_scripts.items() for p in patterns
               if _path_pattern_matches(name, p)]
    if not matches:
        return []
    best = max(len(p.replace("*", "")) for _, p in matches)
    return sorted({s for s, p in matches if len(p.replace("*", "")) == best})


_RESOLVER_CALL = re.compile(r"paths\.([a-z_][a-z_0-9]*)\(\)")

# Resolvers that answer with a directory inside the checkout, or that take arguments. They
# are never the reason a script cannot run, so counting them would make every script look
# unrunnable.
_NOT_AN_ASSET = frozenset({"outputs_dir", "data_dir", "pathways_dir", "require",
                           "resolve_or_exit"})


def absent_assets(source: pathlib.Path, paths_module) -> list[str]:
    """Asset resolvers this script calls that return ``None`` on this machine.

    Read off the script's own ``paths.<resolver>()`` calls rather than from a hand-kept
    table of script-to-asset, because a hand-kept table is a second place to forget: a
    script that grows a new data dependency would keep its old, wrong entry and the answer
    below would be confidently stale. The regex is doing what an import does -- naming the
    resolver the script actually calls -- and then the resolver itself is asked.

    Used by :func:`check_output_freshness` to tell two facts apart that it previously
    reported as one.
    """
    text = source.read_text()
    absent = []
    for name in sorted(set(_RESOLVER_CALL.findall(text))):
        if name in _NOT_AN_ASSET:
            continue
        resolver = getattr(paths_module, name, None)
        if resolver is None:
            continue
        try:
            if resolver() is None:
                absent.append(name)
        except TypeError:  # takes arguments, so it is not a bare asset resolver
            continue
    return absent


#: Where the producer digests live. Tracked, because the whole point is that it survives a
#: clone -- an mtime does not, and that is the hole this closes.
PROVENANCE = "outputs/provenance.json"


def producer_digest(producers: list[pathlib.Path], repo: pathlib.Path = REPO) -> str:
    """A digest of everything that produced a table: its script, its library modules, and
    the tables it read.

    Content, not timestamps. `check_output_freshness` compares mtimes and its own docstring
    says why -- "a fresh clone writes every file at once, so mtimes there are equal and
    nothing fires". That was written as a reassurance and is the vulnerability: `git
    checkout` stamps every file with the checkout time, so on any clone the mtime check
    passes **vacuously**, on content that may not match the code at all.

    It was not hypothetical. On 2026-08-31 the mtime gate reported 284/284 while eight
    committed artefacts did not match what their own scripts produced -- `heldout_scores.csv`
    was 179 training rows against the code's 287, and five biosensor tables moved enough to
    break ten pinned prose numbers. Every one had been stale for weeks behind a green gate.

    A digest fires in exactly the case an mtime cannot: the files are all the same age and
    the content disagrees anyway.
    """
    return artifact_producer_digest(repo, producers)


def check_output_provenance(repo: pathlib.Path, from_scripts: dict[str, list[str]],
                            restamp: bool = False) -> list[dict]:
    """Validate exact artifact identities and reproduction evidence without writing.

    The syntax-tree map is discovery support, not provenance authority. Legacy basename
    digests cannot establish inputs, parameters, runtime, or scientific equivalence.
    A changed identity requires investigation, not an assertion that numbers are wrong.
    ``restamp`` is retained only to refuse unsafe callers explicitly.
    """
    if restamp:
        raise ArtifactError("blind restamping is disabled; reproduce in an isolated directory and compare values, missingness, units and row keys")
    return audit_registry(repo)


#: Files under ``outputs/`` that no script writes, and the reason that is correct.
#:
#: Every entry is an fnmatch pattern relative to ``outputs/`` and every entry carries a
#: reason, because an allowlist without one is indistinguishable from a gap somebody got
#: tired of seeing. :func:`check_writerless_allowlist` fails an entry that matches nothing
#: and an entry whose files have since acquired a writer, so the list cannot quietly
#: outlive what it excuses.
#:
#: Only tracked paths belong here. A clone holds exactly the tracked files, so an entry
#: naming a local run's leftovers would fail on every machine but the one that made it.
WRITERLESS_BY_DESIGN: dict[str, str] = {}


def _writerless_reason(relative: str) -> str | None:
    """Legacy callers receive no filename-based exemption; the registry owns lifecycle."""
    return None


def check_writerless_allowlist(repo: pathlib.Path,
                               from_scripts: dict[str, list[str]]) -> list[dict]:
    """Compatibility entrypoint for callers migrating to exact artifact registry checks."""
    return audit_registry(repo)


def check_output_freshness(repo: pathlib.Path, from_scripts: dict[str, list[str]],
                           paths_module=None,
                           tolerance_s: float = FRESHNESS_TOLERANCE_S) -> list[dict]:
    """No table in ``outputs/`` may predate the code that produced it.

    A table older than its writer documents behaviour the code no longer has, and it does
    so in the most convincing possible form -- committed, tabular, citable. This is the
    check that turns "re-run everything before the deadline" from a memory into a gate.

    "The code" means the script *and* the library modules it reaches, because the scripts
    are thin and the behaviour lives under ``src/ystwin``. Which module is newest goes in
    the detail column, so the report says what to re-run and why rather than only that
    something moved. A table cannot be blamed for a change in a module its script never
    imports, which is why the dependency set is resolved per script rather than taken as
    all of ``src``.

    Modification times rather than commit times: a fresh clone writes every file at once,
    so mtimes there are equal and nothing fires, while on the machine where a table was
    actually produced the mtimes are the real record of what was run when.

    A file no script writes gets its own failure. That is not pedantry -- ``outputs/`` here
    contains one table with no reproducer at all, which the guide admits is a gap.
    """
    if (repo / REGISTRY_PATH).is_file():
        return [_row("artifact freshness policy", "content-bound reproduction", "mtime is not evidence of numerical disagreement", "SKIP",
                     "The registry provenance check validates current reproduction receipts and original frozen identities; source timestamps cannot establish values.")]
    outputs = repo / "outputs"
    if not outputs.is_dir():
        return [_row("outputs/ exists", "a directory", "missing", False)]
    src = repo / "src" / "ystwin"
    rows = []
    cache: dict[str, set[pathlib.Path]] = {}
    for path in sorted(p for p in outputs.rglob("*") if p.is_file()):
        relative = path.relative_to(repo).as_posix()
        writers = _writers_of(relative, from_scripts)
        if not writers:
            rows.append(_row(f"output has a writer: {relative}",
                             "an exact artifact registry entry", "investigation_pending", False,
                             "No producer established for this full path; filename discovery is not proof that a reproducer exists or is absent."))
            continue
        produced = path.stat().st_mtime
        for script in writers:
            source = repo / "scripts" / script
            if not source.exists():
                continue
            if script not in cache:
                cache[script] = library_dependencies(source, src)
            # The tables this script READS are producers of the table it writes, exactly as
            # its library modules are. Without this edge a table rebuilt from a superseded
            # table passes: see _READ_CALLS.
            input_names = set().union(*(outputs_read_by(module) for module in [source, *cache[script]]))
            consumed = [repo / name for name in sorted(input_names)
                        if name != relative and (repo / name).is_file()]
            producers = [source, *sorted(cache[script]), *consumed]
            newest = max(producers, key=lambda p: p.stat().st_mtime)
            written = newest.stat().st_mtime
            fresh = produced >= written - tolerance_s
            age_h = (written - produced) / 3600.0
            # STALE AND REBUILDABLE IS NOT THE SAME FACT AS STALE AND UNREBUILDABLE, and
            # this check reported them identically until now. Six of its eight failures
            # were `run_g4.py` and `make_splits.py` tables, both of which need the plate
            # exports that are not ours to redistribute -- and both scripts REFUSE to
            # rebuild from a partial input, deliberately, because a table rebuilt from one
            # of three replicates "would look like a finding about the data rather than
            # about this machine". So the gate was demanding an action it had also, rightly,
            # made impossible, and it did so on every machine on earth except the one
            # holding the raw exports.
            #
            # A SKIP says the check could not run here. It is not a pass: `_row`'s own
            # doctrine is that "a check that never ran has proved nothing", the summary
            # counts skips separately, and the table still shows the staleness in hours.
            # What it stops doing is failing a clone for not having data the repository
            # tells it that it cannot have.
            missing = absent_assets(source, paths_module) if paths_module is not None else []
            status: object = fresh
            actual = "current" if fresh else f"{age_h:.1f} h stale"
            detail = ("" if fresh else
                      f"{newest.relative_to(repo).as_posix()} changed after this table was "
                      f"written; re-run scripts/{script} or the table describes behaviour "
                      "the code no longer has")
            if not fresh and missing:
                status = "SKIP"
                detail = (
                    f"{age_h:.1f} h stale and NOT REBUILDABLE HERE: scripts/{script} needs "
                    f"{missing}, which resolve to nothing on this machine. "
                    f"{newest.relative_to(repo).as_posix()} changed after the table was "
                    "written, so it may describe behaviour the code no longer has -- but "
                    "re-running here would replace it with a table built from partial "
                    "inputs, which the script itself refuses to do. Supply the data, or "
                    "treat the table as the record of the machine that had it")
            rows.append(_row(
                f"output not older than its code: {relative}",
                f"at least as new as scripts/{script} and the {len(producers) - 1} "
                f"ystwin modules it reaches"
                + (f" and the {len(consumed)} outputs/ table(s) it reads" if consumed else ""),
                actual, status, detail))
    return rows


# ---------------------------------------------------------------------------
# does it work for someone with no data at all
# ---------------------------------------------------------------------------

_SUMMARY = re.compile(r"(\d+)\s+(passed|failed|skipped|errors?|xfailed|xpassed|deselected)")


def _pytest_counts(output: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for number, word in _SUMMARY.findall(output):
        key = "error" if word.startswith("error") else word
        counts[key] = counts.get(key, 0) + int(number)
    return counts


def check_offline_suite(repo: pathlib.Path, env_vars: list[str], target: str = "tests",
                        timeout_s: float = 3600.0) -> list[dict]:
    """The suite must pass with every data-locating variable pointed at nothing.

    This is the check that proves the repository works for someone who has no plate
    reader. The claim is made in ``docs/REPRODUCING.md`` and is easy to break without
    noticing: a fixture that reads a file directly instead of going through
    :mod:`ystwin.paths`, or a module-level constant evaluated at import, turns "you do not
    have the wet-lab exports" into "the tests are broken" for every outside reader.

    Pointing the variables at a directory that does not exist -- rather than unsetting them
    -- is what makes the check meaningful, because ``paths._resolve`` treats a variable
    that is set as authoritative and does *not* fall through to the sibling-directory
    default. Unsetting them would find the author's data again.

    Both runs report their skip count, since the difference between them is the size of
    the real-data surface.

    Args:
        repo: Checkout root.
        env_vars: Variables to point at nothing.
        target: What to hand pytest. Narrow it to a file or two for a smoke test.
        timeout_s: Give up after this long.
    """
    nowhere = f"/nonexistent-ystwin-audit-{uuid.uuid4().hex}"
    # No -q here on purpose. This repository's pytest configuration already passes it, and
    # a second one makes -qq, which suppresses the summary line this function parses.
    command = [sys.executable, "-m", "pytest", target, "--tb=no", "-p", "no:cacheprovider"]

    def run(env: dict) -> tuple[dict[str, int], str]:
        try:
            done = subprocess.run(command, cwd=repo, capture_output=True, text=True,
                                  env=env, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            return {}, "timed out"
        return _pytest_counts(done.stdout + done.stderr), done.stdout[-2000:]

    baseline, _ = run(dict(os.environ))
    offline_env = dict(os.environ)
    offline_env.update({name: nowhere for name in env_vars})
    offline_env["YSTWIN_OUTPUTS"] = f"{nowhere}-outputs"
    offline, tail = run(offline_env)

    if not offline:
        return [_row("suite passes with no real data", "0 failed, 0 errors",
                     "could not run pytest", False, tail[-300:])]

    broken = offline.get("failed", 0) + offline.get("error", 0)
    rows = [_row(
        "suite passes with no real data", "0 failed, 0 errors",
        f"{offline.get('failed', 0)} failed, {offline.get('error', 0)} errors",
        broken == 0,
        "" if broken == 0 else
        "a machine with no wet-lab data cannot run this repository; the failures above "
        "are what an outside reader sees first")]
    rows.append(_row(
        "skips with data present", "reported",
        f"{baseline.get('skipped', 0)} skipped of "
        f"{sum(baseline.get(k, 0) for k in ('passed', 'skipped', 'failed'))}",
        "SKIP", "informational: the real-data surface on this machine"))
    rows.append(_row(
        "skips with data absent", "at least as many as with data present",
        f"{offline.get('skipped', 0)} skipped of "
        f"{sum(offline.get(k, 0) for k in ('passed', 'skipped', 'failed'))}",
        offline.get("skipped", 0) >= baseline.get("skipped", 0),
        "" if offline.get("skipped", 0) >= baseline.get("skipped", 0) else
        "fewer skips without the data than with it, which cannot be right"))
    return rows


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def collect(repo: pathlib.Path, offline_suite: bool = False,
            suite_target: str = "tests", restamp: bool = False) -> list[dict]:
    """Run every check and return the rows, so a caller can report them how it likes."""
    from ystwin import paths

    pyproject = (repo / "pyproject.toml").read_text(encoding="utf-8")
    reproducing = repo / "docs" / "REPRODUCING.md"
    doc_text = reproducing.read_text(encoding="utf-8") if reproducing.exists() else ""

    records = describe_resolvers(paths)
    from_scripts = writers_from_scripts(repo)

    rows: list[dict] = []
    rows += check_asset_resolvers(records, repo)
    rows += check_env_vars_documented(records, doc_text)
    rows += check_python_version(pyproject)
    rows += check_dependency_pins(pyproject, _installed_version)
    rows += check_solver(*_solver())
    documented_writers = dict(from_scripts)
    if (repo / REGISTRY_PATH).is_file():
        try:
            documented_writers.update(_artifacts.registered_writers(repo))
        except (OSError, ValueError, KeyError, TypeError):
            pass
    rows += check_writer_map_matches_doc(documented_writers, writers_from_doc(doc_text))
    rows += check_outputs_tracked(repo)
    rows += check_output_tables(repo)
    rows += check_output_freshness(repo, from_scripts, paths_module=paths)
    rows += check_output_provenance(repo, from_scripts, restamp=restamp)
    if offline_suite:
        env_vars = [r["env_var"] for r in records if r["env_var"]]
        rows += check_offline_suite(repo, env_vars, target=suite_target)
    else:
        rows.append(_row("suite passes with no real data", "0 failed, 0 errors",
                         "not run", "SKIP",
                         "costs a full suite run twice; pass --offline-suite"))
    return rows


def _installed_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _solver() -> tuple[str | None, list[str] | None]:
    try:
        import cobra
        from cobra.util.solver import solvers
    except Exception:
        return None, None
    try:
        return cobra.Configuration().solver.__name__, sorted(solvers)
    except Exception:
        return None, sorted(solvers)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Read-only reproducibility audit using exact artifact identities and verified lineage.")
    parser.add_argument("--quiet", action="store_true", help="only print failures")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="explicit temporary directory outside the repository for CSV and JSON reports")
    parser.add_argument("--restamp", action="store_true", help="disabled: blind restamping is not evidence of reproduction")
    parser.add_argument("--offline-suite", action="store_true",
                        help="also run the suite twice, once with the data hidden (slow)")
    parser.add_argument("--suite-target", default="tests",
                        help="what to hand pytest for --offline-suite; narrow it to smoke-test")
    args = parser.parse_args(argv)
    if args.restamp:
        parser.error("--restamp is disabled; use scripts/regenerate_current_artifacts.py with --run and a temporary --output-dir to reproduce and compare before any reviewed provenance update")
    out = args.output_dir.resolve() if args.output_dir is not None else None
    if out is not None and out.is_relative_to(REPO.resolve()):
        parser.error("--output-dir must be a temporary directory outside the repository")
    if out is not None and any((out / f"audit_reproducibility.{suffix}").exists() for suffix in ("csv", "json")):
        parser.error("audit reports are write-once; choose a fresh --output-dir")
    rows = collect(REPO, offline_suite=args.offline_suite,
                   suite_target=args.suite_target)
    frame = pd.DataFrame(rows, columns=["check", "expected", "actual", "status", "detail"])

    failures = frame[frame.status == "FAIL"]
    skips = frame[frame.status == "SKIP"]
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)
        with (out / "audit_reproducibility.csv").open("x", encoding="utf-8", newline="") as stream:
            frame.to_csv(stream, index=False)
        with (out / "audit_reproducibility.json").open("x", encoding="utf-8") as stream:
            json.dump({"schema_version": 1, "read_only": True, "rows": rows,
                       "summary": {"passed": len(frame) - len(failures) - len(skips), "failed": len(failures), "skipped": len(skips)},
                       "investigation_queue": [row for row in rows if row["status"] == "FAIL"]},
                      stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    width = 72
    print("=" * width)
    print("Reproducibility audit: could a stranger get these numbers?")
    print("=" * width)
    for _, r in frame.iterrows():
        if args.quiet and r.status == "PASS":
            continue
        mark = {"PASS": "ok  ", "SKIP": "skip"}.get(r.status, "FAIL")
        print(f"  [{mark}] {r.check}")
        if r.status != "PASS":
            print(f"         expected {r.expected}, got {r.actual}")
            if r.detail:
                print(f"         {r.detail}")

    print()
    report = str(out / "audit_reproducibility.json") if out is not None else "read-only; no report files written"
    print(f"  {len(frame) - len(failures) - len(skips)}/{len(frame) - len(skips)} checks "
          f"pass ({len(skips)} skipped) -> {report}")
    if len(failures):
        print(f"  {len(failures)} FAILING. These results are not reproducible as they stand.")
    return 1 if len(failures) else 0


if __name__ == "__main__":
    raise SystemExit(main())
