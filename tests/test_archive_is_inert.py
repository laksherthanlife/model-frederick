"""Nothing in `archive/` may be imported, and nothing outside it may name what it holds.

An archive that the live code still reaches into is not an archive, it is a subdirectory
with a discouraging name. The property that makes archiving mean something is that the main
tree runs without it -- so this file checks the boundary rather than trusting it.

The rule for what belongs there is in `archive/README.md` and is narrower than "unused":
superseded, or orphaned with no plan. Parked code is not archived, and code that is only
reachable from its own tests is not archived either, because being called by nothing *yet*
is different from being replaced.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_ARCHIVE = _REPO / "archive"
_LIVE = ("src", "scripts", "tests")


def _archived_module_names() -> set[str]:
    return {p.stem for p in (_ARCHIVE / "code").glob("*.py")}


@pytest.mark.skipif(not _ARCHIVE.is_dir(), reason="no archive/ in this checkout")
class TestTheBoundaryHolds:
    def test_no_live_module_imports_anything_archived(self):
        archived = _archived_module_names()
        offenders = []
        for directory in _LIVE:
            for path in (_REPO / directory).rglob("*.py"):
                try:
                    tree = ast.parse(path.read_text(encoding="utf-8"))
                except SyntaxError:
                    continue
                for node in ast.walk(tree):
                    names: list[str] = []
                    if isinstance(node, ast.ImportFrom):
                        names = [(node.module or "").split(".")[-1]]
                        names += [a.name for a in node.names]
                    elif isinstance(node, ast.Import):
                        names = [a.name.split(".")[-1] for a in node.names]
                    for name in names:
                        if name in archived:
                            offenders.append(
                                f"{path.relative_to(_REPO)}:{node.lineno} imports {name!r}")

        assert offenders == []

    def test_archived_code_imports_nothing_from_the_live_tree_either(self):
        """The other direction, and the one that decides whether archived code can be read
        as a historical record. A module that still imports today's `ystwin` documents
        today's package, not the one it was written against, and will break the moment the
        live signature it depends on changes.

        `archive/tests/` is deliberately exempt. A test is inseparable from the thing it
        tested -- rewriting its imports to reach nothing would leave a file that documents
        neither -- and nothing collects them, so they cannot break the suite.
        """
        offenders = []
        for path in (_ARCHIVE / "code").rglob("*.py"):
            try:
                tree = ast.parse(path.read_text(encoding="utf-8"))
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                module = ""
                if isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                elif isinstance(node, ast.Import):
                    module = ",".join(a.name for a in node.names)
                if "ystwin" in module:
                    offenders.append(
                        f"{path.relative_to(_REPO)}:{node.lineno} imports {module!r}")

        assert offenders == []

    def test_pytest_does_not_collect_the_archive(self):
        """Archived tests exercise archived code. Collecting them would make the suite
        depend on the thing the archive exists to detach from."""
        import tomllib

        config = tomllib.loads((_REPO / "pyproject.toml").read_text())
        testpaths = config.get("tool", {}).get("pytest", {}).get("ini_options", {}).get(
            "testpaths", [])

        assert "archive" not in testpaths

    def test_the_archive_says_why_each_thing_is_there(self):
        """An archive without reasons is a graveyard, and a reader who cannot tell why
        something was retired cannot tell whether it should come back."""
        readme = (_ARCHIVE / "README.md").read_text()

        for module in _archived_module_names():
            assert module in readme, f"{module} is archived and unexplained"
