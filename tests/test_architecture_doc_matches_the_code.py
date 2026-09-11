"""`docs/ARCHITECTURE.md` must name the modules that exist, and all of them.

A 1,500-line hand-maintained map of a 88-module package drifts. It had drifted: it named two
modules that had just been archived, and it had **zero** mentions of the `pathway/` layer
that had become the thing the whole product prediction runs through.

The failure mode is not that the prose is wrong in an interesting way. It is that a reader
uses the map to decide what to read, and a map missing a layer sends them somewhere else.
That is worth a test rather than a habit, for the same reason
`scripts/audit_claims.py` pins numbers to the cells they came from: prose drifts and nobody
notices, and the noticing is the expensive part.

Two directions, and the second is the one that catches a new layer going unmentioned:

    every module the doc names          ->  must exist
    every module in src/ystwin          ->  must be named somewhere in the doc

This does NOT check that what the doc says about a module is true, which is not mechanically
checkable and is what review is for. It checks that the inventory is current, which is the
part that goes wrong silently.
"""

from __future__ import annotations

import pathlib
import re

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_DOC = _REPO / "docs" / "ARCHITECTURE.md"
_SRC = _REPO / "src" / "ystwin"

# Package-relative paths the doc is allowed to omit. Each needs a reason, because an
# exemption list is how a check like this quietly stops checking anything.
_EXEMPT = {
    "__init__.py",          # namespace, not a layer
}

# Paths the doc names ON PURPOSE while saying they are absent. Naming a file in order to
# record that it does not exist is good documentation, not drift, and a check that could
# not tell the two apart would push an author towards deleting the warning instead of the
# staleness. Each entry needs a reason.
_DELIBERATELY_ABSENT = {
    # §6 records that README once listed this and no such file was ever written.
    "generator/augment.py",
    # §7 item 10 records that the parked script's own usage line gives the pre-move path.
    "scripts/run_d1.py",
}


def _modules_on_disk() -> set[str]:
    return {str(p.relative_to(_SRC)) for p in _SRC.rglob("*.py")
            if p.name not in _EXEMPT}


def _modules_the_doc_names() -> set[str]:
    """Every `path/to/module.py` the doc mentions, in prose, tables or mermaid blocks."""
    text = _DOC.read_text(encoding="utf-8")
    return set(re.findall(r"\b((?:[a-z_][a-z0-9_]*/)*[a-z_][a-z0-9_]*\.py)\b", text))


@pytest.mark.skipif(not _DOC.exists(), reason="no ARCHITECTURE.md in this checkout")
class TestTheInventoryIsCurrent:
    def test_the_doc_names_no_module_that_has_been_removed(self):
        """The direction that fires on an archive or a rename. It caught
        `analysis/deconvolve.py` and `gates/module_admission.py` still being drawn into two
        dependency diagrams after both had been archived."""
        named = _modules_the_doc_names()
        on_disk = _modules_on_disk()

        ghosts = sorted(
            name for name in named
            if name not in on_disk
            and name not in _DELIBERATELY_ABSENT
            and not (_REPO / "scripts" / name).exists()
            and not (_REPO / "tests" / name).exists()
            and not (_REPO / name).exists()
            and "/" in name)                      # bare filenames are usually prose

        assert ghosts == []

    def test_every_module_in_the_package_is_named_somewhere(self):
        """The direction that fires on a NEW layer nobody wrote up -- which is how the
        product architecture came to be documented as three layers while running five."""
        named = {n.split("/")[-1] for n in _modules_the_doc_names()}
        missing = sorted(m for m in _modules_on_disk()
                         if m.split("/")[-1] not in named)

        assert missing == [], (
            f"{len(missing)} module(s) exist and the architecture map does not mention "
            f"them: {missing}")


@pytest.mark.skipif(not _DOC.exists(), reason="no ARCHITECTURE.md in this checkout")
class TestTheCheckCanActuallyFail:
    """An inventory check that cannot fire is worse than none, because it reads as
    assurance."""

    def test_it_finds_the_modules_at_all(self):
        assert len(_modules_the_doc_names()) > 30

    def test_a_removed_module_would_be_caught(self):
        """The archived pair, which is what the check found on its first run."""
        on_disk = _modules_on_disk()

        assert "analysis/deconvolve.py" not in on_disk
        assert "gates/module_admission.py" not in on_disk

    def test_every_deliberate_absence_is_still_absent(self):
        """An exemption list is how a check like this quietly stops checking. If one of
        these files is ever written, the entry becomes a lie and should go."""
        on_disk = _modules_on_disk()

        for name in _DELIBERATELY_ABSENT:
            assert name not in on_disk, f"{name} now exists; drop its exemption"
