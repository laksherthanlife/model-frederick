"""No call site may read ``slim_optimize()`` without handling nan.

`fba/solver.py` explains the trap and `growth_or_none` is the fix, but a fix in one module
does not hold a codebase: after three sites in `fba/` were corrected, **seven more were
still live** in `bridge/` and `scripts/`, and one of them was observed this session printing
``max growth nan /h  (  nan% of max)`` from `scripts/parked/run_d1.py`.

Every one of them failed in the same direction -- OPEN. An infeasible LP returns nan, and:

    nan or 0.0        is nan     a fallback that never falls back
    nan <= 0.0        is False   a "cannot grow" refusal that never refuses
    nan < threshold   is False   a bound check that always passes
    nan is None       is False   an "infeasible" branch that is dead code
    nan * 0.9         is nan     a bound that FVA accepts and answers around

So the failure mode is never a crash. It is a confident number computed from a solve that
did not happen, which is the one outcome this repository is built to prevent. The worst was
`bridge/regulation.py`'s `InfeasibleRegulation`, whose entire purpose is refusing a
regulation layer that kills growth and which could not fire.

This test is a grep with a rule, deliberately: the property is syntactic, it holds over
files nothing imports (parked scripts included), and it costs no solver.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]

# The one module allowed to call it raw -- it is where the nan check lives.
EXEMPT = {REPO / "src" / "ystwin" / "fba" / "solver.py"}


def _python_files() -> list[pathlib.Path]:
    roots = (REPO / "src", REPO / "scripts")
    return sorted(p for root in roots for p in root.rglob("*.py")
                  if p not in EXEMPT and "__pycache__" not in p.parts)


def _guarded_by_isnan(tree: ast.Module, node: ast.Call) -> bool:
    """True when an ``isnan`` appears within a few lines of the call, in the same scope.

    Crude on purpose. The point is not to prove a guard correct -- it is to make an
    UNGUARDED read impossible to add without noticing, and the exact reads that were wrong
    all had no ``isnan`` anywhere near them.
    """
    line = node.lineno
    return any(
        isinstance(other, ast.Attribute) and other.attr == "isnan"
        and abs(getattr(other, "lineno", 0) - line) <= 4
        for other in ast.walk(tree))


def _unguarded_calls(path: pathlib.Path) -> list[int]:
    tree = ast.parse(path.read_text())
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        function = node.func
        if isinstance(function, ast.Attribute) and function.attr == "slim_optimize":
            if not _guarded_by_isnan(tree, node):
                out.append(node.lineno)
    return out


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: str(p.relative_to(REPO)))
def test_no_module_reads_slim_optimize_without_an_isnan_check(path):
    """Use ``growth_or_none`` from ``ystwin.fba.solver``, or check ``math.isnan`` yourself."""
    unguarded = _unguarded_calls(path)
    assert not unguarded, (
        f"{path.relative_to(REPO)} calls slim_optimize() at line(s) {unguarded} with no "
        "isnan check nearby. An infeasible LP returns nan, not None, and every ordinary "
        "guard against it is False -- so the value flows on and becomes a reported number. "
        "Call `growth_or_none(model)` from ystwin.fba.solver instead")


def test_the_exemption_is_the_module_that_defines_the_guard():
    """If solver.py stops checking for nan, nothing else in the repository does."""
    source = (REPO / "src" / "ystwin" / "fba" / "solver.py").read_text()
    assert "isnan" in source and "def growth_or_none" in source


def test_every_guard_the_wrong_way_round_is_actually_wrong():
    """The premise, pinned: these are the five expressions that were in the code."""
    nan = float("nan")
    assert str(nan or 0.0) == "nan"        # the fallback that never falls back
    assert (nan <= 0.0) is False
    assert (nan < 0.1 - 1e-6) is False
    assert (nan is None) is False
    assert str(nan * 0.9) == "nan"
