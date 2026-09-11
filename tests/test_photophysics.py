"""The photophysics leaf, and the proof that lifting it out of the generator moved nothing.

`citrine_ph_response` and `CITRINE_PKA` used to live in `generator/context.py`. They now live
in `ystwin/photophysics.py`, a depth-0 leaf that imports nothing from this package, so that
`observation.py` can reach the measurement channel without importing the simulator.

This file tests the *refactor*, not the curve. What the curve does -- half-maximal at the pKa,
0.95 at pH 7, monotone, 30% lost over one unit of acidification -- is already pinned in
`tests/test_environment_sweep.py` against the name `context.py` still exports, and that test
passing unchanged is itself half the evidence the move was clean. The other half is here:
the two names are the *same objects* in both places, and `observe_rfu` returns the same bits
it returned before the file existed.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

from ystwin import observation, photophysics
from ystwin.generator import context
from ystwin.observation import ReporterOptics, observe_rfu

_SRC = pathlib.Path(observation.__file__).resolve().parent


# --------------------------------------------------------------------------- #
# the refactor: one definition, two names
# --------------------------------------------------------------------------- #


def test_the_generator_re_exports_the_leaf_rather_than_carrying_a_second_copy():
    """Identity, not equality. A duplicated definition would pass an equality check.

    `scripts/run_environment_sweep.py` and `tests/test_environment_sweep.py` both import
    these from `generator.context`, so the names have to stay there. What must not happen is
    two definitions drifting apart -- so the test is `is`, which no copy can satisfy.
    """
    assert context.citrine_ph_response is photophysics.citrine_ph_response
    assert context.CITRINE_PKA == photophysics.CITRINE_PKA
    assert "citrine_ph_response" in context.__all__
    assert "CITRINE_PKA" in context.__all__


def test_photophysics_imports_nothing_from_this_package():
    """It is depth 0 or it is not a leaf, and the whole point of the move was the depth.

    Read off the file's own syntax tree rather than off `sys.modules`, because an import
    that happened to be satisfied by another test's import would not show up there.
    """
    tree = ast.parse(pathlib.Path(photophysics.__file__).read_text(encoding="utf-8"))
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            imported.append("." * node.level + (node.module or ""))
        elif isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)

    assert imported == ["__future__"]


def test_the_measurement_model_does_not_import_the_generator():
    """The coupling this refactor exists to prevent, asserted as a property of the file.

    `observation.py` and `generator/context.py` sit at the same depth, so importing the
    generator from the measurement model would have worked -- and would have dragged
    `stress_panel`, `culture` and `chemostat_physiology` behind every predicted RFU.
    """
    tree = ast.parse((_SRC / "observation.py").read_text(encoding="utf-8"))
    modules = [("." * node.level + (node.module or ""))
               for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]

    assert not [m for m in modules if "generator" in m]


# --------------------------------------------------------------------------- #
# bit-for-bit: `observe_rfu` before and after the move
# --------------------------------------------------------------------------- #

_PLAIN = dict(gain=1234.5, background=41.0, autofluorescence=3.25, inner_filter_coeff=None)
_PIGMENTED = dict(gain=1234.5, background=41.0, autofluorescence=3.25,
                  inner_filter_coeff=0.87, detector_max=260000.0)

_BEFORE_THE_MOVE = [
    # (optics, reporter_mature, biomass, carotenoid, reading) -- readings captured by
    # running `observe_rfu` on the commit before photophysics.py existed, as exact hex.
    (_PLAIN, 0.0, 0.011, 0.0, float.fromhex("0x1.4849374bc6a7fp+5")),
    (_PLAIN, 1e-09, 1.4, 0.0, float.fromhex("0x1.6c66674e5e560p+5")),
    (_PLAIN, 0.37, 1.4, 3.0, float.fromhex("0x1.5682b020c49bap+9")),
    (_PIGMENTED, 0.37, 1.4, 3.0, float.fromhex("0x1.616e609be2f38p+6")),
    (_PIGMENTED, 5.5, 0.011, 0.02, float.fromhex("0x1.c9bc7bcf2f22ap+6")),
    (_PIGMENTED, 10000.0, 1.4, 0.0, float.fromhex("0x1.fbd0000000000p+17")),
]


@pytest.mark.parametrize("optics,reporter,biomass,carotenoid,expected", _BEFORE_THE_MOVE)
def test_observe_rfu_is_bit_for_bit_what_it_was_before_the_refactor(
    optics, reporter, biomass, carotenoid, expected
):
    """Exact equality, not ``approx``. A pure refactor that moves a last bit is not pure.

    The expected values are hex float literals rather than decimals so that the comparison
    is over the bit pattern and no decimal round-trip can hide a one-ulp move. Six cases,
    covering an unpigmented well, a pigmented one, the attenuation branch on and off, and
    the detector ceiling.
    """
    assert observe_rfu(reporter, biomass, carotenoid, ReporterOptics(**optics)) == expected


def test_the_four_argument_call_stays_the_historical_path():
    """What makes the pin above survive Phase 1, stated as the invariant it rests on.

    Phase 0 moved the function; it did not wire it in. Multiplying `gain` by
    `citrine_ph_response(pH_c) / citrine_ph_response(pH_ref)` is Phase 1's B3 and belongs to
    whoever owns `observation.py`. When it lands, the four-argument call must still take the
    branch it takes today -- so every parameter past the four has to carry a default, and the
    default has to be the one that leaves the arithmetic alone. Today there are exactly four
    and no fifth; tomorrow there may be a fifth, and this still holds.
    """
    import inspect

    parameters = list(inspect.signature(observe_rfu).parameters.values())

    assert [p.name for p in parameters[:4]] == [
        "reporter_mature", "biomass", "carotenoid", "optics"]
    assert all(p.default is not inspect.Parameter.empty for p in parameters[4:])


# --------------------------------------------------------------------------- #
# the number that makes the seam worth opening
# --------------------------------------------------------------------------- #


def test_one_unit_of_acidification_costs_a_third_of_the_signal():
    """30% against induction folds of ~1.5: the artifact is the size of the effect.

    This is the module's reason to exist, so it is asserted where the module is, on the
    leaf's own name. Griesbeck 2001, PMID 11387331 for the pKa; the 1.5-fold is this
    project's own measured induction range.
    """
    loss = 1.0 - photophysics.citrine_ph_response(6.0) / photophysics.citrine_ph_response(7.0)

    assert loss == pytest.approx(0.3005, abs=1e-4)
    assert photophysics.CITRINE_PKA == 5.7
