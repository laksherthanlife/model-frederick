"""Every function that takes a measured channel must say which correction state it wants.

`docs/ARCHITECTURE.md` §6.3 carried three tables listing, function by function, which ones
expected a raw reading and which expected a blank-corrected one. Every entry was true. None
was enforced, so the tables were a promise about what callers would remember, and §6.3 said
so in its own heading: "conventions two modules must agree on, **unenforced**".

`ystwin.readings` made the state a type. This is what stops it drifting back: a new function
taking `optical_density` or `rfu` as a bare array would restore the old hazard silently, and
the only thing that would notice is a reviewer who had read a forty-line table.

The test is deliberately about the SIGNATURE and not about behaviour. Whether a function
handles its readings correctly is what the rest of the suite is for; this asks only that the
question "raw or corrected?" has an answer written down where the compiler can see it.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

_SRC = pathlib.Path(__file__).resolve().parents[1] / "src" / "ystwin"

#: Parameter names that carry a measured optical channel. Anything named one of these is
#: a reading off a plate and has a correction state, whether or not anyone stated it.
_MEASURED = {
    "optical_density", "measured_od", "od", "expected_od", "measured_above_blank",
    "rfu", "reporter", "total_signal", "biomass",
}

_READING_TYPES = {"RawOD", "CorrectedOD", "RawRFU", "CorrectedRFU", "SpecificFluorescence"}

#: `(function, parameter)` pairs that carry a measured-channel NAME without being a reading.
#: Each needs a reason, because an exemption list is how a check like this quietly stops
#: checking anything -- the same rule `test_architecture_doc_matches_the_code.py` applies to
#: its own list.
_NOT_A_READING = {
    # A dilution factor times a stock density: what the reader SHOULD report. Nothing
    # measured it, so no blank has or has not come off it.
    ("calib/od_linearity.py", "fit_linear_range", "expected_od"): "computed from dilution factors, not measured",
    # Biomass in g/L is a state variable on the simulation side, not a channel. It is the
    # quantity `gdcw_per_od` would convert a density INTO, and nothing has measured that
    # factor -- see SpecificFluorescence's docstring.
    ("observation.py", "observe_rfu", "biomass"): "g/L biomass, a state, not a reading",
    ("observation.py", "observe_absorbance", "biomass"): "g/L biomass, a state, not a reading",
    ("observation.py", "product_from_absorbance", "biomass"): "g/L biomass, a state, not a reading",
    ("observation.py", "inner_filter_coeff_from_extinction", "biomass"): "g/L biomass, a state, not a reading",
    ("mech/contracts.py", "prior", "reporter"): "reporter presence flag or gene identifier, not a measured channel",
    ("mech/heat_zheng2016.py", "initial_state", "reporter"): "initial simulated YFP state in source-native arbitrary abundance units, not measured RFU",
}


def _public_functions():
    """Every public function under ``src/ystwin``, with its module path."""
    for path in sorted(_SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and not node.name.startswith("_"):
                yield path, node


def _measured_parameters():
    """``(module, function, parameter, annotation)`` for every measured-channel argument."""
    for path, node in _public_functions():
        for arg in node.args.args + node.args.kwonlyargs:
            if arg.arg not in _MEASURED:
                continue
            annotation = ast.unparse(arg.annotation) if arg.annotation else ""
            yield str(path.relative_to(_SRC)), node.name, arg.arg, annotation


class TestTheConventionIsEnforcedRatherThanDocumented:
    def test_every_measured_channel_is_annotated_with_a_reading_type(self):
        """The check §6.3 could not make. A bare array here is the old hazard returning."""
        unstated = [
            f"{module}::{function}({parameter}: {annotation or 'no annotation'})"
            for module, function, parameter, annotation in _measured_parameters()
            if not any(kind in annotation for kind in _READING_TYPES)
            and (module, function, parameter) not in _NOT_A_READING
        ]

        assert unstated == [], (
            f"{len(unstated)} measured-channel parameter(s) do not say which correction "
            f"state they want: {unstated}. Annotate with a ystwin.readings type, or add an "
            f"entry to _NOT_A_READING with a reason.")

    def test_the_check_covers_the_functions_it_is_meant_to(self):
        """A contract test that matched nothing would pass forever."""
        found = list(_measured_parameters())

        assert len(found) > 30
        assert len({module for module, *_ in found}) > 8

    def test_every_exemption_still_names_a_parameter_that_exists(self):
        """An exemption for a function that has been renamed or deleted is a lie that makes
        the list longer and the check weaker."""
        live = {(module, function, parameter)
                for module, function, parameter, _ in _measured_parameters()}

        assert set(_NOT_A_READING) <= live, (
            f"stale exemption(s): {sorted(set(_NOT_A_READING) - live)}")

    def test_every_exemption_states_a_reason(self):
        for key, reason in _NOT_A_READING.items():
            assert reason.strip(), f"{key} is exempt with no reason given"

    def test_reporter_selector_exemption_only_covers_flags_and_identifiers(self):
        selectors = [annotation for module, function, parameter, annotation in _measured_parameters()
                     if (module, function, parameter) == ("mech/contracts.py", "prior", "reporter")]
        assert sorted(selectors) == ["bool", "str"]


class TestTheTwoNamesThatUsedToCollide:
    """`calib/od.py` and `calib/od_linearity.py` both defined `fit_od_linearity`, taking a
    parameter called `measured_od` that meant RAW in one and BLANK-CORRECTED in the other,
    and returning different types. Nothing stopped a caller importing the wrong one."""

    def test_neither_module_defines_the_old_ambiguous_name(self):
        from ystwin.calib import od, od_linearity

        assert not hasattr(od, "fit_od_linearity")
        assert not hasattr(od_linearity, "fit_od_linearity")

    def test_each_is_named_for_what_it_returns(self):
        """Strings, not types: the package uses `from __future__ import annotations`."""
        from ystwin.calib.od import fit_od_calibration
        from ystwin.calib.od_linearity import fit_linear_range

        assert fit_od_calibration.__annotations__["return"] == "ODCalibration"
        assert fit_linear_range.__annotations__["return"] == "LinearityFit"

    def test_their_shared_parameter_name_now_carries_opposite_types(self):
        """Which is the whole point: `measured_od` still appears in both signatures, and
        handing either one the other's reading is now a TypeError rather than a number."""
        from ystwin.calib.od import fit_od_calibration
        from ystwin.calib.od_linearity import fit_linear_range

        assert fit_od_calibration.__annotations__["measured_od"] == "RawOD"
        assert fit_linear_range.__annotations__["measured_od"] == "CorrectedOD"


class TestABareArrayIsRefusedByTheFunctionsThatMatterMost:
    """Spot checks at the three entry points a real analysis actually goes through, so the
    contract is tested by behaviour and not only by reading annotations."""

    def test_the_growth_rate_refuses_a_bare_array(self):
        import numpy as np

        from ystwin.growth import specific_growth_rate

        with pytest.raises(TypeError, match="CorrectedOD"):
            specific_growth_rate(np.linspace(0, 10, 20), np.linspace(0.1, 1.0, 20))

    def test_the_activity_inversion_refuses_a_bare_array(self):
        import numpy as np

        from ystwin.reporter import promoter_activity

        with pytest.raises(TypeError, match="SpecificFluorescence"):
            promoter_activity(np.linspace(0, 10, 20), np.ones(20), np.full(20, 0.3))

    def test_the_optical_gate_refuses_a_corrected_density_where_it_wants_raw(self):
        """The direction a reader is least likely to expect: G1 subtracts the blank itself,
        so handing it a corrected trace would subtract it twice."""
        import numpy as np

        from ystwin.gates.g1_optical import assess_well
        from ystwin.readings import CorrectedOD

        with pytest.raises(TypeError, match="before the blank came off"):
            assess_well(np.linspace(0, 10, 20),
                        CorrectedOD(np.linspace(0.01, 0.5, 20)), od_blank=0.09)
