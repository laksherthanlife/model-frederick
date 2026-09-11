"""The ATP-sensor script: how it names a condition, and what it does without the plates.

``scripts/analyse_atp_sensor.py`` reproduces the numbers in ``docs/ATP_SENSOR.md``. Its one
piece of real logic outside the printing is ``conditions``, which turns a Synergy layout
into ``(row, carbon ratio)`` keys -- and the ratio is derived from the well's column number
by arithmetic, so a plate laid out differently is mapped silently rather than refused.

``read_export`` is stubbed throughout, so these run on any machine. The last class covers
the two refusals, which is what a machine without the exports actually hits.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "analyse_atp_sensor.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("analyse_atp_sensor", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _export(layout, wells, hours=(0.0, 6.0, 12.0)):
    """A stand-in for one Synergy export: a layout and the two blanked channels."""
    index = pd.Index(list(hours), name="hour")
    od = pd.DataFrame({well: np.linspace(0.1, 1.5, len(hours)) * scale
                       for well, scale in wells.items()}, index=index)
    fluorescence = pd.DataFrame({well: np.linspace(100.0, 900.0, len(hours)) * scale
                                 for well, scale in wells.items()}, index=index)
    return layout, {"Blank OD600:600": od, "Blank mCitrine:480,530": fluorescence}


@pytest.fixture
def stub_export(script, monkeypatch):
    def _install(layout, wells, hours=(0.0, 6.0, 12.0)):
        payload = _export(layout, wells, hours)
        monkeypatch.setattr(script, "read_export", lambda path: payload)
        return payload

    return _install


class TestACarbonRatioIsReadOffTheColumnNumber:
    """Three wells per condition, four conditions per row, so the ratio is
    ``(column - 1) // 3 + 1``. The mapping is arithmetic over the layout rather than
    anything the export states, which makes it exactly the kind of thing that goes wrong
    quietly when a plate is laid out differently."""

    @staticmethod
    def _row_of_four(script):
        layout, wells = {}, {}
        for column in range(1, 13):
            layout[f"A{column}"] = f"sample{(column - 1) // 3}"
            wells[f"A{column}"] = 1.0
        return layout, wells

    def test_each_triplet_of_columns_is_one_carbon_ratio(self, script, stub_export):
        layout, wells = self._row_of_four(script)
        stub_export(layout, wells)

        got = script.conditions("ignored.xlsx", None)

        assert set(got) == {("A", ratio) for ratio in script.RATIOS.values()}

    def test_the_first_three_columns_are_the_glucose_only_condition(self, script,
                                                                    stub_export):
        """The ratio table is 1-indexed and the columns are too, so an off-by-one here
        would relabel every condition on the plate as its neighbour."""
        stub_export({"A1": "s", "A2": "s", "A3": "s"}, {"A1": 1.0, "A2": 1.0, "A3": 1.0})

        assert list(script.conditions("ignored.xlsx", None)) == [("A", "glu 10:0")]

    def test_the_group_is_taken_from_the_lowest_numbered_well_not_the_first_string(
            self, script, stub_export):
        """Wells are sorted by ``int(well[1:])``. Sorted as strings, "A10" would come
        before "A9" and a condition straddling the boundary would be filed under the next
        ratio along."""
        wells = {"A9": 1.0, "A10": 1.0, "A11": 1.0}
        stub_export({well: "straddles" for well in wells}, wells)

        assert list(script.conditions("ignored.xlsx", None)) == [("A", "gal 6:4")]

    def test_a_plate_wider_than_the_grid_is_refused_rather_than_folded(self, script,
                                                                       stub_export):
        """Column 13 maps to a fifth group the ratio table does not have. Raising is the
        right outcome: silently folding it into an existing condition would average two
        carbon sources together."""
        stub_export({"A13": "s"}, {"A13": 1.0})

        with pytest.raises(KeyError):
            script.conditions("ignored.xlsx", None)


class TestWhichWellsAreRead:
    def test_blanks_are_not_a_condition(self, script, stub_export):
        """"BLK" is the medium-only control. Read as a condition it would enter the fold
        range as a well with no cells in it."""
        wells = {"A1": 1.0, "A2": 1.0, "A3": 1.0}
        stub_export({"A1": "BLK", "A2": "sample", "A3": "sample"}, wells)

        assert list(script.conditions("ignored.xlsx", None)) == [("A", "glu 10:0")]

    def test_a_row_filter_keeps_only_the_named_rows(self, script, stub_export):
        """The 2026-08-09 plate holds two different builds, rows A-C and rows D-F, on one
        export. Reading both as one plate would compare two strains as if they were
        conditions of the same one."""
        wells = {"A1": 1.0, "A2": 1.0, "D1": 2.0, "D2": 2.0}
        stub_export({"A1": "top", "A2": "top", "D1": "bottom", "D2": "bottom"}, wells)

        assert list(script.conditions("ignored.xlsx", set("ABC"))) == [("A", "glu 10:0")]

    def test_no_filter_keeps_every_row(self, script, stub_export):
        wells = {"A1": 1.0, "D1": 2.0}
        stub_export({"A1": "top", "D1": "bottom"}, wells)

        assert set(script.conditions("ignored.xlsx", None)) == {("A", "glu 10:0"),
                                                                ("D", "glu 10:0")}

    def test_replicate_wells_are_averaged_rather_than_taken_one_at_a_time(self, script,
                                                                          stub_export):
        """Three technical replicates per condition, and the specific fluorescence is
        computed once from their mean rather than three times and averaged after."""
        wells = {"A1": 1.0, "A2": 3.0}
        stub_export({"A1": "s", "A2": "s"}, wells)

        _, od, fluorescence = script.conditions("ignored.xlsx", None)[("A", "glu 10:0")]

        assert od[0] == pytest.approx(0.1 * (1.0 + 3.0) / 2)
        assert fluorescence[0] == pytest.approx(100.0 * (1.0 + 3.0) / 2)

    def test_the_time_axis_comes_from_the_optical_density_table(self, script,
                                                                stub_export):
        """The two channels are read minutes apart and the script pairs them by row. The
        OD index is the one carried forward, so a fluorescence table with its own clock
        does not silently become the time axis."""
        stub_export({"A1": "s"}, {"A1": 1.0}, hours=(0.0, 4.0, 8.0))

        times, _, _ = script.conditions("ignored.xlsx", None)[("A", "glu 10:0")]

        assert list(times) == [0.0, 4.0, 8.0]


class TestNoPlatesMeansNoNumbers:
    """The script prints a comparison and writes nothing, so the only way it can mislead is
    by printing a table from part of the data. Both refusals happen before any reading."""

    def test_an_unset_location_exits_naming_the_variable(self, script, monkeypatch,
                                                         tmp_path):
        monkeypatch.setenv("YSTWIN_ATP_SENSOR", str(tmp_path / "nowhere"))

        with pytest.raises(SystemExit) as refusal:
            script.main()

        assert "YSTWIN_ATP_SENSOR" in str(refusal.value)

    def test_a_directory_missing_an_export_names_the_file_rather_than_reading_the_rest(
            self, script, monkeypatch, tmp_path):
        """Two exports, three plate blocks. Analysing whichever file happens to be present
        would print a "what the plates measured" table covering one plate under a heading
        that claims two."""
        (tmp_path / "20260805_ICL.xlsx").write_bytes(b"")
        monkeypatch.setenv("YSTWIN_ATP_SENSOR", str(tmp_path))

        with pytest.raises(SystemExit) as refusal:
            script.main()

        assert "20260809_ACS_ICL.xlsx" in str(refusal.value)
        assert "20260805_ICL.xlsx" not in str(refusal.value)


class TestThePlateTableDescribesTheRealExports:
    """``PLATES`` carries filenames rather than paths, so the directory is resolved at run
    time -- and the row sets are what split one export into two builds."""

    def test_no_entry_carries_a_path(self, script):
        """A path here would put one machine's directory layout back into the script, which
        is what ``paths.py`` exists to have removed."""
        for _, filename, _, _ in script.PLATES:
            assert "/" not in filename and "\\" not in filename

    def test_the_two_builds_on_one_export_do_not_share_a_row(self, script):
        """Rows A-C and D-F are different strains on the same plate. An overlap would put
        one strain's wells into both blocks and report the difference as a condition."""
        by_file = {}
        for _, filename, rows, _ in script.PLATES:
            if rows is not None:
                assert not by_file.get(filename, set()) & rows
                by_file[filename] = by_file.get(filename, set()) | rows

    def test_every_filtered_block_labels_exactly_the_rows_it_keeps(self, script):
        """The sugar table and the row filter are written separately, so a row kept with no
        label -- or labelled but filtered out -- is possible and would be invisible."""
        for _, _, rows, sugar in script.PLATES:
            if rows is not None:
                assert set(sugar) == rows

    def test_the_ratio_table_covers_the_whole_twelve_column_grid(self, script):
        """Four conditions of three wells. A gap would make ``conditions`` raise on a
        perfectly ordinary plate."""
        assert set(script.RATIOS) == {1, 2, 3, 4}
