"""Recover an unannotated qPCR export using the layout the logbook records.

The 2026-07-24 run was exported without Target or Sample filled in, so the analysis
sheet every other replicate has does not exist for it. The logbook records the plate
layout, which is identical across all three runs, so the annotation can be restored
from well position instead of re-running or re-exporting anything.

Guarding the assumption: the layout is only applied where the data agree with it --
rows G and H empty, reference and target rows separated, and the third column of each
block behaving like a no-RT control.
"""

import numpy as np
import pytest

from ystwin.qpcr import QPCR_PLATE, read_positional_cq


def _export(tmp_path, name="raw.xlsx", grid=None):
    """A Bio-Rad style export: Well / Fluor / Target / Content / Sample / Cq."""
    from openpyxl import Workbook

    wb = Workbook()
    ws = wb.active
    ws.append([None, "Well", "Fluor", "Target", "Content", "Sample", "Biological Set", "Cq"])
    for row_i, row in enumerate("ABCDEFGH"):
        for col in range(1, 13):
            cq = None if grid is None else grid.get((row, col))
            ws.append([None, f"{row}{col:02d}", "SYBR", None, "Unkn", None, None, cq])
    path = tmp_path / name
    wb.save(path)
    return path


def _plausible_grid():
    """Reference rows near Cq 19, target rows split by gene, RT- column raised."""
    grid = {}
    for row in "ABC":                       # UBC
        for col in range(1, 13):
            grid[(row, col)] = 19.0 + (2.5 if col % 3 == 0 else 0.0)
    for row in "DEF":                       # Hac1 (cols 1-6) / TRX2 (cols 7-12)
        for col in range(1, 13):
            base = 24.0 if col <= 6 else 16.0
            grid[(row, col)] = base + (7.0 if col % 3 == 0 else 0.0)
    return grid


def test_the_recorded_layout_matches_the_logbook(): 
    assert QPCR_PLATE.reference_rows == {"A": 0.0, "B": 0.2, "C": 0.5}
    assert QPCR_PLATE.target_rows == {"D": 0.0, "E": 0.2, "F": 0.5}
    assert QPCR_PLATE.construct_columns["UPRE1"] == (1, 2, 3)
    assert QPCR_PLATE.roles == ("tech1", "tech2", "rt_minus")


def test_it_produces_the_same_shape_as_the_annotated_reader(tmp_path):
    tidy = read_positional_cq(_export(tmp_path, grid=_plausible_grid()))

    assert set(tidy.columns) >= {"construct", "target", "dose_mM", "tech_rep", "cq", "rt_minus_cq"}
    assert set(tidy.construct) == {"UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"}
    assert sorted(tidy.dose_mM.unique()) == [0.0, 0.2, 0.5]


def test_the_target_gene_follows_the_construct(tmp_path):
    tidy = read_positional_cq(_export(tmp_path, grid=_plausible_grid()))
    by_construct = tidy.groupby("construct").target.unique().apply(set).to_dict()

    assert by_construct["UPRE1"] == {"UBC", "Hac1"}
    assert by_construct["AlteredYap1"] == {"UBC", "TRX2"}


def test_the_third_column_of_each_block_becomes_the_no_rt_control(tmp_path):
    tidy = read_positional_cq(_export(tmp_path, grid=_plausible_grid()))
    row = tidy.query("construct == 'UPRE1' and target == 'Hac1' and dose_mM == 0.0").iloc[0]

    assert row.cq == pytest.approx(24.0)
    assert row.rt_minus_cq == pytest.approx(31.0)


def test_both_technical_replicates_are_kept(tmp_path):
    tidy = read_positional_cq(_export(tmp_path, grid=_plausible_grid()))
    one = tidy.query("construct == 'UPRE2' and target == 'UBC' and dose_mM == 0.5")

    assert sorted(one.tech_rep) == [1, 2]


def test_it_feeds_delta_delta_cq_directly(tmp_path):
    from ystwin.qpcr import delta_delta_cq

    tidy = read_positional_cq(_export(tmp_path, grid=_plausible_grid()))
    out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)

    assert np.isfinite(out.fold_change).all()


class TestItRefusesWhenTheDataContradictTheLayout:
    def test_data_in_the_rows_the_layout_says_are_empty_is_refused(self, tmp_path):
        grid = _plausible_grid()
        for col in range(1, 13):
            grid[("G", col)] = 20.0

        with pytest.raises(ValueError, match="rows G/H"):
            read_positional_cq(_export(tmp_path, grid=grid))

    def test_a_reference_row_that_is_not_lower_than_its_target_is_refused(self, tmp_path):
        grid = _plausible_grid()
        for row in "ABC":
            for col in range(1, 7):
                grid[(row, col)] = 30.0     # reference above the Hac1 target

        with pytest.raises(ValueError, match="reference"):
            read_positional_cq(_export(tmp_path, grid=grid))

    def test_an_export_that_already_carries_annotation_is_refused(self, tmp_path):
        from openpyxl import Workbook

        wb = Workbook()
        ws = wb.active
        ws.append([None, "Well", "Fluor", "Target", "Content", "Sample", "Biological Set", "Cq"])
        ws.append([None, "A01", "SYBR", "UBC", "Unkn", "UPRE1 0mM", None, 19.0])
        path = tmp_path / "annotated.xlsx"
        wb.save(path)

        with pytest.raises(ValueError, match="already annotated"):
            read_positional_cq(path)
