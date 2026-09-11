"""Reader for the hand-annotated Cq sheets (2026-08 biological replicates).

Layout, located by content because the row offset differs between replicates:

    row k    :                      UPRE1                UPRE2            ...
    row k+1  :                      TechRep1 TechRep2 Average RT-Control  ...
    row k+2  : 0 mM   UBC        A  20.09    19.75    19.92   22.45       ...
    row k+3  : 0.2 mM            B  ...
    row k+5  : 0 mM   Hac1/TRX2  D  ...

The target cell reads "Hac1 / TRX2" because which gene it is depends on the column
block, so it is resolved per construct rather than taken literally.
"""

import numpy as np
import pandas as pd
import pytest

from ystwin.qpcr import read_annotated_cq

CONSTRUCTS = ["UPRE1", "UPRE2", "Native Yap1", "Altered Yap1"]


@pytest.fixture
def annotated_sheet(tmp_path):
    def _make(name="cq.xlsx", lead_rows=1, doses=("0 mM", "0.2 mM", "0.5 mM")):
        from openpyxl import Workbook

        wb = Workbook()
        ws = wb.active
        ws.title = "Results"
        for _ in range(lead_rows):
            ws.append([])
        header = [None] * 5
        for i, c in enumerate(CONSTRUCTS):
            header += [c, None, None, None]
        ws.append(header)
        sub = [None] * 5
        for _ in CONSTRUCTS:
            sub += ["Technical Replicate 1", "Technical Replicate 2", "Average", "RT- Control"]
        ws.append(sub)

        letters = iter("ABCDEF")
        for target in ("UBC", "Hac1 / TRX2"):
            for j, dose in enumerate(doses):
                row = [None, None, dose, target if j == 0 else None, next(letters)]
                for k, _ in enumerate(CONSTRUCTS):
                    base = 20.0 + k + (0.0 if target == "UBC" else 6.0) - j
                    row += [base, base + 0.2, base + 0.1, base + 8.0]
                ws.append(row)
        path = tmp_path / name
        wb.save(path)
        return path

    return _make


def test_every_construct_dose_and_target_combination_is_recovered(annotated_sheet):
    tidy = read_annotated_cq(annotated_sheet())

    assert set(tidy.construct) == {"UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"}
    assert sorted(tidy.dose_mM.unique()) == [0.0, 0.2, 0.5]


def test_construct_names_are_normalised_to_one_spelling(annotated_sheet):
    tidy = read_annotated_cq(annotated_sheet())

    assert "Native Yap1" not in set(tidy.construct)
    assert "NativeYap1" in set(tidy.construct)


def test_the_shared_target_cell_resolves_per_construct(annotated_sheet):
    """'Hac1 / TRX2' means Hac1 under the UPRE blocks and TRX2 under the Yap1 blocks."""
    tidy = read_annotated_cq(annotated_sheet())
    by_construct = tidy.groupby("construct").target.unique().apply(set).to_dict()

    assert by_construct["UPRE1"] == {"UBC", "Hac1"}
    assert by_construct["NativeYap1"] == {"UBC", "TRX2"}


def test_both_technical_replicates_are_kept_and_the_average_is_not_a_third(annotated_sheet):
    tidy = read_annotated_cq(annotated_sheet())
    one = tidy.query("construct == 'UPRE1' and target == 'UBC' and dose_mM == 0.0")

    assert sorted(one.tech_rep) == [1, 2]


def test_the_no_rt_control_travels_with_each_reading(annotated_sheet):
    tidy = read_annotated_cq(annotated_sheet())

    assert tidy.rt_minus_cq.notna().all()
    assert (tidy.rt_minus_cq > tidy.cq).all()


def test_a_varying_preamble_offset_does_not_shift_the_parse(annotated_sheet):
    shallow = read_annotated_cq(annotated_sheet(lead_rows=1))
    deep = read_annotated_cq(annotated_sheet(name="deep.xlsx", lead_rows=4))

    pd.testing.assert_frame_equal(
        shallow.sort_values(list(shallow.columns)).reset_index(drop=True),
        deep.sort_values(list(deep.columns)).reset_index(drop=True),
    )


def test_doses_are_parsed_as_numbers_not_labels(annotated_sheet):
    tidy = read_annotated_cq(annotated_sheet())

    assert tidy.dose_mM.dtype.kind == "f"


def test_a_sheet_with_no_recognisable_construct_header_is_refused(tmp_path):
    from openpyxl import Workbook

    wb = Workbook()
    wb.active.append(["Well", "Fluor", "Target", "Cq"])
    wb.active.append(["A01", "SYBR", None, 18.9])
    path = tmp_path / "raw_export.xlsx"
    wb.save(path)

    with pytest.raises(ValueError, match="no construct header"):
        read_annotated_cq(path)


def test_the_result_feeds_delta_delta_cq_directly(annotated_sheet):
    from ystwin.qpcr import delta_delta_cq

    tidy = read_annotated_cq(annotated_sheet())
    out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)

    assert set(out.construct) == {"UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"}
    assert np.isfinite(out.fold_change).all()
