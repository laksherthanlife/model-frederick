"""Some exports split the plate down the page, not across it.

The Synergy exporter usually writes every well as a column beside one ``Time`` column.
Two of the four biosensor replicates instead write one block per group of plate columns,
stacked vertically, and only the first block carries ``Time``. Reading the first block
alone returns 21 wells of 87 -- and because the blank wells sit in a later block, the
plate then looks like one with no blank and is refused. That is how two of four
biological replicates went missing, which halved n on every reported fold.
"""

import pandas as pd
import pytest

from ystwin.plate.synergy import _continuations, _extract_blocks


def _sheet(rows):
    return pd.DataFrame(rows)


def _stacked(n_times=4, second_wells=("A4", "A5")):
    """Two blocks: the first with Time, the second with neither Time nor times."""
    rows = [["Time", "A1", "A2"]]
    for i in range(n_times):
        rows.append([f"0:{i:02d}:00", 0.1 + i, 0.2 + i])
    rows.append([None, *second_wells])
    for i in range(n_times):
        rows.append([None, 1.1 + i, 1.2 + i])
    return _sheet(rows)


class TestItReattachesStackedColumnGroups:
    def test_both_blocks_are_read(self):
        blocks = _extract_blocks(_stacked(), "OD600")

        assert len(blocks) == 1, "the continuation joins the block above, not a new one"
        assert list(blocks[0].data.columns) == ["A1", "A2", "A4", "A5"]

    def test_the_continuation_gets_the_time_axis_above_it(self):
        block = _extract_blocks(_stacked(), "OD600")[0]

        assert len(block.data.index) == 4
        assert block.data["A4"].tolist() == [1.1, 2.1, 3.1, 4.1]

    def test_a_wide_export_is_unchanged(self):
        """The format that already worked must read identically."""
        rows = [["Time", "A1", "A2"]]
        for i in range(4):
            rows.append([f"0:{i:02d}:00", 0.1 + i, 0.2 + i])

        blocks = _extract_blocks(_sheet(rows), "OD600")

        assert len(blocks) == 1
        assert list(blocks[0].data.columns) == ["A1", "A2"]


class TestItRefusesRatherThanMisalign:
    def test_a_short_continuation_is_skipped(self):
        """Merging a group with fewer rows would pair readings with the wrong
        timepoints, and nothing downstream would catch it."""
        rows = [["Time", "A1"]]
        for i in range(4):
            rows.append([f"0:{i:02d}:00", 0.1 + i])
        rows.append([None, "A4"])
        rows.append([None, 9.9])            # one row against the block's four

        block = _extract_blocks(_sheet(rows), "OD600")[0]

        assert "A4" not in block.data.columns

    def test_an_overlapping_group_is_left_alone(self):
        """Repeating the same wells is the raw/blank-subtracted pairing, a different
        thing, and must stay a separate block rather than be concatenated."""
        got = _continuations(_stacked(second_wells=("A1", "A2")), start=0,
                             n_times=4, taken={"A1", "A2"})

        assert got == []

    def test_a_group_above_the_block_is_not_pulled_up(self):
        got = _continuations(_stacked(), start=99, n_times=4, taken={"A1", "A2"})

        assert got == []


class TestTheBlankWellsAreTheReasonThisMatters:
    def test_a_blank_only_group_still_attaches(self):
        """On 20260804 the H1-H3 blanks are their own narrow block, far to the right and
        below the culture blocks. Losing it is what made the plate look unblanked."""
        rows = [["Time", "A1", "A2"]]
        for i in range(4):
            rows.append([f"0:{i:02d}:00", 0.5 + i, 0.6 + i])
        rows.append([None, "H1", "H2"])
        for i in range(4):
            rows.append([None, 0.09, 0.10])

        block = _extract_blocks(_sheet(rows), "OD600")[0]

        assert {"H1", "H2"} <= set(block.data.columns)
        assert block.data["H1"].max() == pytest.approx(0.09)
