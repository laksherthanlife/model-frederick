"""Refusing to replace a tracked table with one built from fewer inputs.

Two scripts were found doing this within an hour of each other, and neither looked wrong.
``make_splits.py`` without the plate exports wrote a manifest missing 714 of 1650 rows and
exited 0; ``run_g4.py`` without two of three qPCR exports turned a measured anchor effect
into ``NaN`` and a replicate requirement into ``0``, which in the table is
indistinguishable from "there is no anchor at all".

Neither failed. That is the point: a partial rebuild produces a well-formed artefact whose
smallness is a fact about the machine that ran it, and nothing in the file says so.
"""

from __future__ import annotations

import pytest

from ystwin.gates.artefact import refuse_partial_rebuild


class TestItRefusesWhereWritingWouldDestroySomething:
    def test_fewer_inputs_over_an_existing_artefact_stops_the_run(self, tmp_path):
        table = tmp_path / "verdicts.csv"
        table.write_text("construct,effect\nUPRE1,0.4\n")

        with pytest.raises(SystemExit, match="1 of 3"):
            refuse_partial_rebuild(table, found=1, expected=3,
                                   what="qPCR replicates", remedy="Set YSTWIN_QPCR_RAW.")

    def test_the_message_carries_the_remedy_and_not_only_the_complaint(self, tmp_path):
        """A refusal a reader cannot act on is an obstacle. The variable that would supply
        the missing inputs is the whole content of the message."""
        table = tmp_path / "verdicts.csv"
        table.write_text("a\n1\n")

        with pytest.raises(SystemExit, match="YSTWIN_QPCR_RAW"):
            refuse_partial_rebuild(table, found=1, expected=3,
                                   what="qPCR replicates", remedy="Set YSTWIN_QPCR_RAW.")

    def test_it_says_how_many_are_missing_not_only_how_many_were_found(self, tmp_path):
        table = tmp_path / "t.csv"
        table.write_text("a\n1\n")

        with pytest.raises(SystemExit, match="missing 2"):
            refuse_partial_rebuild(table, found=1, expected=3, what="plates", remedy="x")


class TestItStaysOutOfTheWayWhereThereIsNothingToLose:
    def test_a_first_run_with_no_artefact_yet_is_allowed(self, tmp_path, capsys):
        """Refusing here would strand the very run that creates the file. It proceeds, and
        announces that what it is building is thin, so the artefact is not later mistaken
        for a complete one."""
        refuse_partial_rebuild(tmp_path / "absent.csv", found=1, expected=3,
                               what="qPCR replicates", remedy="Set YSTWIN_QPCR_RAW.")

        assert "1 of 3 qPCR replicates" in capsys.readouterr().out

    def test_a_complete_run_over_an_existing_artefact_is_allowed(self, tmp_path):
        table = tmp_path / "t.csv"
        table.write_text("a\n1\n")

        refuse_partial_rebuild(table, found=3, expected=3, what="plates", remedy="x")

    def test_more_inputs_than_expected_is_not_treated_as_a_shortfall(self, tmp_path):
        """A machine that finds a fourth replicate is not degrading anything, and a
        strict equality check would refuse the best case this gate could see."""
        table = tmp_path / "t.csv"
        table.write_text("a\n1\n")

        refuse_partial_rebuild(table, found=4, expected=3, what="plates", remedy="x")

    def test_a_complete_run_says_nothing_at_all(self, tmp_path, capsys):
        refuse_partial_rebuild(tmp_path / "absent.csv", found=3, expected=3,
                               what="plates", remedy="x")

        assert capsys.readouterr().out == ""


class TestTheGateIsWiredIntoTheScriptThatNeededIt:
    def test_run_g4_checks_before_it_writes_anything(self):
        """The order is the whole protection. A check after the first ``to_csv`` would
        leave three of the four tables replaced and refuse on the fourth."""
        import pathlib

        source = (pathlib.Path(__file__).resolve().parents[1]
                  / "scripts" / "run_g4.py").read_text()

        assert source.index("refuse_partial_rebuild(") < source.index("to_csv")
