"""What the G4 anchor script does when a replicate is not on this machine.

The qPCR replicates live in two places and the split is not arbitrary: replicate 1's archive
entries were lowercased by the instrument, so it exists only as the copy
``scripts/repair_qpcr_export.py`` rebuilds into the repository, while replicates 2 and 3 open
as exported and stay wherever the lab keeps them. Either location can be absent, and the
column that says which replicate a fold change came from has to keep meaning the same date
on every machine.

That is the property pinned here. ``replicate`` is the index into ``_REPLICATES``, not the
position in whatever was found, so a machine with only the annotated exports reports them as
replicates 1 and 2 -- never renumbered to 0 and 1, which would silently relabel 11 August as
24 July in ``outputs/g4_anchor_fold_change.csv``.

Every test forces the absence with the documented environment variables rather than relying
on what this machine happens to hold, so it says the same thing on a machine with all three
replicates and on a fresh clone with none.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_g4.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("run_g4", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def locations(tmp_path, monkeypatch):
    """Put each qPCR location either at a real empty directory or at nothing.

    ``paths._resolve`` treats a variable that is set as authoritative even when it points at
    nothing, which is what makes this the documented way to test an absent dataset -- and
    forcing both states here is what keeps these tests saying the same thing on a machine
    that holds all three replicates and on a fresh clone.
    """
    def place(*, repaired: bool, lab: bool) -> None:
        for name, present in (("YSTWIN_QPCR", repaired), ("YSTWIN_QPCR_RAW", lab)):
            target = tmp_path / name
            if present:
                target.mkdir(exist_ok=True)
            monkeypatch.setenv(name, str(target))
    return place


class TestAReplicateKeepsItsNumberWhereverItLives:
    def test_each_replicate_reports_the_index_it_has_in_the_canonical_order(self, script):
        """Whatever this machine holds, index and tag have to agree with ``_REPLICATES``:
        that is what makes ``replicate`` mean a date rather than a position."""
        found = script.qpcr_files()

        assert all(script._REPLICATES[index][0] == tag for index, tag, _, _ in found)
        indices = [index for index, _, _, _ in found]
        assert indices == sorted(indices)

    def test_a_missing_repaired_replicate_does_not_renumber_the_others(self, script,
                                                                      locations):
        """The failure this guards: with replicate 1 absent, a comprehension that enumerated
        what it found would report 11 August as replicate 0, and the ``replicate`` column in
        the committed table would mean a different date per machine."""
        locations(repaired=False, lab=True)

        found = script.qpcr_files()

        assert [index for index, _, _, _ in found] == [1, 2]
        assert [tag for _, tag, _, _ in found] == ["qpcr_11Aug", "qpcr_13Aug"]

    def test_and_a_missing_lab_directory_leaves_the_repaired_one_as_replicate_zero(
            self, script, locations):
        locations(repaired=True, lab=False)

        found = script.qpcr_files()

        assert [(index, tag) for index, tag, _, _ in found] == [(0, "qpcr_24Jul")]

    def test_every_replicate_names_a_reader_the_script_actually_has(self, script):
        """``directories[how]`` is a dict lookup on a literal in ``_REPLICATES``, so a typo
        there is a KeyError at the top of the run rather than a wrong number."""
        assert {how for _, how, _ in script._REPLICATES} <= {"positional", "annotated"}


class TestNoAnchorMeansNoTableRatherThanAnEmptyOne:
    def test_both_locations_absent_refuses_and_names_both_variables(self, script,
                                                                   locations):
        """A gate with no anchor is not a gate. The message has to name both variables
        because the two locations hold different replicates for different reasons."""
        locations(repaired=False, lab=False)

        with pytest.raises(SystemExit, match="YSTWIN_QPCR_RAW"):
            script.load_anchor()

    def test_a_directory_that_resolves_but_holds_no_export_refuses_too(self, script,
                                                                      locations, capsys):
        """The nastier shape of the same absence: the variable points somewhere real and the
        file is not in it. The replicate is skipped by name, and with none left the run
        stops instead of writing a fold-change table with no rows."""
        locations(repaired=True, lab=False)

        with pytest.raises(SystemExit, match="no qPCR replicate could be read"):
            script.load_anchor()
        assert "not present; skipped" in capsys.readouterr().out


class TestNoReporterMeansNoTableEither:
    def test_absent_workbooks_replay_the_committed_text_rather_than_refusing(
            self, script, tmp_path, monkeypatch):
        """This used to exit, and exiting was the wrong answer. The Synergy workbooks are
        not committed -- their document properties name a private individual -- so the
        anchor could be checked only by the person who ran the plates. The same numbers are
        committed as text under `data/plates`, and the reader reaches them unchanged."""
        monkeypatch.setenv("YSTWIN_PLATES", str(tmp_path / "absent"))

        reporter = script.load_reporter()

        assert not reporter.empty
        assert set(reporter.construct) >= {"UPRE1", "UPRE2"}

    def test_an_empty_plate_directory_is_refused_rather_than_summarised(self, script,
                                                                       tmp_path,
                                                                       monkeypatch):
        """A resolvable but empty directory is the failure ``paths.py`` warns about -- being
        "left to discover an empty directory". It must not read as a plate set with no
        response in it."""
        empty = tmp_path / "plates"
        empty.mkdir()
        monkeypatch.setenv("YSTWIN_PLATES", str(empty))

        with pytest.raises(SystemExit, match="readable dose-response sheet"):
            script.load_reporter()


class TestTheConstructsTheReportWalksAreOnesTheProjectKnows:
    def test_every_construct_has_a_stressor(self, script):
        """The report prints ``STRESSOR_FOR_CONSTRUCT[construct]`` after two tables are
        already written, so a construct the map does not know fails the run halfway through
        with a bare KeyError and a half-written outputs directory."""
        assert set(script.CONSTRUCTS) <= set(STRESSOR_FOR_CONSTRUCT)

    def test_the_er_and_oxidative_constructs_are_both_represented(self, script):
        """G4 compares two anchors, Hac1 for the ER pair and TRX2 for the oxidative pair.
        A list that lost one would still run and would report one anchor as the gate."""
        agents = {STRESSOR_FOR_CONSTRUCT[c] for c in script.CONSTRUCTS}

        assert len(agents) == 2
