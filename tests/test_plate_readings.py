"""The committed measurement matrix, and whether it is the workbook or only close to it.

``scripts/plate_readings.py`` exists so a fresh clone can re-derive a real-data number
without the four Synergy workbooks, which cannot be committed because each carries a named
private individual in its metadata. That makes the committed text the only copy of the
measurements a stranger sees, and the whole design rests on one claim: it round-trips *bit
for bit*, not to within a tolerance.

So the tests here are about exactness rather than agreement. Elapsed times are stored as
``H:MM:SS`` and rebuilt with the same ``h + m/60 + s/3600`` expression ``plate/synergy.py``
uses; a different expression for the same real number lands on a different float for
26 percent of the second-resolution timestamps in a day, which the round-trip check would
catch and no reader of the CSV would. Well ordering is per block because one export reads
its two channels in two different orders and that ordering reaches the row order of
``sensor_characterisation.csv``.

Nothing here loads ``run_sensor_characterisation.py``: this script's job is to supply that
module's readers, and the readers are what is checked.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pandas as pd
import pytest

from ystwin.plate import replay

from ystwin.plate.synergy import KineticBlock, SynergyRun

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "plate_readings.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("plate_readings", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _block(script, channel="OD600", wells=("A1", "A2"), values=((0.10, 0.20),
                                                                (0.12, 0.24)),
           times=("0:08:18", "1:08:18"), temperatures=(30.1, 30.2), **kw):
    """One kinetic block whose times are the floats the committed text decodes to."""
    data = pd.DataFrame(list(values), columns=list(wells))
    data.index = pd.Index([replay._from_hms(t) for t in times], name="time_h")
    return KineticBlock(
        channel=channel, fluorophore=kw.get("fluorophore", "OD600"),
        optics=kw.get("optics", "600"), sheet=kw.get("sheet", "OD600 (Raw)"),
        derived=kw.get("derived", False),
        blank_subtracted=kw.get("blank_subtracted", False),
        data=data, temperature_c=list(temperatures),
    )


class TestElapsedTimesSurviveAsFloatsAndNotAsRoundedStrings:
    @pytest.mark.parametrize("text", ["0:00:00", "0:01:01", "0:08:18", "1:07:25",
                                      "12:34:56", "25:00:03", "100:59:59"])
    def test_a_timestamp_decodes_and_re_encodes_to_the_same_string(self, script, text):
        assert script._hms(replay._from_hms(text)) == text

    @pytest.mark.parametrize("text", ["0:01:01", "0:08:18", "12:34:56"])
    def test_and_the_float_comes_back_equal_to_the_last_bit(self, script, text):
        hours = replay._from_hms(text)

        assert replay._from_hms(script._hms(hours)) == hours

    def test_the_decoding_expression_is_the_one_synergy_uses_and_not_an_equivalent(
            self, script):
        """``h + m/60 + s/3600`` and ``round(h*3600 + m*60 + s)/3600`` are the same real
        number and different floats for 22536 of the 86400 timestamps in a day. If this
        ever holds by accident the round-trip above stops being a proof of anything."""
        h, m, s = 0, 1, 1

        assert replay._from_hms("0:01:01") == h + m / 60.0 + s / 3600.0
        assert replay._from_hms("0:01:01") != round(h * 3600 + m * 60 + s) / 3600.0

    def test_a_time_that_is_not_a_whole_number_of_seconds_is_refused(self, script):
        """H:MM:SS cannot store it, so storing it would make the committed text lossy and
        the round-trip check downstream would fail without saying why."""
        with pytest.raises(ValueError, match="H:MM:SS"):
            script._hms(0.5 / 3600.0 + 1.0)


class TestFilenamesKeepWhatTheyHaveToDistinguish:
    def test_the_export_slug_matches_the_convention_the_other_scripts_use(self, script):
        assert script._slug("20260722_ER&OxidativeStress New Protocol") == \
            "20260722_ERandOxidativeStress_New_Protocol"

    def test_a_repeat_read_of_one_channel_keeps_its_own_file(self, script):
        """``mCitrine`` and ``mCitrine[2]`` are two reads of one fluorophore. A slug that
        dropped the index would put both blocks in one file and lose one of them."""
        assert script._channel_slug("mCitrine[2]") == "mCitrine_2"
        assert script._channel_slug("mCitrine") != script._channel_slug("mCitrine[2]")


class TestExactEqualityIsWhatIsChecked:
    @staticmethod
    def _frame(columns=("A1", "A2"), values=((1.0, 2.0), (3.0, 4.0)),
               index=(0.5, 1.5), name="time_h"):
        frame = pd.DataFrame(list(values), columns=list(columns))
        frame.index = pd.Index(list(index), name=name)
        return frame

    def test_two_identical_blocks_report_no_difference(self, script):
        assert script._identical(self._frame(), self._frame()) is None

    def test_a_difference_of_one_bit_is_a_difference(self, script):
        """``allclose`` would pass this, and the committed text is reproduced byte for byte
        or not at all -- a float that is nearly the reading is not the reading."""
        nudged = self._frame(values=((1.0, 2.0), (3.0, 4.0 + 2 ** -50)))

        assert "value differs" in script._identical(self._frame(), nudged)

    def test_the_same_values_in_a_different_well_order_is_a_difference(self, script):
        """The reason blocks are committed one file each: 20260804 reads its two channels
        in two different well orders, and that ordering sets the row order of the
        characterisation table."""
        swapped = self._frame(columns=("A2", "A1"))

        assert "well ordering differs" in script._identical(self._frame(), swapped)

    def test_a_shifted_timestamp_is_a_difference_even_where_the_readings_match(self, script):
        moved = self._frame(index=(0.5, 1.6))

        assert script._identical(self._frame(), moved) == "elapsed times differ"

    def test_the_index_name_synergy_gives_the_frame_is_part_of_the_comparison(self, script):
        renamed = self._frame(name="hours")

        assert "index name differs" in script._identical(self._frame(), renamed)

    def test_a_block_with_a_row_missing_is_a_difference(self, script):
        short = self._frame(values=((1.0, 2.0),), index=(0.5,))

        assert "elapsed times differ" in script._identical(self._frame(), short)


class TestAnAbsentOrIncompleteCommittedSetRefuses:
    def test_a_missing_manifest_names_both_ways_of_supplying_the_plates(self, script,
                                                                       tmp_path):
        """The message is the whole value of this refusal: a clone with neither the
        workbooks nor the text has to be told which of the two to go and get.

        It raises ``FileNotFoundError`` and not ``SystemExit``, which it did while this
        lived inside a script. A library that calls ``SystemExit`` takes the decision to
        stop away from every caller, and the replay is now reached from several."""
        with pytest.raises(FileNotFoundError, match="YSTWIN_PLATES"):
            replay.load_manifest(tmp_path)

        with pytest.raises(FileNotFoundError, match="plate_readings.py --export"):
            replay.load_manifest(tmp_path)

    def test_an_export_with_no_committed_blocks_is_refused_rather_than_returned_empty(
            self, script, tmp_path):
        manifest = pd.DataFrame([], columns=script.MANIFEST_COLUMNS)

        with pytest.raises(ValueError, match="absent.xlsx"):
            replay.load_run("absent.xlsx", tmp_path, manifest)

    def test_a_plate_with_no_dose_ladder_of_its_own_raises_rather_than_returning_nothing(
            self, script, tmp_path):
        """``collect()`` tells "this plate has no dose ladder" apart from "this plate is
        unusable" by catching ValueError from the reader. Returning an empty frame instead
        would send a plate with no ladder down the usable path."""
        manifest = pd.DataFrame([{"export": "p.xlsx", "dose_response_file": ""}])

        with pytest.raises(ValueError, match="no per-construct dose-response sheet"):
            replay.load_doses("p.xlsx", tmp_path, manifest)


class TestTheCommittedTextRebuildsTheRunItCameFrom:
    @pytest.fixture
    def exported(self, script, tmp_path):
        """Two exports written out through the real ``export``, with injected readers.

        The second export's two channels are in different well orders, which is the case
        the per-block layout exists for.
        """
        source = tmp_path / "raw"
        source.mkdir()
        runs = {}
        for name, blocks in (
            ("one.xlsx", [_block(script), _block(script, channel="mCitrine",
                                                 fluorophore="mCitrine",
                                                 optics="480,530", sheet="mCitrine")]),
            ("two.xlsx", [_block(script, wells=("A1", "A2")),
                          _block(script, channel="mCitrine[2]", wells=("A2", "A1"),
                                 fluorophore="mCitrine", optics="480,530",
                                 sheet="mCitrine")]),
        ):
            path = source / name
            path.write_bytes(name.encode())
            runs[path.name] = SynergyRun(source=path, blocks=tuple(blocks))
        doses = pd.DataFrame({"construct": ["UPRE1", "UPRE1"], "dose_mM": [0.0, 0.5],
                              "time_h": [1.0, 1.0], "signal": [100.0, 250.0]})

        dest = tmp_path / "committed"
        paths = sorted(source.glob("*.xlsx"))
        manifest = script.export(paths, dest, lambda p: runs[p.name], lambda p: doses)
        return paths, dest, runs, manifest

    def test_every_block_of_every_export_is_reported_in_the_manifest(self, script,
                                                                    exported):
        _, _, runs, manifest = exported

        assert len(manifest) == sum(len(run.blocks) for run in runs.values())
        assert list(manifest.columns) == script.MANIFEST_COLUMNS

    def test_the_rebuilt_run_equals_the_workbook_value_by_value(self, script, exported):
        paths, dest, runs, _ = exported

        assert script.verify(paths, dest, lambda p: runs[p.name]) is True

    def test_a_single_edited_reading_makes_the_check_refuse(self, script, exported):
        """The check has to fail on the one thing it cannot see any other way: text that
        parses, has the right shape, and holds a number the reader never wrote."""
        paths, dest, runs, manifest = exported
        target = dest / manifest.readings_file.iloc[0]
        frame = pd.read_csv(target)
        frame.iloc[0, 2] = float(frame.iloc[0, 2]) + 1e-12
        frame.to_csv(target, index=False)

        assert script.verify(paths, dest, lambda p: runs[p.name]) is False

    def test_each_block_keeps_its_own_well_order_rather_than_a_union_header(self, script,
                                                                           exported):
        paths, dest, runs, manifest = exported

        rebuilt = replay.load_run("two.xlsx", dest, manifest)

        assert [list(b.data.columns) for b in rebuilt.blocks] == [["A1", "A2"],
                                                                 ["A2", "A1"]]

    def test_and_the_blocks_come_back_in_the_order_the_reader_produced_them(self, script,
                                                                           exported):
        """Manifest order preserves exported channel identities, not a raw-selection tie break."""
        paths, dest, runs, manifest = exported

        rebuilt = replay.load_run("one.xlsx", dest, manifest)

        assert [b.channel for b in rebuilt.blocks] == \
            [b.channel for b in runs["one.xlsx"].blocks]

    def test_the_block_metadata_survives_and_is_not_reconstructed_from_the_filename(
            self, script, exported):
        """``derived`` and ``blank_subtracted`` decide whether a block may be used at all,
        and both are booleans written as text -- a truthiness bug on ``"False"`` would turn
        every block into a derived one."""
        _, dest, runs, manifest = exported

        rebuilt = replay.load_run("one.xlsx", dest, manifest)

        original = runs["one.xlsx"].blocks[1]
        assert rebuilt.blocks[1].fluorophore == original.fluorophore
        assert rebuilt.blocks[1].optics == original.optics
        assert rebuilt.blocks[1].sheet == original.sheet
        assert rebuilt.blocks[1].derived is False
        assert rebuilt.blocks[1].blank_subtracted is False
        assert rebuilt.blocks[1].temperature_c == list(original.temperature_c)

    def test_the_dose_ladder_committed_beside_the_readings_reads_back(self, script,
                                                                     exported):
        """Kept because dropping it changes the answer rather than the provenance: without
        it ``collect()`` takes the ladder from the logbook registry instead of the plate."""
        _, dest, _, manifest = exported

        doses = replay.load_doses("one.xlsx", dest, manifest)

        assert list(doses.dose_mM) == [0.0, 0.5]
        assert list(doses.time_h) == [1.0, 1.0]


class TestReExportingCleansUpAfterItselfAndNothingElse:
    @pytest.fixture
    def dest(self, script, tmp_path):
        source = tmp_path / "raw"
        source.mkdir()
        (source / "first.xlsx").write_bytes(b"first")
        dest = tmp_path / "committed"
        script.export([source / "first.xlsx"], dest,
                      lambda p: SynergyRun(source=p, blocks=(_block(script),)),
                      lambda p: (_ for _ in ()).throw(ValueError("no ladder")))
        return source, dest

    def test_a_renamed_export_does_not_leave_its_old_table_behind(self, script, dest):
        """An orphan readings file is a measurement nothing points at, and the next reader
        cannot tell it from a current one.

        A rename is **the same bytes under a new name**, and the manifest already records
        the only thing that can say so: the workbook's sha256. Matching on the name alone
        would keep the old table; matching on the digest retires it."""
        source, committed = dest
        (source / "renamed.xlsx").write_bytes(b"first")   # same content, new name

        script.export([source / "renamed.xlsx"], committed,
                      lambda p: SynergyRun(source=p, blocks=(_block(script),)),
                      lambda p: (_ for _ in ()).throw(ValueError("no ladder")))

        assert not (committed / "first__OD600.csv").exists()
        assert (committed / "renamed__OD600.csv").exists()

    def test_a_different_export_leaves_the_other_plate_set_alone(self, script, dest):
        """The change this class exists to pin, and the opposite of a rename.

        `export` used to clear every file the manifest named, which made it all-or-nothing:
        the two wet-lab plate sets live in two directories, usually only one is configured,
        and exporting from one deleted the other's committed text. It did that here once,
        for real, and the tracked measurements are the only wet-lab data in this repository.

        Different bytes under a different name is a different export. It cannot supersede
        anything, so nothing else is touched."""
        source, committed = dest
        (source / "second.xlsx").write_bytes(b"second")

        script.export([source / "second.xlsx"], committed,
                      lambda p: SynergyRun(source=p, blocks=(_block(script),)),
                      lambda p: (_ for _ in ()).throw(ValueError("no ladder")))

        assert (committed / "first__OD600.csv").exists()
        assert (committed / "second__OD600.csv").exists()

    def test_and_the_manifest_still_names_both(self, script, dest):
        """A kept file with no manifest row would be the orphan this class started with."""
        import pandas as pd

        source, committed = dest
        (source / "second.xlsx").write_bytes(b"second")

        script.export([source / "second.xlsx"], committed,
                      lambda p: SynergyRun(source=p, blocks=(_block(script),)),
                      lambda p: (_ for _ in ()).throw(ValueError("no ladder")))

        manifest = pd.read_csv(committed / script.MANIFEST, dtype=str)
        assert set(manifest.export) == {"first.xlsx", "second.xlsx"}

    def test_a_failed_read_leaves_the_committed_text_exactly_as_it_was(self, script, dest):
        """The bug that produced all of this. `export` cleared the old files at the top and
        then read the workbooks, so one unreadable .xlsx in the directory -- a plate map, a
        summary -- deleted every committed measurement and left no manifest. Reads happen
        first now, and nothing is removed until all of them have succeeded."""
        source, committed = dest
        before = sorted(p.name for p in committed.iterdir())
        (source / "broken.xlsx").write_bytes(b"broken")

        def refuse(path):
            raise ValueError("no kinetic channel blocks found")

        with pytest.raises(SystemExit, match="unchanged"):
            script.export([source / "broken.xlsx"], committed, refuse,
                          lambda p: (_ for _ in ()).throw(ValueError("no ladder")))

        assert sorted(p.name for p in committed.iterdir()) == before

    def test_undeclared_processing_is_refused_before_any_existing_export_is_changed(self, script, dest):
        source, committed = dest
        before = {p.name: p.read_bytes() for p in committed.iterdir()}
        path = source / "first.xlsx"
        unknown = SynergyRun(path, (_block(script, blank_subtracted=None),))

        with pytest.raises(ValueError, match="correction-state.*undeclared"):
            script.export([path], committed, lambda p: unknown,
                          lambda p: pytest.fail("processing must be checked before dose-sheet access"))

        assert {p.name: p.read_bytes() for p in committed.iterdir()} == before

    def test_but_a_table_no_manifest_ever_named_is_left_alone(self, script, dest):
        """Stale files go by name from the previous manifest and not by globbing ``*.csv``,
        so a mistyped ``--dir`` cannot delete somebody's unrelated tables."""
        source, committed = dest
        bystander = committed / "somebody_elses.csv"
        bystander.write_text("keep,me\n1,2\n")

        script.export([source / "first.xlsx"], committed,
                      lambda p: SynergyRun(source=p, blocks=(_block(script),)),
                      lambda p: (_ for _ in ()).throw(ValueError("no ladder")))

        assert bystander.read_text() == "keep,me\n1,2\n"

    def test_an_export_with_no_readable_dose_sheet_still_commits_its_readings(self, script,
                                                                             dest):
        """A plate can be usable without a ladder of its own, so the reader's refusal is
        recorded as an empty cell rather than losing the plate."""
        _, committed = dest

        manifest = pd.read_csv(committed / script.MANIFEST, dtype=str,
                               keep_default_na=False)

        assert list(manifest.dose_response_file) == [""]
        assert (committed / manifest.readings_file.iloc[0]).exists()


@pytest.mark.integration
class TestTheCommittedTextThatIsActuallyTracked:
    """The 17 files in ``data/plates``, which are what a clone reproduces from."""

    @pytest.fixture(scope="class")
    def committed(self, script):
        if not (script.PLATES_DIR / script.MANIFEST).exists():
            pytest.skip(f"no committed plate text at {script.PLATES_DIR}")
        return script.PLATES_DIR, replay.load_manifest(script.PLATES_DIR)

    def test_every_manifest_row_points_at_a_file_that_is_there(self, script, committed):
        directory, manifest = committed

        missing = [name for name in manifest.readings_file
                   if not (directory / name).exists()]
        missing += [name for name in set(manifest.dose_response_file) - {""}
                    if not (directory / name).exists()]
        assert missing == []

    def test_every_committed_timestamp_round_trips_through_the_stored_string(self, script,
                                                                            committed):
        """The claim the whole file format rests on, on the real timestamps rather than on
        chosen ones."""
        directory, manifest = committed

        for name in manifest.readings_file:
            for text in pd.read_csv(directory / name).elapsed_hms:
                assert script._hms(replay._from_hms(text)) == text

    def test_each_export_rebuilds_into_a_run_with_the_blocks_the_manifest_lists(
            self, script, committed):
        directory, manifest = committed

        for export in manifest.export.drop_duplicates():
            run = replay.load_run(export, directory, manifest)
            expected = manifest[manifest.export == export]
            assert [b.channel for b in run.blocks] == list(expected.channel)
            assert [len(b.data) for b in run.blocks] == [int(n) for n in expected.n_times]
            assert [len(b.data.columns) for b in run.blocks] == \
                [int(n) for n in expected.n_wells]
