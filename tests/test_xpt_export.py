"""The fourth biosensor replicate: the instrument file as committed text, and its arithmetic.

Plate ``20260804`` sat outside the panel for one reason. Every ``.xlsx`` of it went through
Gen5's blank-subtraction transform, so its media blanks read within 0.002 of zero and there
was no background left to correct against -- and the biosensor panel therefore rested on
three biological replicates rather than four.

The archive still has the background, and this file is about the two claims that had to hold
before it could be used:

1. **The two files describe one plate.** Subtracting the per-timepoint mean of H1:H3 from the
   archive reproduces the committed subtracted export. That check is now runnable **from the
   committed text alone**, with no ``.xpt`` and no workbook, because both halves are tracked:
   the raw matrix and the transformed one, from the same run, one directory apart. It is the
   strongest evidence in the repository that ``plate/gen5.py`` reads what it says it reads.
2. **One plate is one biological replicate.** Both halves stay committed on purpose -- the
   subtracted export is the evidence of what the transform removed -- so something has to
   decide which of the two is the replicate. ``collect()`` does, by the plate they share, and
   the alternative is a panel that reports n=5 for four cultures.

The tests here are about that seam. ``tests/test_gen5_xpt.py`` is about the binary reader.
"""
from __future__ import annotations

import importlib.util
import os
import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.plate import gen5, replay
from ystwin.plate.layout import RECORDED_PLATES, plate_key

_REPO = pathlib.Path(__file__).resolve().parents[1]

ARCHIVE_EXPORT = "20260804_ER_Oxidative_Replicate4.xpt"
"""The instrument file's row in ``data/plates/manifest.csv``."""

SUBTRACTED_EXPORT = "20260804_ER&OxidativeStress_Replicate4.xlsx"
"""The blank-subtracted workbook's row, kept as the evidence of what Gen5's transform did."""

BLANKS = ("H1", "H2", "H3")
"""The media wells Gen5 subtracted, which the tests below re-derive rather than assume."""

EXPORT_HALF_STEP = 5e-4
"""Half the last digit a three-decimal absorbance export can carry; see test_gen5_xpt.py."""

FLOAT_SLACK = 1e-9


def _script(name: str):
    """Load one script as a module, the way every script test in this repository does."""
    spec = importlib.util.spec_from_file_location(name, _REPO / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def characterisation():
    return _script("run_sensor_characterisation")


@pytest.fixture(scope="module")
def readings():
    return _script("plate_readings")


def _frame(export: str, fluorophore: str) -> pd.DataFrame:
    """One committed block as a frame, found through the manifest rather than by filename."""
    run = replay.load_run(export)
    return run.channel(fluorophore).data


class TestTheArchiveAndTheExportAreOnePlate:
    """Re-derivable from the checkout: no ``.xpt``, no workbook, no environment variable."""

    @pytest.mark.parametrize("fluorophore", ["OD600", "mCitrine"])
    def test_the_superseded_corrected_export_cannot_be_selected_as_raw(self, fluorophore):
        run = replay.load_run(SUBTRACTED_EXPORT)

        assert run.channel(fluorophore).blank_subtracted is True
        with pytest.raises(KeyError, match="already blank-subtracted"):
            run.raw_channel(fluorophore)
        assert replay.load_run(ARCHIVE_EXPORT).raw_channel(fluorophore).blank_subtracted is False

    @pytest.mark.parametrize("fluorophore,tolerance", [
        ("mCitrine", 1 / 3 + FLOAT_SLACK),
        ("OD600", EXPORT_HALF_STEP + FLOAT_SLACK),
    ])
    def test_raw_minus_the_mean_of_h1_h3_is_the_subtracted_export(self, fluorophore,
                                                                  tolerance):
        """The whole justification for committing the archive, checked from committed text.

        Fluorescence is the sharper half. Gen5 stores RFU as integers and rounds the
        subtracted value back to an integer, and the blank is the mean of three integers,
        so every disagreement must be zero or exactly a third -- which is asserted below
        rather than absorbed into a tolerance. Absorbance cannot be exact at all: the
        archive holds four decimals and the export three, so half of the export's last
        digit is the arithmetic maximum.
        """
        raw = _frame(ARCHIVE_EXPORT, fluorophore)
        subtracted = _frame(SUBTRACTED_EXPORT, fluorophore)
        wells = list(subtracted.columns)

        blank = raw[list(BLANKS)].mean(axis=1)
        residual = raw[wells].sub(blank, axis=0).to_numpy() - subtracted.to_numpy()

        assert float(np.max(np.abs(residual))) <= tolerance
        if fluorophore == "mCitrine":
            thirds = np.abs(residual) * 3.0
            assert np.allclose(thirds, np.round(thirds), atol=1e-6)

    def test_the_two_blocks_share_a_time_axis_to_the_last_bit(self):
        """Two readers, two file formats, one committed spelling of elapsed time.

        ``elapsed_ms / 3_600_000`` and ``h + m/60 + s/3600`` are the same real number and
        different floats -- they disagree at 16 of this plate's 50 timestamps -- so the
        archive reader decodes ``H:MM:SS`` rather than dividing the milliseconds out. If it
        did not, the export's own round-trip check would refuse to commit the result.
        """
        for fluorophore in ("OD600", "mCitrine"):
            raw = _frame(ARCHIVE_EXPORT, fluorophore).index.to_numpy(dtype=float)
            subtracted = _frame(SUBTRACTED_EXPORT, fluorophore).index.to_numpy(dtype=float)

            assert np.array_equal(raw, subtracted), fluorophore

    def test_the_background_is_in_one_file_and_gone_from_the_other(self):
        """What n=3 rested on. The blanks are medium in the archive and zero in the export."""
        raw_od = _frame(ARCHIVE_EXPORT, "OD600")[list(BLANKS)].to_numpy()
        raw_rfu = _frame(ARCHIVE_EXPORT, "mCitrine")[list(BLANKS)].to_numpy()
        export_od = _frame(SUBTRACTED_EXPORT, "OD600")[list(BLANKS)].to_numpy()

        assert raw_od.min() > 0.05, "media blanks read around 0.09 OD in the archive"
        assert raw_rfu.min() > 100.0, "and hundreds of RFU"
        assert np.abs(export_od).max() < 0.005, "and within noise of zero in the export"

    def test_the_archive_carries_the_whole_plate_and_the_export_does_not(self):
        """96 wells against 87: row H's water wells were never exported."""
        raw = _frame(ARCHIVE_EXPORT, "OD600")
        subtracted = _frame(SUBTRACTED_EXPORT, "OD600")

        assert len(raw.columns) == 96
        assert set(subtracted.columns) < set(raw.columns)


class TestWhichWellsRowHHolds:
    """The logbook said H4-H6 and the instrument says otherwise, so the entry moved."""

    def test_the_fourth_replicate_is_corrected_against_its_instrument_file(self):
        plate = RECORDED_PLATES["20260804"]

        assert plate.blank_wells == BLANKS
        assert "instrument file settles it" in plate.note

    def test_only_h1_h3_hold_medium_and_the_rest_of_row_h_holds_none(self):
        """Why the correction is a reading rather than a preference.

        Every well of row H sits at the density floor and stays there, so optical density
        cannot tell them apart. The reporter channel can: medium autofluoresces at 480/530
        and water does not, and the two groups differ by roughly a factor of six.
        """
        rfu = _frame(ARCHIVE_EXPORT, "mCitrine")
        medium = rfu[list(BLANKS)].to_numpy().mean()
        rest = rfu[[f"H{i}" for i in range(4, 13)]].to_numpy().mean()

        assert medium > 300.0
        assert rest < 100.0
        assert medium > 5 * rest


class TestOnePlateIsOneBiologicalReplicate:
    """``n_plates`` is what every interval rests on, so a plate counted twice invents one."""

    def test_both_files_resolve_to_the_same_logbook_plate(self):
        assert plate_key(ARCHIVE_EXPORT) == plate_key(SUBTRACTED_EXPORT) == "20260804"

    def test_a_file_outside_the_registry_has_no_plate(self):
        """The July exports and the BY4741 controls pass through untouched."""
        assert plate_key("20260701_ER_preliminary (RAW).xlsx") is None

    def test_the_instrument_file_is_the_one_the_analysis_keeps(self, characterisation):
        kept, set_aside = characterisation.one_file_per_plate(
            [pathlib.Path(SUBTRACTED_EXPORT), pathlib.Path(ARCHIVE_EXPORT)])

        assert [p.name for p in kept] == [ARCHIVE_EXPORT]
        assert set_aside == [SUBTRACTED_EXPORT]

    def test_the_order_the_two_arrive_in_does_not_decide_it(self, characterisation):
        """List order is not evidence; a rule that depended on it would be a coin toss."""
        kept, _ = characterisation.one_file_per_plate(
            [pathlib.Path(ARCHIVE_EXPORT), pathlib.Path(SUBTRACTED_EXPORT)])

        assert [p.name for p in kept] == [ARCHIVE_EXPORT]

    def test_two_workbooks_for_one_plate_are_refused_rather_than_ranked(self,
                                                                       characterisation):
        """Nothing distinguishes them, so picking either would be a guess about which
        culture the replicate is."""
        with pytest.raises(ValueError, match="One plate is one biological replicate"):
            characterisation.one_file_per_plate([
                pathlib.Path("20260722_a.xlsx"), pathlib.Path("20260722_b.xlsx")])

    def test_files_that_are_not_logbook_plates_are_all_kept(self, characterisation):
        names = ["20260701_ER_preliminary (RAW).xlsx", "20260709_ER_stress_1st (RAW).xlsx"]
        kept, set_aside = characterisation.one_file_per_plate(
            [pathlib.Path(n) for n in names])

        assert [p.name for p in kept] == names
        assert set_aside == []

    def test_the_committed_manifest_never_offers_two_files_for_one_plate(self,
                                                                        characterisation):
        """The invariant on ``data/plates`` itself, not on a constructed list.

        This is the check that would fire if a future export re-added a workbook beside an
        archive it cannot improve on, or committed a second copy of a plate under a new
        name.
        """
        kept, _ = characterisation.one_file_per_plate(replay.exports())
        keys = [plate_key(p.name) for p in kept]

        assert len([k for k in keys if k]) == len({k for k in keys if k})

    def test_all_four_current_public_sources_prepare_from_raw_channels_once(self, characterisation, monkeypatch):
        monkeypatch.setattr(characterisation, "read_synergy_kinetic",
                            lambda path: replay.load_run(path.name))
        exports = replay.exports(source_set="newprotocol")
        assert len(exports) == 4
        for path in exports:
            run = replay.load_run(path.name)
            blanks = run.recorded_plate.blank_wells
            prepared = characterisation.prepare(path, recorded_blanks=blanks)
            assert prepared is not None
            od, rfu, normalised, _, od_blank, rfu_blank = prepared
            aligned = run.aligned(run.raw_channel("OD600").channel)
            expected_od = aligned[run.raw_channel("OD600").channel]
            expected_rfu = aligned[run.raw_channel("mCitrine").channel]
            assert set(od.columns) == set(run.recorded_plate.culture_wells)
            pd.testing.assert_frame_equal(od, expected_od[od.columns])
            pd.testing.assert_frame_equal(rfu, expected_rfu[rfu.columns])
            assert od_blank == float(expected_od[list(blanks)].to_numpy().mean())
            assert rfu_blank == float(expected_rfu[list(blanks)].to_numpy().mean())
            pd.testing.assert_frame_equal(normalised, (rfu - rfu_blank) / (od - od_blank))

    def test_the_panel_is_four_plates_and_they_are_four_different_cultures(self):
        """The result this track exists for, asserted against the shipped table."""
        table = pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")

        assert table.plate.nunique() == 4
        assert sorted(table.groupby("construct").plate.nunique().unique()) == [4]


class TestTheReaderDispatchesOnFormatAndRefusesToInvent:
    def test_a_dose_sheet_is_refused_for_an_archive_the_way_collect_listens_for(
            self, readings):
        """``collect()`` tells "no ladder of its own" from "unusable" by catching
        ``ValueError``, so the absence has to be raised and not returned empty."""
        with pytest.raises(ValueError, match="no per-construct derived sheet"):
            readings._no_dose_sheet(pathlib.Path(ARCHIVE_EXPORT))

    def test_the_dispatch_is_on_the_file_format_and_nothing_else(self, readings):
        read = readings._by_suffix(lambda p: ("workbook", p), lambda p: ("archive", p))

        assert read("plate.xlsx")[0] == "workbook"
        assert read("plate.XPT")[0] == "archive"

    def test_the_instrument_file_this_panel_needs_is_named_and_committed(self, readings):
        """The list is short on purpose -- see ``INSTRUMENT_FILES`` -- and this is the one
        entry, so a rename that orphaned it would show up here rather than as a silent
        return to n=3."""
        assert ARCHIVE_EXPORT in readings.INSTRUMENT_FILES
        assert ARCHIVE_EXPORT in set(replay.load_manifest().export)


class TestTheArchiveReaderProducesTheShapeThePipelineSpeaks:
    """Unit-level: no ``.xpt`` needed, because the failure modes are about the container."""

    @staticmethod
    def _channel(elapsed_ms, values=None, name="OD600:600", temperature=None):
        n = 1 if elapsed_ms is None else len(elapsed_ms)
        return gen5.Gen5Channel(
            name=name, fluorophore=name.partition(":")[0], optics=name.partition(":")[2],
            repeat=1, wells=("A1", "A2"),
            values=np.array(values if values is not None else [[0.1, 0.2]] * n, dtype=float),
            elapsed_ms=None if elapsed_ms is None else np.array(elapsed_ms, dtype=np.int64),
            temperature_c=temperature, started_at=None, rows=1, cols=2)

    def _read(self, *channels):
        return gen5.Gen5Read(source=pathlib.Path("synthetic.xpt"), protocol=None,
                             channels=tuple(channels))

    def test_the_time_axis_is_the_one_the_committed_text_decodes_to(self):
        """``548000 / 3_600_000`` is a different float from ``9/60 + 8/3600``. The committed
        text stores ``0:09:08``, so the reader has to land on the second one or the export's
        exactness check refuses the result."""
        run = gen5.as_synergy_run(self._read(self._channel([548_000, 4_148_000])))

        assert list(run.blocks[0].data.index) == [9 / 60 + 8 / 3600, 1 + 9 / 60 + 8 / 3600]
        assert list(run.blocks[0].data.index)[0] != 548_000 / 3_600_000

    def test_a_stamp_that_is_not_a_whole_number_of_seconds_is_refused(self):
        """``H:MM:SS`` cannot carry it, so committing it would make the text lossy in a way
        the round-trip check would report without saying why."""
        with pytest.raises(gen5.Gen5FormatError, match="whole number of seconds"):
            gen5.as_synergy_run(self._read(self._channel([548_400])))

    def test_an_endpoint_read_is_refused_rather_than_put_on_a_read_number_axis(self):
        with pytest.raises(gen5.Gen5FormatError, match="endpoint read"):
            gen5.as_synergy_run(self._read(self._channel(None)))

    def test_archive_processing_state_comes_from_the_format_contract_not_the_sign(self):
        raw = gen5.as_synergy_run(self._read(self._channel([548_000])))
        negative = gen5.as_synergy_run(
            self._read(self._channel([548_000], values=[[-0.01, 0.2]])))

        assert raw.blocks[0].blank_subtracted is False
        assert negative.blocks[0].blank_subtracted is False
        assert negative.raw_channel("OD600").data.iloc[0, 0] == -0.01

    def test_a_repeat_read_is_labelled_the_way_the_workbook_reader_labels_it(self):
        """One convention for ``mCitrine[1]``/``[2]``, shared with ``plate/synergy.py``,
        because the committed text of one source has to compare against the other."""
        run = gen5.as_synergy_run(self._read(
            self._channel([548_000], name="mCitrine:480,530"),
            self._channel([548_000], name="mCitrine:480,530[2]")))

        assert [b.channel for b in run.blocks] == ["mCitrine[1]", "mCitrine[2]"]
        assert [b.fluorophore for b in run.blocks] == ["mCitrine", "mCitrine"]

    def test_the_sheet_field_carries_the_archives_own_data_set_name(self):
        """A ``.xpt`` has no sheets. The locator that plays the same role is the name the
        instrument gave the read, and it goes in as itself rather than as a fabricated
        sheet name."""
        run = gen5.as_synergy_run(self._read(self._channel([548_000])))

        assert run.blocks[0].sheet == "OD600:600"
        assert run.blocks[0].derived is False


@pytest.mark.integration
class TestTheCommittedTextIsTheInstrumentFile:
    """Needs the ``.xpt``; skipped where it is absent, as ``test_gen5_xpt.py`` does."""

    def test_the_archive_reproduces_its_committed_text_exactly(self):
        directory = os.environ.get("YSTWIN_GEN5_XPT")
        if not directory:
            pytest.skip("set YSTWIN_GEN5_XPT to exercise the binary reader")
        path = pathlib.Path(directory).expanduser() / ARCHIVE_EXPORT
        if not path.exists():
            pytest.skip(f"not present: {path}")

        run = gen5.read_xpt_run(path)
        committed = replay.load_run(ARCHIVE_EXPORT)

        assert [b.channel for b in run.blocks] == [b.channel for b in committed.blocks]
        for live, text in zip(run.blocks, committed.blocks):
            assert list(live.data.columns) == list(text.data.columns)
            assert np.array_equal(live.data.index.to_numpy(dtype=float),
                                  text.data.index.to_numpy(dtype=float))
            assert np.array_equal(live.data.to_numpy(), text.data.to_numpy())
