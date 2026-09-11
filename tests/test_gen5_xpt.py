"""The Gen5 ``.xpt`` reader, checked against text this repository already trusts.

``src/ystwin/plate/gen5.py`` decodes an undocumented binary. Nothing in it can be
believed because it looks right; the only reason to believe any of it is that reading a
``.xpt`` reproduces the plate the ``.xlsx`` export of the *same run* produced, and that
export is committed under ``data/plates`` and was verified value-for-value against the
workbook by ``plate_readings.py --export``.

Four plates have both halves on this machine, and they are deliberately unalike -- 25, 110
and 123 timepoints; 45, 87 and 93 exported wells; two runs that ran to the end and two
that were stopped early:

    20260728_ER_Oxidative_Replicate2.xpt  -> ..._Replicate2__{OD600,mCitrine}.csv
    20260803_ER_Oxidative_Replicate3.xpt  -> ..._Replicate3__{OD600_1,mCitrine_2}.csv
    20260708_ER_stress_1st.xpt            -> 20260709_ER_stress_1st_(RAW)__*.csv
    24h_30min-interval_mCitrine-2.xpt     -> 20260701_ER_preliminary_(RAW)__*.csv

Fluorescence has to match **exactly**: Gen5 records RFU as integers and exports the same
integers. Absorbance cannot, and the size of the gap is itself the finding -- the ``.xpt``
holds four decimals and the ``.xlsx`` export holds three, so the two can differ by at most
half of the export's last digit, 0.0005, and never by more. That bound is asserted from
both sides: no value further than 0.0005 away, and at least one value that is not equal,
because an ``.xpt`` that agreed to the last bit would mean the fourth decimal was never
there and the whole point of reading this format would be gone.

Everything that needs a ``.xpt`` skips when ``YSTWIN_GEN5_XPT`` does not name a directory
holding them, in the way ``conftest._real_export`` skips for the workbooks: these files
are wet-lab exports on one machine, and the suite has to pass on a machine that has never
seen a plate reader. The pure-decoding tests below need no data at all.
"""
from __future__ import annotations

import datetime as _dt
import os
import pathlib

import numpy as np
import pytest

from ystwin import paths
from ystwin.plate import gen5

# Half of the last digit an .xlsx export of absorbance can carry. Gen5 writes OD to three
# decimals in the export and stores four in the .xpt, so this is the largest disagreement
# the two can have while describing one reading, and any larger one is a parser bug.
EXPORT_HALF_STEP = 5e-4

# Float slack on top of it. The worst observed gap is 0.000500000000000056: 0.0005 is not
# exactly representable, so the comparison needs room for the last bit and nothing more.
FLOAT_SLACK = 1e-9

# (.xpt name, channel in the archive, committed CSV under data/plates).
# The Replicate3 workbook has four sheets -- raw and blank-subtracted for each channel --
# and the raw ones are OD600[1] and mCitrine[2]; see data/plates/manifest.csv.
PROVEN = [
    ("20260728_ER_Oxidative_Replicate2.xpt", "OD600:600",
     "20260728_ERandOxidativeStress_NewProtocol_Replicate2__OD600.csv"),
    ("20260728_ER_Oxidative_Replicate2.xpt", "mCitrine:480,530",
     "20260728_ERandOxidativeStress_NewProtocol_Replicate2__mCitrine.csv"),
    ("20260803_ER_Oxidative_Replicate3.xpt", "OD600:600",
     "20260803_ERandoxidativestress_Replicate3__OD600_1.csv"),
    ("20260803_ER_Oxidative_Replicate3.xpt", "mCitrine:480,530",
     "20260803_ERandoxidativestress_Replicate3__mCitrine_2.csv"),
    ("20260708_ER_stress_1st.xpt", "OD600:600",
     "20260709_ER_stress_1st_(RAW)__OD600.csv"),
    ("20260708_ER_stress_1st.xpt", "mCitrine:480,530",
     "20260709_ER_stress_1st_(RAW)__mCitrine.csv"),
    ("24h_30min-interval_mCitrine-2.xpt", "OD600:600",
     "20260701_ER_preliminary_(RAW)__OD600.csv"),
    ("24h_30min-interval_mCitrine-2.xpt", "mCitrine:480,530",
     "20260701_ER_preliminary_(RAW)__mCitrine.csv"),
    # Replicate 1 of the biosensor panel. The file is named for its contents rather than
    # its date, which is why it went unrecognised until the .xpt sweep: it is the raw
    # instrument file behind the 2026-07-22 ANALYSED workbook. mCitrine matches the
    # committed text exactly, 2175 of 2175.
    ("UPRE_Yap1_StressTest_Replicate1.xpt", "OD600:600",
     "20260722_ERandOxidativeStress_NewProtocol_ANALYSED__OD600_1.csv"),
    ("UPRE_Yap1_StressTest_Replicate1.xpt", "mCitrine:480,530",
     "20260722_ERandOxidativeStress_NewProtocol_ANALYSED__mCitrine_1.csv"),
]

# The plate whose only surviving .xlsx is blank-subtracted, and the wells it was
# subtracted against. docs/DATA_INVENTORY.md calls recovering its background the cheapest
# remaining gain in the project, which is what this file is here to settle.
SUBTRACTED = "20260804_ER_Oxidative_Replicate4.xpt"
SUBTRACTED_BLANKS = ("H1", "H2", "H3")


def _xpt(name: str) -> pathlib.Path:
    """One named Gen5 experiment file, or a skip naming the variable that would find it.

    ``YSTWIN_GEN5_XPT`` is the only way to reach these, and the absence of a default is
    deliberate: they live in the lab's own folder on one machine, and a sibling-directory
    convention is what ``audit_reproducibility.py`` exists to reject. What is lost by
    skipping is coverage of the binary reader; what is not lost is any number in this
    repository, because ``data/plates`` holds the measurements as text.
    """
    directory = os.environ.get("YSTWIN_GEN5_XPT")
    if not directory:
        pytest.skip("set YSTWIN_GEN5_XPT to the directory holding the Gen5 .xpt files "
                    "to exercise the binary reader")
    path = pathlib.Path(directory).expanduser() / name
    if not path.exists():
        pytest.skip(f"Gen5 experiment file not present: {path}")
    return path


def _committed(name: str) -> tuple[list[str], list[float] | None, list[str], np.ndarray]:
    """One ``data/plates`` CSV -> ``(elapsed_hms, temperature_c, wells, values)``.

    ``temperature_c`` is ``None`` for the sheets exported without it -- the blank-subtracted
    ones -- which is how the second column is told from the first well.
    """
    path = paths.data_dir() / "plates" / name
    rows = [line.split(",") for line in path.read_text(encoding="utf-8").splitlines() if line]
    header = rows[0]
    has_temperature = header[1] == "temperature_c"
    first = 2 if has_temperature else 1
    wells = header[first:]
    times = [row[0] for row in rows[1:]]
    temperatures = [float(row[1]) for row in rows[1:]] if has_temperature else None
    values = np.array([[float(cell) for cell in row[first:]] for row in rows[1:]])
    return times, temperatures, wells, values


def _aligned(channel: gen5.Gen5Channel, wells: list[str], n_times: int) -> np.ndarray:
    """The channel's readings cut to the wells and timepoints one committed CSV holds."""
    index = {well: i for i, well in enumerate(channel.wells)}
    return channel.values[:n_times][:, [index[w] for w in wells]]


# ---------------------------------------------------------------------------
# decoding, with no wet-lab file in sight
# ---------------------------------------------------------------------------


class TestDecoding:
    """The parts that are pure arithmetic on bytes, and so run everywhere."""

    def test_wells_are_labelled_row_major_from_a1(self):
        """The archive stores a timepoint as A1..A12, B1..B12, ... H12, and so must we.

        Column-major here would put every reading in the wrong well while leaving the
        array the right shape, which is the one way to be wrong that no summary catches.
        """
        wells = gen5.well_names(8, 12)

        assert len(wells) == 96
        assert wells[:3] == ("A1", "A2", "A3")
        assert wells[11:14] == ("A12", "B1", "B2")
        assert wells[-1] == "H12"

    def test_a_plate_with_no_letter_for_its_last_row_is_refused(self):
        with pytest.raises(gen5.Gen5FormatError):
            gen5.well_names(27, 12)

    def test_names_are_latin1_because_the_degree_sign_is(self):
        """``T\N{DEGREE SIGN} OD600:600`` is 0xB0 in Latin-1 and invalid UTF-8.

        Decoding as UTF-8 raises rather than mis-reads, so a temperature series would be
        reported absent on every plate that has one.
        """
        buffer = b"\x0cT\xb0 OD600:600rest"

        text, after = gen5._cstring(buffer, 0)

        assert text == "T\N{DEGREE SIGN} OD600:600"
        assert buffer[after:] == b"rest"

    def test_the_protocol_path_is_reported_as_a_bare_file_name(self):
        """The archive stores a path on the instrument PC; only the file name is surfaced."""
        path = b"C:\\Users\\LAB\\Desktop\\Protocols\\4h-10min_mCitrine_OD600.prt"

        name = gen5._protocol_name(b"junk" + bytes([len(path)]) + path + b"junk")

        assert name == "4h-10min_mCitrine_OD600.prt"
        assert gen5._protocol_name(b"no protocol here") is None

    def test_a_file_that_is_not_a_compound_document_is_refused_by_name(self):
        with pytest.raises(gen5.Gen5FormatError, match="signature"):
            gen5._CompoundFile(b"PK\x03\x04" + b"\x00" * 1024)


class TestCompletedReads:
    """A block is sized by the protocol; only some of it may have happened."""

    def _channel(self, n_times: int, elapsed_ms) -> gen5.Gen5Channel:
        return gen5.Gen5Channel(
            name="OD600:600", fluorophore="OD600", optics="600", repeat=1,
            wells=("A1", "A2"), values=np.ones((n_times, 2)),
            elapsed_ms=None if elapsed_ms is None else np.asarray(elapsed_ms),
            temperature_c=None, started_at=None, rows=1, cols=2)

    def test_slots_a_stopped_run_never_filled_are_dropped(self):
        """22 rows of exact zeros are 22 slots, not 22 readings of OD 0.000."""
        values = np.vstack([np.full((3, 2), 0.5), np.zeros((2, 2))])
        elapsed = np.array([1000, 2000, 3000, 0, 0])
        temperature = np.array([30.0, 30.0, 30.1, 0.0, 0.0])

        cut_values, cut_elapsed, cut_temperature = gen5._completed_reads(
            "OD600:600", values, elapsed, temperature)

        assert cut_values.shape == (3, 2)
        assert cut_elapsed.tolist() == [1000, 2000, 3000]
        assert cut_temperature.tolist() == [30.0, 30.0, 30.1]

    def test_readings_and_timestamps_disagreeing_is_a_refusal(self):
        """Four rows carry a reading and three carry a time. Nothing here can choose."""
        values = np.vstack([np.full((4, 2), 0.5), np.zeros((1, 2))])
        elapsed = np.array([1000, 2000, 3000, 0, 0])

        with pytest.raises(gen5.Gen5FormatError, match="will not choose"):
            gen5._completed_reads("OD600:600", values, elapsed, None)

    def test_an_all_zero_table_means_no_time_axis_rather_than_no_reads(self):
        """Two 2026-06 files write a full-length table of zeros over readings that exist.

        A 49-point run did not happen at t=0 forty-nine times.
        """
        values = np.full((4, 2), 0.5)

        cut_values, cut_elapsed, _ = gen5._completed_reads(
            "RFP:530,580", values, np.zeros(4, dtype=np.int64), None)

        assert cut_values.shape == (4, 2)
        assert cut_elapsed is None

    def test_times_that_do_not_advance_are_refused(self):
        values = np.full((3, 2), 0.5)

        with pytest.raises(gen5.Gen5FormatError, match="strictly increasing"):
            gen5._completed_reads("OD600:600", values, np.array([1000, 1000, 3000]), None)

    def test_an_endpoint_read_refuses_an_elapsed_axis_rather_than_inventing_one(self):
        channel = self._channel(1, None)

        with pytest.raises(gen5.Gen5FormatError, match="endpoint"):
            channel.elapsed_hms()

    def test_elapsed_times_are_spelled_the_way_data_plates_spells_them(self):
        channel = self._channel(3, [498_000, 1_098_000, 14_898_000])

        assert channel.elapsed_hms() == ("0:08:18", "0:18:18", "4:08:18")
        assert channel.elapsed_hours()[0] == pytest.approx(0.13833333, abs=1e-7)


# ---------------------------------------------------------------------------
# the proof: does reading the .xpt give back the committed plate?
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestAgainstCommittedPlates:

    @pytest.mark.parametrize("xpt,channel_name,csv_name", PROVEN)
    def test_the_archive_reproduces_the_export(self, xpt, channel_name, csv_name):
        """Every reading the export carries, back out of the binary, to the export's precision.

        Fluorescence exactly; absorbance to the 0.0005 that three exported decimals can
        say about four recorded ones.
        """
        channel = gen5.read_xpt(_xpt(xpt)).channel(channel_name)
        times, temperatures, wells, expected = _committed(csv_name)
        n_times = len(times)

        got = _aligned(channel, wells, n_times)
        worst = float(np.max(np.abs(got - expected)))

        if channel.fluorophore == "OD600":
            assert worst <= EXPORT_HALF_STEP + FLOAT_SLACK, (
                f"{xpt} {channel_name}: worst disagreement {worst} exceeds half of the "
                "export's last decimal")
            assert worst > 0, (
                "an .xpt that agreed with a three-decimal export to the last bit would "
                "mean the fourth decimal is not in the file, and reading it would be "
                "pointless")
        else:
            assert worst == 0.0, f"{xpt} {channel_name}: fluorescence must match exactly"

    @pytest.mark.parametrize("xpt,channel_name,csv_name", PROVEN)
    def test_the_time_and_temperature_axes_come_back_too(self, xpt, channel_name, csv_name):
        """A right-looking matrix on a wrong time axis is still the wrong plate."""
        channel = gen5.read_xpt(_xpt(xpt)).channel(channel_name)
        times, temperatures, _, _ = _committed(csv_name)
        n_times = len(times)

        assert list(channel.elapsed_hms())[:n_times] == times
        assert temperatures is not None
        assert [round(float(t), 1) for t in channel.temperature_c[:n_times]] == temperatures

    @pytest.mark.parametrize("xpt,channel_name,csv_name", PROVEN)
    def test_the_archive_holds_every_well_the_export_dropped(self, xpt, channel_name, csv_name):
        """The exports carry 45-93 wells; the plate has 96, and all of them are in the file."""
        channel = gen5.read_xpt(_xpt(xpt)).channel(channel_name)
        _, _, wells, _ = _committed(csv_name)

        assert len(channel.wells) == 96
        assert set(wells) <= set(channel.wells)
        assert len(wells) < 96

    @pytest.mark.parametrize("xpt,channel_name,csv_name",
                             [row for row in PROVEN if row[1] == "OD600:600"])
    def test_absorbance_carries_the_fourth_decimal_the_export_dropped(
            self, xpt, channel_name, csv_name):
        """Every stored absorbance is an exact multiple of 0.0001, and most are not of 0.001.

        This is the gain: the export rounds to three decimals, and the reading it rounded
        is still here.
        """
        values = gen5.read_xpt(_xpt(xpt)).channel(channel_name).values

        assert np.allclose(values * 1e4, np.round(values * 1e4), atol=1e-6)
        finer_than_the_export = ~np.isclose(values * 1e3, np.round(values * 1e3), atol=1e-6)
        assert finer_than_the_export.mean() > 0.5


@pytest.mark.integration
class TestBlankSubtractedExports:
    """What the export lost, and whether the archive still has it."""

    def test_a_subtracted_sheet_is_the_archive_minus_its_blank_wells(self):
        """Replicate3 exported both halves, so the arithmetic can be checked before it is used.

        ``BLANK (OLD)OD`` and ``Blank mCitrine`` are the raw sheets with the mean of
        H1:H3 taken off at each timepoint. Reproducing them from the ``.xpt`` is what
        licenses reading the 20260804 export the same way, since that plate kept only its
        subtracted half.
        """
        read = gen5.read_xpt(_xpt("20260803_ER_Oxidative_Replicate3.xpt"))
        for channel_name, csv_name, tolerance in (
                ("OD600:600", "20260803_ERandoxidativestress_Replicate3__OD600_2.csv",
                 EXPORT_HALF_STEP + FLOAT_SLACK),
                ("mCitrine:480,530", "20260803_ERandoxidativestress_Replicate3__mCitrine_1.csv",
                 1 / 3 + FLOAT_SLACK)):
            channel = read.channel(channel_name)
            times, _, wells, expected = _committed(csv_name)
            index = {well: i for i, well in enumerate(channel.wells)}
            blank = channel.values[:, [index[w] for w in SUBTRACTED_BLANKS]].mean(axis=1)

            got = _aligned(channel, wells, len(times)) - blank[:len(times), None]
            worst = float(np.max(np.abs(got - expected)))

            assert worst <= tolerance, f"{channel_name}: {worst}"

    def test_plate_20260804_still_has_the_background_its_export_lost(self):
        """The whole point. H1-H3 read medium, not zero, in the file Gen5 wrote.

        ``data/plates/20260804_...`` has those wells within 0.002 of zero because the
        operator exported through a blank-subtraction transform. The archive is untouched:
        the media blanks carry real optical density and real autofluorescence, and both are
        far outside the noise the committed text shows.
        """
        read = gen5.read_xpt(_xpt(SUBTRACTED))
        od = read.channel("OD600:600")
        fluorescence = read.channel("mCitrine:480,530")
        index = {well: i for i, well in enumerate(od.wells)}
        columns = [index[w] for w in SUBTRACTED_BLANKS]

        raw_od = od.values[:, columns]
        raw_rfu = fluorescence.values[:, columns]

        assert raw_od.min() > 0.05, "media blanks should read around 0.09 OD, not zero"
        assert raw_rfu.min() > 100.0, "media blanks should read hundreds of RFU, not zero"
        _, _, committed_wells, committed = _committed(
            "20260804_ERandOxidativeStress_Replicate4__OD600.csv")
        subtracted = committed[:, [committed_wells.index(w) for w in SUBTRACTED_BLANKS]]
        assert np.abs(subtracted).max() < 0.005, "the committed export is the subtracted one"

    def test_the_subtracted_export_of_20260804_is_reproduced_from_the_archive(self):
        """Not just non-zero -- the right non-zero: raw minus mean(H1:H3) is the export.

        Fluorescence is the sharper half. Gen5 rounds the subtracted RFU to an integer,
        and the blank is the mean of three integers, so every disagreement must be 0 or
        exactly a third. Anything else would mean these are not the wells it subtracted.
        """
        read = gen5.read_xpt(_xpt(SUBTRACTED))
        for channel_name, csv_name, tolerance in (
                ("OD600:600", "20260804_ERandOxidativeStress_Replicate4__OD600.csv",
                 EXPORT_HALF_STEP + FLOAT_SLACK),
                ("mCitrine:480,530", "20260804_ERandOxidativeStress_Replicate4__mCitrine.csv",
                 1 / 3 + FLOAT_SLACK)):
            channel = read.channel(channel_name)
            times, _, wells, expected = _committed(csv_name)
            index = {well: i for i, well in enumerate(channel.wells)}
            blank = channel.values[:, [index[w] for w in SUBTRACTED_BLANKS]].mean(axis=1)

            residual = (_aligned(channel, wells, len(times)) - blank[:len(times), None]) - expected

            assert float(np.max(np.abs(residual))) <= tolerance, channel_name
            if channel.fluorophore != "OD600":
                thirds = np.abs(residual) * 3.0
                assert np.allclose(thirds, np.round(thirds), atol=1e-6)


@pytest.mark.integration
class TestTheWholeCorpus:
    """Every ``.xpt`` in the directory, not only the four that can be checked."""

    def test_every_experiment_file_parses_or_says_why(self):
        directory = _xpt(PROVEN[0][0]).parent
        files = sorted(directory.glob("*.xpt"))
        assert files, f"no .xpt files under {directory}"

        for path in files:
            read = gen5.read_xpt(path)

            assert read.channels, f"{path.name}: no channels"
            for channel in read.channels:
                assert channel.values.shape == (channel.n_times, len(channel.wells))
                assert np.isfinite(channel.values).all(), f"{path.name} {channel.name}"
                if channel.elapsed_ms is not None:
                    assert len(channel.elapsed_ms) == channel.n_times
                    assert np.all(np.diff(channel.elapsed_ms) > 0)
                if channel.temperature_c is not None:
                    assert len(channel.temperature_c) == channel.n_times
                if channel.started_at is not None:
                    assert channel.started_at.tzinfo is _dt.UTC

    def test_the_container_reader_agrees_with_olefile(self):
        """The hand-written MS-CFB reader, against the library it replaces.

        ``olefile`` is not a dependency of this project, and adding one to open a
        container whose whole grammar is in the specification would be the wrong trade.
        This is how that decision is kept honest: where the library is installed, every
        stream of every file must come out byte-identical.
        """
        olefile = pytest.importorskip("olefile", reason="cross-check needs olefile")
        directory = _xpt(PROVEN[0][0]).parent

        for path in sorted(directory.glob("*.xpt")):
            mine = gen5._CompoundFile(path.read_bytes())
            with olefile.OleFileIO(str(path)) as reference:
                theirs = {"/".join(entry) for entry in reference.listdir(streams=True)}

                assert set(mine.paths) == theirs, path.name
                for stream in sorted(theirs):
                    assert mine.read(stream) == reference.openstream(stream).read(), \
                        f"{path.name}:{stream}"


def _channel_with(values):
    """A minimal channel around a values array, for the guards that need no file."""
    import numpy as np

    return gen5.Gen5Channel(
        name="test:480,530", fluorophore="test", optics="480,530", repeat=1,
        wells=tuple(f"A{i + 1}" for i in range(len(values[0]))),
        values=np.asarray(values, dtype=float),
        elapsed_ms=None, temperature_c=None, started_at=None,
        rows=1, cols=len(values[0]),
    )

class TestTheOverflowFlagCannotReachAnAverageBySilence:
    """`OVERFLOW` is what Gen5 writes above the top of scale, and it is not a reading.

    11,149 of the 303,552 fluorescence readings in the 22-file corpus are exactly it. No
    caller in this repository reads a channel that contains one, which is precisely why the
    guard is worth having: the failure it prevents is a future one, and it is silent. A
    -99999 averaged into a reporter trace is not a wrong number, it is a number with no
    relationship to the plate.
    """

    def test_the_sentinel_is_defined_once(self):
        """`calib/gain_linearity.py` measured it; `plate/gen5.py` owns it, because it is a
        property of the file format and not of the analysis."""
        from ystwin.calib.gain_linearity import RFU_OVERFLOW
        from ystwin.plate.gen5 import OVERFLOW

        assert RFU_OVERFLOW is OVERFLOW or RFU_OVERFLOW == OVERFLOW

    def test_a_clean_channel_is_unaffected(self):
        channel = _channel_with([[1.0, 2.0], [3.0, 4.0]])

        assert not channel.overflowed().any()
        assert channel.frame().to_numpy().tolist() == [[1.0, 2.0], [3.0, 4.0]]

    def test_a_channel_carrying_the_flag_refuses_by_default(self):
        channel = _channel_with([[1.0, gen5.OVERFLOW], [3.0, 4.0]])

        with pytest.raises(gen5.Gen5FormatError, match="overflow"):
            channel.frame()

    def test_the_refusal_says_how_many_and_what_to_do(self):
        channel = _channel_with([[1.0, gen5.OVERFLOW], [gen5.OVERFLOW, 4.0]])

        with pytest.raises(gen5.Gen5FormatError) as raised:
            channel.frame()

        assert "2 of 4" in str(raised.value)
        assert "on_overflow" in str(raised.value)

    def test_nan_masks_it(self):
        import math

        channel = _channel_with([[1.0, gen5.OVERFLOW], [3.0, 4.0]])

        values = channel.frame(on_overflow="nan").to_numpy()

        assert math.isnan(values[0][1])
        assert values[0][0] == 1.0

    def test_keep_returns_what_the_file_holds(self):
        channel = _channel_with([[1.0, gen5.OVERFLOW], [3.0, 4.0]])

        assert channel.frame(on_overflow="keep").to_numpy()[0][1] == gen5.OVERFLOW

    def test_an_unknown_policy_is_refused(self):
        channel = _channel_with([[1.0, 2.0], [3.0, 4.0]])

        with pytest.raises(ValueError, match="on_overflow"):
            channel.frame(on_overflow="ignore")


def _read_with_overflow(extra=2.0, all_flagged=False):
    """A two-timepoint kinetic read whose channel carries the over-range flag."""
    import numpy as np

    values = ([[gen5.OVERFLOW, gen5.OVERFLOW], [gen5.OVERFLOW, gen5.OVERFLOW]]
              if all_flagged else [[1.0, gen5.OVERFLOW], [3.0, extra]])
    channel = gen5.Gen5Channel(
        name="mCitrine:480,530[2]", fluorophore="mCitrine", optics="480,530", repeat=2,
        wells=("A1", "A2"), values=np.asarray(values, dtype=float),
        elapsed_ms=np.asarray([0, 600_000]), temperature_c=None, started_at=None,
        rows=1, cols=2,
    )
    return gen5.Gen5Read(source=pathlib.Path("synthetic.xpt"), protocol=None,
                         channels=(channel,))

class TestTheSentinelDoesNotReachASynergyRun:
    """The guard on `frame()` protected a door the export pipeline does not use.

    `as_synergy_run` must mask the format's -99999 overflow flag without masking finite
    negative readings or using either one's sign to classify processing.
    """

    def test_a_flagged_channel_becomes_nan_not_minus_ninety_nine_thousand(self):
        run = gen5.as_synergy_run(_read_with_overflow())
        values = run.blocks[0].data.to_numpy()

        assert not (values <= gen5.OVERFLOW + 1.0).any()
        assert np.isnan(values).sum() == 1

    def test_a_flagged_raw_channel_is_not_called_blank_subtracted(self):
        """Overflow affects measurement availability, not the raw archive's processing state."""
        run = gen5.as_synergy_run(_read_with_overflow())

        assert run.blocks[0].blank_subtracted is False

    def test_a_finite_negative_archive_reading_stays_raw_and_is_preserved(self):
        run = gen5.as_synergy_run(_read_with_overflow(extra=-12.5))

        assert run.blocks[0].blank_subtracted is False
        assert run.raw_channel("mCitrine").data.iloc[1, 1] == -12.5
        assert np.isnan(run.blocks[0].data.iloc[0, 1])

    def test_a_channel_that_is_entirely_flagged_claims_nothing(self):
        run = gen5.as_synergy_run(_read_with_overflow(all_flagged=True))

        assert run.blocks[0].blank_subtracted is False
