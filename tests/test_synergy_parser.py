import numpy as np
import pytest

from ystwin.plate.synergy import read_synergy_kinetic


def test_parses_elapsed_clock_strings_into_hours(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file())

    od = run.channel("OD600")
    assert od.times_h == pytest.approx([7.417 / 60, 17.417 / 60, 67.417 / 60], abs=1e-3)


def test_wells_become_columns_with_values_preserved(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file())

    od = run.channel("OD600")
    assert list(od.data.columns) == ["A1", "A2", "B1"]
    assert od.data["A2"].to_numpy() == pytest.approx([0.20, 0.24, 0.40])


def test_discovers_both_channels_and_names_them_from_the_header(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file())

    assert set(run.channel_names) == {"OD600", "mCitrine"}
    assert run.channel("mCitrine").data["B1"].to_numpy() == pytest.approx([3000, 3300, 5000])


def test_header_is_found_even_when_the_export_has_a_blank_preamble(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(lead_blank_rows=2))

    assert list(run.channel("OD600").data.columns) == ["A1", "A2", "B1"]
    assert run.channel("OD600").data.iloc[0, 0] == pytest.approx(0.10)


def test_temperature_trace_is_kept_alongside_the_measurements(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file())

    assert run.channel("OD600").temperature_c == pytest.approx([30.0, 30.1, 30.0])


def test_channels_are_time_aligned_onto_a_common_grid(make_synergy_file):
    """OD and fluorescence are read ~50 s apart; the twin needs them on one clock."""
    run = read_synergy_kinetic(make_synergy_file())

    aligned = run.aligned(reference="OD600")
    assert list(aligned.columns.names) == ["channel", "well"]
    assert aligned[("mCitrine", "A1")].to_numpy() == pytest.approx([1000, 1200, 2000], rel=0.02)
    assert np.all(np.diff(aligned.index.to_numpy()) > 0)


def test_trailing_padding_rows_are_not_mistaken_for_measurements(make_synergy_file):
    """Exports pad to the protocol length with 00:00:00 and blank wells."""
    run = read_synergy_kinetic(make_synergy_file(pad_rows=5))

    od = run.channel("OD600")
    assert len(od.times_h) == 3
    assert not od.data.isna().to_numpy().any()
    assert np.all(np.diff(np.asarray(od.times_h)) > 0)


def test_a_legitimate_zero_timestamp_first_reading_is_still_kept(make_synergy_file):
    """Guard against 'time == 0 means padding' — a run may genuinely start at 00:00:00."""
    path = make_synergy_file(od_times=("00:00:00", "00:17:25", "01:07:25"), pad_rows=4)

    od = read_synergy_kinetic(path).channel("OD600")
    assert od.times_h[0] == pytest.approx(0.0)
    assert len(od.times_h) == 3
