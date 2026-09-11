import numpy as np
import pandas as pd
import pytest

from ystwin.diagnostics.dilution_confound import plate_dilution_report
from ystwin.reporter import ReporterKinetics, simulate_reporter


def _plate(n_growers=4, n_blanks=2, n_stalled=1, seed=0):
    """Synthetic plate: growing cultures, media blanks, and one stalled well."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 24, 300)
    cols, wells = {}, []
    kin = ReporterKinetics(k_deg=0.0)
    for i in range(n_growers):
        w = f"A{i + 1}"
        mu = 0.30 + 0.02 * i
        od = 0.05 * np.exp(mu * t) / (1 + 0.05 * (np.exp(mu * t) - 1) / 1.5)
        mu_t = np.clip(np.gradient(np.log(od), t), 1e-3, None)
        per_cell = simulate_reporter(t, np.ones_like(t), mu_t, kin, r0=1 / mu_t[0])
        cols[("OD600", w)] = od + 0.09
        cols[("mCitrine", w)] = per_cell * od * 500 + 250
        wells.append(w)
    for i in range(n_blanks):
        w = f"H{11 + i}"
        cols[("OD600", w)] = 0.09 + rng.normal(0, 5e-4, t.size)
        cols[("mCitrine", w)] = 250 + rng.normal(0, 3, t.size)
    for i in range(n_stalled):
        w = f"G{i + 1}"
        cols[("OD600", w)] = 0.09 + 0.30 + rng.normal(0, 5e-4, t.size)
        cols[("mCitrine", w)] = 900 + rng.normal(0, 5, t.size)
    frame = pd.DataFrame(cols, index=pd.Index(t, name="time_h"))
    frame.columns = pd.MultiIndex.from_tuples(frame.columns, names=["channel", "well"])
    return frame, wells


def test_media_blanks_are_detected_and_used_as_the_correction():
    frame, _ = _plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert set(rep.blank_wells) == {"H11", "H12"}
    assert rep.od_blank == pytest.approx(0.09, abs=0.005)
    assert rep.rfu_background == pytest.approx(250, abs=5)


def test_one_row_per_analysed_culture_and_blanks_are_not_cultures():
    frame, growers = _plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert set(rep.per_well.well) == set(growers)


def test_excluded_wells_are_recorded_with_the_reason_not_dropped_silently():
    frame, _ = _plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert "G1" in set(rep.excluded.well)
    assert rep.excluded.reason.str.contains("did not grow").all()


def test_a_plate_of_pure_dilution_artefacts_is_summarised_as_such():
    frame, _ = _plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert rep.median_dilution_r2 > 0.8
    assert (rep.per_well.verdict == "dilution-dominated").mean() > 0.5


def test_requires_the_named_channels_to_exist():
    frame, _ = _plate()

    with pytest.raises(KeyError, match="mScarlet"):
        plate_dilution_report(frame, od_channel="OD600", reporter_channel="mScarlet")


def test_refuses_to_guess_a_blank_when_no_well_looks_like_media():
    frame, _ = _plate(n_blanks=0)

    with pytest.raises(ValueError, match="no media-blank"):
        plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")


def test_explicit_blank_wells_override_detection():
    frame, _ = _plate(n_blanks=0)

    rep = plate_dilution_report(
        frame, od_channel="OD600", reporter_channel="mCitrine",
        blank_wells=["G1"],
    )

    assert rep.blank_wells == ["G1"]


def _mismatched_plate():
    """OD read across the plate, fluorescence on a subset -- as in the 2026-07 replicates."""
    frame, growers = _plate(n_growers=4, n_blanks=2)
    od_only = frame[[c for c in frame.columns if c[0] == "OD600"]].copy()
    extra = od_only.rename(columns={"A1": "D1"}, level=1)[[("OD600", "D1")]]
    frame = pd.concat([frame, extra], axis=1)
    frame.columns = pd.MultiIndex.from_tuples(frame.columns, names=["channel", "well"])
    return frame, growers


def test_wells_missing_from_the_reporter_channel_are_excluded_not_crashed_on():
    frame, growers = _mismatched_plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert set(rep.per_well.well) == set(growers)


def test_the_excluded_wells_are_named_with_the_channel_they_were_missing_from():
    frame, _ = _mismatched_plate()

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    row = rep.excluded.set_index("well").loc["D1"]
    assert "mCitrine" in row.reason


def test_a_blank_missing_from_one_channel_is_dropped_from_the_correction():
    frame, _ = _mismatched_plate()
    frame = frame.drop(columns=[("mCitrine", "H12")])

    rep = plate_dilution_report(frame, od_channel="OD600", reporter_channel="mCitrine")

    assert "H12" not in rep.blank_wells
    assert rep.blank_wells == ["H11"]


def test_no_well_shared_between_the_channels_is_refused_by_name():
    frame, _ = _plate()
    od_only = frame[[c for c in frame.columns if c[0] == "OD600"]]
    renamed = frame[[c for c in frame.columns if c[0] == "mCitrine"]].rename(
        columns=lambda w: f"Z{w[1:]}", level=1
    )
    merged = pd.concat([od_only, renamed], axis=1)
    merged.columns = pd.MultiIndex.from_tuples(merged.columns, names=["channel", "well"])

    with pytest.raises(ValueError, match="no well appears in both"):
        plate_dilution_report(merged, od_channel="OD600", reporter_channel="mCitrine")
