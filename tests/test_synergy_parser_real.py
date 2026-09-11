"""The parser must survive the actual exports, not just the synthetic fixture."""
import numpy as np
import pytest

from ystwin.plate.synergy import WELL_RE, read_synergy_kinetic

pytestmark = pytest.mark.integration


def test_real_export_yields_paired_od_and_reporter_channels(real_er_prelim):
    run = read_synergy_kinetic(real_er_prelim)

    assert set(run.channel_names) == {"OD600", "mCitrine"}


def test_real_export_well_labels_are_all_valid_plate_positions(real_er_prelim):
    run = read_synergy_kinetic(real_er_prelim)

    for block in run.blocks:
        bad = [w for w in block.wells if not WELL_RE.match(w)]
        assert bad == [], f"{block.channel} has non-well columns {bad}"


def test_real_export_time_is_strictly_increasing_in_every_channel(real_er_prelim):
    run = read_synergy_kinetic(real_er_prelim)

    for block in run.blocks:
        t = np.asarray(block.times_h)
        assert np.all(np.diff(t) > 0), f"{block.channel} time not monotonic"


def test_real_od_values_are_in_a_physically_plausible_range(real_er_prelim):
    od = read_synergy_kinetic(real_er_prelim).channel("OD600").data

    assert np.nanmin(od.to_numpy()) > 0.0
    assert np.nanmax(od.to_numpy()) < 4.0


def test_aligning_real_channels_loses_no_wells_and_introduces_no_gaps(real_er_prelim):
    run = read_synergy_kinetic(real_er_prelim)

    aligned = run.aligned(reference="OD600")
    od_wells = set(run.channel("OD600").wells)
    fl_wells = set(run.channel("mCitrine").wells)
    assert set(aligned["OD600"].columns) == od_wells
    assert set(aligned["mCitrine"].columns) == fl_wells
    assert not aligned.isna().to_numpy().any()
