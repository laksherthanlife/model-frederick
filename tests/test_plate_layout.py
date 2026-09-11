"""The plate map, recovered by matching numbers rather than assumed from row order.

The team's per-construct sheets hold (RFU - blank)/(OD - blank) per dose. Searching
every construct-to-column-block assignment for the one that reproduces those values
from the raw wells recovers the map, and reports how far ahead of the runner-up it
is -- so a weak recovery is visible rather than silent.
"""

import numpy as np
import pandas as pd
import pytest

from ystwin.plate.layout import (
    NEWPROTOCOL_LAYOUT,
    PlateLayout,
    recover_layout,
    wells_for,
)


def test_the_recorded_layout_covers_all_four_constructs():
    assert set(NEWPROTOCOL_LAYOUT.construct_columns) == {
        "UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"
    }


def test_each_construct_occupies_three_adjacent_columns():
    for cols in NEWPROTOCOL_LAYOUT.construct_columns.values():
        assert len(cols) == 3
        assert list(cols) == list(range(min(cols), min(cols) + 3))


def test_no_two_constructs_share_a_column():
    used = [c for cols in NEWPROTOCOL_LAYOUT.construct_columns.values() for c in cols]

    assert len(used) == len(set(used))


def test_wells_for_a_condition_are_the_row_crossed_with_the_block():
    wells = wells_for("UPRE1", dose_index=0, layout=NEWPROTOCOL_LAYOUT)

    assert wells == ["A1", "A2", "A3"]


def test_later_doses_sit_in_later_rows():
    assert wells_for("AlteredYap1", 3, NEWPROTOCOL_LAYOUT) == ["D10", "D11", "D12"]


def test_a_dose_index_past_the_plate_is_refused():
    with pytest.raises(IndexError, match="dose index"):
        wells_for("UPRE1", 99, NEWPROTOCOL_LAYOUT)


def test_an_unknown_construct_is_refused_by_name():
    with pytest.raises(KeyError, match="mScarletSensor"):
        wells_for("mScarletSensor", 0, NEWPROTOCOL_LAYOUT)


class TestRecovery:
    def _synthetic(self):
        """Build wells from a known layout, then check recovery finds it back."""
        rng = np.random.default_rng(0)
        t = np.linspace(0, 4, 25)
        layout = PlateLayout(
            construct_columns={"UPRE1": (1, 2, 3), "UPRE2": (4, 5, 6)},
            dose_rows=tuple("ABC"),
        )
        wells, derived = {}, []
        for ci, (con, cols) in enumerate(layout.construct_columns.items()):
            for di, row in enumerate(layout.dose_rows):
                trace = 1000 * (1 + ci) + 300 * di * t
                for c in cols:
                    wells[f"{row}{c}"] = trace + rng.normal(0, 5, t.size)
                for tt, v in zip(t, trace):
                    derived.append({"construct": con, "dose_mM": float(di) / 10,
                                    "time_h": tt, "signal": v})
        return pd.DataFrame(wells, index=pd.Index(t, name="time_h")), pd.DataFrame(derived), layout

    def test_it_recovers_the_layout_that_generated_the_wells(self):
        normalised, derived, truth = self._synthetic()

        found = recover_layout(normalised, derived)

        assert found.layout.construct_columns == truth.construct_columns

    def test_it_reports_the_fit_error_and_the_margin_over_the_runner_up(self):
        normalised, derived, _ = self._synthetic()

        found = recover_layout(normalised, derived)

        assert found.median_relative_error < 0.05
        assert found.margin_over_runner_up > 1.0

    def test_indistinguishable_constructs_yield_a_small_margin(self):
        """Two constructs with identical traces cannot be told apart, and it shows."""
        normalised, derived, _ = self._synthetic()
        derived.loc[derived.construct == "UPRE2", "signal"] = derived.loc[
            derived.construct == "UPRE1", "signal"
        ].to_numpy()
        for row in "ABC":
            for c in (4, 5, 6):
                normalised[f"{row}{c}"] = normalised[f"{row}{c - 3}"]

        found = recover_layout(normalised, derived)

        assert found.margin_over_runner_up == pytest.approx(1.0, abs=0.15)


class TestRecordedPlateMaps:
    """Layouts transcribed from the wet-lab logbook, which is the authority.

    The numerical recovery agrees with it on construct columns and dose rows, which
    is the check that matters. Blank positions it could not have known: they moved
    from H1-H3 on 2026-07-22 to H4-H6 on every later plate, and on one plate the
    detector picked an empty well instead.
    """

    def test_every_recorded_plate_has_the_same_construct_blocks(self):
        from ystwin.plate.layout import RECORDED_PLATES

        for name, plate in RECORDED_PLATES.items():
            cols = {c: v[0] for c, v in plate.layout.construct_columns.items()}
            assert cols == {"UPRE1": 1, "UPRE2": 4, "NativeYap1": 7, "AlteredYap1": 10}, name

    def test_the_recovered_map_agrees_with_the_logbook(self):
        from ystwin.plate.layout import NEWPROTOCOL_LAYOUT, RECORDED_PLATES

        recorded = RECORDED_PLATES["20260722"].layout

        assert recorded.construct_columns == NEWPROTOCOL_LAYOUT.construct_columns
        assert recorded.dose_rows == NEWPROTOCOL_LAYOUT.dose_rows

    def test_the_first_plate_blanks_the_first_three_wells_of_row_h(self):
        from ystwin.plate.layout import RECORDED_PLATES

        assert RECORDED_PLATES["20260722"].blank_wells == ("H1", "H2", "H3")

    def test_one_plate_moved_the_blanks_along_row_h(self):
        """20260728 is the only one that really did, and it is corroborated twice over:
        detection finds H4-H6 in its export, and in its ``.xpt`` H1-H3 are cultures --
        they grow from 0.117 to 0.642 OD while a blank does not grow at all."""
        from ystwin.plate.layout import RECORDED_PLATES

        assert RECORDED_PLATES["20260728"].blank_wells == ("H4", "H5", "H6")

    def test_the_third_replicate_is_corrected_against_its_export(self):
        """Its logbook entry says H4-H6, but no such wells exist in the export and
        its H1-H3 are flat in the reporter channel. The entry looks templated."""
        from ystwin.plate.layout import RECORDED_PLATES

        plate = RECORDED_PLATES["20260803"]
        assert plate.blank_wells == ("H1", "H2", "H3")
        assert "Logbook records H4-H6" in plate.note

    def test_the_fourth_replicate_is_corrected_against_its_instrument_file(self):
        """The same templated H4-H6 entry, and here it mattered rather than being moot.

        20260803's export has no H4-H6 wells at all, so the stale entry fell through to
        detection and the right blank was used anyway. 20260804 is committed from its
        ``.xpt``, which carries all 96 wells including H4-H6 -- so the entry would have
        been *used*, and it would have subtracted 56 RFU of instrument dark count where
        341 RFU of medium autofluorescence was right.
        """
        from ystwin.plate.layout import RECORDED_PLATES

        plate = RECORDED_PLATES["20260804"]
        assert plate.blank_wells == ("H1", "H2", "H3")
        assert "instrument file settles it" in plate.note

    def test_dtt_runs_to_five_millimolar_and_peroxide_to_four(self):
        """The logbook's stock prep tops out at 4 mM H2O2 but its map says 5 mM."""
        from ystwin.plate.layout import RECORDED_PLATES

        plate = RECORDED_PLATES["20260722"]
        assert plate.doses_mM["DTT"][-1] == 5.0
        assert plate.doses_mM["H2O2"][-1] == 4.0

    def test_a_plate_can_report_the_blank_wells_for_a_given_export(self):
        from ystwin.plate.layout import blanks_for_export

        assert blanks_for_export("20260803_ER&oxidativestress_Replicate3.xlsx") == ("H1", "H2", "H3")
        assert blanks_for_export("20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx") == (
            "H1", "H2", "H3"
        )

    def test_an_unrecognised_export_returns_no_recorded_blanks(self):
        from ystwin.plate.layout import blanks_for_export

        assert blanks_for_export("something_else.xlsx") is None


def test_recorded_well_roles_distinguish_stressed_cultures_from_media():
    from ystwin.plate.layout import RECORDED_PLATES

    recorded = RECORDED_PLATES["20260722"]
    assert recorded.well_roles["G1"] == "culture"
    assert recorded.well_roles["H1"] == "blank"
    assert "H12" not in recorded.well_roles
    assert set(recorded.culture_wells) == {
        f"{row}{column}" for row in "ABCDEFG" for column in range(1, 13)
    }


def test_a_recorded_well_cannot_be_both_culture_and_blank():
    from ystwin.plate.layout import RecordedPlate

    with pytest.raises(ValueError, match="A1"):
        RecordedPlate(PlateLayout({"UPRE1": (1,)}, ("A",)), ("A1",), {"DTT": (0.0,)})


def test_explicit_plate_identity_survives_raw_channel_and_alignment_on_a_renamed_export():
    from pathlib import Path

    from ystwin.plate.layout import RECORDED_PLATES
    from ystwin.plate.synergy import KineticBlock, SynergyRun

    recorded = RECORDED_PLATES["20260722"]
    t = np.linspace(0.0, 4.0, 25)
    od = pd.DataFrame({"G1": np.full_like(t, 0.08), "H1": np.full_like(t, 0.09)}, index=t)
    block = KineticBlock("OD600", "OD600", "600", "OD", False, False, od, [])
    run = SynergyRun(Path("renamed.xlsx"), (block,), recorded_plate=recorded)

    assert run.raw_channel("OD600").data.attrs["well_roles"] == recorded.well_roles
    assert run.aligned("OD600").attrs["well_roles"] == recorded.well_roles
    assert run.aligned("OD600")["OD600"].attrs["well_roles"] == recorded.well_roles
    assert "well_roles" not in od.attrs


def test_committed_replay_carries_verified_recorded_roles():
    from ystwin.plate import replay
    from ystwin.plate.layout import RECORDED_PLATES

    run = replay.load_run("20260804_ER_Oxidative_Replicate4.xpt")
    recorded = RECORDED_PLATES["20260804"]

    assert run.recorded_plate == recorded
    assert run.raw_channel("OD600").data.attrs["well_roles"]["G1"] == "culture"
    assert run.raw_channel("mCitrine").data.attrs["well_roles"]["H1"] == "blank"


def test_public_plate_prepare_retains_a_low_od_well_recorded_as_a_culture(monkeypatch):
    import importlib.util
    from pathlib import Path

    from ystwin.plate.layout import RECORDED_PLATES
    from ystwin.plate.synergy import KineticBlock, SynergyRun

    root = Path(__file__).resolve().parents[1]
    loader = importlib.util.spec_from_file_location(
        "sensor_correctness", root / "scripts" / "run_sensor_characterisation.py")
    module = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(module)
    path = Path("20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx")
    t = np.linspace(0.0, 4.0, 25)
    od = pd.DataFrame({"G1": np.full_like(t, 0.105), "A1": 0.2 * np.exp(0.2 * t),
                       **{well: np.full_like(t, 0.09) for well in ("H1", "H2", "H3")}}, index=t)
    rfu = pd.DataFrame({"G1": np.full_like(t, 700.0), "A1": 1000.0 * np.exp(0.2 * t),
                        **{well: np.full_like(t, 300.0) for well in ("H1", "H2", "H3")}}, index=t)
    blocks = tuple(KineticBlock(channel, channel, optics, channel, False, False, frame, [])
                   for channel, optics, frame in (("OD600", "600", od),
                                                  ("mCitrine", "480,530", rfu)))
    monkeypatch.setattr(module, "read_synergy_kinetic", lambda _: SynergyRun(path, blocks))

    prepared = module.prepare(path, recorded_blanks=RECORDED_PLATES["20260722"].blank_wells)

    assert prepared is not None
    assert "G1" in prepared[0].columns
    assert "G1" in prepared[1].columns
    assert prepared[4] == pytest.approx(0.09)
    assert prepared[5] == pytest.approx(300.0)
    with pytest.raises(ValueError, match="recorded.*culture"):
        module.prepare(path, recorded_blanks=("G1",))
    assert module.prepare(path, recorded_blanks=()) is None
    missing = tuple(KineticBlock(channel, channel, optics, channel, False, False,
                                 frame[["G1", "A1"]], [])
                    for channel, optics, frame in (("OD600", "600", od),
                                                   ("mCitrine", "480,530", rfu)))
    monkeypatch.setattr(module, "read_synergy_kinetic", lambda _: SynergyRun(path, missing))
    assert module.prepare(path, recorded_blanks=("H1", "H2", "H3")) is None
    assert module.prepare(path) is None


def test_workbook_reader_accepts_explicit_recorded_identity_without_using_the_filename(monkeypatch, tmp_path):
    from ystwin.plate.layout import RECORDED_PLATES
    from ystwin.plate.synergy import read_synergy_kinetic

    sheet = pd.DataFrame([["Time", "T° OD600:600", "G1", "H1"],
                          ["00:00:00", 30.0, 0.08, 0.09],
                          ["00:10:00", 30.0, 0.08, 0.09]])
    monkeypatch.setattr(pd, "read_excel", lambda *a, **k: {"OD600": sheet})
    recorded = RECORDED_PLATES["20260722"]

    path = tmp_path / "renamed.xlsx"
    path.write_bytes(b"synthetic workbook supplied by the read_excel stub")
    run = read_synergy_kinetic(path, recorded_plate=recorded, processing_states={"OD600": False})

    assert run.recorded_plate is recorded
    assert run.raw_channel("OD600").data.attrs["well_roles"]["G1"] == "culture"


def test_replay_dispatch_preserves_explicit_plate_metadata_and_reader_options():
    from ystwin.plate import replay
    from ystwin.plate.layout import RecordedPlate

    namespace = {"read_synergy_kinetic": None}
    replay.install(namespace)
    recorded = RecordedPlate(PlateLayout({"custom": (1,)}, ("A",)), ("H1",), {"DTT": (0.0,)})
    export = "20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx"

    run = namespace["read_synergy_kinetic"](export, True, recorded_plate=recorded)

    assert run.recorded_plate is recorded
    assert run.aligned("OD600[1]").attrs["well_roles"] == {"A1": "culture", "H1": "blank"}
    with pytest.raises(ValueError, match="prefer_raw=False"):
        namespace["read_synergy_kinetic"](export, False)


def test_replay_refuses_conflicting_dose_files_instead_of_picking_one():
    from ystwin.plate import replay

    manifest = replay.load_manifest()
    export = manifest.export.iloc[0]
    original = manifest.dose_response_file.iloc[0]
    alternative = next(name for name in manifest.dose_response_file if name and name != original)
    manifest.loc[manifest.index[0], "dose_response_file"] = alternative

    with pytest.raises(ValueError, match="conflicting"):
        replay.load_doses(export, manifest=manifest)
