"""The human-readable surfaces are part of the deliverable and must not crash."""

import numpy as np
import pandas as pd
import pytest

from ystwin.diagnostics.dilution_confound import detect_blank_wells, dilution_confound_report
from ystwin.gates.g1_optical import assess_well
from ystwin.reporter import ReporterKinetics
from ystwin.readings import RawOD, RawRFU


def _trace():
    t = np.linspace(0, 20, 400)
    od = 0.05 * np.exp(0.25 * t)
    return t, od, od * 4000.0


def test_blank_detection_picks_wells_pinned_at_the_plate_floor():
    t = np.linspace(0, 10, 50)
    od = pd.DataFrame({
        "A1": 0.05 * np.exp(0.3 * t) + 0.09,
        "H11": np.full(50, 0.091),
        "H12": np.full(50, 0.089),
    })

    assert detect_blank_wells(RawOD(od)) == ["H11", "H12"]


def test_blank_detection_tolerance_is_adjustable():
    od = pd.DataFrame({"A1": np.full(50, 0.30), "H12": np.full(50, 0.09)})

    assert detect_blank_wells(RawOD(od), tolerance=0.5) == ["A1", "H12"]


def test_dilution_report_summary_states_both_folds_and_the_verdict():
    t, od, rfu = _trace()

    text = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics()).summary()

    assert "naive RFU/OD fold" in text
    assert "implied min k_deg" in text


def test_well_quality_summary_names_the_failing_checks():
    t = np.linspace(0, 18, 110)
    flat = np.full_like(t, 0.5)

    text = assess_well(t, RawOD(flat), od_blank=0.09).summary()

    assert text.startswith("well: FAIL")
    assert "dynamic_range" in text


def test_well_quality_summary_is_terse_when_everything_passes():
    t = np.linspace(0, 18, 110)
    od = 0.08 * np.exp(0.2 * t) / (1 + 0.08 * (np.exp(0.2 * t) - 1) / 0.8) + 0.09

    assert assess_well(t, RawOD(od), od_blank=0.09).summary() == "well: PASS"


def test_flux_range_summary_reports_the_interval_and_relative_width():
    from ystwin.fba.fva import FluxRange

    text = FluxRange("DM_x", 0.0, 0.02, 0.38, 0.34, 10.0, None).summary()

    assert "relative width 1.000" in text


def test_physiology_report_summary_lists_each_flux_against_its_reference():
    from ystwin.fba.physiology import REFERENCE_AEROBIC_BATCH, PhysiologyReport

    report = PhysiologyReport(
        reference=REFERENCE_AEROBIC_BATCH,
        observed={"growth_rate": 0.38}, expected={"growth_rate": 0.40},
        tolerance=0.35, relative_error={"growth_rate": 0.05},
    )

    assert "PASS" in report.summary()
    assert report.max_relative_error == pytest.approx(0.05)


class TestBlankDetectionRequiresNoGrowth:
    """Low OD alone is not enough: a barely-growing culture at a high stressor dose
    sits close to the floor too. A media blank has no biomass, so the definitive
    test is that it does not grow. Including one poorly-growing well inflated the
    blank by 0.014 OD on a real plate, which pushed every derived growth rate above
    the physiological maximum for yeast."""

    def _plate(self, floor=0.095):
        t = np.linspace(0, 4.14, 25)
        rng = np.random.default_rng(0)
        return pd.DataFrame({
            "H1": floor + rng.normal(0, 5e-4, t.size),
            "H2": floor + 0.005 + rng.normal(0, 5e-4, t.size),
            "F9": np.linspace(0.100, 0.112, t.size),      # near the floor but grows
            "A1": 0.12 * np.exp(0.30 * t),
        }, index=pd.Index(t, name="time_h"))

    def test_a_flat_well_at_the_floor_is_a_blank(self):
        assert set(detect_blank_wells(RawOD(self._plate()))) >= {"H1", "H2"}

    def test_a_slowly_growing_well_near_the_floor_is_not_a_blank(self):
        assert "F9" not in detect_blank_wells(RawOD(self._plate()))

    def test_a_clearly_growing_culture_is_not_a_blank(self):
        assert "A1" not in detect_blank_wells(RawOD(self._plate()))

    def test_the_growth_tolerance_is_adjustable(self):
        plate = self._plate()

        strict = detect_blank_wells(RawOD(plate), max_growth_fold=1.005)
        loose = detect_blank_wells(RawOD(plate), max_growth_fold=1.30)

        assert "F9" not in strict
        assert "F9" in loose

    def test_including_a_growing_well_would_raise_the_blank(self):
        """Guards the specific failure: folding a slow grower into the blank."""
        plate = self._plate()
        correct = float(plate[detect_blank_wells(RawOD(plate))].to_numpy().mean())
        contaminated = float(plate[detect_blank_wells(RawOD(plate)) + ["F9"]].to_numpy().mean())

        assert contaminated > correct


class TestBlankDetectionUsesBothChannels:
    """OD alone mistakes a drifting medium well for a culture.

    On the 2026-08-03 plate two medium wells drifted up 1.10x in OD -- evaporation
    or settling debris -- which an OD-only test reads as growth. Their fluorescence
    was flat (433 -> 426) while a real culture went 639 -> 1705. No cells means
    neither channel grows, so requiring both to be flat is the stronger test.
    """

    def _plate(self):
        t = np.linspace(0, 4.14, 25)
        od = pd.DataFrame({
            "H1": np.full(t.size, 0.094),
            "H2": np.linspace(0.090, 0.099, t.size),   # medium, OD drifts up
            "A1": np.linspace(0.168, 0.500, t.size),   # culture
        }, index=pd.Index(t, name="time_h"))
        rfu = pd.DataFrame({
            "H1": np.linspace(360, 357, t.size),
            "H2": np.linspace(433, 426, t.size),       # flat: no cells
            "A1": np.linspace(639, 1705, t.size),      # rises: cells
        }, index=od.index)
        return od, rfu

    def test_a_drifting_medium_well_is_rejected_on_od_alone(self):
        od, _ = self._plate()

        assert "H2" not in detect_blank_wells(RawOD(od))

    def test_adding_the_reporter_channel_recovers_it(self):
        od, rfu = self._plate()

        assert set(detect_blank_wells(RawOD(od), rfu=RawRFU(rfu))) == {"H1", "H2"}

    def test_a_real_culture_is_still_rejected_with_both_channels(self):
        od, rfu = self._plate()

        assert "A1" not in detect_blank_wells(RawOD(od), rfu=RawRFU(rfu))

    def test_a_well_growing_in_both_channels_is_not_a_blank(self):
        od, rfu = self._plate()
        od["H3"] = np.linspace(0.092, 0.101, len(od))
        rfu["H3"] = np.linspace(430, 900, len(od))     # fluorescence climbs: cells

        assert "H3" not in detect_blank_wells(RawOD(od), rfu=RawRFU(rfu))


def test_recorded_well_roles_outrank_low_od_blank_detection():
    t = np.linspace(0.0, 4.0, 25)
    od = pd.DataFrame({"G1": np.full_like(t, 0.08), "H1": np.full_like(t, 0.09)}, index=t)
    od.attrs["well_roles"] = {"G1": "culture", "H1": "blank"}

    assert detect_blank_wells(RawOD(od)) == ["H1"]


def test_recorded_blanks_are_not_replaced_by_a_culture_when_absent_from_the_channel():
    t = np.linspace(0.0, 4.0, 25)
    od = pd.DataFrame({"G1": np.full_like(t, 0.08)}, index=t)
    od.attrs["well_roles"] = {"G1": "culture", "H1": "blank"}

    assert detect_blank_wells(RawOD(od)) == []


def test_recorded_identity_in_the_reporter_channel_is_not_ignored():
    od = pd.DataFrame({"G1": [0.08, 0.08], "H1": [0.3, 0.5]})
    rfu = pd.DataFrame({"G1": [400.0, 400.0], "H1": [300.0, 300.0]})
    rfu.attrs["well_roles"] = {"G1": "culture", "H1": "blank"}

    assert detect_blank_wells(RawOD(od), rfu=RawRFU(rfu)) == ["H1"]


def test_conflicting_channel_role_metadata_is_refused():
    od = pd.DataFrame({"G1": [0.08, 0.08], "H1": [0.09, 0.09]})
    rfu = od * 1000.0
    od.attrs["well_roles"] = {"G1": "culture", "H1": "blank"}
    rfu.attrs["well_roles"] = {"G1": "blank", "H1": "culture"}

    with pytest.raises(ValueError, match="conflicting.*G1"):
        detect_blank_wells(RawOD(od), rfu=RawRFU(rfu))


def test_unknown_recorded_roles_are_not_reinterpreted_as_missing_metadata():
    od = pd.DataFrame({"G1": [0.08, 0.08]})
    od.attrs["well_roles"] = {"G1": "cultrue"}

    with pytest.raises(ValueError, match="role"):
        detect_blank_wells(RawOD(od))


def test_public_plate_dilution_report_uses_recorded_blank_identity():
    from ystwin.diagnostics.dilution_confound import plate_dilution_report

    t = np.linspace(0.0, 4.0, 25)
    od = pd.DataFrame({"G1": np.full_like(t, 0.08), "H1": np.full_like(t, 0.09),
                       "A1": 0.09 + 0.2 * np.exp(0.2 * t)}, index=t)
    rfu = pd.DataFrame({"G1": np.full_like(t, 500.0), "H1": np.full_like(t, 300.0),
                        "A1": 300.0 + 1000.0 * np.exp(0.2 * t)}, index=t)
    aligned = pd.concat({"OD600": od, "mCitrine": rfu}, axis=1)
    aligned.columns.names = ["channel", "well"]
    aligned.attrs["well_roles"] = {"G1": "culture", "H1": "blank", "A1": "culture"}

    result = plate_dilution_report(aligned, "OD600", "mCitrine")

    assert result.blank_wells == ["H1"]
    assert result.od_blank == pytest.approx(0.09)
    assert result.rfu_background == pytest.approx(300.0)
    with pytest.raises(ValueError, match="recorded.*culture"):
        plate_dilution_report(aligned, "OD600", "mCitrine", blank_wells=["G1"])
