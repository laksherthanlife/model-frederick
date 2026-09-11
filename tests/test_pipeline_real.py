"""Real-plate numerical pipeline: committed measurements -> G1 -> D2.

These assertions encode what the actual wet-lab data showed. If a future change to
the estimators silently moves these numbers, that is the signal to re-derive the
conclusions rather than to relax the test. All measurements are public text; the
optional private-workbook parser coverage stays in test_synergy_parser_real.py.
"""

import numpy as np
import pytest

from ystwin.diagnostics.dilution_confound import detect_blank_wells, plate_dilution_report
from ystwin.gates.g1_optical import OpticalQualityGate, assess_plate
from ystwin.plate import replay
from ystwin.reporter import ReporterKinetics
from ystwin.readings import RawOD

pytestmark = pytest.mark.integration

GATE = OpticalQualityGate(od_linear_max=1.0, max_decline_fraction=0.15)


@pytest.fixture(scope="module")
def public_er_prelim():
    return replay.load_run("20260701_ER_preliminary (RAW).xlsx")


@pytest.fixture(scope="module")
def public_er_stress():
    return replay.load_run("20260709_ER_stress_1st (RAW).xlsx")


@pytest.fixture(scope="module")
def prepared(public_er_prelim):
    aligned = public_er_prelim.aligned("OD600")
    od, rfu = aligned["OD600"], aligned["mCitrine"]
    blanks = detect_blank_wells(RawOD(od))
    cultures = [w for w in od.columns if w not in blanks]
    return aligned, od, rfu, blanks, cultures


def test_media_blanks_are_recovered_from_the_real_plate(prepared):
    _, _, _, blanks, _ = prepared

    assert blanks == ["F11", "G11", "H11"]


def test_the_real_plate_reached_optical_densities_outside_the_linear_range(prepared):
    _, od, _, _, cultures = prepared

    assert od[cultures].to_numpy().max() > 1.5


def test_a_substantial_minority_of_real_wells_fail_optical_quality(prepared):
    _, od, rfu, blanks, cultures = prepared
    od_blank = float(od[blanks].to_numpy().mean())

    table = assess_plate(RawOD(od[cultures]), od_blank=od_blank, gate=GATE)

    assert 0.3 < table.passed.mean() < 0.8


def test_the_dominant_real_failure_modes_are_decline_and_linear_range(prepared):
    _, od, _, blanks, cultures = prepared
    od_blank = float(od[blanks].to_numpy().mean())

    table = assess_plate(RawOD(od[cultures]), od_blank=od_blank, gate=GATE)
    reasons = set(table[~table.passed].failures.str.split(",").explode())

    assert {"sustained_decline", "linear_range"} <= reasons


def test_every_surviving_real_culture_falsifies_the_stable_reporter_assumption(prepared):
    """The finding that blocks latent-state work on this data as it stands."""
    aligned, od, rfu, blanks, cultures = prepared
    od_blank = float(od[blanks].to_numpy().mean())
    survivors = assess_plate(RawOD(od[cultures]), od_blank=od_blank, gate=GATE).query("passed").well.tolist()

    report = plate_dilution_report(
        aligned.loc[:, (slice(None), survivors + blanks)],
        od_channel="OD600", reporter_channel="mCitrine",
        kinetics=ReporterKinetics(k_deg=0.0), blank_wells=blanks,
    )

    assert (report.per_well.verdict == "model-inadequate").all()
    assert report.per_well.implied_min_k_deg.median() > 0.01


def test_the_implied_reporter_half_life_is_hours_not_days(prepared):
    """A stable YFP should be diluted away, not degraded on this timescale."""
    aligned, od, _, blanks, cultures = prepared

    report = plate_dilution_report(
        aligned, od_channel="OD600", reporter_channel="mCitrine",
        kinetics=ReporterKinetics(k_deg=0.0), blank_wells=blanks,
    )
    half_life = np.log(2) / report.per_well.implied_min_k_deg.median()

    assert 2.0 < half_life < 48.0


def _analyse(run, gated):
    aligned = run.aligned("OD600")
    od = aligned["OD600"]
    blanks = detect_blank_wells(RawOD(od))
    cultures = [w for w in od.columns if w not in blanks]
    od_blank = float(od[blanks].to_numpy().mean())
    wells = cultures
    if gated:
        wells = assess_plate(RawOD(od[cultures]), od_blank=od_blank, gate=GATE).query("passed").well.tolist()
    return plate_dilution_report(
        aligned.loc[:, (slice(None), wells + blanks)],
        od_channel="OD600", reporter_channel="mCitrine",
        kinetics=ReporterKinetics(k_deg=0.0), blank_wells=blanks,
    )


def test_the_fitted_slope_has_the_sign_the_dilution_mechanism_predicts(public_er_stress):
    """log S = log k_synth - log(mu + k_deg) predicts slope +1 against -log(mu + k_deg).

    Every well where the statistic is estimable comes out positive, at roughly half
    to two thirds of the predicted magnitude -- consistent with dilution setting the
    gain, with promoter activity supplying the rest.
    """
    slopes = _analyse(public_er_stress, gated=True).per_well.dilution_slope.dropna()

    assert len(slopes) > 0
    assert (slopes > 0).all()
    assert 0.3 < slopes.median() < 1.2


def test_the_preliminary_plate_supports_no_dilution_conclusion_at_all(public_er_prelim):
    """Its statistic is estimable only in wells that fail optical QC.

    Quasi-steady state needs the reporter to relax within the run. On that plate it
    does so only where OD is high and falling -- exactly what G1 rejects. So the
    experiment carries no admissible evidence either way about the confound.
    """
    gated = _analyse(public_er_prelim, gated=True).per_well
    ungated = _analyse(public_er_prelim, gated=False).per_well

    assert gated.dilution_r2.notna().sum() == 0
    assert ungated.dilution_r2.notna().sum() > 0


def test_the_quasi_steady_state_limit_holds_for_only_part_of_each_real_trace(prepared):
    """Qualifies every R^2 reported from this data: it rests on a minority of points."""
    aligned, od, _, blanks, cultures = prepared
    od_blank = float(od[blanks].to_numpy().mean())
    survivors = assess_plate(RawOD(od[cultures]), od_blank=od_blank, gate=GATE).query("passed").well.tolist()

    report = plate_dilution_report(
        aligned.loc[:, (slice(None), survivors + blanks)],
        od_channel="OD600", reporter_channel="mCitrine",
        kinetics=ReporterKinetics(k_deg=0.0), blank_wells=blanks,
    )

    assert report.per_well.qss_fraction.median() < 0.5
