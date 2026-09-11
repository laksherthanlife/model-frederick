"""G1: is an optical channel quantitative for this well, before it informs biology?

Driven by what the real plates showed: cultures inoculated above the reader's
linear range, and OD traces that fall for a third of the run. A growth rate
estimated from a declining OD is not a growth rate, and every downstream
quantity -- dilution correction, promoter activity, latent state -- inherits it.
"""

import numpy as np

from ystwin.gates.g1_optical import OpticalQualityGate, assess_plate, assess_well
from ystwin.readings import RawOD, RawRFU


def _clean(hours=18.0, n=110, od0=0.08, mu=0.30, blank=0.09, cap=0.8):
    t = np.linspace(0, hours, n)
    od = od0 * np.exp(mu * t) / (1 + od0 * (np.exp(mu * t) - 1) / cap) + blank
    return t, od


def test_a_clean_exponential_culture_passes_every_check():
    t, od = _clean()

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert result.passed, result.failures


def test_a_culture_above_the_readers_linear_range_is_failed_by_name():
    # The real preliminary plate reached raw OD 1.9 in a 96-well format.
    t, od = _clean(od0=0.5, mu=0.25, cap=1.8)

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert not result.passed
    assert "linear_range" in result.failures


def test_a_persistently_declining_trace_is_failed_by_name():
    t, _ = _clean()
    od = 0.09 + 0.4 * np.exp(-0.05 * t)

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert "sustained_decline" in result.failures


def test_a_well_barely_above_the_blank_is_failed_by_name():
    t = np.linspace(0, 18, 110)
    od = 0.09 + 0.004 * np.exp(0.05 * t)

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert "above_blank" in result.failures


def test_a_well_with_no_dynamic_range_is_failed_by_name():
    t = np.linspace(0, 18, 110)
    od = np.full_like(t, 0.5)

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert "dynamic_range" in result.failures


def test_a_reporter_channel_at_background_is_failed_by_name():
    t, od = _clean()
    rfu = np.full_like(t, 260.0)

    result = assess_well(t, RawOD(od), od_blank=0.09, rfu=RawRFU(rfu), rfu_background=258.0)

    assert "reporter_above_background" in result.failures


def test_the_metrics_behind_each_verdict_are_returned_for_audit():
    t, od = _clean()

    result = assess_well(t, RawOD(od), od_blank=0.09)

    assert set(result.metrics) >= {"max_raw_od", "decline_fraction", "od_fold_change", "min_od_above_blank"}


def test_the_thresholds_used_are_recorded_with_the_result():
    t, od = _clean()
    gate = OpticalQualityGate(od_linear_max=0.75)

    result = assess_well(t, RawOD(od), od_blank=0.09, gate=gate)

    assert result.gate.od_linear_max == 0.75


def test_a_stricter_linear_limit_can_fail_a_well_that_otherwise_passes():
    t, od = _clean()

    lenient = assess_well(t, RawOD(od), od_blank=0.09, gate=OpticalQualityGate(od_linear_max=1.0))
    strict = assess_well(t, RawOD(od), od_blank=0.09, gate=OpticalQualityGate(od_linear_max=0.3))

    assert lenient.passed and not strict.passed


def test_plate_assessment_returns_one_row_per_well_with_its_failures():
    import pandas as pd

    t, good = _clean()
    bad = 0.09 + 0.4 * np.exp(-0.05 * t)
    frame = pd.DataFrame({"A1": good, "A2": bad}, index=pd.Index(t, name="time_h"))

    table = assess_plate(RawOD(frame), od_blank=0.09)

    assert list(table.well) == ["A1", "A2"]
    assert bool(table.set_index("well").loc["A1", "passed"])
    assert not bool(table.set_index("well").loc["A2", "passed"])
    assert "sustained_decline" in table.set_index("well").loc["A2", "failures"]
