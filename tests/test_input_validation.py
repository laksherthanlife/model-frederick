"""Malformed input must fail by name. An unverified error path is an unverified path.

The time-format branches matter most: a plate reader that exports elapsed time as a
``datetime.time`` rather than a string must not silently produce a wrong clock.
"""

import datetime as dt

import numpy as np
import pandas as pd
import pytest

from ystwin.growth import specific_growth_rate
from ystwin.plate.synergy import _parse_channel_header, _to_hours, read_synergy_kinetic
from ystwin.reporter import (
    ReporterKinetics,
    naive_specific_fluorescence,
    promoter_activity,
    simulate_reporter,
)
from ystwin.readings import CorrectedOD, CorrectedRFU, SpecificFluorescence


class TestElapsedTimeParsing:
    def test_clock_string_with_seconds(self):
        assert _to_hours("02:30:00") == pytest.approx(2.5)

    def test_clock_string_without_hours_is_read_as_minutes_and_seconds(self):
        assert _to_hours("30:00") == pytest.approx(0.5)

    def test_timedelta_cell(self):
        assert _to_hours(dt.timedelta(hours=1, minutes=30)) == pytest.approx(1.5)

    def test_time_cell(self):
        assert _to_hours(dt.time(3, 15, 0)) == pytest.approx(3.25)

    def test_datetime_cell(self):
        assert _to_hours(dt.datetime(2026, 7, 1, 6, 30, 0)) == pytest.approx(6.5)

    def test_excel_serial_day_fraction(self):
        assert _to_hours(0.25) == pytest.approx(6.0)

    def test_a_boolean_is_not_a_timestamp(self):
        with pytest.raises(ValueError, match="unparseable"):
            _to_hours(True)

    def test_unparseable_text_is_rejected_by_name(self):
        with pytest.raises(ValueError, match="unparseable"):
            _to_hours("not a time")


class TestChannelHeaderParsing:
    def test_a_non_string_header_is_not_a_channel(self):
        assert _parse_channel_header(42) is None

    def test_a_header_with_no_name_is_not_a_channel(self):
        assert _parse_channel_header("T° :600") is None

    def test_the_repeat_read_suffix_is_stripped_from_the_optics(self):
        assert _parse_channel_header("T° mCitrine:480,530[2]") == ("mCitrine", "480,530")


class TestGrowthRateValidation:
    def test_mismatched_shapes_are_rejected(self):
        with pytest.raises(ValueError, match="same shape"):
            specific_growth_rate(np.linspace(0, 5, 10), CorrectedOD(np.ones(9)))

    def test_too_few_timepoints_are_rejected(self):
        with pytest.raises(ValueError, match="at least"):
            specific_growth_rate(np.linspace(0, 1, 4), CorrectedOD(np.ones(4)))

    def test_non_increasing_time_is_rejected(self):
        t = np.array([0.0, 1.0, 1.0, 2.0, 3.0, 4.0])
        with pytest.raises(ValueError, match="strictly increasing"):
            specific_growth_rate(t, CorrectedOD(np.ones_like(t)))

    def test_a_very_short_window_is_widened_rather_than_failing(self):
        t = np.linspace(0, 10, 12)
        mu = specific_growth_rate(t, CorrectedOD(0.05 * np.exp(0.3 * t)), window_h=0.01)

        assert np.all(np.isfinite(mu))


class TestReporterValidation:
    def test_negative_degradation_is_rejected(self):
        with pytest.raises(ValueError, match="k_deg"):
            ReporterKinetics(k_deg=-0.1)

    def test_nonpositive_optical_density_is_rejected_in_the_ratio(self):
        with pytest.raises(ValueError, match="positive"):
            naive_specific_fluorescence(CorrectedRFU(np.array([1.0, 2.0])), CorrectedOD(np.array([0.1, 0.0])))

    def test_too_few_timepoints_are_rejected(self):
        t = np.linspace(0, 1, 4)
        with pytest.raises(ValueError, match="at least 5"):
            promoter_activity(t, SpecificFluorescence(np.ones_like(t)), np.ones_like(t))

    def test_a_series_off_the_time_grid_is_rejected(self):
        t = np.linspace(0, 10, 20)
        with pytest.raises(ValueError, match="share the time grid"):
            promoter_activity(t, SpecificFluorescence(np.ones(19)), np.ones(20))

    def test_non_increasing_time_is_rejected(self):
        t = np.array([0.0, 1.0, 0.5, 2.0, 3.0, 4.0])
        with pytest.raises(ValueError, match="strictly increasing"):
            simulate_reporter(t, np.ones_like(t), np.ones_like(t))

    def test_supplying_explicit_initial_conditions_for_a_maturing_reporter(self):
        t = np.linspace(0, 10, 100)
        kin = ReporterKinetics(k_deg=0.1, k_mat=2.0)

        traced = simulate_reporter(t, np.ones_like(t), np.full_like(t, 0.2), kin, r0=5.0, i0=1.0)

        assert traced[0] == pytest.approx(5.0)


class TestWorkbookParsing:
    def test_a_workbook_with_no_kinetic_block_is_rejected_by_path(self, tmp_path):
        path = tmp_path / "plain.xlsx"
        pd.DataFrame({"a": [1, 2], "b": [3, 4]}).to_excel(path, index=False)

        with pytest.raises(ValueError, match="no kinetic channel blocks"):
            read_synergy_kinetic(path)

    def test_an_unknown_channel_lists_what_is_available(self, make_synergy_file):
        run = read_synergy_kinetic(make_synergy_file())

        with pytest.raises(KeyError, match="OD600"):
            run.channel("mScarlet")


@pytest.mark.parametrize("field", ["k_deg", "k_mat", "k_deg_immature"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_reporter_kinetics_require_finite_rates(field, value):
    with pytest.raises(ValueError, match=field):
        ReporterKinetics(**{"k_mat": 2.0, field: value})


def test_immature_degradation_requires_a_maturing_reporter():
    with pytest.raises(ValueError, match="k_deg_immature"):
        ReporterKinetics(k_deg_immature=0.2)


def test_an_unused_immature_initial_condition_is_rejected():
    t = np.linspace(0.0, 4.0, 25)
    with pytest.raises(ValueError, match="i0"):
        simulate_reporter(t, np.ones_like(t), np.full_like(t, 0.2), i0=1.0)


@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_reporter_derivatives_refuse_nonfinite_timestamps(value):
    t = np.linspace(0.0, 4.0, 25)
    t[12] = value
    with pytest.raises(ValueError, match="finite"):
        promoter_activity(t, SpecificFluorescence(np.ones_like(t)), np.full_like(t, 0.2))


@pytest.mark.parametrize("series", ["signal", "growth"])
def test_reporter_derivatives_refuse_nonfinite_series(series):
    t = np.linspace(0.0, 4.0, 25)
    signal, mu = np.ones_like(t), np.full_like(t, 0.2)
    (signal if series == "signal" else mu)[12] = np.nan
    with pytest.raises(ValueError, match="finite"):
        promoter_activity(t, SpecificFluorescence(signal), mu)


def test_reporter_time_grid_must_be_one_dimensional():
    t = np.linspace(0.0, 4.0, 25).reshape(5, 5)
    with pytest.raises(ValueError, match="one-dimensional"):
        promoter_activity(t, SpecificFluorescence(np.ones_like(t)), np.full_like(t, 0.2))


@pytest.mark.parametrize("window", [0.0, -1.0, np.nan, np.inf])
def test_reporter_smoothing_window_must_be_finite_and_positive(window):
    t = np.linspace(0.0, 4.0, 25)
    with pytest.raises(ValueError, match="window_h"):
        promoter_activity(t, SpecificFluorescence(np.ones_like(t)), np.full_like(t, 0.2),
                          window_h=window)


@pytest.mark.parametrize("field", ["entry_expression", "cassette", "dose", "setpoint"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_prediction_inputs_reject_nonfinite_values(field, value):
    from ystwin.predict import Environment, Genotype

    with pytest.raises(ValueError, match="finite"):
        if field == "entry_expression":
            Genotype(value)
        elif field == "cassette":
            Genotype(1.0, cassette={"crtYB": value})
        elif field == "dose":
            Environment(stressor="DTT", dose=value)
        else:
            Environment(growth_rate_setpoint_per_h=value)


def test_a_negative_dose_cannot_become_zero_stress():
    from ystwin.predict import Environment

    with pytest.raises(ValueError, match="dose"):
        Environment(stressor="DTT", dose=-1.0)
