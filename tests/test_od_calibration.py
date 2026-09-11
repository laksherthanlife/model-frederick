"""Move 3: turn a dilution series into the OD linear range, instead of a placeholder.

In a 96-well plate the optical path is short and scattering is multiple, so measured
absorbance saturates against true cell concentration. A one-parameter hyperbolic
saturation captures it and inverts in closed form:

    A_meas = A_true / (1 + k * A_true)      A_true = A_meas / (1 - k * A_meas)

The number the project actually needs from this is the measured OD at which the
deviation from linearity first exceeds a stated tolerance. That value configures
the G1 gate, which currently carries a guess.
"""

import numpy as np
import pytest

from ystwin.calib.od import fit_od_calibration
from ystwin.readings import CorrectedOD, RawOD


def _series(k=0.35, blank=0.09, top_od=2.4, n=12, noise=0.0, seed=0):
    """Two-fold dilution series of a dense culture, read through a saturating detector."""
    rng = np.random.default_rng(seed)
    dilution = 2.0 ** -np.arange(n)
    true_od = top_od * dilution
    measured = true_od / (1 + k * true_od) + blank
    return dilution, measured + rng.normal(0, noise, measured.shape), true_od


def test_recovers_a_known_saturation_constant_from_a_dilution_series():
    dilution, measured, _ = _series(k=0.35)

    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.saturation_k == pytest.approx(0.35, rel=0.02)


def test_linearisation_inverts_the_saturation_back_to_true_density():
    dilution, measured, true_od = _series(k=0.35)

    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.linearise(CorrectedOD(measured - 0.09)) == pytest.approx(true_od, rel=0.02)


def test_it_reports_the_measured_od_where_deviation_first_exceeds_tolerance():
    """The number that configures G1. It reduces exactly to tolerance / k."""
    dilution, measured, _ = _series(k=0.35)

    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.linear_range_max(tolerance=0.05) == pytest.approx(0.05 / 0.35, rel=0.02)


@pytest.mark.parametrize("k,tolerance", [(0.1, 0.05), (0.35, 0.05), (0.35, 0.02), (0.8, 0.10)])
def test_the_usable_ceiling_is_the_tolerance_divided_by_the_saturation_constant(k, tolerance):
    """A reader that saturates twice as hard halves the usable range, exactly."""
    dilution, measured, _ = _series(k=k)
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.linear_range_max(tolerance) == pytest.approx(tolerance / k, rel=0.03)


def test_a_reading_past_the_detector_asymptote_is_refused_not_extrapolated():
    dilution, measured, _ = _series(k=0.35)
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    with pytest.raises(ValueError, match="asymptote"):
        calib.linearise(CorrectedOD(np.array([1.0 / 0.35])))


def test_a_tighter_tolerance_gives_a_lower_usable_ceiling():
    dilution, measured, _ = _series(k=0.35)
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.linear_range_max(tolerance=0.02) < calib.linear_range_max(tolerance=0.10)


def test_a_detector_with_no_saturation_places_no_ceiling_on_the_range():
    dilution, measured, _ = _series(k=0.0)

    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.saturation_k == pytest.approx(0.0, abs=1e-3)
    assert calib.linear_range_max(tolerance=0.05) > 10.0


def test_the_fit_survives_realistic_reader_noise():
    dilution, measured, _ = _series(k=0.35, noise=0.004, seed=3)

    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    assert calib.saturation_k == pytest.approx(0.35, rel=0.20)


def test_a_series_with_no_dilution_variation_is_refused():
    with pytest.raises(ValueError, match="dilution"):
        fit_od_calibration(np.ones(6), RawOD(np.full(6, 0.5)), blank=0.09)


def test_dry_weight_conversion_uses_the_linearised_density():
    dilution, measured, _ = _series(k=0.35)
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09, gdcw_per_od=0.42)

    # 2.4 true OD units of biomass, read as a saturated value
    assert calib.to_dcw(CorrectedOD(measured[0] - 0.09)) == pytest.approx(2.4 * 0.42, rel=0.02)


def test_without_a_dry_weight_factor_the_conversion_is_refused_not_guessed():
    dilution, measured, _ = _series()
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    with pytest.raises(ValueError, match="gdcw_per_od"):
        calib.to_dcw(CorrectedOD(np.array([0.5])))


def test_the_calibration_configures_the_g1_gate_directly():
    """Closes the loop: the measured limit replaces the placeholder."""
    dilution, measured, _ = _series(k=0.35)
    calib = fit_od_calibration(dilution, RawOD(measured), blank=0.09)

    gate = calib.optical_gate(tolerance=0.05)

    assert gate.od_linear_max == pytest.approx(calib.linear_range_max(0.05))
    assert gate.od_linear_max < 1.0, "a real 96-well limit sits well below the cuvette figure"
