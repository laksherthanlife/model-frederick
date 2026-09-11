"""S1: fit the generator to an uploaded plate.

The test that matters is the round trip -- generate from known parameters, calibrate
against the result, recover them. Anything the plate cannot identify must come back
flagged rather than fitted to a plausible-looking number.

The clearest example of that is the reporter gain. The observable is
``(RFU-bg)/(OD-blank) = (autofluorescence + gain*R) * gdcw_per_od``, in which ``gain``
and ``gdcw_per_od`` only ever appear multiplied. A plate alone cannot separate them,
and a fit that reported both would be inventing one.
"""

import numpy as np
import pytest

from ystwin.generator.calibrate import calibrate_from_plate
from ystwin.generator.plate import DEFAULT_PANEL, PlateConditions, generate_plate
from ystwin.readings import RawOD, RawRFU


@pytest.fixture(scope="module")
def calibrated():
    conditions = PlateConditions(seed=11)
    plate = generate_plate(DEFAULT_PANEL, conditions)
    result = calibrate_from_plate(
        RawOD(plate.od), RawRFU(plate.rfu), blank_wells=plate.blank_wells,
        layout=conditions.layout, doses=conditions.doses,
    )
    return plate, result


class TestPlateLevelTerms:
    def test_the_blank_is_recovered(self, calibrated):
        plate, result = calibrated

        assert result.conditions.od_blank == pytest.approx(plate.conditions.od_blank, rel=0.03)

    def test_the_reporter_background_is_recovered(self, calibrated):
        plate, result = calibrated

        assert result.conditions.optics.background == pytest.approx(
            plate.conditions.optics.background, rel=0.05
        )

    def test_the_reader_noise_is_estimated(self, calibrated):
        plate, result = calibrated

        assert result.conditions.reader_cv == pytest.approx(plate.conditions.reader_cv, abs=0.01)


class TestPerConstructKinetics:
    def test_every_construct_in_the_layout_gets_parameters(self, calibrated):
        _, result = calibrated

        assert set(result.panel) == {"UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"}

    def test_the_unstressed_growth_rate_is_recovered(self, calibrated):
        plate, result = calibrated

        for construct, truth in DEFAULT_PANEL.items():
            fitted = result.panel[construct].mu_max
            assert fitted == pytest.approx(truth.mu_max, rel=0.25), construct

    def test_the_growth_inhibition_ranks_the_constructs_correctly(self, calibrated):
        _, result = calibrated

        truth_order = sorted(DEFAULT_PANEL, key=lambda c: DEFAULT_PANEL[c].growth_ic50)
        fitted_order = sorted(result.panel, key=lambda c: result.panel[c].growth_ic50)
        assert fitted_order[0] == truth_order[0]

    def test_the_promoter_induction_is_recovered_in_relative_terms(self, calibrated):
        """Absolute activity is confounded with gain; the fold change is not."""
        _, result = calibrated

        for construct, truth in DEFAULT_PANEL.items():
            expected = truth.promoter_peak / truth.promoter_basal
            fitted = result.panel[construct].promoter_peak / result.panel[construct].promoter_basal
            assert fitted == pytest.approx(expected, rel=0.35), construct


class TestWhatThePlateCannotIdentify:
    def test_the_gain_and_dry_weight_factor_are_reported_only_as_their_product(self, calibrated):
        _, result = calibrated

        assert "gain" in result.unidentifiable
        assert result.gain_times_gdcw > 0

    def test_the_recovered_product_matches_the_one_that_generated_the_plate(self, calibrated):
        plate, result = calibrated
        truth = plate.conditions.optics.gain * plate.conditions.gdcw_per_od

        assert result.gain_times_gdcw == pytest.approx(truth, rel=0.25)

    def test_supplying_a_measured_dry_weight_factor_resolves_the_gain(self, calibrated):
        plate, _ = calibrated
        conditions = plate.conditions

        resolved = calibrate_from_plate(
            RawOD(plate.od), RawRFU(plate.rfu), blank_wells=plate.blank_wells,
            layout=conditions.layout, doses=conditions.doses,
            gdcw_per_od=conditions.gdcw_per_od,
        )

        assert "gain" not in resolved.unidentifiable
        assert resolved.conditions.optics.gain == pytest.approx(conditions.optics.gain, rel=0.25)


def test_the_calibration_reports_its_own_fit_quality(calibrated):
    _, result = calibrated

    assert set(result.fit_quality) == set(result.panel)
    assert all(0.0 <= v <= 1.0 for v in result.fit_quality.values())


def test_a_plate_with_no_blank_well_is_refused(calibrated):
    plate, _ = calibrated

    with pytest.raises(ValueError, match="blank"):
        calibrate_from_plate(
            RawOD(plate.od), RawRFU(plate.rfu), blank_wells=(),
            layout=plate.conditions.layout, doses=plate.conditions.doses,
        )


def test_calibrating_then_regenerating_reproduces_the_original_plate(calibrated):
    """The end of S1: the fitted generator stands in for the plate it was fitted to."""
    plate, result = calibrated

    regenerated = generate_plate(result.panel, result.conditions)

    original = plate.od.drop(columns=list(plate.blank_wells)).to_numpy()
    fresh = regenerated.od.drop(columns=list(regenerated.blank_wells)).to_numpy()
    assert np.median(np.abs(fresh - original) / original) < 0.15


class TestNoiseSourcesAreSeparated:
    """Replicate scatter mixes reader noise with well-to-well inoculum variation.

    They have different signatures: reader noise is independent per reading, while a
    different starting biomass is a fixed offset for the whole trace. Detrending each
    well separates them; using replicate scatter alone reports their sum as reader
    noise and overstates it several fold.
    """

    def _fit(self, reader_cv, well_cv):
        conditions = PlateConditions(seed=3, reader_cv=reader_cv, well_cv=well_cv)
        plate = generate_plate(DEFAULT_PANEL, conditions)
        return calibrate_from_plate(
            RawOD(plate.od), RawRFU(plate.rfu), blank_wells=plate.blank_wells,
            layout=conditions.layout, doses=conditions.doses,
        )

    def test_reader_noise_is_recovered_despite_large_well_variation(self):
        result = self._fit(reader_cv=0.01, well_cv=0.08)

        assert result.conditions.reader_cv == pytest.approx(0.01, abs=0.006)

    def test_well_variation_is_recovered_as_the_remainder(self):
        result = self._fit(reader_cv=0.01, well_cv=0.08)

        assert result.conditions.well_cv == pytest.approx(0.08, abs=0.04)

    def test_a_quiet_reader_is_not_credited_with_the_wells_spread(self):
        noisy_wells = self._fit(reader_cv=0.002, well_cv=0.10)

        assert noisy_wells.conditions.reader_cv < 0.02
