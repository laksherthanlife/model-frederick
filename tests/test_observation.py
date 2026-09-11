"""Move 1: the measurement layer of §8.1, with the reporter as a state.

Two things this makes explicit that the plan's `y_t ~ p(y | z_t, ...)` hides:

**The reporter is a state.** RFU reads the mature protein pool, which has its own
balance. The observation model maps state to signal; it does not map stress to signal.

**Beta-carotene attenuates the reporter.** The product absorbs hard at 450 nm, right
on mCitrine's excitation. So RFU falls as titre rises -- a coupling from the quantity
being predicted back into the channel predicting it. A shuffled-channel control
cannot catch this, because it is real signal. Only a non-producing isogenic strain,
or an explicit correction, separates them.
"""

import numpy as np
import pytest

from ystwin.observation import (
    ReporterOptics,
    correct_inner_filter,
    observe_od,
    observe_rfu,
)
from ystwin.readings import RawRFU

OPTICS = ReporterOptics(
    gain=2.0e5, background=250.0, autofluorescence=400.0, inner_filter_coeff=180.0
)
NO_PRODUCER = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                             inner_filter_coeff=0.0)


class TestOpticalDensity:
    def test_reads_total_biomass_through_the_calibrated_saturation(self):
        from ystwin.calib.od import ODCalibration

        calib = ODCalibration(saturation_k=0.35, top_true_od=2.4, blank=0.09)
        signal = observe_od(viable=0.8, dead=0.2, calibration=calib, gdcw_per_od=1.0)

        assert signal == pytest.approx(1.0 / (1 + 0.35 * 1.0) + 0.09, rel=1e-9)

    def test_dead_biomass_still_scatters_light(self):
        from ystwin.calib.od import ODCalibration

        calib = ODCalibration(saturation_k=0.0, top_true_od=2.0, blank=0.0)
        with_dead = observe_od(0.5, 0.5, calib, gdcw_per_od=1.0)
        viable_only = observe_od(0.5, 0.0, calib, gdcw_per_od=1.0)

        assert with_dead > viable_only


class TestReporterSignal:
    def test_scales_with_the_mature_pool_and_with_biomass(self):
        single = observe_rfu(reporter_mature=1e-4, biomass=0.5, carotenoid=0.0, optics=NO_PRODUCER)
        double = observe_rfu(reporter_mature=2e-4, biomass=0.5, carotenoid=0.0, optics=NO_PRODUCER)

        assert (double - NO_PRODUCER.background) == pytest.approx(
            2 * (single - NO_PRODUCER.background) - NO_PRODUCER.autofluorescence * 0.5, rel=1e-6
        )

    def test_autofluorescence_tracks_biomass_and_not_the_reporter(self):
        dark = observe_rfu(reporter_mature=0.0, biomass=1.0, carotenoid=0.0, optics=NO_PRODUCER)

        assert dark == pytest.approx(250.0 + 400.0, rel=1e-9)

    def test_an_empty_well_reads_only_the_media_background(self):
        blank = observe_rfu(reporter_mature=0.0, biomass=0.0, carotenoid=0.0, optics=NO_PRODUCER)

        assert blank == pytest.approx(250.0)

    def test_the_detector_saturates_when_a_ceiling_is_declared(self):
        capped = ReporterOptics(gain=2e5, background=250.0, autofluorescence=400.0,
                                inner_filter_coeff=0.0, detector_max=5000.0)

        assert observe_rfu(1e-1, 1.0, 0.0, capped) == pytest.approx(5000.0)

    def test_a_reading_below_the_ceiling_is_left_alone(self):
        capped = ReporterOptics(gain=2e5, background=250.0, autofluorescence=400.0,
                                inner_filter_coeff=0.0, detector_max=5000.0)
        uncapped = ReporterOptics(gain=2e5, background=250.0, autofluorescence=400.0,
                                  inner_filter_coeff=0.0)

        assert observe_rfu(1e-3, 1.0, 0.0, capped) == observe_rfu(1e-3, 1.0, 0.0, uncapped)


class TestInnerFilter:
    def test_beta_carotene_attenuates_the_reporter_signal(self):
        clean = observe_rfu(1e-4, 0.5, carotenoid=0.0, optics=OPTICS)
        pigmented = observe_rfu(1e-4, 0.5, carotenoid=4e-3, optics=OPTICS)

        assert pigmented < clean

    def test_a_non_producing_strain_shows_no_attenuation(self):
        """Which is why the isogenic non-producer is the control that detects it."""
        a = observe_rfu(1e-4, 0.5, carotenoid=0.0, optics=OPTICS)
        b = observe_rfu(1e-4, 0.5, carotenoid=4e-3, optics=NO_PRODUCER)

        assert a == pytest.approx(b)

    def test_it_manufactures_a_correlation_between_reporter_and_product(self):
        """The structural confound: RFU tracks titre with promoter activity fixed."""
        titre = np.linspace(0.0, 8e-3, 40)
        rfu = np.array([observe_rfu(1e-4, 0.5, c, OPTICS) for c in titre])

        assert np.corrcoef(titre, rfu)[0, 1] < -0.95

    def test_correcting_for_it_removes_the_manufactured_correlation(self):
        titre = np.linspace(0.0, 8e-3, 40)
        rfu = np.array([observe_rfu(1e-4, 0.5, c, OPTICS) for c in titre])

        corrected = correct_inner_filter(RawRFU(rfu), titre, OPTICS)

        assert np.std(corrected) / np.mean(corrected) < 1e-9

    def test_the_correction_is_the_exact_inverse_of_the_attenuation(self):
        raw = observe_rfu(1e-4, 0.5, carotenoid=5e-3, optics=OPTICS)
        unpigmented = observe_rfu(1e-4, 0.5, carotenoid=0.0, optics=OPTICS)

        corrected = correct_inner_filter(RawRFU(np.array([raw])), np.array([5e-3]), OPTICS)

        assert corrected.values[0] == pytest.approx(unpigmented, rel=1e-9)

    def test_correcting_with_an_uncalibrated_coefficient_is_refused(self):
        unknown = ReporterOptics(gain=2e5, background=250.0, autofluorescence=400.0,
                                 inner_filter_coeff=None)

        with pytest.raises(ValueError, match="inner_filter_coeff"):
            correct_inner_filter(RawRFU(np.array([1000.0])), np.array([1e-3]), unknown)


def test_negative_optical_parameters_are_refused():
    with pytest.raises(ValueError, match="gain"):
        ReporterOptics(gain=-1.0, background=250.0, autofluorescence=400.0, inner_filter_coeff=0.0)
