"""Fitting the panel's dose parameters to real plates, and saying what survives contact.

The panel's EC50s are literature estimates. The uploaded plates carry dose ladders for two
of them, so those two are testable: fit the induction EC50 and the lethal dose from the
measured activity and see whether the encoded values sit inside the fit.

What is fitted has to match what is identifiable. One reporter under one agent gives the
shape of the curve, so the two doses are recoverable; the amplitude is not, because
promoter strength, reporter gain and the target weight enter as one product and no single
ladder separates them. Reporting an amplitude as if it were a target weight would be
inventing a number, so the result says outright that it is confounded.

Wells that stopped growing are excluded rather than fitted. At 2 mM H2O2 blank-corrected
fluorescence is negative and growth is near zero -- that is an absent measurement, and
treating it as low activity would pull the lethal dose toward whatever the noise did.
"""

import numpy as np
import pytest

from ystwin.generator.panel_calibration import (
    calibrate_stressor,
    fit_dose_response,
    usable_points,
)


def curve(doses, ec50, lethal, amplitude=1.0, basal=0.0, hill_lethal=2.5):
    """Basal plus induction, the whole scaled by viability -- a cell that has stopped
    transcribing has stopped transcribing its basal too."""
    d = np.asarray(doses, dtype=float)
    return (basal + amplitude * (d / (ec50 + d))) / (1.0 + (d / lethal) ** hill_lethal)


class TestUsablePoints:
    def test_it_keeps_growing_wells(self):
        keep = usable_points(doses=[0, 0.5, 1.0], activity=[1.0, 1.2, 1.4],
                             growth=[0.35, 0.33, 0.22])

        assert list(keep) == [True, True, True]

    def test_it_drops_a_well_that_stopped_growing(self):
        keep = usable_points(doses=[0, 1.0, 2.0], activity=[1.0, 1.4, -0.02],
                             growth=[0.35, 0.22, 0.019])

        assert list(keep) == [True, True, False]

    def test_it_drops_negative_activity_whatever_the_growth(self):
        keep = usable_points(doses=[0, 1.0], activity=[1.0, -5.0], growth=[0.35, 0.30])

        assert list(keep) == [True, False]

    def test_it_refuses_mismatched_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            usable_points(doses=[0, 1], activity=[1.0], growth=[0.3, 0.3])


class TestFittingACurve:
    def test_it_recovers_the_doses_it_was_generated_from(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0])
        fit = fit_dose_response(doses, curve(doses, ec50=0.4, lethal=2.5))

        assert fit.ec50 == pytest.approx(0.4, rel=0.25)
        assert fit.lethal_dose == pytest.approx(2.5, rel=0.35)

    def test_it_reports_the_amplitude_as_confounded(self):
        doses = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 5.0])
        fit = fit_dose_response(doses, curve(doses, 0.4, 2.5, amplitude=7.0))

        assert fit.amplitude_is_confounded

    def test_it_survives_noise(self):
        rng = np.random.default_rng(0)
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0])
        clean = curve(doses, 0.4, 2.5)
        fit = fit_dose_response(doses, clean * rng.normal(1.0, 0.05, size=clean.shape))

        assert fit.ec50 == pytest.approx(0.4, rel=0.5)

    def test_it_refuses_a_ladder_too_short_to_fit(self):
        with pytest.raises(ValueError, match="at least"):
            fit_dose_response([0.0, 1.0], [0.0, 1.0])

    def test_it_reports_how_well_the_curve_matched(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0])
        fit = fit_dose_response(doses, curve(doses, 0.4, 2.5))

        assert fit.r_squared > 0.99


class TestVerdictAgainstTheEncodedValue:
    def test_a_matching_ladder_confirms_the_encoded_ec50(self):
        """The ladder has to match whatever is encoded, and that value moved.

        H2O2's EC50 was 0.5 mM by assertion and is now 0.15 mM from Goulev 2017. A
        ladder centred on 0.5 therefore refutes the encoded value rather than confirming
        it -- which is the gate working, not failing.
        """
        from ystwin.generator.stress_panel import STRESSORS

        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0])
        encoded = STRESSORS["H2O2"].ec50
        result = calibrate_stressor("H2O2", doses, curve(doses, encoded, 1.6),
                                    growth=np.full(len(doses), 0.3))

        assert result.verdict == "CONFIRMED"

    def test_a_ladder_that_disagrees_refutes_it(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0])
        result = calibrate_stressor("H2O2", doses, curve(doses, 0.02, 1.6),
                                    growth=np.full(len(doses), 0.3))

        assert result.verdict == "REFUTED"

    def test_too_few_usable_points_is_inconclusive_not_a_verdict(self):
        doses = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
        result = calibrate_stressor("H2O2", doses, [1.0, 1.2, -0.1, -0.2, -0.3],
                                    growth=[0.3, 0.3, 0.01, 0.01, 0.01])

        assert result.verdict == "INCONCLUSIVE"

    def test_it_carries_the_encoded_value_it_was_tested_against(self):
        from ystwin.generator.stress_panel import STRESSORS

        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0])
        result = calibrate_stressor("H2O2", doses, curve(doses, 0.5, 1.6),
                                    growth=np.full(len(doses), 0.3))

        assert result.encoded_ec50 == STRESSORS["H2O2"].ec50

    def test_an_unknown_stressor_is_refused(self):
        with pytest.raises(KeyError, match="unobtainium"):
            calibrate_stressor("unobtainium", [0, 1, 2, 3, 4], [0, 1, 2, 1, 0],
                               growth=[0.3] * 5)


class TestTheCurveNeedsABasalTerm:
    """An uninduced promoter is not silent, and a curve through the origin cannot say so.

    On the uploaded plates activity at zero dose is a large fraction of the peak, so a model
    forced through zero fits worse than a horizontal line -- R2 came back at -5. A negative
    R2 is the check that catches a wrong model form, and any verdict issued on top of one is
    meaningless, so the fit reports it and the verdict refuses to lean on it.
    """

    def test_it_fits_a_ladder_with_a_large_basal(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
        observed = curve(doses, ec50=0.4, lethal=2.5, amplitude=600.0, basal=1000.0)
        fit = fit_dose_response(doses, observed)

        assert fit.r_squared > 0.95
        assert fit.ec50 == pytest.approx(0.4, rel=0.4)

    def test_it_recovers_the_basal_it_was_given(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
        fit = fit_dose_response(doses, curve(doses, 0.4, 2.5, amplitude=400.0, basal=800.0))

        assert fit.basal == pytest.approx(800.0, rel=0.2)

    def test_a_zero_basal_ladder_still_fits(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
        fit = fit_dose_response(doses, curve(doses, 0.4, 2.5))

        assert fit.r_squared > 0.99

    def test_a_verdict_refuses_to_rest_on_a_bad_fit(self):
        """Whatever the fitted EC50 says, a curve that fits worse than a flat line says
        nothing about the encoded one."""
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0])
        noise = np.array([5.0, 1.0, 9.0, 2.0, 8.0, 3.0, 7.0])
        result = calibrate_stressor("H2O2", doses, noise, growth=np.full(len(doses), 0.3))

        assert result.verdict == "INCONCLUSIVE"
        assert result.note


class TestADegenerateFitIsNotAMeasurement:
    """A fitter always returns numbers; only some of them mean anything.

    Where a ladder has no falling limb the lethal dose is unconstrained, and an optimiser
    walks both parameters to wherever the bounds stop it. That produced an EC50 of 50 and a
    lethal dose of 500 on the H2O2 plates -- both exactly at their limits -- and a separate
    fit where the lethal dose came out below the inducing one, which no cell can do.

    Neither is a refutation of anything. Both must come back inconclusive, or the panel
    would be corrected against numbers the data never contained.
    """

    def test_a_parameter_pinned_at_its_bound_is_reported_as_unidentified(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0])
        rising = curve(doses, ec50=0.4, lethal=1e6, amplitude=500.0, basal=900.0)
        fit = fit_dose_response(doses, rising)

        assert not fit.identifiable

    def test_a_ladder_with_both_limbs_is_identified(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        fit = fit_dose_response(doses, curve(doses, 0.4, 2.5, amplitude=500.0, basal=900.0))

        assert fit.identifiable

    def test_a_lethal_dose_below_the_inducing_one_is_not_identified(self):
        doses = np.array([0.0, 0.05, 0.1, 0.2, 0.4])
        fit = fit_dose_response(doses, curve(doses, ec50=2.0, lethal=0.3, amplitude=500.0))

        assert not fit.identifiable

    def test_an_unidentified_fit_makes_the_verdict_inconclusive(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0])
        rising = curve(doses, ec50=0.4, lethal=1e6, amplitude=500.0, basal=900.0)
        result = calibrate_stressor("H2O2", doses, rising, growth=np.full(len(doses), 0.3))

        assert result.verdict == "INCONCLUSIVE"

    def test_it_says_why_the_ladder_could_not_answer(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0])
        rising = curve(doses, ec50=0.4, lethal=1e6, amplitude=500.0, basal=900.0)
        result = calibrate_stressor("H2O2", doses, rising, growth=np.full(len(doses), 0.3))

        assert result.verdict == "INCONCLUSIVE"
        assert "climbing" in result.note or "falling" in result.note


class TestAModelFreeCrossCheck:
    """Where a four-parameter fit is weak, the raw ladder can still say where half-max is.

    On the DTT plates the biphasic fit reaches R2 of only 0.53 on seven points, which is
    too little to set a parameter from on its own. Interpolating the dose at which activity
    crosses halfway from basal to peak assumes nothing about the curve's shape, and both
    UPRE replicates land on the same value -- which is what makes the refutation of the
    encoded EC50 safe to act on.
    """

    def test_it_finds_the_crossing_on_a_clean_ladder(self):
        from ystwin.generator.panel_calibration import half_maximal_dose

        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0])
        activity = np.array([100.0, 110.0, 125.0, 150.0, 200.0, 120.0])

        assert half_maximal_dose(doses, activity) == pytest.approx(0.5, abs=0.3)

    def test_it_is_unaffected_by_the_shape_beyond_the_peak(self):
        from ystwin.generator.panel_calibration import half_maximal_dose

        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
        rising = [100.0, 110.0, 125.0, 150.0, 200.0]

        gentle = half_maximal_dose(doses, np.array(rising + [180.0, 160.0]))
        steep = half_maximal_dose(doses, np.array(rising + [110.0, 100.0]))

        assert gentle == pytest.approx(steep)

    def test_a_ladder_that_never_rises_has_no_crossing(self):
        from ystwin.generator.panel_calibration import half_maximal_dose

        assert half_maximal_dose([0.0, 1.0, 2.0], [100.0, 100.0, 100.0]) is None

    def test_a_peak_at_the_first_dose_has_no_crossing(self):
        from ystwin.generator.panel_calibration import half_maximal_dose

        assert half_maximal_dose([0.0, 1.0, 2.0], [200.0, 150.0, 100.0]) is None

    def test_the_verdict_carries_it_beside_the_fitted_value(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        result = calibrate_stressor(
            "DTT", doses, curve(doses, 0.4, 2.5, amplitude=500.0, basal=900.0),
            growth=np.full(len(doses), 0.3))

        assert result.half_maximal_dose == pytest.approx(0.4, abs=0.3)


class TestWhatTheUploadedPlatesSettled:
    """What the plates could measure, and what they could not -- corrected.

    The first pass reported DTT's EC50 as measured at 0.304 mM and refuted the literature
    1.0. That was wrong three times over. The shape-free crossing it rested on is biased
    low by up to 87%; the biphasic fit that disagreed with it was misspecified, leaving
    systematic residuals of +12% at the peak and -22% past it; and the peak anchoring both
    sat at a dose where growth was already down to 71% of control, with the peroxide peak
    at 23%. Restricted to doses where the culture is healthy, both ladders are still
    climbing, so neither EC50 was ever identifiable and both keep their literature values.

    What the plates do measure is lethality, because growth measures it directly instead of
    inferring it from a reporter through a model of transcription. That route gave 7.2 mM
    for an agent whose cultures are at 24% of control by 2 mM.
    """

    def test_the_dtt_ec50_was_not_changed_by_the_plates(self):
        from ystwin.generator.stress_panel import STRESSORS

        assert STRESSORS["DTT"].ec50 == pytest.approx(1.0)

    def test_its_source_says_the_ladder_could_not_test_it(self):
        """The DTT EC50 must be marked as not established by the data.

        Asserting a phrase broke when the source was rewritten to record that 1.0 mM is
        unsourced and that the cited paper contains no such value. Asserting the meaning
        survives rewording and still fails if the caveat is ever dropped.
        """
        from ystwin.generator.stress_panel import STRESSORS

        source = STRESSORS["DTT"].source.lower()
        assert any(mark in source for mark in ("not tested", "unsourced", "cannot")), (
            "the DTT EC50 is not established by this data and the source must say so"
        )

    def test_the_peroxide_ec50_now_follows_goulev_rather_than_assertion(self):
        """0.5 mM was asserted; 0.15 mM comes from a measured nuclear-entry curve.

        Goulev 2017 has Yap1 nuclear entry partial at 0.1 mM and saturating above
        0.2 mM in S288C. At the old 0.5 the encoded ladder put its top four rungs where
        the response is already flat.
        """
        from ystwin.generator.stress_panel import STRESSORS

        assert STRESSORS["H2O2"].ec50 == pytest.approx(0.15)
        assert STRESSORS["H2O2"].ec50 < 0.2, "the response saturates above 0.2 mM"

    def test_both_lethal_doses_come_from_growth(self):
        """Measured where growth halves, not fitted to a reporter."""
        from ystwin.generator.stress_panel import STRESSORS

        assert STRESSORS["DTT"].lethal_dose == pytest.approx(1.55, abs=0.15)
        assert STRESSORS["H2O2"].lethal_dose == pytest.approx(1.0, abs=0.3)

    def test_dtt_has_a_narrow_window_which_is_why_the_ladder_failed(self):
        """Lethal at 1.55 against an inducing 1.0 leaves almost no room to saturate in."""
        from ystwin.generator.stress_panel import STRESSORS

        spec = STRESSORS["DTT"]
        assert spec.lethal_dose / spec.ec50 < 2.0

    def test_neither_ladder_saturates_before_the_culture_is_compromised(self):
        """The finding, stated as the generator now carries it."""
        from ystwin.generator.stress_panel import STRESSORS, module_response, viability

        for name in ("DTT", "H2O2"):
            spec = STRESSORS[name]
            healthy_edge = spec.lethal_dose * 0.7
            assert module_response(name, healthy_edge)[spec.targets and list(spec.targets)[0]] > 0
            assert viability(name, healthy_edge) > 0.5
