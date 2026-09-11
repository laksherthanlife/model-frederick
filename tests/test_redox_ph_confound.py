"""Reconciling the two E_GSH literatures, and the confound that reconciliation exposes.

Two published ranges for the same compartment: -289 mV from rxYFP and -306 from
Grx1-roGFP2, both assuming pH 7.0, against -350 from roGFP2 read alongside a pHluorin that
measured 7.5. They are not in conflict. The glutathione couple is GSSG + 2H+ + 2e- -> 2GSH,
one proton per electron, so its potential moves 60 mV per pH unit, and half a unit accounts
for 30 of the 44 mV between Kojer and Ayer.

The consequence is larger than the reconciliation. Ayer measured roughly one pH unit of
acidification in the cytosol under H2O2 -- 60 mV of apparent shift -- against a genuine
E_GSH response to 1 mM H2O2 of 40 to 50 mV. **The artifact is bigger than the signal**, and
their own text says the pH change "contributed significantly to the observed changes in
redox state".

So a glutathione sensor read without a pH measurement cannot be interpreted under any
stressor that acidifies, which includes most of them. That is a property of the sensor the
panel has to carry, not a caveat for a discussion section.
"""

import pytest

from ystwin.generator.redox import (
    NERNST_MV_PER_PH,
    ROGFP2_MIDPOINT_MV,
    apparent_shift_from_ph,
    e_gsh_from_oxd,
)


class TestTheNernstConversion:
    def test_a_half_oxidised_probe_reads_its_own_midpoint(self):
        assert e_gsh_from_oxd(0.5, ph=7.0) == pytest.approx(ROGFP2_MIDPOINT_MV, abs=1.0)

    def test_a_more_oxidised_probe_reads_less_negative(self):
        assert e_gsh_from_oxd(0.6, ph=7.0) > e_gsh_from_oxd(0.2, ph=7.0)

    def test_a_more_alkaline_cytosol_reads_more_negative(self):
        assert e_gsh_from_oxd(0.12, ph=7.5) < e_gsh_from_oxd(0.12, ph=7.0)

    def test_the_ph_slope_is_sixty_millivolts_per_unit(self):
        one_unit = e_gsh_from_oxd(0.12, ph=7.0) - e_gsh_from_oxd(0.12, ph=8.0)

        assert one_unit == pytest.approx(NERNST_MV_PER_PH, rel=0.02)

    def test_a_fully_reduced_or_oxidised_probe_is_refused(self):
        for degree in (0.0, 1.0):
            with pytest.raises(ValueError, match="between"):
                e_gsh_from_oxd(degree, ph=7.0)


class TestTheTwoLiteraturesReconcile:
    def test_kojers_oxidation_degree_gives_kojers_number_at_the_ph_he_assumed(self):
        """Kojer 2012 reports OxD 0.12 and -306 mV, assuming pH 7.0."""
        assert e_gsh_from_oxd(0.12, ph=7.0) == pytest.approx(-306.0, abs=12.0)

    def test_the_same_measurement_at_the_measured_ph_moves_toward_ayer(self):
        assumed = e_gsh_from_oxd(0.12, ph=7.0)
        measured = e_gsh_from_oxd(0.12, ph=7.5)

        assert measured < assumed
        assert abs(measured - assumed) == pytest.approx(0.5 * NERNST_MV_PER_PH, rel=0.05)

    def test_ph_accounts_for_most_of_the_published_gap(self):
        """44 mV between Kojer and Ayer, of which half a pH unit is 30."""
        assert 0.5 * NERNST_MV_PER_PH > 0.6 * 44.0


class TestTheConfoundIsBiggerThanTheSignal:
    def test_one_ph_unit_of_acidification_looks_like_a_large_oxidation(self):
        assert apparent_shift_from_ph(from_ph=7.5, to_ph=6.5) == pytest.approx(
            NERNST_MV_PER_PH, rel=0.02)

    def test_it_exceeds_the_measured_response_to_peroxide(self):
        """Ayer 2013: 1 mM H2O2 shifts E_GSH by 40-50 mV and acidifies by about one unit."""
        artifact = abs(apparent_shift_from_ph(from_ph=7.5, to_ph=6.5))

        assert artifact > 50.0

    def test_the_sign_is_the_same_as_an_oxidation_which_is_why_it_deceives(self):
        """Acidifying makes the probe read less negative, exactly as oxidising does."""
        assert apparent_shift_from_ph(from_ph=7.5, to_ph=6.5) > 0

    def test_a_small_ph_drift_is_still_worth_a_reading(self):
        assert abs(apparent_shift_from_ph(from_ph=7.5, to_ph=7.3)) > 10.0


class TestThePanelCarriesThis:
    def test_the_glutathione_sensor_is_marked_ph_sensitive(self):
        from ystwin.generator.stress_panel import REPORTERS

        assert REPORTERS["roGFP2-Grx1"].ph_sensitive

    def test_the_peroxide_sensor_is_not(self):
        """HyPer7 was engineered pH-insensitive, which is much of the point of it."""
        from ystwin.generator.stress_panel import REPORTERS

        assert not REPORTERS["HyPer7"].ph_sensitive

    def test_a_promoter_fusion_is_not_marked_this_way(self):
        from ystwin.generator.stress_panel import Kind, REPORTERS

        for name, reporter in REPORTERS.items():
            if reporter.kind is Kind.TRANSCRIPTIONAL:
                assert not reporter.ph_sensitive, name

    def test_every_ph_sensitive_sensor_says_so_in_its_source(self):
        from ystwin.generator.stress_panel import REPORTERS

        for name, reporter in REPORTERS.items():
            if reporter.ph_sensitive:
                assert "pH" in reporter.source, name


class TestABuildMustBeAbleToInterpretItsOwnReadings:
    """A pH-sensitive sensor with no pH reading beside it produces uninterpretable data.

    The recommended build carried roGFP2 and no way to measure pH, and could not have added
    pHluorin if it wanted to -- both are green excitation-ratio sensors and share the slot.
    Under any acidifying stressor the redox channel would have returned a number nobody
    could attribute, and the model had no way to say so.

    Three ways out, and the selector should take one rather than leaving it to be noticed
    later: drop the sensor, swap it for a pH-insensitive one reading a related quantity, or
    measure pH in a parallel strain and accept that the strains must behave alike.
    """

    def test_the_flaw_is_detectable(self):
        from ystwin.analysis.sensor_selection import interpretable

        assert not interpretable(["STRE-general", "UPRE-ER", "roGFP2-Grx1"])

    def test_adding_a_ph_reading_fixes_it(self):
        from ystwin.analysis.sensor_selection import interpretable

        assert interpretable(["STRE-general", "roGFP2-Grx1", "pHluorin"])

    def test_but_that_pairing_cannot_be_built(self):
        """Which is the whole difficulty: the fix conflicts spectrally with the problem."""
        from ystwin.analysis.sensor_selection import spectral_conflict

        assert spectral_conflict(["roGFP2-Grx1", "pHluorin"])

    def test_a_ph_insensitive_sensor_needs_no_partner(self):
        from ystwin.analysis.sensor_selection import interpretable

        assert interpretable(["STRE-general", "UPRE-ER", "HyPer7"])

    def test_promoter_fusions_alone_are_always_interpretable(self):
        from ystwin.analysis.sensor_selection import interpretable
        from ystwin.generator.stress_panel import transcriptional_reporters

        assert interpretable(transcriptional_reporters()[:5])

    def test_the_recommended_build_is_now_one_that_can_be_read(self):
        from ystwin.analysis.sensor_selection import (
            RECOMMENDED_BUILD, interpretable, spectral_conflict)

        assert interpretable(RECOMMENDED_BUILD.stress_reporters)
        assert not spectral_conflict(RECOMMENDED_BUILD.stress_reporters)
