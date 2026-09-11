"""The yeast stress landscape, not just ER and oxidative.

Modules and their crosstalk are transcribed from the literature. The structure that
matters is that they are *not* orthogonal: Msn2/4 sits under most of them, Rpn4 is
induced by Msn2/4, Yap1 and Hsf1 together, and Hog1 routes its transcriptional output
through Msn2. A panel that assumed orthogonal sensors would over-count the states it
can identify.
"""

import pytest

from ystwin.generator.stress_panel import (
    MODULES,
    REPORTERS,
    STRESSORS,
    integration_hours,
    module_ec50,
    module_response,
    reporter_loadings,
)


class TestTheModules:
    def test_the_major_yeast_stress_regulons_are_present(self):
        assert {"ESR", "UPR", "oxidative", "heat", "osmotic", "proteasome", "iron"} <= set(MODULES)

    def test_each_names_its_transcription_factor_and_element(self):
        for name, module in MODULES.items():
            assert module.transcription_factor, name
            assert module.response_element, name
            assert module.source, name

    def test_the_general_stress_module_is_msn2_msn4(self):
        assert "Msn2" in MODULES["ESR"].transcription_factor
        assert MODULES["ESR"].response_element == "STRE"


class TestCrosstalkIsEncoded:
    """Cascade edges are sparse on purpose; the wiring is checked in test_network_wiring."""

    def test_the_proteasome_module_answers_to_the_two_factors_shown_at_its_promoter(self):
        """Yap1 and Hsf1 sites in the RPN4 promoter, both confirmed by mutagenesis."""
        assert set(MODULES["proteasome"].driven_by) == {"oxidative", "heat"}

    def test_osmotic_signalling_is_not_encoded_as_a_transcriptional_edge(self):
        """Hog1 acts on Msn2/4 activity, which is signalling and points the other way."""
        assert MODULES["osmotic"].driven_by == {}

    def test_the_iron_regulon_is_largely_independent(self):
        assert "ESR" not in MODULES["iron"].driven_by

    def test_no_module_drives_itself(self):
        for name, module in MODULES.items():
            assert name not in module.driven_by, name


class TestStressors:
    def test_the_library_covers_more_than_two_axes(self):
        assert len(STRESSORS) >= 6

    def test_dtt_hits_the_er_module_hardest(self):
        response = module_response("DTT", dose=1.0)
        assert max(response, key=response.get) == "UPR"

    def test_the_peroxide_pool_leads_its_regulon_threefold(self):
        """A threefold lead, asserted twice before on two wrong grounds and now measured.

        It was 3.3-fold from a blanket multiplier on a 0.5 mM oxidative EC50. Correcting
        that EC50 to 0.15 mM left the pool arm pinned and the gap looked closed -- this
        test asserted the two were co-timed within 20%. Neither reading was a measurement
        of the pool. HyPer7 detects cytosolic H2O2 at ~20 uM in BY4742 (Kritsiligkou 2021,
        PMID 34118234) against 100 uM for the first oxidised Yap1 (Delaunay 2000, PMID
        11013218), which puts the lead back at threefold on evidence.
        """
        pool = module_ec50("H2O2", "peroxide")
        regulon = module_ec50("H2O2", "oxidative")

        assert regulon / pool == pytest.approx(3.0)

    def test_a_pool_leads_only_where_it_needs_no_enzyme(self):
        """The ordering is not a rule about pools, it is a rule about routes. Acid entering
        the cytosol leads; a thiol pool behind a peroxidase does not, and DTT's own
        glutathione arm asserts no separation because none is measured."""
        assert module_ec50("acetic_acid", "ph") < module_ec50("acetic_acid", "ESR")
        assert module_ec50("DTT", "redox") == pytest.approx(module_ec50("DTT", "UPR"))

    def test_every_stressor_also_moves_the_general_response(self):
        """Which is why a general-stress reporter alone cannot say what happened."""
        for name in STRESSORS:
            assert module_response(name, dose=1.0)["ESR"] > 0.0

    def test_zero_dose_moves_nothing(self):
        assert all(v == pytest.approx(0.0) for v in module_response("DTT", dose=0.0).values())

    def test_response_grows_along_the_rising_limb(self):
        """Only along the rising limb: 2 mM peroxide killed every well on the plates, and
        the generator now reproduces that rather than extrapolating a ceiling."""
        low = module_response("H2O2", dose=0.2)["oxidative"]
        high = module_response("H2O2", dose=0.5)["oxidative"]
        assert high > low

    def test_response_falls_again_past_the_lethal_dose(self):
        peak = module_response("H2O2", dose=0.5)["oxidative"]
        assert module_response("H2O2", dose=4.0)["oxidative"] < peak

    def test_an_unknown_stressor_is_refused(self):
        with pytest.raises(KeyError, match="unobtainium"):
            module_response("unobtainium", dose=1.0)


class TestReporterLoadings:
    def test_every_reporter_names_the_module_it_reads(self):
        for name, reporter in REPORTERS.items():
            assert reporter.module in MODULES, name

    def test_a_reporter_picks_up_the_extra_elements_in_its_promoter(self):
        """HSP104 carries an STRE beside its HSE, so it reads the general response too."""
        row = reporter_loadings(["HSE-heat"])[0]

        assert (row > 0).sum() > 1

    def test_a_single_element_promoter_reads_only_its_own_module(self):
        """PACE alone: everything upstream reaches RPT1 through Rpn4, not at this promoter."""
        row = reporter_loadings(["PACE-proteasome"])[0]

        assert (row > 0).sum() == 1

    def test_a_clean_reporter_loads_mostly_on_one_module(self):
        loadings = reporter_loadings(["FeRE-iron"])
        assert (loadings[0] > 0.01).sum() == 1

    def test_the_matrix_has_one_row_per_reporter_and_column_per_module(self):
        names = ["STRE-general", "UPRE-ER", "TRX2-oxidative"]
        assert reporter_loadings(names).shape == (3, len(MODULES))

    def test_an_unknown_reporter_is_refused(self):
        with pytest.raises(KeyError, match="mystery"):
            reporter_loadings(["mystery"])


class TestTheTwoTopologyCorrections:
    """Both were errors about *which module* an arm belongs to, not about its size.

    A magnitude error scales a column of the loading matrix. These changed its shape: one
    deleted an off-diagonal entry that was never off-diagonal, the other moved a reporter's
    dominant arm off its own module. Pinned because both survived prose review -- the
    sentences beside them were true, and named the wrong module.
    """

    def test_trx2_carries_no_crosstalk_because_its_second_arm_is_its_own_module(self):
        """Skn7 binds TRX2 directly (Morgan 1997) and is already inside `oxidative`.

        The residual induction in yap1 is therefore on-diagonal. Encoding it as a 0.30 term
        to the ESR attributed a Skn7 effect to Msn2/4 and gave `L` an off-diagonal the
        promoter does not have.
        """
        assert REPORTERS["TRX2-oxidative"].also_reads == {}
        assert (reporter_loadings(["TRX2-oxidative"])[0] > 0).sum() == 1

    def test_trx2_cites_skn7_for_the_yap1_independent_arm(self):
        source = REPORTERS["TRX2-oxidative"].source

        assert "Skn7" in source
        assert "9118942" in source

    def test_the_ena1_reporter_reads_calcium_harder_than_its_own_module(self):
        """Crz1 supplies about 60% of the early alkaline response, Rim101 and Snf1 the rest.

        At 0.30 the model had the 60% arm as the minor one. The weight is relative to a
        self-loading of 1.0, so 60:40 is 1.5.
        """
        weights = dict(zip(list(MODULES), reporter_loadings(["RIM101-alkaline"])[0]))

        assert weights["calcium"] > weights["alkaline_ph"]
        assert weights["calcium"] / weights["alkaline_ph"] == pytest.approx(1.5)

    def test_the_ena1_reporter_cites_the_paper_that_measured_the_split(self):
        assert "27362362" in REPORTERS["RIM101-alkaline"].source


class TestTheCorrectedDoseParameters:
    def test_the_peroxide_oxidative_ec50_sits_where_yap1_actually_moves(self):
        """Goulev 2017: nuclear entry partial at 0.1 mM, saturating above 0.2 mM."""
        assert 0.1 <= module_ec50("H2O2", "oxidative") <= 0.2

    def test_the_pool_arms_are_measured_not_pinned_to_the_old_derivation(self):
        """0.175 was defended as a value the correction "has no claim on". It was not a
        measurement: it was 0.5 x the old blanket pool factor, so it carried the very EC50
        being corrected. Both arms now come from papers that measured pools.
        """
        # Kritsiligkou 2021: HyPer7 detects ~20 uM, so the pool leads; 0.05 brackets
        # that floor and the 0.15 regulon EC50.
        assert module_ec50("H2O2", "peroxide") == pytest.approx(0.05)
        # Ayer 2013: half the roGFP2 probe oxidised at 2 mM -- a real half-activation.
        assert module_ec50("H2O2", "redox") == pytest.approx(2.0)
        assert module_ec50("H2O2", "peroxide") != pytest.approx(0.175)
        assert module_ec50("H2O2", "redox") != pytest.approx(0.175)

    def test_the_derived_esr_ec50_no_longer_exceeds_the_lethal_dose(self):
        """It did, at 1.25 mM against a measured 1.0 mM: the model claimed the general
        response half-activates at a dose that has already halved transcription."""
        assert module_ec50("H2O2", "ESR") < STRESSORS["H2O2"].lethal_dose

    def test_the_dtt_ec50_no_longer_claims_a_source_that_does_not_contain_it(self):
        """PMC7646510 is MacGilvray 2020, which doses 2.5 mM and reports no 1.0 mM value.

        The identifier stays in the string deliberately: a retracted attribution that is
        named and corrected is provenance, and one that is quietly deleted looks like it was
        never claimed.
        """
        source = STRESSORS["DTT"].source

        assert "contains no 1.0 mM value" in source
        assert "32597660" in source

    def test_the_dtt_ec50_records_the_disagreement_instead_of_resolving_it(self):
        """Three measurements differ by more than twofold and nothing here settles them.

        Moving the value to Pincus's 2.2 mM would place half-induction above this project's
        own measured lethal dose, so the value stays and the contradiction is written down
        beside it with the experiments that would decide it.
        """
        source = STRESSORS["DTT"].source

        assert STRESSORS["DTT"].ec50 == 1.0
        assert "20625545" in source
        assert "2.2" in source
        assert "UNSOURCED" in source


class TestTheIntegrationWindowIsDerived:
    def test_it_is_one_over_the_growth_rate_rather_than_a_literal(self):
        """A promoter fusion's time constant is 1/(mu + k_deg) and nothing else. Fixed at
        3.0 h it encoded mu = 0.33 invisibly; derived, it follows whatever mu is."""
        from ystwin.generator.stress_panel import _MU_MAX

        assert REPORTERS["STRE-general"].integrates_hours == pytest.approx(1.0 / _MU_MAX)

    def test_it_transfers_to_a_host_by_supplying_that_host_s_growth_rate(self):
        """E. coli at mu ~ 1/h integrates over a third of the window yeast does, which is a
        design change rather than a re-parameterisation."""
        assert integration_hours(1.0) == pytest.approx(1.0)
        assert integration_hours(0.30) > integration_hours(0.46)

    def test_a_culture_that_never_dilutes_is_refused_rather_than_returning_infinity(self):
        with pytest.raises(ValueError, match="never relaxes"):
            integration_hours(0.0)

    def test_a_ratiometric_sensor_still_integrates_over_nothing(self):
        """It reports an equilibrium between two forms of one molecule, so growth-driven
        dilution never enters and the derivation does not apply to it."""
        assert REPORTERS["roGFP2-Grx1"].integrates_hours == 0.0
