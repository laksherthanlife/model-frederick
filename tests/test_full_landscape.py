"""The whole yeast stress landscape, and the two kinds of sensor that read it.

Seven regulons was never the landscape, only the part the first build reached. Yeast
answers osmotic, oxidative, ER, heat, cell-wall, DNA-damage, calcium, carbon, nitrogen,
hypoxia, metal, sulfur, xenobiotic, retrograde and pH stress through largely separate
transcription factors, and a latent state fitted to a fraction of them can only ever be
general over that fraction.

The second thing the landscape contains is a different sort of sensor. A promoter fusion
reports a regulon, accumulates over hours, matures slowly and is diluted by growth. A
ratiometric biochemical sensor -- roGFP2 for glutathione potential, HyPer for peroxide,
QUEEN for ATP, pHluorin for pH -- reports a metabolite pool directly, equilibrates in
seconds, is reversible, and its ratio cancels concentration, so growth does not touch it.
Treating those as one kind of measurement would put a dilution correction on a sensor that
has no dilution to correct, and would claim an hours-long integration for a reading that is
instantaneous.
"""

from ystwin.generator.stress_panel import (
    MODULES,
    REPORTERS,
    STRESSORS,
    Kind,
    module_response,
    ratiometric_reporters,
    reporter_loadings,
    transcriptional_reporters,
)


class TestTheLandscapeIsBroad:
    def test_it_covers_far_more_than_the_first_seven_regulons(self):
        assert len(MODULES) >= 15

    def test_it_names_the_stress_axes_a_yeast_actually_has(self):
        expected = {"ESR", "UPR", "oxidative", "heat", "osmotic", "proteasome", "iron",
                    "cell_wall", "dna_damage", "calcium", "carbon", "nitrogen",
                    "hypoxia", "copper", "zinc", "sulfur", "xenobiotic", "retrograde"}

        assert expected <= set(MODULES)

    def test_every_module_cites_a_source(self):
        assert all(m.source for m in MODULES.values())

    def test_every_module_names_its_factor_and_element(self):
        for name, module in MODULES.items():
            assert module.transcription_factor, name
            assert module.response_element, name

    def test_there_is_a_stressor_for_every_module(self):
        excited = {m for s in STRESSORS.values() for m in s.targets}
        driven = {d for m in MODULES.values() for d in m.driven_by}

        assert set(MODULES) <= excited | driven

    def test_every_stressor_cites_a_source_and_a_dose_scale(self):
        for name, stressor in STRESSORS.items():
            assert stressor.source, name
            assert stressor.ec50 > 0, name

    def test_no_stressor_targets_a_module_that_does_not_exist(self):
        for name, stressor in STRESSORS.items():
            for module in stressor.targets:
                assert module in MODULES, f"{name} targets unknown {module}"

    def test_the_cascade_still_has_no_cycle(self):
        resolved = set()
        for name in MODULES:
            for driver in MODULES[name].driven_by:
                assert driver != name
        for name, module in MODULES.items():
            if not module.driven_by:
                resolved.add(name)
        assert resolved


class TestTwoKindsOfSensor:
    def test_a_promoter_fusion_is_transcriptional(self):
        assert REPORTERS["STRE-general"].kind is Kind.TRANSCRIPTIONAL

    def test_a_glutathione_sensor_is_ratiometric(self):
        assert REPORTERS["roGFP2-Grx1"].kind is Kind.RATIOMETRIC

    def test_the_panel_offers_both(self):
        assert transcriptional_reporters()
        assert ratiometric_reporters()

    def test_the_two_sets_do_not_overlap(self):
        assert not set(transcriptional_reporters()) & set(ratiometric_reporters())

    def test_together_they_are_the_whole_panel(self):
        assert set(transcriptional_reporters()) | set(ratiometric_reporters()) == set(REPORTERS)

    def test_ratiometric_sensors_cover_the_pools_the_user_asked_for(self):
        read = {REPORTERS[r].module for r in ratiometric_reporters()}

        assert {"redox", "peroxide", "atp", "ph"} <= read

    def test_a_ratiometric_sensor_reads_a_pool_not_a_regulon(self):
        """Its module must be a metabolite state, which no promoter can report."""
        assert MODULES["redox"].response_element == "none"

    def test_every_reporter_declares_a_kind(self):
        assert all(r.kind in Kind for r in REPORTERS.values())


class TestRatiometricSensorsAreImmuneToGrowth:
    def test_a_ratiometric_sensor_declares_no_maturation_delay(self):
        """Its readout is an equilibrium ratio, so there is nothing to wait for."""
        assert REPORTERS["roGFP2-Grx1"].integrates_hours == 0.0

    def test_a_promoter_fusion_integrates_over_hours(self):
        assert REPORTERS["STRE-general"].integrates_hours > 0.0

    def test_the_two_kinds_are_distinguished_by_that_and_not_by_name(self):
        for name in ratiometric_reporters():
            assert REPORTERS[name].integrates_hours == 0.0
        for name in transcriptional_reporters():
            assert REPORTERS[name].integrates_hours > 0.0


class TestLoadingsAcrossTheWholePanel:
    def test_every_reporter_loads_on_its_own_module(self):
        names = list(REPORTERS)
        loadings = reporter_loadings(names)
        order = list(MODULES)

        for row, name in enumerate(names):
            assert loadings[row, order.index(REPORTERS[name].module)] > 0

    def test_the_matrix_is_reporters_by_modules(self):
        assert reporter_loadings(list(REPORTERS)).shape == (len(REPORTERS), len(MODULES))

    def test_a_stressor_moves_the_module_it_targets(self):
        assert module_response("NaCl", 0.5)["osmotic"] > 0

    def test_a_dna_damaging_agent_moves_the_damage_module(self):
        assert module_response("MMS", 0.02)["dna_damage"] > 0

    def test_a_cell_wall_agent_moves_the_wall_module(self):
        assert module_response("congo_red", 50.0)["cell_wall"] > 0
