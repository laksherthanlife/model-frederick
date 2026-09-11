"""The regulatory wiring, corrected against the literature rather than assumed.

Two different mechanisms were being carried by one parameter. A module-to-module cascade
is one transcription factor driving the GENE for another, so the downstream regulon really
does rise; promoter crosstalk is a reporter's own promoter carrying a second response
element, so it reads a module its regulon has nothing to do with. Using one weight for
both double-counted the first and misplaced the second.

Most of the cascade edges assumed here do not exist. UPR activation is Ire1-mediated HAC1
splicing and HAC1 has no STRE, so a general-stress arm into it is not a gene-expression
edge at all. Hsf1 and Msn2/4 run in parallel with compensatory, negatively signed
coupling. Hog1 acts on Msn2/4 activity, which points the opposite way to the arrow drawn.
What survives is oxidative and heat driving the proteasome regulon through RPN4, which is
documented by promoter mutagenesis and ChIP.
"""

from ystwin.generator.stress_panel import MODULES, REPORTERS, module_response


class TestCascadeEdgesMatchTheLiterature:
    def test_the_general_stress_module_does_not_drive_the_upr(self):
        """HAC1 is spliced, not transcribed by Msn2/4, and carries no STRE."""
        assert "ESR" not in MODULES["UPR"].driven_by

    def test_the_general_stress_module_does_not_drive_the_heat_module(self):
        """Hsf1 and Msn2/4 act in parallel, and removing one strengthens the other."""
        assert "ESR" not in MODULES["heat"].driven_by

    def test_the_general_stress_module_does_not_drive_the_oxidative_module(self):
        assert "ESR" not in MODULES["oxidative"].driven_by

    def test_the_osmotic_arrow_is_not_drawn_from_the_general_stress_module(self):
        """Hog1 acts on Msn2/4 activity, so any edge runs the other way and is signalling."""
        assert "ESR" not in MODULES["osmotic"].driven_by

    def test_the_oxidative_module_drives_the_proteasome(self):
        """Yap1 -> RPN4 through a YRE, confirmed by promoter mutagenesis."""
        assert "oxidative" in MODULES["proteasome"].driven_by

    def test_the_heat_module_drives_the_proteasome(self):
        """Hsf1 -> RPN4 through an HSE, confirmed by mutagenesis and ChIP."""
        assert "heat" in MODULES["proteasome"].driven_by

    def test_the_iron_regulon_is_independent_of_general_stress(self):
        assert "ESR" not in MODULES["iron"].driven_by

    def test_every_cascade_edge_names_a_real_module(self):
        for name, module in MODULES.items():
            for driver in module.driven_by:
                assert driver in MODULES, f"{name} is driven by unknown {driver}"


class TestTheCascadeActuallyPropagates:
    def test_an_oxidative_stressor_raises_the_proteasome_beyond_its_direct_target(self):
        """H2O2 lists no proteasome target of its own, so what appears comes via Yap1."""
        from ystwin.generator.stress_panel import STRESSORS

        response = module_response("H2O2", 1.0)

        assert "proteasome" not in STRESSORS["H2O2"].targets
        assert response["proteasome"] > 0.0

    def test_a_purely_general_stressor_leaves_the_upr_alone(self):
        response = module_response("glucose_starvation", 1.0)

        assert response["ESR"] > 0.0
        assert response["UPR"] == 0.0

    def test_it_leaves_the_iron_regulon_alone(self):
        assert module_response("glucose_starvation", 1.0)["iron"] == 0.0

    def test_an_unstressed_culture_has_no_activity_anywhere(self):
        assert all(v == 0.0 for v in module_response("DTT", 0.0).values())

    def test_the_cascade_cannot_run_away(self):
        response = module_response("H2O2", 1e6)

        assert all(v < 100.0 for v in response.values())


class TestPromoterCrosstalkIsSeparateFromTheCascade:
    def test_a_reporter_declares_the_elements_in_its_promoter(self):
        assert REPORTERS["HSE-heat"].also_reads

    def test_the_heat_reporter_also_carries_a_general_stress_element(self):
        """HSP104 is documented to carry both HSE and STRE."""
        assert "ESR" in REPORTERS["HSE-heat"].also_reads

    def test_a_single_element_reporter_reads_only_its_own_module(self):
        assert not REPORTERS["FeRE-iron"].also_reads

    def test_promoter_crosstalk_never_names_the_reporters_own_module(self):
        for name, reporter in REPORTERS.items():
            assert reporter.module not in reporter.also_reads, name

    def test_every_crosstalk_target_is_a_real_module(self):
        for name, reporter in REPORTERS.items():
            for module in reporter.also_reads:
                assert module in MODULES, f"{name} reads unknown {module}"
