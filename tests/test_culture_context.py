"""The physiological baseline a stress sits on, which every simulated culture shared.

Training varied the stressor and nothing else, so every well was exponential-phase glucose
at 30 degrees. A stress state learned there describes one corner of the space a real plate
occupies, and two of the twenty-four modules are set by the baseline before any agent is
added: the carbon regulon is derepressed the moment glucose runs out, and the general
stress response is high in stationary phase whatever else is happening.

The consequence is not a small loss of accuracy. A model trained only on exponential
glucose has never seen the ESR high without a stressor, so it reads late-plate ESR as
stress -- which is exactly the reading a plate gives at the end of a run.

Context shifts the baseline, not the response. A stressor's effect is still its own; what
changes is what it is added to, and how fast the culture is growing while it is measured.
"""

import pytest

from ystwin.generator.context import (
    CultureContext,
    baseline_activity,
    context_growth_rate,
)
from ystwin.generator.stress_panel import MODULES


class TestTheContextItself:
    def test_the_default_is_exponential_glucose(self):
        context = CultureContext()

        assert context.carbon_source == "glucose"
        assert context.growth_phase == "exponential"

    def test_it_refuses_a_carbon_source_it_has_no_physiology_for(self):
        with pytest.raises(ValueError, match="carbon source"):
            CultureContext(carbon_source="unobtainium")

    def test_it_refuses_a_growth_phase_it_has_no_physiology_for(self):
        with pytest.raises(ValueError, match="growth phase"):
            CultureContext(growth_phase="interpretive")

    def test_oxygen_above_one_is_not_a_mole_fraction(self):
        with pytest.raises(ValueError, match="cannot exceed 1.0"):
            CultureContext(oxygen=2.0)

    def test_enrichment_above_air_is_refused_rather_than_answered(self):
        """The subtle one, and the reason it is a refusal and not a clamp.

        The aeration term is ``1 - respiring * max(0, 1 - O2/0.21)``, so it saturates at air:
        pure oxygen and air return the *same* growth rate. That number is not wrong --
        respiration is oxygen-saturated well below air -- but returning it silently asserts
        that enrichment changes nothing, when raised pO2 raises ROS and the oxidative burden
        the stress panel exists to carry. A sweep that dialled oxygen up would see a flat line
        and read it as evidence.
        """
        with pytest.raises(ValueError, match="enriched above air"):
            CultureContext(oxygen=1.0)

    def test_air_itself_is_still_allowed(self):
        """The boundary is inclusive, so the default context is not caught by its own guard."""
        assert CultureContext(oxygen=0.21).oxygen == 0.21
        assert CultureContext().oxygen == 0.21

    def test_it_carries_a_temperature_and_an_oxygen_level(self):
        context = CultureContext(temperature_c=37.0, oxygen=0.1)

        assert context.temperature_c == 37.0
        assert context.oxygen == 0.1


class TestGrowthRateFollowsContext:
    def test_glucose_grows_fastest(self):
        glucose = context_growth_rate(CultureContext(carbon_source="glucose"))
        galactose = context_growth_rate(CultureContext(carbon_source="galactose"))
        ethanol = context_growth_rate(CultureContext(carbon_source="ethanol"))

        assert glucose > galactose > ethanol

    def test_stationary_phase_barely_grows(self):
        rate = context_growth_rate(CultureContext(growth_phase="stationary"))

        assert rate < 0.05

    def test_low_oxygen_slows_a_respiring_culture_most(self):
        """Ethanol has to be respired; glucose can be fermented."""
        on_ethanol = (context_growth_rate(CultureContext(carbon_source="ethanol", oxygen=0.05))
                      / context_growth_rate(CultureContext(carbon_source="ethanol")))
        on_glucose = (context_growth_rate(CultureContext(carbon_source="glucose", oxygen=0.05))
                      / context_growth_rate(CultureContext(carbon_source="glucose")))

        assert on_ethanol < on_glucose

    def test_growth_is_never_negative(self):
        harsh = CultureContext(carbon_source="ethanol", growth_phase="stationary", oxygen=0.0)

        assert context_growth_rate(harsh) >= 0.0


class TestBaselineActivityFollowsContext:
    def test_an_exponential_glucose_culture_has_no_baseline_stress(self):
        baseline = baseline_activity(CultureContext())

        assert all(abs(v) < 1e-9 for v in baseline.values())

    def test_stationary_phase_raises_the_general_stress_response(self):
        assert baseline_activity(CultureContext(growth_phase="stationary"))["ESR"] > 0.3

    def test_losing_glucose_derepresses_the_carbon_regulon(self):
        assert baseline_activity(CultureContext(carbon_source="ethanol"))["carbon"] > 0.3

    def test_low_oxygen_raises_the_hypoxic_regulon(self):
        assert baseline_activity(CultureContext(oxygen=0.02))["hypoxia"] > 0.3

    def test_a_respiratory_carbon_source_raises_the_retrograde_response(self):
        assert baseline_activity(CultureContext(carbon_source="ethanol"))["retrograde"] > 0.0

    def test_a_warm_culture_carries_some_heat_response(self):
        assert baseline_activity(CultureContext(temperature_c=37.0))["heat"] > 0.0

    def test_every_named_module_is_real(self):
        for context in (CultureContext(), CultureContext(growth_phase="stationary"),
                        CultureContext(carbon_source="ethanol", oxygen=0.02)):
            assert set(baseline_activity(context)) <= set(MODULES)

    def test_a_baseline_never_exceeds_a_full_response(self):
        for phase in ("exponential", "diauxic", "stationary"):
            for carbon in ("glucose", "galactose", "ethanol"):
                context = CultureContext(carbon_source=carbon, growth_phase=phase, oxygen=0.0)
                assert all(0.0 <= v <= 1.0 for v in baseline_activity(context).values())
