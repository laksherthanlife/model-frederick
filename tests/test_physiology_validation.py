"""Validating the metabolic model against published yeast physiology it was not fitted to.

The chemostat series shipped with yeast-GEM is the data its own growth-associated
maintenance was fitted to, so reproducing it to 5% confirms the pipeline runs and validates
nothing. These are measurements the model never saw.

Two findings came out of it. The maintenance prediction is right to within the measurement
error of an experiment done after the parameter was set, which is as clean a held-out test
as this model gets. And the model does not ferment at all -- no oxygen bound means it will
respire without limit, so the Crabtree effect is simply absent until a measured respiratory
ceiling is imposed.
"""

import pytest

from ystwin import paths
from ystwin.bridge.thermodynamic import _default_model

pytestmark = pytest.mark.skipif(
    paths.yeast_gem() is None,
    reason="needs yeast-GEM; see docs/REPRODUCING.md",
)

GLUCOSE, OXYGEN, ETHANOL, BIOMASS = "r_1714", "r_1992", "r_1761", "r_2111"

MAX_OXYGEN_UPTAKE = 12.0
"""Highest oxygen uptake measured in S. cerevisiae, Postma 1989 PMID 2566299, CBS 8066 at
D = 0.38 /h. yeast-GEM leaves oxygen unbounded, which is why it never ferments."""


@pytest.fixture(scope="module")
def model():
    return _default_model()


def _at(model, glucose, oxygen_cap=1000.0, growth=None):
    with model:
        model.reactions.get_by_id(GLUCOSE).bounds = (-glucose, -glucose)
        model.reactions.get_by_id(OXYGEN).lower_bound = -oxygen_cap
        if growth is not None:
            model.reactions.get_by_id(BIOMASS).bounds = (growth, growth)
        solution = model.optimize()
        return {"growth": solution.objective_value,
                "ethanol": solution.fluxes[ETHANOL],
                "oxygen": -solution.fluxes[OXYGEN]}


class TestMaintenanceAgainstAnExperimentDoneAfterItWasSet:
    """Vos 2016, PMID 27317316: an aerobic retentostat held at near-zero growth measured
    glucose consumption of 0.039 +/- 0.003 mmol/gDW/h. yeast-GEM's maintenance was set in
    2016 from other work, so this is genuinely held out."""

    def test_the_model_predicts_the_measured_maintenance_demand(self, model):
        with model:
            model.reactions.get_by_id(BIOMASS).bounds = (0.0, 0.0)
            model.objective = model.problem.Objective(
                model.reactions.get_by_id(GLUCOSE).flux_expression, direction="max")
            predicted = -model.optimize().fluxes[GLUCOSE]

        assert predicted == pytest.approx(0.039, abs=0.006)

    def test_it_is_not_merely_close_to_zero(self):
        """A model predicting nothing would also look close to a small number."""
        assert 0.039 > 0.006


class TestTheModelDoesNotFermentOnItsOwn:
    """The finding worth keeping. Nothing in yeast-GEM limits respiration, so glucose is
    respired however fast it arrives and the Crabtree effect never appears."""

    def test_oxygen_uptake_is_unbounded_as_shipped(self, model):
        assert model.reactions.get_by_id(OXYGEN).lower_bound <= -1000

    def test_so_no_ethanol_appears_at_any_glucose_uptake(self, model):
        for glucose in (4.0, 8.0, 12.0):
            assert _at(model, glucose)["ethanol"] == pytest.approx(0.0, abs=1e-6)

    def test_and_respiration_runs_past_anything_ever_measured(self, model):
        """Postma's ceiling is 12; the model asks for more than double that."""
        assert _at(model, 12.0)["oxygen"] > 2 * MAX_OXYGEN_UPTAKE


class TestAMeasuredCeilingRestoresIt:
    def test_capping_oxygen_produces_ethanol(self, model):
        assert _at(model, 8.0, MAX_OXYGEN_UPTAKE)["ethanol"] > 1.0

    def test_low_glucose_still_respires_fully(self, model):
        assert _at(model, 3.0, MAX_OXYGEN_UPTAKE)["ethanol"] == pytest.approx(0.0, abs=1e-6)

    def test_the_switch_sits_between_the_two(self, model):
        below = _at(model, 5.0, MAX_OXYGEN_UPTAKE)["ethanol"]
        above = _at(model, 6.0, MAX_OXYGEN_UPTAKE)["ethanol"]

        assert below == pytest.approx(0.0, abs=1e-6)
        assert above > 0.5

    def test_the_cap_binds_once_it_switches(self, model):
        assert _at(model, 8.0, MAX_OXYGEN_UPTAKE)["oxygen"] == pytest.approx(
            MAX_OXYGEN_UPTAKE, rel=1e-3)

    def test_the_switch_is_still_at_too_high_a_growth_rate(self, model):
        """Honest about what it does not fix. van Hoek 1998 puts the critical rate near
        0.28 /h; an oxygen ceiling alone puts it near 0.45, because the real limit is
        proteome allocation rather than gas exchange -- which is what ecYeastGEM models."""
        switching = _at(model, 6.0, MAX_OXYGEN_UPTAKE)["growth"]

        assert switching > 0.35
