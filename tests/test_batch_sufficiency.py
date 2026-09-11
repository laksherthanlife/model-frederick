"""The chemostat requirement, tested and found unnecessary.

Every recommendation in this repository used to assume a HELD growth rate, and
`docs/DISTANCE_TO_THE_VISION.md` listed "shake flask / batch" as REFUSED. An adversarial
audit attacked that premise on 2026-09-03 and it did not survive: the balance
`pathway/solve.py` closes is per gram of biomass, so the vessel does not appear in it.

These tests pin the finding and the two corrections it forced -- one to a table that was
simply false, and one to a protocol that presented a fed-batch needing a sparged fermenter
as the low-equipment option.
"""

from __future__ import annotations

import math

import pytest

from ystwin.generator.context import CultureContext
from ystwin.pathway import calibrations
from ystwin.pathway.solve import solve_pathway
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, predict_product

KINETICS = calibrations.BETA_CAROTENE_KINETICS
V_IN = calibrations.BETA_CAROTENE_FLUX.alpha * 1.0


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


def _integrate(spec, growth_rate, generations, step=1e-4):
    """Forward-integrate the same balance solve_pathway closes, from an empty cell."""
    branch = KINETICS["lycopene"]
    vmax, km = branch.vmax_per_growth * growth_rate, branch.km
    intermediate = terminal = 0.0
    hours = generations * math.log(2.0) / growth_rate
    elapsed = 0.0
    while elapsed < hours:
        export = vmax * intermediate / (km + intermediate)
        intermediate += (V_IN - export - growth_rate * intermediate) * step
        terminal += (export - growth_rate * terminal) * step
        elapsed += step
    return terminal


class TestBatchReachesTheChemostatAnswer:
    """The finding. A flask at constant mu converges on the steady state exactly."""

    @pytest.mark.parametrize("growth_rate", [0.101, 0.18, 0.2543])
    def test_it_converges_on_the_steady_state_solver(self, spec, growth_rate):
        steady = solve_pathway(spec, V_IN, growth_rate,
                               KINETICS).terminal.content_mmol_per_gdcw
        integrated = _integrate(spec, growth_rate, generations=12.0)

        assert integrated == pytest.approx(steady, rel=1e-3)

    @pytest.mark.parametrize("growth_rate", [0.101, 0.18, 0.2543])
    def test_five_generations_gets_within_five_percent(self, spec, growth_rate):
        """The number that makes a flask experiment viable: about 33 h at mu = 0.101."""
        steady = solve_pathway(spec, V_IN, growth_rate,
                               KINETICS).terminal.content_mmol_per_gdcw
        integrated = _integrate(spec, growth_rate, generations=5.1)

        assert abs(integrated - steady) / steady < 0.05

    def test_two_generations_is_not_enough(self, spec):
        """The other half, and the reason the inoculum is the whole experiment: a heavy
        inoculum leaves only ~2.6 generations of exponential phase and lands ~36% low."""
        steady = solve_pathway(spec, V_IN, 0.18, KINETICS).terminal.content_mmol_per_gdcw
        integrated = _integrate(spec, 0.18, generations=2.0)

        assert abs(integrated - steady) / steady > 0.25

    def test_the_generation_count_barely_moves_with_growth_rate(self, spec):
        """Which is why it is quotable as a protocol constant rather than a lookup."""
        def generations_to_five_percent(growth_rate):
            steady = solve_pathway(spec, V_IN, growth_rate,
                                   KINETICS).terminal.content_mmol_per_gdcw
            low, high = 1.0, 12.0
            while high - low > 0.02:
                mid = (low + high) / 2
                if abs(_integrate(spec, growth_rate, mid) - steady) / steady < 0.05:
                    high = mid
                else:
                    low = mid
            return high

        slow, fast = generations_to_five_percent(0.101), generations_to_five_percent(0.2543)

        assert abs(slow - fast) < 0.5
        assert 4.5 < slow < 5.5


class TestTheFlaskRefusalWasAboutGrowthRateNotTheVessel:
    """`DISTANCE_TO_THE_VISION.md` listed batch culture as refused. It is not."""

    @pytest.mark.parametrize("carbon,temperature", [
        ("glucose", 20.0), ("glucose", 15.0), ("galactose", 25.0),
        ("galactose", 20.0), ("ethanol", 30.0), ("ethanol", 25.0)])
    def test_these_flask_contexts_answer(self, spec, carbon, temperature):
        got = predict_product(
            spec, Genotype(entry_expression=1.0),
            Environment(context=CultureContext(carbon_source=carbon,
                                               temperature_c=temperature)),
            calibrations.BETA_CAROTENE_FLUX, KINETICS)

        assert got.content_mg_per_gdcw > 0

    def test_the_default_flask_context_is_refused_for_its_growth_rate(self, spec):
        """glucose at 30 C gives mu = 0.40, outside the fitted [0.101, 0.2543] window. That
        is a statement about mu, not about batch culture."""
        low, high = KINETICS["lycopene"].growth_rate_range

        with pytest.raises(ValueError):
            predict_product(
                spec, Genotype(entry_expression=1.0),
                Environment(context=CultureContext(carbon_source="glucose",
                                                   temperature_c=30.0)),
                calibrations.BETA_CAROTENE_FLUX, KINETICS)
        assert high < 0.40


class TestTheProtocolCorrectionsAreRecorded:
    def test_p4_exists_and_needs_no_chemostat(self):
        import pathlib

        text = (pathlib.Path(__file__).resolve().parents[1]
                / "docs" / "PROTOCOLS.md").read_text()

        assert "## P4" in text
        assert "does not have to be HELD" in text

    def test_p3_no_longer_claims_it_reproduces_the_calibration_states(self):
        import pathlib

        text = (pathlib.Path(__file__).resolve().parents[1]
                / "docs" / "PROTOCOLS.md").read_text()

        assert "cannot reproduce EITHER Elizondo calibration state" in text
        assert "glucose-**EXCESS**" in text
