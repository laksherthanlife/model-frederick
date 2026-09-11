"""The first stress channel large enough to reach the product.

Six earlier routes from a stress state into the GEM were measured and all failed -- NGAM
maintenance, GECKO protein pool, precursor ceilings, NADPH competition, Crabtree
partitioning, cassette burden. The largest moved titre 0.17-1.73% against an entry-flux law
whose own held-out error is 22.2%. Every one was a SCALAR TAX.

This one is a different shape: flux REQUIRED into the named metabolic effectors of the stress
response. Yeast9 carries those effectors and none of the signalling that drives them, so the
GEM has the machinery and no reason to use it -- FBA maximises growth and trehalose costs
carbon. That missing decision is what a learned latent state supplies.

Growth is HELD throughout, which is what separates this from every earlier attempt: those all
reached the product through mu and were indistinguishable from the growth-rate response the
stress panel already measures directly.

Marked `integration`: every test solves the real GSMM.
"""

from __future__ import annotations

import pytest

from ystwin.bridge.stress_diversion import (
    STRESS_EFFECTORS,
    EffectorUnmeasured,
    diversion_cost,
)
from ystwin.pathway.spec import load_pathway

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def carotene():
    return load_pathway("beta_carotene")


class TestTheChannelIsLargeEnoughToMatter:
    """The finding. Against the 22.2% floor every earlier channel sat far below."""

    @pytest.mark.parametrize("effector", ["trehalose", "glutathione"])
    def test_one_mmol_clears_the_calibration_floor(self, yeast_gem, carotene, effector):
        result = diversion_cost(yeast_gem, carotene, effector, 1.0, 0.18)

        assert result.clears_the_calibration_floor
        assert result.relative_change < -0.222

    def test_it_is_an_order_of_magnitude_above_the_ngam_channel(self, yeast_gem, carotene):
        """NGAM moved titre 0.17-1.73%. This moves precursor supply by more than 10% at half
        that flux, which is what makes it the first stress route worth wiring."""
        result = diversion_cost(yeast_gem, carotene, "trehalose", 0.5, 0.18)

        assert abs(result.relative_change) > 0.10

    def test_the_cost_rises_with_the_flux(self, yeast_gem, carotene):
        light = diversion_cost(yeast_gem, carotene, "glutathione", 0.5, 0.18)
        heavy = diversion_cost(yeast_gem, carotene, "glutathione", 2.0, 0.18)

        assert heavy.relative_change < light.relative_change

    def test_glycerol_is_the_cheap_one(self, yeast_gem, carotene):
        """Not every effector is expensive, and the spread is the useful part: an osmotic
        response costs a quarter of what an oxidative one does at the same flux."""
        glycerol = diversion_cost(yeast_gem, carotene, "glycerol", 1.0, 0.18)
        glutathione = diversion_cost(yeast_gem, carotene, "glutathione", 1.0, 0.18)

        assert abs(glycerol.relative_change) < abs(glutathione.relative_change) / 2


class TestGrowthIsHeldSoThisIsNotTheMuChannel:
    def test_the_growth_rate_is_recorded_and_respected(self, yeast_gem, carotene):
        result = diversion_cost(yeast_gem, carotene, "trehalose", 1.0, 0.12)

        assert result.growth_rate_per_h == 0.12

    def test_the_baseline_is_measured_at_the_same_growth_rate(self, yeast_gem, carotene):
        """So the reported change is diversion alone. If the baseline moved with the
        stress this would silently become the growth channel again."""
        first = diversion_cost(yeast_gem, carotene, "trehalose", 0.5, 0.18)
        second = diversion_cost(yeast_gem, carotene, "glutathione", 2.0, 0.18)

        assert first.baseline_supply_mmol_per_gdcw_h == pytest.approx(
            second.baseline_supply_mmol_per_gdcw_h)


class TestItRefusesToInventTheStressState:
    def test_a_zero_flux_is_refused(self, yeast_gem, carotene):
        """There is no default, because the flux IS the learned latent state and this
        repository has not measured it -- latent_bridge's four constructs all returned
        INCONCLUSIVE at G4."""
        with pytest.raises(EffectorUnmeasured, match="learned latent stress state"):
            diversion_cost(yeast_gem, carotene, "trehalose", 0.0, 0.18)

    def test_an_unknown_effector_names_the_ones_that_exist(self, yeast_gem, carotene):
        with pytest.raises(KeyError, match="no effector"):
            diversion_cost(yeast_gem, carotene, "heat_shock_protein", 1.0, 0.18)

    def test_the_refusal_explains_that_the_gem_has_no_regulators(self, yeast_gem, carotene):
        with pytest.raises(KeyError, match="not in a metabolic model"):
            diversion_cost(yeast_gem, carotene, "hsp104", 1.0, 0.18)

    def test_every_effector_id_resolves_in_the_model(self, yeast_gem):
        """A wrong metabolite id would produce a confident number for the wrong chemistry."""
        ids = {metabolite.id for metabolite in yeast_gem.metabolites}

        for name, identifier in STRESS_EFFECTORS.items():
            assert identifier in ids, f"{name} -> {identifier} absent from the model"


class TestItWorksForAProductWithNoCalibration:
    def test_phb_gets_a_cost_despite_having_no_flux_law(self, yeast_gem):
        """It reads only precursor_metabolite and precursor_stoichiometry, so a pathway with
        no calibration at all still gets a diversion cost. That is what makes the transfer
        function a host property rather than a per-product fit."""
        result = diversion_cost(yeast_gem, load_pathway("phb"), "trehalose", 1.0, 0.18)

        assert result.baseline_supply_mmol_per_gdcw_h > 0
        assert result.relative_change < 0


class TestTheCostsAdd:
    """The measurement that makes a multi-module stress state computable at all.

    For each pair of effectors the joint infeasibility point along the ray (t*W_a, t*W_b) was
    found. A SHARED budget predicts t* = 0.5; INDEPENDENT budgets predict 1.0. All 45 pairs
    came back 0.5000-0.5137, and ten effectors at once gave sum(f/W) = 1.011 against 1.0
    shared and 10.0 independent. One budget, to 1.1%.
    """

    def test_two_effectors_cost_about_the_sum_of_their_parts(self):
        from ystwin.bridge.stress_diversion import (
            INFEASIBILITY_WALL_AT_MU_018,
            combined_cost,
        )

        half_a = 0.25 * INFEASIBILITY_WALL_AT_MU_018["trehalose"]
        half_b = 0.25 * INFEASIBILITY_WALL_AT_MU_018["glycerol"]
        solo_a = combined_cost({"trehalose": half_a})
        solo_b = combined_cost({"glycerol": half_b})
        joint = combined_cost({"trehalose": half_a, "glycerol": half_b})

        assert joint == pytest.approx(solo_a + solo_b, rel=1e-9)

    def test_sharing_an_intermediate_does_not_make_them_cheaper(self):
        """The obvious hypothesis, and it is false. Trehalose and glycogen both come off
        UDP-glucose and still share the budget exactly rather than sub-adding."""
        from ystwin.bridge.stress_diversion import combined_cost

        joint = combined_cost({"trehalose": 0.5, "glycogen": 0.5})
        separate = (combined_cost({"trehalose": 0.5})
                    + combined_cost({"glycogen": 0.5}))

        assert joint == pytest.approx(separate, rel=1e-9)

    def test_a_bundle_past_the_budget_is_refused_with_the_reason(self):
        """And the reason matters: the boundary is the HELD GROWTH RATE running out of
        carbon, not a metabolic ceiling. Lower mu and the same bundle becomes affordable."""
        from ystwin.bridge.stress_diversion import combined_cost

        with pytest.raises(ValueError, match="not a metabolic ceiling"):
            combined_cost({"trehalose": 5.0})

    def test_the_budget_widens_as_growth_slows(self):
        from ystwin.bridge.stress_diversion import combined_cost

        fast = combined_cost({"glycerol": 2.0}, 0.30)
        slow = combined_cost({"glycerol": 2.0}, 0.05)

        assert abs(slow) < abs(fast)

    def test_a_growth_rate_past_the_budget_zero_is_refused(self):
        from ystwin.bridge.stress_diversion import (
            BUDGET_ZERO_GROWTH_RATE,
            combined_cost,
        )

        with pytest.raises(ValueError, match="no diversion budget"):
            combined_cost({"trehalose": 0.1}, BUDGET_ZERO_GROWTH_RATE + 0.01)

    def test_the_most_expensive_real_stress_still_misses_the_floor(self):
        """Acetic acid is the panel's costliest state once the sign correction is applied,
        and at the stressed growth rate the panel itself predicts it reaches about -15%
        against a 22.2% floor. Close, and still not checkable."""
        from ystwin.bridge.stress_diversion import combined_cost

        cost = abs(combined_cost({"atp": 17.98}, 0.327))

        assert 0.10 < cost < 0.222

    def test_an_effector_with_no_measured_wall_is_refused(self):
        from ystwin.bridge.stress_diversion import combined_cost

        with pytest.raises(KeyError, match="no measured wall"):
            combined_cost({"heat_shock_protein": 1.0})


class TestTheChannelCollapsesAtMeasuredFlux:
    """The end of this line, and it is a measurement rather than an argument.

    The module's contract says the diversion flux must be supplied from a measurement. A
    search of ~53 papers found exactly one aerobic dataset that measures a stress effector in
    absolute units on one strain across graded conditions -- Hakkaart 2020, PMID 31654410,
    aerobic glucose-limited chemostat at D = 0.025 /h, glycogen and trehalose in mg/gDCW.
    Supplying it gives -0.1%, not the -15% the module first advertised.
    """

    #: q = mu * content at D = 0.025 /h, from the paper's Table, mmol/gDCW/h.
    HAKKAART = {
        "pH 5 reference": {"glycogen": 0.00545, "trehalose": 0.00142},
        "pH 3": {"glycogen": 0.00716, "trehalose": 0.00092},
        "pH 3 + 50% CO2": {"glycogen": 0.00469, "trehalose": 0.00071},
    }

    @pytest.mark.parametrize("condition", list(HAKKAART))
    def test_every_measured_condition_costs_about_a_tenth_of_a_percent(self, condition):
        from ystwin.bridge.stress_diversion import combined_cost

        cost = abs(combined_cost(self.HAKKAART[condition], 0.025))

        assert cost < 0.002

    def test_the_stress_increment_is_thousands_of_times_below_the_floor(self):
        """0.008 percentage points, against a 22.2% calibration floor."""
        from ystwin.bridge.stress_diversion import combined_cost

        reference = combined_cost(self.HAKKAART["pH 5 reference"], 0.025)
        stressed = combined_cost(self.HAKKAART["pH 3"], 0.025)
        increment = abs(stressed - reference)

        assert increment < 0.0002
        assert 0.222 / increment > 1000

    def test_the_module_records_the_retraction_of_the_glycerol_number(self):
        """1.49 mmol/gDCW/h is close to the UNSTRESSED ANAEROBIC baseline, and aerobically
        q_glycerol is 0.00 +/- 0.00 at every oxygen level (PMID 18613954). The largest
        surviving diversion does not exist under the conditions this repo predicts in."""
        import ystwin.bridge.stress_diversion as module

        assert "no\nglycerol flux to divert at all" in module.__doc__

    def test_the_module_records_that_the_sign_is_the_other_way(self):
        """Every published dataset that dosed a stressor and measured absolute content found
        it RISING, and where mechanism was measured it is transcriptional."""
        import ystwin.bridge.stress_diversion as module

        assert "models the wrong mechanism" in module.__doc__
        assert "TRANSCRIPTIONAL" in module.__doc__
