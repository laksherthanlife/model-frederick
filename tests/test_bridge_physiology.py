"""Branch A: measured physiology to FBA constraints to flux. Load-bearing.

This is the uncontroversial half of the plan's g_psi. Nothing latent enters it: the
inputs are quantities the plate actually measures -- biomass and growth rate -- and
the constraint they set is the one dynamic FBA has always used. Fluorescence is not
involved, which is precisely why this branch can carry weight while Branch B cannot.

Every prediction carries `branch="physiology"` so a report can attribute each number
to the evidence that produced it.
"""

import numpy as np
import pytest

from ystwin.bridge.physiology_bridge import (
    MeasuredState,
    physiology_constraints,
    predict_flux,
)

pytestmark = pytest.mark.integration


def _state(**over):
    base = dict(biomass_gl=0.35, growth_rate=0.28, glucose_mM=20.0, time_h=2.0)
    base.update(over)
    return MeasuredState(**base)


class TestConstraints:
    def test_the_measured_growth_rate_becomes_a_growth_constraint(self):
        constraints = physiology_constraints(_state(growth_rate=0.28))

        assert constraints.growth_lower_bound == pytest.approx(0.28, rel=1e-9)

    def test_the_uptake_bound_does_not_depend_on_the_measured_growth_rate(self):
        """Otherwise the model could never call a measurement impossible.

        A bound derived from demand hands the solver exactly as much substrate as
        the observed rate needs, so the feasibility check always passes and the
        flux prediction restates the measurement instead of testing it.
        """
        slow = physiology_constraints(_state(growth_rate=0.10))
        fast = physiology_constraints(_state(growth_rate=0.35))

        assert fast.glucose_uptake == pytest.approx(slow.glucose_uptake)

    def test_substrate_uptake_saturates_as_the_medium_is_replete(self):
        lean = physiology_constraints(_state(glucose_mM=0.05))
        replete = physiology_constraints(_state(glucose_mM=100.0))

        assert lean.glucose_uptake < replete.glucose_uptake

    def test_a_negative_growth_rate_is_carried_not_clipped(self):
        """A dying culture is a real state; silently flooring it hides death."""
        constraints = physiology_constraints(_state(growth_rate=-0.05))

        assert constraints.growth_lower_bound < 0

    def test_the_constraints_record_which_branch_produced_them(self):
        assert physiology_constraints(_state()).branch == "physiology"

    def test_nothing_latent_appears_in_the_inputs(self):
        """The guarantee that lets this branch carry weight."""
        assert not hasattr(MeasuredState(0.3, 0.2, 20.0, 1.0), "promoter_activity")


class TestFluxPrediction:
    def test_it_returns_a_flux_for_every_requested_reaction(self, yeast_gem):
        prediction = predict_flux(
            yeast_gem, physiology_constraints(_state()),
            reactions=["r_2111", "r_1761"],
        )

        assert set(prediction.fluxes) == {"r_2111", "r_1761"}

    def test_the_solved_growth_matches_the_measured_growth(self, yeast_gem):
        prediction = predict_flux(
            yeast_gem, physiology_constraints(_state(growth_rate=0.25)),
            reactions=["r_2111"],
        )

        assert prediction.fluxes["r_2111"] == pytest.approx(0.25, rel=0.05)

    def test_a_growth_rate_the_model_cannot_support_is_refused_by_name(self, yeast_gem):
        with pytest.raises(ValueError, match="cannot support"):
            predict_flux(
                yeast_gem, physiology_constraints(_state(growth_rate=5.0)),
                reactions=["r_2111"],
            )

    def test_the_prediction_carries_its_branch_and_its_inputs(self, yeast_gem):
        prediction = predict_flux(
            yeast_gem, physiology_constraints(_state()), reactions=["r_2111"]
        )

        assert prediction.branch == "physiology"
        assert prediction.constraints.growth_lower_bound == pytest.approx(0.28)

    def test_it_reports_the_range_not_only_a_point(self, yeast_gem):
        prediction = predict_flux(
            yeast_gem, physiology_constraints(_state()),
            reactions=["r_1761"], with_ranges=True,
        )

        low, high = prediction.ranges["r_1761"]
        assert low <= prediction.fluxes["r_1761"] <= high

    def test_the_source_model_is_not_mutated(self, yeast_gem):
        before = yeast_gem.reactions.get_by_id("r_1714").bounds
        predict_flux(yeast_gem, physiology_constraints(_state()), reactions=["r_2111"])

        assert yeast_gem.reactions.get_by_id("r_1714").bounds == before


def test_a_measured_trajectory_yields_a_flux_per_timepoint(yeast_gem):
    from ystwin.bridge.physiology_bridge import predict_flux_trajectory

    times = np.linspace(0, 4, 5)
    states = [
        MeasuredState(biomass_gl=0.1 * np.exp(0.25 * t), growth_rate=0.25,
                      glucose_mM=20.0, time_h=float(t))
        for t in times
    ]

    frame = predict_flux_trajectory(yeast_gem, states, reactions=["r_2111", "r_1761"])

    assert len(frame) == len(times)
    assert (frame.branch == "physiology").all()
