"""Exponential fed-batch, the vessel that replaces the chemostat.

Everything load-bearing here is TESTED rather than asserted, because the module's whole claim
is that a vessel nobody has built can stand in for one nobody has either. Three things carry
the claim and each has its own class below:

  * the setpoint is a stable fixed point for ANY monotone uptake->growth map, so the yield
    barely matters -- checked against a Crabtree-kinked map with the yield wrong by 10x;
  * the quasi-steady state is the CHEMOSTAT'S steady state, so `pathway/solve.py` needs no
    change -- checked by driving an intracellular pool to v_in/mu;
  * mu is MEASURED, not asserted from the pump -- which is the one thing a chemostat cannot
    do, and it comes with a standard error.

TWO CORRECTIONS ARE PINNED HERE because both were wrong in this file's own module first, and
simulation is what caught them. The settling algebra says an over-estimated yield costs more
generations than an under-estimated one; the simulator says the opposite happens, and the
reason is the UPTAKE CEILING rather than the algebra -- see
:class:`TestTheOpeningUptakeIsTheRealGuard`. And `vessel_full` read a volume that was never
recorded, so it returned False on every run that filled.
"""

from __future__ import annotations

import math

import pytest

from ystwin.fba.fedbatch import (
    FedBatchDesign,
    SetpointAboveCapacity,
    design_fedbatch,
    estimate_growth_rate,
    reanchor_feed,
    simulate_fedbatch,
)

# A Crabtree kink: respiratory below q_crit, fermentative above, capped. Deliberately NOT the
# linear map the fixed-point argument is easiest to see in -- the claim is about ANY monotone
# map, so testing it on a linear one would be testing the easy case.
Q_CRIT, Q_MAX = 1.1, 2.8
RESPIRATORY_YIELD, FERMENTATIVE_YIELD = 0.50, 0.088


def crabtree(uptake):
    uptake = min(uptake, Q_MAX)
    if uptake <= Q_CRIT:
        return RESPIRATORY_YIELD * uptake
    return RESPIRATORY_YIELD * Q_CRIT + FERMENTATIVE_YIELD * (uptake - Q_CRIT)


MU_MAX = crabtree(Q_MAX)


def build(mu_set=0.15, assumed_yield=0.088, max_volume=2.0, initial_biomass=0.5):
    return design_fedbatch(
        mu_set, MU_MAX,
        initial_volume_l=1.0, max_volume_l=max_volume,
        initial_biomass_g_per_l=initial_biomass, feed_substrate_g_per_l=500.0,
        assumed_yield_g_per_g=assumed_yield, max_uptake_g_per_gdcw_h=Q_MAX)


def run(design, true_yield=0.088, hours=400.0, steps=60000):
    return simulate_fedbatch(design, crabtree, true_yield_g_per_g=true_yield,
                             max_uptake_g_per_gdcw_h=Q_MAX, hours=hours, steps=steps)


class TestTheSetpointIsAStableFixedPoint:
    """The claim the replacement rests on."""

    @pytest.mark.parametrize("mu_set", [0.05, 0.10, 0.15, 0.20])
    def test_mu_converges_on_the_setpoint(self, mu_set):
        result = run(build(mu_set=mu_set))

        assert abs(result["setpoint_error_relative"]) < 0.02

    @pytest.mark.parametrize("assumed", [0.088, 0.15, 0.30, 0.88])
    def test_it_converges_even_when_the_assumed_yield_is_badly_wrong(self, assumed):
        """Up to a 10x over-estimate. The yield enters only through F0, which is the loop's
        initial condition and not its fixed point."""
        result = run(build(assumed_yield=assumed), true_yield=0.088)

        assert abs(result["setpoint_error_relative"]) < 0.02

    def test_the_map_used_here_really_is_non_linear(self):
        """Guards the guard: if `crabtree` ever became linear, every convergence test above
        would still pass while testing nothing about arbitrary monotone maps."""
        below = crabtree(0.5) / 0.5
        above = crabtree(2.0) / 2.0

        assert below > 1.5 * above

    def test_a_richer_feed_and_a_leaner_one_reach_the_same_place(self):
        lean = run(build(assumed_yield=0.88))
        near = run(build(assumed_yield=0.088))

        assert lean["settled_growth_rate_per_h"] == pytest.approx(
            near["settled_growth_rate_per_h"], rel=0.03)


class TestTheOpeningUptakeIsTheRealGuard:
    """The correction. The settling algebra points one way and the vessel behaves the other.

    ``q_S(0) = mu_set / Y_assumed`` is fixed before inoculation and does not involve the
    vessel at all, so a feed that opens above the uptake ceiling can never be rescued by a
    bigger tank or a heavier inoculum. That is what actually breaks an under-estimated yield,
    NOT the extra generations the settling law attributes to it.
    """

    def test_a_feed_that_opens_above_the_ceiling_is_refused(self):
        with pytest.raises(SetpointAboveCapacity, match="uptake ceiling"):
            build(assumed_yield=0.0088)

    def test_the_refusal_names_the_yield_that_would_be_safe(self):
        with pytest.raises(SetpointAboveCapacity, match=r"at least 0\.05357"):
            build(assumed_yield=0.0088)

    def test_the_opening_uptake_does_not_depend_on_the_vessel(self):
        """Why a bigger tank cannot rescue it: F0 scales with X0*V0 and q_S(0) divides it out."""
        small = build(max_volume=2.0, initial_biomass=0.5)
        large = build(max_volume=50.0, initial_biomass=5.0)

        assert (small.initial_uptake_g_per_gdcw_h
                == pytest.approx(large.initial_uptake_g_per_gdcw_h))

    def test_the_refused_design_really_would_have_failed(self):
        """The guard has to be justified by the behaviour, not merely by the inequality.
        Constructing the refused design directly and running it must actually miss."""
        refused = FedBatchDesign(
            mu_setpoint_per_h=0.15, initial_volume_l=1.0, max_volume_l=2.0,
            initial_biomass_g_per_l=0.5, feed_substrate_g_per_l=500.0,
            assumed_yield_g_per_g=0.0088, mu_max_per_h=MU_MAX)

        assert refused.initial_uptake_g_per_gdcw_h > Q_MAX
        assert abs(run(refused)["setpoint_error_relative"]) > 0.10

    def test_the_settling_law_alone_would_have_allowed_it(self):
        """Pins the reason this correction was needed: on generations alone the rich feed
        looks CHEAPER than the lean one, which is the opposite of how they behave."""
        design = build()
        rich = design.settling_generations(10.0)     # yield under-estimated 10x
        lean = design.settling_generations(0.1)      # yield over-estimated 10x

        assert rich < lean

    def test_the_lowest_safe_yield_is_reported(self):
        design = build()

        assert design.worst_tolerable_yield_overestimate(Q_MAX) == pytest.approx(
            0.15 / Q_MAX)


class TestItRefusesTheUpperElizondoState:
    """This repository's own calibration set contains a state fed-batch cannot reproduce."""

    def test_a_setpoint_at_mu_max_is_refused(self):
        with pytest.raises(SetpointAboveCapacity, match="margin"):
            design_fedbatch(0.254, 0.254, initial_volume_l=1.0, max_volume_l=2.0,
                            initial_biomass_g_per_l=0.5, feed_substrate_g_per_l=500.0,
                            assumed_yield_g_per_g=0.15, max_uptake_g_per_gdcw_h=2.59)

    def test_the_state_is_refused_by_the_uptake_guard_too(self):
        """Both guards independently reject it, and at the published yield the uptake one
        fires first: 0.254/0.098 = 2.59 g/gDCW/h is exactly that strain's measured uptake."""
        with pytest.raises(SetpointAboveCapacity, match="uptake ceiling"):
            design_fedbatch(0.254, 0.254, initial_volume_l=1.0, max_volume_l=2.0,
                            initial_biomass_g_per_l=0.5, feed_substrate_g_per_l=500.0,
                            assumed_yield_g_per_g=0.098, max_uptake_g_per_gdcw_h=2.59)

    def test_the_lower_elizondo_state_is_accepted(self):
        """The asymmetry that matters for the calibration set: one state transfers, one does not."""
        design = design_fedbatch(0.101, 0.254, initial_volume_l=1.0, max_volume_l=2.0,
                                 initial_biomass_g_per_l=0.5, feed_substrate_g_per_l=500.0,
                                 assumed_yield_g_per_g=0.098, max_uptake_g_per_gdcw_h=2.59)

        assert design.capacity_margin > 0.5

    def test_the_refusal_explains_itself_rather_than_just_failing(self):
        from ystwin.fba.fedbatch import REFUSES_THE_UPPER_ELIZONDO_STATE

        assert "0.254" in REFUSES_THE_UPPER_ELIZONDO_STATE
        assert "mu_max" in REFUSES_THE_UPPER_ELIZONDO_STATE


class TestThePoolsLawSurvivesTheVesselSwap:
    """`pathway/solve.py` divides by mu. If the settled fed-batch is the chemostat's steady
    state, it divides by mu_set and nothing in that module changes."""

    @pytest.mark.parametrize("mu_set", [0.05, 0.10, 0.20])
    def test_a_terminal_intracellular_pool_reaches_v_in_over_mu(self, mu_set):
        """d[P]/dt = v_in - mu*[P], integrated along the actual trajectory."""
        design = build(mu_set=mu_set)
        result = run(design)
        v_in = 0.004

        content = 0.0
        times = result["times_h"]
        for previous, current, mu in zip(times, times[1:], result["growth_rate_per_h"]):
            content += (v_in - mu * content) * (current - previous)

        assert content == pytest.approx(v_in / mu_set, rel=0.05)

    def test_the_settled_uptake_is_the_one_that_supports_the_setpoint(self):
        design = build(mu_set=0.10)
        result = run(design)

        assert crabtree(result["uptake_g_per_gdcw_h"][-1]) == pytest.approx(0.10, rel=0.02)

    def test_carbon_does_not_accumulate_once_settled(self):
        """Balanced growth: the residual is stationary, which is what the pools law needs."""
        result = run(build())

        assert result["residual_substrate_g_per_l"] < 1e-6


class TestMuIsMeasuredNotAsserted:
    """The thing a chemostat cannot do. D has no error bar because D is an input."""

    def test_it_recovers_a_known_rate_from_a_clean_series(self):
        mu = 0.12
        times = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0]
        volumes = [1.0 + 0.01 * t for t in times]
        biomass = [0.5 * math.exp(mu * t) / v for t, v in zip(times, volumes)]

        estimate = estimate_growth_rate(times, biomass, volumes)

        assert estimate.mu_per_h == pytest.approx(mu, rel=1e-9)
        assert estimate.standard_error_per_h < 1e-9

    def test_ignoring_the_volume_understates_mu(self):
        """Total biomass X*V, not concentration -- the vessel is filling, and X alone loses
        exactly the volume's own growth."""
        mu = 0.12
        times = [0.0, 4.0, 8.0, 12.0]
        volumes = [1.0 + 0.02 * t for t in times]
        biomass = [0.5 * math.exp(mu * t) / v for t, v in zip(times, volumes)]

        with_volume = estimate_growth_rate(times, biomass, volumes).mu_per_h
        without = estimate_growth_rate(times, biomass).mu_per_h

        assert with_volume == pytest.approx(mu, rel=1e-9)
        assert without < mu

    def test_noise_shows_up_as_a_standard_error(self):
        mu, noise = 0.12, [1.018, 0.982, 1.021, 0.977, 1.012, 0.990]
        times = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0]
        biomass = [0.5 * math.exp(mu * t) * n for t, n in zip(times, noise)]

        estimate = estimate_growth_rate(times, biomass)

        assert estimate.standard_error_per_h > 0
        assert estimate.relative_standard_error < 0.05

    def test_more_points_over_a_longer_span_tighten_it(self):
        mu = 0.12
        noise = [1.02, 0.98, 1.02, 0.98, 1.02, 0.98, 1.02, 0.98, 1.02]
        short = [0.0, 3.0, 6.0]
        long = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0]

        def series(times):
            return [0.5 * math.exp(mu * t) * n for t, n in zip(times, noise)]

        assert (estimate_growth_rate(long, series(long)).standard_error_per_h
                < estimate_growth_rate(short, series(short)).standard_error_per_h)

    @pytest.mark.parametrize("bad", [
        ([0.0, 1.0], [1.0, 2.0]),
        ([0.0, 1.0, 2.0], [1.0, 2.0]),
        ([0.0, 1.0, 2.0], [1.0, 0.0, 2.0]),
        ([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]),
    ])
    def test_it_refuses_a_series_it_cannot_estimate_from(self, bad):
        with pytest.raises(ValueError):
            estimate_growth_rate(*bad)

    def test_two_points_are_refused_because_they_carry_no_uncertainty(self):
        with pytest.raises(ValueError, match="standard error"):
            estimate_growth_rate([0.0, 1.0], [1.0, 2.0])


class TestSamplingIsADisturbance:
    """A chemostat is sampled from its overflow for free; this vessel is not."""

    def test_reanchoring_returns_the_feed_that_holds_the_setpoint(self):
        design = build(mu_set=0.10)
        feed = reanchor_feed(design, 2.0, 1.5)

        assert feed == pytest.approx(0.10 * 2.0 * 1.5 / (0.088 * 500.0))

    def test_removing_biomass_lowers_the_feed_the_vessel_needs(self):
        design = build()

        assert reanchor_feed(design, 1.8, 1.5) < reanchor_feed(design, 2.0, 1.5)

    def test_the_reanchored_feed_matches_the_open_loop_one_when_nothing_was_removed(self):
        """The two laws are the same law; re-anchoring only closes the loop by hand."""
        design = build()

        assert reanchor_feed(design, design.initial_biomass_g_per_l,
                             design.initial_volume_l) == pytest.approx(
            design.initial_feed_rate_l_per_h)

    def test_it_refuses_an_empty_vessel(self):
        with pytest.raises(ValueError, match="must be positive"):
            reanchor_feed(build(), 0.0, 1.5)


class TestTheGenerationBudget:
    def test_it_does_not_depend_on_the_setpoint(self):
        """The setpoint sets how fast the budget is spent, never how large it is."""
        slow = build(mu_set=0.05).generation_budget
        fast = build(mu_set=0.20).generation_budget

        assert slow == pytest.approx(fast)

    def test_a_lighter_inoculum_buys_generations(self):
        assert (build(initial_biomass=0.05).generation_budget
                > build(initial_biomass=0.5).generation_budget + 3.0)

    def test_it_is_computed_at_the_true_yield_not_the_assumed_one(self):
        """The correction: the budget is a physical quantity. Reading it off the assumed
        yield makes it appear to shrink when you under-estimate, which is backwards."""
        design = build(assumed_yield=0.30)

        assert design.generation_budget_at(0.088) < design.generation_budget

    def test_the_final_biomass_is_independent_of_the_inoculum(self):
        assert (build(initial_biomass=0.05).final_biomass_g_per_l
                == pytest.approx(build(initial_biomass=0.5).final_biomass_g_per_l))

    def test_settling_is_logarithmic_in_the_yield_error(self):
        """Which is why a factor-of-two error is affordable and a factor of ten still is."""
        design = build()
        doubled = design.settling_generations(2.0)
        tenfold = design.settling_generations(10.0)

        assert tenfold < 2 * doubled

    def test_a_perfect_yield_costs_no_settling(self):
        assert build().settling_generations(1.0) == 0.0


class TestTheSimulatorReportsWhatHappened:
    def test_a_run_that_fills_the_vessel_says_so(self):
        """`vessel_full` was read off a volume the loop never recorded, so it was False on
        every run that filled."""
        result = run(build(max_volume=2.0))

        assert result["vessel_full"]
        assert result["volume_l"][-1] >= 2.0

    def test_a_run_that_stops_early_does_not_claim_to_have_filled(self):
        result = run(build(max_volume=50.0), hours=5.0, steps=2000)

        assert not result["vessel_full"]

    def test_volume_and_biomass_both_grow(self):
        result = run(build())

        assert result["volume_l"][-1] > result["volume_l"][0]
        assert result["biomass_g_per_l"][-1] > result["biomass_g_per_l"][0]

    @pytest.mark.parametrize("kw", [{"hours": 0.0}, {"steps": 1}])
    def test_it_refuses_a_degenerate_run(self, kw):
        with pytest.raises(ValueError):
            simulate_fedbatch(build(), crabtree, true_yield_g_per_g=0.088,
                              max_uptake_g_per_gdcw_h=Q_MAX,
                              **{"hours": 10.0, "steps": 100, **kw})

    def test_it_refuses_a_missing_uptake_ceiling(self):
        with pytest.raises(ValueError, match="max_uptake"):
            simulate_fedbatch(build(), crabtree, true_yield_g_per_g=0.088,
                              max_uptake_g_per_gdcw_h=0.0, hours=10.0, steps=100)


class TestDesignRefusals:
    @pytest.mark.parametrize("kw,match", [
        ({"mu_set": 0.0}, "mu_setpoint must be positive"),
        ({"max_volume": 1.0}, "must exceed"),
        ({"initial_biomass": 0.0}, "initial_biomass_g_per_l must be positive"),
    ])
    def test_a_contradictory_vessel_is_refused(self, kw, match):
        with pytest.raises(ValueError, match=match):
            build(**kw)

    def test_the_margin_default_is_declared_a_judgement(self):
        """It is 10% because the convergence rate falls off smoothly and there is no sharp
        threshold to read off -- so it must not be presented as measured."""
        import inspect

        assert "judgement, not a measurement" in inspect.getdoc(design_fedbatch)
