"""The two-population system, the unit bridge, and the fed-batch ethanol diagnostic.

Three of these classes are calibrations rather than unit tests and say so in their
docstrings. ``TestCriterionA`` re-runs both ablations on the real committed rows and pins the
duration switch that decides whether this block is worth its states;
``TestTheTimescaleAuditDecidesTheVerdict`` re-runs `mech/integrate.py`'s tau/T audit and pins
the four verdicts; and ``TestTheFedBatchCollapse`` runs the system over the vessel
`scripts/design_fedbatch_run.py` actually publishes and pins what segregation costs it.

Everything numeric here is either recomputed from a committed file or derived from a
conversion in the module under test. Nothing is transcribed.
"""

from __future__ import annotations

import csv
import math

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.fba.fedbatch import design_fedbatch
from ystwin.generator.context import CultureContext
from ystwin.generator.plate import PlateConditions
from ystwin.mech import population as pop
from ystwin.mech.ablation import GROWTH_RATE_FLOOR, REPORTER_ACTIVITY_FLOOR
from ystwin.mech.integrate import Verdict, reduce_for_window
from ystwin.mech.params import (
    FREE_TAGS,
    FreeScalarGateFailed,
    Param,
    RefusedValue,
    Tag,
)
from ystwin.mech.state import FEDBATCH_5D, PLATE_READ_4H
from ystwin.predict import LayerState

LN2 = math.log(2.0)


@pytest.fixture(scope="module")
def growth_rows():
    """The 183 committed wells, read once for the whole module."""
    return pop.plate_chord_growth_rates()


@pytest.fixture(scope="module")
def plate_mu(growth_rows):
    """The MEASURED plate rate the reporter channel converts p with."""
    return float(np.median(growth_rows.growth_rate))


class TestCriterionEIsTwoAgainstTwo:
    """The gate that refused this module at import until the miscount was fixed."""

    def test_the_built_system_has_exactly_two_free_scalars(self):
        free = {p.name for p in pop.POPULATION_PARAMS.free_scalars()}
        assert free == {"p_loss_per_division", "b_burden"}

    def test_it_is_scored_against_two_independent_targets(self):
        targets = {t.name for t in pop.POPULATION_PARAMS.independent_targets()}
        assert targets == {GROWTH_RATE_FLOOR.name, REPORTER_ACTIVITY_FLOOR.name}
        assert pop.POPULATION_GATE.passes

    def test_the_expressing_band_is_the_same_scalar_and_sits_outside_the_gate(self):
        assert "p_loss_per_division_expressing" not in pop.POPULATION_PARAMS
        assert "p_loss_per_division_expressing" in pop.BAND_VARIANTS
        assert not pop.BAND_VARIANTS.targets

    def test_counting_it_twice_is_what_the_gate_refuses(self):
        registry = type(pop.POPULATION_PARAMS)("probe")
        registry.add(pop.LOSS_PER_DIVISION)
        registry.add(pop.BURDEN)
        registry.add(pop.LOSS_PER_DIVISION_EXPRESSING)
        registry.add_target(GROWTH_RATE_FLOOR)
        registry.add_target(REPORTER_ACTIVITY_FLOOR)
        with pytest.raises(FreeScalarGateFailed, match="3 free scalars against 2"):
            registry.require_gate()

    def test_the_diagnostic_does_not_pass_its_own_gate_and_that_is_why_it_is_reported(self):
        assert not pop.ETHANOL_GATE.passes
        assert pop.ETHANOL_GATE.free == ("r_ethanol_per_biomass",)
        assert pop.ETHANOL_GATE.targets == ()
        assert pop.REPORTED == LayerState.REPORTED

    def test_every_refusal_carries_a_reason_and_the_measurement_that_would_close_it(self):
        refusals = pop.NOT_BUILT.refusals()
        assert len(refusals) == len(pop.NOT_BUILT)
        for param in refusals:
            assert param.tag == Tag.REFUSED
            assert param.reason.strip() and param.missing.strip()
            with pytest.raises(RefusedValue):
                float(param)

    def test_the_refusals_are_outside_both_gates(self):
        assert not (set(pop.NOT_BUILT._params) & set(pop.POPULATION_PARAMS._params))
        assert not (set(pop.NOT_BUILT._params) & set(pop.ETHANOL_PARAMS._params))


class TestTheUnitBridge:
    """theta = 1 - 2^-p. The 44% the module exists to stop, and the convergence it buys."""

    @pytest.mark.parametrize("theta", [0.010, 0.0175, 0.057, 0.285, 0.5])
    def test_the_two_conversions_are_inverses(self, theta):
        assert pop.loss_per_generation(pop.divisions_per_generation(theta)) == pytest.approx(
            theta, rel=1e-12)

    def test_they_differ_by_exactly_one_over_ln2_in_the_small_limit(self):
        theta = 1e-6
        assert pop.divisions_per_generation(theta) / theta == pytest.approx(1 / LN2, rel=1e-5)

    def test_the_swept_band_is_the_measured_band_converted_and_not_transcribed(self):
        low = pop.divisions_per_generation(min(pop.MEASURED_LOSS_PER_GENERATION))
        high = pop.divisions_per_generation(max(pop.MEASURED_LOSS_PER_GENERATION))
        assert pop.LOSS_PER_DIVISION.bounds == pytest.approx((low, high), rel=1e-12)

    def test_the_expressing_band_is_the_2018_multiplier_applied_to_the_same_ends(self):
        low = pop.divisions_per_generation(0.020)
        high = pop.divisions_per_generation(0.285)
        assert pop.LOSS_PER_DIVISION_EXPRESSING.bounds == pytest.approx((low, high), rel=1e-12)

    def test_chengs_fitted_p_lands_inside_hohnholz_only_after_the_conversion(self):
        result = pop.cheng_hohnholz_convergence()
        assert result["converted_lands_inside"] == 1.0
        assert result["unconverted_lands_inside"] == 0.0
        assert result["cheng_theta_per_generation"] == pytest.approx(0.0539423, rel=1e-5)

    def test_the_convergence_survives_the_discrete_reading_of_chengs_kernel(self):
        """The verifier's second ambiguity in the same symbol, and it changes no verdict."""
        result = pop.cheng_hohnholz_convergence()
        assert result["cheng_theta_discrete_kernel"] == pytest.approx(0.04, rel=1e-12)
        assert result["discrete_reading_lands_inside"] == 1.0
        assert result["converted_lands_inside"] == 1.0
        assert result["unconverted_lands_inside"] == 0.0

    def test_the_two_kernel_readings_differ_by_the_priced_39_percent(self):
        result = pop.cheng_hohnholz_convergence()
        assert result["kernel_readings_differ_by"] == pytest.approx(0.34856, rel=1e-4)
        small = pop.loss_per_generation(1e-9) / (0.5 * 1e-9) - 1.0
        assert small == pytest.approx(2 * LN2 - 1.0, rel=1e-6)

    def test_the_ambiguity_cannot_reach_a_constant_this_module_builds(self):
        """p is defined by this module's ODE, so the swept band is its own inverse and not Cheng's."""
        for theta in pop.MEASURED_LOSS_PER_GENERATION:
            assert pop.loss_per_generation(
                pop.divisions_per_generation(theta)) == pytest.approx(theta, rel=1e-12)

    def test_a_negative_p_and_a_theta_of_one_are_both_refused(self):
        with pytest.raises(ValueError):
            pop.loss_per_generation(-0.01)
        with pytest.raises(ValueError):
            pop.divisions_per_generation(1.0)


class TestThePartitionKernels:
    """A microscope predicting an agar plate, and the bacterial form that fails by 1e19."""

    def test_the_foci_prediction_overlaps_the_replica_plated_band(self):
        result = pop.foci_prediction_overlap()
        assert result["overlaps"] == 1.0
        assert result["overlap_low"] < result["overlap_high"]
        assert result["predicted_high"] == pytest.approx(2.0 ** -6)

    def test_the_kernel_returns_the_per_daughter_fraction_replica_plating_estimates(self):
        """theta = 1 - (F/I)^(1/N) counts plasmid-free DAUGHTERS, not divisions that make one.

        The per-division form is twice this -- the same factor of 2 the discrete reading of
        Cheng's kernel carries -- and applying it here would double every prediction.
        """
        for foci in (1, 2, 3, 5):
            copies = 2 * foci
            per_daughter = math.comb(copies, 0) / 2.0 ** copies
            assert pop.foci_partition_loss_per_generation(foci) == pytest.approx(
                per_daughter, rel=1e-12)
            assert 2.0 * per_daughter == pytest.approx(2.0 ** (1 - copies), rel=1e-12)

    def test_the_kernel_is_in_this_module_s_own_theta_units(self):
        theta = pop.foci_partition_loss_per_generation(3)
        system = pop.TwoPopulation(pop.divisions_per_generation(theta), 0.0)
        assert system.fraction_after_generations(1.0, 1.0) == pytest.approx(
            1.0 - theta, rel=1e-12)

    def test_the_molecule_kernel_is_refuted_by_seven_to_nineteen_orders_of_magnitude(self):
        measured_low = min(pop.MEASURED_LOSS_PER_GENERATION)
        for copies, orders in ((14, 6), (34, 18)):
            predicted = pop.binomial_partition_loss_per_generation(copies)
            assert math.log10(measured_low / predicted) > orders

    def test_a_segregating_unit_count_below_one_is_refused(self):
        with pytest.raises(ValueError):
            pop.foci_partition_loss_per_generation(0)


class TestTheTwoPopulationSystem:
    """dF/dt, and the generation invariance that is the layer's most useful statement."""

    def test_a_probability_outside_zero_one_is_refused_at_construction(self):
        with pytest.raises(ValueError, match="probability per division"):
            pop.TwoPopulation(1.0, 0.1)
        with pytest.raises(ValueError, match="lethal construct"):
            pop.TwoPopulation(0.05, 1.0)

    def test_the_bearing_fraction_only_ever_falls(self):
        system = pop.TwoPopulation(0.05, 0.1)
        for fraction in (0.05, 0.5, 0.95, 1.0):
            assert system.d_fraction_d_time(fraction, 0.3) < 0
        assert system.d_fraction_d_time(0.0, 0.3) == 0.0

    def test_the_generation_derivative_carries_no_growth_rate(self):
        system = pop.TwoPopulation(0.05, 0.15)
        rates = [system.fraction_after_generations(0.95, 5.0) for _ in range(3)]
        assert rates[0] == rates[1] == rates[2]
        for mu in (0.05, 0.101, 0.363):
            assert system.d_fraction_d_time(0.8, mu) / mu == pytest.approx(
                system.d_fraction_d_generation(0.8) * system.population_growth_rate(0.8, 1.0)
                / LN2, rel=1e-12)

    def test_at_zero_burden_the_closed_form_is_the_per_generation_definition(self):
        system = pop.TwoPopulation(0.05, 0.0)
        assert system.fraction_after_generations(0.95, 10.0) == pytest.approx(
            0.95 * 2.0 ** -0.5, rel=1e-12)
        assert system.loss_per_generation == pytest.approx(1 - 2.0 ** -0.05, rel=1e-12)

    def test_generations_to_a_target_inverts_the_trajectory(self):
        for burden in (0.0, 0.15):
            system = pop.TwoPopulation(0.05, burden)
            generations = system.generations_to_fraction(0.95, 0.5)
            assert system.fraction_after_generations(0.95, generations) == pytest.approx(
                0.5, rel=1e-7)

    def test_the_ablated_model_never_reaches_a_target_and_says_so(self):
        with pytest.raises(ValueError, match="ablated model"):
            pop.TwoPopulation(0.0, 0.0).generations_to_fraction(0.95, 0.5)

    def test_burden_makes_the_population_rate_rise_as_the_plasmid_is_lost(self):
        system = pop.TwoPopulation(0.05, 0.2)
        assert system.population_growth_rate(0.2, 0.3) > system.population_growth_rate(1.0, 0.3)
        assert system.population_growth_rate(0.0, 0.3) == pytest.approx(0.3)


class TestTheCommittedRowsCannotSeparateTheTwoModels:
    """A calibration, and the sharpest statement of what the 4.14 h caveat costs.

    The criterion-(a) pass lives entirely at the DECLARED 24 h read. This measures the other
    half: on every row this repository has actually committed, the model WITH the bearing
    fraction and the model without it fit equally well, so nothing measured here chooses
    between them and the 24 h verdict is an extrapolation of two equally good fits.
    """

    def test_the_two_reporter_fits_are_indistinguishable_on_every_committed_block(self,
                                                                                  plate_mu):
        from ystwin.mech.ablation import SaturatingReporter, reporter_block
        worst = 0.0
        for export, construct in pop.committed_reporter_blocks():
            data = reporter_block(export, construct)
            for p in pop.LOSS_PER_DIVISION.bounds:
                full = pop.SegregatingReporter(
                    population=pop.TwoPopulation(p, 0.0), mu_free_per_h=plate_mu)
                ratio = full.refit(data).residual_rms / SaturatingReporter().refit(
                    data).residual_rms
                worst = max(worst, abs(ratio - 1.0))
        assert worst < 0.04
        assert worst / REPORTER_ACTIVITY_FLOOR.noise_floor < 1.0

    def test_the_segregation_decay_is_absorbed_into_the_incumbents_own_decay_constant(
            self, plate_mu):
        from ystwin.mech.ablation import SaturatingReporter, reporter_block
        export, construct = pop.committed_reporter_blocks()[0]
        data = reporter_block(export, construct)
        full = pop.SegregatingReporter(
            population=pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[1], 0.0),
            mu_free_per_h=plate_mu)
        assert full.refit(data).parameters["lambda_"] < \
            SaturatingReporter().refit(data).parameters["lambda_"]

    def test_the_growth_channel_reads_one_scalar_out_of_the_hundred_and_eighty_three_wells(
            self, growth_rows):
        """Both models are constant across rows, so only the mean of the wells can matter."""
        spread = growth_rows.copy()
        mean = spread.growth_rate.mean()
        spread["growth_rate"] = mean + (spread.growth_rate - mean) * 5.0
        assert pop.score_growth_ablation(
            pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1],
            read_h=pop.DECLARED_READ_H, data=spread).effect == pytest.approx(
            pop.score_growth_ablation(
                pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1],
                read_h=pop.DECLARED_READ_H, data=growth_rows).effect, rel=1e-6)


class TestTheGrowthChannelIsConfoundedByASubstrateTermTheRepoAlreadyHas:
    """A calibration, and the limit on how the 24 h growth pass may be read.

    Deleting the segregation weight is a correct ablation of the segregation weight. It does
    not show that a 24 h OD trace could tell segregation from substrate exhaustion, and this
    measures the competing term the repository already ships.
    """

    def test_the_committed_logistic_incumbent_falls_further_than_this_block_moves(
            self, growth_rows):
        from ystwin.generator.culture import CultureParameters, simulate_culture
        conditions = PlateConditions()
        x0 = conditions.inoculum_for("any") * conditions.gdcw_per_od
        params = CultureParameters(
            mu_max=float(np.median(growth_rows.growth_rate)), growth_ic50=10.0,
            growth_hill=1.0, promoter_basal=1.0, promoter_peak=10.0, promoter_ec50=1.0,
            promoter_hill=1.0, carrying_capacity=10.0)
        chords = {}
        for duration in (pop.PLATE_READ_H, pop.DECLARED_READ_H):
            grid = np.linspace(0.0, duration, 1601)
            biomass = simulate_culture(grid, 0.0, params, x0)["biomass"]
            chords[duration] = float(np.log(biomass[-1] / biomass[0]) / duration)
        substrate_drop = chords[pop.PLATE_READ_H] - chords[pop.DECLARED_READ_H]
        block_effect = max(
            pop.score_growth_ablation(p, pop.BURDEN.bounds[1], read_h=pop.DECLARED_READ_H,
                                      data=growth_rows).effect
            for p in pop.LOSS_PER_DIVISION.bounds)
        assert substrate_drop / GROWTH_RATE_FLOOR.noise_floor > 5.0
        assert substrate_drop > 3.0 * block_effect

    def test_the_reporter_channel_does_not_share_the_confound(self, plate_mu):
        """Its competing decay is a FITTED parameter of both models, not an unmodelled term."""
        from ystwin.mech.ablation import SaturatingReporter
        assert "lambda_" in SaturatingReporter.parameter_names
        assert "lambda_" in pop.SegregatingReporter.parameter_names
        assert pop.SegregatingReporter.parameter_names == SaturatingReporter.parameter_names

    def test_the_reporter_pass_shrinks_but_survives_at_the_rate_the_incumbent_predicts(
            self, growth_rows, plate_mu):
        """The other extrapolation in the 24 h verdict: the growth rate driving segregation.

        F(t) is driven at the 4.14 h chord for the whole declared read, and the repository's
        own logistic says a 24 h read does not grow that fast. Re-scored at the rate it does
        predict, the pass survives on all twelve blocks and the margin roughly halves.
        """
        from ystwin.generator.culture import CultureParameters, simulate_culture
        conditions = PlateConditions()
        x0 = conditions.inoculum_for("any") * conditions.gdcw_per_od
        params = CultureParameters(
            mu_max=plate_mu, growth_ic50=10.0, growth_hill=1.0, promoter_basal=1.0,
            promoter_peak=10.0, promoter_ec50=1.0, promoter_hill=1.0, carrying_capacity=10.0)
        grid = np.linspace(0.0, pop.DECLARED_READ_H, 1601)
        biomass = simulate_culture(grid, 0.0, params, x0)["biomass"]
        slow = float(np.log(biomass[-1] / biomass[0]) / pop.DECLARED_READ_H)
        assert slow < plate_mu
        ratios = {}
        for mu in (plate_mu, slow):
            ratios[mu] = [pop.score_reporter_ablation(
                export, construct, pop.LOSS_PER_DIVISION.bounds[1], 0.0, mu_free_per_h=mu,
                read_h=pop.DECLARED_READ_H).ratio
                for export, construct in pop.committed_reporter_blocks()]
        assert min(ratios[slow]) > 1.0
        assert min(ratios[slow]) < 0.8 * min(ratios[plate_mu])
        assert max(ratios[slow]) < max(ratios[plate_mu])


class TestTheCollapseAcceleratesOnlyPerCapita:
    """The word "accelerates" hides a turnover, and this pins which of the two is meant."""

    def test_the_per_capita_rate_rises_monotonically_all_the_way_to_zero(self):
        system = pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1])
        per_capita = [-system.d_fraction_d_time(f, 1.0) / f
                      for f in np.linspace(1.0, 0.01, 40)]
        assert all(b > a for a, b in zip(per_capita, per_capita[1:]))

    def test_but_the_absolute_rate_turns_over_inside_the_range_of_a_fed_batch(self):
        system = pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1])
        grid = np.linspace(1.0, 0.01, 2001)
        absolute = np.array([abs(system.d_fraction_d_time(f, 1.0)) for f in grid])
        assert grid[int(np.argmax(absolute))] == pytest.approx(
            system.peak_loss_fraction, abs=2e-3)
        assert 0.0 < system.peak_loss_fraction < 1.0

    def test_the_verdict_string_says_per_capita_and_names_the_turnover(self):
        detail = pop.TwoPopulation(0.05, 0.234).takeover().summary()
        assert "PER-CAPITA" in detail
        assert "peaks at F" in detail

    def test_a_burdenless_plasmid_never_turns_over(self):
        assert pop.TwoPopulation(0.05, 0.0).peak_loss_fraction == 1.0
        assert "largest at F = 1" in pop.TwoPopulation(0.05, 0.0).takeover().summary()


class TestTheAssertedInoculumFraction:
    """The one number here with no citation. Named, tagged, and its sensitivity measured."""

    def test_it_is_the_default_of_every_entry_point(self):
        import inspect
        defaults = []
        for name in ("simulate_population", "bearing_fraction_at", "mean_bearing_fraction",
                     "score_growth_ablation", "score_reporter_ablation", "fedbatch_collapse"):
            parameter = inspect.signature(getattr(pop, name)).parameters["initial_fraction"]
            defaults.append(parameter.default)
        assert set(defaults) == {pop.INOCULUM_BEARING_FRACTION}
        assert pop.SegregatingCulture(pop.TwoPopulation(0.05, 0.1)).initial_fraction == \
            pop.INOCULUM_BEARING_FRACTION

    def test_it_is_deliberately_outside_the_gate_because_asserted_counts_as_free(self):
        assert Tag.ASSERTED in FREE_TAGS
        assert len(pop.POPULATION_GATE.free) == 2
        assert not any("inocul" in p.name for p in pop.POPULATION_PARAMS.free_scalars())

    def test_it_sits_at_the_optimistic_end_of_the_range_hohnholz_measured(self):
        measured = (0.853, 0.960)
        assert measured[0] < pop.INOCULUM_BEARING_FRACTION < measured[1]
        assert pop.INOCULUM_BEARING_FRACTION > 0.5 * (measured[0] + measured[1])

    def test_the_reporter_channel_is_exactly_blind_to_it(self, plate_mu):
        export, construct = pop.committed_reporter_blocks()[0]
        ratios = [pop.score_reporter_ablation(
            export, construct, pop.LOSS_PER_DIVISION.bounds[1], 0.0, mu_free_per_h=plate_mu,
            read_h=pop.DECLARED_READ_H, initial_fraction=f0).ratio
            for f0 in (0.80, pop.INOCULUM_BEARING_FRACTION, 1.0)]
        assert max(ratios) - min(ratios) < 1e-3

    def test_the_growth_channel_moves_with_it_and_the_verdict_does_not_flip(self, growth_rows):
        ratios = [pop.score_growth_ablation(
            pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1],
            read_h=pop.DECLARED_READ_H, initial_fraction=f0, data=growth_rows).ratio
            for f0 in (0.853, pop.INOCULUM_BEARING_FRACTION, 1.0)]
        assert all(r > 1.0 for r in ratios)
        assert ratios[0] < ratios[1] < ratios[2]
        assert ratios[2] / ratios[0] < 1.15

    def test_half_the_least_severe_fed_batch_headline_is_this_number_and_not_the_run(self):
        model = pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[0], 0.0)
        run = pop.fedbatch_collapse(model, generations=6.6294, mu_free_per_h=0.101)
        lost_in_total = 1.0 - run.biomass_weighted_fraction
        lost_during_run = 1.0 - run.biomass_weighted_fraction / pop.INOCULUM_BEARING_FRACTION
        assert lost_in_total == pytest.approx(0.099, abs=0.002)
        assert lost_during_run == pytest.approx(0.051, abs=0.002)
        assert lost_during_run < 0.6 * lost_in_total


class TestTakeoverIsUnconditionalWithoutSelection:
    """The productivity-collapse mechanism the outcome panel names."""

    def test_no_interior_equilibrium_exists_in_a_non_selective_medium(self):
        for burden in (0.0, 0.1, 0.234):
            verdict = pop.TwoPopulation(0.05, burden).takeover()
            assert verdict.unconditional
            assert verdict.interior_fraction is None
            assert verdict.selective_advantage <= 0.0

    def test_it_is_only_arrested_when_a_refused_death_rate_exceeds_b_times_mu(self):
        system = pop.TwoPopulation(0.05, 0.1)
        assert system.takeover(differential_death_per_h=0.05, mu_free_per_h=0.101
                               ).interior_fraction == pytest.approx(0.88609, rel=1e-4)
        assert system.takeover(differential_death_per_h=0.011,
                               mu_free_per_h=0.101).unconditional

    def test_the_arrest_names_k_d_as_a_refusal_so_it_cannot_be_quoted_as_a_result(self):
        detail = pop.TwoPopulation(0.05, 0.1).takeover(
            differential_death_per_h=0.05, mu_free_per_h=0.101).summary()
        assert "NOT_BUILT" in detail
        assert "k_d_plasmid_free_death" in {p.name for p in pop.NOT_BUILT.refusals()}


class TestTheTimescaleAuditDecidesTheVerdict:
    """A calibration. Re-runs `mech/integrate.py`'s tau/T rule and pins the four verdicts."""

    def test_the_plate_read_freezes_a_non_expressing_plasmid_and_the_fed_batch_keeps_it(self):
        states = pop.population_states(
            PLATE_READ_4H, p_bounds=pop.LOSS_PER_DIVISION.bounds, burden=0.0)
        assert reduce_for_window(states, PLATE_READ_4H).verdict(
            "plasmid_bearing") is Verdict.FREEZE
        states = pop.population_states(
            FEDBATCH_5D, p_bounds=pop.LOSS_PER_DIVISION.bounds, burden=0.0)
        assert reduce_for_window(states, FEDBATCH_5D).verdict(
            "plasmid_bearing") is Verdict.KEEP

    def test_an_expressing_construct_straddles_the_threshold_on_the_plate(self):
        states = pop.population_states(
            PLATE_READ_4H, p_bounds=pop.LOSS_PER_DIVISION_EXPRESSING.bounds, burden=0.0)
        audit = reduce_for_window(states, PLATE_READ_4H)
        assert audit.verdict("plasmid_bearing") is Verdict.STRADDLES
        row = audit.row("plasmid_bearing")
        assert row.ratio_low == pytest.approx(1.375, rel=1e-2)
        assert row.ratio_high == pytest.approx(33.83, rel=1e-2)

    def test_the_generation_clock_is_an_integrator_and_is_kept_in_both_windows(self):
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            states = pop.population_states(
                window, p_bounds=pop.LOSS_PER_DIVISION.bounds, burden=0.0)
            assert reduce_for_window(states, window).verdict("generations") is Verdict.KEEP

    def test_the_time_constant_is_the_per_generation_multiple_not_the_per_e_fold_one(self):
        per_generation = 1.0 / -math.log2(1 - 0.05)
        per_e_fold = 1.0 / -math.log1p(-0.05)
        assert per_generation == pytest.approx(13.5134, rel=1e-5)
        assert per_e_fold / per_generation == pytest.approx(1 / LN2, rel=1e-12)
        p_from_theta = pop.divisions_per_generation(0.05)
        states = pop.population_states(FEDBATCH_5D, p_bounds=(p_from_theta, p_from_theta),
                                       burden=0.0)
        _, tau_high = states["plasmid_bearing"].taus_in(FEDBATCH_5D)
        assert tau_high == pytest.approx(per_generation / FEDBATCH_5D.growth_rate_low_per_h,
                                         rel=1e-12)

    def test_a_band_outside_zero_one_is_refused(self):
        with pytest.raises(ValueError, match="increasing band"):
            pop.population_states(FEDBATCH_5D, p_bounds=(0.2, 0.1))


class TestTheTrajectoryAgreesWithTheGenerationDomain:
    """Two integrations of the same system, in different independent variables."""

    def test_the_stiff_driver_and_the_closed_form_agree(self):
        system = pop.TwoPopulation(0.05, 0.1)
        trajectory = pop.simulate_population(system, FEDBATCH_5D, mu_free_per_h=0.101)
        generations = trajectory.of("generations")[-1]
        assert trajectory.of("plasmid_bearing")[-1] == pytest.approx(
            system.fraction_after_generations(0.95, generations), rel=1e-5)

    def test_the_generation_count_matches_the_population_rate_it_integrated(self):
        system = pop.TwoPopulation(0.0, 0.0)
        trajectory = pop.simulate_population(system, FEDBATCH_5D, mu_free_per_h=0.101)
        assert trajectory.of("generations")[-1] == pytest.approx(
            0.101 * FEDBATCH_5D.duration_h / LN2, rel=1e-6)

    def test_an_asserted_growth_rate_is_refused(self):
        with pytest.raises(ValueError, match="mu_free must be positive"):
            pop.simulate_population(pop.TwoPopulation(0.05, 0.0), FEDBATCH_5D,
                                    mu_free_per_h=0.0)


class TestTheBearingFractionHelpers:
    """The closed form and the integrated branch have to be the same function."""

    def test_the_two_branches_agree_as_the_burden_goes_to_zero(self):
        times = np.linspace(0.0, 24.0, 9)
        closed = pop.bearing_fraction_at(pop.TwoPopulation(0.05, 0.0), times,
                                         mu_free_per_h=0.3227)
        integrated = pop.bearing_fraction_at(pop.TwoPopulation(0.05, 1e-9), times,
                                             mu_free_per_h=0.3227)
        assert integrated == pytest.approx(closed, rel=1e-6)

    def test_the_mean_fraction_is_the_initial_one_when_nothing_is_lost(self):
        assert pop.mean_bearing_fraction(pop.TwoPopulation(0.0, 0.0), mu_free_per_h=0.3227,
                                         duration_h=24.0) == pytest.approx(0.95)

    def test_the_mean_lies_between_the_endpoints_and_falls_with_the_read(self):
        system = pop.TwoPopulation(0.0847, 0.0)
        short = pop.mean_bearing_fraction(system, mu_free_per_h=0.3227, duration_h=4.14)
        long = pop.mean_bearing_fraction(system, mu_free_per_h=0.3227, duration_h=24.0)
        end = float(pop.bearing_fraction_at(system, [24.0], mu_free_per_h=0.3227)[0])
        assert end < long < short < 0.95


class TestCriterionA:
    """A calibration. Both ablations on the real committed rows; the switch is duration."""

    def test_the_growth_channel_is_exactly_blind_on_the_read_that_was_run(self, growth_rows):
        for burden in (0.0, 0.234):
            result = pop.score_growth_ablation(
                pop.LOSS_PER_DIVISION.bounds[1], burden,
                read_h=pop.PLATE_READ_H, data=growth_rows)
            assert result.effect < 1e-12
            assert not result.clears_floor

    def test_it_is_blind_to_p_at_zero_burden_however_fast_the_plasmid_is_lost(self, growth_rows):
        for p in pop.LOSS_PER_DIVISION.bounds:
            result = pop.score_growth_ablation(p, 0.0, read_h=24.0, data=growth_rows)
            assert result.effect == 0.0

    def test_at_twenty_four_hours_it_clears_the_floor_at_the_top_of_the_burden_band(
            self, growth_rows):
        result = pop.score_growth_ablation(
            pop.LOSS_PER_DIVISION.bounds[1], pop.BURDEN.bounds[1], read_h=24.0,
            data=growth_rows)
        assert result.clears_floor
        assert result.ratio == pytest.approx(2.57, rel=0.02)
        assert result.floor is GROWTH_RATE_FLOOR
        assert result.full_n_free == result.reduced_n_free == 1

    def test_the_burden_that_clears_it_is_inside_the_measured_band(self, growth_rows):
        crossing = pop.burden_to_clear_the_growth_floor(
            pop.LOSS_PER_DIVISION.bounds[1], read_h=24.0, data=growth_rows)
        assert crossing is not None
        assert pop.BURDEN.bounds[0] < crossing <= pop.BURDEN.bounds[1]

    def test_the_reporter_channel_clears_no_block_on_a_four_hour_read(self, plate_mu):
        blocks = pop.committed_reporter_blocks()
        assert len(blocks) == 12
        ratios = [pop.score_reporter_ablation(export, construct,
                                              pop.LOSS_PER_DIVISION.bounds[1], 0.0,
                                              mu_free_per_h=plate_mu,
                                              read_h=pop.PLATE_READ_H).ratio
                  for export, construct in blocks]
        assert max(ratios) < 1.0
        assert max(ratios) == pytest.approx(0.141, rel=0.05)

    def test_at_twenty_four_hours_it_clears_every_block_at_the_top_of_the_loss_band(
            self, plate_mu):
        ratios = [pop.score_reporter_ablation(export, construct,
                                              pop.LOSS_PER_DIVISION.bounds[1], 0.0,
                                              mu_free_per_h=plate_mu, read_h=24.0).ratio
                  for export, construct in pop.committed_reporter_blocks()]
        assert min(ratios) > 1.0
        assert min(ratios) == pytest.approx(1.58, rel=0.05)

    def test_and_clears_none_of_them_at_the_bottom_of_it(self, plate_mu):
        ratios = [pop.score_reporter_ablation(export, construct,
                                              pop.LOSS_PER_DIVISION.bounds[0], 0.0,
                                              mu_free_per_h=plate_mu, read_h=24.0).ratio
                  for export, construct in pop.committed_reporter_blocks()]
        assert max(ratios) < 1.0

    def test_the_two_models_are_refit_on_the_same_rows_at_the_same_parameter_count(
            self, plate_mu):
        export, construct = pop.committed_reporter_blocks()[0]
        result = pop.score_reporter_ablation(export, construct, 0.05, 0.0,
                                             mu_free_per_h=plate_mu, read_h=24.0)
        assert result.full_n_free == result.reduced_n_free == 5
        assert result.floor is REPORTER_ACTIVITY_FLOOR
        assert result.n_rows == len(pop.reporter_block(export, construct))

    def test_the_reduced_model_is_the_full_one_at_zero_loss(self, plate_mu):
        export, construct = pop.committed_reporter_blocks()[0]
        result = pop.score_reporter_ablation(export, construct, 0.0, 0.0,
                                             mu_free_per_h=plate_mu, read_h=24.0)
        assert result.effect < 1e-4

    def test_the_effect_rises_monotonically_with_the_loss_rate(self, plate_mu):
        export, construct = pop.committed_reporter_blocks()[0]
        ratios = [pop.score_reporter_ablation(export, construct, p, 0.0,
                                              mu_free_per_h=plate_mu, read_h=24.0).ratio
                  for p in (0.0, 0.01, 0.03, 0.0847)]
        assert ratios == sorted(ratios)
        assert ratios[0] < 1.0 < ratios[-1]

    def test_fold_induction_is_exactly_blind_to_this_block_which_is_why_it_is_not_scored(self):
        parameters = {"R0": 100.0, "k_basal": 50.0, "k_max": 500.0, "K_dose": 1.0,
                      "lambda_": 0.4}
        query = pd.DataFrame({"time_h": [24.0, 24.0], "dose_mM": [4.0, 0.0],
                              "signal": [0.0, 0.0]})
        folds = []
        for p in (0.0, 0.0847, 0.3):
            model = pop.SegregatingReporter(population=pop.TwoPopulation(p, 0.0),
                                            mu_free_per_h=0.3227)
            dosed, control = model.predict(query, parameters)
            folds.append(dosed / control)
        assert folds[0] == pytest.approx(folds[1], rel=1e-12)
        assert folds[0] == pytest.approx(folds[2], rel=1e-12)

    def test_a_reporter_model_without_a_measured_growth_rate_is_refused(self):
        with pytest.raises(ValueError, match="never fitted"):
            pop.SegregatingReporter(population=pop.TwoPopulation(0.05, 0.0))
        with pytest.raises(ValueError, match="does not invent them"):
            pop.SegregatingCulture()

    def test_the_growth_wells_are_the_ones_the_floor_was_measured_on(self, growth_rows):
        assert len(growth_rows) == 183
        assert set(growth_rows.plate.unique()) == {
            "g1_20260722_ERandOxidativeStress_NewProtocol_ANALYSED.csv",
            "g1_20260728_ERandOxidativeStress_NewProtocol_Replicate2.csv",
            "g1_20260803_ERandoxidativestress_Replicate3.csv"}
        assert float(np.median(growth_rows.growth_rate)) == pytest.approx(
            0.32482577175526683, rel=1e-3)

    def test_the_read_length_is_the_committed_one(self):
        assert pop.PLATE_READ_H == PlateConditions().duration_h == PLATE_READ_4H.duration_h


class TestTheFedBatchCollapse:
    """A calibration. What segregation costs the vessel `design_fedbatch_run.py` publishes."""

    @staticmethod
    def _design(inoculum):
        return design_fedbatch(0.101, 0.254, initial_biomass_g_per_l=inoculum,
                               assumed_yield_g_per_g=0.098, max_uptake_g_per_gdcw_h=2.59,
                               initial_volume_l=1.0, max_volume_l=2.0,
                               feed_substrate_g_per_l=500.0)

    def test_the_titre_a_run_keeps_falls_with_the_generation_budget(self):
        system = pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[1], 0.0)
        rich = pop.fedbatch_collapse(
            system, generations=self._design(0.5).generation_budget_at(0.098))
        lean = pop.fedbatch_collapse(
            system, generations=self._design(0.05).generation_budget_at(0.098))
        assert lean.generations > rich.generations
        assert lean.titre_retained < rich.titre_retained
        assert rich.titre_retained == pytest.approx(0.700, rel=0.01)
        assert lean.titre_retained == pytest.approx(0.579, rel=0.01)

    def test_biomass_weighting_is_harsher_than_cell_hours_because_growth_is_exponential(self):
        collapse = pop.fedbatch_collapse(pop.TwoPopulation(0.0847, 0.234), generations=6.63)
        assert (collapse.final_fraction < collapse.biomass_weighted_fraction
                < collapse.cell_hour_fraction < 0.95)

    def test_wall_clock_moves_with_the_feed_and_generations_do_not(self):
        system = pop.TwoPopulation(0.05, 0.1)
        slow = pop.fedbatch_collapse(system, generations=6.0, mu_free_per_h=0.05)
        fast = pop.fedbatch_collapse(system, generations=6.0, mu_free_per_h=0.3)
        assert slow.generations_to_half == pytest.approx(fast.generations_to_half)
        assert slow.final_fraction == pytest.approx(fast.final_fraction)
        assert slow.hours_to_half > 5 * fast.hours_to_half

    def test_nothing_is_lost_when_nothing_is_lost(self):
        collapse = pop.fedbatch_collapse(pop.TwoPopulation(0.0, 0.0), generations=10.0)
        assert collapse.final_fraction == pytest.approx(0.95)
        assert collapse.titre_retained == pytest.approx(0.95, rel=1e-6)

    def test_a_run_of_no_generations_is_refused(self):
        with pytest.raises(ValueError, match="generations must be positive"):
            pop.fedbatch_collapse(pop.TwoPopulation(0.05, 0.0), generations=0.0)

    def test_hours_to_half_does_not_depend_on_whether_the_run_reached_it(self):
        """The same integral either way. The closed form it replaced was 7.6% high at b = 0.234."""
        system = pop.TwoPopulation(pop.LOSS_PER_DIVISION.bounds[0], pop.BURDEN.bounds[1])
        short = pop.fedbatch_collapse(system, generations=2.0, mu_free_per_h=0.101,
                                      n_points=201)
        reached = pop.fedbatch_collapse(system, generations=40.0, mu_free_per_h=0.101,
                                        n_points=201)
        assert short.generations_to_half > short.generations
        assert reached.generations_to_half < reached.generations
        assert short.hours_to_half == pytest.approx(reached.hours_to_half, rel=2e-3)


class TestTheEthanolDiagnostic:
    """REPORTED, and it feeds nothing. The band is recomputed from the committed table."""

    def test_the_swept_band_is_recomputed_from_the_committed_states_not_transcribed(self):
        source = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
        with open(source, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        yields = [float(r["q_ethanol"]) * float(pop.ETHANOL_G_PER_MOL) / 1000.0
                  / float(r["mu_per_h"]) for r in rows]
        assert len(yields) == 6
        assert pop.ethanol_yield_band() == (min(yields), max(yields))
        assert pop.ETHANOL_YIELD.bounds == pytest.approx((min(yields), max(yields)), rel=1e-12)

    def test_the_standing_concentration_r_was_measured_at_is_recomputed_too(self):
        """The scale caveat, as a number: r comes from states holding ~1 g/L, not tens."""
        with open(paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv",
                  newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected = [float(r["q_ethanol"]) * float(pop.ETHANOL_G_PER_MOL) / 1000.0
                    / float(r["mu_per_h"]) * float(r["biomass_g_per_l"]) for r in rows]
        assert pop.standing_ethanol_band() == pytest.approx(
            (min(expected), max(expected)), rel=1e-12)
        assert pop.standing_ethanol_band() == pytest.approx((0.2787, 1.2632), rel=1e-3)

    def test_the_diagnostic_extrapolates_r_by_at_least_twentyfive_fold(self):
        design = design_fedbatch(
            0.101, 0.254, initial_biomass_g_per_l=0.5, assumed_yield_g_per_g=0.098,
            max_uptake_g_per_gdcw_h=2.59, initial_volume_l=1.0, max_volume_l=2.0,
            feed_substrate_g_per_l=500.0)
        harvest = pop.EthanolDiagnostic(*pop.ethanol_yield_band()).at_harvest_g_per_l(design)
        assert harvest[0] / pop.standing_ethanol_band()[1] > 25.0
        assert "q_ethanol_uptake" in {p.name for p in pop.NOT_BUILT.refusals()}

    def test_the_molar_masses_are_the_ones_the_rest_of_the_package_uses(self):
        assert float(pop.ETHANOL_G_PER_MOL) == pytest.approx(46.069, abs=5e-4)
        assert float(pop.GLUCOSE_G_PER_MOL) == pytest.approx(180.156, abs=5e-4)

    def test_the_plate_ceiling_is_parameter_free_gay_lussac_stoichiometry(self):
        diagnostic = pop.EthanolDiagnostic()
        expected = (20.0 / float(pop.GLUCOSE_G_PER_MOL)) * 2.0 * float(pop.ETHANOL_G_PER_MOL)
        assert diagnostic.plate_ceiling_g_per_l(20.0) == pytest.approx(expected)
        assert diagnostic.plate_ceiling_g_per_l() == pytest.approx(10.2287, rel=1e-4)
        assert pop.PLATE_GLUCOSE_G_PER_L == CultureContext().glucose_g_per_L

    def test_the_fed_batch_harvest_sits_several_fold_above_that_ceiling(self):
        design = TestTheFedBatchCollapse._design(0.5)
        low, high = pop.EthanolDiagnostic().at_harvest_g_per_l(design)
        assert (low, high) == pytest.approx((32.1, 70.8), rel=1e-2)
        assert low / pop.EthanolDiagnostic().plate_ceiling_g_per_l() > 3.0

    def test_it_scales_with_the_biomass_made_and_refuses_an_empty_vessel(self):
        diagnostic = pop.EthanolDiagnostic()
        assert diagnostic.accumulated_g_per_l(0.0, 2.0) == (0.0, 0.0)
        one, two = (diagnostic.accumulated_g_per_l(g, 2.0)[0] for g in (10.0, 20.0))
        assert two == pytest.approx(2 * one)
        with pytest.raises(ValueError):
            diagnostic.accumulated_g_per_l(1.0, 0.0)
        with pytest.raises(ValueError):
            diagnostic.accumulated_g_per_l(-1.0, 2.0)

    def test_asking_it_what_the_ethanol_does_raises_and_names_the_five_constants(self):
        with pytest.raises(RefusedValue) as excinfo:
            pop.EthanolDiagnostic().growth_penalty_per_h(50.0)
        message = str(excinfo.value)
        assert "FIVE" in message
        for name in ("i_ethanol_growth_inhibition", "E_theta_inhibition_threshold",
                     "q_ethanol_uptake", "k_strip_ethanol_evaporation"):
            assert name in {p.name for p in pop.NOT_BUILT.refusals()}

    def test_it_cannot_be_promoted_out_of_reported_by_setting_a_field(self):
        with pytest.raises(ValueError, match="cannot be promoted"):
            pop.EthanolDiagnostic(1.0, 2.0, LayerState.IN_CHAIN)

    def test_an_unordered_band_is_refused(self):
        with pytest.raises(ValueError, match="positive and ordered"):
            pop.EthanolDiagnostic(2.0, 1.0)


class TestProvenance:
    """Every constant graded, and the summary that says both gates out loud."""

    def test_every_registered_parameter_carries_a_tag_and_units(self):
        for registry in (pop.POPULATION_PARAMS, pop.ETHANOL_PARAMS, pop.BAND_VARIANTS,
                         pop.NOT_BUILT):
            assert len(registry) > 0
            for param in registry.params:
                assert isinstance(param, Param)
                assert param.tag in Tag.ALL
                assert param.units.strip() and param.source.strip()

    def test_no_pinned_number_in_this_module_is_uncited(self):
        assert pop.POPULATION_PARAMS.uncited() == ()
        assert pop.BAND_VARIANTS.uncited() == ()

    def test_nothing_here_is_borrowed_from_another_organism(self):
        for registry in (pop.POPULATION_PARAMS, pop.ETHANOL_PARAMS, pop.BAND_VARIANTS):
            assert registry.borrowed() == ()

    def test_the_summary_reports_both_gates_and_every_refusal(self):
        summary = pop.provenance_summary()
        assert "PASSES" in summary and "REFUSED" in summary
        for param in pop.NOT_BUILT.refusals():
            assert param.name in summary

    def test_the_criterion_a_table_is_computed_and_names_both_floors(self, plate_mu):
        table = pop.criterion_a_table(mu_free_per_h=plate_mu, read_hours=(24.0,),
                                      burdens=(0.234,))
        assert "12" in table and "183 wells" in table
        assert f"{REPORTER_ACTIVITY_FLOOR.noise_floor:g}" in table
        assert f"{GROWTH_RATE_FLOOR.noise_floor:g}" in table
