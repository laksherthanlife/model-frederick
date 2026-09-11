"""Batch and fed-batch integrated over time, where titre is an integral rather than a balance.

Everything else in this package solves a steady state, which is right for a chemostat and
wrong for the vessel anyone manufactures in. These tests pin the integrator's arithmetic, the
conservation it must not break, and the two failure modes that actually bit while it was
being written:

  * a product rate that is a STEP FUNCTION of substrate permits unlimited production as
    substrate asymptotes toward zero. The first version also had a numerical error, returning
    3.14 mg/L at 1,600 steps and 64.2 at 6,400; convergence alone does not validate the rates.
  * clamping the substrate STATE instead of the RATE lets a culture integrate a negative
    concentration and report a titre computed from it.

Both are pinned below, because a convergence bug that only shows at one step count is exactly
the kind of thing that survives review.
"""

from __future__ import annotations

import math

import cobra
import numpy as np
import pytest

from ystwin.fba.dynamic import BatchResult, FeedProfile, simulate_batch

MU_MAX, KS, Y_XS, Q_P = 0.35, 0.5, 0.09, 2.0e-4
MOLAR_MASS = 536.87


def _rates(substrate):
    """A continuous closure: product formation saturates with the same carbon as growth."""
    f = substrate / (KS + substrate) if substrate > 0 else 0.0
    return MU_MAX * f, MU_MAX * f / Y_XS, Q_P * f


def _run(**kw):
    kw.setdefault("initial_biomass_g_per_l", 0.05)
    kw.setdefault("initial_substrate_mmol_per_l", 111.0)
    kw.setdefault("hours", 72.0)
    kw.setdefault("product_molar_mass_g_per_mol", MOLAR_MASS)
    return simulate_batch(_rates, **kw)


class TestTheIntegratorConverges:
    def test_the_titre_is_the_same_across_a_256_fold_range_of_step_counts(self):
        """The check that caught the step-function bug. A trajectory that moves with the
        step count is a result about the integrator, not about the culture."""
        titres = [_run(steps=n).titre_mg_per_l for n in (100, 400, 1600, 6400, 25600)]

        assert max(titres) - min(titres) < 1e-6 * max(titres)

    def test_a_converged_step_function_product_rate_is_still_physically_invalid(self):
        """The failure mode, reproduced deliberately so the guard above cannot rot.

        Monod substrate remains positive at every finite time. A constant `q_P` therefore
        integrates accurately but keeps making product without measurable carbon uptake.
        """
        def stepwise(s):
            f = s / (KS + s) if s > 0 else 0.0
            return MU_MAX * f, MU_MAX * f / Y_XS, (Q_P if s > 0 else 0.0)
        coarse = simulate_batch(stepwise, 0.05, 111.0, 72.0, MOLAR_MASS, steps=400)
        fine = simulate_batch(stepwise, 0.05, 111.0, 72.0, MOLAR_MASS, steps=6400)

        assert fine.titre_mg_per_l == pytest.approx(coarse.titre_mg_per_l, rel=1e-6)
        midpoint = len(fine.hours) // 2
        assert fine.substrate_mmol_per_l[midpoint] < 1e-9
        extra_product = fine.product_mmol_per_l[-1] - fine.product_mmol_per_l[midpoint]
        assert extra_product == pytest.approx(Q_P * fine.final_biomass_g_per_l * 36.0, rel=1e-6)
        valid = _run(steps=400)
        assert valid.substrate_exhausted
        assert fine.titre_mg_per_l > 10 * valid.titre_mg_per_l


class TestConservation:
    def test_the_biomass_yield_is_recovered_from_the_substrate_consumed(self):
        run = _run(steps=2000)
        consumed = 111.0 - run.substrate_mmol_per_l[-1]

        assert run.final_biomass_g_per_l / consumed == pytest.approx(Y_XS, rel=0.02)

    def test_substrate_never_goes_negative(self):
        run = _run(steps=400)

        assert (run.substrate_mmol_per_l >= 0.0).all()

    def test_growth_stops_when_the_carbon_does(self):
        run = _run(steps=2000)
        exhausted = (run.substrate_mmol_per_l <= 1e-9).argmax()

        assert run.substrate_exhausted
        assert run.growth_rate_per_h[exhausted + 1:].max() < 1e-9

    def test_biomass_is_monotonic(self):
        run = _run(steps=400)

        assert all(b >= a - 1e-12 for a, b in zip(run.biomass_g_per_l,
                                                  run.biomass_g_per_l[1:]))


class TestFedBatchDiffersFromBatchInTheRightWay:
    def test_feeding_raises_biomass_and_titre(self):
        batch = _run(steps=2000)
        fed = _run(steps=2000, feed=FeedProfile(rate_mmol_per_l_h=3.0, start_h=24.0))

        assert fed.final_biomass_g_per_l > batch.final_biomass_g_per_l
        assert fed.titre_mg_per_l > batch.titre_mg_per_l

    def test_a_zero_feed_is_exactly_a_batch(self):
        """One number separates the two modes, and it separates them exactly."""
        batch = _run(steps=400)
        zero_fed = _run(steps=400, feed=FeedProfile(rate_mmol_per_l_h=0.0))

        assert zero_fed.titre_mg_per_l == batch.titre_mg_per_l

    def test_the_feed_window_is_respected(self):
        profile = FeedProfile(rate_mmol_per_l_h=3.0, start_h=24.0, end_h=48.0)

        assert profile.at(12.0) == 0.0
        assert profile.at(30.0) == 3.0
        assert profile.at(60.0) == 0.0


class TestItRefusesBadInputs:
    @pytest.mark.parametrize("kw", [{"hours": 0.0}, {"hours": -1.0}])
    def test_a_non_positive_duration_is_refused(self, kw):
        with pytest.raises(ValueError, match="hours"):
            _run(**kw)

    def test_a_zero_inoculum_is_refused_rather_than_silently_never_growing(self):
        """`dX/dt = mu*X` is zero for all time at X = 0, so this would be a no-growth run
        reported as an ordinary result."""
        with pytest.raises(ValueError, match="stays at nothing"):
            _run(initial_biomass_g_per_l=0.0)

    def test_one_step_is_refused(self):
        with pytest.raises(ValueError, match="at least 2"):
            _run(steps=1)


class TestTheReportedQuantities:
    def test_titre_is_the_product_concentration_times_molar_mass(self):
        run = _run(steps=400)

        assert run.titre_mg_per_l == pytest.approx(
            run.product_mmol_per_l[-1] * MOLAR_MASS)

    def test_content_is_titre_over_final_biomass(self):
        run = _run(steps=400)

        assert run.content_mg_per_gdcw == pytest.approx(
            run.titre_mg_per_l / run.final_biomass_g_per_l)

    def test_the_summary_names_all_three(self):
        text = _run(steps=400).summary()

        assert "mg/L" in text and "gDCW/L" in text and "mg/gDCW" in text


class TestTheRateClosureBuilder:
    """`fba/dynamic_rates.py` -- where the GEM, stress and regulation enter."""

    def test_it_refuses_without_a_kcat_or_an_abundance(self):
        from ystwin.fba.dynamic_rates import build_rates

        with pytest.raises(ValueError, match="no default"):
            build_rates(None, "g", "o", "n", kcat_per_s=0.0, enzyme_mmol_per_gdcw=1e-5)

    def test_stress_outside_zero_to_one_is_refused(self):
        from ystwin.fba.dynamic_rates import build_rates

        with pytest.raises(ValueError, match="0-1"):
            build_rates(None, "g", "o", "n", kcat_per_s=0.1,
                        enzyme_mmol_per_gdcw=1e-5, stress=1.5)

    def test_the_product_ceiling_is_kcat_times_enzyme_times_regulation(self):
        from ystwin.fba.dynamic_rates import RateClosure

        info = RateClosure(kcat_per_s=0.1, enzyme_mmol_per_gdcw=2e-5,
                           regulation_scale=3.0, stress=0.0, ngam=0.7,
                           uptake_vmax_mmol_per_gdcw_h=10.0, uptake_km_mmol_per_l=0.5)

        assert info.q_product_max_mmol_per_gdcw_h == pytest.approx(0.1 * 3600 * 2e-5 * 3.0)

    def test_regulation_scales_the_product_ceiling_linearly(self):
        """The route by which regulation acts: on how much enzyme there is, not on a
        reaction ceiling. E-Flux's ceilings never bind in this model."""
        from ystwin.fba.dynamic_rates import RateClosure

        def q(scale):
            return RateClosure(0.1, 2e-5, scale, 0.0, 0.7, 10.0, 0.5
                               ).q_product_max_mmol_per_gdcw_h

        assert q(2.0) == pytest.approx(2 * q(1.0))
        assert q(0.5) == pytest.approx(0.5 * q(1.0))

    def test_stress_maps_to_the_latent_bridge_maintenance_bound(self):
        from ystwin.fba.dynamic_rates import MAX_STRESS_NGAM, RESTING_NGAM

        assert RESTING_NGAM == 0.7
        assert RESTING_NGAM + MAX_STRESS_NGAM == pytest.approx(7.2)


def _coupled_network():
    model = cobra.Model("carbon_accounting")
    carbon_e = cobra.Metabolite("carbon_e", formula="C", compartment="e")
    carbon_c = cobra.Metabolite("carbon_c", formula="C", compartment="c")
    product = cobra.Metabolite("product_c", formula="C2", compartment="c")
    oxygen = cobra.Metabolite("oxygen_e", formula="O2", compartment="e")
    energy = cobra.Metabolite("energy_c", compartment="c")
    definitions = (
        ("EX_carbon", {carbon_e: -1}, (-10, 1000)),
        ("transport", {carbon_e: -1, carbon_c: 1}, (0, 1000)),
        ("growth", {carbon_c: -1}, (0, 1000)),
        ("synthesis", {carbon_c: -2, product: 1}, (0, 1000)),
        ("DM_product", {product: -1}, (0, 1000)),
        ("EX_oxygen", {oxygen: -1}, (-20, 1000)),
        ("energy", {carbon_c: -1, energy: 1}, (0, 1000)),
        ("maintenance", {energy: -1}, (0, 1000)),
    )
    for name, metabolites, bounds in definitions:
        reaction = cobra.Reaction(name, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites(metabolites)
        model.add_reactions([reaction])
    model.objective = "growth"
    return model


def _coupled_run(model=None, **kwargs):
    from ystwin.fba.dynamic_rates import build_rates

    options = {
        "product_reaction_id": "DM_product",
        "biomass_reaction_id": "growth",
        "growth_fraction": 0.8,
        "grid_points": 9,
    }
    options.update(kwargs)
    return build_rates(
        _coupled_network() if model is None else model,
        "EX_carbon", "EX_oxygen", "maintenance",
        kcat_per_s=1.0, enzyme_mmol_per_gdcw=1.0, **options)


class TestProductAndGrowthShareTheSameCarbon:
    def test_a_product_reaction_is_required_before_a_rate_can_be_reported(self):
        from ystwin.fba.dynamic_rates import build_rates

        with pytest.raises(ValueError, match="product_reaction_id"):
            build_rates(_coupled_network(), "EX_carbon", "EX_oxygen", "maintenance",
                        kcat_per_s=1.0, enzyme_mmol_per_gdcw=1.0)

    def test_product_cannot_be_created_outside_the_metabolic_solve(self):
        rates, _ = _coupled_run()
        mu, uptake, product = rates(100.0)

        assert product > 0
        assert mu + 2 * product + 0.7 <= uptake + 1e-7
        assert uptake <= 10.0

    def test_blocking_the_product_pathway_stops_product(self):
        model = _coupled_network()
        model.reactions.synthesis.bounds = (0.0, 0.0)
        rates, _ = _coupled_run(model)

        assert rates(100.0)[2] == pytest.approx(0.0, abs=1e-9)

    def test_scenarios_do_not_change_the_callers_model(self):
        model = _coupled_network()
        before = {r.id: r.bounds for r in model.reactions}
        objective = str(model.objective.expression)
        _coupled_run(model, stress=0.5)

        assert {r.id: r.bounds for r in model.reactions} == before
        assert str(model.objective.expression) == objective

    def test_a_learned_growth_retention_reaches_the_metabolic_solve(self):
        normal, _ = _coupled_run(growth_retention=1.0)
        stressed, _ = _coupled_run(growth_retention=0.5)

        assert stressed(100.0)[0] < normal(100.0)[0]
        assert stressed(100.0) != normal(100.0)

    def test_starvation_does_not_return_any_formation_or_uptake(self):
        rates, _ = _coupled_run()

        assert rates(0.0) == (0.0, 0.0, 0.0)

    def test_actual_uptake_not_the_allowed_uptake_is_reported(self):
        rates, _ = _coupled_run(regulation_scale=1e-6, growth_fraction=0.5)
        mu, uptake, product = rates(100.0)

        assert uptake == pytest.approx(mu + 2 * product + 0.7, abs=1e-7)
        assert uptake < 9.0

    def test_the_integrated_run_keeps_its_carbon_ledger(self):
        rates, _ = _coupled_run()
        run = simulate_batch(rates, 0.05, 20.0, 12.0, 24.0, steps=400)
        growth_carbon = run.final_biomass_g_per_l - 0.05
        product_carbon = 2 * run.product_mmol_per_l[-1]
        consumed = 20.0 - run.substrate_mmol_per_l[-1]

        assert product_carbon > 0
        assert growth_carbon + product_carbon <= consumed + 1e-6


class TestRateTitreAndYieldHaveExplicitUnits:
    def test_average_productivity_uses_the_elapsed_time(self):
        run = _run()

        assert run.average_productivity_mg_per_l_h == pytest.approx(run.titre_mg_per_l / 72.0)

    def test_yield_uses_consumed_substrate_and_its_own_molar_mass(self):
        run = _run(substrate_molar_mass_g_per_mol=180.156)
        consumed = 111.0 - run.substrate_mmol_per_l[-1]

        assert run.yield_g_per_g == pytest.approx(run.titre_mg_per_l / (consumed * 180.156))

    def test_yield_includes_substrate_that_arrived_in_feed(self):
        run = _run(substrate_molar_mass_g_per_mol=180.156,
                   feed=FeedProfile(rate_mmol_per_l_h=3.0, start_h=24.0, end_h=48.0),
                   steps=720)
        consumed = 111.0 + run.total_feed_mmol_per_l - run.substrate_mmol_per_l[-1]

        assert run.total_feed_mmol_per_l == pytest.approx(72.0, abs=0.31)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(consumed)
        assert run.yield_g_per_g == pytest.approx(run.titre_mg_per_l / (consumed * 180.156))

    def test_a_substrate_molar_mass_is_not_guessed(self):
        run: BatchResult = _run()

        with pytest.raises(ValueError, match="substrate_molar_mass"):
            _ = run.yield_g_per_g

    @pytest.mark.parametrize("name,value", [
        ("initial_substrate_mmol_per_l", -1.0),
        ("initial_substrate_mmol_per_l", math.nan),
        ("hours", math.inf),
        ("product_molar_mass_g_per_mol", -1.0),
        ("substrate_molar_mass_g_per_mol", 0.0),
    ])
    def test_invalid_physical_inputs_are_refused(self, name, value):
        with pytest.raises(ValueError, match=name):
            _run(**{name: value})

    def test_a_negative_product_rate_is_not_integrated(self):
        with pytest.raises(ValueError, match="rates"):
            simulate_batch(lambda s: (0.1, 1.0, -0.1), 0.05, 111.0, 24.0, MOLAR_MASS)


class TestAdaptiveIntegrationAndEvents:
    @pytest.mark.parametrize("steps", [2, 7, 400])
    def test_maintenance_continues_at_its_full_rate_until_depletion(self, steps):
        maintenance, uptake_slope = 0.7, 0.4
        depletion_time = math.log1p(uptake_slope / maintenance) / uptake_slope

        def rates(substrate):
            assert substrate >= 0.0
            if substrate == 0.0:
                return 0.0, 0.0, 0.0
            return 0.0, maintenance + uptake_slope * substrate, uptake_slope * substrate / 2

        run = simulate_batch(rates, 1.0, 1.0, 3.0, 24.0, steps=steps)
        expected_product = (1.0 - maintenance * depletion_time) / 2.0

        assert run.substrate_exhausted
        assert run.product_mmol_per_l[-1] == pytest.approx(expected_product, abs=1e-8)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(1.0, abs=1e-10)
        assert (run.substrate_mmol_per_l >= 0.0).all()
        assert run.product_mmol_per_l[-1] == run.product_mmol_per_l[-2]

    def test_a_depletion_event_does_not_allocate_the_rest_of_a_step_to_production(self):
        def rates(substrate):
            return (0.0, 2.0, 0.5) if substrate > 0.0 else (0.0, 0.0, 0.0)

        run = simulate_batch(rates, 1.0, 1.0, 2.0, 24.0, steps=2)

        np.testing.assert_allclose(run.substrate_mmol_per_l, [1.0, 0.0, 0.0], atol=1e-10)
        np.testing.assert_allclose(run.product_mmol_per_l, [0.0, 0.25, 0.25], atol=1e-9)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(1.0, abs=1e-9)

    @pytest.mark.parametrize("steps", [2, 37, 400])
    def test_off_grid_feed_boundaries_match_the_piecewise_analytic_solution(self, steps):
        rate, uptake, start, end, hours = 3.0, 0.8, 0.37, 1.13, 2.0
        feed = FeedProfile(rate, start, end)

        def rates(substrate):
            return 0.0, uptake * substrate, uptake * substrate / 2.0

        run = simulate_batch(rates, 1.0, 0.0, hours, 24.0, feed=feed, steps=steps)
        expected_substrate = []
        for hour in run.hours:
            if hour <= start:
                value = 0.0
            elif hour <= end:
                value = rate / uptake * (1.0 - math.exp(-uptake * (hour - start)))
            else:
                value = rate / uptake * (1.0 - math.exp(-uptake * (end - start)))
                value *= math.exp(-uptake * (hour - end))
            expected_substrate.append(value)
        supply = np.array([feed.integral(0.0, hour) for hour in run.hours])

        np.testing.assert_allclose(run.substrate_mmol_per_l, expected_substrate, atol=1e-8)
        np.testing.assert_allclose(run.substrate_mmol_per_l + 2 * run.product_mmol_per_l,
                                   supply, atol=1e-12, rtol=0.0)
        assert run.total_feed_mmol_per_l == rate * (end - start)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(
            2 * run.product_mmol_per_l[-1], abs=1e-12)

    def test_feed_is_not_spent_before_it_arrives(self):
        def rates(substrate):
            return (0.5, 1.0, 0.1) if substrate > 0.0 else (0.0, 0.0, 0.0)

        feed = FeedProfile(4.0, 0.75, 1.0)
        run = simulate_batch(rates, 1.0, 0.0, 1.0, 24.0, feed=feed, steps=2)
        exposure = math.expm1(0.5 * 0.25) / 0.5

        assert run.product_mmol_per_l[1] == 0.0
        assert run.biomass_g_per_l[1] == 1.0
        assert run.final_biomass_g_per_l == pytest.approx(math.exp(0.5 * 0.25), rel=1e-7)
        assert run.product_mmol_per_l[-1] == pytest.approx(0.1 * exposure, rel=1e-7)
        assert run.substrate_mmol_per_l[-1] == pytest.approx(1.0 - exposure, abs=1e-8)

    def test_a_feed_window_shorter_than_an_output_interval_is_not_lost(self):
        feed = FeedProfile(3.0, 0.31, 0.310001)
        run = simulate_batch(lambda s: (0.0, 0.0, 0.0), 1.0, 0.0, 2.0, 24.0,
                             feed=feed, steps=2)

        assert run.total_feed_mmol_per_l == feed.integral(0.0, 2.0)
        assert run.substrate_mmol_per_l[-1] == pytest.approx(run.total_feed_mmol_per_l, abs=1e-15)
        assert run.substrate_consumed_mmol_per_l == 0.0
        assert run.product_mmol_per_l[-1] == 0.0

    def test_a_starved_batch_stays_at_rest_without_negative_or_created_carbon(self):
        run = _run(initial_substrate_mmol_per_l=0.0, steps=2)

        np.testing.assert_array_equal(run.biomass_g_per_l, [0.05, 0.05, 0.05])
        assert not run.substrate_mmol_per_l.any()
        assert not run.product_mmol_per_l.any()
        assert not run.growth_rate_per_h.any()
        assert run.substrate_consumed_mmol_per_l == 0.0
        assert run.mean_specific_rate_mmol_per_gdcw_h == 0.0

    def test_a_dormant_positive_substrate_state_keeps_its_unconsumed_carbon(self):
        def rates(substrate):
            return (0.0, 1.0, 0.0) if substrate > 0.25 else (0.0, 0.0, 0.0)

        run = simulate_batch(rates, 1.0, 1.0, 3.0, 24.0, steps=2)

        assert run.substrate_mmol_per_l[-1] == pytest.approx(0.25, abs=1e-7)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(0.75, abs=1e-7)
        assert not run.substrate_exhausted

    @pytest.mark.parametrize("steps", [2, 400])
    def test_mean_specific_rate_uses_integrated_biomass_not_left_endpoint_exposure(self, steps):
        mu, uptake, product = 0.4, 1.0, 0.2
        run = simulate_batch(lambda s: (mu, uptake, product), 1.0, 100.0, 2.0, 24.0,
                             steps=steps)
        exposure = math.expm1(mu * 2.0) / mu

        assert run.final_biomass_g_per_l == pytest.approx(math.exp(mu * 2.0), rel=1e-7)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(uptake * exposure, rel=1e-7)
        assert run.product_mmol_per_l[-1] == pytest.approx(product * exposure, rel=1e-7)
        assert run.mean_specific_rate_mmol_per_gdcw_h == pytest.approx(product, abs=1e-12)
        assert run.growth_rate_per_h[-1] == mu

    def test_actual_rates_are_sampled_at_the_output_state_including_the_endpoint(self):
        run = _run(hours=6.0, initial_substrate_mmol_per_l=1.0, steps=2)
        expected = [_rates(substrate)[0] for substrate in run.substrate_mmol_per_l]

        np.testing.assert_array_equal(run.growth_rate_per_h, expected)
        assert run.growth_rate_per_h[-1] != run.growth_rate_per_h[-2]

    def test_tighter_tolerances_improve_a_coarse_nonlinear_trajectory(self):
        def rates(substrate):
            return 0.0, substrate ** 2, substrate ** 2 / 2.0

        loose = simulate_batch(rates, 1.0, 1.0, 2.0, 24.0, steps=2, rtol=1e-3, atol=1e-6)
        tight = simulate_batch(rates, 1.0, 1.0, 2.0, 24.0, steps=2, rtol=1e-9, atol=1e-12)
        exact = 1.0 / 3.0

        assert abs(tight.substrate_mmol_per_l[-1] - exact) < 1e-4 * abs(
            loose.substrate_mmol_per_l[-1] - exact)
        assert tight.product_mmol_per_l[-1] == pytest.approx((1.0 - exact) / 2, abs=1e-9)

    @pytest.mark.parametrize("rates", [(0.0, 0.0, 0.1), (0.1, 0.0, 0.0), (0.1, 1.0, 0.1)])
    def test_rates_cannot_create_biomass_or_product_from_depleted_substrate(self, rates):
        with pytest.raises(ValueError, match="substrate"):
            simulate_batch(lambda s: rates, 1.0, 0.0, 2.0, 24.0, steps=2)

    @pytest.mark.parametrize("name,value", [
        ("rtol", math.nan), ("rtol", math.inf), ("rtol", 0.0), ("rtol", 1e-20),
        ("atol", math.nan), ("atol", math.inf), ("atol", -1.0),
        ("steps", True), ("steps", 2.5), ("steps", math.inf),
        ("product_molar_mass_g_per_mol", None), ("feed", 0.0),
    ])
    def test_invalid_numerical_inputs_fail_before_requesting_rates(self, name, value):
        def unexpected_rates(substrate):
            pytest.fail("invalid input reached the rate callable")

        options = dict(initial_biomass_g_per_l=1.0, initial_substrate_mmol_per_l=1.0,
                       hours=1.0, product_molar_mass_g_per_mol=24.0, steps=2)
        options[name] = value
        with pytest.raises(ValueError, match=name):
            simulate_batch(unexpected_rates, **options)

    @pytest.mark.parametrize("rates", [None, (0.0, 1.0), (0.0, 1.0, 0.0, 1.0), [[0.0, 1.0, 0.0]]])
    def test_malformed_rates_are_refused(self, rates):
        with pytest.raises(ValueError, match="rates"):
            simulate_batch(lambda s: rates, 1.0, 1.0, 1.0, 24.0, steps=2)

    def test_a_finite_feed_rate_cannot_overflow_the_substrate_ledger(self):
        def unexpected_rates(substrate):
            pytest.fail("overflowed feed reached the rate callable")

        with pytest.raises(ValueError, match="feed"):
            simulate_batch(unexpected_rates, 1.0, 1.0, 2.0, 24.0,
                           feed=FeedProfile(1e308), steps=2)

    def test_finite_rates_cannot_overflow_the_integrated_state(self):
        with pytest.raises(ValueError, match="finite"):
            simulate_batch(lambda s: (1e308, 1.0, 0.1), 1e308, 1.0, 1.0, 24.0, steps=2)

    def test_integrated_uptake_is_not_lost_to_subtraction_of_large_inventories(self):
        run = simulate_batch(lambda s: (0.0, 1.0, 0.5), 1.0, 1e20, 1.0, 24.0, steps=2)

        assert run.substrate_mmol_per_l[0] == run.substrate_mmol_per_l[-1]
        assert run.substrate_consumed_mmol_per_l == pytest.approx(1.0, abs=1e-12)
        assert run.product_mmol_per_l[-1] == pytest.approx(0.5, abs=1e-12)

    def test_an_inactive_feed_window_is_identical_to_no_feed(self):
        batch = _run(steps=2)
        fed = _run(feed=FeedProfile(0.0, 1.23, 4.56), steps=2)

        np.testing.assert_array_equal(fed.biomass_g_per_l, batch.biomass_g_per_l)
        np.testing.assert_array_equal(fed.substrate_mmol_per_l, batch.substrate_mmol_per_l)
        np.testing.assert_array_equal(fed.product_mmol_per_l, batch.product_mmol_per_l)

    def test_a_depleted_culture_restarts_only_when_feed_arrives_and_stops_after_it_ends(self):
        def rates(substrate):
            return (0.0, 2.0, 0.5) if substrate > 0.0 else (0.0, 0.0, 0.0)

        run = simulate_batch(rates, 1.0, 0.5, 2.0, 24.0,
                             feed=FeedProfile(4.0, 0.75, 1.0), steps=4)

        np.testing.assert_allclose(run.substrate_mmol_per_l, [0.5, 0.0, 0.5, 0.0, 0.0], atol=1e-8)
        np.testing.assert_allclose(run.product_mmol_per_l, [0.0, 0.125, 0.25, 0.375, 0.375], atol=1e-8)
        assert run.substrate_consumed_mmol_per_l == pytest.approx(1.5, abs=1e-8)
        assert run.total_feed_mmol_per_l == 1.0

    def test_a_discontinuous_starvation_switch_cannot_hide_submaintenance_feed_in_averaged_rates(self):
        def rates(substrate):
            return (0.0, 1.0, 0.0) if substrate > 0.25 else (0.0, 0.0, 0.0)

        with pytest.raises(RuntimeError, match="starvation"):
            simulate_batch(rates, 1.0, 1.0, 3.0, 24.0, feed=FeedProfile(0.1), steps=2)
