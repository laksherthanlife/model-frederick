"""Anaerobic and hypoxic mode: the diagnosis, the supplements, and Verduyn 1990.

Vendored Yeast9 grows at exactly zero without oxygen and opening ergosterol does not fix
it. What these tests hold in place is the full list of what does, the two numbers the
result is checked against, and the two places the model departs from a real anaerobic
culture -- no glycerol, and a yield 7.2% high.
"""

import pytest

from ystwin.fba import stress_oxygen as so
from ystwin.fba.physiology import cap_uptake
from ystwin.fba.solver import FVA_PROCESSES, configure, growth_or_none

pytestmark = pytest.mark.integration

AEROBIC_GROWTH = 0.887685
ANAEROBIC_GROWTH = 0.199918
HEME_A_FORMULA = {"C": 49, "H": 55, "Fe": 1, "N": 4, "O": 6}
HEME_A_CHARGE = -3
HEME_A_COEFFICIENT = 1e-06
# mu / qO2 in the oxygen-limited regime of the unsupplemented model, /h per mmol/gDCW/h.
HYPOXIC_SLOPE = 1.8246
# Oleic acid C18H34O2, for Verduyn 1990's 35 mg (g biomass)-1.
OLEATE_G_PER_MMOL = 0.28246


@pytest.fixture(scope="module")
def _gem_template(yeast_gem_factory):
    model = yeast_gem_factory()
    configure(model)
    return model


@pytest.fixture
def gem(_gem_template, model_copy):
    return model_copy(_gem_template)


@pytest.fixture(scope="module")
def _anaerobic_template(_gem_template, model_copy):
    return so.anaerobic_model(model_copy(_gem_template), glucose_uptake=10.0)


@pytest.fixture
def anaerobic(_anaerobic_template, model_copy):
    return model_copy(_anaerobic_template)


class TestTheDiagnosis:
    def test_the_vendored_model_cannot_grow_without_oxygen(self, gem):
        with gem as m:
            cap_uptake(m, "r_1714", 10.0)
            assert growth_or_none(m) == pytest.approx(AEROBIC_GROWTH, rel=1e-5)
            m.reactions.get_by_id("r_1992").bounds = (0.0, 0.0)

            assert growth_or_none(m) == pytest.approx(0.0, abs=1e-9)

    def test_opening_ergosterol_alone_does_not_fix_it(self, gem):
        """The 1953 requirement is real and it is not the whole requirement."""
        with gem as m:
            cap_uptake(m, "r_1714", 10.0)
            m.reactions.get_by_id("r_1992").bounds = (0.0, 0.0)
            cap_uptake(m, "r_1757", 1000.0)

            assert growth_or_none(m) == pytest.approx(0.0, abs=1e-9)

    def test_the_recorded_essential_set_is_what_deletion_actually_finds(self, gem):
        measured = {rid for rid, _, _ in so.diagnose_oxygen_essentiality(gem, 10.0)}

        assert measured == {r.reaction_id for r in so.ESSENTIAL_OXYGEN_REACTIONS}

    def test_the_essential_set_does_not_depend_on_the_glucose_bound(self, gem):
        """It is a stoichiometric requirement, so it should survive a change of regime."""
        measured = {rid for rid, _, _ in so.diagnose_oxygen_essentiality(gem, 11.1)}

        assert measured == {r.reaction_id for r in so.ESSENTIAL_OXYGEN_REACTIONS}

    def test_the_essential_set_is_a_small_part_of_the_oxygen_consuming_set(self, gem):
        assert len(so.oxygen_consuming_reactions(gem)) == 74
        assert len(so.ESSENTIAL_OXYGEN_REACTIONS) == 19

    def test_every_requirement_names_what_resolves_it(self):
        for requirement in so.ESSENTIAL_OXYGEN_REACTIONS:
            assert requirement.requirement
            assert requirement.resolved_by


class TestWhichSpeciesCountAsOxygen:
    def test_it_resolves_to_the_recorded_set(self, gem):
        assert so.oxygen_metabolite_ids(gem) == tuple(
            sorted(so.MOLECULAR_OXYGEN_METABOLITES))

    def test_superoxide_shares_the_formula_and_is_excluded_by_charge(self, gem):
        for mid in ("s_3813", "s_3931"):
            superoxide = gem.metabolites.get_by_id(mid)

            assert superoxide.formula == "O2"
            assert superoxide.charge == -1
            assert mid not in so.oxygen_metabolite_ids(gem)

    def test_the_reactions_that_dispose_of_superoxide_are_not_called_oxygen_consumers(
            self, gem):
        """r_4190/r_4270 PRODUCE O2 from superoxide; counting them would invert the sign."""
        consumers = so.oxygen_consuming_reactions(gem)

        assert "r_4190" not in consumers
        assert "r_4270" not in consumers


class TestTheAnaerobicCofactorPseudoreaction:
    def test_it_differs_from_the_vendored_one_by_exactly_heme_a(self, gem, anaerobic):
        vendored = gem.reactions.get_by_id(so.COFACTOR_PSEUDOREACTION).check_mass_balance()
        added = anaerobic.reactions.get_by_id(
            so.ANAEROBIC_COFACTOR_PSEUDOREACTION).check_mass_balance()

        for element, count in HEME_A_FORMULA.items():
            delta = added.get(element, 0.0) - vendored.get(element, 0.0)
            assert delta == pytest.approx(count * HEME_A_COEFFICIENT, rel=1e-6)
        assert added.get("charge") - vendored["charge"] == pytest.approx(
            HEME_A_CHARGE * HEME_A_COEFFICIENT, rel=1e-6)

    def test_it_is_a_pseudoreaction_so_it_does_not_balance_and_neither_does_the_original(
            self, gem, anaerobic):
        """s_4205 has no formula, so both reactions are imbalanced by construction."""
        assert gem.metabolites.get_by_id("s_4205").formula in (None, "")
        assert gem.reactions.get_by_id(so.COFACTOR_PSEUDOREACTION).check_mass_balance() != {}
        assert anaerobic.reactions.get_by_id(
            so.ANAEROBIC_COFACTOR_PSEUDOREACTION).check_mass_balance() != {}

    def test_it_carries_no_iron_and_the_vendored_one_does(self, gem, anaerobic):
        assert "Fe" in gem.reactions.get_by_id(so.COFACTOR_PSEUDOREACTION).check_mass_balance()
        assert "Fe" not in anaerobic.reactions.get_by_id(
            so.ANAEROBIC_COFACTOR_PSEUDOREACTION).check_mass_balance()

    def test_it_carries_flux_and_the_vendored_one_is_closed(self, anaerobic):
        from cobra.flux_analysis import flux_variability_analysis

        span = flux_variability_analysis(
            anaerobic, reaction_list=[so.ANAEROBIC_COFACTOR_PSEUDOREACTION],
            fraction_of_optimum=0.0, processes=FVA_PROCESSES,
        ).loc[so.ANAEROBIC_COFACTOR_PSEUDOREACTION]

        assert span["maximum"] > 1e-6
        assert anaerobic.reactions.get_by_id(so.COFACTOR_PSEUDOREACTION).bounds == (0.0, 0.0)

    def test_installing_it_twice_is_refused(self, anaerobic):
        with pytest.raises(ValueError, match="already"):
            so.add_anaerobic_cofactor_pseudoreaction(anaerobic)


class TestTheAnaerobicMode:
    def test_it_grows(self, anaerobic):
        assert growth_or_none(anaerobic) == pytest.approx(ANAEROBIC_GROWTH, rel=1e-4)

    def test_no_oxygen_is_taken_up(self, anaerobic):
        solution = anaerobic.optimize()

        assert solution.fluxes["r_1992"] == pytest.approx(0.0, abs=1e-9)

    def test_every_supplement_is_individually_necessary(self, gem):
        for dropped in so.ANAEROBIC_SUPPLEMENTS:
            kept = tuple(s for s in so.ANAEROBIC_SUPPLEMENTS if s is not dropped)
            model = so.anaerobic_model(gem, glucose_uptake=10.0, supplements=kept)

            assert growth_or_none(model) == pytest.approx(0.0, abs=1e-9), dropped.component

    def test_the_unsupplemented_medium_does_not_grow(self, gem):
        model = so.anaerobic_model(gem, glucose_uptake=10.0, supplements=())

        assert growth_or_none(model) == pytest.approx(0.0, abs=1e-9)

    def test_every_supplement_carries_a_citation(self):
        for supplement in so.ANAEROBIC_SUPPLEMENTS:
            assert "PMID" in supplement.source or "USP-NF" in supplement.source

    def test_importing_heme_a_gives_the_same_growth_as_taking_it_out_of_the_biomass(self, gem):
        """Both routes are available and they agree, so the default is the honest one."""
        imported = so.anaerobic_model(gem, glucose_uptake=10.0, import_heme_a=True)
        removed = so.anaerobic_model(gem, glucose_uptake=10.0)

        assert growth_or_none(imported) == pytest.approx(growth_or_none(removed), rel=1e-9)

    def test_the_default_does_not_open_the_heme_a_exchange(self, anaerobic):
        assert anaerobic.reactions.get_by_id("r_4780").lower_bound == 0.0

    def test_the_model_it_is_given_is_not_modified(self, gem):
        before = growth_or_none(gem)
        so.anaerobic_model(gem, glucose_uptake=10.0)

        assert growth_or_none(gem) == pytest.approx(before, rel=1e-12)
        assert not gem.reactions.has_id(so.ANAEROBIC_COFACTOR_PSEUDOREACTION)


class TestTheEsteraseDirectionConstraint:
    def test_it_changes_nothing_aerobically(self, gem):
        with gem as m:
            cap_uptake(m, "r_1714", 10.0)
            before = growth_or_none(m)
            so.constrain_esterases_to_hydrolysis(m)

            assert growth_or_none(m) == pytest.approx(before, rel=1e-12)

    def test_without_it_the_anaerobic_optimum_leaks_carbon_into_an_ester(self, gem):
        free = so.anaerobic_model(gem, glucose_uptake=10.0, hydrolysis_only_esterases=False)
        solution = free.optimize()

        assert solution.fluxes["r_4731"] > 0.5

    def test_with_it_the_ester_sink_is_gone_and_the_carbon_leaves_as_ethanol(self, anaerobic):
        solution = anaerobic.optimize()

        assert solution.fluxes["r_4731"] == pytest.approx(0.0, abs=1e-9)
        assert solution.fluxes["r_1761"] > 17.0

    def test_it_constrains_the_reversible_esterases_and_not_the_irreversible_ones(
            self, anaerobic):
        for rid in so.REVERSIBLE_ESTERASES:
            assert anaerobic.reactions.get_by_id(rid).upper_bound == 0.0
        for rid in ("r_0369", "r_0656", "r_0657"):
            assert anaerobic.reactions.get_by_id(rid).bounds == (0.0, 1000.0)


class TestTheMissingGlycerol:
    def test_the_optimum_secretes_none_although_the_route_is_open(self, anaerobic):
        from cobra.flux_analysis import flux_variability_analysis

        at_optimum = flux_variability_analysis(
            anaerobic, reaction_list=[so.GLYCEROL_EXCHANGE], fraction_of_optimum=1.0,
            processes=FVA_PROCESSES).loc[so.GLYCEROL_EXCHANGE]
        unconstrained = flux_variability_analysis(
            anaerobic, reaction_list=[so.GLYCEROL_EXCHANGE], fraction_of_optimum=0.0,
            processes=FVA_PROCESSES).loc[so.GLYCEROL_EXCHANGE]

        assert at_optimum["maximum"] == pytest.approx(0.0, abs=1e-9)
        assert unconstrained["maximum"] > 10.0

    def test_closing_both_nadh_ammonium_routes_does_not_bring_it_back(self, gem):
        """The obvious escape route, checked rather than assumed: it costs nothing."""
        blocked = so.anaerobic_model(gem, glucose_uptake=10.0)
        for rid in ("r_0470", "r_0472"):
            blocked.reactions.get_by_id(rid).bounds = (0.0, 0.0)

        assert growth_or_none(blocked) == pytest.approx(ANAEROBIC_GROWTH, rel=1e-4)
        assert blocked.optimize().fluxes[so.GLYCEROL_EXCHANGE] == pytest.approx(0.0, abs=1e-9)

    def test_an_imposed_glycerol_flux_is_priced_rather_than_invented(self, gem):
        priced = so.anaerobic_glycerol_price(gem, fluxes=(0.5, 1.0, 2.0, 5.0))
        costs = [cost for _, _, cost in priced]

        assert costs == pytest.approx([-0.0187, -0.0438, -0.1073, -0.4203], abs=5e-4)
        assert costs == sorted(costs, reverse=True)

    def test_the_price_list_is_measured_against_the_regime_it_names(self, gem):
        priced = so.anaerobic_glycerol_price(gem, fluxes=(1.0,))
        _, growth, _ = priced[0]

        assert growth == pytest.approx(0.191153, rel=1e-4)


class TestAgainstVerduyn1990:
    @pytest.fixture(scope="class")
    def _report_template(self, _gem_template, model_copy):
        return so.validate_anaerobic(model_copy(_gem_template))

    @pytest.fixture
    def report(self, _report_template):
        return _report_template

    def test_the_reference_matches_the_constant_the_generator_already_carries(self):
        from ystwin.generator.context import ANAEROBIC_MU_MAX_PER_H

        assert so.VERDUYN_1990_ANAEROBIC.mu_max == ANAEROBIC_MU_MAX_PER_H

    def test_the_biomass_yield_at_the_measured_dilution_rate_agrees(self, report):
        assert report.model_biomass_yield == pytest.approx(0.1072, rel=1e-3)
        assert abs(report.yield_relative_error) < 0.15

    def test_reaching_the_measured_mu_max_needs_a_reported_glucose_uptake(self, report):
        assert report.glucose_for_reference_mu_max == pytest.approx(15.314, rel=1e-3)

    def test_the_implied_uptake_is_arithmetic_on_the_papers_own_two_numbers(self):
        reference = so.VERDUYN_1990_ANAEROBIC

        assert reference.implied_glucose_uptake() == pytest.approx(
            reference.mu_max / (reference.biomass_yield * 0.18016), rel=1e-12)
        assert reference.implied_glucose_uptake() == pytest.approx(17.2069, rel=1e-4)

    def test_the_growth_at_that_uptake_overshoots_the_measured_mu_max(self, report):
        assert report.mu_at_implied_uptake == pytest.approx(0.349223, rel=1e-4)
        assert report.mu_relative_error == pytest.approx(0.1265, abs=1e-3)

    def test_the_ethanol_yield_there_is_below_the_gay_lussac_maximum(self, report):
        gay_lussac = 2 * 0.04607 / 0.18016

        assert report.ethanol_yield_g_per_g == pytest.approx(0.4444, rel=1e-3)
        assert report.ethanol_yield_g_per_g < gay_lussac

    def test_the_report_carries_the_missing_glycerol_rather_than_hiding_it(self, report):
        assert report.glycerol_at_optimum == pytest.approx(0.0, abs=1e-9)
        assert "glycerol" in report.summary()

    def test_the_summary_names_its_source(self, report):
        assert "PMID 1975265" in report.summary()

    def test_a_growth_rate_the_regime_cannot_reach_is_refused(self, anaerobic):
        with pytest.raises(RuntimeError, match="infeasible"):
            so.biomass_yield(anaerobic, 5.0)


class TestTheMaintenanceDiagnostic:
    def test_the_vendored_maintenance_is_what_the_module_says_it_is(self, gem):
        assert gem.reactions.get_by_id(so.NGAM_REACTION).bounds == (
            so.VENDORED_NGAM, so.VENDORED_NGAM)

    def test_a_larger_maintenance_moves_the_yield_towards_the_measured_one(self, gem):
        rows = so.maintenance_sensitivity(gem)
        yields = [y for _, y, _ in rows]

        assert yields == pytest.approx([0.1072, 0.1042, 0.0995, 0.0953], abs=5e-4)
        assert yields == sorted(yields, reverse=True)

    def test_it_does_not_adopt_anything(self, gem, anaerobic):
        """A diagnostic that changed the model would be fitting, so it must not."""
        so.maintenance_sensitivity(gem)

        assert gem.reactions.get_by_id(so.NGAM_REACTION).bounds == (0.7, 0.7)
        assert anaerobic.reactions.get_by_id(so.NGAM_REACTION).bounds == (0.7, 0.7)


class TestHypoxia:
    def test_growth_is_continuous_in_oxygen_with_no_supplements_at_all(self, gem):
        rates = [growth_or_none(so.hypoxic_model(gem, oxygen_uptake=q, glucose_uptake=10.0))
                 for q in (0.001, 0.1, 1.0, 5.0)]

        assert rates == pytest.approx([0.001825, 0.182462, 0.229329, 0.362420], rel=1e-3)
        assert rates == sorted(rates)

    def test_below_the_knee_growth_is_exactly_proportional_to_oxygen(self, gem):
        """Oxygen there is a biosynthetic reagent, not an energy source."""
        slopes = [growth_or_none(so.hypoxic_model(gem, oxygen_uptake=q, glucose_uptake=10.0)) / q
                  for q in (0.001, 0.01, 0.05, 0.1)]

        assert slopes == pytest.approx([HYPOXIC_SLOPE] * 4, rel=1e-3)

    def test_above_the_knee_the_proportionality_is_gone(self, gem):
        slope = growth_or_none(so.hypoxic_model(gem, oxygen_uptake=0.5,
                                                glucose_uptake=10.0)) / 0.5

        assert slope < 0.5 * HYPOXIC_SLOPE

    def test_supplementing_lifts_the_most_oxygen_limited_point_a_hundredfold(self, gem):
        bare, supplemented, ratio = so.hypoxic_supplement_effect(gem, 0.001)

        assert bare == pytest.approx(0.001825, rel=1e-3)
        assert supplemented == pytest.approx(0.202429, rel=1e-3)
        assert ratio == pytest.approx(110.9, rel=1e-2)

    def test_the_supplements_stop_mattering_once_oxygen_is_not_the_limit(self, gem):
        _, _, ratio = so.hypoxic_supplement_effect(gem, 1.0)

        assert ratio == pytest.approx(1.089, rel=1e-2)

    def test_zero_oxygen_is_refused_rather_than_returning_zero_growth(self, gem):
        with pytest.raises(ValueError, match="anaerobic, not hypoxic"):
            so.hypoxic_model(gem, oxygen_uptake=0.0)

    def test_the_model_it_is_given_is_not_modified(self, gem):
        before = growth_or_none(gem)
        so.hypoxic_model(gem, oxygen_uptake=1.0)

        assert growth_or_none(gem) == pytest.approx(before, rel=1e-12)


class TestWhatTheVerifierMeasured:
    """Three facts found while trying to break this module, pinned so they cannot be lost."""

    def test_the_missing_glycerol_is_a_transhydrogenase_cycle_that_can_be_named(self, gem):
        """GPD -> r_0489 phosphatase -> GCY1 r_0487 makes the glycerol and re-eats it."""
        base = so.anaerobic_model(gem, glucose_uptake=10.0)
        blocked = so.anaerobic_model(gem, glucose_uptake=10.0)
        for rid in ("r_0487", "r_0470", "r_0472"):
            blocked.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        growth = growth_or_none(blocked)

        assert base.optimize().fluxes["r_0487"] > 0.3
        assert blocked.optimize().fluxes[so.GLYCEROL_EXCHANGE] == pytest.approx(0.418, rel=1e-2)
        assert growth == pytest.approx(0.185951, rel=1e-4)
        assert growth / growth_or_none(base) - 1.0 == pytest.approx(-0.0699, abs=1e-3)

    def test_the_unsaturated_fatty_acid_draw_matches_verduyns_measured_requirement(self, gem):
        """Verduyn 1990 PMID 1975265: ~35 mg oleic acid per g biomass for optimal growth."""
        model = so.anaerobic_model(gem, glucose_uptake=10.0)
        growth = growth_or_none(model)
        solution = model.optimize()
        ufa = -(solution.fluxes["r_2189"] + solution.fluxes["r_1994"])

        assert ufa / growth * OLEATE_G_PER_MMOL * 1000 == pytest.approx(35.0, rel=0.05)

    def test_the_supplements_are_a_requirement_and_not_fuel_in_the_range_this_module_uses(
            self, gem):
        """Unbounded lipid uptake becomes food once oxygen is free; here it must not be."""
        carbon = {"r_1757": 28, "r_2189": 18, "r_1994": 16, "r_1967": 6, "r_1548": 9}
        for oxygen in (0.001, 5.0):
            model = so.hypoxic_model(gem, oxygen_uptake=oxygen, glucose_uptake=10.0,
                                     supplements=so.ANAEROBIC_SUPPLEMENTS)
            fluxes = model.optimize().fluxes
            supplement_carbon = sum(-fluxes[e] * n for e, n in carbon.items())

            assert supplement_carbon / (-fluxes["r_1714"] * 6) < 0.03, oxygen
