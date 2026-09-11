"""The environment axes, and the dataset that spans them.

Everything the project has fitted was generated in one corner: exponential-phase glucose, 30
degrees, air, one pH. Both channels a plate reader gives -- optical density and reporter -- move
with growth rate, and in that corner the only thing that moves the growth rate is the stressor.
So the fitted model has no way to distinguish "slowed" from "stressed", and there is no sample
size at which it learns to.

These tests cover the three parts of the fix: the environment terms that were wrong or missing,
the refusals that stop a wide sweep from inventing what it cannot ground, and the sweep design
itself, which has to record every condition on every row or a downstream fit cannot condition
on any of them.
"""

import numpy as np
import pytest

from ystwin.generator.context import (
    ANAEROBIC_MU_MAX_PER_H,
    CARDINAL_TEMPERATURES_C,
    CITRINE_PKA,
    CultureContext,
    citrine_ph_response,
    context_growth_rate,
    context_physiology,
    ph_growth_factor,
    temperature_factor,
)
from ystwin.generator.design import (
    collinearity,
    context_grid,
    environment_grid,
    stressor_only_series,
)
from ystwin.generator.literature import parameters_from_literature

PARAMS = parameters_from_literature("NativeYap1")


class TestTemperatureIsNoLongerASymmetricLine:
    """`1 - 0.04 * |T - 30|` said 20 and 40 degrees cost a culture the same 40%.

    Yeast is not symmetric about its reference temperature and does not decline linearly. The
    cardinal temperature model of Rosso 1993 (PMID 8412234) at Salvado 2011's measured cardinals
    for *S. cerevisiae* (PMID 21317255, T_opt 32.3, T_max 45.4) says 20 degrees is a mild
    slowdown and 40 degrees is most of the way to a dead culture. The two disagree most at
    exactly the temperatures a heat-shock experiment is run at.
    """

    def test_thirty_degrees_is_still_exactly_the_reference(self):
        """Every carbon-source rate in this module was read from a 30-degree culture.

        The temperature term has to be a shape normalised at 30, not the model's own mu_opt.
        Taking mu_opt would rescale every existing result in the repository silently.
        """
        assert temperature_factor(30.0) == pytest.approx(1.0)

    def test_the_optimum_is_warmer_than_the_reference(self):
        assert temperature_factor(CARDINAL_TEMPERATURES_C["opt"]) > temperature_factor(30.0)

    def test_growth_stops_at_and_above_the_maximum_temperature(self):
        for temperature in (CARDINAL_TEMPERATURES_C["max"], 46.0, 50.0):
            assert temperature_factor(temperature) == 0.0

    def test_growth_stops_at_and_below_the_minimum_temperature(self):
        for temperature in (CARDINAL_TEMPERATURES_C["min"], 1.0, 0.0):
            assert temperature_factor(temperature) == 0.0

    def test_the_curve_is_asymmetric_about_the_reference(self):
        """The specific defect the old form had.

        Ten degrees warmer than 30 is nearer the optimum than ten degrees cooler, so 40 must
        cost less than 20 does. The old linear form scored them identically at 0.60.
        """
        assert temperature_factor(40.0) > temperature_factor(20.0)

    def test_thirty_seven_degrees_is_a_mild_slowdown_not_a_third_of_growth(self):
        """The number this changes most, and the one the project actually uses.

        37 degrees is the standard heat-stress condition. The old linear form put a culture
        there at 0.72 of maximum; the measured cardinals put it at 0.92, because 37 is still on
        the near side of a 45.4-degree ceiling.
        """
        assert temperature_factor(37.0) == pytest.approx(0.92, abs=0.02)

    def test_forty_five_degrees_is_almost_no_growth(self):
        """The other end of the same error, and the larger one.

        The old form gave a culture 0.40 of maximum growth at 45 degrees. The measured maximum
        is 45.4, so the true answer is a culture on the point of not growing at all.
        """
        assert temperature_factor(45.0) < 0.1


class TestAnoxiaIsAMeasurementNotAnExtrapolation:
    """The old model took an unsupplemented anaerobic glucose culture to 0.28 /h.

    It got there by multiplying 0.40 by one minus an asserted `respiring` weight of 0.3. The
    measured answer is zero: sterol and unsaturated fatty acid synthesis both need molecular
    oxygen, so a defined medium without them supports no anaerobic growth at all (Andreasen &
    Stier 1953, PMID 13034889). With them, Verduyn 1990 (PMID 1975265) measured mu_max 0.31 /h.
    """

    def test_an_unsupplemented_anaerobic_culture_does_not_grow(self):
        assert context_growth_rate(CultureContext(oxygen=0.0)) == 0.0

    def test_a_supplemented_anaerobic_culture_grows_at_the_measured_rate(self):
        rate = context_growth_rate(CultureContext(oxygen=0.0, anaerobic_supplements=True))

        assert rate == pytest.approx(ANAEROBIC_MU_MAX_PER_H)

    def test_supplements_do_not_rescue_a_non_fermentable_substrate(self):
        """Ergosterol replaces what oxygen was needed to *build*, not what it was needed to do.

        A culture on ethanol has to oxidise it, so no supplement makes anaerobic growth
        possible. Treating the supplement flag as a general anaerobic rescue would produce a
        culture respiring without a terminal electron acceptor.
        """
        context = CultureContext(carbon_source="ethanol", oxygen=0.0,
                                 anaerobic_supplements=True)

        assert context_growth_rate(context) == 0.0

    def test_the_anaerobic_ceiling_never_exceeds_the_aerobic_rate(self):
        """0.31 /h is a glucose measurement and the table also holds galactose.

        Applying it as a value rather than a ceiling would make an anaerobic galactose culture
        grow faster than an aerated one, which nothing supports.
        """
        aerobic = context_growth_rate(CultureContext(carbon_source="galactose"))
        anaerobic = context_growth_rate(
            CultureContext(carbon_source="galactose", oxygen=0.0, anaerobic_supplements=True))

        assert anaerobic <= aerobic

    def test_full_air_is_unchanged_by_any_of_this(self):
        """No other lane's numbers may move because anoxia was fixed."""
        assert context_growth_rate(CultureContext()) == pytest.approx(0.40)


class TestPhIsRecordedAndItsGrowthEffectIsRefused:
    """pH is a real condition and this module will not pretend to know what it does to growth.

    No cardinal pH values for *S. cerevisiae* were found, and Arroyo-Lopez 2009 (PMID 19246112)
    found pH *not* significant for *S. cerevisiae* where temperature was. Returning 1.0 would
    assert that pH does not matter in every sweep that varies it; returning an invented curve
    would be worse.
    """

    def test_a_context_carries_both_a_medium_and_a_cytosolic_ph(self):
        """They are different numbers and confusing them predicts a dark cell that is not dark.

        Yeast holds a near-neutral cytosol against a medium at 4. The reporter sees the cytosol.
        """
        context = CultureContext(ph_medium=4.0, ph_cytosolic=6.8)

        assert context.ph_medium == 4.0
        assert context.ph_cytosolic == 6.8

    def test_asking_for_a_growth_response_to_ph_without_measured_cardinals_is_refused(self):
        with pytest.raises(ValueError, match="refuses to invent"):
            ph_growth_factor(4.0)

    def test_supplying_measured_cardinals_gives_the_rosso_cardinal_ph_model(self):
        cardinals = {"min": 2.5, "opt": 5.5, "max": 8.5}

        assert ph_growth_factor(5.5, cardinals) == pytest.approx(1.0)
        assert ph_growth_factor(2.5, cardinals) == 0.0
        assert ph_growth_factor(8.5, cardinals) == 0.0
        assert 0.0 < ph_growth_factor(4.0, cardinals) < 1.0

    def test_growth_rate_does_not_secretly_depend_on_ph(self):
        """The refusal is worth nothing if the rate quietly uses pH anyway."""
        neutral = context_growth_rate(CultureContext(ph_medium=5.5))
        acid = context_growth_rate(CultureContext(ph_medium=3.5))

        assert neutral == acid

    def test_an_impossible_ph_is_refused_at_construction(self):
        with pytest.raises(ValueError, match="ph_medium"):
            CultureContext(ph_medium=-1.0)


class TestThePhEffectThatIsMeasuredIsOnTheReporter:
    """Citrine's pKa is 5.7 (Griesbeck 2001, PMID 11387331), and that is not a rounding error.

    A full unit of cytosolic acidification from 7.0 to 6.0 removes 30% of the fluorescence with
    the promoter untouched. The induction folds this project measures are around 1.5, so the
    artifact is a third of the signal -- and it has dose structure, because the stressors that
    acidify are the ones being dosed.
    """

    def test_the_response_is_one_half_at_the_measured_pka(self):
        assert citrine_ph_response(CITRINE_PKA) == pytest.approx(0.5)

    def test_a_neutral_cytosol_is_nearly_fully_bright(self):
        assert citrine_ph_response(7.0) == pytest.approx(0.95, abs=0.01)

    def test_one_unit_of_acidification_costs_about_a_third_of_the_signal(self):
        loss = 1.0 - citrine_ph_response(6.0) / citrine_ph_response(7.0)

        assert loss == pytest.approx(0.30, abs=0.02)

    def test_the_response_rises_monotonically_with_ph(self):
        values = [citrine_ph_response(p) for p in np.linspace(4.0, 9.0, 40)]

        assert all(a < b for a, b in zip(values, values[1:]))


class TestTheContextKnowsWhatMetabolismItIsIn:
    def test_the_default_glucose_culture_is_respiro_fermentative(self):
        """0.40 /h is above the critical rate, so the corner everything was trained in ferments.

        Worth pinning because it is easy to picture the reference condition as the tidy
        respiratory one. It is not: the untreated well is the *fermenting* one, and dosing it
        is what makes it respire.
        """
        physiology = context_physiology(CultureContext())

        assert physiology.fermentative
        assert physiology.biomass_yield == pytest.approx(0.20)

    def test_a_context_outside_the_measured_curve_refuses_rather_than_guessing(self):
        """A stationary-phase culture grows at 0.02 /h, below anything van Hoek ran."""
        with pytest.raises(ValueError, match="outside the measured chemostat curve"):
            context_physiology(CultureContext(growth_phase="stationary"))


class TestTheSweepDesign:
    def test_the_default_grid_is_the_one_corner_everything_was_trained_in(self):
        """Widening one axis has to be a deliberate act, not the default.

        If the defaults were wide, a caller who varied one thing would be varying six, and no
        result from the sweep would be attributable to the axis they meant to move.
        """
        contexts = context_grid()

        assert contexts == [CultureContext()]

    def test_it_is_a_full_factorial_over_the_axes_given(self):
        contexts = context_grid(carbon_sources=("glucose", "ethanol"),
                                temperatures_c=(25.0, 30.0, 37.0),
                                oxygen_fractions=(0.21, 0.05))

        assert len(contexts) == 12
        assert len({(c.carbon_source, c.temperature_c, c.oxygen) for c in contexts}) == 12

    def test_an_unsupported_axis_value_is_refused_rather_than_dropped(self):
        """Returning a smaller grid than was asked for would hide a design error silently."""
        with pytest.raises(ValueError, match="carbon source"):
            context_grid(carbon_sources=("glucose", "xylose"))

    def test_every_row_records_every_condition_that_produced_it(self):
        """The whole point of the sweep. A condition not written down cannot be conditioned on.

        Growth rate alone is not enough: two environments can produce the same rate for
        different reasons -- a cool aerated culture and a warm starved one -- and a model handed
        only the rate cannot tell them apart.
        """
        contexts = context_grid(temperatures_c=(30.0, 37.0), oxygen_fractions=(0.21, 0.05))
        conditions = environment_grid(PARAMS, contexts, doses=(0.0, 1.0))

        for condition in conditions:
            row = condition.to_row()
            for field in ("carbon_source", "growth_phase", "temperature_c", "oxygen",
                          "ph_medium", "ph_cytosolic", "glucose_g_per_L",
                          "anaerobic_supplements", "dose_mM", "nutrient_factor",
                          "growth_rate", "promoter_activity"):
                assert field in row

    def test_a_refusal_travels_as_data_rather_than_stopping_the_sweep(self):
        """A wide sweep visits corners the measured curve does not cover, and that is a result.

        Raising would mean the sweep stops at the first slow culture and produces nothing. The
        column says which cells are grounded, so the fraction that are is reportable instead of
        being silently the whole thing.
        """
        contexts = context_grid(growth_phases=("exponential", "stationary"))
        conditions = environment_grid(PARAMS, contexts, doses=(0.0,))

        grounded = [c for c in conditions if not c.physiology_refusal]
        refused = [c for c in conditions if c.physiology_refusal]
        assert grounded and refused
        for condition in refused:
            assert condition.physiology is None
            assert condition.to_row()["q_ethanol"] is None
        for condition in grounded:
            assert condition.physiology is not None

    def test_the_environment_axes_break_the_stress_growth_collinearity(self):
        """The reason for doing any of this, measured on the same scale as the nutrient axis.

        A dose ladder in one environment is almost perfectly collinear, so stress and 1/mu are
        one axis. Crossing the same ladder with environments that move the growth rate for
        reasons unrelated to the promoter breaks that.

        Three grids rather than two, because the size of the effect is the finding and it is
        not automatic. A narrow environment grid -- three temperatures and two aeration levels
        -- buys only about 0.19 of correlation, which is real and is nowhere near enough. It
        takes a grid that spans carbon source, temperature, aeration and the nutrient axis
        together to get the correlation down to where a latent state is separable. Asserting
        only "wide beats one" would let a future narrowing of the sweep pass unnoticed.
        """
        doses = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0)
        one_environment = stressor_only_series(PARAMS, doses=doses)
        narrow = environment_grid(
            PARAMS,
            context_grid(temperatures_c=(25.0, 30.0, 37.0), oxygen_fractions=(0.21, 0.05)),
            doses=doses,
        )
        wide = environment_grid(
            PARAMS,
            context_grid(carbon_sources=("glucose", "galactose"),
                         temperatures_c=(20.0, 25.0, 30.0, 34.0, 40.0),
                         oxygen_fractions=(0.21, 0.10, 0.05)),
            doses=doses,
            nutrient_factors=(1.0, 0.6),
        )

        assert collinearity(one_environment) > 0.85
        assert collinearity(narrow) < collinearity(one_environment)
        assert collinearity(wide) < 0.5

    def test_the_dose_acts_as_a_fraction_of_whatever_the_environment_allows(self):
        """The composition rule, asserted and therefore worth pinning where it can be seen.

        The same dose in a cooler incubator has to remove the same *share* of a smaller rate.
        If it removed the same absolute rate it could drive a slow culture negative, and if it
        were ignored the environment axis would be the only one moving growth.
        """
        warm = environment_grid(PARAMS, context_grid(temperatures_c=(30.0,)), doses=(1.0,))[0]
        cool = environment_grid(PARAMS, context_grid(temperatures_c=(25.0,)), doses=(1.0,))[0]
        undosed_warm = environment_grid(
            PARAMS, context_grid(temperatures_c=(30.0,)), doses=(0.0,))[0]
        undosed_cool = environment_grid(
            PARAMS, context_grid(temperatures_c=(25.0,)), doses=(0.0,))[0]

        assert (warm.growth_rate / undosed_warm.growth_rate ==
                pytest.approx(cool.growth_rate / undosed_cool.growth_rate))

    def test_the_promoter_response_does_not_move_with_the_environment(self):
        """Context shifts the baseline and the growth rate, not the response.

        The generator has taken this position since `context.py` was written, and the sweep must
        not quietly break it -- if the environment changed the dose-response too, no stressor
        effect learned in one environment would mean anything in another.
        """
        contexts = context_grid(temperatures_c=(25.0, 37.0), oxygen_fractions=(0.21, 0.05))
        conditions = environment_grid(PARAMS, contexts, doses=(0.5,))

        activities = {condition.promoter_activity for condition in conditions}
        assert len(activities) == 1

    def test_the_context_baseline_travels_with_the_row(self):
        """A stationary well carries a raised general stress response with no agent in it.

        That is the confound the whole context layer exists for, so it has to reach the dataset
        rather than being recomputed by whoever fits it.
        """
        conditions = environment_grid(
            PARAMS, context_grid(growth_phases=("stationary",)), doses=(0.0,))

        assert conditions[0].to_row()["baseline_ESR"] > 0.3
