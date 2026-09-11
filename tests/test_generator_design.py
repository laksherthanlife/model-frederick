"""Experiment design: break the collinearity between stress and growth.

Every perturbation run so far slows growth, so "stress" and "1/mu" are one axis and no
latent stress state is separable from growth rate at any sample size. The fix is a
design, not a model: perturbations that slow growth without stressing (nutrient
limitation) and that stress without slowing growth (sub-inhibitory dose).
"""

import numpy as np
import pytest

from ystwin.generator.culture import CultureParameters, simulate_culture
from ystwin.generator.design import (
    Condition,
    collinearity,
    decoupling_grid,
    stressor_only_series,
)

PARAMS = CultureParameters(
    mu_max=0.40, growth_ic50=0.8, growth_hill=1.5,
    promoter_basal=1.0e-3, promoter_peak=2.5e-3, promoter_ec50=0.7, promoter_hill=2.0,
    carrying_capacity=3.0,
)


class TestNutrientLimitationSlowsGrowthWithoutStress:
    def test_it_reduces_the_growth_rate(self):
        t = np.linspace(0, 6, 60)
        full = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=1.0)
        limited = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=0.4)

        assert limited["growth_rate"][0] < full["growth_rate"][0]

    def test_it_leaves_promoter_activity_untouched(self):
        """The whole point: a growth knob that is not a stress knob."""
        t = np.linspace(0, 6, 60)
        full = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=1.0)
        limited = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=0.4)

        assert limited["promoter_activity"][0] == pytest.approx(full["promoter_activity"][0])

    def test_it_still_raises_the_reporter_through_dilution_alone(self):
        t = np.linspace(0, 6, 60)
        full = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=1.0)
        limited = simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=0.35)

        assert limited["reporter"][-1] > full["reporter"][-1]

    def test_a_factor_outside_zero_to_one_is_refused(self):
        t = np.linspace(0, 6, 60)
        with pytest.raises(ValueError, match="nutrient_factor"):
            simulate_culture(t, 0.0, PARAMS, 0.05, nutrient_factor=1.4)


class TestAnEmptyAxisIsRefused:
    """A cross product with one empty factor is empty, and an empty design is the one
    failure that looks like success: the sweep writes a valid table with no rows, exits 0,
    and reads as "swept, found nothing" instead of "never swept". Every builder guards every
    axis, so the refusal lands where the mistake was made."""

    def test_a_dose_ladder_with_no_doses(self):
        with pytest.raises(ValueError, match="'doses' is empty"):
            stressor_only_series(PARAMS, ())

    def test_a_decoupling_grid_missing_either_axis(self):
        with pytest.raises(ValueError, match="'doses' is empty"):
            decoupling_grid(PARAMS, (), (1.0, 0.5))
        with pytest.raises(ValueError, match="'nutrient_factors' is empty"):
            decoupling_grid(PARAMS, (0.0, 1.0), ())

    @pytest.mark.parametrize("axis", ["carbon_sources", "growth_phases", "temperatures_c",
                                      "oxygen_fractions", "ph_medium", "ph_cytosolic",
                                      "glucose_g_per_L", "anaerobic_supplements"])
    def test_every_environment_axis_is_guarded_not_just_the_first(self, axis):
        """Guarding only the first factor would pass a one-axis test and still let a later
        empty axis through, because ``product`` short-circuits on the first empty one."""
        from ystwin.generator.design import context_grid

        with pytest.raises(ValueError, match=f"{axis!r} is empty"):
            context_grid(**{axis: ()})

    def test_an_environment_grid_with_no_contexts(self):
        from ystwin.generator.design import environment_grid

        with pytest.raises(ValueError, match="'contexts' is empty"):
            environment_grid(PARAMS, (), (0.0, 1.0))


class TestCollinearity:
    def test_a_stressor_only_series_is_almost_perfectly_collinear(self):
        conditions = stressor_only_series(PARAMS, doses=(0.0, 0.2, 0.5, 1.0, 2.0))

        assert collinearity(conditions) > 0.9

    def test_the_decoupling_grid_breaks_it(self):
        conditions = decoupling_grid(
            PARAMS, doses=(0.0, 0.2, 0.5, 1.0), nutrient_factors=(1.0, 0.6, 0.35)
        )

        assert collinearity(conditions) < 0.6

    def test_the_grid_covers_all_four_quadrants(self):
        conditions = decoupling_grid(
            PARAMS, doses=(0.0, 0.2, 1.0), nutrient_factors=(1.0, 0.4)
        )
        stress = np.array([c.promoter_activity for c in conditions])
        growth = np.array([c.growth_rate for c in conditions])
        mid_s, mid_g = np.median(stress), np.median(growth)

        quadrants = {
            (bool(s > mid_s), bool(g > mid_g)) for s, g in zip(stress, growth)
        }
        assert len(quadrants) == 4

    def test_every_condition_records_what_produced_it(self):
        conditions = decoupling_grid(PARAMS, doses=(0.0, 1.0), nutrient_factors=(1.0, 0.5))

        for c in conditions:
            assert isinstance(c, Condition)
            assert c.nutrient_factor in (1.0, 0.5)
            assert c.dose_mM in (0.0, 1.0)

    def test_a_single_condition_has_no_collinearity_to_measure(self):
        with pytest.raises(ValueError, match="at least"):
            collinearity(decoupling_grid(PARAMS, doses=(0.0,), nutrient_factors=(1.0,)))


class TestTheDesignIsWhatMakesRecoveryPossible:
    """The claim the argument rests on, tested where truth is known exactly.

    At steady state R = k_synth/mu, so log R = log k_synth - log mu and the true
    coefficients are exactly (+1, -1). Collinearity does not bias that fit on clean
    data -- it destabilises it, because the two regressors carry the same information.

    Run at fixed growth rather than in batch. A batch culture cannot reach steady state
    at the slow end: at 0.35 nutrients and 2 mM, mu falls to 0.028/h and the reporter
    needs 35 h to relax, while saturation changes mu long before that. Fixing growth is
    what a chemostat does, and is why one is worth having for this experiment.
    """

    FIXED = CultureParameters(
        mu_max=0.40, growth_ic50=0.8, growth_hill=1.5,
        promoter_basal=1.0e-3, promoter_peak=2.5e-3, promoter_ec50=0.7,
        promoter_hill=2.0, carrying_capacity=50.0,
    )

    def _fit(self, conditions, seed):
        from ystwin.reporter import ReporterKinetics, simulate_reporter

        rng = np.random.default_rng(seed)
        t = np.linspace(0, 60, 300)
        rows = []
        for c in conditions:
            mu = np.full_like(t, c.growth_rate)
            reporter = simulate_reporter(
                t, np.full_like(t, c.promoter_activity), mu,
                ReporterKinetics(k_deg=0.0), r0=c.promoter_activity / c.growth_rate,
            )
            rows.append((c.growth_rate, c.promoter_activity,
                         reporter[-1] * (1 + rng.normal(0, 0.02))))
        growth, stress, signal = (np.array(x) for x in zip(*rows))
        design = np.column_stack([np.ones_like(growth), np.log(stress), np.log(growth)])
        coefficients, *_ = np.linalg.lstsq(design, np.log(signal), rcond=None)
        return coefficients[1], coefficients[2]

    def _spread(self, conditions, n=25):
        fits = np.array([self._fit(conditions, seed) for seed in range(n)])
        return fits.mean(axis=0), fits.std(axis=0)

    def _collinear(self):
        return stressor_only_series(self.FIXED, doses=(0.05, 0.1, 0.2, 0.5, 1.0, 2.0))

    def _decoupled(self):
        return decoupling_grid(
            self.FIXED, doses=(0.05, 0.1, 0.2, 0.5, 1.0, 2.0),
            nutrient_factors=(1.0, 0.6, 0.35),
        )

    def test_the_decoupled_design_recovers_the_true_coefficients(self):
        mean, _ = self._spread(self._decoupled())

        assert mean[0] == pytest.approx(1.0, abs=0.1)
        assert mean[1] == pytest.approx(-1.0, abs=0.1)

    def test_collinearity_costs_stability_not_accuracy(self):
        _, collinear_sd = self._spread(self._collinear())
        _, decoupled_sd = self._spread(self._decoupled())

        assert collinear_sd[0] > decoupled_sd[0] * 3

    def test_adding_nutrient_levels_keeps_improving_the_conditioning(self):
        one = decoupling_grid(self.FIXED, doses=(0.05, 0.2, 1.0), nutrient_factors=(1.0,))
        three = decoupling_grid(
            self.FIXED, doses=(0.05, 0.2, 1.0), nutrient_factors=(1.0, 0.6, 0.35)
        )

        assert collinearity(three) < collinearity(one)

    def test_a_batch_run_cannot_reach_steady_state_at_the_slow_end(self):
        """Why the fixed-growth design above is not a convenience."""
        t = np.linspace(0, 12, 120)
        slow = simulate_culture(t, 2.0, self.FIXED, 0.02, nutrient_factor=0.35)
        expected = self.FIXED.promoter_activity_at(2.0) / slow["growth_rate"][-1]

        assert slow["reporter"][-1] < 0.6 * expected
