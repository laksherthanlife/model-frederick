"""Does a co-expressed reference reporter rescue the design that was already run?

The dose-only design has 4% power to attribute a 1.5x response to dose rather than to
the growth the dose caused. A co-expressed reference changes the observable from R to
R_stress/R_reference, in which mu cancels exactly -- so the collinearity that made
attribution impossible should stop mattering.

What can still break it: mismatched maturation between the two fluorophores, and error
introduced by spectral unmixing.
"""


from ystwin.analysis.power import ratio_discrimination_power
from ystwin.generator.design import decoupling_grid, stressor_only_series
from ystwin.generator.literature import parameters_from_literature
from ystwin.generator.unmixing import PANEL, FluorophoreSpec

PARAMS = parameters_from_literature("UPRE2")
DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0)
COLLINEAR = stressor_only_series(PARAMS, doses=DOSES)
DECOUPLED = decoupling_grid(PARAMS, doses=DOSES, nutrient_factors=(1.0, 0.6, 0.35))
MATCHED = PANEL["YFP"]


class TestTheReferenceRescuesTheCollinearDesign:
    def test_a_matched_reference_attributes_on_the_dose_only_design(self):
        """The design already run becomes usable."""
        power = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, n_simulations=60, seed=0,
        )

        assert power > 0.8

    def test_that_is_a_large_gain_over_a_single_reporter(self):
        from ystwin.analysis.power import discrimination_power

        single = discrimination_power(PARAMS, COLLINEAR, induction_fold=1.5,
                                      n_replicates=6, n_simulations=60, seed=0)
        paired = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, n_simulations=60, seed=0,
        )

        assert paired > single + 0.4

    def test_it_also_works_on_the_decoupled_design(self):
        power = ratio_discrimination_power(
            PARAMS, DECOUPLED, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, n_simulations=60, seed=0,
        )

        assert power > 0.8

    def test_no_induction_is_still_not_detected(self):
        power = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.0, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, n_simulations=60, seed=0,
        )

        assert power < 0.2


class TestGrowthMustBeLeftOutOfTheModel:
    """The ratio cancels mu structurally, so the term has no true effect.

    Putting it back does not make the analysis more careful; on a collinear design it
    inflates the dose standard error and throws away the power the ratio just bought.
    """

    def test_including_growth_costs_most_of_the_power(self):
        without = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, n_simulations=60, seed=0,
        )
        with_growth = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED, include_growth=True,
            n_simulations=60, seed=0,
        )

        assert without > with_growth + 0.4


class TestTheCancellationIsCheckableOnce:
    def test_a_matched_pair_is_exactly_invariant(self):
        from ystwin.analysis.power import ratio_growth_invariance

        departure = ratio_growth_invariance(2e-3, 1e-3, MATCHED, MATCHED)

        assert departure < 1e-9

    def test_a_mismatched_pair_is_not(self):
        from ystwin.analysis.power import ratio_growth_invariance

        slow = FluorophoreSpec("slow", 549, 565, maturation_h=2.3)
        departure = ratio_growth_invariance(2e-3, 1e-3, slow, MATCHED)

        assert departure > 0.05

    def test_the_departure_grows_with_the_mismatch(self):
        from ystwin.analysis.power import ratio_growth_invariance

        small = FluorophoreSpec("a", 549, 565, maturation_h=0.6)
        large = FluorophoreSpec("b", 549, 565, maturation_h=4.0)

        assert (ratio_growth_invariance(2e-3, 1e-3, large, MATCHED)
                > ratio_growth_invariance(2e-3, 1e-3, small, MATCHED))


class TestWhatBreaksIt:
    def test_mismatched_maturation_breaks_the_invariance_it_relies_on(self):
        from ystwin.analysis.power import ratio_growth_invariance

        slow = FluorophoreSpec("mOrange2-like", 549, 565, maturation_h=2.3)

        assert ratio_growth_invariance(2e-3, 1e-3, MATCHED, MATCHED) < 1e-9
        assert ratio_growth_invariance(2e-3, 1e-3, slow, MATCHED) > 0.05

    def test_unmixing_error_costs_power(self):
        clean = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.2, n_replicates=4,
            stress_spec=MATCHED, reference_spec=MATCHED,
            unmixing_cv=0.0, n_simulations=60, seed=2,
        )
        noisy = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.2, n_replicates=4,
            stress_spec=MATCHED, reference_spec=MATCHED,
            unmixing_cv=0.15, n_simulations=60, seed=2,
        )

        assert clean >= noisy

    def test_a_reference_that_is_itself_stress_responsive_destroys_it(self):
        """The reference must be constitutive, or the ratio divides out the signal."""
        power = ratio_discrimination_power(
            PARAMS, COLLINEAR, induction_fold=1.5, n_replicates=6,
            stress_spec=MATCHED, reference_spec=MATCHED,
            reference_responds=True, n_simulations=60, seed=0,
        )

        assert power < 0.3
