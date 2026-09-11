"""Full-family multiplicity, explicit refusals, and method-specific resolution.

Current results are recomputed from the sensor readings without writing tracked artifacts.
The per-dose corrections use exact sign-test p-values, not empirical bootstrap tails.
Sign tests concern median usable-plate folds and have their own discrete p-value floor;
that is not a universal resolution bound for bootstrap or permutation inference.

Two-stage plate/well intervals are not subject to the cluster-only atom policy, although
finite resampling and the shared estimator's refusal rules still matter. Nominal coverage
is approximate, and an unavailable interval is not evidence of equivalence. Sensor-level
permutation tests form a separate family and use actual Holm and BH adjustments.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.analysis.multiplicity import (
    adjusted_interval,
    benjamini_hochberg,
    dose_response_permutation,
    finest_resolvable_alpha,
    holm,
    sign_test_floor,
)


@pytest.fixture(scope="module")
def multiplicity_table(readings):
    """Recompute in memory; historical output snapshots do not validate the current runner."""
    script = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fold_multiplicity.py"
    spec = importlib.util.spec_from_file_location("multiplicity_current_analysis", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.analyse(readings, n_resamples=1000)


@pytest.fixture(scope="module")
def per_dose(multiplicity_table):
    return multiplicity_table[multiplicity_table.measurement == "per_dose_fold"]


@pytest.fixture(scope="module")
def readings():
    return pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")


class TestExactSignResolutionIsMethodSpecific:
    """The chosen sign test has a discrete floor, not every test on the same design."""

    def test_three_plates_cannot_reach_five_percent_by_an_exact_sign_test(self):
        assert sign_test_floor(3) == pytest.approx(0.25)

    def test_nor_can_four(self):
        assert sign_test_floor(4) == pytest.approx(0.125)

    def test_six_plates_is_the_first_count_that_can(self):
        """Doubling the panel is not enough. 2*(1/2)**5 = 0.0625, still above 0.05."""
        assert sign_test_floor(6) < 0.05
        assert sign_test_floor(5) > 0.05

    def test_and_ten_is_the_first_that_survives_correcting_for_twenty(self):
        assert sign_test_floor(10) < 0.05 / 20
        assert sign_test_floor(9) > 0.05 / 20

    def test_the_floor_is_a_probability(self):
        for n in range(1, 12):
            assert 0.0 < sign_test_floor(n) <= 1.0

    def test_it_refuses_a_design_with_no_clusters(self):
        with pytest.raises(ValueError, match="at least one cluster"):
            sign_test_floor(0)


class TestTheBootstrapTailIsNotAPValue:
    """The confirmatory-p columns must not turn an empirical bootstrap tail into a test."""

    def test_adjustments_use_exact_sign_tests_instead_of_bootstrap_tails(self, per_dose):
        family = int(per_dose.family_size.iloc[0])

        assert "p_bootstrap" not in per_dose
        assert per_dose.p_value_method.str.contains("exact sign test").all()
        np.testing.assert_allclose(per_dose.p_holm,
                                   holm(list(per_dose.p_exact_sign), family_size=family),
                                   equal_nan=True)
        np.testing.assert_allclose(per_dose.p_bh,
                                   benjamini_hochberg(list(per_dose.p_exact_sign), family_size=family),
                                   equal_nan=True)

    def test_agreeing_non_tied_plates_reach_their_own_sign_test_floor(self, per_dose):
        agreeing = per_dose[(per_dose.n_plates_sign > 0)
                            & ((per_dose.plates_above_one == per_dose.n_plates_sign)
                               | (per_dose.plates_below_one == per_dose.n_plates_sign))]

        assert len(agreeing) > 0
        expected = agreeing.n_plates_sign.map(lambda n: sign_test_floor(int(n)))
        np.testing.assert_allclose(agreeing.p_exact_sign, expected)
        np.testing.assert_allclose(agreeing.exact_sign_floor, expected)

    def test_no_estimable_sign_p_is_below_its_method_specific_floor(self, per_dose):
        signed = per_dose[per_dose.n_plates_sign > 0]
        expected = signed.n_plates_sign.map(lambda n: sign_test_floor(int(n)))

        assert (signed.p_exact_sign >= expected).all()
        assert per_dose.loc[per_dose.n_plates_contributing == 0, "p_exact_sign"].isna().all()


class TestFamilyIntervalsRespectTheResamplingMethod:
    """Use the planned divisor, contributors and actual resampling method for each row."""

    def test_unadjusted_calls_require_an_estimable_interval_excluding_unity(self, per_dose):
        expected = per_dose.estimable & ((per_dose.low > 1) | (per_dose.high < 1))

        assert per_dose.clears.equals(expected)

    def test_the_family_is_the_folds_asked_about_not_the_ones_estimable(self, per_dose, readings):
        family = readings[readings.dose_mM > 0].groupby(["construct", "dose_mM"]).ngroups

        assert set(per_dose.family_size) == {family}
        assert len(per_dose) == family
        np.testing.assert_allclose(per_dose.family_alpha, 0.05 / family)

    def test_two_stage_family_intervals_are_computable_without_claiming_exact_coverage(self, per_dose):
        computed = per_dose[per_dose.family_wise_resolvable]

        assert not computed.empty
        assert set(computed.resampling) == {"two_stage"}
        assert np.isfinite(computed[["fw_low", "fw_high"]]).all().all()
        assert (computed.fw_n_resamples * computed.family_alpha / 2 >= 1).all()
        assert computed.coverage_status.str.contains("approximate").all()

    def test_refused_family_bounds_are_empty_and_carry_the_reason(self, per_dose):
        refused = per_dose[~per_dose.family_wise_resolvable]

        assert refused.fw_low.isna().all()
        assert refused.fw_high.isna().all()
        assert not refused.clears_family_wise.any()
        assert refused.family_wise_method.str.contains("refused").all()

    def test_the_interval_api_retains_its_cluster_only_tail_policy(self):
        draws = np.log(np.random.default_rng(0).lognormal(0.3, 0.2, 4000))

        with pytest.raises(ValueError, match="finer than"):
            adjusted_interval(draws, 1.35, 4, 0.05 / 24, resampling="cluster_only")

    def test_cluster_only_policy_first_allows_this_level_at_five_plates(self):
        assert finest_resolvable_alpha(4) > 0.05 / 24
        assert finest_resolvable_alpha(5) < 0.05 / 24

    def test_because_the_t_quantile_explodes_at_two_degrees_of_freedom(self):
        """The t/normal inflation ratio grows at a stricter simultaneous level.
        This arithmetic is not a coverage guarantee or an impossibility result."""
        from scipy import stats

        nominal = stats.t.ppf(0.975, df=2) / stats.norm.ppf(0.975)
        family = stats.t.ppf(1 - 0.05 / 40, df=2) / stats.norm.ppf(1 - 0.05 / 40)

        assert family > 3 * nominal

    def test_and_costs_a_third_less_at_three_degrees_of_freedom(self):
        """The corrected-to-nominal t-quantile ratio falls as degrees of freedom rise.
        This isolates the inflation arithmetic, not the empirical interval's coverage.
        """
        from scipy import stats

        def price(df):
            return stats.t.ppf(1 - 0.05 / 40, df) / stats.t.ppf(0.975, df)

        assert price(2) == pytest.approx(4.64, abs=0.01)
        assert price(3) == pytest.approx(2.97, abs=0.01)
        assert price(3) < 0.7 * price(2)

    def test_a_corrected_interval_still_contains_the_nominal_one(self):
        """For centred draws, a stricter level widens this approximate interval.
        The example stays within the explicitly chosen cluster-only tail policy.
        """
        draws = np.log(np.random.default_rng(0).lognormal(0.3, 0.2, 4000))

        low, high = adjusted_interval(draws, 1.35, 4, 0.05)
        wider_low, wider_high = adjusted_interval(draws, 1.35, 4, 0.05 / 6)

        assert wider_low < low
        assert wider_high > high


class TestTheCorrectionsThemselves:
    def test_holm_is_monotone_and_never_shrinks_a_p_value(self):
        raw = [0.001, 0.01, 0.04, 0.3]
        adjusted = holm(raw)

        assert all(a >= r for a, r in zip(adjusted, raw))
        assert adjusted == sorted(adjusted)

    def test_holm_multiplies_the_smallest_by_the_family_size(self):
        assert holm([0.001, 0.5, 0.5, 0.5])[0] == pytest.approx(0.004)

    def test_benjamini_hochberg_is_never_more_conservative_than_holm(self):
        raw = [0.001, 0.008, 0.02, 0.04, 0.3]

        assert all(bh <= h + 1e-12 for bh, h in zip(benjamini_hochberg(raw), holm(raw)))

    def test_both_cap_at_one(self):
        assert max(holm([0.9, 0.9, 0.9])) <= 1.0
        assert max(benjamini_hochberg([0.9, 0.9, 0.9])) <= 1.0

    def test_a_single_test_is_unchanged_by_either(self):
        assert holm([0.03]) == pytest.approx([0.03])
        assert benjamini_hochberg([0.03]) == pytest.approx([0.03])

    def test_an_interval_refuses_a_level_outside_zero_to_one(self):
        with pytest.raises(ValueError, match="alpha must be in"):
            adjusted_interval(np.log(np.array([0.9, 1.0, 1.1])), 1.0, 3, 1.5)


class TestTheQuestionThatSurvives:
    """One test per sensor, not six. The plate effect is held by permuting inside it."""

    @pytest.mark.parametrize("construct",
                             ["UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"])
    def test_every_sensor_shows_a_dose_response(self, multiplicity_table, construct):
        row = multiplicity_table[
            (multiplicity_table.measurement == "dose_response_permutation")
            & (multiplicity_table.construct == construct)]

        assert len(row) == 1
        assert float(row.p_permutation.iloc[0]) < 0.001

    def test_oxidative_trends_are_not_per_dose_effect_size_claims(self, multiplicity_table,
                                                                per_dose):
        """Dose-response and median per-dose fold tests ask distinct questions."""
        oxidative = per_dose[per_dose.construct.isin(["NativeYap1", "AlteredYap1"])]
        trends = multiplicity_table[
            (multiplicity_table.measurement == "dose_response_permutation")
            & (multiplicity_table.construct.isin(["NativeYap1", "AlteredYap1"]))]

        assert (trends.p_holm < 0.05).all()
        assert not (oxidative.p_holm <= 0.05).any()
        assert trends.p_value_method.str.contains("within-plate").all()
        assert oxidative.target.str.contains("geometric mean").all()

    def test_and_all_four_survive_correcting_for_all_four(self, multiplicity_table):
        trends = multiplicity_table[
            multiplicity_table.measurement == "dose_response_permutation"]

        assert len(trends) == 4
        assert (trends.p_holm < 0.05).all()

    def test_the_statistic_is_positive_for_every_sensor(self, multiplicity_table):
        trends = multiplicity_table[
            multiplicity_table.measurement == "dose_response_permutation"]

        assert (trends.statistic > 0).all()

    def test_a_permutation_p_is_never_exactly_zero(self, readings):
        """The observed assignment is itself one of the permutations, so no finite number of
        draws can report that it is unreachable under the null."""
        got = dose_response_permutation(readings, "UPRE1", max_dose=1.0,
                                        n_permutations=200)

        assert got.p_value >= 1 / 201

    def test_shuffled_activity_gives_no_trend(self, readings):
        """The control. Break the dose-activity association and the test must stop firing."""
        rng = np.random.default_rng(0)
        scrambled = readings.copy()
        scrambled["activity_late"] = rng.permutation(scrambled.activity_late.to_numpy())

        got = dose_response_permutation(scrambled, "UPRE1", max_dose=1.0,
                                        n_permutations=2000)

        assert got.p_value > 0.05

    def test_it_refuses_a_construct_with_one_dose_left_after_the_cap(self, readings):
        with pytest.raises(ValueError, match="a trend needs at least two"):
            dose_response_permutation(readings, "UPRE1", max_dose=0.0)

    def test_it_refuses_a_construct_that_is_not_there(self, readings):
        with pytest.raises(ValueError, match="no rows for construct"):
            dose_response_permutation(readings, "NotASensor")


class TestThePermutationTestRefusesRatherThanReportingSignificance:
    """Its failure mode was `p = 1/(n+1)`, the most significant answer it can give.

    `statistic` returns NaN when no plate carries more than one dose, and
    `np.abs(null) >= np.nan` is False for every draw, so the exceedance count was zero and
    the smallest reachable p came back for a question that had not been asked. Found by the
    2026-08-30 audit's `new-code` lens.
    """

    def _one_dose_per_plate(self):
        return pd.DataFrame({
            "plate": ["p1", "p1", "p2", "p2"],
            "construct": ["UPRE1"] * 4,
            "stressor": ["DTT"] * 4,
            "dose_mM": [0.0, 0.0, 0.5, 0.5],
            "activity_late": [1.0, 1.1, 1.4, 1.5],
        })

    def test_it_refuses_when_no_plate_carries_a_ladder(self):
        with pytest.raises(ValueError, match="undefined"):
            dose_response_permutation(self._one_dose_per_plate(), "UPRE1")

    def test_the_refusal_says_why_rather_than_naming_a_number(self):
        with pytest.raises(ValueError) as raised:
            dose_response_permutation(self._one_dose_per_plate(), "UPRE1")

        assert "more than one distinct dose" in str(raised.value)

    def test_it_does_not_return_the_most_significant_p_it_can(self):
        """The regression, stated as the number it used to give."""
        try:
            got = dose_response_permutation(self._one_dose_per_plate(), "UPRE1",
                                            n_permutations=200)
        except ValueError:
            return
        assert got.p_value != pytest.approx(1 / 201)
        pytest.fail("it returned a p-value for an undefined statistic")

    def test_a_real_within_plate_ladder_still_works(self):
        """The fix must not refuse the case the function exists for."""
        frame = pd.DataFrame({
            "plate": ["p1"] * 4 + ["p2"] * 4,
            "construct": ["UPRE1"] * 8,
            "stressor": ["DTT"] * 8,
            "dose_mM": [0.0, 0.2, 0.5, 1.0] * 2,
            "activity_late": [1.0, 1.2, 1.4, 1.7, 1.1, 1.25, 1.5, 1.8],
        })

        got = dose_response_permutation(frame, "UPRE1", n_permutations=500)

        assert np.isfinite(got.statistic) and got.statistic > 0
        assert 0 < got.p_value <= 1


def test_holm_and_bh_have_distinct_correct_arithmetic():
    values = [0.01, 0.04, 0.03, 0.002]
    assert holm(values) == pytest.approx([0.03, 0.06, 0.06, 0.008])
    assert benjamini_hochberg(values) == pytest.approx([0.02, 0.04, 0.04, 0.008])


def test_unestimated_tests_remain_in_a_prespecified_family():
    assert holm([0.001, 0.01], family_size=24) == pytest.approx([0.024, 0.23])
    assert benjamini_hochberg([0.001, 0.01], family_size=24) == pytest.approx([0.024, 0.12])
    with pytest.raises(ValueError, match="family_size"):
        holm([0.1, 0.2], family_size=1)


@pytest.mark.parametrize("values", [[-0.01], [1.01], [np.inf]])
def test_invalid_probabilities_are_not_silently_corrected(values):
    for correct in (holm, benjamini_hochberg):
        with pytest.raises(ValueError, match="probabilities"):
            correct(values)


def test_sign_test_resolution_does_not_bound_parametric_inference():
    from scipy import stats
    from ystwin.analysis.multiplicity import MultipleTests

    p = stats.ttest_1samp([9.9, 10.0, 10.1], 0.0).pvalue
    assert p < 0.001 < sign_test_floor(3)
    result = MultipleTests(1, 3, 0, 0, sign_test_floor(3), pd.DataFrame())
    assert "exact sign test" in result.summary()
    assert "not a resolution bound" in result.summary()


def test_two_stage_bootstrap_is_not_subject_to_cluster_only_atom_policy():
    draws = np.random.default_rng(1).normal(0.1, 0.1, 10000)
    low, high = adjusted_interval(draws, np.exp(0.1), 4, 0.05 / 24, resampling="two_stage")
    assert 0 < low < np.exp(0.1) < high
    assert finest_resolvable_alpha(1) <= 1


def test_constant_response_does_not_become_a_significant_permutation_result():
    frame = pd.DataFrame({"plate": ["a"] * 4 + ["b"] * 4, "construct": "R",
                          "dose_mM": [0.0, 0.2, 0.5, 1.0] * 2, "activity_late": 1.0})
    result = dose_response_permutation(frame, "R", n_permutations=20)
    assert result.p_value == 1.0
    assert result.n_informative_plates == 2
    assert "exchangeable" in result.exchangeability
    assert "permutation" in result.method
