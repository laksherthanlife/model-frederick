"""The window rule is measured where it runs, and the numbers justifying it are produced.

`reporter.py::default_activity_window_h` sets the smoothing window every dilution-corrected
number in this repository depends on. It carried a five-cell accuracy table in its own
docstring that **no script produced and no test read**, and the only test of the rule itself
checked it at 12 h and 48 h -- six and twelve times the density of the committed plates,
which are 25 points over 4.00 h. The rule was tested only where it works.

These tests close both halves. They are slow enough to keep the replicate counts small:
the point is that the shape of the answer holds, and `scripts/estimator_accuracy.py` is what
produces the numbers at full precision.
"""

from __future__ import annotations

import numpy as np
import pytest

from ystwin.analysis.estimator_accuracy import (
    DOCSTRING_REGIME,
    MATURATION_HALF_LIFE_H,
    PLATE,
    PLATE_DURATION_H,
    PLATE_POINTS,
    Geometry,
    NisArm,
    accuracy_curve,
    fold_recovery,
    maturation_bias,
    measurement_inflation_for_nis,
    nis_error_decomposition,
    recovery_error,
    true_activity,
)
from ystwin.reporter import TARGET_ACTIVITY_WINDOW_H, default_activity_window_h
from ystwin.readings import CorrectedOD, CorrectedRFU

REPLICATES = 40


class TestTheGeometryIsThePlatesOwn:
    def test_the_committed_exports_are_twenty_five_points_over_four_hours(self, plate_text):
        """Read off the committed text rather than taken from this module's constants,
        so the constants cannot drift away from the plates they claim to describe."""
        times = plate_text
        assert len(times) == PLATE_POINTS
        assert float(times[-1] - times[0]) == pytest.approx(PLATE_DURATION_H, abs=0.01)

    def test_the_rule_picks_a_third_of_the_run_there(self):
        assert PLATE.auto_window_h == pytest.approx(PLATE_DURATION_H / 3.0)

    def test_which_is_wider_than_the_sixth_it_used_to_pick(self):
        """The whole change, in one assertion. A sixth of a four-hour run is 0.667 h."""
        assert PLATE.auto_window_h > PLATE_DURATION_H / 6.0


class TestTheErrorIsSetByDurationNotByPointCount:
    """The finding that makes the old rule wrong rather than merely unlucky.

    A window's accuracy is a property of how many HOURS it spans, not of how many samples
    fall in it and not of what fraction of the run it is. So a rule expressed as a fraction
    of the run is correct only at the run length it was tuned on -- which is exactly how a
    table measured over 24 h came to justify a window used over 4 h.
    """

    def test_two_geometries_agree_when_their_windows_match_in_hours(self):
        dense_long = Geometry("dense long", 145, 24.0)
        sparse_short = Geometry("sparse short", 25, 4.0)

        a, _ = recovery_error(dense_long, 1.17, n_replicates=REPLICATES)
        b, _ = recovery_error(sparse_short, 1.17, n_replicates=REPLICATES)

        assert a == pytest.approx(b, abs=0.03)

    def test_and_disagree_when_their_windows_match_in_points(self):
        """Three points spans 20 minutes on one geometry and 3 hours on the other. If point
        count governed the error these would agree; they differ by more than fourfold."""
        dense_long = Geometry("dense long", 145, 24.0)
        sparse_long = Geometry("sparse long", 25, 24.0)

        a, _ = recovery_error(dense_long, 3 * dense_long.dt_h, n_replicates=REPLICATES)
        b, _ = recovery_error(sparse_long, 3 * sparse_long.dt_h, n_replicates=REPLICATES)

        assert a > 3 * b

    def test_a_wider_window_is_more_accurate_at_the_plate_geometry(self):
        narrow, _ = recovery_error(PLATE, PLATE_DURATION_H / 6.0, n_replicates=REPLICATES)
        chosen, _ = recovery_error(PLATE, PLATE.auto_window_h, n_replicates=REPLICATES)

        assert chosen < narrow


class TestWhatTheWindowBuysIsPrecisionNotAccuracy:
    """The reason the change did not move the point estimates, stated as a test.

    Every committed result is a FOLD -- a ratio of a dosed well to its own plate's control
    -- and a smoother treats numerator and denominator alike, so whatever it does to one it
    largely does to the other. That is why widening the window moved twenty-one folds by at
    most a few thousandths while changing which of them clear 1.0.
    """

    @pytest.mark.parametrize("window_h", [0.667, 1.333, 2.0])
    def test_the_fold_is_accurate_at_every_window(self, window_h):
        got = fold_recovery(PLATE, window_h=window_h, n_replicates=REPLICATES)

        assert abs(got["relative_error"]) < 0.02

    def test_but_a_wider_window_gives_a_tighter_spread(self):
        narrow = fold_recovery(PLATE, window_h=0.667, n_replicates=REPLICATES)
        wide = fold_recovery(PLATE, window_h=2.0, n_replicates=REPLICATES)

        assert wide["spread_of_recovered_fold"] < narrow["spread_of_recovered_fold"]

    def test_the_naive_fold_overstates_it_because_the_dose_slows_growth(self):
        """The confound the whole package exists for, reproduced end to end: a dosed well
        keeping 75% of the control's growth reads high on RFU/OD even with the induction
        exactly as simulated."""
        got = fold_recovery(PLATE, n_replicates=REPLICATES)

        assert got["naive_fold"] > got["target_fold"] * 1.1


class TestIgnoringMaturationIsPricedRatherThanAssumed:
    """`promoter_activity_from_total` refuses maturation outright, so the pipeline runs a
    non-maturing model against a reporter that matures. That was never costed."""

    def test_it_biases_a_single_wells_activity_by_about_ten_percent(self):
        got = maturation_bias(PLATE, n_replicates=REPLICATES)

        assert got["added_by_maturation"] > 0.05

    def test_and_a_fold_by_far_less_because_a_ratio_cancels_it(self):
        clean = fold_recovery(PLATE, n_replicates=REPLICATES)
        maturing = fold_recovery(PLATE, half_life_h=MATURATION_HALF_LIFE_H,
                                 n_replicates=REPLICATES)
        single_well = maturation_bias(PLATE, n_replicates=REPLICATES)

        on_fold = abs(maturing["relative_error"] - clean["relative_error"])

        assert on_fold < single_well["added_by_maturation"] / 2.0

    def test_the_pipeline_cannot_correct_for_it_even_if_it_wanted_to(self):
        """Recorded as a test because it is the reason the bias is a bias and not an error:
        there is no code path by which a maturation term could reach a committed number."""
        from ystwin.reporter import ReporterKinetics, promoter_activity_from_total

        t = PLATE.times()
        with pytest.raises(ValueError, match="maturation is not recoverable"):
            promoter_activity_from_total(t, CorrectedRFU(np.ones_like(t)), CorrectedOD(np.ones_like(t)),
                                         ReporterKinetics(k_deg=0.0, k_mat=2.0))


class TestTheSimulatedTruthIsShapedLikeAnInduction:
    def test_it_starts_at_basal_and_rises(self):
        t = PLATE.times()
        activity = true_activity(t)

        assert activity[-1] > activity[0]

    def test_a_flat_truth_is_recovered_at_any_window(self):
        """The control on the whole method: if the estimator could not recover a constant
        the error table would be measuring the smoother rather than the transient."""
        flat, _ = recovery_error(PLATE, 0.667, n_replicates=REPLICATES)
        curve = accuracy_curve(PLATE, windows_h=(0.667,), n_replicates=REPLICATES)

        assert curve[0].median_relative_error == pytest.approx(flat, abs=1e-9)


class TestTheAccuracyCurveIsAProducedTableNotADocstring:
    def test_every_window_reported_fits_inside_the_run(self):
        for point in accuracy_curve(PLATE, n_replicates=8):
            assert point.window_h <= PLATE.duration_h

    def test_the_auto_window_is_always_in_the_curve_and_flagged(self):
        flagged = [p for p in accuracy_curve(PLATE, n_replicates=8) if p.is_auto_window]

        assert len(flagged) == 1
        assert flagged[0].window_h == pytest.approx(PLATE.auto_window_h)

    def test_the_long_regime_still_picks_a_sixth_of_its_run(self):
        """The constant only decides short runs. A 24 h run's sixth is 4 h, which already
        outranks it, so nothing about the older geometry moved."""
        assert DOCSTRING_REGIME.auto_window_h == pytest.approx(4.0)
        assert DOCSTRING_REGIME.auto_window_h > TARGET_ACTIVITY_WINDOW_H


class TestPredictionCalibrationMetrics:
    def test_coverage_and_nis_are_for_predictions_not_refitted_residuals(self):
        from ystwin.analysis.estimator_accuracy import prediction_metrics

        result = prediction_metrics([0.0, 1.0, 3.0, np.nan], np.zeros(4), np.ones(4))
        assert result["n_predictions"] == 3
        assert result["n_excluded"] == 1
        assert result["coverage"] == pytest.approx(2 / 3)
        assert result["mean_nis"] == pytest.approx(10 / 3)
        assert result["nominal_coverage"] == 0.95

    def test_autocorrelation_never_bridges_a_missing_timepoint(self):
        from ystwin.analysis.estimator_accuracy import prediction_metrics

        result = prediction_metrics([1.0, -1.0, np.nan, 1.0, -1.0], np.zeros(5), np.ones(5))
        assert result["lag1_autocorr"] == pytest.approx(-0.5)
        assert result["lag1_pairs"] == 2

    def test_omitted_rows_can_be_identified_by_their_original_step_indices(self):
        from ystwin.analysis.estimator_accuracy import prediction_metrics

        result = prediction_metrics([1.0, -1.0, 1.0, -1.0], np.zeros(4), np.ones(4),
                                    step_indices=[0, 1, 3, 4])
        assert result["lag1_autocorr"] == pytest.approx(-0.5)
        assert result["lag1_pairs"] == 2

    def test_a_constant_sequence_has_no_estimable_whiteness(self):
        from ystwin.analysis.estimator_accuracy import prediction_metrics

        result = prediction_metrics(np.ones(5), np.zeros(5), np.ones(5))
        assert np.isnan(result["lag1_autocorr"])

    def test_invalid_shapes_and_coverage_are_refused(self):
        from ystwin.analysis.estimator_accuracy import prediction_metrics

        with pytest.raises(ValueError, match="shape"):
            prediction_metrics([1.0, 2.0], [1.0], [1.0])
        with pytest.raises(ValueError, match="coverage"):
            prediction_metrics([1.0], [1.0], [1.0], nominal_coverage=1.0)


class TestCustomLossIsNotSilentlyReplaced:
    def test_maturation_bias_preserves_the_supplied_degradation(self, monkeypatch):
        from ystwin.analysis import estimator_accuracy as analysis
        from ystwin.reporter import ReporterKinetics

        calls = []

        def record(*args, **kwargs):
            calls.append(kwargs)
            return 0.0, 0.0

        monkeypatch.setattr(analysis, "recovery_error", record)
        kinetics = ReporterKinetics(k_deg=0.7)
        analysis.maturation_bias(PLATE, kinetics=kinetics)
        assert calls[1]["generate_with"].k_deg == kinetics.k_deg
        assert calls[1]["generate_with"].immature_loss == kinetics.immature_loss


class TestTheOlderGeometriesAreUnchanged:
    """The change is minimal or it is not defensible. These are the two cases the previous
    test suite asserted, and both must still hold."""

    def test_twelve_hours_at_ten_minute_sampling_still_gives_two_hours(self):
        assert default_activity_window_h(np.linspace(0, 12, 73)) == pytest.approx(2.0, rel=0.01)

    def test_forty_eight_hours_still_gives_eight(self):
        assert default_activity_window_h(np.linspace(0, 48, 289)) == pytest.approx(8.0, rel=0.01)

    def test_and_the_sampling_floor_still_outranks_everything(self):
        """Hourly sampling over six hours: four intervals is 4 h, which exceeds both the
        target and the third-of-the-run cap, and must still win."""
        assert default_activity_window_h(np.linspace(0, 6, 7)) >= 4.0


class TestTheMeasurementInflationANisWouldDemand:
    """`measurement_inflation_for_nis` is the algebra the measurement arms then test.

    An inflated NIS says the squared innovation is large against the declared predictive
    variance. Blaming the measurement model means claiming the declared relative sigma is
    too small; this asks by how much, exactly, and the answer is the cheapest possible
    measurement-only story. If the empirical arm does not reproduce it, the assumption it
    holds -- that inflating sigma leaves the prediction and the state variance alone -- is
    what failed, and that is itself a finding.
    """

    def test_a_calibrated_filter_needs_no_inflation(self):
        assert measurement_inflation_for_nis(1.0, 0.5) == pytest.approx(1.0)

    def test_an_over_dispersed_filter_is_told_to_shrink_rather_than_inflate(self):
        """NIS below one is a filter claiming more uncertainty than it has. The
        multiplier is below one, which is a real instruction and not a refusal."""
        assert 0.0 < measurement_inflation_for_nis(0.25, 1.0) < 1.0

    def test_when_measurement_is_the_whole_variance_the_answer_is_the_square_root(self):
        """The floor case, and the one worth quoting: NIS 45 buys a factor of 6.7 at
        best, so a 2% declared sigma would have to be 13.5%."""
        assert measurement_inflation_for_nis(45.3, 1.0) == pytest.approx(np.sqrt(45.3))

    def test_a_smaller_measurement_share_demands_a_larger_multiplier(self):
        """Measurement variance is what the multiplier scales. The less of the predictive
        variance it is, the more it has to be scaled to cover the same innovation, so the
        full-share answer is a lower bound on every measurement-only account."""
        shares = [1.0, 0.5, 0.2, 0.05]
        needed = [measurement_inflation_for_nis(45.3, share) for share in shares]
        assert needed == sorted(needed)
        assert all(value >= np.sqrt(45.3) for value in needed)

    def test_no_declared_measurement_variance_makes_the_story_unavailable(self):
        """Scaling zero is zero. An excess cannot be bought at any price here, and the
        infinity says so rather than a finite number implying it could be."""
        assert measurement_inflation_for_nis(45.3, 0.0) == float("inf")
        assert measurement_inflation_for_nis(0.5, 0.0) == 0.0

    def test_an_unestimable_input_stays_unestimable(self):
        assert np.isnan(measurement_inflation_for_nis(float("nan"), 0.5))
        assert np.isnan(measurement_inflation_for_nis(45.3, float("nan")))

    def test_a_share_that_is_not_a_share_is_refused(self):
        with pytest.raises(ValueError):
            measurement_inflation_for_nis(45.3, 1.4)
        with pytest.raises(ValueError):
            measurement_inflation_for_nis(45.3, 0.5, target=0.0)


def _arm(label, family, nis, *, particles=500, multiplier=1.0, share=0.5, min_ess=1.0):
    return NisArm(label=label, family=family, n_particles=particles,
                  sigma_multiplier=multiplier, mean_nis=nis, measurement_variance_share=share,
                  min_ess=min_ess, min_unique_ancestors=1.0, coverage=0.23, well_channels=319)


class TestTheNisErrorDecomposition:
    """What each arm removes of the excess, as a share of the excess and not of the NIS.

    The three families are not independently estimable from the NIS alone -- that is why
    process error is a residual here rather than an arm. What the table has to get right
    is the bookkeeping: shares against the baseline excess, unclipped signs, and a refusal
    when the excess it would divide by does not exist.
    """

    def test_an_arm_that_reaches_the_target_removes_the_whole_excess(self):
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                        _arm("sigma_x7", "measurement", 1.0, multiplier=7.0)])
        assert rows[1]["excess_removed"] == pytest.approx(1.0)

    def test_an_arm_that_changes_nothing_removes_nothing(self):
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                        _arm("particles=8000", "particles", 45.0, particles=8000)])
        assert rows[1]["excess_removed"] == pytest.approx(0.0)

    def test_an_arm_that_makes_it_worse_reports_a_negative_share(self):
        """Clipping this at zero would hide an arm that hurt, which is the one result
        that would change what the sweep concludes."""
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                        _arm("particles=8000", "particles", 60.0, particles=8000)])
        assert rows[1]["excess_removed"] < 0.0

    def test_the_shares_are_against_the_excess_and_not_against_the_nis(self):
        """A baseline of 45 that falls to 23 has removed half its excess above one, not
        the 49% of its NIS that a ratio of the levels would report."""
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                        _arm("sigma_x2", "measurement", 23.0, multiplier=2.0)])
        assert rows[1]["excess_removed"] == pytest.approx(0.5)

    def test_a_baseline_already_at_the_target_apportions_nothing(self):
        """There is no excess to divide, and a share of a negative excess would invert
        the sign of every arm in the table."""
        rows = nis_error_decomposition([_arm("baseline", "baseline", 0.22),
                                        _arm("particles=8000", "particles", 1.4, particles=8000)])
        assert all(np.isnan(row["excess_removed"]) for row in rows)
        assert all(row["baseline_has_excess"] is False for row in rows)

    def test_the_baseline_row_comes_first_whatever_order_it_arrives_in(self):
        rows = nis_error_decomposition([_arm("sigma_x2", "measurement", 23.0, multiplier=2.0),
                                        _arm("baseline", "baseline", 45.0)])
        assert rows[0]["label"] == "baseline"
        assert [row["label"] for row in rows[1:]] == ["sigma_x2"]

    def test_every_row_carries_the_baseline_it_is_a_share_of(self):
        """A decomposition row read on its own must not be mistakable for an absolute."""
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                        _arm("sigma_x2", "measurement", 23.0, multiplier=2.0)])
        assert {row["baseline_mean_nis"] for row in rows} == {45.0}
        assert {row["target_nis"] for row in rows} == {1.0}

    def test_a_measurement_arm_that_moved_the_ensemble_is_shown_to_have_moved_it(self):
        """The confound that decides whether the elimination means anything. Inflating the
        declared sigma loosens the likelihood, which keeps particles alive; an arm that
        bought its NIS that way has not isolated the measurement model."""
        rows = nis_error_decomposition([
            _arm("baseline", "baseline", 45.3, min_ess=1.029),
            _arm("sigma_x8", "measurement", 0.42, multiplier=8.0, min_ess=101.8)])
        assert rows[0]["min_ess_ratio"] == pytest.approx(1.0)
        assert rows[1]["min_ess_ratio"] == pytest.approx(101.8 / 1.029)

    def test_an_arm_that_left_the_ensemble_alone_reports_a_ratio_of_one(self):
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0, min_ess=1.03),
                                        _arm("sigma_x2", "measurement", 23.0, multiplier=2.0,
                                             min_ess=1.03)])
        assert rows[1]["min_ess_ratio"] == pytest.approx(1.0)

    def test_no_baseline_ensemble_gives_no_ratio_rather_than_an_infinity(self):
        """A ratio against zero would read as an enormous confound where there is simply
        nothing recorded to compare against."""
        rows = nis_error_decomposition([_arm("baseline", "baseline", 45.0, min_ess=0.0),
                                        _arm("sigma_x2", "measurement", 23.0, multiplier=2.0,
                                             min_ess=30.0)])
        assert all(np.isnan(row["min_ess_ratio"]) for row in rows)

    @pytest.mark.parametrize("arms", [
        [],
        [_arm("a", "particles", 4.0)],
        [_arm("a", "baseline", 4.0), _arm("b", "baseline", 4.0)],
    ])
    def test_exactly_one_baseline_or_nothing_to_compare_against(self, arms):
        with pytest.raises(ValueError, match="exactly one baseline"):
            nis_error_decomposition(arms)

    def test_a_repeated_label_is_refused_rather_than_silently_collapsed(self):
        with pytest.raises(ValueError, match="unique"):
            nis_error_decomposition([_arm("same", "baseline", 45.0),
                                     _arm("same", "particles", 40.0)])

    def test_an_unknown_family_is_refused(self):
        """Process error is a residual, not an arm. A row claiming to *be* the process
        arm would assert an identification the design cannot make."""
        with pytest.raises(ValueError, match="unknown decomposition families"):
            nis_error_decomposition([_arm("baseline", "baseline", 45.0),
                                     _arm("process", "process", 40.0)])
