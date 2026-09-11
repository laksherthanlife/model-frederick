"""Power: how many biological replicates to detect a real induction after correction.

The question a design has to answer before it is run. Detection means the recovered
promoter activity varies with dose *after* the dilution correction -- not that the raw
signal does, which growth inhibition alone guarantees.
"""

import numpy as np
import pytest

from ystwin.analysis.power import (
    detection_power,
    discrimination_power,
    power_curve,
    replicates_needed,
)
from ystwin.generator.design import decoupling_grid, stressor_only_series
from ystwin.generator.literature import parameters_from_literature

PARAMS = parameters_from_literature("UPRE2")
DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0)
COLLINEAR = stressor_only_series(PARAMS, doses=DOSES)
DECOUPLED = decoupling_grid(PARAMS, doses=DOSES, nutrient_factors=(1.0, 0.6, 0.35))


class TestDetectionPower:
    def test_a_real_induction_is_detected_with_enough_replicates(self):
        power = detection_power(PARAMS, DECOUPLED, induction_fold=2.0,
                                n_replicates=8, n_simulations=60, seed=0)

        assert power > 0.8

    def test_no_induction_gives_power_near_the_false_positive_rate(self):
        """A flat promoter must not be detected, or the test is measuring growth."""
        power = detection_power(PARAMS, DECOUPLED, induction_fold=1.0,
                                n_replicates=8, n_simulations=60, seed=0)

        assert power < 0.20

    def test_power_rises_with_replicate_count(self):
        low = detection_power(PARAMS, DECOUPLED, induction_fold=1.3,
                              n_replicates=2, n_simulations=60, seed=1)
        high = detection_power(PARAMS, DECOUPLED, induction_fold=1.3,
                               n_replicates=12, n_simulations=60, seed=1)

        assert high > low

    def test_power_rises_with_effect_size(self):
        small = detection_power(PARAMS, DECOUPLED, induction_fold=1.15,
                                n_replicates=4, n_simulations=60, seed=2)
        large = detection_power(PARAMS, DECOUPLED, induction_fold=2.5,
                                n_replicates=4, n_simulations=60, seed=2)

        assert large > small

    def test_power_is_a_probability(self):
        power = detection_power(PARAMS, DECOUPLED, induction_fold=1.5,
                                n_replicates=3, n_simulations=40, seed=3)

        assert 0.0 <= power <= 1.0


class TestTheDesignChangesWhatItCosts:
    """Detection is easy on any design; attribution is what collinearity destroys."""

    def test_both_designs_detect_a_dose_slope_easily(self):
        decoupled = replicates_needed(PARAMS, DECOUPLED, induction_fold=1.5,
                                      target_power=0.8, n_simulations=40, seed=0)
        collinear = replicates_needed(PARAMS, COLLINEAR, induction_fold=1.5,
                                      target_power=0.8, n_simulations=40, seed=0)

        assert decoupled is not None and collinear is not None

    def test_only_the_decoupled_design_can_attribute_it_to_dose(self):
        """With growth in the model too, the collinear design cannot separate them."""
        decoupled = discrimination_power(PARAMS, DECOUPLED, induction_fold=2.0,
                                         n_replicates=6, n_simulations=50, seed=0)
        collinear = discrimination_power(PARAMS, COLLINEAR, induction_fold=2.0,
                                         n_replicates=6, n_simulations=50, seed=0)

        assert decoupled > collinear + 0.3

    def test_a_chemostat_beats_a_plate_at_the_same_replicate_count(self):
        """Knowing mu exactly removes the error that propagates into the correction."""
        plate = discrimination_power(PARAMS, DECOUPLED, induction_fold=1.4,
                                     n_replicates=4, n_simulations=50, seed=1)
        chemostat = discrimination_power(PARAMS, DECOUPLED, induction_fold=1.4,
                                         n_replicates=4, n_simulations=50, seed=1,
                                         growth_known=True)

        assert chemostat >= plate

    def test_it_reports_the_cap_rather_than_a_fabricated_number(self):
        needed = replicates_needed(PARAMS, DECOUPLED, induction_fold=1.001,
                                   target_power=0.99, n_simulations=25,
                                   max_replicates=4, seed=0)

        assert needed is None


class TestPowerCurve:
    def test_it_returns_one_row_per_combination(self):
        curve = power_curve(PARAMS, DECOUPLED, induction_folds=(1.2, 2.0),
                            replicate_counts=(2, 4, 8), n_simulations=25, seed=0)

        assert len(curve) == 6
        assert set(curve.columns) >= {"induction_fold", "n_replicates", "power"}

    def test_power_increases_along_both_axes(self):
        curve = power_curve(PARAMS, DECOUPLED, induction_folds=(1.2, 2.5),
                            replicate_counts=(2, 10), n_simulations=40, seed=0)
        grid = curve.pivot_table(index="induction_fold", columns="n_replicates", values="power")

        assert grid.loc[2.5, 10] >= grid.loc[1.2, 2]

    def test_the_curve_records_which_design_produced_it(self):
        curve = power_curve(PARAMS, DECOUPLED, induction_folds=(2.0,),
                            replicate_counts=(4,), n_simulations=20, seed=0)

        assert curve.n_conditions.iloc[0] == len(DECOUPLED)
        assert 0.0 <= curve.collinearity.iloc[0] <= 1.0


class TestTheClosedFormMatchesTheIntegrator:
    """The power sweep uses a closed form; it must agree with the tested integrator."""

    @pytest.mark.parametrize("growth,k_deg,t", [(0.4, 0.0, 20.0), (0.1, 0.3, 5.0), (0.05, 0.0, 60.0)])
    def test_it_agrees_at_a_range_of_rates(self, growth, k_deg, t):
        from ystwin.analysis.power import _reporter_at
        from ystwin.reporter import ReporterKinetics, simulate_reporter

        grid = np.linspace(0.0, t, 400)
        activity = 1.2e-3
        integrated = simulate_reporter(
            grid, np.full_like(grid, activity), np.full_like(grid, growth),
            ReporterKinetics(k_deg=k_deg), r0=2.0e-3,
        )
        closed = _reporter_at(activity, growth, k_deg, t, initial=2.0e-3)

        assert closed == pytest.approx(float(integrated[-1]), rel=1e-4)

    def test_with_no_initial_value_it_starts_at_steady_state(self):
        from ystwin.analysis.power import _reporter_at

        steady = _reporter_at(1.2e-3, 0.3, 0.0, 0.0)

        assert steady == pytest.approx(1.2e-3 / 0.3)


class TestEstimatingGrowthDoesNotBiasTheAttribution:
    """Checked because it would invalidate the plate result if true.

    Error in a regressor attenuates its coefficient toward zero, which would leave
    spurious variance for dose and make a plate look better than a chemostat. It does
    not happen here: mu is estimated from 50 readings of a clean exponential, so its
    error is small against the spread of mu across the design. The consequence is
    practical -- the nutrient axis is what matters, not the chemostat.
    """

    def _coefficients(self, growth_known, fold=1.2, n_plates=6, trials=60, seed=0):
        from dataclasses import replace

        from ystwin.analysis.power import _observe, _scaled

        rng = np.random.default_rng(seed)
        scaled = _scaled(PARAMS, fold)
        out = []
        for _ in range(trials):
            rows = []
            for _ in range(n_plates):
                batch = replace(scaled, promoter_basal=scaled.promoter_basal
                                * float(np.exp(rng.normal(0, 0.12))))
                for condition in DECOUPLED:
                    signal, growth = _observe(condition, batch, rng, 0.02, growth_known)
                    rows.append((condition.dose_mM, growth, signal))
            dose, growth, signal = (np.asarray(x, float) for x in zip(*rows))
            design = np.column_stack([
                np.ones_like(dose), np.log1p(dose), np.log(np.clip(growth, 1e-9, None))
            ])
            beta, *_ = np.linalg.lstsq(design, np.log(np.clip(signal, 1e-30, None)), rcond=None)
            out.append(beta[1:])
        return np.asarray(out)

    def test_the_growth_coefficient_recovers_its_true_value_of_minus_one(self):
        for known in (True, False):
            coefficients = self._coefficients(known)
            assert coefficients[:, 1].mean() == pytest.approx(-1.0, abs=0.02)

    def test_estimating_growth_does_not_attenuate_its_coefficient(self):
        exact = self._coefficients(True)[:, 1].mean()
        estimated = self._coefficients(False)[:, 1].mean()

        assert abs(estimated) > abs(exact) - 0.02

    def test_the_dose_coefficient_is_the_same_either_way(self):
        exact = self._coefficients(True)[:, 0].mean()
        estimated = self._coefficients(False)[:, 0].mean()

        assert estimated == pytest.approx(exact, rel=0.15)


class TestDilutionInversionMatchesTheReporterModule:
    """The power sweep and the real inversion must agree, including for a degron.

    power.py inverts the reporter at steady state to recover activity, and
    reporter.promoter_activity does the same job on a real trace. If the two disagree,
    the replicate counts are for a different estimator than the one that will be used.

    They agreed for a long time by coincidence: the sweep dropped k_deg from the loss
    term, and every construct in the repository has k_deg = 0.0.
    """

    def test_a_stable_reporter_inverts_to_its_true_activity(self):
        from dataclasses import replace

        from ystwin.analysis.power import _reporter_at
        from ystwin.generator.design import stressor_only_series

        params = replace(PARAMS, k_deg=0.0)
        condition = stressor_only_series(params, doses=(1.0,))[0]
        activity = params.promoter_activity_at(condition.dose_mM)
        signal = _reporter_at(activity, condition.growth_rate, params.k_deg, 60.0)
        recovered = signal * (condition.growth_rate + params.k_deg)
        assert recovered == pytest.approx(activity, rel=1e-6)

    def test_a_degron_reporter_also_inverts_to_its_true_activity(self):
        """The case that was silently wrong. Fails if k_deg is dropped again."""
        from dataclasses import replace

        from ystwin.analysis.power import _reporter_at
        from ystwin.generator.design import stressor_only_series

        params = replace(PARAMS, k_deg=0.8)  # degron-tagged: loss dominates dilution
        condition = stressor_only_series(params, doses=(1.0,))[0]
        activity = params.promoter_activity_at(condition.dose_mM)
        signal = _reporter_at(activity, condition.growth_rate, params.k_deg, 60.0)

        correct = signal * (condition.growth_rate + params.k_deg)
        dropped = signal * condition.growth_rate  # the old expression

        assert correct == pytest.approx(activity, rel=1e-6)
        assert dropped != pytest.approx(activity, rel=0.1), (
            "dropping k_deg should be visibly wrong for a degron reporter"
        )

    def test_detection_power_still_runs_with_a_degron_reporter(self):
        from dataclasses import replace

        from ystwin.analysis.power import detection_power

        params = replace(PARAMS, k_deg=0.8)
        power = detection_power(params, DECOUPLED, induction_fold=2.0,
                                n_replicates=6, n_simulations=40, seed=0)
        assert 0.0 <= float(power) <= 1.0


def test_noiseless_flat_promoter_stays_null_despite_batch_variation():
    estimate = detection_power(PARAMS, COLLINEAR, induction_fold=1.0,
                               n_replicates=4, n_simulations=20, seed=3,
                               reader_cv=0, biological_cv=0.4, well_cv=0, growth_known=True)
    assert estimate.hits == 0
    assert "plate-slope" in estimate.method


def test_discrimination_uses_one_coefficient_per_biological_plate(monkeypatch):
    import ystwin.analysis.power as module

    seen = []
    monkeypatch.setattr(module, "_positive_plate_test", lambda values: seen.extend(values) or False)
    module._discrimination_trial(PARAMS, DECOUPLED, 1.4, 4,
                                 np.random.default_rng(3), 0.02, 0.12, True)
    assert len(seen) == 4


def test_duplicating_rows_cannot_increase_the_plate_coefficient_count():
    from ystwin.analysis.power import _plate_coefficients

    dose = np.tile([0.0, 0.2, 0.5, 1.0], 3)
    growth = np.tile([0.1, 0.3, 0.2, 0.25], 3)
    plate = np.repeat(np.arange(3), 4)
    slope = np.repeat([0.2, 0.4, 0.1], 4)
    signal = np.exp(slope * np.log1p(dose) - np.log(growth) + plate)
    original = _plate_coefficients(dose, growth, signal, plate)
    duplicated = _plate_coefficients(*(np.tile(values, 10) for values in (dose, growth, signal, plate)))
    assert len(original) == len(duplicated) == 3
    assert duplicated == pytest.approx(original)


def test_replication_search_adjusts_monte_carlo_confidence_for_its_search(monkeypatch):
    import ystwin.analysis.power as module

    levels = []

    def estimate(*args, confidence=0.95, **kwargs):
        levels.append(confidence)
        return module.PowerEstimate(0, 20, confidence)

    monkeypatch.setattr(module, "detection_power", estimate)
    assert module.replicates_needed(PARAMS, COLLINEAR, 1.2, max_replicates=5) is None
    assert levels == pytest.approx([0.99] * 5)


def test_adaptive_power_does_not_exceed_its_budget_when_initial_batch_is_larger(monkeypatch):
    import ystwin.analysis.power as module

    calls = []
    monkeypatch.setattr(module, "_one_trial", lambda *args: calls.append(1) or False)
    result = module.detection_power(PARAMS, COLLINEAR, 1.2, 3, n_simulations=50,
                                    max_simulations=10, tolerance=0.01, seed=1)
    assert len(calls) == result.n_simulations == 10


@pytest.mark.parametrize("hits,n,confidence", [(-1, 10, 0.95), (11, 10, 0.95), (1, 10, 1.0)])
def test_impossible_power_estimates_are_refused(hits, n, confidence):
    from ystwin.analysis.power import PowerEstimate

    with pytest.raises(ValueError):
        PowerEstimate(hits, n, confidence)


def test_discrimination_controls_total_reporter_loss_not_only_dilution(monkeypatch):
    from dataclasses import replace
    import ystwin.analysis.power as module

    captured = []

    def coefficients(dose, loss, signal, plate):
        captured.extend(loss)
        return [0.0] * 3

    monkeypatch.setattr(module, "_plate_coefficients", coefficients)
    params = replace(PARAMS, k_deg=0.8)
    module._discrimination_trial(params, DECOUPLED, 1.0, 3,
                                 np.random.default_rng(4), 0.0, 0.1, True)
    expected = np.tile([condition.growth_rate + params.k_deg for condition in DECOUPLED], 3)
    assert captured == pytest.approx(expected)
