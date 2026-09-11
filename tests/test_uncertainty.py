"""The interval has to be right, and it has to be absent when it cannot be right.

Two properties carry these tests. A cluster bootstrap must be *wider* than a naive
well-level one on clustered data, because that is the entire reason to use it; and the
module must refuse an interval at two plates rather than emit a narrow one, because a
narrow interval there would be the most damaging output the module could produce.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.uncertainty import (
    MIN_PLATES_FOR_INTERVAL,
    activity_uncertainty,
    equivalence,
    fold_change,
    window_sensitivity,
)
from ystwin.reporter import ReporterKinetics
from ystwin.readings import CorrectedOD, SpecificFluorescence


@pytest.fixture(scope="module")
def sweep_table():
    """The committed late-window sweep. Tracked, so this never skips."""
    from ystwin import paths

    return pd.read_csv(paths.outputs_dir() / "late_window_sensitivity.csv")


def _clustered_plate(
    n_plates: int,
    n_wells: int,
    true_fold: float,
    plate_cv: float,
    well_cv: float,
    seed: int,
    fold_cv: float = 0.0,
    control_level: float = 1000.0,
) -> pd.DataFrame:
    """Synthetic readings with a known fold and three separable variance components.

    ``plate_cv`` shifts a whole plate's level, hitting the dosed and control wells
    equally. ``fold_cv`` perturbs the plate's *fold* -- a plate-by-dose interaction,
    which is what differing inoculum density and therefore differing growth actually
    produce. The distinction matters: a shared multiplicative shift cancels in a ratio
    and an interaction does not, so only ``fold_cv`` should widen a fold's interval.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for p in range(n_plates):
        batch = float(np.exp(rng.normal(0.0, plate_cv)))
        fold = true_fold * float(np.exp(rng.normal(0.0, fold_cv)))
        for dose, level in ((0.0, control_level), (1.0, control_level * fold)):
            for w in range(n_wells):
                value = level * batch * float(np.exp(rng.normal(0.0, well_cv)))
                rows.append({
                    "plate": f"P{p}", "well": f"w{w}", "rep": w, "construct": "R",
                    "dose_mM": dose, "activity_late": value, "naive_late": value * 1.5,
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# refusal at low cluster counts -- the property that protects the real result
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_plates", [1, 2])
def test_interval_is_refused_below_the_minimum(n_plates):
    readings = _clustered_plate(n_plates, 3, true_fold=1.4, plate_cv=0.2, well_cv=0.1, seed=0)
    fold = fold_change(readings, "R", 1.0)
    assert not fold.estimable
    assert np.isnan(fold.low) and np.isnan(fold.high)
    assert str(n_plates) in fold.method
    # The point estimate is still reported; only the interval is withheld.
    assert np.isfinite(fold.point)


def test_refused_interval_never_claims_an_effect():
    """A missing interval must not read as evidence, in either direction."""
    readings = _clustered_plate(2, 3, true_fold=3.0, plate_cv=0.2, well_cv=0.1, seed=1)
    fold = fold_change(readings, "R", 1.0)
    assert fold.point > 2.0          # a large apparent effect
    assert not fold.excludes_unity   # and still no claim
    assert equivalence(fold).verdict == "INCONCLUSIVE"


def test_interval_appears_at_the_minimum():
    readings = _clustered_plate(MIN_PLATES_FOR_INTERVAL, 3, true_fold=1.4,
                                plate_cv=0.15, well_cv=0.1, seed=2)
    fold = fold_change(readings, "R", 1.0)
    assert fold.estimable
    assert fold.low < fold.point < fold.high


# ---------------------------------------------------------------------------
# clustering
# ---------------------------------------------------------------------------


def test_a_shared_plate_shift_cancels_in_the_ratio():
    """A multiplicative effect hitting both doses equally divides out.

    This is the same cancellation the co-expression design relies on, and it is worth
    locking in: it means plate-level *level* variation is not what threatens a fold
    change, so an interval that grew with it would be measuring the wrong thing.
    """
    shifted = _clustered_plate(6, 8, 1.3, plate_cv=0.40, well_cv=0.05, seed=3)
    flat = _clustered_plate(6, 8, 1.3, plate_cv=0.0, well_cv=0.05, seed=3)
    wide = fold_change(shifted, "R", 1.0, seed=11)
    narrow = fold_change(flat, "R", 1.0, seed=11)
    assert wide.point == pytest.approx(narrow.point, rel=1e-9)
    assert (wide.high - wide.low) == pytest.approx(narrow.high - narrow.low, rel=0.05)


def test_cluster_bootstrap_is_wider_when_the_fold_itself_varies_by_plate():
    """Plate-by-dose interaction is the variance component that must reach the interval."""
    interacting = _clustered_plate(6, 8, 1.3, plate_cv=0.05, well_cv=0.05,
                                   seed=3, fold_cv=0.30)
    flat = _clustered_plate(6, 8, 1.3, plate_cv=0.05, well_cv=0.05, seed=3, fold_cv=0.0)
    wide = fold_change(interacting, "R", 1.0, seed=11)
    narrow = fold_change(flat, "R", 1.0, seed=11)
    assert (wide.high - wide.low) > 2.0 * (narrow.high - narrow.low)


def test_ignoring_clusters_would_understate_the_interval():
    """Treating each well-pair as its own biological replicate narrows the interval.

    That is the error being avoided: it turns 6 plates into 48 pseudo-replicates and
    buys precision that was never measured.
    """
    readings = _clustered_plate(6, 8, 1.3, plate_cv=0.05, well_cv=0.05, seed=4, fold_cv=0.30)
    honest = fold_change(readings, "R", 1.0, seed=12)

    naive_frame = readings.copy()
    naive_frame["plate"] = naive_frame.plate + "-" + naive_frame.rep.astype(str)
    naive = fold_change(naive_frame, "R", 1.0, seed=12)
    assert (naive.high - naive.low) < (honest.high - honest.low)


def test_coverage_is_near_nominal_on_clustered_data():
    """Across repeated worlds the interval should contain the truth about 95% of the time."""
    true_fold, hits, trials = 1.35, 0, 60
    for seed in range(trials):
        readings = _clustered_plate(6, 6, true_fold=true_fold, plate_cv=0.2,
                                    well_cv=0.08, seed=100 + seed, fold_cv=0.15)
        fold = fold_change(readings, "R", 1.0, n_resamples=600, seed=seed)
        if fold.estimable and fold.low <= true_fold <= fold.high:
            hits += 1
    # Loose bound: this is a small-cluster bootstrap and undercoverage is expected.
    assert 0.80 <= hits / trials <= 1.0


# ---------------------------------------------------------------------------
# equivalence
# ---------------------------------------------------------------------------


def test_equivalence_detects_a_genuine_null_at_adequate_replication():
    readings = _clustered_plate(12, 8, true_fold=1.0, plate_cv=0.03, well_cv=0.03, seed=5)
    verdict = equivalence(fold_change(readings, "R", 1.0, seed=13), margin=(0.9, 1.1))
    assert verdict.verdict == "EQUIVALENT"


def test_equivalence_detects_a_real_effect():
    readings = _clustered_plate(12, 8, true_fold=1.8, plate_cv=0.05, well_cv=0.05, seed=6)
    verdict = equivalence(fold_change(readings, "R", 1.0, seed=14), margin=(0.9, 1.1))
    assert verdict.verdict == "DIFFERENT"


def test_equivalence_is_inconclusive_when_the_interval_straddles_the_margin():
    readings = _clustered_plate(3, 3, true_fold=1.05, plate_cv=0.35, well_cv=0.25,
                                seed=7, fold_cv=0.35)
    verdict = equivalence(fold_change(readings, "R", 1.0, seed=15), margin=(0.9, 1.1))
    assert verdict.verdict == "INCONCLUSIVE"


def test_margin_must_bracket_no_induction():
    readings = _clustered_plate(4, 4, true_fold=1.0, plate_cv=0.1, well_cv=0.1, seed=8)
    fold = fold_change(readings, "R", 1.0)
    with pytest.raises(ValueError, match="must bracket 1.0"):
        equivalence(fold, margin=(1.2, 1.5))


# ---------------------------------------------------------------------------
# input validation -- refuse rather than guess
# ---------------------------------------------------------------------------


def test_missing_columns_are_refused():
    with pytest.raises(ValueError, match="missing required columns"):
        fold_change(pd.DataFrame({"plate": ["P0"], "construct": ["R"]}), "R", 1.0)


def test_absent_construct_is_refused():
    readings = _clustered_plate(4, 3, 1.2, 0.1, 0.1, seed=9)
    with pytest.raises(ValueError, match="no rows for construct"):
        fold_change(readings, "NotAConstruct", 1.0)


def test_absent_dose_is_refused():
    readings = _clustered_plate(4, 3, 1.2, 0.1, 0.1, seed=10)
    with pytest.raises(ValueError, match="no rows at dose"):
        fold_change(readings, "R", 99.0)


def test_dose_and_control_on_different_plates_is_refused():
    """A fold across plates confounds the batch effect with the dose."""
    readings = _clustered_plate(4, 3, 1.2, 0.1, 0.1, seed=11)
    split = readings.copy()
    split.loc[split.dose_mM == 0.0, "plate"] = "OTHER"
    with pytest.raises(ValueError, match="no plate carrying both"):
        fold_change(split, "R", 1.0)


def test_is_deterministic_given_a_seed():
    readings = _clustered_plate(5, 4, 1.3, 0.2, 0.1, seed=12)
    a = fold_change(readings, "R", 1.0, seed=42)
    b = fold_change(readings, "R", 1.0, seed=42)
    assert (a.low, a.high) == (b.low, b.high)


# ---------------------------------------------------------------------------
# propagating input error into the recovered activity
# ---------------------------------------------------------------------------


def _trace(mu: float, activity: float, noise: float, seed: int, n: int = 120):
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 24.0, n)
    od = 0.1 * np.exp(mu * t) * (1.0 + rng.normal(0.0, noise, n))
    od = np.clip(od, 1e-6, None)
    reporter = np.full(n, activity / mu) * (1.0 + rng.normal(0.0, noise, n))
    # Typed here rather than at the four call sites, which unpack this with `*`.
    return t, SpecificFluorescence(reporter), CorrectedOD(od), np.full(n, mu)


def test_activity_uncertainty_grows_with_noisier_input():
    quiet = activity_uncertainty(*_trace(0.3, 100.0, 0.005, seed=1), n_draws=120)
    loud = activity_uncertainty(*_trace(0.3, 100.0, 0.05, seed=1), n_draws=120)
    assert np.median(loud.sigma) > np.median(quiet.sigma)
    assert loud.reporter_cv > quiet.reporter_cv


def test_growth_rate_standard_error_is_reported_and_finite():
    result = activity_uncertainty(*_trace(0.3, 100.0, 0.02, seed=2), n_draws=100)
    assert np.isfinite(result.growth_rate_se)
    assert result.growth_rate_se > 0


def test_delta_and_monte_carlo_agree_in_the_quiet_regime():
    args = _trace(0.3, 100.0, 0.004, seed=3)
    delta = activity_uncertainty(*args, method="delta")
    mc = activity_uncertainty(*args, method="monte-carlo", n_draws=300)
    ratio = np.median(delta.sigma) / max(np.median(mc.sigma), 1e-12)
    assert 0.2 < ratio < 5.0  # same order; they are different approximations


def test_band_widens_with_z():
    result = activity_uncertainty(*_trace(0.3, 100.0, 0.02, seed=4), n_draws=100)
    lo1, hi1 = result.band(z=1.0)
    lo2, hi2 = result.band(z=2.0)
    assert np.all(hi2 - lo2 >= hi1 - lo1)


def test_non_positive_optical_density_is_refused():
    t, reporter, od, mu = _trace(0.3, 100.0, 0.01, seed=5)
    # `_trace` returns readings, so the bad value goes into the array the reading holds.
    od.values[3] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        activity_uncertainty(t, reporter, od, mu)


def test_shape_mismatch_is_refused():
    t, reporter, od, mu = _trace(0.3, 100.0, 0.01, seed=6)
    with pytest.raises(ValueError, match="share a shape"):
        activity_uncertainty(t, reporter, CorrectedOD(od.values[:-1]), mu)


def test_unknown_method_is_refused():
    t, reporter, od, mu = _trace(0.3, 100.0, 0.01, seed=7)
    with pytest.raises(ValueError, match="unknown method"):
        activity_uncertainty(t, reporter, od, mu, method="bootstrap")


def test_maturation_kinetics_are_accepted():
    t, reporter, od, mu = _trace(0.3, 100.0, 0.01, seed=8)
    result = activity_uncertainty(t, reporter, od, mu,
                                  kinetics=ReporterKinetics(k_mat=2.0), n_draws=60)
    assert result.sigma.shape == t.shape


# ---------------------------------------------------------------------------
# window sensitivity
# ---------------------------------------------------------------------------


def test_window_sensitivity_reports_both_columns():
    readings = _clustered_plate(4, 4, 1.4, 0.15, 0.1, seed=13)
    table = window_sensitivity(readings, "R", 1.0)
    assert set(table.column) == {"activity_late", "naive_late"}
    assert len(table) == 2


def test_window_sensitivity_refuses_when_no_column_present():
    readings = _clustered_plate(4, 4, 1.4, 0.15, 0.1, seed=14)
    with pytest.raises(ValueError, match="none of"):
        window_sensitivity(readings, "R", 1.0, value_columns=("absent_column",))


class TestAFoldWithNoPointEstimateIsRefused:
    """A dying culture is not a repressed reporter, and it used to come back as one.

    At 2.0 mM H2O2 the AlteredYap1 wells lose total fluorescence faster than growth dilutes
    it, so the recovered activity is negative and `_combine`'s geometric mean over positive
    per-plate folds does not exist. The bootstrap's own guard only requires half the
    RESAMPLES to be positive, so it passed, and `fold_change` returned `estimable=True` with
    `point=nan` and an interval of [0.003, 0.585].

    `excludes_unity` reads `low > 1.0 or high < 1.0` and knows nothing about the point, so
    that row reported as a confident finding of strong repression -- from a well whose naive
    fold is -0.22.
    """

    @staticmethod
    def _readings(dosed_values, control_value=100.0):
        rows = []
        for plate, dosed in zip(("p1", "p2", "p3"), dosed_values):
            for well, value in (("A1", dosed), ("A2", dosed)):
                rows.append({"plate": plate, "construct": "X", "dose_mM": 1.0,
                             "well": well, "activity_late": value})
            for well in ("B1", "B2"):
                rows.append({"plate": plate, "construct": "X", "dose_mM": 0.0,
                             "well": well, "activity_late": control_value})
        return pd.DataFrame(rows)

    def test_every_plate_negative_is_refused_rather_than_bounded(self):
        got = fold_change(self._readings([-30.0, -25.0, -40.0]), "X", 1.0)

        assert not got.estimable
        assert "no positive per-plate fold" in got.method

    def test_and_carries_no_interval_to_be_misread(self):
        got = fold_change(self._readings([-30.0, -25.0, -40.0]), "X", 1.0)

        assert np.isnan(got.low) and np.isnan(got.high)

    def test_so_it_cannot_be_reported_as_a_call(self):
        """The assertion that matters: the old behaviour made `excludes_unity` True."""
        got = fold_change(self._readings([-30.0, -25.0, -40.0]), "X", 1.0)

        assert not got.excludes_unity

    def test_an_ordinary_fold_is_untouched_by_the_guard(self):
        got = fold_change(self._readings([150.0, 145.0, 160.0]), "X", 1.0)

        assert got.estimable
        assert got.point == pytest.approx(1.51, abs=0.02)

    def test_no_committed_row_has_an_interval_without_a_point(self, sweep_table):
        """The invariant, over the whole shipped table rather than the one row that
        exposed it. An estimable fold must have a fold; the defect was 144 rows in which
        one of them did not, and nothing looked.

        Stated this way on purpose. The row that exposed it -- AlteredYap1 at 2.0 mM H2O2 --
        is still estimable at the two widest windows, and correctly so: there one plate does
        return a positive fold, so a point estimate exists. It is a preposterous point
        estimate, 0.03 and 0.01, and both intervals straddle 1.0, so it claims nothing. The
        guard is for the case where NO plate is positive and the geometric mean therefore
        does not exist at all.
        """
        estimable = sweep_table[sweep_table.estimable]

        assert len(estimable) > 0
        assert estimable.fold.notna().all()

    def test_and_the_row_that_exposed_it_claims_nothing_at_any_window(self, sweep_table):
        row = sweep_table[(sweep_table.construct == "AlteredYap1")
                          & np.isclose(sweep_table.dose_mM, 2.0)]
        calls = row[row.estimable & ((row.low > 1.0) | (row.high < 1.0))]

        assert len(row) == 6
        assert calls.empty, "a culture dying under peroxide is being reported as a call"


def test_interval_records_target_level_and_two_stage_assumptions():
    readings = _clustered_plate(6, 4, 1.2, 0.1, 0.1, seed=12)
    fold = fold_change(readings, "R", 1.0, confidence=0.90, n_resamples=200)
    assert fold.confidence == 0.90
    assert fold.resampling == "two_stage"
    assert "geometric mean" in fold.target
    assert "approximate" in fold.coverage_status
    assert any("exchangeable wells" in assumption for assumption in fold.assumptions)
    assert equivalence(fold, margin=(0.9, 1.1)).one_sided_alpha == pytest.approx(0.05)


def test_duplicate_summaries_of_a_well_do_not_create_replicates():
    readings = _clustered_plate(6, 4, 1.2, 0.1, 0.1, seed=3)
    original = fold_change(readings, "R", 1.0, n_resamples=200)
    duplicated = fold_change(pd.concat([readings, readings.iloc[[0]]]), "R", 1.0, n_resamples=200)
    assert duplicated.n_wells == original.n_wells
    assert duplicated.point == original.point
    assert (duplicated.low, duplicated.high) == (original.low, original.high)


def test_invalid_plates_cannot_reenter_bootstrap_after_being_excluded_from_the_point():
    readings = _clustered_plate(4, 3, 1.2, 0.1, 0.1, seed=3)
    invalid = readings[readings.plate == "P0"].copy()
    invalid["plate"] = "invalid"
    invalid.loc[invalid.dose_mM == 1.0, "activity_late"] = [-1000, -1000, 100]
    original = fold_change(readings, "R", 1.0, n_resamples=200)
    combined = fold_change(pd.concat([readings, invalid]), "R", 1.0, n_resamples=200)
    assert combined.n_plates == 5
    assert combined.n_plates_contributing == 4
    assert combined.point == original.point
    assert (combined.low, combined.high) == (original.low, original.high)


def test_two_negative_activities_do_not_make_a_positive_fold():
    readings = _clustered_plate(4, 3, 1.2, 0.1, 0.1, seed=3)
    readings["activity_late"] *= -1
    result = fold_change(readings, "R", 1.0, n_resamples=200)
    assert not result.estimable
    assert result.n_plates_contributing == 0


def test_including_unity_is_not_a_no_response_claim():
    from ystwin.analysis.uncertainty import FoldChange

    fold = FoldChange("R", 1.0, 1.0, 0.5, 2.0, "constructed", 6, 6, 12, 200, True)
    assert not fold.excludes_unity
    assert equivalence(fold, margin=(0.9, 1.1)).verdict == "INCONCLUSIVE"


def test_reference_identity_is_not_independent_equivalence_evidence():
    readings = _clustered_plate(4, 3, 1.0, 0.1, 0.1, seed=3)
    fold = fold_change(readings, "R", 0.0, n_resamples=200)
    assert fold.is_reference
    assert equivalence(fold).verdict == "INCONCLUSIVE"


def test_bootstrap_coverage_check_counts_refusals_instead_of_calling_them_covered():
    from ystwin.analysis.uncertainty import fold_change_coverage

    result = fold_change_coverage(2, n_trials=4, n_resamples=50)
    assert result["n_issued"] == result["n_covered"] == 0
    assert result["n_refused"] == result["n_trials"] == 4
    assert np.isnan(result["actual_coverage"])
    assert "not real-data calibration" in result["coverage_scope"]


def test_bootstrap_coverage_returns_actual_frequency_and_monte_carlo_error():
    from ystwin.analysis.uncertainty import fold_change_coverage

    result = fold_change_coverage(6, n_trials=12, n_resamples=100, seed=8)
    assert result["n_issued"] == 12
    assert result["actual_coverage"] == result["n_covered"] / 12
    assert result["coverage_mc_low"] < result["actual_coverage"] <= result["coverage_mc_high"]
    assert result["nominal_coverage"] == 0.95


def test_activity_noise_is_estimated_on_the_numerator_before_od_division(monkeypatch):
    import ystwin.analysis.uncertainty as module

    t, reporter, od, mu = _trace(0.3, 100.0, 0.01, seed=2)
    seen = []
    original = module._reporter_noise_cv

    def record(times, values, window):
        seen.append(np.asarray(values).copy())
        return original(times, values, window)

    monkeypatch.setattr(module, "_reporter_noise_cv", record)
    module.activity_uncertainty(t, reporter, od, mu, n_draws=20)
    assert any(np.allclose(values, reporter.array * od.array) for values in seen)
