import numpy as np
import pytest

from ystwin.diagnostics.dilution_confound import dilution_confound_report
from ystwin.reporter import ReporterKinetics, simulate_reporter
from ystwin.readings import RawOD, RawRFU


def _logistic_od(t, od0=0.05, cap=2.0, rate=0.45):
    return cap / (1 + (cap / od0 - 1) * np.exp(-rate * t))


def test_flags_a_pure_dilution_artefact_as_almost_entirely_growth_explained():
    """Constant promoter, saturating culture: the apparent induction is all dilution."""
    t = np.linspace(0, 30, 600)
    od = _logistic_od(t)
    mu = np.gradient(np.log(od), t)
    kin = ReporterKinetics(k_deg=0.0)
    per_cell = simulate_reporter(t, np.ones_like(t), np.clip(mu, 1e-3, None), kin, r0=1 / mu[0])
    rfu = per_cell * od

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=kin)

    assert rep.dilution_explained_r2 > 0.9
    assert rep.naive_fold_change > 3.0
    assert rep.activity_fold_change < 1.5
    assert rep.verdict == "dilution-dominated"


def test_does_not_flag_a_genuine_induction_at_constant_growth_rate():
    t = np.linspace(0, 30, 600)
    od = 0.05 * np.exp(0.2 * t)
    mu = np.full_like(t, 0.2)
    kin = ReporterKinetics(k_deg=0.0)
    k_syn = 1.0 + 4.0 / (1 + np.exp(-(t - 15.0)))
    rfu = simulate_reporter(t, k_syn, mu, kin, r0=k_syn[0] / mu[0]) * od

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=kin)

    assert rep.activity_fold_change > 3.0
    assert rep.verdict == "promoter-driven"


def test_separates_a_real_induction_that_coincides_with_a_growth_slowdown():
    """The case that matters: true induction AND slowing growth at the same time."""
    t = np.linspace(0, 30, 900)
    od = _logistic_od(t)
    mu = np.clip(np.gradient(np.log(od), t), 1e-3, None)
    kin = ReporterKinetics(k_deg=0.0)
    k_syn = 1.0 + 2.0 / (1 + np.exp(-(t - 14.0)))
    rfu = simulate_reporter(t, k_syn, mu, kin, r0=k_syn[0] / mu[0]) * od

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=kin)

    assert rep.activity_fold_change == pytest.approx(3.0, rel=0.20)
    assert rep.naive_fold_change > rep.activity_fold_change * 2


def test_background_and_blank_corrections_are_applied_before_any_ratio():
    t = np.linspace(0, 20, 400)
    od = 0.05 * np.exp(0.25 * t)
    kin = ReporterKinetics(k_deg=0.0)
    per_cell = simulate_reporter(t, np.ones_like(t), np.full_like(t, 0.25), kin)
    clean = dilution_confound_report(t, RawOD(od), RawRFU(per_cell * od), kinetics=kin)
    offset = dilution_confound_report(
        t, RawOD(od + 0.04), RawRFU(per_cell * od + 350.0), kinetics=kin, od_blank=0.04, rfu_background=350.0
    )

    assert offset.activity_fold_change == pytest.approx(clean.activity_fold_change, rel=1e-6)


def test_report_carries_the_per_timepoint_series_for_plotting_and_audit():
    t = np.linspace(0, 20, 400)
    od = 0.05 * np.exp(0.25 * t)
    rfu = od * 4.0

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics())

    for series in (rep.growth_rate, rep.naive_specific, rep.promoter_activity):
        assert series.shape == t.shape


def test_wells_that_never_grow_are_rejected_rather_than_reported_as_stressed():
    t = np.linspace(0, 20, 400)
    od = np.full_like(t, 0.041) + np.random.default_rng(0).normal(0, 1e-4, t.shape)

    with pytest.raises(ValueError, match="did not grow"):
        dilution_confound_report(t, RawOD(od), RawRFU(od * 100), kinetics=ReporterKinetics(), od_blank=0.04)


def _shutoff_trace(k_deg_true, mu_const=0.25, t_off=8.0, floor=0.0):
    """Promoter switches off; reporter then decays at (mu + k_deg_true).

    ``floor`` sets the residual activity. With floor=0 the late decay is a clean
    exponential, which is what makes the loss-rate bound recoverable; with a
    positive floor the true activity never reaches zero, which is what makes the
    negative-activity fraction meaningful.
    """
    t = np.linspace(0, 24, 900)
    od = 0.05 * np.exp(mu_const * t)
    mu = np.full_like(t, mu_const)
    k_syn = floor + (1.0 - floor) / (1 + np.exp((t - t_off) / 0.35))
    kin_true = ReporterKinetics(k_deg=k_deg_true)
    per_cell = simulate_reporter(t, k_syn, mu, kin_true, r0=1.0 / (mu_const + k_deg_true))
    return t, od, per_cell * od


def test_activity_below_zero_is_reported_as_model_inadequacy_not_clipped_away():
    """Signal falling faster than dilution means the assumed k_deg is wrong."""
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert rep.negative_activity_fraction > 0.4
    assert rep.verdict == "model-inadequate"


def test_reports_the_minimum_reporter_loss_rate_the_data_demand():
    """k_synth >= 0 implies k_deg >= -(dR/dt)/R - mu; recover a known true value."""
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert rep.implied_min_k_deg == pytest.approx(0.30, abs=0.03)


def test_the_implied_loss_bound_is_a_property_of_the_data_not_of_the_assumption():
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    assuming_none = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))
    assuming_true = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.30))

    assert assuming_none.implied_min_k_deg == pytest.approx(assuming_true.implied_min_k_deg, rel=1e-6)


def test_excess_loss_rate_measures_how_far_the_assumed_kinetics_fall_short():
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    assuming_none = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))
    assuming_true = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.30))

    assert assuming_none.excess_loss_rate == pytest.approx(0.30, abs=0.05)
    assert assuming_true.excess_loss_rate <= 0.03


def test_the_right_kinetics_leave_far_less_of_the_trace_implying_negative_activity():
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    assuming_none = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))
    assuming_true = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.30))

    assert assuming_true.negative_activity_fraction < assuming_none.negative_activity_fraction / 3


def test_activity_that_stays_genuinely_positive_is_never_flagged_as_negative():
    """With the right kinetics and a promoter that never fully shuts off."""
    t, od, rfu = _shutoff_trace(k_deg_true=0.30, floor=0.25)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.30))

    assert rep.negative_activity_fraction == pytest.approx(0.0, abs=0.01)
    assert rep.verdict != "model-inadequate"


def test_fold_change_is_undefined_rather_than_astronomical_when_activity_is_negative():
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert np.isnan(rep.activity_fold_change)


def test_quasi_steady_state_points_are_counted_and_exposed():
    t = np.linspace(0, 20, 400)
    od = 0.05 * np.exp(0.25 * t)
    rfu = od * 4.0

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics())

    assert rep.qss_fraction == pytest.approx(1.0)


def test_dilution_r2_is_withheld_when_relaxation_is_slower_than_the_experiment():
    """S -> k_synth/(mu+k_deg) only holds if the reporter can relax within the run.

    A stalled culture with a stable FP never reaches that limit, so the ratio
    carries no information about dilution and R^2 must not be reported.
    """
    t = np.linspace(0, 12, 400)
    od = 0.4 * np.exp(0.02 * t)  # crawls: mu = 0.02 /h, relaxation 50 h vs a 12 h run
    rfu = od * np.linspace(4000, 1500, t.size)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert rep.qss_fraction < 0.2
    assert np.isnan(rep.dilution_explained_r2)


def test_r2_is_estimated_only_from_the_points_where_the_limit_is_valid():
    """Growth phase (valid) followed by a stalled tail (invalid) must not be pooled."""
    t = np.linspace(0, 40, 800)
    mu_true = np.where(t < 12, 0.40, 0.004)
    od = 0.05 * np.exp(np.concatenate([[0], np.cumsum(mu_true[:-1] * np.diff(t))]))
    kin = ReporterKinetics(k_deg=0.0)
    per_cell = simulate_reporter(t, np.ones_like(t), np.clip(mu_true, 1e-4, None), kin, r0=1 / 0.40)

    rep = dilution_confound_report(t, RawOD(od), RawRFU(per_cell * od), kinetics=kin)

    assert 0.0 < rep.qss_fraction < 0.6
    assert np.isfinite(rep.dilution_explained_r2)


@pytest.mark.parametrize("window_h", [0.6, 1.0, 2.0, 4.0, 8.0])
def test_the_implied_loss_bound_is_stable_across_smoothing_windows(window_h):
    """It must estimate a decay rate, not the curvature of the smoother.

    Fitting a polynomial to an exponentially decaying signal and dividing by the
    signal is wildly window-dependent; differentiating the logarithm is not.
    """
    t, od, rfu = _shutoff_trace(k_deg_true=0.30)

    rep = dilution_confound_report(
        t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0), activity_window_h=window_h
    )

    assert rep.implied_min_k_deg == pytest.approx(0.30, abs=0.05)


class TestNotEstimableIsNotMixed:
    """A statistic that could not be computed is not an ambiguous result.

    The dilution R^2 needs quasi-steady state, mu + k_deg > 2/duration. Below that no
    timepoint qualifies, the R^2 is NaN, and the verdict used to fall through to
    "mixed" -- reporting a finding where there was no measurement. A 48 h run needs
    mu > 0.042 /h, so slow cultures and most mammalian lines fall through it.
    """

    def _trace(self, mu, duration=48.0, n=49):
        """Unpacked with `*` into `dilution_confound_report`, so the readings are typed
        here rather than at the four call sites."""
        t = np.linspace(0.0, duration, n)
        od = 0.05 * np.exp(mu * t)
        rfu = 1000 * od * (1 + 0.01 * np.random.default_rng(0).normal(size=t.size))
        return t, RawOD(od), RawRFU(rfu)

    def test_a_slow_culture_reports_not_estimable(self):
        report = dilution_confound_report(*self._trace(0.03),
                                          kinetics=ReporterKinetics(k_deg=0.0))
        assert report.qss_fraction == 0.0
        assert not np.isfinite(report.dilution_explained_r2)
        assert report.verdict == "not-estimable"

    def test_a_fast_culture_still_gets_a_real_verdict(self):
        report = dilution_confound_report(*self._trace(0.30),
                                          kinetics=ReporterKinetics(k_deg=0.0))
        assert report.qss_fraction > 0
        assert report.verdict != "not-estimable"

    def test_a_degron_reporter_rescues_a_slow_culture(self):
        """The threshold is on mu + k_deg, so active degradation restores the test."""
        slow = self._trace(0.03)
        stable = dilution_confound_report(*slow, kinetics=ReporterKinetics(k_deg=0.0))
        degron = dilution_confound_report(*slow, kinetics=ReporterKinetics(k_deg=0.5))
        assert stable.verdict == "not-estimable"
        assert degron.qss_fraction > 0
        assert degron.verdict != "not-estimable"


def test_model_inadequate_outranks_not_estimable():
    """Negative activity refutes the kinetics whether or not the regression ran.

    negative_activity_fraction comes from promoter_activity and reads neither the R^2
    nor the QSS fraction, so it is determinable exactly when they are not. Placing the
    not-estimable guard above it shadowed the strongest real-data finding here: on the
    2026-07-01 plate all 53 G1-passing wells have 57-100% negative inferred activity
    and zero QSS fraction, and were reported as unmeasurable rather than as refuting
    the stable-reporter assumption.
    """
    t = np.linspace(0.0, 48.0, 49)
    od = 0.05 * np.exp(0.03 * t)                    # too slow for quasi-steady state
    rfu = 1000 * od * np.exp(-0.08 * t)             # falling signal -> negative activity
    report = dilution_confound_report(t, RawOD(od), RawRFU(rfu), kinetics=ReporterKinetics(k_deg=0.0))

    assert report.qss_fraction == 0.0
    assert not np.isfinite(report.dilution_explained_r2)
    assert report.negative_activity_fraction > 0.25
    assert report.verdict == "model-inadequate"
