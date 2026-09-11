"""State estimation over the joint biomass / reporter / activity state.

This is §8.3 made concrete, on the corrected §8.1 state. Promoter activity is
inferred as a state with uncertainty rather than read off RFU/OD, and growth rate is
inferred jointly rather than assumed, so the dilution term is not conditioned on a
point estimate that is itself in doubt.

Two requirements from the plan the filter has to honour literally:
  * a missing sensor is skipped, never substituted with zero
  * only data at or before the update time may enter
"""

from dataclasses import replace

import numpy as np
import pytest

from ystwin.estimator import Observation, ParticleFilter, TwinPriors
from ystwin.observation import ReporterOptics

OPTICS = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                        inner_filter_coeff=0.0)
PRIORS = TwinPriors(
    biomass=(0.05, 0.01), reporter=(4e-3, 1e-3),
    promoter_activity=(1.0e-3, 5e-4), growth_rate=(0.30, 0.08),
    k_deg=0.05,
)


def _truth(hours=18.0, n=110, mu=0.28, k_syn=1.2e-3, k_deg=0.05, seed=0):
    """Simulate a culture, then observe it through the calibrated optics."""
    from ystwin.reporter import ReporterKinetics, simulate_reporter

    rng = np.random.default_rng(seed)
    t = np.linspace(0, hours, n)
    biomass = 0.05 * np.exp(mu * t)
    reporter = simulate_reporter(
        t, np.full_like(t, k_syn), np.full_like(t, mu),
        ReporterKinetics(k_deg=k_deg), r0=k_syn / (mu + k_deg),
    )
    od = biomass / 0.42 + 0.09
    rfu = np.array([
        OPTICS.background + (OPTICS.autofluorescence + OPTICS.gain * r) * b
        for r, b in zip(reporter, biomass)
    ])
    return t, biomass, reporter, od * (1 + rng.normal(0, 0.01, n)), rfu * (1 + rng.normal(0, 0.01, n))


def _run(t, od, rfu, n_particles=800, seed=1, rfu_mask=None, priors=PRIORS):
    filt = ParticleFilter(priors=priors, optics=OPTICS, gdcw_per_od=0.42, od_blank=0.09,
                          n_particles=n_particles, seed=seed)
    posteriors = []
    for i, time_h in enumerate(t):
        obs = Observation(
            time_h=float(time_h),
            od=float(od[i]),
            rfu=None if (rfu_mask is not None and not rfu_mask[i]) else float(rfu[i]),
        )
        posteriors.append(filt.update(obs))
    return filt, posteriors


def test_it_recovers_the_true_biomass_trajectory():
    t, biomass, _, od, rfu = _truth()

    _, posteriors = _run(t, od, rfu)

    estimated = np.array([p.mean["biomass"] for p in posteriors])
    assert np.median(np.abs(estimated - biomass) / biomass) < 0.10


def test_it_recovers_the_promoter_activity_that_generated_the_signal():
    t, _, _, od, rfu = _truth(k_syn=1.2e-3)

    _, posteriors = _run(t, od, rfu)

    late = np.median([p.mean["promoter_activity"] for p in posteriors[-30:]])
    assert late == pytest.approx(1.2e-3, rel=0.35)


def test_it_recovers_the_growth_rate_jointly_with_everything_else():
    t, _, _, od, rfu = _truth(mu=0.28)

    _, posteriors = _run(t, od, rfu)

    late = np.median([p.mean["growth_rate"] for p in posteriors[-30:]])
    assert late == pytest.approx(0.28, abs=0.06)


def test_the_posterior_narrows_as_observations_accumulate():
    t, _, _, od, rfu = _truth()

    _, posteriors = _run(t, od, rfu)

    assert posteriors[-1].std["biomass"] / posteriors[-1].mean["biomass"] < \
           posteriors[2].std["biomass"] / posteriors[2].mean["biomass"]


class TestMissingModality:
    def test_a_missing_reporter_reading_is_skipped_not_read_as_zero(self):
        """The plan's rule, stated as a test. Zero-filling would imply a dark cell."""
        t, _, _, od, rfu = _truth()
        mask = np.ones(len(t), dtype=bool)
        mask[40:70] = False

        _, gapped = _run(t, od, rfu, rfu_mask=mask)
        _, complete = _run(t, od, rfu)

        gapped_late = gapped[-1].mean["promoter_activity"]
        complete_late = complete[-1].mean["promoter_activity"]
        assert gapped_late == pytest.approx(complete_late, rel=0.6)

    def test_uncertainty_grows_across_a_gap_in_the_reporter_channel(self):
        t, _, _, od, rfu = _truth()
        mask = np.ones(len(t), dtype=bool)
        mask[40:70] = False

        _, posteriors = _run(t, od, rfu, rfu_mask=mask)

        assert posteriors[68].std["promoter_activity"] > posteriors[38].std["promoter_activity"]

    def test_an_update_with_no_channels_at_all_only_propagates(self):
        filt = ParticleFilter(priors=PRIORS, optics=OPTICS, gdcw_per_od=0.42, od_blank=0.09,
                              n_particles=400, seed=0)
        filt.update(Observation(time_h=0.0, od=0.21, rfu=1000.0))
        before = filt.effective_sample_size

        posterior = filt.update(Observation(time_h=0.5, od=None, rfu=None))

        assert posterior.effective_sample_size >= before * 0.99


class TestFirewall:
    def test_an_observation_before_the_last_update_is_refused(self):
        filt = ParticleFilter(priors=PRIORS, optics=OPTICS, gdcw_per_od=0.42, od_blank=0.09,
                              n_particles=200, seed=0)
        filt.update(Observation(time_h=5.0, od=0.3, rfu=2000.0))

        with pytest.raises(ValueError, match="not run backwards"):
            filt.update(Observation(time_h=4.0, od=0.3, rfu=2000.0))

    def test_forecasting_uses_no_observation_at_all(self):
        t, biomass, _, od, rfu = _truth()
        filt, _ = _run(t[:60], od[:60], rfu[:60])

        forecast = filt.forecast(horizon_h=3.0)

        assert forecast.mean["biomass"] > filt.posterior().mean["biomass"]


class TestAForecastDoesNotDisturbTheFilter:
    """A forecast reads the state. Until 2026-08-31 it also spent the random stream.

    ``forecast`` copies the particle array, which says plainly that it was meant to leave
    the filter alone -- but ``_propagate`` draws from the shared generator, and those
    draws were not put back. The method was therefore a hidden coupling between a
    diagnostic and the estimate it was diagnosing, and nothing could see it: no module
    under ``src/`` forecasts mid-run, and ``audit_determinism.py`` calls ``forecast``
    once, last, which is the one ordering under which the fault is invisible.

    These three tests are the orderings that ordering hides.
    """

    def test_it_leaves_the_generator_where_it_found_it(self):
        """The direct statement. Everything below is a consequence of this line."""
        t, _, _, od, rfu = _truth()
        filt, _ = _run(t[:60], od[:60], rfu[:60])
        before = filt._rng.bit_generator.state

        filt.forecast(horizon_h=2.0)

        assert filt._rng.bit_generator.state == before

    def test_forecasting_twice_gives_the_same_answer_twice(self):
        """A method whose contract is 'read the state' cannot answer differently on the
        second reading. This returned two different biomasses."""
        t, _, _, od, rfu = _truth()
        filt, _ = _run(t[:60], od[:60], rfu[:60])

        first = filt.forecast(horizon_h=2.0).mean
        second = filt.forecast(horizon_h=2.0).mean

        assert first == second

    def test_a_forecast_midway_does_not_move_the_posterior_that_follows_it(self):
        """The consequence that would have been read as a result rather than as a bug: two
        runs over identical readings, differing only in whether a diagnostic was called."""
        t, _, _, od, rfu = _truth()
        plain, _ = _run(t[:60], od[:60], rfu[:60])

        # Built by hand rather than through `_run`, because the point is the interleaving
        # -- so it must match `_run`'s construction argument for argument, and then differ
        # only in the forecast call.
        instrumented = ParticleFilter(priors=PRIORS, optics=OPTICS, gdcw_per_od=0.42,
                                      od_blank=0.09, n_particles=800, seed=1)
        for k in range(60):
            instrumented.update(Observation(time_h=float(t[k]), od=float(od[k]),
                                            rfu=float(rfu[k])))
            if k == 30:
                instrumented.forecast(horizon_h=2.0)

        assert instrumented.posterior().mean == plain.posterior().mean


class TestTheReportedESSCanSeeACollapse:
    """The ESS a posterior carries has to be the one that could be small.

    Until 2026-09-08 it could not be. ``update`` reweighted, resampled when the ensemble
    had degenerated, and only then built the ``Posterior`` -- off weights ``_resample``
    had just set back to uniform. So the field read ``n_particles`` at exactly the steps
    it existed to flag, and the test guarding it asserted ``effective_sample_size > 1.0``
    on a quantity that cannot be below ``resample_threshold * n_particles``: a tautology
    standing in for a diagnostic.

    The three tests below are the ones that tautology admitted. Each fails on the
    pre-fix code.
    """

    def test_every_resample_has_a_reported_ess_below_the_threshold_that_caused_it(self):
        """The exact invariant, not an inequality that happens to hold.

        A resample fires when and only when the post-reweight ESS is under the
        threshold, so the count of reported ESS values under the threshold must equal
        the number of resamples. Pre-fix the left side was 0 against 41 resamples.
        """
        t, _, _, od, rfu = _truth()

        filt, posteriors = _run(t, od, rfu)

        floor = filt.resample_threshold * filt.n_particles
        collapsed = [p for p in posteriors if p.effective_sample_size < floor]
        assert filt.n_resamples > 0
        assert len(collapsed) == filt.n_resamples

    def test_the_reported_ess_reaches_a_small_fraction_of_the_ensemble(self):
        """Not just "below the threshold": the ensemble really does collapse.

        On this trace the worst step leaves a few per cent of the particles carrying the
        weight. A diagnostic that reported 800 of 800 there is not a diagnostic.
        """
        t, _, _, od, rfu = _truth()

        filt, posteriors = _run(t, od, rfu)

        worst = min(p.effective_sample_size for p in posteriors)
        assert worst < 0.25 * filt.n_particles
        assert worst >= 1.0

    def test_the_innovation_and_the_posterior_report_the_same_step(self):
        """One number, two carriers. They are read at the same point in the step, so a
        residual and the belief it produced cannot disagree about how much of the
        ensemble was behind them."""
        filt = ParticleFilter(priors=PRIORS, optics=OPTICS, gdcw_per_od=0.42,
                              od_blank=0.09, n_particles=400, seed=3)
        t, _, _, od, rfu = _truth(n=40)

        seen = 0
        for i, time_h in enumerate(t):
            posterior, innovations = filt.update_with_innovations(Observation(
                time_h=float(time_h), od=float(od[i]), rfu=float(rfu[i])))
            for innovation in innovations:
                assert innovation.effective_sample_size == posterior.effective_sample_size
                seen += 1
        assert seen == 2 * len(t)


def test_the_live_effective_sample_size_is_still_what_the_resample_rule_reads():
    """The filter's own property keeps its meaning: the weights as they stand.

    Both quantities are wanted, and the fix was to report the right one rather than to
    redefine the other. Immediately after a resample the live property is n_particles by
    construction -- which is precisely why it is not the one a ``Posterior`` carries.
    """
    filt = ParticleFilter(priors=PRIORS, optics=OPTICS, gdcw_per_od=0.42, od_blank=0.09,
                          n_particles=200, seed=0)

    posterior = filt.update(Observation(time_h=0.0, od=0.21, rfu=1000.0))

    assert filt.n_resamples == 1
    assert filt.effective_sample_size == pytest.approx(200.0)
    assert posterior.effective_sample_size < 0.5 * 200


class TestTheDiagnosticFixChangedNoDynamics:
    """Keep the archived diagnostic-only invariant without freezing a known dynamics bug.

    The literals below describe commit 84561f5 and are not restamped. These two tests
    explicitly install its historical propagator and both legacy noise conventions,
    so diagnostic collection and resampling are checked without freezing the corrected
    observation likelihood, physical-noise or reporter-loss equations. Current dynamics
    and observation models are tested against analytic moments and normalized densities.
    """

    # The final posterior of `_run(*_truth())`, as the pre-fix code produced it.
    PRE_FIX_MEAN = {
        "biomass": "0x1.ec98d50212e39p+2",
        "reporter": "0x1.e44f67ed72997p-9",
        "promoter_activity": "0x1.4643c682ba140p-10",
        "growth_rate": "0x1.2102255374efcp-2",
    }
    PRE_FIX_STD = {
        "biomass": "0x1.bc1f4dea94ecep-4",
        "reporter": "0x1.ee5a4150b26d6p-14",
        "promoter_activity": "0x1.71e61b0be4369p-13",
        "growth_rate": "0x1.4dac544aed4efp-5",
    }
    PRE_FIX_RESAMPLES = 41
    PRE_FIX_BIT_GENERATOR_STATE = 86995385033515949729767634416769399672
    HISTORICAL_PRIORS = replace(PRIORS, positive_noise_mode="median_preserving_legacy",
                                observation_noise_mode="observed_scale_legacy")

    @staticmethod
    def _historical_propagator(filt, x, dt):
        if dt <= 0:
            return x
        biomass, reporter, activity, mu = x.T
        scale = np.sqrt(dt)
        mu_next = mu + filt._rng.normal(0, filt.priors.growth_walk * scale, mu.shape)
        activity_next = activity * np.exp(
            filt._rng.normal(0, filt.priors.activity_walk * scale, activity.shape))
        loss = np.clip(mu_next + filt.priors.k_deg, 1e-6, None)
        steady = activity_next / loss
        reporter_next = steady + (reporter - steady) * np.exp(-loss * dt)
        biomass_next = biomass * np.exp(mu_next * dt) * np.exp(
            filt._rng.normal(0, filt.priors.biomass_walk * scale, biomass.shape))
        return np.column_stack([np.clip(biomass_next, 1e-9, None),
                                np.clip(reporter_next, 0.0, None),
                                np.clip(activity_next, 0.0, None), mu_next])

    def test_the_filtered_posterior_is_bit_for_bit_what_it_was(self, monkeypatch):
        monkeypatch.setattr(ParticleFilter, "_propagate", self._historical_propagator)
        t, _, _, od, rfu = _truth()

        _, posteriors = _run(t, od, rfu, priors=self.HISTORICAL_PRIORS)

        last = posteriors[-1]
        for state, literal in self.PRE_FIX_MEAN.items():
            assert last.mean[state] == float.fromhex(literal)
        for state, literal in self.PRE_FIX_STD.items():
            assert last.std[state] == float.fromhex(literal)

    def test_the_generator_is_spent_in_the_same_order_and_no_more(self, monkeypatch):
        """Diagnostic changes consume no draws under the explicitly frozen dynamics."""
        monkeypatch.setattr(ParticleFilter, "_propagate", self._historical_propagator)
        t, _, _, od, rfu = _truth()

        filt, _ = _run(t, od, rfu, priors=self.HISTORICAL_PRIORS)

        assert filt.n_resamples == self.PRE_FIX_RESAMPLES
        assert (filt._rng.bit_generator.state["state"]["state"]
                == self.PRE_FIX_BIT_GENERATOR_STATE)


def test_more_particles_give_a_more_stable_estimate():
    t, _, _, od, rfu = _truth()

    spread = []
    for n in (150, 2000):
        runs = [_run(t, od, rfu, n_particles=n, seed=s)[1][-1].mean["growth_rate"]
                for s in range(4)]
        spread.append(np.std(runs))

    assert spread[1] < spread[0]


class TestInnovations:
    """The residual is the only thing a filter produces that is checkable without truth.

    So the instrumentation has to be exactly that -- instrumentation. If collecting
    innovations changed the filtering even slightly, every calibration verdict would be
    about a filter nobody runs.
    """

    def _filter(self, seed=7, n=300):
        from ystwin.estimator import ParticleFilter, TwinPriors
        from ystwin.observation import ReporterOptics

        priors = TwinPriors(biomass=(0.05, 0.01), reporter=(1000.0, 200.0),
                            promoter_activity=(300.0, 100.0), growth_rate=(0.3, 0.05))
        optics = ReporterOptics(gain=1.0, background=50.0, autofluorescence=0.0,
                                inner_filter_coeff=None)
        return ParticleFilter(priors, optics, gdcw_per_od=0.3, od_blank=0.0,
                              n_particles=n, seed=seed)

    def _observations(self, n=12, mu=0.3):
        from ystwin.estimator import Observation

        return [Observation(time_h=float(t), od=float(0.05 * np.exp(mu * t)),
                            rfu=float(50 + 1000 * 0.05 * np.exp(mu * t)))
                for t in np.linspace(0.0, 8.0, n)]

    def test_collecting_innovations_does_not_change_the_posterior(self):
        plain, instrumented = self._filter(), self._filter()
        for obs in self._observations():
            a = plain.update(obs)
            b, _ = instrumented.update_with_innovations(obs)
            for state in ("biomass", "reporter", "promoter_activity", "growth_rate"):
                assert a.mean[state] == pytest.approx(b.mean[state], rel=1e-12)
                assert a.std[state] == pytest.approx(b.std[state], rel=1e-12)

    def test_one_innovation_per_measured_channel(self):
        from ystwin.estimator import Observation

        filt = self._filter()
        _, both = filt.update_with_innovations(Observation(0.0, od=0.05, rfu=1000.0))
        assert sorted(i.channel for i in both) == ["od", "rfu"]

        _, od_only = filt.update_with_innovations(Observation(1.0, od=0.06))
        assert [i.channel for i in od_only] == ["od"]

    def test_an_unmeasured_step_yields_no_fabricated_residual(self):
        from ystwin.estimator import Observation

        filt = self._filter()
        _, none = filt.update_with_innovations(Observation(0.0))
        assert none == []

    def test_nis_and_standardised_agree(self):
        filt = self._filter()
        for obs in self._observations(n=6):
            _, innovations = filt.update_with_innovations(obs)
            for i in innovations:
                assert i.nis == pytest.approx(i.standardised**2, rel=1e-9)

    def test_total_variance_sums_its_parts(self):
        from ystwin.estimator import Innovation

        i = Innovation(0.0, "od", 1.0, 0.9, state_variance=0.01,
                       measurement_variance=0.04, effective_sample_size=500.0)
        assert i.total_variance == pytest.approx(0.05)
        assert i.nis == pytest.approx(0.01 / 0.05)

    def test_a_degenerate_variance_is_refused_not_divided_by(self):
        from ystwin.estimator import Innovation

        i = Innovation(0.0, "od", 1.0, 0.9, state_variance=0.0,
                       measurement_variance=0.0, effective_sample_size=1.0)
        assert np.isnan(i.nis)
        assert np.isnan(i.standardised)

    def test_a_well_specified_filter_lands_near_the_expected_nis(self):
        """Truth generated from the filter's own model: NIS should average near one.

        This is the positive control. Without it a uniformly over-confident verdict on
        real data could just as easily mean the statistic is computed wrongly.
        """
        from ystwin.estimator import Observation, ParticleFilter, TwinPriors
        from ystwin.observation import ReporterOptics, observe_od, observe_rfu

        rng = np.random.default_rng(3)
        mu, gain, background, gdcw = 0.30, 1.0, 50.0, 0.3
        priors = TwinPriors(biomass=(0.015, 0.004), reporter=(3000.0, 900.0),
                            promoter_activity=(900.0, 400.0), growth_rate=(mu, 0.06))
        optics = ReporterOptics(gain=gain, background=background,
                                autofluorescence=0.0, inner_filter_coeff=None)
        cal = None
        filt = ParticleFilter(priors, optics, gdcw_per_od=gdcw, od_blank=0.0,
                              n_particles=4000, seed=11, od_calibration=cal)

        collected = []
        biomass, reporter = 0.015, 3000.0
        times = np.linspace(0.0, 10.0, 30)
        for k in range(1, times.size):
            dt = float(times[k] - times[k - 1])
            biomass *= float(np.exp(mu * dt))
            reporter = 900.0 / mu + (reporter - 900.0 / mu) * float(np.exp(-mu * dt))
            od = observe_od(biomass, 0.0, filt.calibration, gdcw)
            rfu = observe_rfu(reporter, biomass, 0.0, optics)
            obs = Observation(time_h=float(times[k]),
                              od=float(od * (1 + rng.normal(0, priors.od_rel_sigma))),
                              rfu=float(rfu * (1 + rng.normal(0, priors.rfu_rel_sigma))))
            _, innovations = filt.update_with_innovations(obs)
            collected += [i.nis for i in innovations if np.isfinite(i.nis)]

        # Generous band: the filter's process noise is not exactly the truth's, so this
        # asserts the statistic is not wildly mis-scaled rather than pinning a value.
        assert 0.01 < float(np.median(collected)) < 20.0
