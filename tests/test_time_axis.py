"""When to read the plate, which the panel could not express at all.

Every simulated reading was a steady state. A real plate is read on a timetable, and a
promoter fusion does not arrive at steady state on demand: reporter accumulates as
dR/dt = k - (mu + k_deg) R, so it approaches k/(mu + k_deg) with a time constant of
1/(mu + k_deg). Reading before that is reading a fraction of the response.

The consequence is worse than a lag, and it runs the wrong way. The relaxation rate *is*
the growth rate, so a stressed culture -- the one that grows slowly, the one the experiment
is about -- takes longest to reach steady state. An early read therefore understates stress
most in exactly the wells where stress is greatest, which is a bias with dose structure and
not a delay that averages out.

A ratiometric sensor has no such term. It reports an equilibrium between two forms of one
molecule and settles in seconds, so it is at its final value whenever the plate is read.
That is a second, independent reason to carry one, and the time axis is what makes it
visible.
"""

import numpy as np
import pytest

from ystwin.generator.kinetics import reporter_at_time, time_to_fraction
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import REPORTERS


class TestTheReporterBalance:
    def test_it_starts_where_it_started(self):
        assert reporter_at_time(0.0, steady_state=5.0, growth_rate=0.3, initial=1.0) == pytest.approx(1.0)

    def test_it_approaches_the_steady_state(self):
        late = reporter_at_time(100.0, steady_state=5.0, growth_rate=0.3, initial=1.0)

        assert late == pytest.approx(5.0, rel=1e-6)

    def test_it_rises_monotonically_toward_a_higher_steady_state(self):
        values = [reporter_at_time(t, 5.0, 0.3, initial=1.0) for t in (0, 1, 2, 4, 8)]

        assert values == sorted(values)

    def test_it_falls_toward_a_lower_steady_state(self):
        values = [reporter_at_time(t, 1.0, 0.3, initial=5.0) for t in (0, 1, 2, 4, 8)]

        assert values == sorted(values, reverse=True)

    def test_a_slow_culture_takes_longer_to_get_there(self):
        """The bias that runs the wrong way: slow means stressed."""
        fast = reporter_at_time(2.0, 5.0, growth_rate=0.40, initial=1.0)
        slow = reporter_at_time(2.0, 5.0, growth_rate=0.10, initial=1.0)

        assert (fast - 1.0) / 4.0 > (slow - 1.0) / 4.0

    def test_degradation_speeds_it_up(self):
        """An unstable reporter reaches its own steady state sooner, at a lower level."""
        stable = reporter_at_time(2.0, 5.0, 0.1, initial=1.0, k_deg=0.0)
        unstable = reporter_at_time(2.0, 5.0, 0.1, initial=1.0, k_deg=0.5)

        assert unstable != pytest.approx(stable)

    def test_a_non_growing_stable_reporter_never_relaxes(self):
        """Nothing dilutes it and nothing degrades it, so it holds whatever it had."""
        assert reporter_at_time(10.0, 5.0, growth_rate=0.0, initial=1.0, k_deg=0.0) == pytest.approx(1.0)

    def test_it_refuses_a_negative_time(self):
        with pytest.raises(ValueError, match="negative"):
            reporter_at_time(-1.0, 5.0, 0.3, initial=1.0)


class TestHowLongToWait:
    def test_it_reports_the_time_to_reach_a_fraction(self):
        t = time_to_fraction(0.9, growth_rate=0.35)

        assert reporter_at_time(t, 1.0, 0.35, initial=0.0) == pytest.approx(0.9, rel=1e-6)

    def test_a_slower_culture_needs_longer(self):
        assert time_to_fraction(0.9, growth_rate=0.10) > time_to_fraction(0.9, growth_rate=0.40)

    def test_a_non_growing_stable_culture_never_arrives(self):
        assert not np.isfinite(time_to_fraction(0.9, growth_rate=0.0))

    def test_it_refuses_a_fraction_outside_the_approach(self):
        with pytest.raises(ValueError, match="between"):
            time_to_fraction(1.5, growth_rate=0.3)


class TestTheDatasetCanBeReadOnATimetable:
    def test_it_can_produce_several_read_times(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, read_times_h=(2.0, 6.0, 24.0))

        assert len(set(data.read_times)) == 3

    def test_an_early_read_understates_the_response(self):
        early = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                              read_times_h=(1.0,))
        late = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                             read_times_h=(48.0,))
        floor = REPORTERS["UPRE-ER"].basal

        assert (early.readings.max() - floor) < (late.readings.max() - floor)

    def test_a_long_read_reproduces_the_steady_state_the_panel_used_to_assume(self):
        timed = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                              read_times_h=(500.0,))
        steady = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0)

        assert timed.readings == pytest.approx(steady.readings, rel=1e-3)

    def test_a_ratiometric_sensor_is_already_there_at_any_read_time(self):
        early = panel_dataset(reporters=["roGFP2-Grx1"], stressors=["H2O2"], noise_cv=0.0,
                              read_times_h=(0.5,))
        late = panel_dataset(reporters=["roGFP2-Grx1"], stressors=["H2O2"], noise_cv=0.0,
                             read_times_h=(48.0,))

        assert early.readings == pytest.approx(late.readings, rel=1e-6)

    def test_the_understatement_is_worst_where_the_stress_is_worst(self):
        """A dose that slows growth also slows the reporter reaching its new level, so the
        early read loses most exactly where the response is largest."""
        early = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                              read_times_h=(1.0,), replicates=1)
        late = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                             read_times_h=(500.0,), replicates=1)
        floor = REPORTERS["UPRE-ER"].basal
        captured = (early.readings[:, 0] - floor) / np.maximum(late.readings[:, 0] - floor, 1e-9)

        assert captured[-1] < captured[0]

    def test_the_true_module_activities_do_not_depend_on_when_it_was_read(self):
        early = panel_dataset(stressors=["DTT"], noise_cv=0.0, read_times_h=(1.0,))
        late = panel_dataset(stressors=["DTT"], noise_cv=0.0, read_times_h=(48.0,))

        assert early.modules == pytest.approx(late.modules)

    def test_no_read_time_still_means_steady_state(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0)

        assert data.read_times == [] or len(set(data.read_times)) == 1
