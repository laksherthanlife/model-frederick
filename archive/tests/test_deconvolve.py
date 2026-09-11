"""Joint deconvolution of N co-expressed reporters, with maturation as information.

Steady state with a maturation step is

    R_i = k_i * m_i / [(mu + m_i + kdeg_i)(mu + kdeg_i)]

Two consequences, and the second is the one that changes the design.

A *matched* pair cancels mu, so its ratio reads relative promoter activity. A *mismatched*
pair does not -- its ratio is (m_A/m_B)(mu+m_B)/(mu+m_A), monotonic in mu from 1 to
m_A/m_B. So a deliberately mismatched pair of constitutive reporters is a growth-rate
sensor, and the wider the mismatch the more sensitive it is.

mOrange2 is therefore not a liability. It is the best mu sensor in the panel.
"""

import numpy as np
import pytest

from ystwin.analysis.deconvolve import (
    growth_from_pair,
    growth_from_reference,
    growth_sensitivity,
    solve_activities,
    steady_state_gain,
)
from ystwin.generator.unmixing import PANEL, FluorophoreSpec

YFP, ORANGE = PANEL["YFP"], PANEL["mOrange2"]


class TestTheGain:
    def test_it_matches_the_closed_form(self):
        """maturation_h is a time; the rate entering the balance is its reciprocal."""
        spec = FluorophoreSpec("x", 500, 520, maturation_h=0.5, k_deg=0.0)
        mu, rate = 0.3, 1 / 0.5
        expected = rate / ((mu + rate) * mu)

        assert steady_state_gain(mu, spec) == pytest.approx(expected)

    def test_slower_maturation_gives_a_lower_gain(self):
        fast = FluorophoreSpec("f", 500, 520, maturation_h=0.2)
        slow = FluorophoreSpec("s", 500, 520, maturation_h=4.0)

        assert steady_state_gain(0.3, fast) > steady_state_gain(0.3, slow)

    def test_the_gain_falls_as_growth_rises(self):
        assert steady_state_gain(0.1, YFP) > steady_state_gain(0.4, YFP)


class TestMismatchIsInformation:
    def test_a_mismatched_pair_ratio_is_monotonic_in_growth(self):
        rates = np.linspace(0.02, 0.6, 30)
        ratios = [steady_state_gain(m, ORANGE) / steady_state_gain(m, YFP) for m in rates]

        assert all(a < b for a, b in zip(ratios, ratios[1:])) or \
               all(a > b for a, b in zip(ratios, ratios[1:]))

    def test_growth_is_recoverable_from_a_mismatched_constitutive_pair(self):
        for truth in (0.05, 0.15, 0.30, 0.45):
            observed = steady_state_gain(truth, ORANGE) / steady_state_gain(truth, YFP)
            assert growth_from_pair(observed, ORANGE, YFP) == pytest.approx(truth, rel=1e-4)

    def test_a_matched_pair_carries_no_growth_information(self):
        """Which is exactly why a matched pair is the right reference for activity."""
        with pytest.raises(ValueError, match="maturation"):
            growth_from_pair(1.0, YFP, YFP)

    def test_a_wider_mismatch_is_a_more_sensitive_growth_sensor(self):
        mild = FluorophoreSpec("mild", 500, 520, maturation_h=0.6)

        assert growth_sensitivity(0.25, ORANGE, YFP) > growth_sensitivity(0.25, mild, YFP)

    def test_the_best_growth_pair_in_the_panel_involves_the_slowest_fluorophore(self):
        names = list(PANEL)
        best = max(
            ((a, b) for a in names for b in names if a != b),
            key=lambda pair: growth_sensitivity(0.25, PANEL[pair[0]], PANEL[pair[1]]),
        )

        assert "mOrange2" in best


class TestTwoRoutesToGrowth:
    """Absolute calibration of one channel, or relative calibration of a mismatched pair."""

    def test_a_calibrated_reference_gives_growth_whatever_its_maturation(self):
        for truth in (0.05, 0.22, 0.45):
            activity = 2.0e-3
            observed = activity * steady_state_gain(truth, YFP)

            assert growth_from_reference(observed, activity, YFP) == pytest.approx(truth, rel=1e-6)

    def test_that_route_works_even_when_every_channel_matches(self):
        """Contrast with the pair route, which needs a mismatch."""
        activity = 2.0e-3
        observed = activity * steady_state_gain(0.3, YFP)

        assert growth_from_reference(observed, activity, YFP) == pytest.approx(0.3, rel=1e-6)

    def test_the_pair_route_needs_only_a_relative_calibration(self):
        truth, k_slow, k_fast = 0.2, 1.0e-3, 2.0e-3
        ratio = ((k_slow * steady_state_gain(truth, ORANGE))
                 / (k_fast * steady_state_gain(truth, YFP)))

        recovered = growth_from_pair(ratio, ORANGE, YFP, relative_activity=k_slow / k_fast)

        assert recovered == pytest.approx(truth, rel=1e-4)


class TestJointSolve:
    def _panel(self):
        return [PANEL["YFP"], PANEL["eGFP"], PANEL["mApple"], PANEL["mOrange2"]]

    def test_activities_are_recovered_once_growth_is_known(self):
        specs = self._panel()
        truth_mu = 0.22
        truth_k = np.array([2.0e-3, 1.2e-3, 1.5e-3, 1.0e-3])
        observed = np.array([k * steady_state_gain(truth_mu, s) for k, s in zip(truth_k, specs)])

        assert solve_activities(observed, specs, truth_mu) == pytest.approx(truth_k, rel=1e-9)

    def test_maturation_is_calibrated_out_rather_than_avoided(self):
        """The answer to 'can we calibrate mOrange2': yes, it is a known parameter.

        Its slow maturation lowers its gain, which the solve divides out exactly. Nothing
        about a slow fluorophore prevents its activity being read.
        """
        specs = [PANEL["YFP"], PANEL["mOrange2"]]
        truth_mu, truth_k = 0.18, np.array([2.0e-3, 3.0e-3])
        observed = np.array([k * steady_state_gain(truth_mu, s) for k, s in zip(truth_k, specs)])

        recovered = solve_activities(observed, specs, truth_mu)

        assert recovered[1] == pytest.approx(3.0e-3, rel=1e-9)

    def test_end_to_end_growth_then_activities(self):
        specs = self._panel()
        truth_mu = 0.22
        truth_k = np.array([2.0e-3, 1.2e-3, 1.5e-3, 1.0e-3])
        observed = np.array([k * steady_state_gain(truth_mu, s) for k, s in zip(truth_k, specs)])

        mu = growth_from_reference(observed[0], truth_k[0], specs[0])
        activities = solve_activities(observed, specs, mu)

        assert mu == pytest.approx(truth_mu, rel=1e-6)
        assert activities == pytest.approx(truth_k, rel=1e-6)

    def test_a_reading_per_channel_is_required(self):
        with pytest.raises(ValueError, match="one reading"):
            solve_activities(np.array([1.0, 2.0]), self._panel(), 0.2)

    def test_noise_degrades_growth_recovery_gracefully(self):
        specs = self._panel()
        rng = np.random.default_rng(0)
        truth_mu, truth_k = 0.22, np.array([2.0e-3, 1.2e-3, 1.5e-3, 1.0e-3])
        errors = []
        for _ in range(40):
            observed = truth_k[0] * steady_state_gain(truth_mu, specs[0])
            noisy = observed * (1 + rng.normal(0, 0.02))
            errors.append(abs(growth_from_reference(noisy, truth_k[0], specs[0]) - truth_mu))

        assert np.median(errors) / truth_mu < 0.15
