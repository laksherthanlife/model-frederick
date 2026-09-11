"""A verdict needs an interval, and the earlier one was issued from a single pooled fit.

The first pass pooled every replicate into one curve, fitted it, and refuted the encoded
EC50 because the point estimate was threefold away. That was not sound. A point estimate
from a misspecified model carries no information about whether the encoded value is
excluded, and the pooled fit hid the fact that independent replicates disagreed with each
other by more than twofold.

Fitting each replicate separately gives a spread, and the spread is what decides. The
encoded value is refuted only when it falls outside the range the replicates cover -- which
is a weaker and more honest claim than a fold-change between two numbers.

The shape-free crossing keeps its place but not its old role. It is biased low, badly:
against known parameters it returns anywhere from 0.89 to 0.13 of the truth, because the
peak it measures against is pulled down by the same toxicity that turns the curve over. It
is a lower bound on the EC50, and calling it an estimate was the mistake.
"""

import numpy as np
import pytest

from ystwin.generator.panel_calibration import (
    calibrate_replicates,
    half_maximal_dose,
)


def truth(doses, ec50, lethal, amp=600.0, basal=1000.0, hill=2.5):
    d = np.asarray(doses, dtype=float)
    return (basal + amp * (d / (ec50 + d))) / (1.0 + (d / lethal) ** hill)


LADDER = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])


class TestTheCrossingIsALowerBoundNotAnEstimate:
    def test_it_lands_below_the_ec50_that_generated_it(self):
        for ec50 in (0.2, 0.5, 1.0, 2.0):
            crossing = half_maximal_dose(LADDER, truth(LADDER, ec50, lethal=7.0))
            assert crossing < ec50, ec50

    def test_the_bias_worsens_as_toxicity_bites_earlier(self):
        gentle = half_maximal_dose(LADDER, truth(LADDER, 1.0, lethal=20.0))
        harsh = half_maximal_dose(LADDER, truth(LADDER, 1.0, lethal=3.0))

        assert harsh < gentle

    def test_it_is_nearly_unbiased_when_toxicity_is_far_away(self):
        crossing = half_maximal_dose(LADDER, truth(LADDER, 0.3, lethal=1e6))

        assert crossing == pytest.approx(0.3, rel=0.2)


class TestAVerdictFromReplicates:
    def test_it_returns_one_estimate_per_replicate(self):
        replicates = {f"r{i}": truth(LADDER, 0.4, 5.0) for i in range(3)}
        result = calibrate_replicates(
            "DTT", LADDER, replicates,
            {k: np.full(len(LADDER), 0.3) for k in replicates})

        assert len(result.estimates) == 3

    def test_agreeing_replicates_confirm_a_value_inside_their_range(self):
        rng = np.random.default_rng(0)
        replicates = {f"r{i}": truth(LADDER, 0.4, 5.0) * rng.normal(1, 0.01, len(LADDER))
                      for i in range(4)}
        growth = {k: np.full(len(LADDER), 0.3) for k in replicates}
        result = calibrate_replicates("DTT", LADDER, replicates, growth, encoded=0.4)

        assert result.verdict == "CONFIRMED"

    def test_a_value_far_outside_the_replicate_range_is_refuted(self):
        rng = np.random.default_rng(1)
        replicates = {f"r{i}": truth(LADDER, 0.4, 5.0) * rng.normal(1, 0.01, len(LADDER))
                      for i in range(4)}
        growth = {k: np.full(len(LADDER), 0.3) for k in replicates}
        result = calibrate_replicates("DTT", LADDER, replicates, growth, encoded=20.0)

        assert result.verdict == "REFUTED"

    def test_replicates_that_disagree_with_each_other_cannot_refute_anything(self):
        """Two replicates differing by fourfold say nothing about a third number."""
        replicates = {"a": truth(LADDER, 0.2, 5.0), "b": truth(LADDER, 0.8, 5.0)}
        growth = {k: np.full(len(LADDER), 0.3) for k in replicates}
        result = calibrate_replicates("DTT", LADDER, replicates, growth, encoded=1.0)

        assert result.verdict == "INCONCLUSIVE"

    def test_it_reports_the_interval_it_decided_on(self):
        replicates = {f"r{i}": truth(LADDER, 0.4, 5.0) for i in range(3)}
        growth = {k: np.full(len(LADDER), 0.3) for k in replicates}
        result = calibrate_replicates("DTT", LADDER, replicates, growth, encoded=0.4)

        assert result.low <= result.high
        assert result.low <= np.median(result.estimates) <= result.high

    def test_a_single_replicate_cannot_produce_a_verdict(self):
        replicates = {"only": truth(LADDER, 0.4, 5.0)}
        growth = {"only": np.full(len(LADDER), 0.3)}
        result = calibrate_replicates("DTT", LADDER, replicates, growth, encoded=1.0)

        assert result.verdict == "INCONCLUSIVE"
        assert "one replicate" in result.note
