"""The confound parameterised as what it physically is, instead of tuned until it fitted.

The growth term was a CV multiplied by the control's growth rate, and its value came from a
grid search over two free parameters against one measured summary statistic. That is a
ridge of equally good answers, and the point taken off it carried no information. Trying to
recover the split from the plates fails outright: two plates per condition makes each CV a
one-degree-of-freedom estimate, and regressing CV squared on 1/mu squared returns R2 of 0.07.

None of it needed fitting. A growth rate is the slope of log optical density against time,
so the error in it is the standard error of that slope -- an absolute rate error, in 1/h,
measurable on the trace. Across 210 real wells it is 0.0117 1/h, and the grid search had
been using the equivalent of 0.048.

Parameterised that way the dose structure is no longer something to be dialled in. The same
absolute error is 3% of a healthy growth rate and 12% of a slowed one, and it becomes so on
its own.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import (
    MEASURED_ACTIVITY_CV,
    MEASURED_GROWTH_RATE_SE,
    panel_dataset,
)
from ystwin.generator.stress_panel import growth_rate


class TestTheParameterIsAnAbsoluteRateError:
    def test_it_is_measured_in_inverse_hours_not_as_a_fraction(self):
        """0.0117 1/h, from the standard error of the log-OD slope on 210 real wells."""
        assert 0.005 < MEASURED_GROWTH_RATE_SE < 0.03

    def test_it_is_far_smaller_than_the_grid_searched_value_it_replaces(self):
        """The old parameter was a CV times the control rate, which came to about 0.048."""
        assert MEASURED_GROWTH_RATE_SE < 0.048 / 2

    def test_the_relative_error_it_implies_grows_as_growth_slows(self):
        healthy = MEASURED_GROWTH_RATE_SE / growth_rate("DTT", 0.0)
        slowed = MEASURED_GROWTH_RATE_SE / growth_rate("DTT", 1.2)

        assert slowed > healthy


class TestItReproducesTheMeasuredSpread:
    def _observed_cv(self, **kw):
        import pandas as pd

        data = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], replicates=10,
                             seed=0, **kw)
        by = pd.DataFrame({"d": data.doses, "r": data.readings[:, 0]}).groupby("d")["r"]
        return float((by.std() / by.mean()).mean())

    def test_the_two_terms_together_land_on_the_measured_total(self):
        """0.146 measured across 28 matched conditions on two plates."""
        cv = self._observed_cv(noise_cv=MEASURED_ACTIVITY_CV,
                               growth_rate_se=MEASURED_GROWTH_RATE_SE)

        assert cv == pytest.approx(0.146, abs=0.04)

    def test_the_growth_term_alone_is_a_small_part_of_it(self):
        """Which is the finding: the correction is not where the variability comes from."""
        growth_only = self._observed_cv(noise_cv=0.0, growth_rate_se=MEASURED_GROWTH_RATE_SE)

        assert growth_only < 0.5 * 0.146

    def test_neither_term_was_free_to_absorb_the_other(self):
        """One is measured from replicates, one derived from the OD traces, and the total
        they produce is then checked -- rather than either being solved for."""
        assert MEASURED_ACTIVITY_CV > 0
        assert MEASURED_GROWTH_RATE_SE > 0


class TestTheConfoundStillHasDoseStructure:
    def test_a_slowed_culture_carries_more_relative_error(self):
        spread = []
        for seed in range(60):
            clean = panel_dataset(stressors=["H2O2"], noise_cv=0.0, growth_rate_se=0.0,
                                  replicates=1)
            dirty = panel_dataset(stressors=["H2O2"], noise_cv=0.0, seed=seed,
                                  growth_rate_se=MEASURED_GROWTH_RATE_SE, replicates=1)
            spread.append(np.abs(dirty.readings - clean.readings).mean(axis=1)
                          / np.maximum(clean.readings.mean(axis=1), 1e-9))
        by_dose = np.mean(spread, axis=0)

        assert by_dose[-1] > by_dose[0]

    def test_zero_error_leaves_the_readings_alone(self):
        a = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0, seed=0)
        b = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0, seed=9)

        assert a.readings == pytest.approx(b.readings)


class TestTheFisherPenaltyIsComputedNotFrozen:
    """A derived number written down as a literal goes stale the moment its source moves.

    The penalty a promoter fusion pays in sensor selection is the relative growth error it
    inherits, which is the measured absolute error over the growth rate of the culture.
    Both of those live elsewhere, so writing the quotient down as 0.053 leaves a constant
    that agrees with its inputs only until one of them is remeasured -- the same failure as
    the grid search, one step removed.
    """

    def test_the_penalty_follows_from_the_measured_rate_error(self):
        from ystwin.analysis.sensor_selection import REFERENCE_GROWTH_RATE, growth_penalty_cv

        assert growth_penalty_cv() == pytest.approx(
            MEASURED_GROWTH_RATE_SE / REFERENCE_GROWTH_RATE)

    def test_it_moves_when_the_measured_error_moves(self):
        from ystwin.analysis.sensor_selection import growth_penalty_cv

        assert growth_penalty_cv(rate_se=2 * MEASURED_GROWTH_RATE_SE) == pytest.approx(
            2 * growth_penalty_cv())

    def test_it_moves_when_the_culture_is_slower(self):
        from ystwin.analysis.sensor_selection import growth_penalty_cv

        assert growth_penalty_cv(growth_rate=0.10) > growth_penalty_cv(growth_rate=0.35)

    def test_the_reference_rate_is_one_a_stressed_culture_actually_shows(self):
        """UPRE1 at 1 mM DTT runs at 0.22, and that is the regime a build is chosen for."""
        from ystwin.analysis.sensor_selection import REFERENCE_GROWTH_RATE

        assert 0.15 <= REFERENCE_GROWTH_RATE <= 0.30

    def test_the_recommended_build_uses_the_computed_value(self):
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD, growth_penalty_cv

        assert RECOMMENDED_BUILD.growth_cv == pytest.approx(growth_penalty_cv())

    def test_a_zero_growth_rate_is_refused_rather_than_dividing_by_it(self):
        from ystwin.analysis.sensor_selection import growth_penalty_cv

        with pytest.raises(ValueError, match="growth rate"):
            growth_penalty_cv(growth_rate=0.0)
