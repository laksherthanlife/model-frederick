"""Putting the growth confound into the panel generator, which did not have one.

`culture.py` says why it matters in its own docstring: growth sits in the reporter's
balance, so a culture whose promoter never moves still shows an apparent dose response, and
data without that term validates the pipeline against a strawman. The panel generator was
that strawman -- readings were loadings times activity times a noise draw, with no growth,
no dilution and no maturation anywhere in it.

The confound is derivable rather than invented. Activity is recovered as k = dR/dt + mu*R,
so an error in the estimated growth rate propagates as k_hat = k * (1 + eps / (mu + k_deg)).
Growth falls with dose, so the same absolute error in mu becomes a larger relative error in
activity at every rung up the ladder -- which is why this is a bias with dose structure and
not another noise term. The plates put the growth-rate CV at 0.33, so the size of it is
measured too.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import STRESSORS, growth_rate


class TestGrowthRate:
    def test_an_undosed_culture_grows_at_its_maximum(self):
        assert growth_rate("DTT", 0.0) == pytest.approx(growth_rate("DTT", 0.0))
        assert growth_rate("DTT", 0.0) > growth_rate("DTT", 1.0)

    def test_it_halves_at_the_measured_lethal_dose(self):
        """That is what the lethal dose was measured as: where growth falls to half."""
        for name in ("DTT", "H2O2"):
            lethal = STRESSORS[name].lethal_dose
            assert growth_rate(name, lethal) == pytest.approx(0.5 * growth_rate(name, 0.0), rel=0.05)

    def test_it_never_goes_negative(self):
        assert growth_rate("H2O2", 1e6) >= 0.0

    def test_an_unknown_stressor_is_refused(self):
        with pytest.raises(KeyError, match="unobtainium"):
            growth_rate("unobtainium", 1.0)


class TestTheConfoundHasDoseStructure:
    def test_a_dataset_can_be_generated_with_it(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117, seed=0)

        assert data.readings.shape[0] > 0

    def test_it_does_nothing_when_growth_is_estimated_perfectly(self):
        clean = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0, seed=0)
        also = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0, seed=1)

        assert clean.readings == pytest.approx(also.readings)

    def test_it_perturbs_the_readings_when_growth_is_uncertain(self):
        clean = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0, seed=0)
        confounded = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117, seed=0)

        assert not np.allclose(clean.readings, confounded.readings)

    def test_the_error_grows_with_dose_which_plain_noise_would_not(self):
        """The whole point: the same absolute error in mu bites harder as growth falls."""
        errors = []
        for seed in range(40):
            clean = panel_dataset(stressors=["H2O2"], noise_cv=0.0, growth_rate_se=0.0, replicates=1)
            dirty = panel_dataset(stressors=["H2O2"], noise_cv=0.0, growth_rate_se=0.0117,
                                  replicates=1, seed=seed)
            relative = np.abs(dirty.readings - clean.readings) / np.maximum(clean.readings, 1e-9)
            errors.append(relative.mean(axis=1))
        mean_by_dose = np.mean(errors, axis=0)

        assert mean_by_dose[-1] > mean_by_dose[0]

    def test_it_is_reproducible_for_a_seed(self):
        first = panel_dataset(stressors=["DTT"], growth_rate_se=0.0117, seed=3).readings
        second = panel_dataset(stressors=["DTT"], growth_rate_se=0.0117, seed=3).readings

        assert first == pytest.approx(second)

    def test_the_true_module_activities_are_left_untouched(self):
        """The confound is a measurement error, not a change in the biology."""
        clean = panel_dataset(stressors=["DTT"], growth_rate_se=0.0, seed=0)
        dirty = panel_dataset(stressors=["DTT"], growth_rate_se=0.0117, seed=0)

        assert dirty.modules == pytest.approx(clean.modules)
