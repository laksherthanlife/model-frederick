"""The constitutive floor every real promoter sits on, which the generator did not have.

Simulated readings rose from zero, so a fully induced module looked like an infinite fold
change. Real promoters are not silent: the calibration fit puts UPRE's basal at 1847 against
an induced amplitude of 493, and the measured induction across the whole DTT ladder is only
1.38 to 1.47 fold. Starting from zero inflates the contrast by roughly an order of magnitude,
and every identifiability, transfer and recovery number was computed at that inflated
contrast.

The posterior-predictive gate is what caught it. Simulated spread came out at 1.5 to 1.8
against 0.17 measured, and no amount of adding noise closes that -- the difference is not
noise, it is a missing floor. This is the gate doing the only job it has: comparing the
generator against data rather than against itself.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.panel_gates import resembles_reference
from ystwin.generator.stress_panel import REPORTERS


class TestReportersHaveABasal:
    def test_every_reporter_declares_one(self):
        assert all(r.basal > 0 for r in REPORTERS.values())

    def test_the_floor_is_larger_than_what_induction_adds_to_it(self):
        """Measured: basal 1847 against an induced 713 at peak on the UPRE ladder, so the
        floor dominates. Compared here in generator units against the actual induced
        reading, not against an RFU amplitude that lives on a different scale."""
        import numpy as np
        from ystwin.generator.stress_panel import (
            MODULES, healthy_ladder, module_response, reporter_loadings)

        row = reporter_loadings(["UPRE-ER"])[0]
        induced = max(float(row @ np.array([module_response("DTT", d)[m] for m in MODULES]))
                      for d in healthy_ladder("DTT"))

        assert REPORTERS["UPRE-ER"].basal > induced


class TestTheFloorAppearsInTheReadings:
    def test_an_undosed_well_reads_its_basal_not_zero(self):
        data = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], doses=(0.0,),
                             noise_cv=0.0, replicates=1)

        assert data.readings[0, 0] == pytest.approx(REPORTERS["UPRE-ER"].basal)

    def test_induction_is_a_modest_fold_not_an_infinite_one(self):
        data = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                             replicates=1)
        fold = data.readings.max() / data.readings.min()

        assert 1.1 < fold < 3.0

    def test_the_measured_fold_is_reproduced_within_reason(self):
        """UPRE1 and UPRE2 peak at 1.38 and 1.47 fold over their own zero-dose wells."""
        data = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], noise_cv=0.0,
                             replicates=1)
        fold = data.readings.max() / data.readings.min()

        assert fold == pytest.approx(1.42, abs=0.5)

    def test_a_ratiometric_sensor_keeps_its_own_floor(self):
        """A pool sensor reports displacement, so its floor is its set point."""
        data = panel_dataset(reporters=["roGFP2-Grx1"], stressors=["H2O2"], doses=(0.0,),
                             noise_cv=0.0, replicates=1)

        assert data.readings[0, 0] == pytest.approx(REPORTERS["roGFP2-Grx1"].basal)

    def test_the_true_module_activities_are_untouched(self):
        """The floor is expression, not biology: the modules are still what they were."""
        data = panel_dataset(reporters=["UPRE-ER"], stressors=["DTT"], doses=(0.0,),
                             noise_cv=0.0, replicates=1)

        assert np.allclose(data.modules, 0.0)


class TestItNowResemblesTheRealPlates:
    def test_the_simulated_spread_lands_near_the_measured_one(self):
        import pandas as pd

        real = pd.read_csv("outputs/sensor_characterisation.csv")
        reference = real[real.construct.isin(["UPRE1", "UPRE2"])].pivot_table(
            index=["plate", "dose_mM"], columns="construct",
            values="activity_late").dropna().to_numpy()
        sim = panel_dataset(reporters=["UPRE-ER", "STRE-general"], stressors=["DTT"],
                            replicates=7, noise_cv=0.146, seed=0)

        assert resembles_reference(sim.readings, reference, tolerance=1.0)
