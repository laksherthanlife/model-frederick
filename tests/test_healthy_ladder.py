"""Dosing a simulated plate the way the real one should have been dosed.

The default ladder ran to four times each agent's EC50, which for DTT is well past the
1.55 mM where growth halves. That is the same mistake the uploaded plates made: doses were
spent above the point where the culture stops reporting anything, and the ladder never
resolved the induction it was there to measure.

Spanning up to a fraction of the measured lethal dose instead keeps every well in the range
where activity can be read, and spends the wells that were being wasted on resolving the
rising limb. Where an agent has a wide window the ladder still reaches saturation; where it
has a narrow one, the ladder says so by stopping early rather than by returning dead wells.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import STRESSORS, healthy_ladder, viability


class TestTheLadder:
    def test_every_dose_leaves_the_culture_healthy(self):
        for name in STRESSORS:
            for dose in healthy_ladder(name):
                assert viability(name, dose) > 0.6, (name, dose)

    def test_it_reaches_at_least_the_inducing_dose_where_the_window_allows(self):
        """Peroxide halves growth at twice its EC50, so the ladder should still get there."""
        assert max(healthy_ladder("H2O2")) >= STRESSORS["H2O2"].ec50 * 0.5

    def test_a_narrow_window_gives_a_shorter_reach(self):
        """DTT is lethal at 1.55 times its inducing dose; BPS has room to spare."""
        dtt = max(healthy_ladder("DTT")) / STRESSORS["DTT"].ec50
        bps = max(healthy_ladder("BPS")) / STRESSORS["BPS"].ec50

        assert dtt < bps

    def test_doses_are_distinct_and_ordered(self):
        rungs = healthy_ladder("DTT")

        assert list(rungs) == sorted(rungs)
        assert len(set(np.round(rungs, 9))) == len(rungs)

    def test_it_gives_the_requested_number_of_rungs(self):
        assert len(healthy_ladder("DTT", n_doses=6)) == 6

    def test_it_refuses_an_unknown_stressor(self):
        with pytest.raises(KeyError, match="unobtainium"):
            healthy_ladder("unobtainium")


class TestTheSimulatedPlateUsesIt:
    def test_no_simulated_well_is_dosed_past_health(self):
        data = panel_dataset(stressors=["DTT", "H2O2", "BPS"], noise_cv=0.0)

        for name in ("DTT", "H2O2", "BPS"):
            for dose in np.unique(data.doses[data.labels == name]):
                assert viability(name, float(dose)) > 0.6, (name, dose)

    def test_activity_rises_across_the_whole_ladder(self):
        """Nothing in a healthy ladder should be on the falling limb."""
        data = panel_dataset(stressors=["BPS"], noise_cv=0.0, replicates=1)
        by_dose = {}
        for dose, row in zip(data.doses, data.readings):
            by_dose[float(dose)] = row.sum()
        rising = [by_dose[d] for d in sorted(by_dose)]

        assert rising == sorted(rising)

    def test_an_explicit_dose_list_is_still_honoured(self):
        data = panel_dataset(stressors=["DTT"], doses=(0.5, 1.0), noise_cv=0.0, replicates=1)

        assert len(np.unique(data.doses)) == 2
