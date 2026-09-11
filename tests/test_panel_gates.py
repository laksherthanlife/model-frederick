"""Acceptance gates on the panel path, which had none.

Simulated data was never checked against anything before being trained on. Three gates
now stand between a generator and a model:

physical -- a reading that could not have come off a plate. Activity below zero is not a
weak promoter, it is a blank subtraction gone wrong; a reading past the reader's ceiling is
not a strong one.

posterior-predictive -- summary statistics outside the spread of a real plate. This is the
gate worth having most, because it is the only one that compares simulated data against
measured data rather than against its own assumptions. A generator can be internally
consistent and still produce activity distributions no yeast ever made.

diversity -- a batch whose members are near-identical, which passes every other check while
carrying the information of a single plate.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.panel_gates import (
    GateReport,
    accept_panel,
    diversity_ok,
    physically_possible,
    resembles_reference,
)


@pytest.fixture(scope="module")
def dataset():
    from ystwin.generator.stress_panel import transcriptional_reporters

    return panel_dataset(reporters=transcriptional_reporters()[:6],
                         stressors=["DTT", "H2O2", "heat"], noise_cv=0.05, seed=0)


class TestPhysical:
    def test_a_clean_dataset_passes(self, dataset):
        assert physically_possible(dataset.readings) is None

    def test_a_non_finite_reading_is_refused(self, dataset):
        broken = dataset.readings.copy()
        broken[0, 0] = np.nan

        assert physically_possible(broken) is not None

    def test_negative_activity_is_refused(self, dataset):
        broken = dataset.readings.copy()
        broken[0, 0] = -5.0

        assert "below zero" in physically_possible(broken)

    def test_a_reading_past_the_ceiling_is_refused(self, dataset):
        broken = dataset.readings.copy()
        broken[0, 0] = 1e9

        assert "ceiling" in physically_possible(broken, ceiling=1e3)

    def test_a_dataset_with_no_induction_anywhere_is_refused(self):
        assert "induction" in physically_possible(np.zeros((10, 4)))

    def test_a_reductant_stays_above_the_floor_it_sits_on(self):
        """DTT drives the glutathione pool the reducing way, but roGFP2 reports a ratio
        around a set point it never reaches, so the reading falls without going negative."""
        data = panel_dataset(reporters=["STRE-general", "roGFP2-Grx1"], stressors=["DTT"],
                             noise_cv=0.0, seed=0)
        redox = data.readings[:, 1]

        assert redox.min() > 0
        assert redox.min() < redox.max()
        assert physically_possible(data.readings, data.reporters) is None


class TestPosteriorPredictive:
    def test_a_dataset_like_its_reference_passes(self, dataset):
        assert resembles_reference(dataset.readings, dataset.readings)

    def test_a_rescaled_dataset_passes_because_gain_is_arbitrary(self):
        """Simulated readings are module activity, not RFU. Holding them to the reference's
        absolute scale would refuse correct data for using different units."""
        from ystwin.generator.panel_experiment import panel_dataset
        from ystwin.generator.stress_panel import transcriptional_reporters

        data = panel_dataset(reporters=transcriptional_reporters()[:6],
                             stressors=["DTT", "H2O2"], noise_cv=0.05, seed=0)

        assert resembles_reference(data.readings * 100.0, data.readings)

    def test_a_dataset_with_the_wrong_dynamic_range_fails(self):
        """What it must catch: a spread no plate produces, whatever the units."""
        rng = np.random.default_rng(0)
        reference = rng.lognormal(0.0, 0.3, size=(60, 5))
        too_wide = rng.lognormal(0.0, 3.0, size=(60, 5))

        assert not resembles_reference(too_wide, reference)

    def test_a_dataset_with_no_spread_fails(self, dataset):
        flat = np.full_like(dataset.readings, float(dataset.readings.mean()))

        assert not resembles_reference(flat, dataset.readings)

    def test_a_tighter_tolerance_refuses_more(self):
        rng = np.random.default_rng(1)
        reference = rng.lognormal(0.0, 0.3, size=(60, 5))
        somewhat_wider = rng.lognormal(0.0, 0.55, size=(60, 5))

        assert resembles_reference(somewhat_wider, reference, tolerance=2.0)
        assert not resembles_reference(somewhat_wider, reference, tolerance=0.05)

    def test_it_compares_shape_not_absolute_units(self, dataset):
        """Reporter gain is arbitrary, so a rescaled reference must still be recognisable."""
        assert resembles_reference(dataset.readings * 3.0, dataset.readings * 3.0)


class TestDiversity:
    def test_a_varied_batch_passes(self):
        rng = np.random.default_rng(0)
        batch = [rng.normal(1.0, 0.2, size=(20, 4)) for _ in range(4)]

        assert diversity_ok(batch)

    def test_a_batch_of_copies_fails(self, dataset):
        assert not diversity_ok([dataset.readings] * 4)

    def test_a_single_member_batch_cannot_be_judged(self, dataset):
        assert not diversity_ok([dataset.readings])


class TestTheGatesTogether:
    def test_a_good_dataset_is_accepted(self, dataset):
        report = accept_panel(dataset.readings, reference=dataset.readings)

        assert report.accepted
        assert isinstance(report, GateReport)

    def test_the_report_names_the_gate_that_refused(self, dataset):
        broken = dataset.readings.copy()
        broken[0, 0] = -1.0
        report = accept_panel(broken, reference=dataset.readings)

        assert not report.accepted
        assert report.refused_by == "physical"

    def test_a_missing_reference_skips_that_gate_and_says_so(self, dataset):
        report = accept_panel(dataset.readings)

        assert report.accepted
        assert "posterior-predictive" in report.skipped

    def test_a_batch_of_copies_is_refused_for_diversity(self, dataset):
        report = accept_panel(dataset.readings, batch=[dataset.readings] * 3)

        assert not report.accepted
        assert report.refused_by == "diversity"

    def test_every_gate_run_is_recorded_whether_it_passed_or_not(self, dataset):
        report = accept_panel(dataset.readings, reference=dataset.readings)

        assert set(report.checked) == {"physical", "posterior-predictive"}
