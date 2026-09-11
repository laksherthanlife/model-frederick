"""Making the sweep affordable without changing a single number it reports.

A leave-one-stressor-out sweep fits two latent models per stressor: one on the training
set, which genuinely differs each time, and one oracle on the whole dataset, which does
not. Refitting that identical oracle once per stressor is pure waste, and it is half the
cost of every design evaluation -- which is what makes a search over 300 candidate pairs
unaffordable.

Caching it is exact, so the guard is an equality: the swept table must match what
independent calls produce, entry for entry. A speedup that moves a result is a bug, and
these tests are here to make that impossible to miss.
"""

import pytest

from ystwin.analysis.transfer import leave_one_stressor_out, transfer_test
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import transcriptional_reporters

SUBSET = ["DTT", "H2O2", "heat", "NaCl", "BPS"]
READERS = transcriptional_reporters()[:6]


@pytest.fixture(scope="module")
def dataset():
    return panel_dataset(reporters=READERS, stressors=SUBSET, doses=(0.5, 1.0, 2.0),
                         replicates=2, noise_cv=0.02, seed=0)


class TestTheSweepMatchesIndividualCalls:
    def test_transfer_scores_are_identical(self, dataset):
        swept = leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2))

        for _, row in swept.iterrows():
            one = transfer_test(dataset, row["held_out"], n_states=2, observed=(0, 1, 2))
            assert row["r2"] == pytest.approx(one.r2, rel=1e-12)

    def test_oracle_scores_are_identical(self, dataset):
        swept = leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2))

        for _, row in swept.iterrows():
            one = transfer_test(dataset, row["held_out"], n_states=2, observed=(0, 1, 2))
            assert row["oracle_r2"] == pytest.approx(one.oracle_r2, rel=1e-12)

    def test_alignments_are_identical(self, dataset):
        swept = leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2))

        for _, row in swept.iterrows():
            one = transfer_test(dataset, row["held_out"], n_states=2, observed=(0, 1, 2))
            assert row["alignment"] == pytest.approx(one.alignment, rel=1e-12)

    def test_narrowing_the_sweep_does_not_change_what_is_left(self, dataset):
        everything = leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2))
        narrowed = leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2),
                                          only=["DTT", "BPS"])

        wanted = everything[everything["held_out"].isin(["DTT", "BPS"])]
        assert narrowed["r2"].to_numpy() == pytest.approx(wanted["r2"].to_numpy(), rel=1e-12)

    def test_the_sweep_is_cheaper_than_the_calls_it_replaces(self, dataset):
        """The oracle is one model, however many stressors are held out against it."""
        import ystwin.analysis.transfer as module

        calls = []
        original = module.fit_latent
        module.fit_latent = lambda *a, **k: (calls.append(1), original(*a, **k))[1]
        try:
            leave_one_stressor_out(dataset, n_states=2, observed=(0, 1, 2))
        finally:
            module.fit_latent = original

        assert len(calls) < 2 * len(SUBSET)


class TestParallelSelectionAgrees:
    def test_it_picks_what_the_serial_search_picks(self):
        from ystwin.analysis.experiment_design import select_combinations

        scoring = dict(n_states=2, doses=(0.5, 1.0, 2.0), replicates=2, noise_cv=0.02,
                       seed=0, stressors=SUBSET, reporters=READERS, observed=(0, 1, 2))

        serial = select_combinations(2, library=SUBSET, workers=1, **scoring)
        parallel = select_combinations(2, library=SUBSET, workers=4, **scoring)

        assert serial == parallel
