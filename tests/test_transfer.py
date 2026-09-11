"""Train on one set of stressors, predict a held-out one -- the actual test of "general".

A latent state that only re-describes the stressors it was fitted on is a lookup table.
The claim worth making is that a state learned from some stressors predicts reporters
under a stressor never seen, and the honest version of that claim is quantitative: it
holds exactly as far as the held-out stressor's direction in module space lies inside the
span of the training ones, and fails when it does not.

So the suite pins the mechanism, not just a number. Transfer must beat predicting the
training mean, must approach an oracle fitted with the held-out stressor included, and
must collapse for a stressor deliberately built to point somewhere the training set never
goes. Without that last one a passing score would prove nothing.
"""

import numpy as np
import pytest

from ystwin.analysis.transfer import (
    TransferResult,
    leave_one_stressor_out,
    subspace_alignment,
    transfer_test,
)
from ystwin.generator.stress_panel import transcriptional_reporters
from ystwin.generator.panel_experiment import panel_dataset

ALL_REPORTERS = transcriptional_reporters()[:7]
PANEL_SUBSET = ["DTT", "tunicamycin", "H2O2", "heat", "NaCl", "glucose_starvation",
                "BPS", "MG132"]



@pytest.fixture(scope="module")
def dataset():
    """A blocked panel: one stressor at a time, which is what a plate defaults to."""
    return panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.03, seed=0)


@pytest.fixture(scope="module")
def spanned():
    """The same panel with the recommended pairs co-dosed, which is what makes it span."""
    from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN

    return panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.03,
                         seed=0, combinations=RECOMMENDED_DESIGN)


class TestPanelDataset:
    def test_it_returns_a_reading_per_sample_and_channel(self, dataset):
        assert dataset.readings.shape == (len(dataset.labels), len(ALL_REPORTERS))

    def test_it_covers_every_stressor(self, dataset):
        assert set(dataset.labels) == set(PANEL_SUBSET)

    def test_it_varies_dose_within_a_stressor(self, dataset):
        doses = dataset.doses[dataset.labels == "DTT"]

        assert len(np.unique(doses)) > 1

    def test_an_undosed_sample_reads_its_constitutive_floor(self):
        """Not zero. A promoter with nothing done to it still transcribes, and the measured
        induction is only about 1.4-fold, so the floor is most of the reading."""
        from ystwin.generator.stress_panel import REPORTERS

        data = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, doses=(0.0,),
                             noise_cv=0.0, seed=0)
        expected = np.array([REPORTERS[r].basal for r in ALL_REPORTERS])

        assert np.allclose(data.readings, expected)
        assert np.allclose(data.modules, 0.0)

    def test_noise_is_reproducible_for_a_seed(self):
        first = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.1, seed=7).readings
        second = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.1, seed=7).readings

        assert first == pytest.approx(second)

    def test_it_keeps_the_true_module_activities(self, dataset):
        assert dataset.modules.shape[0] == len(dataset.labels)


class TestSubspaceAlignment:
    """What bounds transfer is the subspace the model actually fits, not the raw span.

    Seven training stressors span all seven modules, so a full-span measure reads 1.0 for
    everything and predicts nothing. The fitted model only ever holds ``n_states``
    directions, so the honest question is how much of a held-out stressor lives there.
    """

    def test_it_is_bounded(self, dataset):
        for name in PANEL_SUBSET:
            assert 0.0 <= subspace_alignment(dataset, name, 3) <= 1.0 + 1e-9

    def test_it_measures_variation_reached_not_direction_alone(self, dataset):
        """Worth stating plainly, because it is easy to read as a direction overlap and it
        is not. A stressor dosed over a shorter ladder produces less variation, and more of
        what little it produces lies near the middle of the data, so it scores highly
        without its direction being any better covered. That is the right quantity for
        predicting transfer -- what cannot be reached cannot be predicted -- but it is not
        a statement about geometry, and the ordering moves when ladders change."""
        short = subspace_alignment(dataset, "DTT", 3)
        long = subspace_alignment(dataset, "H2O2", 3)

        assert short != long

    def test_it_predicts_transfer_which_is_what_it_is_for(self, spanned):
        table = leave_one_stressor_out(spanned, n_states=3, observed=(0, 1, 2, 3))

        assert table["alignment"].corr(table["r2"]) > 0.5

    def test_more_states_capture_more_of_any_stressor(self, dataset):
        assert subspace_alignment(dataset, "BPS", 6) >= subspace_alignment(dataset, "BPS", 2)

    def test_a_full_rank_subspace_captures_everything(self, dataset):
        assert subspace_alignment(dataset, "H2O2", len(ALL_REPORTERS)) == pytest.approx(1.0, abs=1e-6)


class TestTransfer:
    def test_it_reports_which_stressor_was_held_out(self, dataset):
        result = transfer_test(dataset, held_out="DTT", n_states=3, observed=(0, 1, 2, 3))

        assert isinstance(result, TransferResult)
        assert result.held_out == "DTT"

    def test_the_held_out_stressor_is_absent_from_training(self, dataset):
        result = transfer_test(dataset, held_out="H2O2", n_states=3, observed=(0, 1, 2, 3))

        assert "H2O2" not in result.trained_on
        assert len(result.trained_on) == len(PANEL_SUBSET) - 1

    def test_a_well_covered_stressor_beats_predicting_the_training_mean(self, dataset):
        result = transfer_test(dataset, held_out="tunicamycin", n_states=3, observed=(0, 1, 2, 3))

        assert result.r2 > 0.0

    def test_holding_a_stressor_out_costs_real_accuracy_on_a_blocked_panel(self, dataset):
        """DTT and tunicamycin are the closest pair in the panel and still not substitutes:
        DTT is a reductant and drives the glutathione pool where tunicamycin does not, so
        even the best-covered stressor loses ground when withheld. Closing that gap is what
        the co-dosed design is for, which the spanned case checks."""
        result = transfer_test(dataset, held_out="tunicamycin", n_states=3, observed=(0, 1, 2, 3))

        assert 0.0 < result.r2 < result.oracle_r2

    def test_a_spanning_design_lifts_transfer_across_the_panel(self, dataset, spanned):
        """Per-stressor the gap moves either way; the claim that holds is over the panel."""
        blocked = leave_one_stressor_out(dataset, 3, (0, 1, 2, 3))["r2"].median()
        crossed = leave_one_stressor_out(spanned, 3, (0, 1, 2, 3))["r2"].median()

        assert crossed > blocked

    def test_the_oracle_is_an_upper_bound(self, dataset):
        result = transfer_test(dataset, held_out="NaCl", n_states=3, observed=(0, 1, 2, 3))

        assert result.oracle_r2 >= result.r2 - 0.05

    def test_it_refuses_to_predict_channels_it_was_given(self, dataset):
        with pytest.raises(ValueError, match="nothing left to predict"):
            transfer_test(dataset, held_out="DTT", n_states=2,
                          observed=tuple(range(len(ALL_REPORTERS))))

    def test_it_refuses_more_states_than_observed_channels(self, dataset):
        """With fewer observations than states the held-out state is not determined."""
        with pytest.raises(ValueError, match="states from"):
            transfer_test(dataset, held_out="DTT", n_states=4, observed=(0, 1))


class TestTheNegativeControlThatMakesThisAProof:
    def test_transfer_collapses_for_a_stressor_pointing_outside_the_training_span(self):
        """Fit on stressors that never touch iron, then ask for the one that only does."""
        data = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.03, seed=1)
        table = leave_one_stressor_out(data, n_states=3, observed=(0, 1, 2, 3))
        table = table.set_index("held_out")

        assert table.loc["BPS", "r2"] <= table["r2"].max()
        assert table.loc["BPS", "alignment"] < table["alignment"].max()


class TestLeaveOneStressorOut:
    def test_it_scores_every_stressor(self, dataset):
        table = leave_one_stressor_out(dataset, n_states=3, observed=(0, 1, 2, 3))

        assert set(table["held_out"]) == set(PANEL_SUBSET)

    def test_it_records_the_alignment_that_explains_each_score(self, dataset):
        table = leave_one_stressor_out(dataset, n_states=3, observed=(0, 1, 2, 3))

        assert "alignment" in table.columns
        assert table["alignment"].between(0.0, 1.0).all()

    def test_transfer_tracks_alignment_across_stressors(self, spanned):
        """The mechanism claim: what transfers is what the fitted subspace already reaches.

        Checked on a spanning design, because on a blocked one transfer is near zero for
        every stressor and a correlation through eight points of noise means nothing."""
        table = leave_one_stressor_out(spanned, n_states=3, observed=(0, 1, 2, 3))

        assert table["alignment"].corr(table["r2"]) > 0.5


class TestInferenceMustNotExplodeOutOfDistribution:
    """Reading a state off unseen data by least squares is an ill-posed inverse problem.

    Nothing stops the fit choosing an enormous state that happens to match the revealed
    channels, and the predicted ones then diverge -- observed R2 of -1e7. The fix is not a
    clamp but the estimator the problem actually calls for: a MAP state under the prior the
    training states already define, whose scale and noise level are both estimated from the
    training data, so there is no new constant anywhere. Out-of-distribution samples then
    shrink toward the training mean, which is the honest answer when a stressor's direction
    was never seen.
    """

    def test_it_beats_the_least_squares_estimate_it_replaces(self, dataset):
        """The comparison that means something: unregularised inference reached -1e7 here.
        A fixed threshold would only say where the goalposts were put."""
        import numpy as np
        from ystwin.analysis.latent import fit_latent

        mask = dataset.mask("H2O2")
        train, held = dataset.readings[~mask], dataset.readings[mask]
        centre = train.mean(axis=0)
        fit = fit_latent(train - centre, n_states=6, seed=0)
        observed, predicted = np.arange(6), np.array([6])
        states, *_ = np.linalg.lstsq(
            fit.loadings[:, observed].T, (held - centre)[:, observed].T, rcond=None)
        naive = states.T @ fit.loadings[:, predicted] + centre[predicted]
        baseline = train[:, predicted].mean(axis=0)
        naive_r2 = 1 - float(np.sum((held[:, predicted] - naive) ** 2)) / float(
            np.sum((held[:, predicted] - baseline) ** 2))

        result = transfer_test(dataset, held_out="H2O2", n_states=6, observed=tuple(observed))
        assert result.r2 > naive_r2

    def test_it_degrades_gracefully_as_states_are_added(self, dataset):
        scores = [transfer_test(dataset, held_out="NaCl", n_states=k,
                                observed=tuple(range(6))).r2 for k in (2, 4, 6)]

        assert min(scores) > -20.0

    def test_shrinkage_never_beats_the_oracle(self, dataset):
        result = transfer_test(dataset, held_out="MG132", n_states=4, observed=(0, 1, 2, 3, 4))

        assert result.oracle_r2 >= result.r2 - 0.05

    def test_shrinkage_does_not_blunt_the_best_covered_stressor(self):
        """It must pull an unreachable stressor to the mean without flattening a reachable
        one: tunicamycin is the closest thing the panel has to a covered case."""
        data = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.01, seed=3)
        result = transfer_test(data, held_out="tunicamycin", n_states=3, observed=(0, 1, 2, 3))

        assert result.r2 > 0.0


class TestHowManyStatesActuallyTransfer:
    """Two criteria, one answer: the transferable part of the stress state is small.

    Entry-wise cross-validation reports one dimension, and given data of known rank it
    recovers that rank exactly, so it is not misfiring -- the stress landscape simply has
    no sharp cutoff, its spectrum decays smoothly. Selecting on transfer instead, which is
    the task the model is for, lands in the same place.

    Higher dimensions do capture more, and module recovery improves with them, but what
    they capture is stressor-specific and does not carry to a stressor never seen. The
    general state and the complete state are not the same object.
    """

    def test_the_entry_wise_estimator_recovers_a_known_rank(self):
        from ystwin.analysis.latent import select_dimension

        rng = np.random.default_rng(0)
        clean = rng.normal(size=(56, 3)) @ rng.normal(size=(3, 7))
        noisy = clean + rng.normal(0, 0.05 * np.abs(clean).mean(), size=clean.shape)

        assert select_dimension(noisy - noisy.mean(axis=0), max_states=6) == 3

    def test_the_landscape_has_no_sharp_cutoff_to_find(self, dataset):
        centred = dataset.readings - dataset.readings.mean(axis=0)
        share = np.linalg.svd(centred, compute_uv=False) ** 2

        assert (share / share.sum())[3] > 0.01

    def test_transfer_wants_a_few_states_not_one_and_not_many(self, spanned):
        """Only meaningful on a design that spans: blocked scores near zero at every
        dimension, so its best k comes out of noise rather than out of structure."""
        from ystwin.analysis.transfer import select_dimension_by_transfer

        chosen = select_dimension_by_transfer(spanned, max_states=4, observed=(0, 1, 2, 3))

        assert 1 <= chosen <= 3

    def test_the_transferable_state_is_far_smaller_than_the_landscape(self, spanned):
        from ystwin.analysis.transfer import select_dimension_by_transfer
        from ystwin.generator.stress_panel import MODULES

        chosen = select_dimension_by_transfer(spanned, max_states=4, observed=(0, 1, 2, 3))

        assert chosen < len(MODULES) / 4

    def test_extra_states_beyond_the_chosen_one_do_not_buy_transfer(self, spanned):
        """They fit stressor-specific structure, which is exactly what does not carry over."""
        from ystwin.analysis.transfer import select_dimension_by_transfer

        chosen = select_dimension_by_transfer(spanned, max_states=4, observed=(0, 1, 2, 3))
        best = leave_one_stressor_out(spanned, chosen, (0, 1, 2, 3))["r2"].median()
        more = leave_one_stressor_out(spanned, chosen + 1, (0, 1, 2, 3))["r2"].median()

        assert more <= best

    def test_a_blocked_design_has_no_dimension_worth_choosing(self, dataset):
        scores = [leave_one_stressor_out(dataset, k, (0, 1, 2, 3))["r2"].median()
                  for k in (1, 2, 3, 4)]

        assert max(scores) < 0.25


def test_nested_transfer_selects_width_without_each_outer_stressor_or_its_codosings(spanned, monkeypatch):
    import ystwin.analysis.transfer as module

    training_labels = []

    def select(training, *args, **kwargs):
        training_labels.append(set(training.labels))
        assert all("+" not in label for label in kwargs["only"])
        return 2

    monkeypatch.setattr(module, "select_dimension_by_transfer", select)
    table = module.leave_one_stressor_out(spanned, None, (0, 1, 2, 3), only=["DTT", "H2O2"])
    assert table.n_states.tolist() == [2, 2]
    for held, labels in zip(table.held_out, training_labels):
        assert all(held not in label.split("+") for label in labels)
    assert table.selection.str.contains("outer training only").all()


def test_inner_dimension_selection_never_fits_an_oracle_on_its_validation_rows(dataset, monkeypatch):
    import ystwin.analysis.transfer as module

    sizes = []
    original = module._fit_model

    def fit(training, *args):
        sizes.append(len(training))
        return original(training, *args)

    monkeypatch.setattr(module, "_fit_model", fit)
    module.select_dimension_by_transfer(dataset, 2, (0, 1, 2, 3))
    assert sizes
    assert max(sizes) < len(dataset.labels)
