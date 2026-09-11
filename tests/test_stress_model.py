"""The trained artefact: readings in, module activities out.

Everything up to here measures whether a stress state can be learned. This is the thing
that gets used -- a fitted model that takes a plate's reporter readings and returns an
estimate of what each stress module was doing, with the channels it expects named so it
cannot be handed the wrong columns silently.

Two pieces are learned and they are learned differently. The latent basis comes from the
readings alone, unsupervised. The map from that basis to named modules is supervised, and
in simulation the truth is available to supervise it. That is exactly why simulation earns
its place: it is the only setting where the readout can be fitted at all, and once fitted
it applies to real readings that carry no labels.
"""

import numpy as np
import pytest

from ystwin.analysis.stress_model import StressModel, train_stress_model
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import MODULES, STRESSORS, transcriptional_reporters

READERS = transcriptional_reporters()[:7]
SUBSET = ["DTT", "H2O2", "heat", "NaCl", "glucose_starvation", "BPS", "MG132"]


@pytest.fixture(scope="module")
def dataset():
    return panel_dataset(reporters=READERS, stressors=SUBSET, doses=(0.25, 0.5, 1.0, 2.0),
                         replicates=3, noise_cv=0.02, seed=0)


@pytest.fixture(scope="module")
def model(dataset):
    return train_stress_model(dataset, n_states=3)


class TestTraining:
    def test_it_returns_a_model(self, model):
        assert isinstance(model, StressModel)

    def test_it_remembers_which_channels_it_expects(self, model):
        assert model.reporters == READERS

    def test_it_names_every_module_it_reports(self, model):
        assert model.modules == list(MODULES)

    def test_the_latent_basis_has_the_requested_width(self, model):
        assert model.loadings.shape == (3, len(READERS))

    def test_it_records_how_well_each_module_was_recovered(self, model):
        assert set(model.recovery) == set(MODULES)
        assert all(0.0 <= v <= 1.0 for v in model.recovery.values())

    def test_it_refuses_more_states_than_channels(self, dataset):
        with pytest.raises(ValueError, match="cannot exceed channels"):
            train_stress_model(dataset, n_states=len(READERS) + 1)


class TestInference:
    def test_it_returns_one_state_per_reading(self, model, dataset):
        states = model.infer(dataset.readings)

        assert states.shape == (len(dataset.readings), 3)

    def test_it_accepts_a_single_reading(self, model, dataset):
        assert model.infer(dataset.readings[0]).shape == (1, 3)

    def test_it_refuses_the_wrong_number_of_channels(self, model, dataset):
        with pytest.raises(ValueError, match="expects 7 channels"):
            model.infer(dataset.readings[:, :4])

    def test_it_recovers_the_modules_it_was_trained_on(self, model, dataset):
        predicted = model.predict_modules(dataset.readings)

        assert predicted["ESR"].corr(dataset.modules_frame()["ESR"]) > 0.5

    def test_an_unstressed_reading_reports_little_activity(self, model):
        quiet = panel_dataset(reporters=READERS, stressors=SUBSET, doses=(0.0,),
                              noise_cv=0.0, seed=0)
        predicted = model.predict_modules(quiet.readings)

        assert abs(predicted["ESR"]).max() < 0.3

    def test_the_prior_pulls_every_state_below_the_least_squares_answer(self, model, dataset):
        """MAP is linear, so it caps nothing outright -- what it does is shrink, and the
        comparison that shows it is against the unregularised fit on the same readings."""
        extreme = dataset.readings * 50.0
        mapped = model.infer(extreme)
        least_squares, *_ = np.linalg.lstsq(
            model.loadings.T, (extreme - model.centre).T, rcond=None)

        assert np.abs(mapped).max() < np.abs(least_squares.T).max()


class TestPersistence:
    def test_it_round_trips_through_a_file(self, model, dataset, tmp_path):
        path = tmp_path / "model.npz"
        model.save(path)
        loaded = StressModel.load(path)

        assert loaded.reporters == model.reporters
        assert loaded.infer(dataset.readings) == pytest.approx(model.infer(dataset.readings))

    def test_a_loaded_model_predicts_identically(self, model, dataset, tmp_path):
        path = tmp_path / "model.npz"
        model.save(path)
        loaded = StressModel.load(path)

        assert loaded.predict_modules(dataset.readings).to_numpy() == pytest.approx(
            model.predict_modules(dataset.readings).to_numpy())

    def test_it_keeps_the_recovery_scores(self, model, tmp_path):
        path = tmp_path / "model.npz"
        model.save(path)

        assert StressModel.load(path).recovery == pytest.approx(model.recovery)


class TestPredictingModulesForAStressorNeverSeen:
    """The end-to-end claim, in the units anyone actually cares about.

    Channel-level transfer says the model can fill in a reading. This says something
    stronger and more useful: trained without a stressor, does it report the right
    biology when that stressor arrives? The readout was fitted on other stressors'
    labels, so nothing about the held-out one entered either stage.

    A module the training stressors never excited cannot come back, and the suite says so
    rather than averaging that failure away.
    """

    def test_it_scores_every_module(self, dataset):
        from ystwin.analysis.stress_model import module_transfer

        scored = module_transfer(dataset, held_out="H2O2", n_states=3)

        assert set(scored) == set(MODULES)

    def test_a_module_every_training_stressor_drives_is_among_the_best_recovered(self, dataset):
        """Every stressor in the panel raises the ESR, so it is the axis with most support
        and should rank near the top -- a sharper claim than a threshold picked by hand."""
        from ystwin.analysis.stress_model import module_transfer

        scored = module_transfer(dataset, held_out="H2O2", n_states=3)
        ranked = sorted(scored, key=scored.get, reverse=True)

        assert "ESR" in ranked[:3]

    def test_only_the_target_with_peer_support_comes_back_at_all(self, dataset):
        """Sharper than a threshold, and it is the mechanism. H2O2 drives four modules
        here. Six of the other stressors raise the ESR, one raises redox, and nothing else
        in this subset touches oxidative or peroxide -- so only the ESR is recovered, even
        though oxidative is the arm H2O2 drives hardest."""
        from ystwin.analysis.stress_model import module_transfer

        scored = module_transfer(dataset, held_out="H2O2", n_states=3)
        peers = {m: [s for s in SUBSET if s != "H2O2" and m in STRESSORS[s].targets]
                 for m in STRESSORS["H2O2"].targets}

        assert scored["ESR"] > 0.0
        assert len(peers["ESR"]) == 6
        for module in ("oxidative", "peroxide"):
            assert peers[module] == []
            assert scored[module] == 0.0

    def test_a_module_only_the_held_out_stressor_drives_does_not(self, dataset):
        """Nothing but BPS excites iron here, so withholding BPS withholds the axis."""
        from ystwin.analysis.stress_model import module_transfer

        scored = module_transfer(dataset, held_out="BPS", n_states=3)

        assert scored["iron"] < scored["ESR"]

    def test_it_refuses_a_stressor_the_dataset_does_not_contain(self, dataset):
        from ystwin.analysis.stress_model import module_transfer

        with pytest.raises(KeyError, match="unobtainium"):
            module_transfer(dataset, held_out="unobtainium", n_states=3)

    def test_scores_are_bounded_fractions(self, dataset):
        from ystwin.analysis.stress_model import module_transfer

        scored = module_transfer(dataset, held_out="heat", n_states=3)

        assert all(0.0 <= v <= 1.0 for v in scored.values())


class TestChoosingTheLatentWidthForTheJobItWillDo:
    """Latent width has to be selected on the task the model is deployed for.

    Selecting on channel prediction picked one state for a four-channel build, which then
    reported almost nothing for the glutathione pool even though the build carries a
    glutathione sensor: one state cannot hold the general response and a redox pool at
    once. Selecting on module recovery picks a wider state and roughly doubles what the
    model reports.

    On a narrow panel the answer is often the full width, and that is not a failure. There
    is nothing to compress out of four channels; what the latent layer buys there is the
    readout that names modules, not the reduction. Compression earns its place only when
    the channels outnumber the structure behind them.
    """

    def test_it_never_exceeds_the_channels_available(self, dataset):
        from ystwin.analysis.stress_model import select_width_for_modules

        chosen = select_width_for_modules(dataset, max_states=99)

        assert chosen <= len(dataset.reporters)

    def test_it_beats_the_channel_prediction_choice_at_reporting_modules(self, dataset):
        from ystwin.analysis.stress_model import module_transfer, select_width_for_modules

        chosen = select_width_for_modules(dataset, max_states=4)

        def driven_mean(k):
            scores = [module_transfer(dataset, s, k) for s in ("DTT", "H2O2", "heat")]
            return np.mean([v for s in scores for v in s.values()])

        assert driven_mean(chosen) >= driven_mean(1)

    def test_it_returns_at_least_one_state(self, dataset):
        from ystwin.analysis.stress_model import select_width_for_modules

        assert select_width_for_modules(dataset, max_states=3) >= 1

    def test_a_narrow_panel_wants_most_of_its_width(self, dataset):
        """Nothing to compress out of a handful of channels."""
        from ystwin.analysis.stress_model import select_width_for_modules

        assert select_width_for_modules(dataset, max_states=4) >= 2


class TestMissingReadingsSurviveCentring:
    """`fit_latent` masks missing readings and says so. The centring above it must not
    throw that away: a plain column mean turns one NaN into a NaN centre, which makes the
    whole channel useless. On the Gasch arrays that cost 15 arrays and 2 channels."""

    @staticmethod
    def _with_a_hole(dataset):
        import numpy as np

        readings = dataset.readings.copy()
        readings[0, 0] = np.nan
        from ystwin.generator.panel_experiment import PanelDataset

        return PanelDataset(readings=readings, labels=dataset.labels, doses=dataset.doses,
                            modules=dataset.modules, reporters=list(dataset.reporters))

    def test_one_missing_reading_does_not_void_its_channel(self, dataset):
        holed = self._with_a_hole(dataset)

        model = train_stress_model(holed, n_states=2)

        assert np.isfinite(model.centre).all()

    def test_and_the_readout_still_predicts(self, dataset):
        holed = self._with_a_hole(dataset)
        model = train_stress_model(holed, n_states=2)

        predicted = model.predict_modules(dataset.readings).to_numpy()

        assert np.isfinite(predicted).all()
