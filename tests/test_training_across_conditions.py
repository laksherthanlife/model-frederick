"""Training across conditions, not just across stressors.

A dataset generated in one context teaches a model that every raised ESR is a stressor,
because in exponential glucose it always is. On a real plate read late it is not, and the
model has no way to tell -- it never saw the case.

Two things have to hold for the fix to be worth anything. The context must reach the
readings, so a stationary-phase well genuinely looks different from an exponential one. And
a model trained across contexts must do better on a context it was not trained in than one
trained in a single context does, which is the same held-out argument used for stressors and
the only evidence that the breadth bought anything.
"""

import numpy as np

from ystwin.analysis.stress_model import module_transfer, train_stress_model
from ystwin.generator.context import CultureContext
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import transcriptional_reporters

READERS = transcriptional_reporters()[:8]
SUB = ["DTT", "H2O2", "heat", "NaCl", "MG132"]


class TestContextReachesTheReadings:
    def test_a_dataset_records_the_context_of_each_sample(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0,
                             contexts=[CultureContext(), CultureContext(growth_phase="stationary")])

        assert len(data.contexts) == len(data.labels)
        assert len(set(data.contexts)) == 2

    def test_a_stationary_well_reads_differently_from_an_exponential_one(self):
        data = panel_dataset(reporters=["STRE-general"], stressors=["DTT"], noise_cv=0.0,
                             contexts=[CultureContext(), CultureContext(growth_phase="stationary")])
        exponential = data.readings[np.array(data.contexts) == "exponential/glucose"]
        stationary = data.readings[np.array(data.contexts) == "stationary/glucose"]

        assert stationary.mean() > exponential.mean()

    def test_the_baseline_appears_in_the_true_module_activities(self):
        """Not just in the reading: a stationary culture really has a raised ESR."""
        data = panel_dataset(stressors=["DTT"], doses=(0.0,), noise_cv=0.0, replicates=1,
                             contexts=[CultureContext(growth_phase="stationary")])
        from ystwin.generator.stress_panel import MODULES

        assert data.modules[0, list(MODULES).index("ESR")] > 0.3

    def test_a_slower_context_carries_a_larger_growth_confound(self):
        """Stationary cultures barely grow, so dividing growth out is far less certain."""
        def spread(context):
            data = panel_dataset(reporters=["STRE-general"], stressors=["DTT"], noise_cv=0.0,
                                 growth_rate_se=0.0117, replicates=12, seed=0,
                                 contexts=[context])
            return float(np.std(data.readings) / np.mean(data.readings))

        assert spread(CultureContext(growth_phase="stationary")) > spread(CultureContext())

    def test_one_context_is_still_the_default(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0)

        assert set(data.contexts) == {"exponential/glucose"}


class TestBreadthBuysGeneralisation:
    def test_a_model_trained_in_one_context_misreads_another(self):
        """The failure the breadth is for: raised ESR with no stressor behind it."""
        narrow = panel_dataset(reporters=READERS, stressors=SUB, noise_cv=0.05, seed=0,
                               contexts=[CultureContext()])
        model = train_stress_model(narrow, n_states=4)

        held_out = panel_dataset(reporters=READERS, stressors=SUB, doses=(0.0,), noise_cv=0.05,
                                 seed=1, contexts=[CultureContext(growth_phase="stationary")])
        predicted = model.predict_modules(held_out.readings)["ESR"].mean()
        truth = held_out.modules_frame()["ESR"].mean()

        assert abs(predicted - truth) > 0.1

    def test_training_across_contexts_narrows_that_error(self):
        contexts = [CultureContext(), CultureContext(growth_phase="diauxic"),
                    CultureContext(carbon_source="galactose")]
        broad = panel_dataset(reporters=READERS, stressors=SUB, noise_cv=0.05, seed=0,
                              contexts=contexts)
        narrow = panel_dataset(reporters=READERS, stressors=SUB, noise_cv=0.05, seed=0,
                               contexts=[CultureContext()])

        held_out = panel_dataset(reporters=READERS, stressors=SUB, doses=(0.0,), noise_cv=0.05,
                                 seed=1, contexts=[CultureContext(growth_phase="stationary")])
        truth = held_out.modules_frame()["ESR"].mean()

        def error(data):
            model = train_stress_model(data, n_states=4)
            return abs(model.predict_modules(held_out.readings)["ESR"].mean() - truth)

        assert error(broad) < error(narrow)

    def test_stressor_transfer_still_works_across_contexts(self):
        contexts = [CultureContext(), CultureContext(growth_phase="diauxic")]
        data = panel_dataset(reporters=READERS, stressors=SUB, noise_cv=0.05, seed=0,
                             contexts=contexts)
        scored = module_transfer(data, held_out="H2O2", n_states=4)

        assert scored["ESR"] > 0.0


class TestAWellTooSlowToCorrectIsUnusable:
    """Dividing growth out of a culture that is not growing does not give a bad number.

    Activity is recovered as k = dR/dt + mu*R, so the correction divides by mu. In a
    stationary culture mu is near zero and the estimate is not merely noisy -- a growth
    estimate that lands slightly negative sends the recovered activity below zero, which is
    not a weak promoter but an arithmetic artifact. The real plates showed exactly this at
    2 mM peroxide, where blank-corrected signal came back negative.

    An experimenter drops those wells, and so does the generator. Producing them and letting
    a model train on them would teach it that promoters run backwards in stationary phase.
    """

    def test_a_stationary_context_yields_fewer_usable_wells_than_an_exponential_one(self):
        def usable(phase):
            data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117,
                                 replicates=20, seed=0,
                                 contexts=[CultureContext(growth_phase=phase)])
            return len(data.labels)

        assert usable("stationary") < usable("exponential")

    def test_no_reading_that_survives_is_negative(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117,
                             replicates=30, seed=0,
                             contexts=[CultureContext(growth_phase="stationary")])

        assert data.readings.min() >= 0

    def test_it_records_how_many_wells_it_dropped(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117,
                             replicates=30, seed=0,
                             contexts=[CultureContext(growth_phase="stationary")])

        assert data.n_unusable > 0

    def test_an_exponential_culture_loses_almost_nothing(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0117,
                             replicates=30, seed=0, contexts=[CultureContext()])

        assert data.n_unusable == 0

    def test_a_dataset_with_no_confound_drops_nothing(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, growth_rate_se=0.0,
                             replicates=10, contexts=[CultureContext(growth_phase="stationary")])

        assert data.n_unusable == 0
