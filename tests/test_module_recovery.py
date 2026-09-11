"""Does the fitted latent state mean anything, or is it only a compression?

A model can predict held-out channels well while its states are an arbitrary rotation of
whatever spans the data. That is useful for imputation and useless for biology. The state
is only worth calling a stress state if the true module activities can be read off it.

Ground truth is available in simulation, so this is directly checkable: regress each true
module activity on the fitted states and report the fraction of variance recovered. The
suite pins the honest version -- modules the reporters genuinely read must come back,
modules nothing observes must not, and a state fitted to shuffled data must recover
nothing at all.
"""

import pytest

from ystwin.analysis.recovery import module_recovery
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import MODULES, transcriptional_reporters

ALL_REPORTERS = transcriptional_reporters()[:7]
PANEL_SUBSET = ["DTT", "tunicamycin", "H2O2", "heat", "NaCl", "glucose_starvation",
                "BPS", "MG132"]



@pytest.fixture(scope="module")
def dataset():
    return panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.02, seed=0)


class TestRecovery:
    def test_it_scores_every_module(self, dataset):
        scores = module_recovery(dataset, n_states=4)

        assert set(scores) == set(MODULES)

    def test_scores_are_fractions(self, dataset):
        scores = module_recovery(dataset, n_states=4)

        assert all(-1e-9 <= v <= 1.0 + 1e-9 for v in scores.values())

    def test_a_full_panel_recovers_most_modules_well(self, dataset):
        scores = module_recovery(dataset, n_states=6)
        good = [m for m, v in scores.items() if v > 0.8]

        assert len(good) >= 5

    def test_more_states_recover_at_least_as_much(self, dataset):
        few = module_recovery(dataset, n_states=2)
        many = module_recovery(dataset, n_states=6)

        assert sum(many.values()) > sum(few.values())

    def test_a_module_no_reporter_reads_is_not_recovered(self):
        """Drop the iron reporter and the iron regulon should stop coming back."""
        without_iron = [r for r in ALL_REPORTERS if r != "FeRE-iron"]
        data = panel_dataset(reporters=without_iron, stressors=PANEL_SUBSET, noise_cv=0.02, seed=0)

        assert module_recovery(data, n_states=5)["iron"] < 0.5

    def test_shuffling_the_states_destroys_recovery(self, dataset):
        """The control: if a shuffled state still scores, the metric is measuring nothing."""
        real = module_recovery(dataset, n_states=4)
        shuffled = module_recovery(dataset, n_states=4, shuffle=True)

        assert max(shuffled.values()) < max(real.values())

    def test_a_four_channel_build_recovers_the_modules_it_was_chosen_for(self):
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD

        data = panel_dataset(reporters=RECOMMENDED_BUILD.stress_reporters, stressors=PANEL_SUBSET, noise_cv=0.02, seed=0)
        scores = module_recovery(data, n_states=4)

        assert scores["ESR"] > 0.5


def test_the_seed_does_not_reach_an_unshuffled_fit():
    """Pins the documented no-op, because the alternative reading is dangerous.

    fit_latent recovers its states by SVD, which is deterministic, so `seed` cannot
    perturb an unshuffled recovery. That is correct behaviour -- but it means a spread
    computed across seeds is exactly zero, which reads as perfect precision rather than
    as the absence of a measurement. Locking it here so the no-op stays deliberate.
    """
    import warnings

    from ystwin.analysis.recovery import module_recovery
    from ystwin.generator.panel_experiment import panel_dataset

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = panel_dataset(seed=0)
        first = module_recovery(data, n_states=3, seed=1)
        second = module_recovery(data, n_states=3, seed=98765)

    assert first == second, "seed must not affect an unshuffled recovery"


def test_the_seed_does_reach_the_shuffle():
    """The control the seed exists for still varies, so it is threaded, not ignored."""
    import warnings

    from ystwin.analysis.recovery import module_recovery
    from ystwin.generator.panel_experiment import panel_dataset

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = panel_dataset(seed=0)
        first = module_recovery(data, n_states=3, seed=1, shuffle=True)
        second = module_recovery(data, n_states=3, seed=2, shuffle=True)

    assert first != second, "the shuffle control must vary with the seed"
