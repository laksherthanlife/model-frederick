"""Multi-channel latent stress state, with its dimension chosen by evidence.

The plan's warning is that a latent state must not be a reporter rename. The check is
mechanical: if two channels are driven by one underlying process the model should find
one factor, and if they are driven by two it should find two -- regardless of how many
sensors were plugged in.

Inputs are dilution-corrected promoter activities, not raw RFU/OD, so the growth
confound is removed before anything is called latent.
"""

import numpy as np
import pytest

from ystwin.analysis.latent import (
    fit_latent,
    select_dimension,
)


def _activities(n_conditions=40, loadings=None, n_latent=1, noise=0.05, seed=0):
    """Channels = latent states x loadings + noise."""
    rng = np.random.default_rng(seed)
    states = rng.normal(size=(n_conditions, n_latent))
    if loadings is None:
        loadings = np.ones((n_latent, 2))
    signal = states @ np.asarray(loadings, dtype=float)
    return signal + rng.normal(0, noise, signal.shape), states


class TestDimensionIsChosenByEvidence:
    def test_two_channels_driven_by_one_process_give_one_factor(self):
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.8]])

        assert select_dimension(data, max_states=3) == 1

    def test_two_channels_driven_by_two_processes_give_two(self):
        data, _ = _activities(n_latent=2, loadings=[[1.0, 0.0], [0.0, 1.0]])

        assert select_dimension(data, max_states=3) == 2

    def test_adding_a_pure_noise_channel_does_not_add_a_state(self):
        """A sensor that reports nothing must not create a latent state."""
        rng = np.random.default_rng(1)
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.9]])
        with_noise = np.column_stack([data, rng.normal(0, 0.05, len(data))])

        assert select_dimension(with_noise, max_states=3) == 1

    def test_a_duplicated_channel_does_not_add_a_state(self):
        """Two copies of one sensor are one sensor."""
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.9]])
        duplicated = np.column_stack([data, data[:, 0]])

        assert select_dimension(duplicated, max_states=3) == 1

    def test_the_dimension_never_exceeds_the_channel_count(self):
        data, _ = _activities(n_latent=1)

        assert select_dimension(data, max_states=8) <= data.shape[1]


class TestTheFit:
    def test_it_recovers_a_single_shared_state(self):
        data, truth = _activities(n_latent=1, loadings=[[1.0, 0.85]], noise=0.02)

        fit = fit_latent(data, n_states=1)
        r = abs(np.corrcoef(fit.states[:, 0], truth[:, 0])[0, 1])

        assert r > 0.95

    def test_it_reports_loadings_per_channel(self):
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.5]], noise=0.02)

        fit = fit_latent(data, n_states=1)

        assert fit.loadings.shape == (1, 2)
        assert abs(fit.loadings[0, 0]) > abs(fit.loadings[0, 1])

    def test_it_reports_how_much_each_channel_is_explained(self):
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.85]], noise=0.02)

        fit = fit_latent(data, n_states=1)

        assert all(0.0 <= v <= 1.0 for v in fit.explained_per_channel)
        assert min(fit.explained_per_channel) > 0.8

    def test_an_unexplained_channel_is_visible_as_such(self):
        rng = np.random.default_rng(2)
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.9]], noise=0.02)
        with_noise = np.column_stack([data, rng.normal(0, 1.0, len(data))])

        fit = fit_latent(with_noise, n_states=1)

        assert fit.explained_per_channel[2] < 0.3


class TestMissingChannels:
    def test_a_missing_reading_is_not_zero_filled(self):
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.85]], noise=0.02)
        gapped = data.copy()
        gapped[5:15, 1] = np.nan

        fit = fit_latent(gapped, n_states=1)

        assert np.isfinite(fit.states).all()

    def test_a_gap_does_not_shift_the_recovered_state(self):
        data, truth = _activities(n_latent=1, loadings=[[1.0, 0.85]], noise=0.02)
        gapped = data.copy()
        gapped[5:15, 1] = np.nan

        complete = fit_latent(data, n_states=1)
        partial = fit_latent(gapped, n_states=1)
        r = abs(np.corrcoef(complete.states[:, 0], partial.states[:, 0])[0, 1])

        assert r > 0.95

    def test_zero_filling_would_have_shifted_it(self):
        """Why the mask matters, demonstrated rather than asserted."""
        data, _ = _activities(n_latent=1, loadings=[[1.0, 0.85]], noise=0.02)
        gapped, zeroed = data.copy(), data.copy()
        gapped[5:15, 1] = np.nan
        zeroed[5:15, 1] = 0.0

        masked = fit_latent(gapped, n_states=1)
        naive = fit_latent(zeroed, n_states=1)
        complete = fit_latent(data, n_states=1)

        r_masked = abs(np.corrcoef(masked.states[:, 0], complete.states[:, 0])[0, 1])
        r_naive = abs(np.corrcoef(naive.states[:, 0], complete.states[:, 0])[0, 1])
        assert r_masked > r_naive

    def test_a_channel_that_is_entirely_missing_is_refused(self):
        data, _ = _activities(n_latent=1)
        data[:, 1] = np.nan

        with pytest.raises(ValueError, match="no observations"):
            fit_latent(data, n_states=1)


def test_more_states_than_channels_is_refused():
    data, _ = _activities(n_latent=1)

    with pytest.raises(ValueError, match="cannot exceed"):
        fit_latent(data, n_states=5)


def test_the_fit_reports_evidence_so_dimensions_can_be_compared():
    data, _ = _activities(n_latent=1, loadings=[[1.0, 0.85]])

    one = fit_latent(data, n_states=1)
    two = fit_latent(data, n_states=2)

    assert np.isfinite(one.heldout_error) and np.isfinite(two.heldout_error)
    assert one.heldout_error < two.heldout_error


def test_an_in_sample_criterion_would_have_picked_the_wrong_dimension():
    """Why selection is cross-validated: an exact fit wins on any in-sample score."""
    data, _ = _activities(n_latent=1, loadings=[[1.0, 0.85]])

    one = fit_latent(data, n_states=1)
    two = fit_latent(data, n_states=2)

    assert two.bic < one.bic
    assert select_dimension(data, max_states=2) == 1


class TestTwoStressorsAreNotOneLatentState:
    """The demonstration that matters for a *general* stress state.

    UPRE responds to DTT, Yap1 to peroxide. Driving all four channels off one dose axis
    would collapse them to one factor trivially -- an artefact of the simulation, not a
    property of the biology. Run on the real design, where the two stressors vary
    independently, the model must find two states unless something genuinely couples them.
    """

    def _channels(self, general_stress_weight, seed=0, n_per_axis=9):
        rng = np.random.default_rng(seed)
        dtt = rng.uniform(0, 2, n_per_axis)
        peroxide = rng.uniform(0, 2, n_per_axis)
        # Each construct answers its own stressor, plus a shared burden term.
        dtt_axis = np.concatenate([dtt, np.zeros(n_per_axis)])
        h2o2_axis = np.concatenate([np.zeros(n_per_axis), peroxide])
        general = general_stress_weight * (dtt_axis + h2o2_axis)
        columns = [
            dtt_axis * 1.0 + general, dtt_axis * 0.85 + general,
            h2o2_axis * 1.0 + general, h2o2_axis * 0.9 + general,
        ]
        data = np.column_stack(columns)
        return data + rng.normal(0, 0.03, data.shape)

    def _crossed(self, general_stress_weight, seed=0, n=24):
        """Both stressors present at independently varied levels."""
        rng = np.random.default_rng(seed)
        dtt_axis = rng.uniform(0, 2, n)
        h2o2_axis = rng.uniform(0, 2, n)
        general = general_stress_weight * (dtt_axis + h2o2_axis)
        data = np.column_stack([
            dtt_axis * 1.0 + general, dtt_axis * 0.85 + general,
            h2o2_axis * 1.0 + general, h2o2_axis * 0.9 + general,
        ])
        return data + rng.normal(0, 0.03, data.shape)

    def test_two_independent_stressors_give_two_states(self):
        data = self._channels(general_stress_weight=0.0)

        assert select_dimension(data, max_states=4) == 2

    def test_the_channels_split_along_the_stressor_they_sense(self):
        data = self._channels(general_stress_weight=0.0)

        fit = fit_latent(data, n_states=2)
        er = fit.loadings[:, :2]
        oxidative = fit.loadings[:, 2:]
        dominant_er = int(np.argmax(np.abs(er).sum(axis=1)))
        dominant_ox = int(np.argmax(np.abs(oxidative).sum(axis=1)))

        assert dominant_er != dominant_ox

    def test_a_strong_shared_burden_term_collapses_them_toward_one(self):
        """What a genuinely general stress state would look like in the data."""
        specific = self._channels(general_stress_weight=0.0)
        shared = self._channels(general_stress_weight=4.0)

        one_state_specific = fit_latent(specific, n_states=1).explained_per_channel.min()
        one_state_shared = fit_latent(shared, n_states=1).explained_per_channel.min()

        assert one_state_shared > one_state_specific

    def test_a_second_state_explains_substantially_more(self):
        data = self._channels(general_stress_weight=0.0)

        one = fit_latent(data, n_states=1).explained_per_channel.min()
        two = fit_latent(data, n_states=2).explained_per_channel.min()

        assert two > one + 0.1

    def test_applying_stressors_one_at_a_time_leaves_them_anti_correlated(self):
        """A design consequence: alone-only dosing is not the same as independent.

        With one stressor on at a time the other is off, so a single factor captures
        which block a sample came from and explains most of both channels without any
        shared biology. Crossing the two stressors removes that.
        """
        blocked = self._channels(general_stress_weight=0.0)
        crossed = self._crossed(general_stress_weight=0.0)

        one_blocked = fit_latent(blocked, n_states=1).explained_per_channel.min()
        one_crossed = fit_latent(crossed, n_states=1).explained_per_channel.min()

        assert one_blocked > one_crossed

    def test_the_crossed_design_still_finds_two_states(self):
        crossed = self._crossed(general_stress_weight=0.0)

        assert select_dimension(crossed, max_states=4) == 2

    def test_a_single_stressor_axis_cannot_reveal_this_at_all(self):
        """Why the existing dataset could never have answered the question."""
        rng = np.random.default_rng(3)
        dose = rng.uniform(0, 2, 18)
        data = np.column_stack([dose * w for w in (1.0, 0.85, 0.95, 0.9)])
        data = data + rng.normal(0, 0.03, data.shape)

        assert select_dimension(data, max_states=4) == 1
