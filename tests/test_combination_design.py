"""Dosing stressors together, and whether it is what makes a stress state transferable.

One stressor at a time produces one ray per stressor, so a held-out stressor's direction
is reachable only if some training stressor happens to point the same way. That is a
property of the experiment, not of the biology, and it is the reason leave-one-out
transfer is weak on a blocked panel.

Combinations span the directions between stressors, so the training set covers ground no
single-agent panel reaches. If that is the real constraint then adding combinations must
lift transfer for stressors whose modules are excited by nothing else on its own -- and
must still fail for a module no combination touches at all, which is the control that
keeps the claim honest.
"""

import numpy as np
import pytest

from ystwin.analysis.transfer import leave_one_stressor_out, subspace_alignment, transfer_test
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.stress_panel import (
    combination_response,
    transcriptional_reporters,
)

ALL_REPORTERS = transcriptional_reporters()[:7]
PANEL_SUBSET = ["DTT", "tunicamycin", "H2O2", "heat", "NaCl", "glucose_starvation",
                "BPS", "MG132"]



class TestCombinationResponse:
    def test_two_stressors_raise_a_shared_module_above_either_alone(self):
        alone = combination_response({"DTT": 1.0})
        together = combination_response({"DTT": 1.0, "heat": 6.0})

        assert together["ESR"] > alone["ESR"]

    def test_it_reduces_to_a_single_stressor(self):
        from ystwin.generator.stress_panel import module_response

        assert combination_response({"H2O2": 0.5}) == pytest.approx(module_response("H2O2", 0.5))

    def test_a_module_no_applied_stressor_touches_stays_at_zero(self):
        response = combination_response({"DTT": 1.0, "heat": 6.0})

        assert response["iron"] == 0.0

    def test_it_refuses_an_unknown_stressor(self):
        with pytest.raises(KeyError, match="vibes"):
            combination_response({"vibes": 1.0})

    def test_an_empty_combination_is_unstressed(self):
        assert all(v == 0.0 for v in combination_response({}).values())


class TestCombinationsInADataset:
    def test_a_combination_is_labelled_by_its_members(self):
        data = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, combinations=[("DTT", "heat")],
                             noise_cv=0.0, seed=0)

        assert "DTT+heat" in set(data.labels)

    def test_combinations_add_samples_without_removing_singles(self):
        singles = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.0, seed=0)
        mixed = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, combinations=[("DTT", "heat")],
                              noise_cv=0.0, seed=0)

        assert set(PANEL_SUBSET) <= set(mixed.labels)
        assert len(mixed.labels) > len(singles.labels)

    def test_a_combination_lies_off_both_of_its_single_stressor_rays(self):
        """Measured on the induced component. Every reading sits on a constitutive floor
        that is larger than what induction adds to it, so raw readings all point nearly the
        same way whatever the treatment -- the geometry that matters is what the dose moved,
        not where it started."""
        from ystwin.generator.stress_panel import REPORTERS

        data = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                             combinations=[("NaCl", "BPS")], noise_cv=0.0, seed=0)
        floor = np.array([REPORTERS[r].basal for r in ALL_REPORTERS])
        induced = data.readings - floor
        mixed = induced[data.labels == "NaCl+BPS"]
        salt = induced[data.labels == "NaCl"]

        cosines = (mixed @ salt.T) / np.outer(
            np.linalg.norm(mixed, axis=1), np.linalg.norm(salt, axis=1))

        assert np.nanmax(cosines) < 0.999


class TestCombinationsAreWhatMakeTransferWork:
    def test_they_cannot_rescue_a_module_nothing_else_touches(self):
        """The boundary of the method, and it is worth pinning rather than hoping otherwise.

        Combinations work by visiting the directions between agents, so they need at least
        two agents that reach a module. Nothing but BPS excites the iron regulon, and every
        treatment containing BPS is withdrawn along with BPS when it is held out -- so its
        pairs leave with it and the iron axis is as absent as it ever was. Pairing it more
        aggressively makes things worse, not better, by spending wells on directions that
        vanish at test time.
        """
        blocked = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                                noise_cv=0.02, seed=5)
        paired = panel_dataset(
            reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.02, seed=5,
            combinations=[("BPS", "MG132"), ("BPS", "DTT"), ("BPS", "heat")])

        before = transfer_test(blocked, held_out="BPS", n_states=3, observed=(0, 1, 2, 3))
        after = transfer_test(paired, held_out="BPS", n_states=3, observed=(0, 1, 2, 3))

        assert after.r2 <= before.r2

    def test_the_fix_for_a_private_module_is_a_second_agent_that_reaches_it(self):
        """Not a combination -- another way in. Cobalt perturbs iron handling too, so with
        it in the panel the axis survives BPS being withdrawn."""
        without = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                                noise_cv=0.02, seed=5)
        with_second = panel_dataset(reporters=ALL_REPORTERS,
                                    stressors=PANEL_SUBSET + ["cobalt_chloride"],
                                    noise_cv=0.02, seed=5)

        before = transfer_test(without, held_out="BPS", n_states=3, observed=(0, 1, 2, 3))
        after = transfer_test(with_second, held_out="BPS", n_states=3, observed=(0, 1, 2, 3))

        assert after.r2 > before.r2

    def test_they_raise_the_subspace_alignment_that_explains_it(self):
        blocked = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.02, seed=5)
        spanned = panel_dataset(
            reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.02, seed=5,
            combinations=[("BPS", "heat"), ("BPS", "NaCl"), ("DTT", "NaCl"), ("H2O2", "heat")])

        assert subspace_alignment(spanned, "BPS", 3) > subspace_alignment(blocked, "BPS", 3)

    def test_a_chosen_design_lifts_median_transfer_across_the_whole_panel(self):
        from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN

        blocked = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                                noise_cv=0.02, seed=5)
        spanned = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                                noise_cv=0.02, seed=5, combinations=RECOMMENDED_DESIGN)

        before = leave_one_stressor_out(blocked, n_states=3, observed=(0, 1, 2, 3))
        after = leave_one_stressor_out(spanned, n_states=3, observed=(0, 1, 2, 3))

        assert after["r2"].median() > before["r2"].median()

    def test_combinations_help_generally_once_the_ladder_stays_healthy(self):
        """This one reversed when the dose ladder was corrected, which is worth recording.

        On a ladder running past the lethal dose, arbitrary pairs scored below no pairs at
        all -- the extra wells were spent in the range where nothing reports. Dosed inside
        the healthy range they help, and the chosen design still helps more. The lesson is
        not that combinations are risky; it is that a design cannot be judged on a ladder
        that wastes half its wells on dead cultures.
        """
        blocked = panel_dataset(reporters=ALL_REPORTERS, stressors=PANEL_SUBSET,
                                noise_cv=0.02, seed=5)
        arbitrary = panel_dataset(
            reporters=ALL_REPORTERS, stressors=PANEL_SUBSET, noise_cv=0.02, seed=5,
            combinations=[("DTT", "heat"), ("H2O2", "NaCl"), ("BPS", "MG132"),
                          ("heat", "NaCl"), ("DTT", "BPS"), ("H2O2", "MG132")])

        before = leave_one_stressor_out(blocked, n_states=3, observed=(0, 1, 2, 3))
        after = leave_one_stressor_out(arbitrary, n_states=3, observed=(0, 1, 2, 3))

        assert after["r2"].median() > before["r2"].median()
