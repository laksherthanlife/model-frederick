"""Choosing which stressor pairs to run, scored on the transfer they actually deliver.

Three analytical proxies for "does this design cover the space" failed in the same way:
eight single-agent rays already span all seven modules, so any full-span measure ranks
every design equally, and restricting to the fitted subspace is not monotonic in adding
treatments. The proxy was never the objective.

Transfer is simulable, so it is scored directly -- run the design, hold each stressor out,
take the median R2. Slower per candidate and correct by construction, which is the trade
worth making when the shortcut has been wrong three times.
"""

import pytest

from ystwin.analysis.experiment_design import (
    RECOMMENDED_DESIGN,
    design_transfer,
    select_combinations,
)
from ystwin.generator.stress_panel import STRESSORS, transcriptional_reporters

PANEL_SUBSET = ["DTT", "tunicamycin", "H2O2", "heat", "NaCl", "glucose_starvation",
                "BPS", "MG132"]
SMALL = transcriptional_reporters()[:7]
FAST = dict(doses=(0.5, 1.0, 2.0), replicates=2, n_states=3, stressors=PANEL_SUBSET,
            reporters=SMALL)


class TestScoringADesign:
    def test_a_blocked_panel_scores_worse_than_a_fully_crossed_one(self):
        from itertools import combinations

        wide = dict(FAST, stressors=["DTT", "H2O2", "heat", "NaCl", "glucose_starvation",
                                     "BPS", "MG132", "diamide", "rapamycin", "caffeine",
                                     "cobalt_chloride", "MMS"])
        blocked = design_transfer([], **wide)
        crossed = design_transfer(list(combinations(wide["stressors"], 2)), **wide)

        assert crossed > blocked

    def test_it_is_reproducible(self):
        first = design_transfer([("BPS", "heat")], **FAST)
        second = design_transfer([("BPS", "heat")], **FAST)

        assert first == pytest.approx(second)

    def test_it_refuses_an_unknown_stressor(self):
        with pytest.raises(KeyError, match="vibes"):
            design_transfer([("vibes", "heat")], **FAST)


class TestSelection:
    def test_it_returns_the_requested_number_of_pairs(self):
        assert len(select_combinations(2, library=PANEL_SUBSET, **FAST)) == 2

    def test_every_pair_is_a_real_distinct_stressor_pair(self):
        for pair in select_combinations(2, library=PANEL_SUBSET, **FAST):
            assert len(set(pair)) == 2
            assert all(s in STRESSORS for s in pair)

    def test_it_never_repeats_a_pair(self):
        chosen = select_combinations(3, library=PANEL_SUBSET, **FAST)

        assert len({frozenset(p) for p in chosen}) == len(chosen)

    def test_a_selected_design_beats_the_blocked_panel_it_started_from(self):
        chosen = select_combinations(3, library=PANEL_SUBSET, **FAST)

        assert design_transfer(chosen, **FAST) > design_transfer([], **FAST)

    def test_it_refuses_a_budget_larger_than_the_pairs_available(self):
        with pytest.raises(ValueError, match="only"):
            select_combinations(999, library=PANEL_SUBSET, **FAST)


class TestTheRecommendedDesign:
    """Six chosen pairs beat all twenty-eight, which is worth pinning down.

    Adding every pair does not add information indefinitely: the fitted subspace is
    dominated by whatever the design samples most, so a plate full of redundant mixtures
    pulls those directions toward an average stress and away from the ones that separate
    stressors. The recommendation is therefore a specific six, not "as many as fit".
    """

    def test_it_is_a_handful_of_pairs_that_fit_beside_the_singles(self):
        assert 4 <= len(RECOMMENDED_DESIGN) <= 8

    def test_every_pair_is_two_distinct_real_stressors(self):
        for pair in RECOMMENDED_DESIGN:
            assert len(set(pair)) == 2
            assert all(s in STRESSORS for s in pair)

    def test_it_beats_the_blocked_panel_it_replaces(self):
        assert design_transfer(RECOMMENDED_DESIGN, **FAST) > design_transfer([], **FAST)

    def test_it_beats_six_arbitrary_pairs(self):
        arbitrary = [("DTT", "heat"), ("H2O2", "NaCl"), ("BPS", "MG132"),
                     ("heat", "NaCl"), ("DTT", "BPS"), ("H2O2", "MG132")]

        assert design_transfer(RECOMMENDED_DESIGN, **FAST) > design_transfer(arbitrary, **FAST)

    def test_it_beats_dosing_one_agent_at_a_time_wherever_it_is_run(self):
        """The robust claim. Beating every pair holds where the design was derived -- 0.295
        against 0.121 over 25 stressors on 12 channels -- but not on an arbitrary narrower
        panel, because a design is only optimal for the readout it was optimised against."""
        wide = dict(FAST, stressors=["DTT", "H2O2", "heat", "NaCl", "glucose_starvation",
                                     "BPS", "MG132", "diamide", "rapamycin", "caffeine",
                                     "cobalt_chloride", "MMS"])

        assert design_transfer(RECOMMENDED_DESIGN, **wide) > design_transfer([], **wide)


    def test_the_advantage_is_not_an_artifact_of_one_noise_draw(self):
        scores = [design_transfer(RECOMMENDED_DESIGN, seed=s, **FAST) for s in (1, 2, 3)]
        blocked = [design_transfer([], seed=s, **FAST) for s in (1, 2, 3)]

        assert min(scores) > max(blocked)
