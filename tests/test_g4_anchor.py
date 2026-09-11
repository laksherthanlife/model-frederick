"""G4: does the biosensor track an independent anchor, or does it track growth?

The test has three possible answers, and keeping them apart is the whole point:

  PASS          the anchor moved and the reporter tracked it
  REFUTED       the anchor moved and the reporter did not, or moved with growth instead
  INCONCLUSIVE  the anchor did not move detectably, so nothing can be concluded

An underpowered anchor returning INCONCLUSIVE must never be reported as PASS. The
usual failure in practice is the opposite: a flat anchor and a moving reporter get
written up as agreement because nobody checked whether the anchor had any signal.
"""

import numpy as np
import pandas as pd

from ystwin.gates.g4_anchor import anchor_agreement


def _frame(doses, anchor_by_rep, reporter, growth=None):
    """anchor_by_rep: list of per-replicate fold-change lists."""
    rows = []
    for rep, folds in enumerate(anchor_by_rep):
        for dose, fold in zip(doses, folds):
            rows.append({
                "replicate": rep,
                "dose_mM": dose,
                "anchor_fold_change": fold,
                "reporter_response": reporter[doses.index(dose)],
                "growth_rate": (growth or [0.30] * len(doses))[doses.index(dose)],
            })
    return pd.DataFrame(rows)


DOSES = [0.0, 0.2, 0.5]


def test_a_reporter_tracking_a_strong_consistent_anchor_passes():
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]],
        reporter=[1.0, 4.0, 9.2],
    )

    result = anchor_agreement(frame)

    assert result.verdict == "PASS"


def test_an_anchor_smaller_than_its_own_replicate_spread_is_inconclusive():
    """The real UPRE1 / Hac1 values: one replicate falls, the other rises 2.4 fold."""
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 1.29, 0.90], [1.0, 0.12, 2.44]],
        reporter=[1.0, 3.0, 6.0],
    )

    result = anchor_agreement(frame)

    assert result.verdict == "INCONCLUSIVE"
    assert "anchor" in result.reason


def test_an_inconclusive_anchor_is_never_reported_as_agreement():
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 1.05, 0.95], [1.0, 0.98, 1.02]],
        reporter=[1.0, 1.04, 0.96],
    )

    result = anchor_agreement(frame)

    assert result.verdict == "INCONCLUSIVE"


def test_a_reporter_that_ignores_a_responding_anchor_is_refuted():
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 5.0, 12.0], [1.0, 4.8, 11.5]],
        reporter=[1.0, 1.02, 0.99],
    )

    result = anchor_agreement(frame)

    assert result.verdict == "REFUTED"


def test_a_reporter_explained_by_growth_rate_is_refuted_even_if_it_correlates():
    """Both anchor and reporter rise with dose; only the reporter tracks growth."""
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 4.0, 9.0], [1.0, 4.2, 8.8]],
        reporter=[1.0, 4.0, 9.0],
        growth=[0.30, 0.15, 0.066],
    )

    result = anchor_agreement(frame, growth_confound_r2=0.9)

    assert result.verdict == "REFUTED"
    assert "growth" in result.reason


def test_the_effect_size_and_spread_behind_the_verdict_are_reported():
    frame = _frame(DOSES, [[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]], [1.0, 4.0, 9.2])

    result = anchor_agreement(frame)

    assert result.anchor_effect > 1.0
    assert result.anchor_spread >= 0.0
    assert np.isfinite(result.reporter_anchor_r)


def test_a_single_replicate_cannot_establish_the_anchors_reliability():
    frame = _frame(DOSES, [[1.0, 4.0, 9.0]], [1.0, 4.0, 9.0])

    result = anchor_agreement(frame)

    assert result.verdict == "INCONCLUSIVE"
    assert "replicate" in result.reason


def test_a_single_dose_cannot_establish_a_dose_response():
    frame = _frame([0.0], [[1.0], [1.0]], [1.0])

    result = anchor_agreement(frame)

    assert result.verdict == "INCONCLUSIVE"
    assert "dose" in result.reason


def test_the_minimum_detectable_effect_is_reported_for_planning():
    """Tells the team how many replicates would be needed to settle it."""
    frame = _frame(
        DOSES,
        anchor_by_rep=[[1.0, 1.29, 0.90], [1.0, 0.12, 2.44]],
        reporter=[1.0, 3.0, 6.0],
    )

    result = anchor_agreement(frame)

    assert result.replicates_needed > 2


class TestAnchorQualityGatesTheVerdict:
    """A verdict computed on an invalid anchor is not a verdict.

    If the no-RT control shows the anchor signal is largely genomic DNA, then
    'the anchor responded and the reporter did not track it' is not evidence
    against the reporter -- it is evidence against the measurement.
    """

    def test_a_failing_anchor_qc_forces_inconclusive_not_refuted(self):
        frame = _frame(
            DOSES,
            anchor_by_rep=[[1.0, 5.0, 12.0], [1.0, 4.8, 11.5]],
            reporter=[1.0, 1.02, 0.99],
        )

        without_qc = anchor_agreement(frame)
        with_qc = anchor_agreement(frame, anchor_qc_pass_rate=0.0)

        assert without_qc.verdict == "REFUTED"
        assert with_qc.verdict == "INCONCLUSIVE"
        assert "control" in with_qc.reason

    def test_a_failing_anchor_qc_also_blocks_a_pass(self):
        frame = _frame(
            DOSES,
            anchor_by_rep=[[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]],
            reporter=[1.0, 4.0, 9.2],
        )

        assert anchor_agreement(frame, anchor_qc_pass_rate=0.2).verdict == "INCONCLUSIVE"

    def test_a_passing_anchor_qc_leaves_the_verdict_alone(self):
        frame = _frame(
            DOSES,
            anchor_by_rep=[[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]],
            reporter=[1.0, 4.0, 9.2],
        )

        assert anchor_agreement(frame, anchor_qc_pass_rate=1.0).verdict == "PASS"


class TestGrowthIsNotEvidenceWhenItDidNotVary:
    def test_a_constant_growth_rate_is_reported_as_not_evaluated(self):
        """R^2 = 0.00 would read as 'checked and cleared'. It was not checked."""
        frame = _frame(DOSES, [[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]], [1.0, 4.0, 9.2])

        result = anchor_agreement(frame)

        assert np.isnan(result.growth_r2)
        assert "growth rate was not measured" in result.reason

    def test_a_varying_growth_rate_is_actually_evaluated(self):
        frame = _frame(
            DOSES, [[1.0, 4.1, 9.0], [1.0, 3.9, 9.4]], [1.0, 4.0, 9.2],
            growth=[0.30, 0.28, 0.31],
        )

        result = anchor_agreement(frame)

        assert np.isfinite(result.growth_r2)
        assert "not measured" not in result.reason
