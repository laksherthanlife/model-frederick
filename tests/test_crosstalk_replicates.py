"""The crosstalk tables hold three biological replicates, and used to hold one.

`data/crosstalk/erox_2026-08_endpoint.tsv` was the workbook's first `Summary` block alone,
while `data/crosstalk/SOURCE.md` described it as "the workbook's own biological-replicate
mean across the three plates". Nothing checked, because nothing could: a derived table with
a wrong provenance note reads exactly like a right one.

These tests read the committed tables, so they run on any clone and never touch the
workbook. What they cannot check is the extraction itself -- that is
`scripts/extract_crosstalk_workbook.py`'s own gate, which recomputes every replicate-1
per-cell value from the blanked plate in the instrument export and refuses to write on a
mismatch.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ystwin import paths

CONSTRUCTS = ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1")
#: The stressor each sensor was given in this experiment: the wrong one, by design.
OFF_TARGET = {"UPRE1": "H2O2", "UPRE2": "H2O2",
              "NativeYap1": "DTT", "AlteredYap1": "DTT"}


@pytest.fixture(scope="module")
def endpoint():
    return pd.read_csv(paths.data_dir() / "crosstalk" / "erox_2026-08_endpoint.tsv",
                       sep="\t")


@pytest.fixture(scope="module")
def on_vs_off():
    return pd.read_csv(
        paths.data_dir() / "crosstalk" / "erox_2026-08_on_vs_off_target.tsv", sep="\t")


@pytest.fixture(scope="module")
def disagreement():
    return pd.read_csv(
        paths.data_dir() / "crosstalk" / "erox_2026-08_workbook_fold_disagreement.tsv",
        sep="\t")


class TestItIsThreeReplicates:
    def test_the_table_carries_three(self, endpoint):
        assert sorted(endpoint.replicate.unique()) == [1, 2, 3]

    def test_every_construct_and_dose_has_all_three(self, endpoint):
        counts = endpoint.groupby(["construct", "dose_mM"]).replicate.nunique()

        assert counts.eq(3).all()

    def test_the_controls_differ_between_plates(self, endpoint):
        """The check that would have caught the original error. If the table were one
        replicate copied three times, or one replicate alone, these would not differ --
        and they differ by more than 60%."""
        controls = endpoint[(endpoint.construct == "UPRE1") & (endpoint.dose_mM == 0.0)]

        assert controls.od600_mean.max() / controls.od600_mean.min() > 1.5

    def test_replicate_one_is_still_the_values_that_were_there_before(self, endpoint):
        """Nothing was silently re-derived. The old table's numbers survive as replicate 1,
        so the change is an addition rather than a substitution."""
        row = endpoint[(endpoint.construct == "UPRE1") & (endpoint.dose_mM == 0.0)
                       & (endpoint.replicate == 1)].iloc[0]

        assert row.od600_mean == pytest.approx(0.320333, abs=1e-6)
        assert row.rfu_mean == pytest.approx(1774.0, abs=1e-6)


class TestTheViabilityGateSeesAllThreePlates:
    """Why the sample size mattered here beyond the sample size.

    On replicate 1 UPRE1 keeps 90% of its control OD at 1 mM H2O2. On the other two it
    keeps 40% and 31%. The scorer passed that dose because the only plate it could see was
    the one where the cells were alive.
    """

    def test_upre1_viability_at_one_millimolar_differs_threefold_between_plates(
            self, endpoint):
        subset = endpoint[endpoint.construct == "UPRE1"]
        control = subset[subset.dose_mM == 0.0].set_index("replicate").od600_mean
        dosed = subset[subset.dose_mM == 1.0].set_index("replicate").od600_mean
        viability = (dosed / control).sort_index()

        assert viability.loc[1] > 0.85
        assert viability.loc[2] < 0.45
        assert viability.loc[3] < 0.35

    def test_the_scored_table_takes_the_worst_plate_not_the_mean(self):
        scored = pd.read_csv(paths.outputs_dir() / "sensor_crosstalk.csv")
        row = scored[(scored.construct == "UPRE1") & np.isclose(scored.dose_mM, 1.0)].iloc[0]

        assert row.viability_min < row.viability
        assert not row.culture_healthy

    def test_every_scored_row_carries_its_replicate_count(self):
        scored = pd.read_csv(paths.outputs_dir() / "sensor_crosstalk.csv")

        assert scored.n_replicates.eq(3).all()


class TestTheWorkbooksOwnFoldRowsDisagreeWithItsOwnData:
    """Found by checking the extraction against the instrument export rather than the sheet.

    The workbook's fold rows are typed constants. For replicate 1 the two oxidative sensors
    disagree with the per-cell values directly above them, and the per-cell values are the
    ones that recompute exactly from the raw blanked plate.
    """

    def test_the_disagreement_is_recorded_rather_than_silently_resolved(self, disagreement):
        assert {"derived_fold", "workbook_fold", "difference"} <= set(disagreement.columns)

    def test_it_is_confined_to_the_two_oxidative_sensors(self, disagreement):
        moved = disagreement[disagreement.difference.abs() > 1e-6]

        assert len(moved) > 0
        assert set(moved.construct) == {"NativeYap1", "AlteredYap1"}

    def test_the_er_sensors_agree_exactly(self, disagreement):
        er = disagreement[disagreement.construct.isin(["UPRE1", "UPRE2"])]

        assert er.difference.abs().max() < 1e-9

    def test_the_worst_difference_is_large_enough_to_matter(self, disagreement):
        assert disagreement.difference.abs().max() > 0.1


class TestTheComparisonThatMakesTheNegativeMeanSomething:
    """A sensor staying quiet for the wrong stressor says nothing until you know it shouts
    for the right one. A dead sensor also stays quiet."""

    def test_both_stressors_are_present_for_every_sensor(self, on_vs_off):
        assert set(on_vs_off.construct) == set(CONSTRUCTS)
        for construct in CONSTRUCTS:
            row = on_vs_off[on_vs_off.construct == construct].iloc[0]
            assert row.off_target_stressor == OFF_TARGET[construct]
            assert row.on_target_stressor != row.off_target_stressor

    @pytest.mark.parametrize("construct", CONSTRUCTS)
    def test_every_sensor_responds_more_to_its_own_stressor(self, on_vs_off, construct):
        """The result, in one assertion per sensor. Compared at the highest dose where the
        on-target culture is still alive -- 1.0 mM, since above it the oxidative sensors'
        cultures collapse and their folds go negative, which is death rather than
        repression."""
        row = on_vs_off[(on_vs_off.construct == construct)
                        & np.isclose(on_vs_off.on_target_dose_mM, 1.0)].iloc[0]

        assert row.on_target_fold > row.off_target_fold_derived

    def test_the_er_sensors_roughly_double_for_their_own_stressor(self, on_vs_off):
        peak = on_vs_off[on_vs_off.construct.isin(["UPRE1", "UPRE2"])
                         & np.isclose(on_vs_off.on_target_dose_mM, 2.0)]

        assert (peak.on_target_fold > 1.9).all()

    def test_and_fall_below_one_for_the_wrong_one(self, on_vs_off):
        same = on_vs_off[on_vs_off.construct.isin(["UPRE1", "UPRE2"])
                         & np.isclose(on_vs_off.off_target_dose_mM, 2.0)]

        assert (same.off_target_fold_derived < 0.6).all()

    def test_the_derived_off_target_column_is_carried_beside_the_workbooks(self, on_vs_off):
        """Both are in the committed table, so preferring the derived one is auditable
        rather than asserted."""
        assert "off_target_fold" in on_vs_off.columns
        assert "off_target_fold_derived" in on_vs_off.columns

    def test_the_negative_folds_are_only_where_the_culture_is_dead(self, on_vs_off):
        negative = on_vs_off[on_vs_off.on_target_fold < 0]

        assert set(negative.construct) <= {"NativeYap1", "AlteredYap1"}
        assert (negative.on_target_dose_mM >= 2.0).all()
