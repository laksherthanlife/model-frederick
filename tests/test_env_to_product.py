"""Retained PHB retrospective scores and all six fixed-gene beta-carotene requests.

`scripts/env_to_product.py` scores PHB within one source cohort, not on independent biology.
The median and the bad tail both belong here: a four-parameter fit on eleven condition
summaries is not a validated environment law. Beta-carotene has four answers and two
setpoint refusals; missing answers cannot establish carbon-source invariance. These tests
read the retained artifact, so stale outputs must fail until governed regeneration rather
than silently borrowing a fresh candidate or dropping its refusals.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_OUT = _REPO / "outputs" / "env_to_product.csv"
_STATES = _REPO / "data" / "phb" / "kocharin2013_chemostat_states.tsv"
_BETA_REQUESTS = [(carbon, rate) for carbon in ("glucose", "ethanol")
                  for rate in (0.101, 0.15, 0.2543)]
_BETA_FIELDS = {
    "refused", "refusal_reason", "growth_rate_per_h", "content_mg_per_gdw",
    "environment_layer_state", "environment_sets_the_number", "prediction_mode",
    "prediction_scope",
}


@pytest.fixture(scope="module")
def comparison():
    assert _OUT.is_file(), "retained env_to_product.csv is missing; governed regeneration required"
    return pd.read_csv(_OUT)


@pytest.fixture(scope="module")
def scored(comparison):
    return comparison[comparison["product"] == "phb"]


@pytest.fixture(scope="module")
def carotene(comparison):
    return comparison[comparison["product"] == "beta_carotene"]


def _fold_errors(frame, column):
    values = frame[[column, "measured_mg_per_gdw"]].to_numpy(dtype=float)
    assert len(values) > 0
    assert np.isfinite(values).all(), "every PHB measurement and prediction must be present"
    assert (values > 0).all(), "log-fold scoring requires positive values"
    return np.exp(np.abs(np.log(frame[column] / frame.measured_mg_per_gdw)))


def test_the_complete_product_population_is_retained(comparison):
    assert comparison["product"].value_counts(dropna=False).to_dict() == {
        "phb": 11, "beta_carotene": 6,
    }


class TestThePHBEnvironmentFitImprovesTheMedianWithinItsSourceCohort:
    def test_it_beats_the_growth_rate_alone(self, scored):
        """1.320x against 1.505x, retrospective leave-one-state-out on eleven summaries."""
        with_env = np.median(_fold_errors(scored, "predicted_with_environment"))
        mu_only = np.median(_fold_errors(scored, "predicted_mu_only"))
        assert with_env == pytest.approx(1.3201, abs=5e-5)
        assert mu_only == pytest.approx(1.5052, abs=5e-5)
        assert with_env < mu_only

    def test_it_removes_roughly_a_third_of_the_error(self, scored):
        mu_only = np.median(_fold_errors(scored, "predicted_mu_only"))
        with_env = np.median(_fold_errors(scored, "predicted_with_environment"))

        assert 0.25 < (mu_only - with_env) / (mu_only - 1.0) < 0.55

    def test_all_eleven_measured_states_are_scored(self, scored):
        states = pd.read_csv(_STATES, sep="\t").rename(
            columns={"phb_mg_per_gdw": "measured_mg_per_gdw"})
        columns = ["carbon_source", "mu_per_h", "measured_mg_per_gdw"]
        assert len(scored) == len(states) == 11
        assert scored.state_id.is_unique
        pd.testing.assert_frame_equal(
            scored.set_index("state_id")[columns].sort_index(),
            states.set_index("state_id")[columns].sort_index(), check_exact=True)
        for column in ("predicted_mu_only", "predicted_with_environment"):
            assert len(_fold_errors(scored, column)) == 11

    def test_the_scope_does_not_promote_state_holdouts_to_independent_validation(self, scored):
        assert "prediction_scope" in scored, "stale artifact lacks the retrospective scope"
        assert scored.prediction_scope.eq(
            "retrospective leave-one-state-out within one Kocharin 2013 source cohort; "
            "not independent external validation").all()


class TestAndTheTailIsStillBad:
    """The half that must not be dropped when the median is quoted. A 4-parameter law on
    11 states predicts held-out states poorly even where it improves the middle."""

    def test_the_worst_held_out_state_is_more_than_double(self, scored):
        assert _fold_errors(scored, "predicted_with_environment").max() > 2.0

    def test_more_states_are_badly_wrong_than_under_mu_alone(self, scored):
        """Adding the environment terms tightens the median and WIDENS the count of states
        off by more than 2x, from one to two. Both are true."""
        bad_env = _fold_errors(scored, "predicted_with_environment") > 2.0
        bad_mu = _fold_errors(scored, "predicted_mu_only") > 2.0

        assert set(scored.loc[bad_env, "state_id"]) == {"glc_D005", "etoh_D005"}
        assert set(scored.loc[bad_mu, "state_id"]) == {"glc_D005"}

    def test_the_separate_entry_flux_error_is_only_a_scale_comparison(self, scored):
        """22.2% is the entry-flux law's own leave-one-strain-out error, on a different
        product. PHB's median excess fold error is 32%, worse but the same order. The other
        product's error is not a PHB noise floor or independent support for this fit, and
        neither summary removes the two PHB states wrong by more than 2x."""
        from ystwin.pathway import calibrations

        median = np.median(_fold_errors(scored, "predicted_with_environment"))
        reference = calibrations.BETA_CAROTENE_FLUX.typical_fold_error

        assert reference < median < 2 * reference


class TestBetaCaroteneRetainsReachabilityRefusalsRatherThanInferringInvariance:
    def test_all_six_requested_carbon_rate_pairs_remain(self, carotene):
        assert list(zip(carotene.carbon_source, carotene.mu_per_h)) == _BETA_REQUESTS

    def test_four_answers_and_two_typed_refusals_have_different_layer_states(self, carotene):
        assert _BETA_FIELDS <= set(carotene), "stale artifact lacks refusal/scope fields"
        answered = carotene[carotene.refused.isna()]
        refused = carotene[carotene.refused.notna()]

        assert len(answered) == 4
        assert list(zip(refused.carbon_source, refused.mu_per_h)) == [
            ("ethanol", 0.15), ("ethanol", 0.2543),
        ]
        assert refused.refused.eq("SetpointUnreachable").all()
        assert refused.refusal_reason.notna().all()
        assert refused.refusal_reason.str.len().gt(0).all()
        assert refused[["growth_rate_per_h", "content_mg_per_gdw"]].isna().all().all()
        assert refused.environment_layer_state.eq("not-run").all()
        assert refused.environment_sets_the_number.isna().all()

        assert answered.refusal_reason.isna().all()
        assert np.isfinite(answered[["growth_rate_per_h", "content_mg_per_gdw"]]).all().all()
        assert answered.content_mg_per_gdw.gt(0).all()
        np.testing.assert_array_equal(answered.growth_rate_per_h, answered.mu_per_h)
        assert answered.environment_layer_state.eq("reported").all()
        assert answered.environment_sets_the_number.eq(False).all()

    def test_only_the_complete_low_rate_pair_has_equal_mu_only_answers(self, carotene):
        """Equality is a property of two returned model answers, not biological invariance."""
        assert _BETA_FIELDS <= set(carotene), "stale artifact lacks refusal/scope fields"
        pair = carotene[carotene.mu_per_h == 0.101].set_index("carbon_source")
        assert len(pair) == 2
        assert set(pair.index) == {"glucose", "ethanol"}
        assert pair.refused.isna().all()
        assert np.isfinite(pair.content_mg_per_gdw).all()
        assert pair.content_mg_per_gdw.gt(0).all()
        assert pair.loc["glucose", "content_mg_per_gdw"] == pair.loc[
            "ethanol", "content_mg_per_gdw"]

    @pytest.mark.parametrize("rate", [0.15, 0.2543])
    def test_higher_rate_carbon_comparisons_are_unavailable(self, carotene, rate):
        """A single answer must never pass a two-carbon invariance check via nunique."""
        assert _BETA_FIELDS <= set(carotene), "stale artifact lacks refusal/scope fields"
        pair = carotene[carotene.mu_per_h == rate].set_index("carbon_source")
        assert len(pair) == 2
        assert set(pair.index) == {"glucose", "ethanol"}
        assert pd.isna(pair.loc["glucose", "refused"])
        assert np.isfinite(pair.loc["glucose", "content_mg_per_gdw"])
        assert pair.loc["ethanol", "refused"] == "SetpointUnreachable"
        assert pd.isna(pair.loc["ethanol", "content_mg_per_gdw"])
        assert pair.loc["ethanol", "environment_layer_state"] == "not-run"

    def test_fixed_gene_scope_does_not_claim_environment_validation(self, carotene):
        assert _BETA_FIELDS <= set(carotene), "stale artifact lacks refusal/scope fields"
        assert carotene.prediction_mode.eq("empirical").all()
        assert carotene.prediction_scope.eq(
            "fixed-gene comparison, not environment validation").all()
