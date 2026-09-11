from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from ystwin.analysis.hybrid_validation import score_external_targets, summarize_scores


KEYS = ["study_id", "case_id", "product", "metric"]
PREDICTION_FIELDS = [
    "predicted",
    "prediction_units",
    "prediction_status",
    "comparable",
    "prediction_basis",
    "reason",
]
ERROR_COLUMNS = [
    "abs_error",
    "relative_error_pct",
    "absolute_percent_error",
    "fold_error",
]
SUMMARY_COLUMNS = [
    "product",
    "metric",
    "units",
    "n_total",
    "n_scored",
    "n_censored",
    "n_unsupported",
    "n_infeasible",
    "n_missing",
    "n_invalid",
    "n_percent_scored",
    "n_fold_scored",
    "mae",
    "median_abs_percent_error",
    "median_fold_error",
]


def _targets(observed=(10.0,)):
    return pd.DataFrame(
        [
            {
                "study_id": "independent_study",
                "case_id": f"case_{i}",
                "product": "glycogen",
                "metric": "intracellular_content",
                "observed": value,
                "units": "mg/g_total_dry_biomass",
                "observation_kind": "quantified",
                "detection_limit": np.nan,
                "biomass_basis": "total_dry_biomass",
                "uncertainty": 0.5,
                "uncertainty_kind": "reported_pm_undefined",
                "source_url": "https://example.org/primary-source",
                "source_location": f"Table 1 row {i + 1}",
                "local_source": f"primary.csv:{i + 2}",
                "limitations": "assay_basis_uncertain;retrospective",
            }
            for i, value in enumerate(observed)
        ]
    )


def _predictions(targets, values=None):
    result = targets[KEYS].copy()
    result["predicted"] = 10.0 if values is None else values
    result["prediction_units"] = targets["units"].to_numpy(copy=True)
    result["prediction_status"] = "predicted"
    result["comparable"] = True
    result["prediction_basis"] = "physiology_conditioned:model_v1"
    result["reason"] = "independent_inputs_only"
    return result


@pytest.fixture
def benchmark_targets():
    directory = Path(__file__).resolve().parents[1] / "data" / "hybrid_benchmark"
    metadata = json.loads((directory / "sources.json").read_text(encoding="utf-8"))
    targets = pd.read_csv(directory / "external_targets.csv")
    assert metadata["csv_schema"]["primary_key"] == KEYS
    assert len(targets) == metadata["coverage"]["total_rows"] == 58
    assert targets["observation_kind"].eq("below_detection_limit").sum() == 8
    assert targets["detection_limit"].isna().all()
    return targets


@pytest.mark.parametrize("frame_name", ["targets", "predictions"])
def test_duplicate_keys_are_rejected(frame_name):
    targets = _targets()
    predictions = _predictions(targets)
    if frame_name == "targets":
        targets = pd.concat([targets, targets], ignore_index=True)
    else:
        predictions = pd.concat([predictions, predictions], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate.*key"):
        score_external_targets(targets, predictions)


@pytest.mark.parametrize("key", KEYS)
def test_prediction_keys_must_identify_known_targets(key):
    targets = _targets()
    predictions = _predictions(targets)
    predictions.loc[0, key] = "unknown"
    with pytest.raises(ValueError, match="unknown.*target"):
        score_external_targets(targets, predictions)


@pytest.mark.parametrize("frame_name", ["targets", "predictions"])
@pytest.mark.parametrize("value", [None, pd.NA, "", " "])
def test_missing_key_values_are_rejected(frame_name, value):
    targets = _targets()
    predictions = _predictions(targets)
    frame = targets if frame_name == "targets" else predictions
    frame.loc[0, "case_id"] = value
    with pytest.raises(ValueError, match="key"):
        score_external_targets(targets, predictions)


@pytest.mark.parametrize("field", PREDICTION_FIELDS)
def test_prediction_schema_is_required(field):
    targets = _targets()
    predictions = _predictions(targets).drop(columns=field)
    with pytest.raises(ValueError, match=field):
        score_external_targets(targets, predictions)


@pytest.mark.parametrize("field", ["observed", "units", "observation_kind", "detection_limit"])
def test_target_measurement_schema_is_required(field):
    targets = _targets()
    predictions = _predictions(targets)
    with pytest.raises(ValueError, match=field):
        score_external_targets(targets.drop(columns=field), predictions)


def test_duplicate_column_names_are_rejected():
    targets = _targets()
    predictions = _predictions(targets)
    predictions = pd.concat([predictions, predictions[["comparable"]]], axis=1)
    with pytest.raises(ValueError, match="duplicate.*column"):
        score_external_targets(targets, predictions)


@pytest.mark.parametrize(
    "units",
    ["mmol/g_total_dry_biomass/h", "mg/g_viable_dry_biomass", "mg/g_total_dry_biomass ", "", None],
)
def test_unit_mismatch_is_retained_without_scoring_or_conversion(units):
    targets = _targets()
    predictions = _predictions(targets, [12.0]).assign(prediction_units=units)
    scores = score_external_targets(targets, predictions)
    assert len(scores) == 1
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "units_mismatch"
    assert scores.loc[0, "predicted"] == 12.0
    assert scores.loc[0, ERROR_COLUMNS].isna().all()
    assert scores.loc[0, "units"] == targets.loc[0, "units"]
    assert summarize_scores(scores).loc[0, "n_invalid"] == 1


def test_wrong_basis_uses_explicit_noncomparability_even_when_units_match():
    targets = _targets()
    predictions = _predictions(targets, [12.0]).assign(
        comparable=False,
        prediction_basis="viable_dry_biomass_PI",
        reason="viable_and_total_biomass_are_not_interchangeable",
    )
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "not_comparable"
    assert scores.loc[0, "prediction_basis"] == "viable_dry_biomass_PI"
    assert scores.loc[0, "reason"] == predictions.loc[0, "reason"]
    assert scores.loc[0, ERROR_COLUMNS].isna().all()


@pytest.mark.parametrize("value", ["True", "False", 1, 0, None, np.nan, pd.NA])
def test_comparable_must_be_an_explicit_boolean_not_truthiness(value):
    targets = _targets()
    predictions = _predictions(targets)
    predictions["comparable"] = pd.Series([value], dtype=object)
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "comparable_must_be_boolean"
    assert scores.loc[0, ERROR_COLUMNS].isna().all()


@pytest.mark.parametrize("value", [True, np.bool_(True)])
def test_boolean_comparability_accepts_opaque_declared_provenance(value):
    targets = _targets()
    predictions = _predictions(targets, [12.0])
    predictions["comparable"] = pd.Series([value], dtype=object)
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "scored"
    assert scores.loc[0, "prediction_basis"] == "physiology_conditioned:model_v1"


@pytest.mark.parametrize("basis", ["", " ", None, pd.NA, 42])
def test_a_comparable_prediction_still_requires_basis_provenance(basis):
    targets = _targets()
    predictions = _predictions(targets).assign(prediction_basis=basis)
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "prediction_basis_required"
    assert scores.loc[0, ERROR_COLUMNS].isna().all()


@pytest.mark.parametrize("value", [-1.0, np.nan, np.inf, -np.inf, None, pd.NA, "bad", True, 1 + 2j])
def test_invalid_predictions_are_preserved_but_never_scored(value):
    targets = _targets()
    predictions = _predictions(targets)
    predictions["predicted"] = pd.Series([value], dtype=object)
    before = predictions.copy(deep=True)
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "invalid"
    assert "prediction" in scores.loc[0, "score_reason"]
    assert scores.loc[0, ERROR_COLUMNS].isna().all()
    assert_frame_equal(predictions, before)
    assert_frame_equal(scores[PREDICTION_FIELDS], predictions[PREDICTION_FIELDS])


@pytest.mark.parametrize("status", ["optimal", "missing", "", None, pd.NA])
def test_unknown_prediction_status_is_invalid_not_a_missing_prediction(status):
    targets = _targets()
    predictions = _predictions(targets).assign(prediction_status=status)
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "invalid_prediction_status"
    assert scores.loc[0, ERROR_COLUMNS].isna().all()
    summary = summarize_scores(scores).iloc[0]
    assert summary.n_invalid == 1
    assert summary.n_missing == 0


@pytest.mark.parametrize("status", ["unsupported", "infeasible"])
def test_unavailable_predictions_remain_named_not_zero_or_invalid(status):
    targets = _targets()
    predictions = _predictions(targets, [np.nan]).assign(
        prediction_status=status,
        comparable=False,
        prediction_basis="not_available",
        reason="missing_independent_capacity" if status == "unsupported" else "solver_infeasible",
    )
    scores = score_external_targets(targets, predictions)
    assert scores.loc[0, "score_status"] == status
    assert scores.loc[0, "reason"] == predictions.loc[0, "reason"]
    assert pd.isna(scores.loc[0, "predicted"])
    assert scores.loc[0, ERROR_COLUMNS].isna().all()
    summary = summarize_scores(scores).iloc[0]
    assert summary[f"n_{status}"] == 1
    assert summary.n_invalid == 0
    assert summary.n_scored == 0
    assert pd.isna(summary.mae)


def test_positive_quantified_scores_are_signed_absolute_and_symmetric_fold():
    targets = _targets([10.0, 10.0, 10.0, 10.0])
    predictions = _predictions(targets, [15.0, 5.0, 10.0, 0.0])
    scores = score_external_targets(targets, predictions)
    assert scores["score_status"].eq("scored").all()
    assert scores["abs_error"].tolist() == pytest.approx([5.0, 5.0, 0.0, 10.0])
    assert scores["relative_error_pct"].tolist() == pytest.approx([50.0, -50.0, 0.0, -100.0])
    assert scores["absolute_percent_error"].tolist() == pytest.approx([50.0, 50.0, 0.0, 100.0])
    assert scores.loc[:2, "fold_error"].tolist() == pytest.approx([1.5, 2.0, 1.0])
    assert pd.isna(scores.loc[3, "fold_error"])


def test_nondetects_are_not_true_zeros_and_never_receive_errors_or_a_pass():
    targets = _targets([0.0, 0.0, 0.0, 0.0])
    targets.loc[:1, "observation_kind"] = "below_detection_limit"
    targets["limitations"] = "nondetect_unknown_limit"
    predictions = _predictions(targets, [0.0, 10000.0, 0.0, 2.0])
    scores = score_external_targets(targets, predictions)
    assert scores["score_status"].tolist() == ["censored", "censored", "scored", "scored"]
    assert scores["is_censored"].tolist() == [True, True, False, False]
    assert scores.loc[:1, "score_reason"].eq("censored_unknown_limit").all()
    assert scores.loc[:1, ERROR_COLUMNS].isna().all().all()
    assert scores.loc[2:, "abs_error"].tolist() == [0.0, 2.0]
    assert scores.loc[2:, ERROR_COLUMNS[1:]].isna().all().all()
    assert "pass" not in scores["score_status"].tolist()
    summary = summarize_scores(scores).iloc[0]
    assert summary.n_scored == 2
    assert summary.n_censored == 2
    assert summary.n_percent_scored == summary.n_fold_scored == 0
    assert summary.mae == 1.0
    assert pd.isna(summary.median_abs_percent_error)
    assert pd.isna(summary.median_fold_error)


def test_known_limits_still_do_not_make_censored_values_point_targets():
    targets = _targets([0.0, 0.0]).assign(
        observation_kind="below_detection_limit", detection_limit=0.1
    )
    scores = score_external_targets(targets, _predictions(targets, [0.05, 0.2]))
    assert scores["score_status"].eq("censored").all()
    assert scores[ERROR_COLUMNS].isna().all().all()
    assert scores["detection_limit"].tolist() == [0.1, 0.1]
    assert scores["predicted"].tolist() == [0.05, 0.2]


@pytest.mark.parametrize("value", [-1.0, np.nan, np.inf, pd.NA, "bad"])
def test_invalid_quantified_observations_are_not_scored(value):
    targets = _targets([value])
    scores = score_external_targets(targets, _predictions(targets))
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "invalid_observed"
    assert scores.loc[0, ERROR_COLUMNS].isna().all()


@pytest.mark.parametrize("kind", ["nondetect", "zero", None, pd.NA])
def test_unknown_observation_kind_does_not_default_to_quantified(kind):
    targets = _targets().assign(observation_kind=kind)
    scores = score_external_targets(targets, _predictions(targets))
    assert scores.loc[0, "score_status"] == "invalid"
    assert scores.loc[0, "score_reason"] == "invalid_observation_kind"
    assert scores.loc[0, ERROR_COLUMNS].isna().all()


def test_all_benchmark_rows_and_provenance_survive_a_partial_prediction_join(benchmark_targets):
    targets = benchmark_targets
    predictions = _predictions(targets.iloc[[0, 1, 2, 3, 8, 9]], [10.0, np.nan, np.nan, -1.0, 0.0, np.nan])
    predictions["prediction_status"] = ["predicted", "unsupported", "infeasible", "predicted", "predicted", "infeasible"]
    predictions["comparable"] = [True, False, False, True, True, False]
    predictions = predictions.iloc[::-1]
    scores = score_external_targets(targets, predictions)
    assert_frame_equal(scores[targets.columns], targets)
    assert len(scores) == 58
    assert scores["score_status"].value_counts().to_dict() == {
        "missing": 52,
        "infeasible": 2,
        "scored": 1,
        "unsupported": 1,
        "invalid": 1,
        "censored": 1,
    }
    assert scores["is_censored"].sum() == 8
    unavailable = scores["score_status"].isin(["missing", "unsupported", "infeasible"])
    assert scores.loc[unavailable, "predicted"].isna().all()
    assert scores.loc[unavailable, ERROR_COLUMNS].isna().all().all()
    summary = summarize_scores(scores)
    assert summary["n_total"].sum() == 58
    assert summary["n_censored"].sum() == 8
    assert summary["n_missing"].sum() == 52
    assert summary["n_infeasible"].sum() == 2
    assert summary["n_unsupported"].sum() == 1
    assert summary["n_invalid"].sum() == 1
    assert summary["n_scored"].sum() == 1
    assert "glutathione" not in summary["product"].tolist()


def test_all_missing_predictions_keep_the_full_censored_and_product_census(benchmark_targets):
    predictions = _predictions(benchmark_targets.iloc[:0])
    scores = score_external_targets(benchmark_targets, predictions)
    assert scores["score_status"].eq("missing").all()
    assert scores["predicted"].isna().all()
    summary = summarize_scores(scores)
    assert summary["n_total"].sum() == summary["n_missing"].sum() == 58
    assert summary["n_censored"].sum() == 8
    assert len(summary) == 5
    assert summary["n_scored"].sum() == 0
    assert summary[["mae", "median_abs_percent_error", "median_fold_error"]].isna().all().all()


def test_prediction_metadata_does_not_overwrite_target_provenance():
    targets = _targets()
    predictions = _predictions(targets).assign(
        source_url="https://example.org/model-inputs",
        limitations="preexisting_model_priors",
        run_id="fixed_run",
    )
    scores = score_external_targets(targets, predictions)
    assert_frame_equal(scores[targets.columns], targets)
    assert scores.loc[0, "source_url_prediction"] == predictions.loc[0, "source_url"]
    assert scores.loc[0, "limitations_prediction"] == "preexisting_model_priors"
    assert scores.loc[0, "run_id"] == "fixed_run"


@pytest.mark.parametrize("frame_name", ["targets", "predictions"])
def test_supplied_score_columns_cannot_overwrite_computed_results(frame_name):
    targets = _targets()
    predictions = _predictions(targets)
    frame = targets if frame_name == "targets" else predictions
    frame["abs_error"] = 0.0
    with pytest.raises(ValueError, match="reserved"):
        score_external_targets(targets, predictions)


def test_scoring_and_summary_do_not_mutate_inputs_or_share_writable_columns():
    targets = _targets([10.0, 20.0])
    targets.index = pd.Index([7, 7], name="source_row")
    predictions = _predictions(targets, [5.0, 15.0])
    before_targets = targets.copy(deep=True)
    before_predictions = predictions.copy(deep=True)
    scores = score_external_targets(targets, predictions)
    assert_frame_equal(scores[targets.columns], targets)
    before_scores = scores.copy(deep=True)
    summarize_scores(scores)
    assert_frame_equal(scores, before_scores)
    scores.loc[:, "observed"] = 999.0
    scores.loc[:, "predicted"] = 999.0
    scores.loc[:, "limitations"] = "changed"
    assert_frame_equal(targets, before_targets)
    assert_frame_equal(predictions, before_predictions)


def test_scoring_is_deterministic_row_local_and_does_not_fit_predictions():
    targets = _targets([10.0, 20.0, 30.0])
    predictions = _predictions(targets, [3.0, 4.0, 5.0])
    before_predictions = predictions.copy(deep=True)
    baseline = score_external_targets(targets, predictions)
    assert_frame_equal(baseline, score_external_targets(targets, predictions))
    assert_frame_equal(baseline, score_external_targets(targets, predictions.iloc[::-1]))
    reverse = score_external_targets(targets.iloc[::-1], predictions)
    assert_frame_equal(reverse.sort_values(KEYS), baseline.sort_values(KEYS))
    one_row = score_external_targets(targets.iloc[:1], predictions.iloc[:1])
    assert_frame_equal(one_row, baseline.iloc[:1])
    changed_targets = targets.copy(deep=True)
    changed_targets["observed"] = [1000.0, 0.0, 1000000.0]
    changed = score_external_targets(changed_targets, predictions)
    assert changed["predicted"].tolist() == baseline["predicted"].tolist() == [3.0, 4.0, 5.0]
    assert changed["abs_error"].tolist() == [997.0, 4.0, 999995.0]
    assert_frame_equal(predictions, before_predictions)


def test_summary_uses_metric_specific_denominators_and_no_failed_row_imputation():
    targets = _targets([10.0, 20.0, 5.0, 0.0, 0.0, 0.0, 9.0, 9.0, 9.0, 9.0])
    targets.loc[5, "observation_kind"] = "below_detection_limit"
    predictions = _predictions(targets.iloc[:9], [15.0, 10.0, 0.0, 2.0, 0.0, 0.0, np.nan, np.nan, -1.0])
    predictions.loc[6, "prediction_status"] = "unsupported"
    predictions.loc[7, "prediction_status"] = "infeasible"
    scores = score_external_targets(targets, predictions)
    summary = summarize_scores(scores)
    assert summary.columns.tolist() == SUMMARY_COLUMNS
    row = summary.iloc[0]
    assert row.n_total == 10
    assert row.n_scored == 5
    assert row.n_censored == row.n_unsupported == row.n_infeasible == row.n_invalid == row.n_missing == 1
    assert row.n_percent_scored == 3
    assert row.n_fold_scored == 2
    assert row.mae == pytest.approx(22.0 / 5.0)
    assert row.median_abs_percent_error == 50.0
    assert row.median_fold_error == 1.75
    assert row.units == "mg/g_total_dry_biomass"
    before = scores.copy(deep=True)
    scores.loc[scores["score_status"].ne("scored"), ERROR_COLUMNS] = 9999.0
    assert_frame_equal(summarize_scores(scores), summarize_scores(before))


def test_summary_keeps_products_and_metrics_separate_and_ignores_unused_categories():
    targets = _targets([10.0, 2.0, 4.0])
    targets.loc[1:, "product"] = "squalene"
    targets.loc[2, "metric"] = "culture_titre"
    targets.loc[2, "units"] = "mg/L_culture"
    targets["product"] = pd.Categorical(targets["product"], categories=["squalene", "glycogen", "glutathione"])
    scores = score_external_targets(targets, _predictions(targets, [15.0, 3.0, 12.0]))
    summary = summarize_scores(scores).set_index(["product", "metric"])
    assert len(summary) == 3
    assert summary["n_total"].eq(1).all()
    assert summary.loc[("glycogen", "intracellular_content"), "mae"] == 5.0
    assert summary.loc[("squalene", "intracellular_content"), "mae"] == 1.0
    assert summary.loc[("squalene", "culture_titre"), "mae"] == 8.0


def test_summary_refuses_to_average_absolute_errors_in_mixed_units():
    targets = _targets([10.0, 10.0])
    targets.loc[1, "units"] = "g/g_total_dry_biomass"
    scores = score_external_targets(targets, _predictions(targets, [15.0, 15.0]))
    with pytest.raises(ValueError, match="mixed units"):
        summarize_scores(scores)


def test_summary_rejects_duplicate_target_keys():
    targets = _targets()
    scores = score_external_targets(targets, _predictions(targets))
    with pytest.raises(ValueError, match="duplicate.*key"):
        summarize_scores(pd.concat([scores, scores], ignore_index=True))


def test_empty_target_tables_return_stable_schemas_without_synthetic_rows():
    targets = _targets().iloc[:0]
    scores = score_external_targets(targets, _predictions(targets))
    assert scores.empty
    assert set([*targets.columns, *PREDICTION_FIELDS, *ERROR_COLUMNS, "score_status", "score_reason", "is_censored"]).issubset(scores.columns)
    summary = summarize_scores(scores)
    assert summary.empty
    assert summary.columns.tolist() == SUMMARY_COLUMNS


@pytest.mark.parametrize("observed,predicted", [(1.0, 1e306), (1e308, 1.0)])
def test_large_finite_errors_have_finite_summary_medians(observed, predicted):
    targets = _targets([observed, observed])
    scores = score_external_targets(targets, _predictions(targets, [predicted, predicted]))
    assert np.isfinite(scores[ERROR_COLUMNS]).all().all()
    summary = summarize_scores(scores).iloc[0]
    assert summary.n_scored == summary.n_percent_scored == summary.n_fold_scored == 2
    assert summary.mae == pytest.approx(scores.loc[0, "abs_error"])
    assert summary.median_abs_percent_error == pytest.approx(scores.loc[0, "absolute_percent_error"])
    assert summary.median_fold_error == pytest.approx(scores.loc[0, "fold_error"])


def test_overflowing_ratios_from_finite_inputs_are_not_dropped_as_easy_to_ignore_errors():
    targets = _targets([1e-300, 10.0])
    scores = score_external_targets(targets, _predictions(targets, [1e300, 20.0]))
    assert scores["score_status"].eq("scored").all()
    assert scores.loc[0, "absolute_percent_error"] == np.inf
    assert scores.loc[0, "fold_error"] == np.inf
    summary = summarize_scores(scores).iloc[0]
    assert summary.n_scored == summary.n_percent_scored == summary.n_fold_scored == 2
    assert summary.mae == pytest.approx(5e299)
    assert summary.median_abs_percent_error == np.inf
    assert summary.median_fold_error == np.inf
