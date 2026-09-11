from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis import native_physiology as native


SOURCE_PATH = (
    Path(__file__).resolve().parents[1]
    / "data/physiology/chemostatData_VanHoek1998.tsv"
)
SOURCE_SHA256 = "4c01b3465ea1410c86709568da44a32b1de652c21a88de623bbf566838506117"
FIELDS = (
    "GlucoseUptake", "O2uptake", "CO2production", "Ethanol", "Glycerol",
    "Acetate", "Pyruvate", "BiomassYield_g_per_g", "CarbonRecovery_pct",
)
PUBLISHED = np.array([
    [.025, .30, .80, .80, 0, 0, 0, 0, .45, 98.9],
    [.05, .60, 1.30, 1.40, 0, 0, 0, 0, .47, 95.0],
    [.10, 1.10, 2.50, 2.70, 0, 0, 0, 0, .48, 96.0],
    [.15, 1.70, 3.90, 4.20, 0, 0, 0, 0, .49, 102.4],
    [.20, 2.30, 5.30, 5.70, 0, 0, 0, 0, .48, 100.9],
    [.25, 2.80, 7.00, 7.50, 0, 0, 0, 0, .48, 102.6],
    [.28, 3.40, 7.40, 8.00, .11, 0, .08, .01, .46, 97.0],
    [.30, 4.50, 6.10, 8.80, 2.30, 0, .41, .01, .37, 99.1],
    [.35, 8.60, 5.10, 14.90, 9.50, .05, .62, .03, .23, 99.4],
    [.40, 11.10, 3.70, 18.90, 13.90, .15, .60, .05, .20, 97.9],
])


@pytest.fixture
def dataset():
    return native.load_native_chemostat_data()


@pytest.fixture
def model(dataset):
    return native.fit_native_exchange_model(dataset.split()["train"])


def prediction(model, growth, field):
    return next(row for row in model.predict([growth]) if row["observable_id"] == field)


def test_exact_published_table_and_structured_provenance(dataset):
    rows = dataset.observations
    assert len(rows) == 90
    assert tuple(rows.observable_id.unique()) == FIELDS
    wide = rows.pivot(index="growth_rate_per_h", columns="observable_id", values="reported_value")
    np.testing.assert_array_equal(wide.index, PUBLISHED[:, 0])
    np.testing.assert_allclose(wide.loc[:, list(FIELDS)], PUBLISHED[:, 1:], rtol=0, atol=0)
    assert rows.source_table_sha256.eq(SOURCE_SHA256).all()
    assert hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest() == SOURCE_SHA256
    assert rows.source_doi.eq("10.1128/AEM.64.11.4226-4233.1998").all()
    assert rows.source_column.equals(rows.observable_id)
    assert set(rows.source_line) == set(range(2, 12))
    source = dataset.metadata["source"]
    assert source["table_sha256"] == SOURCE_SHA256
    assert source["pmid"] == "9797269"
    assert source["table"] == "Table 1"
    assert source["table_page"] == 4227
    assert source["conditions"] == {
        "organism": "Saccharomyces cerevisiae", "strain": "DS28911",
        "cultivation": "aerobic glucose-limited chemostat", "medium": "defined mineral medium",
        "temperature_c": 30.0, "ph": 5.0,
    }
    assert source["biomass_carbon"] == {
        "mass_fraction": .48, "status": "source_assumption_for_carbon_recovery",
    }
    json.dumps(dataset.metadata, allow_nan=False)


def test_custom_path_and_byte_identity_guard(dataset, tmp_path):
    path = tmp_path / "native.tsv"
    path.write_bytes(SOURCE_PATH.read_bytes())
    pd.testing.assert_frame_equal(
        native.load_native_chemostat_data(str(path)).observations, dataset.observations,
    )
    path.write_bytes(path.read_bytes().replace(b"0.30", b"0.31", 1))
    with pytest.raises(ValueError, match="SHA-256"):
        native.load_native_chemostat_data(path)


def test_source_zeros_are_censored_without_a_numeric_detection_limit(dataset):
    rows = dataset.observations
    zero = rows.reported_value.eq(0)
    assert zero.sum() == 26
    assert rows.loc[zero, "observation_status"].eq("below_detection_limit").all()
    assert rows.loc[zero, "value"].isna().all()
    assert rows.loc[~zero, "observation_status"].eq("quantified").all()
    assert rows.detection_limit.isna().all()
    assert rows.loc[rows.observable_id.eq("GlucoseUptake"), "unit"].eq("mmol/gDW/h").all()
    assert rows.loc[rows.observable_id.eq("BiomassYield_g_per_g"), "unit"].eq("gDW/g_glucose").all()
    assert rows.loc[rows.observable_id.eq("CarbonRecovery_pct"), "unit"].eq("%").all()


def test_predeclared_condition_split_is_exclusive_and_complete(dataset):
    splits = dataset.split()
    assert tuple(native.TRAIN_GROWTH_RATES_PER_H) == (.025, .05, .10, .20, .25, .30, .40)
    assert tuple(native.HOLDOUT_GROWTH_RATES_PER_H) == (.15, .28, .35)
    assert set(splits) == {"train", "holdout"}
    assert len(splits["train"]) == 63
    assert len(splits["holdout"]) == 27
    train = set(splits["train"].growth_rate_per_h)
    held = set(splits["holdout"].growth_rate_per_h)
    assert train.isdisjoint(held)
    assert train | held == set(PUBLISHED[:, 0])
    with pytest.raises(ValueError, match="holdout"):
        native.fit_native_exchange_model(dataset.observations)


def test_quantified_interpolation_uses_adjacent_training_knots(model):
    expected = [(.15, "GlucoseUptake", 1.7), (.28, "GlucoseUptake", 3.82),
                (.35, "GlucoseUptake", 7.8), (.28, "O2uptake", 6.46),
                (.35, "Ethanol", 8.1), (.35, "Pyruvate", .03)]
    for growth, field, value in expected:
        row = prediction(model, growth, field)
        assert row["prediction_status"] == "quantified"
        assert row["prediction"] == pytest.approx(value)
        assert row["method"] == "linear_interpolation"
    exact = prediction(model, .40, "Glycerol")
    assert exact["prediction_status"] == "quantified"
    assert exact["prediction"] == .15
    assert exact["method"] == "observed_knot"


@pytest.mark.parametrize("growth,field,support", [
    (.15, "Glycerol", [.10, .20]), (.28, "Ethanol", [.25, .30]),
    (.35, "Glycerol", [.30, .40]), (.30, "Glycerol", [.30]),
])
def test_censored_support_is_unquantified_not_zero(model, growth, field, support):
    row = prediction(model, growth, field)
    assert row["prediction_status"] == "unquantified_censored_support"
    assert row["prediction"] is None
    assert row["bounds"] is None
    assert row["detection_limit"] is None
    assert [point["growth_rate_per_h"] for point in row["support"]] == support
    assert "below_detection_limit" in {point["observation_status"] for point in row["support"]}
    json.dumps(row, allow_nan=False)


def test_no_extrapolation_and_all_requested_conditions_are_retained(model):
    rows = model.predict([-.1, 0, .024, .025, .40, .401])
    assert len(rows) == 54
    for row in rows:
        if row["growth_rate_per_h"] not in (.025, .40):
            assert row["prediction_status"] == "unsupported_growth_rate"
            assert row["prediction"] is None
            assert row["support"] == []
    assert len(model.predict([.15, .15])) == 18
    assert model.predict([]) == []
    json.dumps(rows, allow_nan=False)


@pytest.mark.parametrize("rates", [[float("nan")], [float("inf")], [True], [[.1]], [".1"]])
def test_invalid_prediction_conditions_are_rejected(model, rates):
    with pytest.raises(ValueError, match="growth"):
        model.predict(rates)


def test_fit_never_reads_hidden_rows_and_domain_comes_from_supplied_rows(dataset, monkeypatch):
    rows = dataset.split()["train"]
    subset = rows.loc[rows.growth_rate_per_h.isin([.10, .20])].copy()

    def forbidden(*args, **kwargs):
        raise AssertionError("Fitting/prediction/reload must not reload the source table")

    monkeypatch.setattr(native, "load_native_chemostat_data", forbidden)
    fitted = native.fit_native_exchange_model(subset)
    assert fitted.to_dict()["training_growth_rates_per_h"] == [.10, .20]
    assert prediction(fitted, .15, "GlucoseUptake")["prediction"] == pytest.approx(1.7)
    assert prediction(fitted, .05, "GlucoseUptake")["prediction_status"] == "unsupported_growth_rate"
    restored = native.NativeExchangeModel.from_dict(fitted.to_dict())
    assert restored.predict([.15]) == fitted.predict([.15])


def test_heldout_label_changes_cannot_change_fitted_model(dataset, model):
    changed = deepcopy(dataset)
    mask = changed.observations.split.eq("holdout") & changed.observations.observation_status.eq("quantified")
    changed.observations.loc[mask, ["value", "reported_value"]] *= 2
    refit = native.fit_native_exchange_model(changed.split()["train"])
    assert refit.to_dict() == model.to_dict()


def test_actual_training_values_are_fitted_and_persisted(dataset, model):
    rows = dataset.split()["train"]
    mask = rows.growth_rate_per_h.eq(.20) & rows.observable_id.eq("GlucoseUptake")
    rows.loc[mask, ["value", "reported_value"]] = 2.9
    changed = native.fit_native_exchange_model(rows)
    assert prediction(changed, .15, "GlucoseUptake")["prediction"] == pytest.approx(2.0)
    assert changed.to_dict()["training_rows_sha256"] != model.to_dict()["training_rows_sha256"]
    assert any(k["value"] == 2.9 for k in changed.to_dict()["knots"]["GlucoseUptake"])


def test_row_id_renaming_row_order_and_dataframe_index_do_not_change_fit(dataset, model):
    rows = dataset.split()["train"].sample(frac=1, random_state=7).reset_index(drop=True)
    rows["row_id"] = [f"arbitrary-{i}" for i in range(len(rows))]
    rows.index = range(100, 100 + len(rows))
    assert native.fit_native_exchange_model(rows).to_dict() == model.to_dict()


def test_json_serialization_and_reload_are_complete_and_detached(model):
    payload = json.loads(json.dumps(model.to_dict(), allow_nan=False))
    assert payload["model_kind"] == "empirical_growth_dependent_exchange_law"
    assert payload["input"] == {
        "name": "growth_rate_per_h", "unit": "1/h", "role": "imposed_condition_not_scored",
    }
    assert payload["extrapolation"] == "unsupported"
    assert payload["source"]["table_sha256"] == SOURCE_SHA256
    assert set(payload["knots"]) == set(FIELDS)
    assert all(len(knots) == 7 for knots in payload["knots"].values())
    assert payload["knots"]["Glycerol"][0]["value"] is None
    assert payload["knots"]["Glycerol"][0]["reported_value"] == 0
    restored = native.NativeExchangeModel.from_dict(payload)
    assert restored.to_dict() == payload
    assert restored.predict([.01, .15, .28, .35, .40]) == model.predict([.01, .15, .28, .35, .40])
    payload["knots"]["GlucoseUptake"][0]["value"] = 123
    assert model.to_dict()["knots"]["GlucoseUptake"][0]["value"] == .30
    with pytest.raises(ValueError, match="payload|quantified"):
        native.NativeExchangeModel.from_dict(payload)


def test_validation_keeps_every_readout_and_reports_denominators(dataset, model):
    result = native.evaluate_native_exchange_model(model, dataset.split()["holdout"])
    assert result["n_conditions"] == 3
    assert result["n_readouts"] == len(result["readouts"]) == 27
    assert result["n_scored"] == 18
    assert result["coverage"] == pytest.approx(2 / 3)
    assert set(result["metrics"]) == set(FIELDS)
    assert "growth_rate_per_h" not in result["metrics"]
    glycerol = result["metrics"]["Glycerol"]
    assert glycerol["n_total"] == 3
    assert glycerol["n_observed_quantified"] == 1
    assert glycerol["n_observed_censored"] == 2
    assert glycerol["n_scored"] == 0
    assert glycerol["coverage"] == 0
    assert glycerol["mae"] is None and glycerol["rmse"] is None
    assert glycerol["prediction_status_counts"]["unquantified_censored_support"] == 3
    ethanol = result["metrics"]["Ethanol"]
    assert ethanol["n_observed_quantified"] == 2
    assert ethanol["n_scored"] == 1
    assert ethanol["quantified_observation_coverage"] == .5
    assert ethanol["mae"] == pytest.approx(1.4)
    assert result["metrics"]["GlucoseUptake"]["mae"] == pytest.approx((.42 + .8) / 3)
    assert {row["growth_rate_per_h"] for row in result["readouts"]} == {.15, .28, .35}
    assert all(row["residual"] is None for row in result["readouts"] if not row["scored"])
    json.dumps(result, allow_nan=False)


def test_validation_retains_unsupported_predictions(dataset):
    train = dataset.split()["train"]
    model = native.fit_native_exchange_model(train.loc[train.growth_rate_per_h.isin([.1, .2])])
    result = native.evaluate_native_exchange_model(model, dataset.split()["holdout"])
    assert result["n_readouts"] == 27
    assert result["n_scored"] == 5
    assert sum(row["prediction_status"] == "unsupported_growth_rate" for row in result["readouts"]) == 18


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "status", "zero", "lod", "unit", "source", "nan"])
def test_malformed_training_observations_are_not_silently_dropped(dataset, mutation):
    rows = dataset.split()["train"].copy()
    index = rows.index[0]
    if mutation == "missing":
        rows = rows.drop(index)
    elif mutation == "duplicate":
        rows = pd.concat([rows, rows.iloc[[0]]])
    elif mutation == "status":
        rows.loc[index, "observation_status"] = "unknown"
    elif mutation == "zero":
        rows.loc[index, ["value", "reported_value"]] = 0
    elif mutation == "lod":
        rows.loc[index, "detection_limit"] = .01
    elif mutation == "unit":
        rows.loc[index, "unit"] = "mol/L"
    elif mutation == "source":
        rows.loc[index, "source_table_sha256"] = "unverified"
    else:
        rows.loc[index, "value"] = np.nan
    with pytest.raises(ValueError):
        native.fit_native_exchange_model(rows)
