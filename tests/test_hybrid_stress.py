import json
from dataclasses import FrozenInstanceError

import numpy as np
import pandas as pd
import pytest

from ystwin import paths, reporter
from ystwin.analysis import hybrid_stress
from ystwin.analysis.hybrid_stress import (
    CalibratedSensorModel,
    StressEstimate,
    fit_sensor_model,
)


CONSTRUCTS = ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1")
TRAINING_PLATES = ("train-a", "train-b")
FEATURES = tuple(f"{name}_activity_fold" for name in CONSTRUCTS)


@pytest.fixture(scope="module")
def readings():
    rows = []
    for plate, gain, growth_scale in (
        ("train-a", 1.0, 1.0),
        ("train-b", 1.6, 1.1),
        ("held-out", 3.0, 0.85),
    ):
        for index, construct in enumerate(CONSTRUCTS):
            basal = (800.0 + 200.0 * index) * gain
            control_mu = (0.3 + 0.02 * index) * growth_scale
            for dose in (0.0, 0.1, 0.2, 0.5, 1.0, 2.0):
                retention = 1.0 - 0.35 * dose if dose <= 1 else 0.25
                fold = 1 + (0.7 + index * 0.2) * dose / (0.4 + dose)
                if dose > 1:
                    fold = 0.3
                for replicate, offset in enumerate((0.98, 1.0, 1.02)):
                    activity = basal * fold * offset
                    rows.append({
                        "plate": plate,
                        "construct": construct,
                        "stressor": "DTT" if index < 2 else "H2O2",
                        "dose_mM": dose,
                        "well": f"{construct}-{dose}-{replicate}",
                        "mu_late": control_mu * retention,
                        "mu_max": control_mu * retention * 1.5,
                        "n_points": 25,
                        "naive_late": activity / (control_mu * retention),
                        "activity_late": activity,
                        "activity_late_percell": activity,
                        "mu_late_se": 0.002,
                        "activity_late_sigma": 0.05 * basal,
                        "reporter_cv": 0.01,
                    })
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def model(readings):
    return fit_sensor_model(readings, training_plates=TRAINING_PLATES, seed=7)


def test_public_estimate_is_frozen_and_names_only_supported_axes():
    estimate = StressEstimate(0.2, 0.4, 0.8, "observed", False)
    assert estimate.modules == {"UPR": 0.2, "oxidative": 0.4}
    with pytest.raises(FrozenInstanceError):
        estimate.upr = 0.9


def test_training_and_condition_prediction_are_deterministic(readings, model):
    again = fit_sensor_model(readings, training_plates=TRAINING_PLATES, seed=7)
    assert isinstance(model, CalibratedSensorModel)
    assert again.metadata == model.metadata
    assert again.predict_condition("DTT", 0.5, seed=19) == model.predict_condition(
        "DTT", 0.5, seed=19
    )
    assert model.predict_condition("H2O2", 0.5, seed=11) == model.predict_condition(
        "H2O2", 0.5, seed=11
    )
    json.dumps(model.metadata, allow_nan=False)


def test_held_out_mutation_cannot_change_any_fitted_coefficient(readings, model):
    changed = readings.copy()
    held = changed.plate == "held-out"
    changed.loc[held, ["activity_late", "mu_late", "n_points", "naive_late"]] = np.nan
    changed.loc[held, "dose_mM"] = -1000.0
    changed.loc[held, "construct"] = "not-a-sensor"
    changed.loc[held, "stressor"] = "not-a-stressor"
    fitted = fit_sensor_model(changed, training_plates=TRAINING_PLATES, seed=7)
    assert fitted.metadata == model.metadata
    assert fitted.predict_condition("DTT", 0.8) == model.predict_condition("DTT", 0.8)


def test_only_allowlisted_sensor_columns_are_used_without_external_reads(
    readings, model, monkeypatch
):
    def forbidden(*args, **kwargs):
        raise AssertionError("sensor fitting must not read another data source")

    monkeypatch.setattr(pd, "read_csv", forbidden)
    extended = readings.assign(product_titre=object(), product_identity=object())
    fitted = fit_sensor_model(extended, training_plates=TRAINING_PLATES, seed=7)
    assert fitted.metadata == model.metadata
    assert fitted.predict_condition("H2O2", 0.5) == model.predict_condition("H2O2", 0.5)
    assert "product_titre" not in fitted.metadata["measurement_columns_used"]
    assert fitted.metadata["training_plates"] == list(TRAINING_PLATES)
    assert fitted.metadata["evidence_status"]["product_measurements_used"] is False


def test_schema_is_named_and_channel_order_is_not_positional(model):
    assert model.feature_names == FEATURES
    assert set(model.metadata["observation_schema"]) == set(FEATURES)
    for description in model.metadata["observation_schema"].values():
        assert description["units"] == "dimensionless"
        assert "control" in description["definition"]
    features = dict(zip(FEATURES, (1.1, 1.2, 1.0, 1.0)))
    reversed_features = dict(reversed(list(features.items())))
    assert model.infer(features) == model.infer(reversed_features)
    assert model.infer(features, source="plate measurement").source == "plate measurement"


@pytest.mark.parametrize("change", ["missing", "extra", "negative", "zero", "nan", "inf", "array"])
def test_inference_refuses_wrong_or_invalid_channels(model, change):
    features = dict.fromkeys(FEATURES, 1.0)
    if change == "missing":
        features.pop(FEATURES[0])
    elif change == "extra":
        features["raw_RFU"] = 12.0
    else:
        features[FEATURES[0]] = {
            "negative": -1.0, "zero": 0.0, "nan": np.nan,
            "inf": np.inf, "array": [1.0, 2.0],
        }[change]
    with pytest.raises((ValueError, TypeError), match="feature|channel|positive|finite|scalar"):
        model.infer(features)


def test_inference_refuses_an_unnamed_vector(model):
    with pytest.raises((ValueError, TypeError), match="Mapping|mapping|named"):
        model.infer(np.ones(4))


@pytest.mark.parametrize("stressor,dose", [("NaCl", 0.0), ("DTT", -0.1), ("DTT", np.nan), ("H2O2", np.inf), (None, 0.5)])
def test_invalid_conditions_are_rejected(model, stressor, dose):
    with pytest.raises(ValueError, match="stressor|dose"):
        model.predict_condition(stressor, dose)


def test_zero_stress_is_anchored_to_measured_controls(model):
    for stressor in (None, "DTT", "H2O2"):
        prediction = model.predict_condition(stressor)
        assert prediction.upr == pytest.approx(0.0, abs=1e-7)
        assert prediction.oxidative == pytest.approx(0.0, abs=1e-7)
        assert prediction.growth_retention == pytest.approx(1.0, abs=1e-7)
        assert "SIMULATED" in prediction.source
        assert not prediction.extrapolated
    observed = model.infer(dict.fromkeys(FEATURES, 1.0))
    assert observed.modules == pytest.approx({"UPR": 0.0, "oxidative": 0.0}, abs=1e-7)
    assert observed.growth_retention == pytest.approx(1.0, abs=1e-7)


def test_growth_effect_is_calibrated_from_growth_phenotypes(readings, model):
    stressed = model.predict_condition("DTT", 0.5)
    assert stressed.upr > 0.1
    assert stressed.oxidative < 0.01
    assert 0.0 < stressed.growth_retention < 0.95
    unaffected = readings.copy()
    for _, group in unaffected.groupby(["plate", "construct"]):
        control = group.loc[group.dose_mM == 0, "mu_late"].mean()
        indices = group.index[group.dose_mM <= 1]
        unaffected.loc[indices, "mu_late"] = control
    refit = fit_sensor_model(unaffected, training_plates=TRAINING_PLATES, seed=7)
    assert refit.predict_condition("DTT", 0.5).growth_retention == pytest.approx(1.0)


def test_supported_doses_and_observation_extrapolation_are_explicit(model):
    for stressor in ("DTT", "H2O2"):
        assert model.metadata["calibrated_ranges"][stressor]["dose_mM"] == [0.0, 1.0]
        assert not model.predict_condition(stressor, 0.5).extrapolated
        assert not model.predict_condition(stressor, 1.0).extrapolated
        assert model.predict_condition(stressor, 2.0).extrapolated
    extreme = model.infer(dict.fromkeys(FEATURES, 100.0))
    assert extreme.extrapolated
    assert 0.0 <= extreme.upr <= 1.0
    assert 0.0 <= extreme.oxidative <= 1.0
    assert 0.0 <= extreme.growth_retention <= 1.0


def test_existing_reporter_ode_generates_training_and_prediction_features(readings, monkeypatch):
    original = reporter.simulate_reporter
    calls = []

    def tracked(times_h, k_synth, growth_rate, *args, **kwargs):
        calls.append((np.asarray(k_synth).copy(), np.asarray(growth_rate).copy()))
        return original(times_h, k_synth, growth_rate, *args, **kwargs)

    monkeypatch.setattr(reporter, "simulate_reporter", tracked)
    fitted = fit_sensor_model(readings, training_plates=TRAINING_PLATES, seed=3)
    assert len(calls) >= 4 * fitted.metadata["synthetic_training_size"]
    assert np.ptp([mu[0] for _, mu in calls]) > 0.02
    assert min(activity[0] for activity, _ in calls) > 500.0
    count = len(calls)
    fitted.predict_condition("DTT", 0.5)
    assert len(calls) >= count + 4

    def dark(times_h, *args, **kwargs):
        return np.zeros_like(times_h)

    monkeypatch.setattr(reporter, "simulate_reporter", dark)
    with pytest.raises(ValueError, match="positive|ODE|activity"):
        fitted.predict_condition("DTT", 0.5)


def test_negative_and_missing_measurements_are_recorded_not_zero_filled(readings):
    changed = readings.copy()
    high = (changed.plate == "train-a") & (changed.dose_mM == 2.0)
    changed.loc[high & (changed.construct == "UPRE1"), "activity_late"] = -10.0
    changed.loc[high & (changed.construct == "UPRE2"), "mu_late"] = np.nan
    fitted = fit_sensor_model(changed, training_plates=TRAINING_PLATES)
    exclusions = fitted.metadata["exclusions"]
    assert any(row["construct"] == "UPRE1" and "activity" in str(row["reasons"])
               for row in exclusions)
    assert any(row["construct"] == "UPRE2" and "growth" in str(row["reasons"])
               for row in exclusions)
    assert fitted.metadata["calibrated_ranges"]["DTT"]["dose_mM"] == [0.0, 1.0]
    table = fitted.validation_table(changed[high])
    assert table.activity_late_measured.min() < 0
    assert not table.activity_valid.all()
    assert not table.growth_valid.all()


@pytest.mark.parametrize("count", [np.nan, np.inf])
def test_nonfinite_timepoint_counts_are_recorded_as_exclusions(readings, count):
    changed = readings.copy()
    bad = (changed.plate == "train-a") & (changed.construct == "UPRE1") & (changed.dose_mM == 2)
    changed["n_points"] = changed.n_points.astype(float)
    changed.loc[bad, "n_points"] = count
    fitted = fit_sensor_model(changed, training_plates=TRAINING_PLATES)
    assert any("insufficient_or_missing_timepoints" in row["reasons"]
               and not row["activity_used"] and not row["growth_used"]
               for row in fitted.metadata["exclusions"])
    json.dumps(fitted.metadata, allow_nan=False)


def test_a_negative_well_is_not_hidden_by_a_positive_condition_mean(readings, model):
    changed = readings[readings.plate == "held-out"].copy()
    bad = changed.index[(changed.construct == "UPRE1") & (changed.dose_mM == 0.2)][0]
    changed.loc[bad, "activity_late"] = -1.0
    table = model.validation_table(changed)
    condition = table[(table.construct == "UPRE1") & (table.dose_mM == 0.2)].iloc[0]
    assert condition.activity_late_measured > 0
    assert not condition.activity_valid
    assert np.isnan(condition.activity_fold_measured)
    assert "activity_nonpositive_or_missing" in condition.exclusion_reason


@pytest.mark.parametrize("plates", [("absent",), (), ("train-a", "train-a")])
def test_training_plate_selection_rejects_absent_empty_or_duplicate_ids(readings, plates):
    with pytest.raises(ValueError, match="plate"):
        fit_sensor_model(readings, training_plates=plates)


def test_controls_are_required_and_stressor_construct_pairing_is_checked(readings):
    single = readings[readings.plate == "train-a"].copy()
    single = single[~((single.construct == "UPRE1") & (single.dose_mM == 0))]
    with pytest.raises(ValueError, match="control|paired|UPR"):
        fit_sensor_model(single)
    wrong = readings.copy()
    wrong.loc[wrong.construct == "UPRE1", "stressor"] = "H2O2"
    with pytest.raises(ValueError, match="stressor|pair"):
        fit_sensor_model(wrong)
    with pytest.raises(ValueError, match="column"):
        fit_sensor_model(readings.drop(columns="activity_late"))


def test_plate_names_do_not_decide_health_or_pairing(readings, model):
    renamed = readings.copy()
    renamed["plate"] = renamed.plate.map({
        "train-a": "bad filename", "train-b": "20260101 perfectly healthy", "held-out": "other",
    })
    fitted = fit_sensor_model(
        renamed, training_plates=("bad filename", "20260101 perfectly healthy"), seed=7
    )
    expected = model.predict_condition("DTT", 0.5)
    observed = fitted.predict_condition("DTT", 0.5)
    assert observed.modules == pytest.approx(expected.modules, abs=1e-6)
    assert observed.growth_retention == pytest.approx(expected.growth_retention, abs=1e-6)


def test_flat_induction_does_not_manufacture_an_identifiable_axis(readings):
    flat = readings.copy()
    for (plate, construct), group in flat.groupby(["plate", "construct"]):
        if construct in CONSTRUCTS[:2]:
            baseline = group.loc[group.dose_mM == 0, "activity_late"].mean()
            flat.loc[group.index, "activity_late"] = baseline
    with pytest.raises(ValueError, match="UPR|identifiable|induction"):
        fit_sensor_model(flat, training_plates=TRAINING_PLATES)


def test_validation_never_refits_or_uses_held_out_controls_for_predictions(
    readings, model, monkeypatch
):
    def forbidden(*args, **kwargs):
        raise AssertionError("validation must not refit")

    monkeypatch.setattr(hybrid_stress, "fit_sensor_model", forbidden)
    monkeypatch.setattr(hybrid_stress, "fit_latent", forbidden)
    monkeypatch.setattr(hybrid_stress, "least_squares", forbidden)
    held = readings[readings.plate == "held-out"].copy()
    before = json.dumps(model.metadata, sort_keys=True)
    table = model.validation_table(held)
    assert len(table) == len(held.groupby(["plate", "construct", "stressor", "dose_mM"]))
    required = {
        "plate", "construct", "stressor", "dose_mM", "held_out", "extrapolated",
        "activity_late_measured", "activity_late_predicted", "activity_fold_measured",
        "activity_fold_predicted", "activity_fold_error", "mu_late_measured",
        "mu_late_predicted", "growth_retention_measured", "growth_retention_predicted",
        "growth_retention_error", "activity_valid", "growth_valid",
    }
    assert required <= set(table.columns)
    assert table.held_out.all()
    supported = table[~table.extrapolated]
    assert np.max(np.abs(supported.activity_fold_error)) < 0.02
    assert np.max(np.abs(supported.growth_retention_error)) < 0.01
    assert table.activity_fold_error.to_numpy() == pytest.approx(
        (table.activity_fold_predicted - table.activity_fold_measured).to_numpy()
    )
    held.loc[held.dose_mM == 0, ["activity_late", "mu_late", "naive_late"]] *= 4
    changed = model.validation_table(held)
    prediction_columns = [name for name in table if name.endswith("_predicted")]
    pd.testing.assert_frame_equal(table[prediction_columns], changed[prediction_columns])
    assert json.dumps(model.metadata, sort_keys=True) == before


@pytest.mark.integration
def test_real_sensor_input_reports_quality_without_claiming_latent_ground_truth():
    path = paths.outputs_dir() / "sensor_characterisation.csv"
    if not path.exists():
        pytest.skip("committed sensor characterisation table is absent")
    frame = pd.read_csv(path)
    plates = tuple(sorted(frame.plate.unique()))
    fitted = fit_sensor_model(frame, training_plates=plates[:-1], seed=0)
    table = fitted.validation_table(frame[frame.plate == plates[-1]])
    metrics = []
    for construct, group in table.groupby("construct"):
        supported = group[(~group.extrapolated) & group.activity_valid & group.growth_valid]
        metrics.append({
            "construct": construct,
            "supported_conditions": len(supported),
            "activity_fold_rmse": float(np.sqrt(np.mean(supported.activity_fold_error**2))),
            "growth_retention_rmse": float(np.sqrt(np.mean(supported.growth_retention_error**2))),
        })
    all_plates = fit_sensor_model(frame, seed=0)
    prospective = {
        stressor: all_plates.predict_condition(stressor, 0.5).__dict__
        for stressor in ("DTT", "H2O2")
    }
    print(json.dumps({
        "held_out_plate": plates[-1],
        "training_plates": fitted.metadata["training_plates"],
        "calibrated_ranges": fitted.metadata["calibrated_ranges"],
        "excluded_training_conditions": len(fitted.metadata["exclusions"]),
        "held_out_sensor_quality": metrics,
        "synthetic_module_rmse": fitted.metadata["synthetic_module_rmse"],
        "evidence_status": fitted.metadata["evidence_status"],
        "all_plate_prospective_simulations": prospective,
    }, indent=2, allow_nan=False))
    json.dumps(all_plates.metadata, allow_nan=False)
    assert table.held_out.all()
    assert fitted.metadata["training_plates"] == list(plates[:-1])
    assert set(fitted.metadata["observation_schema"]) == set(FEATURES)
    assert all(0 <= value["growth_retention"] <= 1 for value in prospective.values())
