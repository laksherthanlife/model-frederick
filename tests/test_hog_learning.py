from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.hog_learning import (
    NativeHogFit, evaluate_hog_fit, fit_hog_parameters, fit_native_hog_parameters,
)
from ystwin.mech.hog import HogProtocol
from ystwin.mech.kinetic_sbml import KineticSimulationError


class ResponseModel:
    parameter_values = {"kv16f_1": 0.008, "kv16r_1": 0.003}
    source_metadata = {"publication": {"doi": "fixture"}}

    def simulate(self, times_s, *, protocol, parameters=None, **kwargs):
        parameters = {**self.parameter_values, **(parameters or {})}
        relative = np.maximum(np.asarray(times_s) - protocol.shock_time_s, 0)
        values = (1 - np.exp(-parameters["kv16f_1"] * relative)) * np.exp(-parameters["kv16r_1"] * relative)
        return SimpleNamespace(variables={"Hog1PP_measured": values})


def training_frame():
    times = np.linspace(0, 1800, 31)
    values = 3e-7 * (1 - np.exp(-0.014 * times)) * np.exp(-0.0018 * times)
    return pd.DataFrame({"time_model_s": times, "time_relative_s": times, "value": values,
                         "nacl_molar": 0.4, "observable_id": "Hog1PP_measured", "split": "fit",
                         "genotype": "wild_type", "experiment_id": "training", "unit": "mol/L"})


def test_rates_and_observation_gain_are_fitted_from_values_not_copied_from_source():
    model = ResponseModel()
    before = model.parameter_values.copy()
    fit = fit_hog_parameters(model, training_frame(), parameter_names=("kv16f_1", "kv16r_1"), max_nfev=40)
    assert fit.parameter_values["kv16f_1"] == pytest.approx(0.014, rel=0.002)
    assert fit.parameter_values["kv16r_1"] == pytest.approx(0.0018, rel=0.002)
    assert fit.observation_gain == pytest.approx(3e-7, rel=0.002)
    assert fit.training_nrmse < 1e-4
    assert fit.diagnostics["local_sensitivity_rank"] == 2
    assert model.parameter_values == before


def test_fitting_refuses_heldout_rows_and_unexcited_observations():
    frame = training_frame()
    frame.loc[0, "split"] = "heldout_dose"
    with pytest.raises(ValueError, match="fit"):
        fit_hog_parameters(ResponseModel(), frame)
    frame["split"] = "fit"
    frame["value"] = 0.0
    with pytest.raises(ValueError, match="excitation"):
        fit_hog_parameters(ResponseModel(), frame)


def test_heldout_values_change_scores_not_fitted_parameters():
    model = ResponseModel()
    fit = fit_hog_parameters(model, training_frame(), parameter_names=("kv16f_1", "kv16r_1"), max_nfev=40)
    before = fit.parameter_values.copy()
    heldout = training_frame()
    heldout["value"] /= heldout.value.max()
    heldout["split"] = "heldout_dose"
    heldout["experiment_id"] = "heldout"
    heldout["measurement_block"] = "measurement1"
    heldout["replicate_label"] = "trace1"
    heldout["unit"] = "relative_intensity"
    first, _ = evaluate_hog_fit(model, fit, heldout)
    changed = heldout.copy()
    changed["value"] = changed.value.to_numpy()[::-1]
    second, _ = evaluate_hog_fit(model, fit, changed)
    assert first.iloc[0].fitted_shape_rmse < 0.001
    assert second.iloc[0].fitted_shape_rmse > 0.1
    assert fit.parameter_values == before


def test_parameter_protocol_and_mutant_rows_are_not_silently_accepted():
    frame = training_frame()
    frame.loc[0, "nacl_molar"] = 0.8
    with pytest.raises(ValueError, match="protocol"):
        fit_hog_parameters(ResponseModel(), frame)
    frame = training_frame()
    frame["genotype"] = "hog1_del"
    with pytest.raises(ValueError, match="wild_type"):
        fit_hog_parameters(ResponseModel(), frame)


class NativeResponseModel:
    parameter_values = {
        "kv6_1": 0.01, "kv6b_4": 0.00004, "kv13a_1": 0.002,
        "kv16f_1": 0.008, "kv16r_1": 0.003, "kv13b_1": 0.0001,
        "kv17f_1": 0.002, "kv18f_1": 0.01,
    }

    def __init__(self):
        self.calls = []

    def simulate(self, times_s, *, protocol, parameters=None, **kwargs):
        times = np.asarray(times_s, dtype=float)
        assert times[0] == 0 and np.all(np.diff(times) > 0)
        self.calls.append((protocol, dict(parameters or {}), times.copy()))
        current = {**self.parameter_values, **(parameters or {})}
        relative = np.maximum(times - protocol.shock_time_s, 0)
        dose = protocol.nacl_molar / 0.4
        synthesis = (0.6 * current["kv6_1"] / self.parameter_values["kv6_1"]
                     + 0.4 * current["kv6b_4"] / self.parameter_values["kv6b_4"])
        diffusion = current["kv13a_1"] / self.parameter_values["kv13a_1"]
        uptake = current["kv13b_1"] / self.parameter_values["kv13b_1"]
        hog = (0.03 + dose * (1 - np.exp(-current["kv16f_1"] * relative))
               * np.exp(-current["kv16r_1"] * relative))
        gpd = 0.08 + dose * 0.3 * (1 - np.exp(-relative / 700))
        inside = (0.05 + synthesis * (0.03 * (1 - np.exp(-times / 700))
                  + dose * 0.3 * (1 - np.exp(-relative / 400)))
                  / (1 + diffusion * (1 - np.exp(-relative / 900)))
                  + 1e-10 * uptake * (1 - np.exp(-relative / 300)))
        outside = (0.002 + synthesis * diffusion * (0.0005 * (1 - np.exp(-times / 1000))
                   + dose * 0.002 * (1 - np.exp(-relative / 700))))
        return SimpleNamespace(variables={
            "Hog1PP_measured": hog, "Gpd1_measured": gpd,
            "glycerol_measured": inside, "glycerol_e": outside,
            "glycerol_i": inside * 7, "cellvol": np.full(len(times), 0.2),
            "cellvol_init": np.ones(len(times)),
        })


def native_training_frame(*, reference=False):
    model = NativeResponseModel()
    parameters = {} if reference else {
        "kv6_1": 0.017, "kv6b_4": 0.000068, "kv13a_1": 0.0012, "kv16f_1": 0.012,
    }
    gains = {"Hog1PP_measured": 4e-7, "Gpd1_measured": 9e-7}
    times = np.array([0, 60, 120, 180, 300, 500, 800, 1200, 1800], dtype=float)
    records = []
    for dose in (0.4, 0.0):
        result = model.simulate(times, protocol=HogProtocol(dose, 120), parameters=parameters)
        for observable in ("Hog1PP_measured", "Gpd1_measured", "glycerol_measured", "glycerol_e"):
            protein = observable in gains
            if dose == 0 and protein:
                continue
            for index, time in enumerate(times):
                value = result.variables[observable][index] * gains.get(observable, 1)
                records.append({
                    "time_model_s": time, "time_relative_s": time - 120, "value": value,
                    "nacl_molar": dose, "observable_id": observable, "split": "fit",
                    "genotype": "wild_type", "medium": "YPD", "supplement_id": "s002",
                    "experiment_id": f"native:{dose}:{observable if protein else 'hplc'}",
                    "unit": "mol/L", "sd": value * 0.2 if protein else np.nan,
                    "measurement_provenance": ("literature_scaled_western_blot" if protein
                                               else "processed_hplc_measurement"),
                })
    return pd.DataFrame(records)


def test_native_fit_recovers_shared_kinetics_without_rescaling_hplc_or_volume_twice():
    model = NativeResponseModel()
    source = model.parameter_values.copy()
    frame = native_training_frame()
    original = frame.copy(deep=True)
    fit = fit_native_hog_parameters(model, frame, max_nfev=30)
    assert isinstance(fit, NativeHogFit)
    assert fit.parameter_values == pytest.approx({
        "kv6_1": 0.017, "kv6b_4": 0.000068, "kv13a_1": 0.0012, "kv16f_1": 0.012,
    }, rel=0.002)
    assert fit.observation_gains == pytest.approx({"Hog1PP_measured": 4e-7, "Gpd1_measured": 9e-7}, rel=2e-4)
    assert fit.training_nrmse < 1e-5
    assert fit.diagnostics["kinetic_parameter_count"] == 3
    assert fit.diagnostics["effective_parameter_count"] == 5
    assert fit.diagnostics["local_sensitivity_rank"] == 3
    assert fit.diagnostics["fixed_observation_gains"] == {"glycerol_measured": 1.0, "glycerol_e": 1.0}
    assert fit.diagnostics["fixed_deactivation_parameter"] == source["kv16r_1"]
    assert fit.diagnostics["log_parameter_difference_step"] == 1e-3
    assert len(fit.diagnostics["protocols"]) == 2
    assert all(set(parameters) <= set(fit.parameter_values) | {"kv13b_1"}
               for _, parameters, _ in model.calls)
    assert model.parameter_values == source
    pd.testing.assert_frame_equal(frame, original)
    json.dumps(fit.diagnostics, allow_nan=False)


def test_native_glycerol_disagreement_cannot_be_hidden_in_a_concentration_gain():
    frame = native_training_frame(reference=True)
    frame.loc[frame.observable_id == "glycerol_measured", "value"] *= 2
    fit = fit_native_hog_parameters(
        NativeResponseModel(), frame, parameter_groups={"activation_balance": ("kv16f_1",)},
        max_nfev=2)
    assert fit.diagnostics["per_observable_errors"]["glycerol_measured"]["nrmse"] == pytest.approx(0.5)
    assert fit.training_nrmse == pytest.approx(0.25)
    assert set(fit.observation_gains) == {"Hog1PP_measured", "Gpd1_measured"}


def test_native_loss_uses_training_scales_and_equal_observable_not_equal_row_weight():
    frame = native_training_frame()
    first = fit_native_hog_parameters(
        NativeResponseModel(), frame, observable_weights={"Hog1PP_measured": 2.0}, max_nfev=1)
    scales = frame.groupby("observable_id").value.apply(lambda values: np.sqrt(np.mean(values ** 2)))
    assert first.diagnostics["normalization_scales"] == pytest.approx(scales.to_dict())
    expected_weights = {name: 0.4 if name == "Hog1PP_measured" else 0.2 for name in scales.index}
    assert first.diagnostics["observable_weights"] == pytest.approx(expected_weights)
    errors = first.diagnostics["per_observable_errors"]
    expected_loss = sum(expected_weights[name] * errors[name]["nrmse"] ** 2 for name in scales.index)
    assert first.training_nrmse ** 2 == pytest.approx(expected_loss)
    assert first.diagnostics["sd_treatment"] == "reported SD retained for diagnostics, not inverse-variance weighted"
    changed = frame.copy()
    changed.loc[changed.sd.notna(), "sd"] *= 100
    second = fit_native_hog_parameters(
        NativeResponseModel(), changed, observable_weights={"Hog1PP_measured": 2.0}, max_nfev=1)
    assert second.parameter_values == first.parameter_values
    assert second.observation_gains == first.observation_gains
    assert second.training_nrmse == first.training_nrmse
    assert first.diagnostics["optimization_converged"] is False


@pytest.mark.parametrize("column,value,match", [
    ("split", "unused", "fit"),
    ("split", "heldout_assay", "fit"),
    ("supplement_id", "s004", "s002"),
    ("genotype", "hog1_del", "wild_type"),
    ("medium", "different_medium", "YPD"),
    ("observable_id", "glycerol_i", "observable"),
    ("unit", "relative_intensity", "mol/L"),
    ("measurement_provenance", "unscaled_measurement", "provenance"),
    ("time_model_s", -1.0, "time"),
    ("time_relative_s", np.nan, "time"),
    ("value", np.inf, "finite"),
    ("sd", -1.0, "SD"),
])
def test_native_fit_refuses_unassigned_non_native_or_ambiguous_training(column, value, match):
    frame = native_training_frame()
    frame.loc[0, column] = value
    model = NativeResponseModel()
    with pytest.raises(ValueError, match=match):
        fit_native_hog_parameters(model, frame)
    assert not model.calls


def test_native_fit_validates_protocol_per_experiment_and_does_not_mutate_indices():
    frame = native_training_frame()
    frame["experiment_id"] = "inconsistent"
    with pytest.raises(ValueError, match="protocol"):
        fit_native_hog_parameters(NativeResponseModel(), frame)
    frame = native_training_frame()
    frame.index = np.zeros(len(frame), dtype=int)
    fit = fit_native_hog_parameters(NativeResponseModel(), frame, max_nfev=1)
    assert fit.diagnostics["training_rows"] == len(frame)
    np.testing.assert_array_equal(frame.index, np.zeros(len(frame), dtype=int))


@pytest.mark.parametrize("options,match", [
    ({"parameter_groups": {}}, "groups"),
    ({"parameter_groups": {"empty": ()}}, "groups"),
    ({"parameter_groups": {"one": ("kv13a_1",), "two": ("kv13a_1",)}}, "distinct"),
    ({"parameter_groups": {"one": ("not_a_parameter",)}}, "positive"),
    ({"parameter_groups": {"one": ("kv6_1",), "two": ("kv6b_4",)}}, "shared"),
    ({"parameter_groups": {"deactivation": ("kv16r_1",)}}, "fixed"),
    ({"fold_bounds": (1, 2)}, "bounds"),
    ({"max_nfev": True}, "max_nfev"),
    ({"max_nfev": 0}, "max_nfev"),
    ({"observable_weights": {"glycerol_measured": -1}}, "weights"),
    ({"observable_weights": {"glycerol_measured": "large"}}, "weights"),
    ({"observable_weights": {"glycerol_measured": np.nan}}, "weights"),
    ({"observable_weights": {"unobserved": 1}}, "weights"),
])
def test_native_fit_rejects_invalid_or_confounded_parameterizations(options, match):
    with pytest.raises(ValueError, match=match):
        fit_native_hog_parameters(NativeResponseModel(), native_training_frame(), **options)


def test_native_fitter_records_failed_trials_but_never_uses_them_as_sensitivities():
    class OneSidedModel(NativeResponseModel):
        def simulate(self, times_s, *, protocol, parameters=None, **kwargs):
            if (parameters or {}).get("kv6_1", self.parameter_values["kv6_1"]) > self.parameter_values["kv6_1"]:
                raise KineticSimulationError("fixture invalid positive synthesis perturbation")
            return super().simulate(times_s, protocol=protocol, parameters=parameters, **kwargs)

    fit = fit_native_hog_parameters(
        OneSidedModel(), native_training_frame(reference=True),
        parameter_groups={"glycerol_synthesis": ("kv6_1", "kv6b_4")}, max_nfev=1)
    assert fit.training_nrmse < 1e-12
    assert fit.diagnostics["failed_simulations"] > 0
    assert fit.diagnostics["simulation_failures"]
    assert fit.diagnostics["local_sensitivity_steps"]["glycerol_synthesis"] == -1e-3
    assert fit.diagnostics["local_sensitivity_rank"] == 1


def test_native_finite_differences_use_full_log_steps_near_bounds(monkeypatch):
    import ystwin.analysis.hog_learning as learning

    def near_bound_result(*args, **kwargs):
        return SimpleNamespace(x=np.array([np.log(1.05) - 1e-7]), success=False,
                               nfev=1, message="fixture bounded iterate", active_mask=np.array([0]))

    monkeypatch.setattr(learning, "least_squares", near_bound_result)
    fit = fit_native_hog_parameters(
        NativeResponseModel(), native_training_frame(), fold_bounds=(0.95, 1.05),
        parameter_groups={"activation_balance": ("kv16f_1",)}, max_nfev=1)
    assert fit.diagnostics["local_sensitivity_steps"]["activation_balance"] == -1e-3
    assert fit.diagnostics["bound_proximity_fraction"]["activation_balance"] < 0.05


def test_native_fitter_rejects_a_failed_reference_or_nonfinite_observable():
    class FailedModel(NativeResponseModel):
        def simulate(self, *args, **kwargs):
            raise KineticSimulationError("fixture failed integration")

    class NonfiniteModel(NativeResponseModel):
        def simulate(self, *args, **kwargs):
            result = super().simulate(*args, **kwargs)
            result.variables["glycerol_measured"][0] = np.nan
            return result

    for model in (FailedModel(), NonfiniteModel()):
        with pytest.raises(RuntimeError, match="reference"):
            fit_native_hog_parameters(model, native_training_frame())


def test_native_uptake_stays_a_prior_with_an_explicit_weak_local_sensitivity_diagnostic():
    fit = fit_native_hog_parameters(NativeResponseModel(), native_training_frame(), max_nfev=1)
    assert "kv13b_1" not in fit.parameter_values
    uptake = fit.diagnostics["fixed_parameter_sensitivities"]["kv13b_1"]
    assert uptake["column_norm"] < 1e-7
    assert uptake["orthogonal_column_norm"] < 1e-7
    assert uptake["estimated"] is False
    groups = {**fit.diagnostics["parameter_groups"], "active_uptake": ("kv13b_1",)}
    expanded = fit_native_hog_parameters(
        NativeResponseModel(), native_training_frame(), parameter_groups=groups, max_nfev=1)
    assert expanded.diagnostics["kinetic_parameter_count"] == 4
    assert expanded.diagnostics["effective_parameter_count"] == 6
    assert expanded.diagnostics["local_sensitivity_rank"] < 4
    assert "active_uptake" in expanded.diagnostics["weak_parameter_groups"]
    json.dumps(expanded.diagnostics, allow_nan=False)
