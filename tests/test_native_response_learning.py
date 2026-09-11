from __future__ import annotations

import copy
import hashlib
import json
import runpy
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pytest

from ystwin.analysis.native_response_learning import (
    ContractError,
    NativeExperiment,
    NativeResponseModel,
    PhysicalHistory,
    Readout,
    SelectionRecipe,
    aggregate_observations,
    causal_response_basis,
    check_data_gate,
    fit_native_response,
    load_verified_train_development,
    verify_freeze,
    write_freeze,
)


READOUT = Readout(
    observable_id="synthetic_localization_assay",
    context_id="synthetic_single_population_and_assay_context",
    unit="dimensionless",
    quantity_kind="nuclear_localization_proxy",
    normalization="raw_max5_over_median",
    normalization_reference="synthetic engineering fixture, not biological evidence",
)
RECIPE = SelectionRecipe(
    tau_seconds=(2.0, 7.0),
    ridge=(0.0, 0.01),
    max_state_terms=2,
    time_only_degrees=(1, 2),
    support_policy="flag",
)


def experiment(name, dose=0.4, onset=0.0, partition="train", gain=1.5, adaptive=False):
    times = tuple(float(t) for t in range(-4, 25))
    elapsed = np.maximum(np.asarray(times) - onset, 0.0)
    response = gain * dose * (1.0 - np.exp(-elapsed / 2.0))
    if adaptive:
        response -= 0.8 * gain * dose * (1.0 - np.exp(-elapsed / 7.0))
    return NativeExperiment(
        experiment_id=name,
        biological_replicate_id=f"replicate_{name}",
        source_groups=(f"movie_{name}",),
        partition=partition,
        history=PhysicalHistory(
            quantity="solute_concentration",
            dose_unit="M",
            time_unit="s",
            start_time=-4.0,
            initial_dose=0.0,
            changes=((onset, dose),),
        ),
        times=times,
        values=(tuple(1.25 + response),),
        readout=READOUT,
    )


def panel(adaptive=False):
    train = [
        experiment("low", 0.15, 1.0, adaptive=adaptive),
        experiment("high", 0.6, 0.0, adaptive=adaptive),
        experiment("middle", 0.3, 4.0, adaptive=adaptive),
    ]
    development = [experiment("development", 0.4, 2.0, "development", adaptive=adaptive)]
    return train, development


def learned(adaptive=False):
    return fit_native_response(*panel(adaptive), recipe=RECIPE)


def model(result):
    return NativeResponseModel.from_dict(result["model"])


def test_synthetic_sustained_and_adaptive_responses_do_not_require_fixed_adaptation():
    sustained = learned()
    adaptive = learned(adaptive=True)
    for result in (sustained, adaptive):
        assert result["selected_family"] == "causal_state"
        assert result["selection"]["development_mse"] < 1e-18
        assert result["claim"] == "empirical_response_only_not_a_universal_mechanism"
        assert result["data_status"] == "unverified_or_synthetic"
        assert result["selection_recipe"]["refit_on_development"] is False
    assert len(sustained["model"]["candidate"]["terms"]) == 1
    assert len(adaptive["model"]["candidate"]["terms"]) == 2
    coefficients = adaptive["model"]["coefficients"]
    assert min(coefficients) < 0.0 < max(coefficients)


def test_time_origin_and_unit_equivariance():
    train, development = panel()
    original = model(fit_native_response(train, development, recipe=RECIPE))
    query = development[0]
    expected = original.predict(query.history, query.times)
    shifted = replace(
        query.history,
        start_time=query.history.start_time + 50000.0,
        changes=tuple((t + 50000.0, u) for t, u in query.history.changes),
    )
    np.testing.assert_allclose(
        original.predict(shifted, tuple(t + 50000.0 for t in query.times))["values"],
        expected["values"],
        rtol=0.0,
        atol=1e-12,
    )
    converted = replace(
        query.history,
        time_unit="min",
        dose_unit="mM",
        start_time=query.history.start_time / 60.0,
        initial_dose=query.history.initial_dose * 1000.0,
        changes=tuple((t / 60.0, u * 1000.0) for t, u in query.history.changes),
    )
    actual = original.predict(converted, tuple(t / 60.0 for t in query.times))
    np.testing.assert_allclose(actual["values"], expected["values"], atol=1e-12)
    assert actual["support_flags"] == expected["support_flags"]


def test_labels_are_not_doses_and_physical_dose_changes_prediction():
    train, development = panel()
    renamed = [
        replace(e, experiment_id=f"dose_999_{i}", biological_replicate_id=f"other_{i}", source_groups=(f"other_movie_{i}",))
        for i, e in enumerate(train)
    ]
    a = model(fit_native_response(train, development, recipe=RECIPE))
    b = model(fit_native_response(renamed, development, recipe=RECIPE))
    query = development[0]
    np.testing.assert_allclose(a.predict(query.history, query.times)["values"], b.predict(query.history, query.times)["values"])
    doubled = replace(query.history, changes=((2.0, 0.8),))
    assert a.predict(doubled, query.times)["values"][-1] > a.predict(query.history, query.times)["values"][-1]


def test_causal_features_do_not_use_future_events_or_observations():
    e = experiment("causality")
    times = (-4.0, -1.0, 0.0, 1.0, 2.0)
    base = causal_response_basis(e.history, times, (2.0, 7.0))
    future = replace(e.history, changes=e.history.changes + ((100.0, 999.0),))
    repeated = replace(e.history, changes=((0.0, 0.4), (0.5, 0.4), (100.0, 999.0)))
    np.testing.assert_array_equal(base, causal_response_basis(future, times, (2.0, 7.0)))
    np.testing.assert_allclose(base, causal_response_basis(repeated, times, (2.0, 7.0)), atol=1e-15)
    np.testing.assert_array_equal(base[:2], 0.0)
    assert base[2, 0] == 0.4
    np.testing.assert_array_equal(base[2, 1:], 0.0)
    fitted = model(learned())
    assert fitted.predict(future, times) == fitted.predict(e.history, times)


def test_masks_preserve_missingness_and_ignore_masked_payloads():
    e = replace(
        experiment("masked"),
        times=(-4.0, 0.0, 1.0),
        values=((1.0, None, 3.0), (3.0, 1e100, None)),
        mask=((True, False, True), (True, False, False)),
    )
    summary = aggregate_observations(e)
    assert summary["mean"] == [2.0, None, 3.0]
    assert summary["count"] == [2, 0, 1]
    assert summary["mask"] == [True, False, True]
    with pytest.raises(ContractError, match="observed.*finite"):
        replace(e, values=((1.0, None, float("nan")), e.values[1]))
    with pytest.raises(ContractError, match="mask"):
        replace(e, mask=((True,),))


def test_calibration_is_train_only_and_predictions_never_use_query_readouts():
    train, development = panel()
    first = fit_native_response(train, development, recipe=RECIPE)
    changed = [replace(development[0], values=(tuple(1000.0 for _ in development[0].times),))]
    second = fit_native_response(train, changed, recipe=RECIPE)
    assert first["observation_calibration"] == second["observation_calibration"]
    assert first["baselines"]["constant"]["model"] == second["baselines"]["constant"]["model"]
    assert first["observation_calibration"]["fitted_partition"] == "train"
    assert first["observation_calibration"]["unit"] == "dimensionless"
    assert model(first).predict(development[0].history, development[0].times)["readout"] == asdict(READOUT)


@pytest.mark.parametrize("field", ["source_groups", "biological_replicate_id", "experiment_id"])
def test_whole_source_and_replicate_split_leakage_is_rejected(field):
    train, development = panel()
    development[0] = replace(development[0], **{field: getattr(train[0], field)})
    with pytest.raises(ContractError, match="leakage|duplicate"):
        fit_native_response(train, development, recipe=RECIPE)


def test_source_groups_require_a_structured_sequence_not_a_label_string():
    with pytest.raises(ContractError, match="source groups"):
        replace(experiment("source_contract"), source_groups="abc")


def test_masked_future_inputs_do_not_create_identifying_training_variation():
    train = [experiment(f"repeat_{i}") for i in range(3)]
    extended = []
    for i, e in enumerate(train):
        extended.append(replace(e, times=e.times + (200.0,), values=(e.values[0] + (None,),), history=replace(e.history, changes=e.history.changes + ((100.0, float(i + 2)),))))
    development = [experiment("dev", partition="development")]
    result = fit_native_response(extended, development, recipe=RECIPE)
    assert "single_physical_history_input_vs_time_not_identifiable" in result["failures"]
    assert result["development_input_response_gain"] is False


def test_duplicate_source_alias_within_training_and_reserved_partition_are_rejected():
    train, development = panel()
    alias = replace(train[0], experiment_id="alias", biological_replicate_id="alias_rep")
    with pytest.raises(ContractError, match="duplicate|source"):
        fit_native_response(train + [alias], development, recipe=RECIPE)
    with pytest.raises(ContractError, match="partition"):
        replace(development[0], partition="reserved_test")


def test_readout_quantity_units_and_normalization_cannot_be_silently_reinterpreted():
    for changed in ({"unit": "absolute_enzyme"}, {"quantity_kind": "phosphorylated_fraction"}, {"normalization": "unaudited_published"}):
        with pytest.raises(ContractError):
            replace(READOUT, **changed)
    fitted = model(learned())
    e = experiment("query")
    with pytest.raises(ContractError, match="readout"):
        fitted.predict(e.history, e.times, readout=replace(READOUT, observable_id="different_assay"))
    with pytest.raises(ContractError, match="input.*unit|input.*contract"):
        fitted.predict(replace(e.history, dose_unit="g/L"), e.times)
    with pytest.raises(ContractError, match="unit"):
        replace(e.history, dose_unit="percent")
    with pytest.raises(ContractError, match="prehistory"):
        replace(e.history, prehistory="unknown")
    with pytest.raises(ContractError, match="readout"):
        fitted.predict(e.history, e.times, readout=replace(READOUT, context_id="different_background"))


def test_serialized_causal_family_cannot_hide_a_time_only_term():
    raw = learned()["model"]
    raw["candidate"]["terms"] = [{"kind": "elapsed", "degree": 1}]
    raw["coefficients"] = [1.0]
    with pytest.raises(ContractError, match="candidate"):
        NativeResponseModel.from_dict(raw)


def test_out_of_support_is_flagged_or_refused_not_clipped():
    fitted = model(learned())
    e = experiment("query", dose=0.6)
    large = replace(e.history, changes=((0.0, 6.0),))
    inside = fitted.predict(e.history, e.times)
    outside = fitted.predict(large, e.times)
    assert outside["values"][-1] > inside["values"][-1] + 1.0
    assert "input_level" in outside["support_flags"][-1]
    refused = fitted.predict(large, e.times, support_policy="refuse")
    assert refused["values"][-1] is None
    assert refused["status"] == "partially_refused_out_of_support"
    late = fitted.predict(e.history, (40.0,))
    assert "elapsed_seconds" in late["support_flags"][0]
    switched = replace(e.history, changes=((0.0, 0.4), (2.0, 0.0), (4.0, 0.4)))
    assert "change_count" in fitted.predict(switched, (6.0,))["support_flags"][0]


def test_refusal_does_not_improve_scores_by_dropping_hard_development_points():
    train, development = panel()
    development[0] = replace(development[0], history=replace(development[0].history, changes=((2.0, 10.0),)))
    with pytest.raises(ContractError, match="complete development coverage") as failure:
        fit_native_response(train, development, recipe=replace(RECIPE, support_policy="refuse"))
    assert failure.value.attempts
    assert all("failure" in attempt for attempt in failure.value.attempts)


def test_no_input_and_constant_observation_null_is_not_a_learned_dynamic_law():
    train, development = panel()
    def constant(e):
        return replace(e, history=replace(e.history, changes=()), values=(tuple(2.0 for _ in e.times),))
    result = fit_native_response([constant(e) for e in train], [constant(e) for e in development], recipe=RECIPE)
    assert result["selected_family"] == "constant"
    assert "no_temporal_input_excitation" in result["failures"]
    assert "no_training_observation_variation" in result["failures"]
    assert result["empirical_input_response_supported"] is False
    assert any(row.get("failure") == "rank_deficient" for row in result["candidates"])


def test_time_only_comparator_is_explicitly_nonmechanistic():
    train, development = panel()
    def clock(e):
        return replace(e, values=(tuple(1.0 + 0.02 * (t - e.history.start_time) for t in e.times),))
    result = fit_native_response([clock(e) for e in train], [clock(e) for e in development], recipe=RECIPE)
    assert result["selected_family"] == "time_only"
    assert result["model"]["candidate"]["interpretation"] == "nonmechanistic_time_only_comparator"
    assert result["empirical_input_response_supported"] is False


def test_single_physical_history_does_not_identify_input_vs_clock():
    train = [experiment(f"repeat_{i}") for i in range(3)]
    development = [experiment("dev", partition="development")]
    result = fit_native_response(train, development, recipe=RECIPE)
    assert "single_physical_history_input_vs_time_not_identifiable" in result["failures"]
    assert result["empirical_input_response_supported"] is False


def test_development_scoring_balances_experiments_not_cells_or_timepoint_counts():
    train, development = panel()
    second = experiment("another_dev", 0.25, 1.0, "development")
    second = replace(second, values=tuple(tuple(v + 0.1 for v in second.values[0]) for _ in range(20)))
    result = fit_native_response(train, development + [second], recipe=RECIPE)
    rows = [row for row in result["predictions"] if row["partition"] == "development"]
    assert result["selection"]["development_mse"] == pytest.approx(np.mean([row["mse"] for row in rows]))
    assert result["selection_recipe"]["score_weighting"] == "equal_experiment_then_equal_observed_cell_time"
    assert result["selection"]["development_mae"] == pytest.approx(np.mean([row["mae"] for row in rows]))
    assert result["selection_recipe"]["selection_loss"] == "equal_experiment_mae"


def test_single_cell_noise_cannot_cancel_out_in_reported_prediction_loss():
    train, development = panel()
    e = development[0]
    development[0] = replace(e, values=(tuple(v - 0.2 for v in e.values[0]), tuple(v + 0.2 for v in e.values[0])))
    result = fit_native_response(train, development, recipe=RECIPE)
    assert result["selection"]["development_mae"] == pytest.approx(0.2)
    assert result["selection"]["development_mse"] >= 0.04 - 1e-12
    assert result["predictions"][-1]["n_observed_values"] == 2 * len(e.times)


def test_scoring_scale_uses_training_group_means_not_development_or_cell_count():
    train, development = panel()
    result = fit_native_response(train, development, recipe=RECIPE)
    means = [np.mean(e.values[0]) for e in train]
    expected = float(np.quantile(means, 0.75) - np.quantile(means, 0.25))
    assert result["scoring_scale"]["training_group_mean_iqr"] == pytest.approx(expected)
    changed = [replace(development[0], values=(tuple(v + 100.0 for v in development[0].values[0]),))]
    assert result["scoring_scale"] == fit_native_response(train, changed, recipe=RECIPE)["scoring_scale"]


def test_negative_localization_predictions_are_refused_without_clipping():
    fitted = model(learned())
    raw = fitted.to_dict()
    raw["intercept"] = -100.0
    query = panel()[1][0]
    prediction = NativeResponseModel.from_dict(raw).predict(query.history, query.times)
    assert all(value is None for value in prediction["values"])
    assert all("negative_localization_proxy" in flags for flags in prediction["domain_flags"])
    assert prediction["status"] == "refused_physical_domain"


def test_serialization_and_freeze_are_deterministic_and_tamper_evident(tmp_path):
    result = learned()
    restored = NativeResponseModel.from_dict(json.loads(json.dumps(result["model"], allow_nan=False)))
    query = panel()[1][0]
    assert restored.predict(query.history, query.times) == model(result).predict(query.history, query.times)
    path = tmp_path / "freeze.json"
    envelope = write_freeze(result, path)
    assert verify_freeze(json.loads(path.read_text())) is True
    assert write_freeze(result, tmp_path / "second.json") == envelope
    with pytest.raises(FileExistsError):
        write_freeze(result, path)
    changed = copy.deepcopy(envelope)
    changed["result"]["predictions"][0]["values"][0] += 1.0
    with pytest.raises(ContractError, match="digest"):
        verify_freeze(changed)


def test_recipe_search_is_bounded_and_has_visible_priors():
    with pytest.raises(ContractError, match="bounded|candidate"):
        replace(RECIPE, tau_seconds=tuple(float(i) for i in range(1, 20)))
    with pytest.raises(ContractError, match="tau"):
        replace(RECIPE, tau_seconds=(-1.0,))
    result = learned()
    assert "dz_tau/dt = (u - z_tau) / tau" in result["equations"]["causal_state"]
    assert result["selection_recipe"]["coefficient_sign_constraint"] == "none"
    assert len(result["candidates"]) <= 128


def verified_bundle(tmp_path):
    train, development = panel()
    exports = {}
    for partition, rows in (("train", train), ("development", development)):
        payload = {"schema_version": 1, "partition": partition, "experiments": [asdict(e) for e in rows]}
        path = tmp_path / f"{partition}.json"
        path.write_text(json.dumps(payload, allow_nan=False))
        exports[partition] = {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "contains_only_partition": partition,
            "experiment_ids": [e.experiment_id for e in rows],
        }
    manifest = {
        "schema_version": 1,
        "status": "verified_train_development_exports",
        "observations_origin": "experimental_measurements",
        "source_assets": [{"asset_id": "synthetic_fixture_only", "sha256": "1" * 64}],
        "verification": {
            "input_history_and_units": True,
            "readout_definition_and_normalization": True,
            "whole_experiment_lineage": True,
            "reserved_arrays_excluded": True,
            "no_teacher_or_kinetic_targets": True,
            "normalization_fit_scope": "none",
            "normalization_uses_future_or_development_outcomes": False,
            "positive_raw_denominators_verified": True,
            "assay_context_verified": True,
            "evidence": "Synthetic gate fixture, no source verification is claimed.",
        },
        "readout": asdict(READOUT),
        "experiments": [
            {key: asdict(e)[key] for key in ("experiment_id", "biological_replicate_id", "source_groups", "partition")}
            for e in train + development
        ],
        "approved_exports": exports,
    }
    manifest_path = tmp_path / "granados_manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    protocol = {
        "schema_version": 1,
        "status": "frozen_before_learning",
        "purpose": "native_observed_input_response_train_development",
        "test_access": "forbidden",
        "split_unit": "whole_biological_experiment",
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "selection_recipe": asdict(RECIPE),
    }
    protocol_path = tmp_path / "protocol.json"
    protocol_path.write_text(json.dumps(protocol))
    return manifest_path, protocol_path, manifest, protocol


def rewrite_metadata(manifest_path, protocol_path, manifest, protocol):
    manifest_path.write_text(json.dumps(manifest))
    protocol["manifest_sha256"] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    protocol_path.write_text(json.dumps(protocol))


def test_gate_checks_metadata_without_reading_exports_and_never_opens_reserved(tmp_path, monkeypatch):
    manifest_path, protocol_path, manifest, protocol = verified_bundle(tmp_path)
    manifest["approved_exports"]["reserved_test"] = {"path": "reserved.json", "sha256": "2" * 64}
    manifest["experiments"].append({"experiment_id": "reserved", "biological_replicate_id": "reserved_rep", "source_groups": ["reserved_movie"], "partition": "reserved_test"})
    rewrite_metadata(manifest_path, protocol_path, manifest, protocol)
    opened = []
    original = Path.read_bytes

    def guarded(path):
        opened.append(path.name)
        assert path.name != "reserved.json"
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    assert check_data_gate(manifest_path, protocol_path)["ready"] is True
    assert "train.json" not in opened and "development.json" not in opened
    loaded = load_verified_train_development(manifest_path, protocol_path)
    assert {e.partition for e in loaded["train"]} == {"train"}
    assert {e.partition for e in loaded["development"]} == {"development"}
    assert loaded["recipe"] == RECIPE
    assert "reserved.json" not in opened


@pytest.mark.parametrize("ineligible", [False, True])
def test_custodian_metadata_remains_blocked_until_a_verified_projection_adapter_exists(tmp_path, monkeypatch, ineligible):
    manifest_path = tmp_path / "granados_manifest.json"
    protocol_path = tmp_path / "protocol.json"
    manifest_path.write_text(json.dumps({"$schema": "granados/manifest.schema.json", "status": "metadata_schema_verified_ineligible_current_protocol" if ineligible else "metadata_schema_verified_eligibility_gated", "custody": {"export_authorized": ineligible, "training_allowlist": ["never_open.json"] if ineligible else []}}))
    protocol_path.write_text(json.dumps({"schema_version": 1, "status": "criteria_declared_before_new_fitting_or_measurement_access_not_a_scientific_verdict"}))
    original = Path.read_bytes

    def only_metadata(path):
        assert path.name in {"granados_manifest.json", "protocol.json"}
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", only_metadata)
    status = check_data_gate(manifest_path, protocol_path)
    assert status["ready"] is False
    assert "custodian" in status["failure"]
    if ineligible:
        assert "ineligible" in status["failure"]
    assert status["measurement_arrays_opened"] is False


def test_missing_or_unverified_protocol_blocks_before_any_measurement_read(tmp_path, monkeypatch):
    manifest_path, protocol_path, manifest, protocol = verified_bundle(tmp_path)
    protocol["status"] = "draft"
    protocol_path.write_text(json.dumps(protocol))
    original = Path.read_bytes
    opened = []

    def guarded(path):
        opened.append(path.name)
        assert path.name not in ("train.json", "development.json")
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    with pytest.raises(ContractError, match="protocol"):
        load_verified_train_development(manifest_path, protocol_path)
    assert check_data_gate(manifest_path, tmp_path / "missing_protocol.json")["ready"] is False


@pytest.mark.parametrize("defect", ["digest", "lineage", "path", "normalization", "future_normalization", "denominator", "reserved_alias"])
def test_verified_loader_rejects_tampering_or_unsafe_metadata_before_fit(tmp_path, defect):
    manifest_path, protocol_path, manifest, protocol = verified_bundle(tmp_path)
    if defect == "digest":
        manifest["approved_exports"]["train"]["sha256"] = "0" * 64
    elif defect == "lineage":
        manifest["experiments"][-1]["source_groups"] = manifest["experiments"][0]["source_groups"]
    elif defect == "path":
        manifest["approved_exports"]["train"]["path"] = "../outside.json"
    elif defect == "normalization":
        manifest["verification"]["normalization_fit_scope"] = "all_experiments"
    elif defect == "future_normalization":
        manifest["verification"]["normalization_uses_future_or_development_outcomes"] = True
    elif defect == "denominator":
        manifest["verification"]["positive_raw_denominators_verified"] = False
    else:
        manifest["experiments"].append({"experiment_id": "hidden", "biological_replicate_id": "hidden_rep", "source_groups": manifest["experiments"][0]["source_groups"], "partition": "reserved_test"})
    rewrite_metadata(manifest_path, protocol_path, manifest, protocol)
    with pytest.raises(ContractError):
        load_verified_train_development(manifest_path, protocol_path)


def runner_main():
    path = Path(__file__).resolve().parents[1] / "scripts/run_native_response_learning.py"
    return runpy.run_path(str(path))["main"]


def test_runner_contract_is_engineering_only_and_missing_gate_does_not_fit(tmp_path, capsys):
    main = runner_main()
    assert main(["--contract"]) == 0
    contract = json.loads(capsys.readouterr().out)
    assert contract["default_recipe_status"] == "engineering_example_not_registered_for_native_data"
    assert len(contract["candidate_library"]) <= 128
    assert main(["--manifest", str(tmp_path / "absent.json"), "--protocol", str(tmp_path / "missing.json")]) == 2
    gate = json.loads(capsys.readouterr().out)
    assert gate["measurement_arrays_opened"] is False
    assert gate["scientific_authorization"] == "none"


def test_runner_freezes_only_synthetic_fixture_and_never_overwrites(tmp_path, capsys):
    manifest_path, protocol_path, _, _ = verified_bundle(tmp_path)
    output = tmp_path / "synthetic_freeze.json"
    main = runner_main()
    argv = ["--manifest", str(manifest_path), "--protocol", str(protocol_path), "--fit", "--freeze", str(output)]
    assert main(argv) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "training_development_frozen_not_scientifically_authorized"
    assert summary["final_test_accessed"] is False
    envelope = json.loads(output.read_text())
    assert verify_freeze(envelope)
    assert envelope["result"]["empirical_input_response_supported"] is False
    assert envelope["result"]["provenance"]["learner_source_sha256"]
    assert main(argv) == 2
    refusal = json.loads(capsys.readouterr().out)
    assert refusal["measurement_access_attempted"] is False


def test_training_equivariance_uses_canonical_dose_and_relative_time():
    train, development = panel()

    def converted(e):
        history = replace(e.history, time_unit="min", dose_unit="mM", start_time=e.history.start_time / 60.0 + 100.0, initial_dose=e.history.initial_dose * 1000.0, changes=tuple((t / 60.0 + 100.0, dose * 1000.0) for t, dose in e.history.changes))
        return replace(e, history=history, times=tuple(t / 60.0 + 100.0 for t in e.times))

    original = learned()
    transformed = fit_native_response([converted(e) for e in train], [converted(e) for e in development], recipe=RECIPE)
    assert transformed["selected_family"] == original["selected_family"]
    np.testing.assert_allclose(transformed["model"]["coefficients"], original["model"]["coefficients"], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(transformed["predictions"][-1]["values"], original["predictions"][-1]["values"], atol=1e-10)
