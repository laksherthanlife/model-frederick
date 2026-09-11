"""Independent tests of the retrospective population-development contract.

All numbers constructed here are synthetic engineering fixtures. No reserved
raw source or unreleased development response is read. Tests exercise actual
population aggregation, fitting, prediction and custody interfaces.
"""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
ASSET = "1b38366b-7473-4ddf-9cf9-26e625533686"
FAMILIES = (
    "prefix_persistence", "training_constant", "prefix_offset",
    "prefix_affine_time", "first_order_response", "transient_response",
)


@pytest.fixture(scope="module")
def population():
    from ystwin.analysis import native_population_development

    return native_population_development


@pytest.fixture(scope="module")
def protocol(population):
    return population.load_protocol(ROOT / "data/native_law_v2/development_protocol.json")


def _group(rows=5):
    return {"group_id": f"{ASSET}:/rep1", "source_pointer": "/rep1", "source_rows": rows, "split": "train"}


def _cell(group, index, frame, ratio, median, time):
    pointer = group["source_pointer"]
    return {
        "group_id": group["group_id"], "split": group["split"],
        "window": "prefix" if 28 <= frame <= 47 else "response",
        "cell_index": index, "frame_1based": frame, "relative_time_min": time,
        "max5": ratio * median, "median": median,
        "time_pointer": f"{pointer}/general/times/{index}/{frame - 1}",
        "max5_pointer": f"{pointer}/GFP/max5/{index}/{frame - 1}",
        "median_pointer": f"{pointer}/GFP/median/{index}/{frame - 1}",
    }


def _source_cells():
    group = _group()
    cells = []
    for index in range(group["source_rows"]):
        for frame in range(28, 48):
            cells.append(_cell(group, index, frame, 1.0 + 2.0 * index + 0.1 * (frame - 28), 2.0 ** index, float(frame - 48)))
        for frame in range(50, 74):
            time = 1.0 + 10.0 * (frame - 50) + (0.0, 2.0, 8.0, 26.0, 1000.0)[index]
            cells.append(_cell(group, index, frame, 2.0 + 2.0 * index, 2.0 ** index, time))
    next(row for row in cells if row["cell_index"] == 4 and row["frame_1based"] == 28)["median"] = None
    return group, cells


def _prepared(population, protocol, levels=(1.0, 2.0, 7.0), train_mu=1.0, development_mu=(50.0, 60.0)):
    """Declared aggregate values, not a claim that these are biological data."""
    prefixes, frames, coordinates = [], [], []
    for index, spec in enumerate(population.group_specs(protocol)):
        count = spec["source_rows"]
        prefixes.append({
            "group_id": spec["group_id"], "origin_1based": 49,
            "eligible_cell_indices": list(range(count)),
            "mu": train_mu if index < 3 else development_mu[index - 3], "s": 0.0,
            "input_projection_sha256": "a" * 64, "operator_id": population.PREFIX_OPERATOR,
        })
        for frame in range(50, 74):
            row = {"group_id": spec["group_id"], "frame_1based": frame, "time_min": 3.0 * (frame - 49) + 0.1 * index}
            if index < 3:
                row.update(observed_mean=levels[index], prefix_cells=count, valid_cells=count,
                           source_missing_cells=0, invalid_cells=0, cell_sample_sd=0.0)
                frames.append(row)
            else:
                coordinates.append(row)
    return population.finalize_prepared(prefixes, frames, coordinates, protocol)


def test_population_estimand_keeps_prefix_cohort_and_invalid_response_times(population):
    group, cells = _source_cells()
    prefix, _ = population.derive_prefix(group, cells, "a" * 64)
    assert prefix["eligible_cell_indices"] == [0, 1, 2, 3]
    assert prefix["mu"] == pytest.approx(4.95)
    assert prefix["s"] == pytest.approx(np.sqrt(20.0 / 3.0))
    # Every future entry in two prefix-eligible rows becomes unusable. Neither
    # future completeness nor valid-response counts may redefine the cohort.
    for row in cells:
        if row["window"] == "response":
            if row["cell_index"] == 2 or (row["frame_1based"] == 51 and row["cell_index"] in (0, 1)):
                row["max5"] = None
            if row["cell_index"] == 3:
                row["median"] = 0.0
    after, _ = population.derive_prefix(group, cells, "a" * 64)
    assert after == prefix
    response = [row for row in cells if row["window"] == "response"]
    frames = population.derive_frames(group, response, response, prefix)
    assert [row["frame_1based"] for row in frames] == list(range(50, 74))
    first = frames[0]
    assert first["prefix_cells"] == 4
    assert (first["valid_cells"], first["source_missing_cells"], first["invalid_cells"]) == (2, 1, 1)
    assert first["time_min"] == pytest.approx(10.0)  # includes times for both invalid response cells
    assert first["observed_mean"] == pytest.approx(3.0)  # mean(2/1, 8/2), not 10/3
    assert first["cell_sample_sd"] == pytest.approx(np.sqrt(2.0))
    assert frames[1]["observed_mean"] is None
    assert frames[1]["cell_sample_sd"] is None
    assert frames[1]["prefix_cells"] == 4
    assert (frames[1]["valid_cells"], frames[1]["source_missing_cells"], frames[1]["invalid_cells"]) == (0, 3, 1)
    assert frames[1]["time_min"] == pytest.approx(20.0)


@pytest.mark.parametrize("mutation", ["intensity_as_time", "max5_as_median", "first_stress_frame", "extra_future_frame"])
def test_source_field_ancestry_and_exact_frame_windows_are_not_column_masks(population, mutation):
    group, cells = _source_cells()
    population.validate_group_records(group, cells, ("prefix", "response"))
    altered = deepcopy(cells)
    if mutation == "intensity_as_time":
        altered[0]["time_pointer"] = altered[0]["max5_pointer"]
    elif mutation == "max5_as_median":
        altered[0]["median_pointer"] = altered[0]["max5_pointer"]
    else:
        # All pointer coordinates agree with the rogue frame: refusal must rest
        # on the registered window, not an incidental stale-pointer mismatch.
        index = 0 if mutation == "first_stress_frame" else 20
        frame = 49 if mutation == "first_stress_frame" else 74
        altered[index] = _cell(group, 0, frame, 2.0, 1.0, 1.0)
        altered[index]["window"] = "prefix" if mutation == "first_stress_frame" else "response"
    with pytest.raises(population.DevelopmentError):
        population.validate_group_records(group, altered, ("prefix", "response"))


def test_real_training_fit_weights_cultures_not_cells_and_ignores_development_conditioning(population, protocol):
    prepared = _prepared(population, protocol)
    assert prepared["H_train"] == 7.0
    assert prepared["S_train"] == 6.0  # IQR of 72 residual means, not three trajectory means
    fitted = population.fit_training(prepared, protocol)
    assert tuple(fitted) == FAMILIES
    constant = fitted["training_constant"]
    assert constant["status"] == "ok"
    assert constant["parameters"]["c"] == pytest.approx(10.0 / 3.0, rel=1e-7)
    assert constant["training_mean_mse"] == pytest.approx(62.0 / 9.0, rel=1e-7)
    changed = _prepared(population, protocol, development_mu=(1050.0, 1060.0))
    again = population.fit_training(changed, protocol)
    assert changed["H_train"] == prepared["H_train"] and changed["S_train"] == prepared["S_train"]
    for model_id in FAMILIES:
        assert again[model_id]["parameters"] == fitted[model_id]["parameters"]
        assert again[model_id]["selected_start_index"] == fitted[model_id]["selected_start_index"]
        assert [a["initial_parameters"] for a in again[model_id]["attempts"]] == [a["initial_parameters"] for a in fitted[model_id]["attempts"]]
    predictions = population.predict_development(prepared, fitted, protocol)[0]
    conditioned = population.predict_development(changed, again, protocol)[0]
    assert len(predictions) == len(conditioned) == 288
    for before, after in zip(predictions, conditioned, strict=True):
        if before["model_id"] == "prefix_persistence":
            assert after["mean"] == pytest.approx(before["mean"] + 1000.0)
        elif before["model_id"] == "training_constant":
            assert after["mean"] == before["mean"]


def test_six_actual_families_allow_simple_truth_and_apply_the_registered_link(population, protocol):
    prepared = _prepared(population, protocol, levels=(2.0, 2.0, 2.0), train_mu=2.0, development_mu=(2.0, 2.0))
    assert prepared["S_train"] == 0.0  # undefined scaled diagnostics must not create an epsilon or failure
    fitted = population.fit_training(prepared, protocol)
    predictions, audit = population.predict_development(prepared, fitted, protocol)
    assert tuple(fitted) == FAMILIES
    assert len(predictions) == len(audit) == 288
    scores = {}
    for model_id in FAMILIES:
        rows = [row for row in predictions if row["model_id"] == model_id]
        assert len(rows) == 48 and all(row["status"] == "ok" for row in rows)
        scores[model_id] = {"complete": True, "mean_mse": float(np.mean([(row["mean"] - 2.0) ** 2 for row in rows]))}
    assert population.select_model(scores) == "prefix_persistence"
    # A registered common link is lawful; a later model-specific rescue clip is not.
    raw, mapped, _ = population.evaluate_mean("prefix_affine_time", {"a": -3.0, "b": 4.0}, 1.0, np.array([2.5, 30.0, 60.0]))
    assert raw[0] < 0.0
    np.testing.assert_array_equal(mapped, [0.0, 0.0, 2.0])
    assert population.mean_mse(mapped, np.array([1.0, 1.0, 3.0])) == pytest.approx(1.0)
    scores["transient_response"] = {"complete": False, "mean_mse": None}
    assert population.select_model(scores) is None
    incomplete = dict(fitted)
    incomplete.pop("transient_response")
    with pytest.raises(population.DevelopmentError):
        population.predict_development(prepared, incomplete, protocol)
    aliases = dict(fitted)
    aliases["observed_prefix_constant"] = deepcopy(fitted["prefix_persistence"])
    with pytest.raises(population.DevelopmentError):
        population.predict_development(prepared, aliases, protocol)


@pytest.mark.parametrize("model_id,parameter,value", [
    ("first_order_response", "tau", 0.25),
    ("training_constant", "c", 77.0),
])
def test_registered_prediction_seam_refuses_out_of_domain_checkpoint_parameters(population, protocol, model_id, parameter, value):
    prepared = _prepared(population, protocol)
    fitted = population.fit_training(prepared, protocol)
    assert all(row["status"] == "ok" for row in population.predict_development(prepared, fitted, protocol)[0])
    changed = deepcopy(fitted)
    changed[model_id]["parameters"][parameter] = value
    try:
        predictions, _ = population.predict_development(prepared, changed, protocol)
    except population.DevelopmentError:
        return
    affected = [row for row in predictions if row["model_id"] == model_id]
    assert len(affected) == 48
    assert all(row["status"] == "failed" and row["mean"] is None for row in affected), (
        "The protocol-bound predictor emitted successful predictions for parameters "
        "outside the registered family domain. Validate checkpoints against the "
        "training-derived H_train and fixed time-constant bounds, not merely finite values."
    )


@pytest.fixture(scope="module")
def custody():
    spec = importlib.util.spec_from_file_location(
        "population_adversarial_custody", ROOT / "scripts/export_granados_training.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_custodian_and_population_operators_agree_and_source_failure_is_not_a_schema_error(population, custody):
    # Invoke only the actual in-memory aggregation function. Never invoke
    # dev_selected_source/dev_source_audit, which would open the reserved source.
    group, cells = _source_cells()
    group["source_rows"] = custody.DEV_ROWS["rep1"]
    for index in range(5, group["source_rows"]):
        for frame in (*range(28, 48), *range(50, 74)):
            row = _cell(group, index, frame, 0.0, 1.0, -1.0 if frame < 49 else 1.0)
            row.update(max5=None, median=None, relative_time_min=None)
            cells.append(row)

    class WindowRow(list):
        def __getitem__(self, index):
            assert isinstance(index, int) and index + 1 in (*range(28, 48), *range(50, 74)), (
                "Custodial aggregation read the first stress frame or another unregistered frame"
            )
            return super().__getitem__(index)

    matrices = {field: [WindowRow([None] * 97) for _ in range(group["source_rows"])]
                for field in ("max5", "median", "times")}
    for row in cells:
        for field, column in (("max5", "max5"), ("median", "median"), ("times", "relative_time_min")):
            matrices[field][row["cell_index"]][row["frame_1based"] - 1] = row[column]
    matrices["origin"] = 49
    prefix, _ = population.derive_prefix(group, cells, "a" * 64)
    response = [row for row in cells if row["window"] == "response"]
    frames = population.derive_frames(group, response, response, prefix)
    authoritative = custody.dev_aggregate_group("rep1", matrices)
    assert authoritative["eligible_cell_indices"] == prefix["eligible_cell_indices"]
    assert authoritative["mu"] == pytest.approx(prefix["mu"])
    assert authoritative["s"] == pytest.approx(prefix["s"])
    for original, replayed in zip(authoritative["frames"], frames, strict=True):
        for field in ("frame_1based", "time_min", "observed_mean", "prefix_cells", "valid_cells", "source_missing_cells", "invalid_cells", "cell_sample_sd"):
            assert original[field] == pytest.approx(replayed[field])
    # Correctly shaped source data with an unusable prefix is a real coverage
    # failure, unlike an incorrect projection schema or adapter field name.
    for row in cells:
        if row["window"] == "prefix" and row["cell_index"] != 0:
            row["max5"] = None
            matrices["max5"][row["cell_index"]][row["frame_1based"] - 1] = None
    population.validate_group_records(group, cells, ("prefix", "response"))
    with pytest.raises(population.DevelopmentError) as refusal:
        population.derive_prefix(group, cells, "a" * 64)
    assert refusal.value.code == "prefix_cohort_unusable"
    unavailable = custody.dev_aggregate_group("rep1", matrices)
    assert unavailable["eligible_cell_indices"] == [0]
    assert unavailable["mu"] is None and unavailable["s"] is None


def test_development_success_flags_cannot_authorize_the_original_strong_gate(tmp_path, population, protocol):
    from ystwin.analysis.native_law_validation import ValidationTrust, assess_native_law

    development_report = {
        **population.environment_record(),
        "document_type": "development_diagnostics_not_biological_validation",
        "protocol_id": protocol["protocol_id"],
        "status": "complete_development_diagnostics",
        "strong_claim_authorized": True,
        "ready": True,
        "passed": True,
    }
    result = assess_native_law(development_report, root=tmp_path, trust=ValidationTrust((), {}))
    assert not result.authorized
    assert any(gate.state == "fail" for gate in result.gates)
    assert not list(tmp_path.iterdir())


def _seal_engineering_fit(population, protocol, prepared, fitted, root, name):
    store = population.ArtifactStore(root, f"outputs/native_population_development/{name}")
    code_ref = store.write_bytes("executed_population_code.py", Path(population.__file__).read_bytes(), "execution_code")
    provenance = {
        "approval_sha256": "b" * 64,
        "training_projection_sha256": "a" * 64,
        "development_inputs_sha256": "a" * 64,
        "code_artifacts": [code_ref],
        "operational_event_reference": "synthetic-engineering-only-not-custody-permission",
    }
    sealed = population.seal_training_freeze(protocol, prepared, fitted, provenance, store)
    assert sealed["strong_claim_authorized"] is False
    assert sealed["selected_development_model"] is None
    freeze = population.read_artifact(root, sealed["prediction_freeze"], protocol, "development_predictions")
    return sealed, freeze, store


def _engineering_development_means(population, protocol, prepared):
    counts = {row["group_id"]: row["source_rows"] for row in population.group_specs(protocol)}
    rows = []
    for coordinate in prepared["development_coordinates"]:
        count = counts[coordinate["group_id"]]
        outlier = coordinate["frame_1based"] == 73
        valid = 1 if outlier else count
        rows.append({**coordinate, "observed_mean": 3.8 if outlier else 1.4,
                     "prefix_cells": count, "valid_cells": valid, "source_missing_cells": count - valid,
                     "invalid_cells": 0, "cell_sample_sd": None if outlier else 0.0})
    return rows


def test_sealed_population_report_uses_frame_mse_not_mae_or_survivor_cell_weights(tmp_path, population, protocol):
    prepared = _prepared(population, protocol, levels=(1.4, 1.4, 1.4), train_mu=1.5, development_mu=(1.5, 1.5))
    fitted = population.fit_training(prepared, protocol)
    sealed, freeze, store = _seal_engineering_fit(population, protocol, prepared, fitted, tmp_path, "complete")
    assert set(sealed["selected_checkpoints"]) == set(FAMILIES)
    assert all(sealed["selected_checkpoints"].values())
    assert len(freeze["predictions"]) == 288
    assert all(ref["role"] != "development_responses" for ref in store.references)
    # These synthetic outcome means are constructed only after every family is
    # actually serialized. The one-cell frame must still get its full 1/24 weight.
    observations = _engineering_development_means(population, protocol, prepared)
    report = population.build_development_report(protocol, freeze, observations,
                                                  sealed["prediction_freeze"]["sha256"], "c" * 64)
    assert report["status"] == "complete_development_diagnostics"
    assert report["selected_development_model"] == "prefix_persistence"
    assert report["strong_claim_authorized"] is False and report["final_test_defined"] is False
    groups = report["group_evaluations"]
    assert len(groups) == 12 and all(len(group["frames"]) == 24 for group in groups)
    for group in groups:
        if group["model_id"] == "prefix_persistence":
            assert group["mean_mse"] == pytest.approx(0.23)
        elif group["model_id"] == "training_constant":
            assert group["mean_mse"] == pytest.approx(0.24, abs=1e-6)
            assert group["mae"] == pytest.approx(0.1, abs=1e-6)
    # A missing target mean is not a zero or a removable hard frame.
    missing = deepcopy(observations)
    missing[0].update(observed_mean=None, valid_cells=0,
                      source_missing_cells=missing[0]["prefix_cells"], cell_sample_sd=None)
    incomplete = population.build_development_report(protocol, freeze, missing,
                                                      sealed["prediction_freeze"]["sha256"], "d" * 64)
    assert incomplete["status"] == "incomplete_observation_coverage"
    assert incomplete["selected_development_model"] is None
    assert len(incomplete["group_evaluations"]) == 12
    assert all(len(group["frames"]) == 24 for group in incomplete["group_evaluations"])
    assert sum(row["status"] == "missing_observation"
               for group in incomplete["group_evaluations"] for row in group["frames"]) == 6


def test_failed_family_remains_in_the_real_sealed_comparison(tmp_path, population, protocol):
    prepared = _prepared(population, protocol, levels=(1.4, 1.4, 1.4), train_mu=1.5, development_mu=(1.5, 1.5))
    fitted = population.fit_training(prepared, protocol)
    failed = fitted["transient_response"]
    failed.update(status="failed", parameters=None, selected_start_index=None, training_mean_mse=None,
                  reason="synthetic recorded optimizer failure; no biological inference")
    for attempt in failed["attempts"]:
        attempt.update(status="failed", training_mean_mse=None, reason=failed["reason"])
    sealed, freeze, _ = _seal_engineering_fit(population, protocol, prepared, fitted, tmp_path, "failed")
    assert len([row for row in freeze["attempts"] if row["model_id"] == "transient_response"]) == 8
    affected = [row for row in freeze["predictions"] if row["model_id"] == "transient_response"]
    assert len(affected) == 48 and all(row["status"] == "failed" and row["mean"] is None for row in affected)
    report = population.build_development_report(
        protocol, freeze, _engineering_development_means(population, protocol, prepared),
        sealed["prediction_freeze"]["sha256"], "c" * 64,
    )
    assert report["status"] == "incomplete_model_execution"
    assert report["selected_development_model"] is None
    assert len(report["group_evaluations"]) == 12
    omitted = deepcopy(freeze)
    omitted["predictions"] = [row for row in omitted["predictions"] if row["model_id"] != "transient_response"]
    with pytest.raises(population.DevelopmentError):
        population.validate_freeze(protocol, omitted)


@pytest.fixture(scope="module")
def runner():
    spec = importlib.util.spec_from_file_location(
        "population_adversarial_runner", ROOT / "scripts/run_native_population_development.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("namespace", ["sealed/raw_mixed", "development_exports/granados-sfp1-published-cohort-development-01"])
def test_runner_cannot_read_a_mixed_source_as_allowlist_even_if_named_train(
    tmp_path, monkeypatch, capsys, population, protocol, runner, namespace,
):
    protocol_path = tmp_path / "data/native_law_v2/development_protocol.json"
    protocol_path.parent.mkdir(parents=True)
    protocol_path.write_bytes((ROOT / "data/native_law_v2/development_protocol.json").read_bytes())
    raw_path = tmp_path / "data/native_law_v2/granados" / namespace / "fig1_train_mixed.json"
    raw_path.parent.mkdir(parents=True)
    raw_path.write_text(json.dumps({"rep1": {"GFP": {"max5": [[9.0]], "median": [[2.0]]}}}))
    forbidden = raw_path.resolve()
    original_read = Path.read_bytes
    opened = []

    def guarded_read(path):
        resolved = path.resolve()
        assert resolved != forbidden, (
            "Runner opened an unapproved mixed measurement file as custody metadata. "
            "Metadata admission must reject projection/raw namespaces before reading values; "
            "a filename containing 'train' or the --allowlist option is not custody approval."
        )
        opened.append(resolved)
        return original_read(path)

    monkeypatch.setattr(Path, "read_bytes", guarded_read)
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    assert runner.main(["fit", "--protocol", str(protocol_path), "--allowlist", str(raw_path),
                        "--output", "outputs/native_population_development/blocked"]) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["strong_claim_authorized"] is False
    assert forbidden not in opened
    assert not (tmp_path / "outputs").exists()


def test_real_custodial_projection_adapter_keeps_unreleased_responses_out_of_fitting(tmp_path, population, protocol, custody):
    class PrefixOnlyIntensity(list):
        def __getitem__(self, index):
            assert isinstance(index, int) and 28 <= index + 1 <= 47, (
                "Development conditioning export consulted a response intensity"
            )
            return super().__getitem__(index)

    selected = {}
    for spec in population.group_specs(protocol):
        group = {field: [[None] * 97 for _ in range(spec["source_rows"])] for field in ("max5", "median", "times")}
        group["origin"] = 49
        for cell in range(spec["source_rows"]):
            denominator = 1.0 + cell % 3
            for frame in (*range(28, 48), *range(50, 74)):
                is_prefix = frame < 49
                group["times"][cell][frame - 1] = float(frame - 48) if is_prefix else (frame - 49) * 2.5
                group["median"][cell][frame - 1] = denominator
                group["max5"][cell][frame - 1] = denominator * (1.5 if is_prefix else 1.4)
        group["max5"][2][30] = None  # a real source null, not a missing serialized row
        if spec["split"] == "development":
            for field in ("max5", "median"):
                group[field] = [PrefixOnlyIntensity(row) for row in group[field]]
        selected[spec["source_pointer"].removeprefix("/")] = group
    store = population.ArtifactStore(tmp_path, "outputs/native_population_development/source-adapter")
    projections = {}
    for phase in ("training_values", "development_inputs"):
        payload = custody.dev_projection_payload(selected, phase)
        replay = custody.dev_verify_payload(payload, selected, phase)
        audit = store.write(f"{phase}_engineering_binding.json", {
            "fixture_scope": "synthetic_engineering_not_biological_source_verification", "replay": replay,
        }, "source_binding_audit")
        projection = {
            "document_type": "source_bound_development_projection", "protocol_id": protocol["protocol_id"],
            "phase": phase, "approval_sha256": "b" * 64,
            "source_asset_id": ASSET, "source_sha256": protocol["source"]["source_sha256"],
            **payload, "source_binding_audit": audit,
        }
        population.validate_projection(protocol, projection, phase)
        projections[phase] = projection
    training, conditioning = projections["training_values"], projections["development_inputs"]
    prepared = population.prepare_projections(protocol, training, conditioning,
                                               population.digest(training), population.digest(conditioning))
    fitted = population.fit_training(prepared, protocol)
    predictions = population.predict_development(prepared, fitted, protocol)[0]
    for spec in population.group_specs(protocol)[3:]:
        group = selected[spec["source_pointer"].removeprefix("/")]
        for cell in range(spec["source_rows"]):
            for frame in range(50, 74):
                group["max5"][cell][frame - 1] = 1e12
                group["median"][cell][frame - 1] = 0.0
    poisoned_conditioning = custody.dev_projection_payload(selected, "development_inputs")
    assert poisoned_conditioning == {key: conditioning[key] for key in ("group_ids", "source_cells", "source_times")}
    again = population.prepare_projections(protocol, training, {**conditioning, **poisoned_conditioning},
                                            population.digest(training), population.digest(conditioning))
    repeated = population.fit_training(again, protocol)
    assert again["H_train"] == prepared["H_train"] and again["S_train"] == prepared["S_train"]
    for family in FAMILIES:
        for field in ("parameters", "bounds", "selected_start_index", "training_mean_mse"):
            assert repeated[family][field] == fitted[family][field]
    assert population.predict_development(again, repeated, protocol)[0] == predictions


def test_registered_family_predictions_survive_the_modules_own_canonical_json(population, protocol):
    prepared = _prepared(population, protocol)
    fitted = population.fit_training(prepared, protocol)
    expected = population.predict_development(prepared, fitted, protocol)
    restored = population.strict_json(population.canonical_bytes(fitted))
    assert set(restored) == set(FAMILIES)
    assert population.predict_development(prepared, restored, protocol) == expected, (
        "Canonical JSON storage must not invalidate the same six complete families. "
        "Validate the exact family set and iterate in protocol order; do not depend "
        "on incidental object-key insertion order."
    )
