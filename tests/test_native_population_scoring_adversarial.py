"""Independent phase-C scoring review on explicitly authorized projections.

The numerical oracle below does not call the production aggregation, mean-law,
normalization, optimizer, or scoring helpers. Negative probes mutate memory or
pytest temporary copies only. Mixed/private source files and rep6 stay unopened.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path

import pytest

from ystwin.analysis import native_population_development as frozen
from ystwin.analysis import native_population_scoring as scoring

from test_population_historical_sources import historical_population_root as historical_population_root


ROOT = Path(__file__).resolve().parents[1]
RUN = "outputs/native_population_development/native_training_01_phase_c_01_scoring"
MANIFEST = "data/native_law_v2/granados/development/phase_c_01/release_manifest.json"
MANIFEST_SHA = "0e2d18b1d239f33fd21165219c512b52546635b733ec47d67c85f726a919ec82"
FREEZE_SHA = "e03f6b9407b29e7e76da7a0be566ce9036c5cf27ac65bfaa9a4d3d7e009b5753"
RESPONSE_SHA = "9090f1e8ced8094f4ed445f5285bf5f85a7c457b07c34c28e9fa0cb0e0ea6178"
MODULE_SHA = "654dd5729d204c5963f5e51576a4830b61f736296c877ffc419efa77d009e2a0"
PROTOCOL_SHA = "0f53f40699ebd0db6df62f00e0a0e6bfc5b56b374ec84a6e36846a8488e6dad6"
SUMMARY_SHA = "3b2093a057f5a6819c93138fc22c7b973a4235410f99b916e74acfaa434dd7ca"
REPORT_SHA = "a98fd9b755c0e8a21c43c05a07b96f76b7c780e25c3686901fedbb6924db72ef"
ASSET = "1b38366b-7473-4ddf-9cf9-26e625533686"
PARAMETERS = {
    "prefix_persistence": (), "training_constant": ("c",), "prefix_offset": ("a",),
    "prefix_affine_time": ("a", "b"), "first_order_response": ("A", "tau"),
    "transient_response": ("A", "tau_off", "tau_on"),
}
PREFIX = tuple(range(28, 48))
RESPONSE = tuple(range(50, 74))


def _bytes(relative, expected_sha=None, *, root=ROOT):
    relative = Path(relative)
    assert not relative.is_absolute() and ".." not in relative.parts
    path = root / relative
    for forbidden in (root / "data/native_law_v2/granados/sealed", root / "data/native_law_v2/granados/development/private"):
        assert not path.resolve().is_relative_to(forbidden.resolve()), "Independent review cannot read raw/private source data"
    raw = path.read_bytes()
    if expected_sha is not None:
        assert hashlib.sha256(raw).hexdigest() == expected_sha, str(relative)
    return raw


def _json(relative, expected_sha=None):
    return json.loads(_bytes(relative, expected_sha))


def _bound(reference):
    raw = _bytes(reference["path"], reference["sha256"])
    assert len(raw) == reference["bytes"]
    return json.loads(raw)


def _finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def _mean(values):
    return math.fsum(values) / len(values) if values else None


def _sd(values):
    if len(values) < 2:
        return None
    center = _mean(values)
    return math.sqrt(math.fsum((value - center) ** 2 for value in values) / (len(values) - 1))


def _quantile(values, probability):
    values = sorted(values)
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    fraction = position - lower
    return values[lower] * (1 - fraction) + values[min(lower + 1, len(values) - 1)] * fraction


def _valid(row, prefix=False):
    x, denominator, time = (row[key] for key in ("max5", "median", "relative_time_min"))
    if not all(_finite(value) for value in (x, denominator, time)):
        return False
    return x >= 0 and denominator > 0 and (time < 0 if prefix else time > 0) and math.isfinite(x / denominator)


def _index_cells(projection, specifications, windows):
    specs = {row["group_id"]: row for row in specifications}
    assert projection["group_ids"] == list(specs)
    split = "train" if projection["phase"] == "training_values" else "development"
    expected = {(gid, cell, frame) for gid, spec in specs.items() for cell in range(spec["source_rows"]) for frame in windows}
    indexed = {}
    for row in projection["source_cells"]:
        assert row["split"] == split
        key = row["group_id"], row["cell_index"], row["frame_1based"]
        assert key in expected and key not in indexed
        gid, cell, frame = key
        parent = specs[gid]["source_pointer"]
        assert gid == f"{ASSET}:{parent}"
        for field, pointer in (("max5", "GFP/max5"), ("median", "GFP/median"), ("relative_time_min", "general/times")):
            pointer_key = "time_pointer" if field == "relative_time_min" else f"{field}_pointer"
            assert row[pointer_key] == f"{parent}/{pointer}/{cell}/{frame - 1}"
            assert row[field] is None or _finite(row[field])
        assert row["window"] == ("prefix" if frame in PREFIX else "response")
        indexed[key] = row
    assert set(indexed) == expected
    return indexed


def _prefix(spec, indexed):
    gid = spec["group_id"]
    cells = [cell for cell in range(spec["source_rows"]) if all(_valid(indexed[gid, cell, frame], prefix=True) for frame in PREFIX)]
    assert len(cells) >= 2
    per_frame = [[indexed[gid, cell, frame]["max5"] / indexed[gid, cell, frame]["median"] for cell in cells] for frame in PREFIX]
    return {"cells": cells, "mu": _mean([_mean(values) for values in per_frame]), "s": _mean([_sd(values) for values in per_frame])}


def _frames(spec, indexed, prefix, timing):
    gid = spec["group_id"]
    frames = []
    for frame in RESPONSE:
        coordinate_times, ratios = [], []
        missing, invalid = 0, 0
        for cell in prefix["cells"]:
            key = gid, cell, frame
            row = indexed[key]
            time = timing[key]
            assert row["relative_time_min"] == time
            if _finite(time) and time > 0:
                coordinate_times.append(time)
            if any(row[field] is None for field in ("max5", "median", "relative_time_min")):
                missing += 1
            elif not _valid(row):
                invalid += 1
            else:
                ratios.append(row["max5"] / row["median"])
        frames.append({
            "group_id": gid, "frame_1based": frame, "time_min": _mean(coordinate_times),
            "observed_mean": _mean(ratios), "prefix_cells": len(prefix["cells"]),
            "valid_cells": len(ratios), "source_missing_cells": missing, "invalid_cells": invalid,
            "cell_sample_sd": _sd(ratios),
        })
        assert len(ratios) + missing + invalid == len(prefix["cells"])
    return frames


def _model_mean(name, parameters, mu, time):
    if name == "prefix_persistence":
        raw = mu
    elif name == "training_constant":
        raw = parameters["c"]
    elif name == "prefix_offset":
        raw = mu + parameters["a"]
    elif name == "prefix_affine_time":
        raw = mu + parameters["a"] + parameters["b"] * time / 60.0
    elif name == "first_order_response":
        raw = mu + parameters["A"] * (-math.expm1(-time / parameters["tau"]))
    else:
        assert name == "transient_response"
        raw = mu + parameters["A"] * (-math.expm1(-time / parameters["tau_on"])) * math.exp(-time / parameters["tau_off"])
    return max(0.0, raw), raw <= 0.0


def _close(actual, expected):
    if expected is None:
        assert actual is None
    else:
        assert actual == pytest.approx(expected, rel=1e-10, abs=1e-12)


@pytest.fixture(scope="module")
def admitted(historical_population_root):
    return scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)


@pytest.fixture(scope="module")
def real():
    manifest = _json(MANIFEST, MANIFEST_SHA)
    protocol = _json("data/native_law_v2/development_protocol.json", PROTOCOL_SHA)
    seal = _bound(manifest["prediction_freeze"])
    assert manifest["prediction_freeze"]["sha256"] == FREEZE_SHA
    assert manifest["response_projection"]["sha256"] == RESPONSE_SHA
    training_ref = seal["environment"]["provenance"]["training_projection"]
    training = _bound(training_ref)
    conditioning = _bound(manifest["previously_released_development_inputs"])
    responses = _bound(manifest["response_projection"])
    training_specs, dev_specs = protocol["partition"]["train"], protocol["partition"]["development"]
    assert [spec["group_id"] for spec in training_specs + dev_specs] == [f"{ASSET}:/rep{i}" for i in range(1, 6)]
    train_index = _index_cells(training, training_specs, PREFIX + RESPONSE)
    prefix_index = _index_cells(conditioning, dev_specs, PREFIX)
    response_index = _index_cells(responses, dev_specs, RESPONSE)
    prefixes = {spec["group_id"]: _prefix(spec, train_index if spec in training_specs else prefix_index) for spec in training_specs + dev_specs}
    times = {}
    for row in conditioning["source_times"]:
        key = row["group_id"], row["cell_index"], row["frame_1based"]
        assert key not in times and row["group_id"] in {spec["group_id"] for spec in dev_specs}
        spec = next(spec for spec in dev_specs if spec["group_id"] == row["group_id"])
        assert row["time_pointer"] == f"{spec['source_pointer']}/general/times/{row['cell_index']}/{row['frame_1based'] - 1}"
        times[key] = row["relative_time_min"]
    assert set(times) == set(response_index)
    train_times = {key: row["relative_time_min"] for key, row in train_index.items() if key[2] in RESPONSE}
    train_frames = [row for spec in training_specs for row in _frames(spec, train_index, prefixes[spec["group_id"]], train_times)]
    dev_frames = [row for spec in dev_specs for row in _frames(spec, response_index, prefixes[spec["group_id"]], times)]
    assert len(train_frames) == 72 and len(dev_frames) == 48
    assert all(row["observed_mean"] is not None and row["time_min"] is not None for row in train_frames + dev_frames)
    return {
        "manifest": manifest, "protocol": protocol, "seal": seal, "training": training,
        "conditioning": conditioning, "responses": responses, "prefixes": prefixes,
        "train_frames": train_frames, "dev_frames": dev_frames,
        "summary": _json(f"{RUN}/scoring_summary.json", SUMMARY_SHA),
        "report": _json(f"{RUN}/development_report.json", REPORT_SHA),
    }


def test_real_authorized_metrics_recompute_from_individual_ratios_and_frozen_predictions(real, admitted, monkeypatch):
    seal, prefixes = real["seal"], real["prefixes"]
    for prefix in seal["prefix_summaries"]:
        expected = prefixes[prefix["group_id"]]
        assert prefix["eligible_cell_indices"] == expected["cells"]
        _close(prefix["mu"], expected["mu"])
        _close(prefix["s"], expected["s"])
    frames = {(row["group_id"], row["frame_1based"]): row for row in real["dev_frames"]}
    predictions = {(row["model_id"], row["group_id"], row["frame_1based"]): row for row in seal["predictions"]}
    expected_keys = {(model, group, frame) for model in PARAMETERS for group, frame in frames}
    assert len(seal["predictions"]) == 288 and set(predictions) == expected_keys
    assert Counter(row["model_id"] for row in seal["predictions"]) == {name: 48 for name in PARAMETERS}
    assert all(row["status"] == "ok" for row in predictions.values())
    evaluations = {(group["model_id"], group["group_id"]): group for group in real["report"]["group_evaluations"]}
    assert len(evaluations) == len(real["report"]["group_evaluations"]) == 12
    metrics = {}
    for model in PARAMETERS:
        group_metrics = []
        for spec in real["protocol"]["partition"]["development"]:
            gid = spec["group_id"]
            group = evaluations[model, gid]
            assert len(group["frames"]) == 24
            assert [row["frame_1based"] for row in group["frames"]] == list(RESPONSE)
            errors = []
            for recorded in group["frames"]:
                frame = frames[gid, recorded["frame_1based"]]
                prediction = predictions[model, gid, recorded["frame_1based"]]
                _close(prediction["time_min"], frame["time_min"])
                for field in ("observed_mean", "time_min", "prefix_cells", "valid_cells", "source_missing_cells", "invalid_cells", "cell_sample_sd"):
                    _close(recorded[field], frame[field])
                _close(recorded["predicted_mean"], prediction["mean"])
                error = prediction["mean"] - frame["observed_mean"]
                _close(recorded["squared_error"], error * error)
                _close(recorded["absolute_error"], abs(error))
                assert recorded["status"] == "quantified"
                errors.append(error)
            mse, mae = _mean([error * error for error in errors]), _mean([abs(error) for error in errors])
            for field, value in (("mean_mse", mse), ("rmse", math.sqrt(mse)), ("mae", mae)):
                _close(group[field], value)
                _close(real["summary"]["models"][model]["groups"][gid][field], value)
            group_metrics.append((mse, mae))
        mse, mae = _mean([row[0] for row in group_metrics]), _mean([row[1] for row in group_metrics])
        metrics[model] = {"mean_mse": mse, "rmse": math.sqrt(mse), "mae": mae}
        for field, value in metrics[model].items():
            _close(real["summary"]["models"][model][field], value)
        _close(real["summary"]["models"][model]["scaled_rmse"], math.sqrt(mse) / seal["S_train"])
        _close(real["summary"]["models"][model]["scaled_mae"], mae / seal["S_train"])
    best = min(row["mean_mse"] for row in metrics.values())
    ties = [name for name in PARAMETERS if metrics[name]["mean_mse"] <= best + 1e-12 * (1 + abs(best))]
    selected = min(ties, key=lambda name: (len(PARAMETERS[name]), name))
    assert real["summary"]["selected_development_model"] == real["report"]["selected_development_model"] == selected

    def no_refit(*args, **kwargs):
        raise AssertionError("Scoring may not refit parameters or retry optimization after response access")

    monkeypatch.setattr(frozen, "fit_training", no_refit)
    monkeypatch.setattr(frozen, "fit_family", no_refit)
    monkeypatch.setattr(frozen, "execute_scoring", no_refit)
    loaded = scoring.load_released_observations(ROOT, admitted)
    assert len(loaded) == 48
    for row in loaded:
        oracle = frames[row["group_id"], row["frame_1based"]]
        for field in ("time_min", "observed_mean", "prefix_cells", "valid_cells", "source_missing_cells", "invalid_cells", "cell_sample_sd"):
            _close(row[field], oracle[field])
    # The API intentionally preserves the exact presealed floating-point time
    # representation; the independent fsum oracle has already checked its value.
    report, summary = scoring.score_frozen_means(admitted, loaded)
    assert report["selected_development_model"] == summary["selected_development_model"] == selected
    for model, values in metrics.items():
        for field, value in values.items():
            _close(summary["models"][model][field], value)
    print(json.dumps({"independent_metrics": metrics, "selected": selected,
                      "cohort_sizes": {key: len(value["cells"]) for key, value in prefixes.items()},
                      "response_valid_cell_range": [min(row["valid_cells"] for row in frames.values()), max(row["valid_cells"] for row in frames.values())]}, sort_keys=True))


def test_presealed_training_selection_normalizers_bounds_and_adapter_disclosure(real, historical_population_root):
    seal = real["seal"]
    prefixes = real["prefixes"]
    train = real["train_frames"]
    training_mu = [prefixes[spec["group_id"]]["mu"] for spec in real["protocol"]["partition"]["train"]]
    H = max(1.0, *map(abs, training_mu), *(abs(row["observed_mean"]) for row in train))
    residual_means = [row["observed_mean"] - prefixes[row["group_id"]]["mu"] for row in train]
    S = _quantile(residual_means, 0.75) - _quantile(residual_means, 0.25)
    _close(seal["H_train"], H)
    _close(seal["S_train"], S)
    _close(real["summary"]["H_train"], H)
    _close(real["summary"]["S_train"], S)
    assert len(seal["attempts"]) == 41
    warnings, selected_starts, map_counts = {}, {}, {}
    predicted = {(row["model_id"], row["group_id"], row["frame_1based"]): row for row in seal["predictions"]}
    for model, names in PARAMETERS.items():
        attempts = [row for row in seal["attempts"] if row["model_id"] == model]
        assert [row["start_index"] for row in attempts] == list(range(8 if names else 1))
        successful = []
        for attempt in attempts:
            assert attempt["seed"] == 201801 + attempt["start_index"] and 0 <= attempt["loss_evaluations"] <= 2000
            if attempt["status"] != "ok":
                assert attempt["reason"] and attempt["checkpoint"] is None
                continue
            checkpoint = _bound(attempt["checkpoint"])
            assert checkpoint["parameters"] == attempt["parameters"] and set(checkpoint["parameters"]) == set(names)
            losses = []
            for row in train:
                value, _ = _model_mean(model, checkpoint["parameters"], prefixes[row["group_id"]]["mu"], row["time_min"])
                losses.append((value - row["observed_mean"]) ** 2)
            _close(attempt["training_mean_mse"], _mean(losses))
            _close(checkpoint["training_mean_mse"], _mean(losses))
            successful.append(attempt)
        best = min(row["training_mean_mse"] for row in successful)
        tied = [row for row in successful if row["training_mean_mse"] <= best + 1e-12 * (1 + abs(best))]
        selected = min(tied, key=lambda row: (tuple(row["parameters"][name] for name in names), row["start_index"]))
        assert seal["environment"]["selected_checkpoints"][model] == selected["checkpoint"]
        selected_starts[model] = selected["start_index"]
        parameters = selected["parameters"]
        hits = {}
        for name, value in parameters.items():
            low, high = (2.5, 480.0) if name in {"tau", "tau_on", "tau_off"} else (0.0 if name == "c" else -10 * H, 10 * H)
            assert low <= value <= high
            if math.isclose(value, low, rel_tol=1e-9, abs_tol=1e-9):
                hits[name] = "lower"
            elif math.isclose(value, high, rel_tol=1e-9, abs_tol=1e-9):
                hits[name] = "upper"
        if hits:
            warnings[model] = hits
        map_counts[model] = {
            name: sum(_model_mean(model, parameters, prefixes[row["group_id"]]["mu"], row["time_min"])[1] for row in rows)
            for name, rows in (("training", train), ("development", real["dev_frames"]))
        }
        for row in real["dev_frames"]:
            value, _ = _model_mean(model, parameters, prefixes[row["group_id"]]["mu"], row["time_min"])
            _close(predicted[model, row["group_id"], row["frame_1based"]]["mean"], value)
    assert real["summary"]["nonnegative_output_map_bindings"] == map_counts
    assert real["summary"]["parameter_bound_warnings"] == warnings
    warning_report = _bound(real["summary"]["warnings"])
    assert warning_report["parameter_bound_warnings"] == warnings
    assert warning_report["bounds_or_predictions_changed"] is False
    assert warnings.get("transient_response", {}).get("A") == "lower"
    _bytes("src/ystwin/analysis/native_population_development.py", MODULE_SHA, root=historical_population_root)
    for reference in real["manifest"]["frozen_executor_artifacts"]:
        raw = _bytes(reference["path"], reference["sha256"])
        assert len(raw) == reference["bytes"]
        if Path(reference["path"]).name == "native_population_development.py":
            assert reference["sha256"] == MODULE_SHA
    intent = _bound(real["summary"]["scoring_access_intent"])
    assert intent["prediction_freeze"] == real["manifest"]["prediction_freeze"]
    assert intent["created_after_prediction_sealing"] is True and intent["created_after_response_availability"] is True
    assert intent["preregistered_new_biology"] is False and intent["scientific_metric_or_prediction_change"] is False
    assert intent["unchanged_frozen_scientific_code"] == real["manifest"]["frozen_executor_artifacts"]
    for reference in intent["adapter_code"]:
        _bytes(reference["path"], reference["sha256"])
        assert Path(reference["path"]).is_relative_to(Path(RUN) / "adapter_code")
    assert real["summary"]["strong_claim_authorized"] is False and real["summary"]["final_test_defined"] is False
    print(json.dumps({"independent_H_train": H, "independent_S_train": S, "selected_training_starts": selected_starts, "retained_bound_warnings": warnings}, sort_keys=True))


def test_valid_distinct_approvals_are_linked_and_rotation_cannot_masquerade_as_training(admitted, real):
    training = admitted.training.data_approval
    executor = admitted.training.executor_approval
    response = admitted.release.manifest["custody_decision"]
    assert len({training["sha256"], executor["sha256"], response["sha256"]}) == 3
    assert real["training"]["approval_sha256"] == real["conditioning"]["approval_sha256"] == real["seal"]["approval_sha256"] == training["sha256"]
    assert real["responses"]["approval_sha256"] == response["sha256"]
    assert scoring.validate_lineage(admitted) is True
    altered = deepcopy(admitted)
    altered = replace(altered, training=replace(altered.training, data_approval=deepcopy(response), executor_approval=deepcopy(training)))
    altered.release.manifest.update(original_training_data_approval=deepcopy(response), executor_refresh_approval=deepcopy(training), custody_decision=deepcopy(executor))
    # All three references remain real, valid and distinct. It is the recorded
    # training seal, not reference shape or inequality alone, that fixes roles.
    with pytest.raises(scoring.ScoringError, match="prediction freeze must retain"):
        scoring.validate_lineage(altered)


@pytest.mark.parametrize("artifact", ["freeze", "response", "current_predictor", "saved_predictor"])
def test_bound_artifact_bytes_cannot_change_even_when_json_values_or_code_meaning_do_not(admitted, artifact, monkeypatch, historical_population_root):
    manifest = admitted.release.manifest
    relative = {
        "freeze": manifest["prediction_freeze"]["path"],
        "response": manifest["response_projection"]["path"],
        "current_predictor": "src/ystwin/analysis/native_population_development.py",
        "saved_predictor": manifest["frozen_executor_artifacts"][0]["path"],
    }[artifact]
    target = (historical_population_root / relative).resolve()
    original_read = Path.read_bytes
    affected = []

    def changed_bytes(path):
        raw = original_read(path)
        if path.resolve() == target:
            affected.append(path)
            return raw + b"\n"  # valid JSON/Python still; must fail actual byte binding
        return raw

    monkeypatch.setattr(Path, "read_bytes", changed_bytes)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="digest|byte"):
        if artifact == "response":
            scoring.load_released_observations(historical_population_root, admitted)
        else:
            scoring.inspect_phase_c(historical_population_root, MANIFEST, MANIFEST_SHA)
    assert affected


def _response_copy_with_rebound_bytes(tmp_path, admitted, response):
    """Use real temporary files and coherent byte hashes to reach semantic gates.

    This is a deliberately altered test copy, not a new authorized data release.
    No repository artifact, approval or reference is written or replaced.
    """
    altered = deepcopy(admitted)
    reference = altered.release.manifest["response_projection"]
    raw = (json.dumps(response, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    path = tmp_path / reference["path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    changed = {**reference, "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    altered.release.manifest["response_projection"] = changed
    altered.release.manifest["response_data_allowlist"] = [changed]
    altered.release.source_replay["response_projection"] = changed
    payload = {key: response[key] for key in ("group_ids", "source_cells", "source_times")}
    payload_bytes = (json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    altered.release.source_binding["ordered_payload_sha256"] = hashlib.sha256(payload_bytes).hexdigest()
    # Keep later real reads available; a missing temporary file must not be the
    # reason the negative test passes if the earlier semantic guard disappears.
    for key in ("previously_released_development_inputs", "coverage_audit"):
        ref = altered.release.manifest[key]
        copy_path = tmp_path / ref["path"]
        copy_path.parent.mkdir(parents=True, exist_ok=True)
        copy_path.write_bytes(_bytes(ref["path"], ref["sha256"]))
    return altered


@pytest.mark.parametrize("wrong_role", ["training", "executor"])
def test_response_projection_requires_its_phase_c_approval_after_valid_file_hash_checks(tmp_path, admitted, real, wrong_role):
    response = deepcopy(real["responses"])
    wrong = admitted.training.data_approval if wrong_role == "training" else admitted.training.executor_approval
    response["approval_sha256"] = wrong["sha256"]
    altered = _response_copy_with_rebound_bytes(tmp_path, admitted, response)
    assert scoring.validate_lineage(altered) is True
    with pytest.raises(scoring.ScoringError, match="response projection must bind its separate phase-C approval"):
        scoring.load_released_observations(tmp_path, altered)


@pytest.mark.parametrize("defect", ["duplicate_prediction", "missing_family", "duplicate_source_cell", "missing_source_cell"])
def test_duplicate_or_missing_scientific_ids_cannot_disappear_in_maps_or_averages(tmp_path, admitted, real, defect):
    if defect in {"duplicate_prediction", "missing_family"}:
        altered = deepcopy(admitted)
        rows = altered.prediction.data["predictions"]
        if defect == "duplicate_prediction":
            rows[-1] = deepcopy(rows[0])
        else:
            altered.prediction.data["predictions"] = [row for row in rows if row["model_id"] != "transient_response"]
        with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="inventory|item count|duplicate|canonical famil"):
            scoring.validate_lineage(altered)
    else:
        response = deepcopy(real["responses"])
        if defect == "duplicate_source_cell":
            response["source_cells"][-1] = deepcopy(response["source_cells"][0])
        else:
            response["source_cells"].pop()
        altered = _response_copy_with_rebound_bytes(tmp_path, admitted, response)
        assert scoring.validate_lineage(altered) is True
        with pytest.raises((scoring.ScoringError, frozen.DevelopmentError), match="source row/frame|item count|duplicate|inventory"):
            scoring.load_released_observations(tmp_path, altered)


def test_finished_score_and_training_seal_cannot_be_overwritten_or_retried_in_place(admitted, historical_population_root):
    for directory in (RUN, str(Path(admitted.prediction.reference["path"]).parent)):
        with pytest.raises(scoring.ScoringError, match="already exists|immutable"):
            scoring.new_score_store(ROOT, admitted, directory)
    # Identity rechecks are read-only; no scoring outputs are regenerated.
    _bytes(f"{RUN}/scoring_summary.json", SUMMARY_SHA)
    _bytes(f"{RUN}/development_report.json", REPORT_SHA)
    _bytes(admitted.prediction.reference["path"], FREEZE_SHA)
    _bytes("src/ystwin/analysis/native_population_development.py", MODULE_SHA, root=historical_population_root)


def test_cached_admission_cannot_rebind_unapproved_response_values_with_recomputed_hashes(tmp_path, admitted, real):
    response = deepcopy(real["responses"])
    row = next(row for row in response["source_cells"]
               if row["frame_1based"] == 50 and row["cell_index"] in real["prefixes"][row["group_id"]]["cells"] and _valid(row))
    changed_key = row["group_id"], row["frame_1based"]
    row["max5"] += 10.0 * row["median"]
    altered = _response_copy_with_rebound_bytes(tmp_path, admitted, response)
    try:
        frames = scoring.load_released_observations(tmp_path, altered)
    except (scoring.ScoringError, frozen.DevelopmentError):
        return
    original = next(frame for frame in real["dev_frames"] if (frame["group_id"], frame["frame_1based"]) == changed_key)
    accepted = next(frame for frame in frames if (frame["group_id"], frame["frame_1based"]) == changed_key)
    assert accepted["observed_mean"] == pytest.approx(original["observed_mean"], abs=1e-12), (
        "The admitted-bundle loader accepted an unapproved numeric response after "
        "its mutable cached manifest/replay fields were given caller-recomputed hashes. "
        "Rebind cached metadata to the authoritative manifest/receipt bytes before "
        "reading projection references, not just to internally consistent reference shapes. "
        f"Authorized mean={original['observed_mean']!r}, accepted mean={accepted['observed_mean']!r}."
    )
