from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace

import pytest

from ystwin.analysis.training_freeze import (
    ArtifactRef,
    InputSpec,
    TrainingContract,
    build_training_manifest,
    freeze_checkpoint,
    prepare_evaluation,
    score_interval_predictions,
    score_predictions,
    seal_interval_predictions,
    seal_predictions,
    validate_checkpoint,
    validate_interval_predictions,
    validate_predictions,
    validate_training_manifest,
)


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def input_spec(root, path, role, scope):
    return InputSpec(
        path=path,
        sha256=hashlib.sha256((root / path).read_bytes()).hexdigest(),
        role=role,
        source={
            "authority": "synthetic fixture",
            "source_id": path,
            "reference": f"fixture://{path}",
            "scope": scope,
        },
    )


def label_bytes(rows):
    return json.dumps({"quantity": "response", "unit": "a.u.", "observations": rows})


@pytest.fixture
def case(tmp_path):
    root = tmp_path / "repo"
    for directory in ("data/hog2013", "src/native", "tasks", "heldout", "artifacts"):
        (root / directory).mkdir(parents=True)
    (root / "data/hog2013/observations.csv").write_text("1\n3\n")
    (root / "src/native/learn.py").write_text("from . import dynamics\n")
    (root / "src/native/dynamics.py").write_text("def rate(x):\n    return x\n")
    (root / "tasks/conditions.json").write_text('{"level": 3}')
    (root / "tasks/structure.json").write_text('{"nodes": ["a", "b"]}')
    (root / "tasks/evaluate.py").write_text(
        'def evaluate(model, inputs):\n    return model["native_rate"] * inputs["level"]\n'
    )
    label_path = root / "heldout/observations.json"
    label_path.write_text(label_bytes([
        {"experiment_id": "exp-a", "value": 1.0},
        {"experiment_id": "exp-b", "value": 3.0},
        {"experiment_id": "exp-c", "value": 5.0},
    ]))
    contract = TrainingContract(
        root=root,
        inputs=(input_spec(
            root, "data/hog2013/observations.csv", "fit", "native_physiology"
        ),),
        code=tuple(input_spec(root, path, "training_code", "native_training_code") for path in (
            "src/native/learn.py", "src/native/dynamics.py", "tasks/evaluate.py",
        )),
    )
    task_inputs = tuple(input_spec(root, path, role, scope) for path, role, scope in (
        ("tasks/conditions.json", "task_conditions", "evaluation_task"),
        ("tasks/structure.json", "task_structure", "evaluation_task"),
        ("tasks/evaluate.py", "evaluation_code", "evaluation_code"),
    ))
    experiment_inputs = [
        {"experiment_id": identifier, "input_paths": [
            "tasks/conditions.json", "tasks/structure.json",
        ]}
        for identifier in ("exp-a", "exp-b", "exp-c")
    ]
    return SimpleNamespace(
        root=root,
        label_path=label_path,
        contract=contract,
        parameters={"native_rate": 2.0},
        config={"optimizer": {"max_steps": 7}, "seed": 11},
        task_inputs=task_inputs,
        experiment_inputs=experiment_inputs,
        predictions=[
            {"experiment_id": "exp-a", "status": "ok", "value": 2.0, "reason": None},
            {"experiment_id": "exp-b", "status": "ok", "value": 4.0, "reason": None},
            {"experiment_id": "exp-c", "status": "infeasible", "value": None,
             "reason": "native resource constraint"},
        ],
    )


def freeze(case, name="checkpoint.json", manifest=None):
    return freeze_checkpoint(
        case.root / "artifacts" / name,
        contract=case.contract,
        training_manifest=manifest or build_training_manifest(case.contract),
        parameters=case.parameters,
        config=case.config,
    )


def live(case, checkpoint):
    return {
        "checkpoint": checkpoint,
        "contract": case.contract,
        "task_inputs": case.task_inputs,
        "parameters": case.parameters,
        "config": case.config,
    }


def prepare(case, checkpoint):
    return prepare_evaluation(
        **live(case, checkpoint),
        experiment_inputs=case.experiment_inputs,
        quantity="response",
        unit="a.u.",
    )


def seal(case, checkpoint, name="predictions.json", manifest=None, predictions=None):
    return seal_predictions(
        case.root / "artifacts" / name,
        **live(case, checkpoint),
        task_manifest=manifest or prepare(case, checkpoint),
        predictions=case.predictions if predictions is None else predictions,
    )


def labels(case):
    return (input_spec(
        case.root, "heldout/observations.json", "outcome_labels", "heldout_outcomes"
    ),)


def test_manifest_records_only_exact_approved_inputs_and_dependencies(case):
    manifest = build_training_manifest(case.contract)
    assert manifest["kind"] == "training_input_manifest"
    assert manifest["dependency_mode"] == "explicit"
    assert [entry["path"] for entry in manifest["inputs"]] == [
        "data/hog2013/observations.csv"
    ]
    assert [entry["path"] for entry in manifest["code"]] == [
        "src/native/dynamics.py", "src/native/learn.py", "tasks/evaluate.py",
    ]
    entry = manifest["inputs"][0]
    assert entry["resolved_path"] == entry["path"]
    assert entry["role"] == "fit"
    assert entry["source"] == case.contract.inputs[0].source
    assert entry["sha256"] == case.contract.inputs[0].sha256
    assert manifest["sha256"] == digest({k: v for k, v in manifest.items() if k != "sha256"})
    assert "heldout" not in json.dumps(manifest)
    assert validate_training_manifest(manifest, case.contract) == manifest


def test_manifests_are_order_independent_and_defensively_copied(case):
    manifest = build_training_manifest(case.contract)
    reverse = replace(case.contract, code=tuple(reversed(case.contract.code)))
    assert build_training_manifest(reverse) == manifest
    case.contract.inputs[0].source["reference"] = "changed"
    assert manifest["inputs"][0]["source"]["reference"] != "changed"
    with pytest.raises(ValueError, match="manifest"):
        validate_training_manifest(manifest, case.contract)


@pytest.mark.parametrize("field", ["inputs", "code"])
def test_manifest_requires_native_inputs_and_explicit_code(case, field):
    with pytest.raises(ValueError, match="nonempty"):
        build_training_manifest(replace(case.contract, **{field: ()}))


@pytest.mark.parametrize("path", ["data/hog2013", "data/hog2013/*.csv", "../outside.csv"])
def test_directory_glob_and_path_escape_are_not_allowlists(case, path):
    spec = replace(case.contract.inputs[0], path=path)
    with pytest.raises((ValueError, FileNotFoundError), match="file|path|escape|allowlist"):
        build_training_manifest(replace(case.contract, inputs=(spec,)))


def test_symlink_escape_is_rejected(case, tmp_path):
    outside = tmp_path / "outside.csv"
    outside.write_text("do not read")
    alias = case.root / "data/hog2013/alias.csv"
    alias.symlink_to(outside)
    spec = replace(case.contract.inputs[0], path="data/hog2013/alias.csv")
    with pytest.raises(ValueError, match="escape"):
        build_training_manifest(replace(case.contract, inputs=(spec,)))


@pytest.mark.parametrize("link_type", ["symlink", "hardlink"])
def test_duplicate_physical_files_are_rejected(case, link_type):
    original = case.root / case.contract.inputs[0].path
    alias = case.root / "data/hog2013/alias.csv"
    if link_type == "symlink":
        alias.symlink_to(original)
    else:
        alias.hardlink_to(original)
    duplicate = replace(case.contract.inputs[0], path="data/hog2013/alias.csv")
    with pytest.raises(ValueError, match="duplicate"):
        build_training_manifest(replace(case.contract, inputs=(*case.contract.inputs, duplicate)))


def test_retargeting_an_internal_symlink_with_identical_bytes_is_detected(case):
    first = case.root / case.contract.inputs[0].path
    second = case.root / "data/hog2013/copy.csv"
    second.write_bytes(first.read_bytes())
    alias = case.root / "data/hog2013/alias.csv"
    alias.symlink_to(first)
    contract = replace(case.contract, inputs=(
        replace(case.contract.inputs[0], path="data/hog2013/alias.csv"),
    ))
    manifest = build_training_manifest(contract)
    assert manifest["inputs"][0]["resolved_path"] == "data/hog2013/observations.csv"
    alias.unlink()
    alias.symlink_to(second)
    with pytest.raises(ValueError, match="manifest"):
        validate_training_manifest(manifest, contract)


@pytest.mark.parametrize("change", [
    {"role": "outcome_labels"},
    {"source": {"authority": "fixture", "source_id": "asset", "reference": "fixture://asset",
                "scope": "heldout_outcomes"}},
    {"source": {}},
    {"sha256": "0" * 64},
])
def test_source_classification_and_pinned_hash_not_names_control_admission(case, change):
    spec = replace(case.contract.inputs[0], **change)
    with pytest.raises(ValueError, match="role|source|scope|SHA-256"):
        build_training_manifest(replace(case.contract, inputs=(spec,)))


def test_a_natively_classified_exact_path_is_not_rejected_by_filename(case):
    path = case.root / "data/hog2013/heldout_native_fit.csv"
    path.write_text("1\n3\n")
    contract = replace(case.contract, inputs=(input_spec(
        case.root, "data/hog2013/heldout_native_fit.csv", "fit", "native_physiology"
    ),))
    assert build_training_manifest(contract)["inputs"][0]["path"].endswith(
        "heldout_native_fit.csv"
    )


def test_recomputed_manifest_digest_cannot_relabel_an_authorized_source(case):
    manifest = build_training_manifest(case.contract)
    manifest["inputs"][0]["source"]["source_id"] = "replacement"
    manifest["sha256"] = digest({k: v for k, v in manifest.items() if k != "sha256"})
    with pytest.raises(ValueError, match="manifest"):
        freeze(case, manifest=manifest)
    assert not (case.root / "artifacts/checkpoint.json").exists()


def test_checkpoint_contains_generic_native_payload_and_external_byte_digest(case):
    checkpoint = freeze(case)
    raw = checkpoint.path.read_bytes()
    assert checkpoint.sha256 == hashlib.sha256(raw).hexdigest()
    saved = validate_checkpoint(checkpoint, contract=case.contract)
    assert saved["kind"] == "native_checkpoint"
    assert saved["model"] == {"parameters": case.parameters, "config": case.config}
    assert "task_manifest" not in saved
    case.parameters["native_rate"] = 100.0
    assert saved["model"]["parameters"]["native_rate"] == 2.0


@pytest.mark.parametrize("parameters", [
    {}, {"native_rate": float("nan")}, {"native_rate": float("inf")},
    {"native_rate": True}, {"native_rate": {"value": 2}}, {"product": "not a parameter"},
])
def test_fit_payload_requires_finite_native_parameter_mapping(case, parameters):
    case.parameters = parameters
    with pytest.raises(ValueError, match="parameter|finite|numeric"):
        freeze(case)


def test_non_json_config_and_ambiguous_mapping_keys_are_rejected(case):
    for config in ({"optimizer": {1: "integer", "1": "string"}}, {"x": {1, 2}},
                   {"x": float("nan")}):
        case.config = config
        with pytest.raises(ValueError, match="JSON|finite|string"):
            freeze(case)


def test_both_artifacts_are_write_once_even_for_identical_payloads(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    before = checkpoint.path.read_bytes(), prediction.path.read_bytes()
    with pytest.raises(FileExistsError):
        freeze(case)
    with pytest.raises(FileExistsError):
        seal(case, checkpoint)
    assert before == (checkpoint.path.read_bytes(), prediction.path.read_bytes())


def test_artifact_path_escape_and_existing_symlink_do_not_write(case, tmp_path):
    outside = tmp_path / "outside.json"
    outside.write_text("untouched")
    alias = case.root / "artifacts/checkpoint.json"
    alias.symlink_to(outside)
    with pytest.raises((FileExistsError, ValueError)):
        freeze(case)
    assert outside.read_text() == "untouched"
    with pytest.raises(ValueError, match="escape"):
        freeze_checkpoint(
            outside, contract=case.contract,
            training_manifest=build_training_manifest(case.contract),
            parameters=case.parameters, config=case.config,
        )


@pytest.mark.parametrize("field", ["parameters", "config"])
def test_live_model_and_configuration_mutation_are_rejected(case, field):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    if field == "parameters":
        case.parameters["native_rate"] += 1
    else:
        case.config["optimizer"]["max_steps"] += 1
    with pytest.raises(ValueError, match="model|config|parameter"):
        validate_checkpoint(
            checkpoint, contract=case.contract, parameters=case.parameters, config=case.config
        )
    with pytest.raises(ValueError, match="model|config|parameter"):
        validate_predictions(prediction, **live(case, checkpoint))


@pytest.mark.parametrize("path", ["data/hog2013/observations.csv", "src/native/dynamics.py"])
def test_native_file_and_dependency_mutation_invalidate_checkpoint_and_score(case, path):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    (case.root / path).write_text("changed")
    case.label_path.write_text("poison: not valid labels")
    with pytest.raises(ValueError, match="SHA-256"):
        validate_checkpoint(checkpoint, contract=case.contract)
    with pytest.raises(ValueError, match="SHA-256"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


@pytest.mark.parametrize("artifact_kind", ["checkpoint", "predictions"])
def test_byte_mutation_is_detected_even_when_json_meaning_is_unchanged(case, artifact_kind):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    artifact = checkpoint if artifact_kind == "checkpoint" else prediction
    artifact.path.write_bytes(artifact.path.read_bytes() + b"\n")
    case.label_path.write_text("poison")
    with pytest.raises(ValueError, match="SHA-256"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_evaluation_requires_valid_checkpoint_before_task_inputs_are_opened(case):
    checkpoint = ArtifactRef(case.root / "artifacts/absent.json", "0" * 64)
    case.task_inputs = (replace(case.task_inputs[0], path="tasks/also_absent.json"),)
    with pytest.raises(FileNotFoundError, match="absent.json"):
        prepare(case, checkpoint)


def test_task_manifest_binds_code_conditions_and_explicit_ids_to_checkpoint(case):
    checkpoint = freeze(case)
    manifest = prepare(case, checkpoint)
    assert manifest["checkpoint_sha256"] == checkpoint.sha256
    assert manifest["quantity"] == "response"
    assert manifest["unit"] == "a.u."
    assert manifest["experiment_inputs"] == case.experiment_inputs
    assert {entry["role"] for entry in manifest["inputs"]} == {
        "task_conditions", "task_structure", "evaluation_code",
    }


@pytest.mark.parametrize("problem", ["duplicate", "unknown_input", "unused_input", "empty"])
def test_experiment_manifest_rejects_ambiguous_or_unbound_inputs(case, problem):
    checkpoint = freeze(case)
    if problem == "duplicate":
        case.experiment_inputs.append(deepcopy(case.experiment_inputs[0]))
    elif problem == "unknown_input":
        case.experiment_inputs[0]["input_paths"].append("heldout/observations.json")
    elif problem == "unused_input":
        for row in case.experiment_inputs:
            row["input_paths"].remove("tasks/structure.json")
    else:
        case.experiment_inputs = []
    with pytest.raises(ValueError, match="experiment|input|duplicate"):
        prepare(case, checkpoint)


def test_task_manifest_requires_explicit_evaluator_code_and_refuses_label_roles(case):
    checkpoint = freeze(case)
    original = case.task_inputs
    case.task_inputs = original[:2]
    with pytest.raises(ValueError, match="evaluation_code"):
        prepare(case, checkpoint)
    case.task_inputs = (*original, *labels(case))
    with pytest.raises(ValueError, match="role|scope"):
        prepare(case, checkpoint)


@pytest.mark.parametrize("problem", ["duplicate", "missing", "extra"])
def test_predictions_require_exactly_one_row_per_experiment(case, problem):
    checkpoint = freeze(case)
    if problem == "duplicate":
        case.predictions.append(deepcopy(case.predictions[0]))
    elif problem == "missing":
        case.predictions.pop()
    else:
        case.predictions.append({"experiment_id": "extra", "status": "ok", "value": 0,
                                 "reason": None})
    with pytest.raises(ValueError, match="duplicate|missing|extra"):
        seal(case, checkpoint)
    assert not (case.root / "artifacts/predictions.json").exists()


@pytest.mark.parametrize("change", [
    {"status": "unknown"}, {"value": float("nan")}, {"value": True},
    {"status": "failed", "value": 0, "reason": "failure"},
    {"status": "infeasible", "value": None, "reason": ""},
    {"value": None}, {"reason": "hidden warning"}, {"observed": 99},
])
def test_prediction_failures_cannot_masquerade_as_zero_filled_successes(case, change):
    checkpoint = freeze(case)
    case.predictions[0].update(change)
    with pytest.raises(ValueError, match="prediction|finite|numeric|status|reason"):
        seal(case, checkpoint)


def test_reordered_predictions_produce_identical_sealed_hashes(case):
    checkpoint = freeze(case)
    first = seal(case, checkpoint)
    case.predictions.reverse()
    second = seal(case, checkpoint, name="reordered.json")
    assert first.sha256 == second.sha256
    assert validate_predictions(first, **live(case, checkpoint))["checkpoint_sha256"] == (
        checkpoint.sha256
    )


def test_different_checkpoint_cannot_be_substituted_for_prediction_model(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    case.parameters["native_rate"] += 1
    other = freeze(case, name="other-checkpoint.json")
    with pytest.raises(ValueError, match="checkpoint"):
        validate_predictions(prediction, **live(case, other))


@pytest.mark.parametrize("stage", ["seal", "score"])
def test_task_file_mutation_invalidates_predictions_before_labels_are_read(case, stage):
    checkpoint = freeze(case)
    task_manifest = prepare(case, checkpoint)
    prediction = seal(case, checkpoint, manifest=task_manifest) if stage == "score" else None
    (case.root / "tasks/conditions.json").write_text('{"level": 99}')
    case.label_path.write_text("poison")
    with pytest.raises(ValueError, match="SHA-256"):
        if stage == "seal":
            seal(case, checkpoint, manifest=task_manifest)
        else:
            score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_poisoned_heldout_labels_do_not_change_training_or_prediction_hashes(case):
    first_manifest = build_training_manifest(case.contract)
    first_checkpoint = freeze(case)
    first_predictions = seal(case, first_checkpoint)
    case.label_path.write_text("POISON: labels must not be parsed by training or prediction")
    second_manifest = build_training_manifest(case.contract)
    second_checkpoint = freeze(case, name="second-checkpoint.json")
    second_predictions = seal(case, second_checkpoint, name="second-predictions.json")
    assert first_manifest == second_manifest
    assert first_checkpoint.sha256 == second_checkpoint.sha256
    assert first_predictions.sha256 == second_predictions.sha256
    assert "heldout/observations.json" not in first_predictions.path.read_text()


def test_scoring_joins_by_id_and_preserves_failed_rows(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    case.label_path.write_text(label_bytes([
        {"experiment_id": "exp-c", "value": 9.0},
        {"experiment_id": "exp-a", "value": 1.5},
        {"experiment_id": "exp-b", "value": 6.0},
    ]))
    result = score_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert result["checkpoint_sha256"] == checkpoint.sha256
    assert result["predictions_sha256"] == prediction.sha256
    assert result["complete"] is False
    assert (result["n_total"], result["n_scored"], result["n_failed"]) == (3, 2, 1)
    assert result["mae_on_successes"] == pytest.approx(1.25)
    assert result["rmse_on_successes"] == pytest.approx((4.25 / 2) ** 0.5)
    rows = result["rows"]
    assert [row["experiment_id"] for row in rows] == ["exp-a", "exp-b", "exp-c"]
    assert [row["observed"] for row in rows] == [1.5, 6.0, 9.0]
    assert [row["residual"] for row in rows] == [0.5, -2.0, None]
    assert rows[-1]["predicted"] is None
    assert rows[-1]["status"] == "infeasible"
    assert rows[-1]["reason"] == "native resource constraint"
    json.dumps(result, allow_nan=False)


def test_all_failed_predictions_report_no_numeric_score_instead_of_zero(case):
    checkpoint = freeze(case)
    case.predictions = [
        {"experiment_id": row["experiment_id"], "status": "failed", "value": None,
         "reason": "solver failed"}
        for row in case.experiment_inputs
    ]
    prediction = seal(case, checkpoint)
    result = score_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert result["n_scored"] == 0
    assert result["n_failed"] == 3
    assert result["mae_on_successes"] is None
    assert result["rmse_on_successes"] is None
    assert len(result["rows"]) == 3


@pytest.mark.parametrize("problem", ["duplicate", "missing", "extra", "nonnumeric", "unit"])
def test_labels_require_exact_ids_finite_values_and_matching_units(case, problem):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    payload = json.loads(case.label_path.read_text())
    if problem == "duplicate":
        payload["observations"].append(payload["observations"][0])
    elif problem == "missing":
        payload["observations"].pop()
    elif problem == "extra":
        payload["observations"].append({"experiment_id": "extra", "value": 1.0})
    elif problem == "nonnumeric":
        payload["observations"][0]["value"] = None
    else:
        payload["unit"] = "different"
    case.label_path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="duplicate|missing|extra|numeric|unit"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_label_file_hash_is_pinned_and_duplicate_ids_across_files_are_rejected(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    expected = labels(case)
    case.label_path.write_text(label_bytes([{"experiment_id": "exp-a", "value": 7}]))
    with pytest.raises(ValueError, match="SHA-256"):
        score_predictions(prediction, **live(case, checkpoint), labels=expected)
    duplicate = case.root / "heldout/duplicate.json"
    duplicate.write_bytes(case.label_path.read_bytes())
    extra = input_spec(case.root, "heldout/duplicate.json", "outcome_labels", "heldout_outcomes")
    with pytest.raises(ValueError, match="duplicate"):
        score_predictions(prediction, **live(case, checkpoint), labels=(*labels(case), extra))


@pytest.mark.parametrize("reuse", ["path", "symlink", "copy", "source_identity"])
def test_labels_cannot_be_reused_as_task_inputs_even_under_another_name(case, reuse):
    task = case.root / "tasks/conditions.json"
    label_path = "heldout/observations.json"
    if reuse == "path":
        case.experiment_inputs = [
            {"experiment_id": row["experiment_id"], "input_paths": [
                label_path, "tasks/structure.json",
            ]}
            for row in case.experiment_inputs
        ]
        changed = input_spec(case.root, label_path, "task_conditions", "evaluation_task")
    elif reuse == "symlink":
        task.unlink()
        task.symlink_to(case.label_path)
        changed = input_spec(case.root, "tasks/conditions.json", "task_conditions", "evaluation_task")
    elif reuse == "copy":
        task.write_bytes(case.label_path.read_bytes())
        changed = input_spec(case.root, "tasks/conditions.json", "task_conditions", "evaluation_task")
    else:
        changed = replace(case.task_inputs[0], source={
            **case.task_inputs[0].source, "source_id": label_path,
        })
    case.task_inputs = (changed, *case.task_inputs[1:])
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    with pytest.raises(ValueError, match="overlap|reus"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_labels_cannot_reuse_native_training_input_bytes(case):
    native = case.root / case.contract.inputs[0].path
    native.write_bytes(case.label_path.read_bytes())
    case.contract = replace(case.contract, inputs=(input_spec(
        case.root, case.contract.inputs[0].path, "fit", "native_physiology"
    ),))
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    with pytest.raises(ValueError, match="overlap|reus"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_scoring_does_not_mutate_or_refit_model_and_sealed_artifacts(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    before = deepcopy(case.parameters), deepcopy(case.config)
    artifacts = checkpoint.path.read_bytes(), prediction.path.read_bytes()
    score_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert before == (case.parameters, case.config)
    assert artifacts == (checkpoint.path.read_bytes(), prediction.path.read_bytes())


@pytest.mark.parametrize("field", ["parameters", "config"])
@pytest.mark.parametrize("stage", ["prepare", "seal", "validate", "score"])
def test_live_payload_cannot_be_omitted_at_evaluation_boundaries(case, field, stage):
    checkpoint = freeze(case)
    task_manifest = prepare(case, checkpoint)
    prediction = seal(case, checkpoint, manifest=task_manifest)
    setattr(case, field, None)
    with pytest.raises(ValueError, match="parameter|config"):
        if stage == "prepare":
            prepare(case, checkpoint)
        elif stage == "seal":
            seal(case, checkpoint, name="other.json", manifest=task_manifest)
        elif stage == "validate":
            validate_predictions(prediction, **live(case, checkpoint))
        else:
            score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_updating_allowlist_hash_cannot_rebind_a_frozen_training_input(case):
    checkpoint = freeze(case)
    (case.root / case.contract.inputs[0].path).write_text("changed")
    case.contract = replace(case.contract, inputs=(input_spec(
        case.root, case.contract.inputs[0].path, "fit", "native_physiology"
    ),))
    with pytest.raises(ValueError, match="manifest"):
        validate_checkpoint(checkpoint, contract=case.contract)


def test_updating_task_allowlist_hash_cannot_rebind_sealed_predictions(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    (case.root / "tasks/conditions.json").write_text('{"level": 17}')
    case.task_inputs = (input_spec(
        case.root, "tasks/conditions.json", "task_conditions", "evaluation_task"
    ), *case.task_inputs[1:])
    with pytest.raises(ValueError, match="manifest"):
        validate_predictions(prediction, **live(case, checkpoint))


def test_scoring_requires_a_prediction_artifact_before_opening_labels(case):
    checkpoint = freeze(case)
    missing = ArtifactRef(case.root / "artifacts/not-sealed.json", "0" * 64)
    case.label_path.write_text("poison")
    with pytest.raises(FileNotFoundError, match="not-sealed.json"):
        score_predictions(missing, **live(case, checkpoint), labels=labels(case))


@pytest.mark.parametrize("identifier", [None, 3, "", " exp-a"])
def test_experiment_identifiers_are_explicit_strings_not_positional_indices(case, identifier):
    checkpoint = freeze(case)
    case.predictions[0]["experiment_id"] = identifier
    with pytest.raises(ValueError, match="experiment_id"):
        seal(case, checkpoint)


def test_scoring_rejects_duplicate_json_fields_instead_of_last_value_wins(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    case.label_path.write_text(
        '{"quantity":"response","unit":"a.u.","observations":[],"observations":[]}'
    )
    with pytest.raises(ValueError, match="duplicate JSON field"):
        score_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_changed_labels_affect_only_score_not_model_or_prediction_artifacts(case):
    checkpoint = freeze(case)
    prediction = seal(case, checkpoint)
    first = score_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    case.label_path.write_text(label_bytes([
        {"experiment_id": "exp-a", "value": 12.0},
        {"experiment_id": "exp-b", "value": 21.0},
        {"experiment_id": "exp-c", "value": 1.0},
    ]))
    second = score_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert first["label_manifest"]["sha256"] != second["label_manifest"]["sha256"]
    assert first["mae_on_successes"] != second["mae_on_successes"]
    assert first["checkpoint_sha256"] == second["checkpoint_sha256"] == checkpoint.sha256
    assert first["predictions_sha256"] == second["predictions_sha256"] == prediction.sha256


@pytest.fixture
def interval_case(case):
    case.config.update({
        "prediction_kind": "feasible_envelope",
        "prediction_scope": "native_resource_feasibility_only",
    })
    case.intervals = [
        {"experiment_id": "exp-a", "status": "ok", "lower": 0.0, "upper": 2.0, "reason": None},
        {"experiment_id": "exp-b", "status": "ok", "lower": 1.0, "upper": 4.0, "reason": None},
        {"experiment_id": "exp-c", "status": "unsupported", "lower": None, "upper": None,
         "reason": "required native constraint absent"},
    ]
    return case


def seal_intervals(case, checkpoint, name="intervals.json", manifest=None):
    return seal_interval_predictions(
        case.root / "artifacts" / name, **live(case, checkpoint),
        task_manifest=manifest or prepare(case, checkpoint), predictions=case.intervals,
    )


def test_interval_artifact_seals_bounds_and_declared_feasibility_scope(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    saved = validate_interval_predictions(prediction, **live(case, checkpoint))
    assert set(saved) == {
        "schema_version", "kind", "checkpoint_sha256", "model_sha256", "task_manifest",
        "prediction_kind", "prediction_scope", "predictions",
    }
    assert saved["kind"] == "sealed_interval_predictions"
    assert saved["prediction_kind"] == "feasible_envelope"
    assert saved["prediction_scope"] == case.config["prediction_scope"]
    assert saved["checkpoint_sha256"] == checkpoint.sha256
    assert saved["predictions"] == case.intervals
    assert saved["task_manifest"]["quantity"] == "response"
    assert saved["task_manifest"]["unit"] == "a.u."
    assert prediction.sha256 == hashlib.sha256(prediction.path.read_bytes()).hexdigest()
    case.intervals[0]["upper"] = 100.0
    assert validate_interval_predictions(prediction, **live(case, checkpoint))["predictions"] == (
        saved["predictions"]
    )


@pytest.mark.parametrize("field,value", [
    ("prediction_kind", None), ("prediction_kind", "point"),
    ("prediction_kind", "confidence_interval"), ("prediction_scope", None),
    ("prediction_scope", ""), ("prediction_scope", 1),
])
def test_interval_kind_and_scope_must_be_declared_in_frozen_config(interval_case, field, value):
    case = interval_case
    if value is None:
        case.config.pop(field)
    else:
        case.config[field] = value
    checkpoint = freeze(case)
    with pytest.raises(ValueError, match="prediction_kind|prediction_scope"):
        seal_intervals(case, checkpoint)
    assert not (case.root / "artifacts/intervals.json").exists()


def test_feasibility_config_cannot_be_used_to_seal_a_point_estimate(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    with pytest.raises(ValueError, match="prediction_kind|point"):
        seal(case, checkpoint)


@pytest.mark.parametrize("change", [
    {"status": "unknown"}, {"lower": 3.0}, {"lower": None}, {"upper": None},
    {"lower": float("nan")}, {"upper": float("nan")},
    {"lower": -float("inf")}, {"upper": float("inf")},
    {"lower": True}, {"upper": "2.0"}, {"reason": "not a failure"},
    {"status": "failed", "lower": 0.0, "upper": None, "reason": "solver failed"},
    {"status": "unsupported", "lower": None, "upper": None, "reason": ""},
    {"status": "infeasible", "lower": None, "upper": None, "reason": None},
    {"value": 2.0}, {"point": 1.0}, {"midpoint": 1.0},
    {"lower": -1e308, "upper": 1e308},
])
def test_interval_schema_requires_finite_ordered_bounds_without_a_point(interval_case, change):
    case = interval_case
    checkpoint = freeze(case)
    case.intervals[0].update(change)
    with pytest.raises(ValueError, match="interval|finite|numeric|status|reason"):
        seal_intervals(case, checkpoint)


@pytest.mark.parametrize("field", ["lower", "upper", "status", "reason"])
def test_interval_schema_does_not_invent_missing_bounds_or_status(interval_case, field):
    case = interval_case
    checkpoint = freeze(case)
    del case.intervals[0][field]
    with pytest.raises(ValueError, match="interval.*fields"):
        seal_intervals(case, checkpoint)


@pytest.mark.parametrize("problem", ["duplicate", "missing", "extra"])
def test_interval_predictions_require_the_exact_condition_id_set(interval_case, problem):
    case = interval_case
    checkpoint = freeze(case)
    if problem == "duplicate":
        case.intervals.append(deepcopy(case.intervals[0]))
    elif problem == "missing":
        case.intervals.pop()
    else:
        case.intervals.append({**case.intervals[0], "experiment_id": "extra"})
    with pytest.raises(ValueError, match="duplicate|missing|extra"):
        seal_intervals(case, checkpoint)


def test_interval_and_scalar_artifacts_cannot_be_substituted_for_each_other(interval_case):
    case = interval_case
    interval_checkpoint = freeze(case)
    intervals = seal_intervals(case, interval_checkpoint)
    with pytest.raises(ValueError, match="kind|point"):
        validate_predictions(intervals, **live(case, interval_checkpoint))
    case.config["prediction_kind"] = "point"
    point_checkpoint = freeze(case, name="point-checkpoint.json")
    point = seal(case, point_checkpoint)
    with pytest.raises(ValueError, match="kind"):
        validate_interval_predictions(point, **live(case, point_checkpoint))
    case.predictions[0].update({"status": "unsupported", "value": None, "reason": "unsupported"})
    with pytest.raises(ValueError, match="status"):
        seal(case, point_checkpoint, name="unsupported-point.json")


def test_reordered_intervals_and_task_ids_have_identical_sealed_digests(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    first = seal_intervals(case, checkpoint)
    case.intervals.reverse()
    case.experiment_inputs.reverse()
    for row in case.experiment_inputs:
        row["input_paths"].reverse()
    second = seal_intervals(case, checkpoint, name="reordered-intervals.json")
    assert first.sha256 == second.sha256


def test_interval_artifacts_are_write_once(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    before = prediction.path.read_bytes()
    with pytest.raises(FileExistsError):
        seal_intervals(case, checkpoint)
    case.intervals[0]["upper"] += 1
    with pytest.raises(FileExistsError):
        seal_intervals(case, checkpoint)
    assert prediction.path.read_bytes() == before


@pytest.mark.parametrize("mutation", [
    "checkpoint", "intervals", "native_input", "training_code", "task_input", "evaluation_code",
    "parameters", "config", "prediction_scope",
])
def test_interval_scoring_rejects_mutation_before_reading_labels(interval_case, mutation):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    files = {
        "checkpoint": checkpoint.path, "intervals": prediction.path,
        "native_input": case.root / case.contract.inputs[0].path,
        "training_code": case.root / case.contract.code[0].path,
        "task_input": case.root / "tasks/conditions.json",
        "evaluation_code": case.root / "tasks/evaluate.py",
    }
    if mutation in files:
        path = files[mutation]
        path.write_bytes(path.read_bytes() + b"\nchanged")
    elif mutation == "parameters":
        case.parameters["native_rate"] += 1
    elif mutation == "config":
        case.config["optimizer"]["max_steps"] += 1
    else:
        case.config["prediction_scope"] = "different scope"
    case.label_path.write_text("POISON: not JSON")
    with pytest.raises(ValueError, match="SHA-256|model|config|parameter"):
        score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_interval_task_inputs_cannot_change_between_preparation_and_seal(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    task_manifest = prepare(case, checkpoint)
    (case.root / "tasks/conditions.json").write_text("changed")
    with pytest.raises(ValueError, match="SHA-256"):
        seal_intervals(case, checkpoint, manifest=task_manifest)


def test_poisoned_labels_do_not_change_frozen_interval_bounds_or_hashes(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    before = validate_interval_predictions(prediction, **live(case, checkpoint))
    case.label_path.write_text("POISON: neither training nor prediction may parse labels")
    second_checkpoint = freeze(case, name="second-checkpoint.json")
    second = seal_intervals(case, second_checkpoint, name="second-intervals.json")
    assert checkpoint.sha256 == second_checkpoint.sha256
    assert prediction.sha256 == second.sha256
    assert validate_interval_predictions(second, **live(case, second_checkpoint)) == before


@pytest.mark.parametrize("value,side", [(-1.0, "below"), (3.0, "above")])
def test_interval_scoring_aligns_ids_and_reports_conditional_containment(interval_case, value, side):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    case.label_path.write_text(label_bytes([
        {"experiment_id": "exp-c", "value": 5.0},
        {"experiment_id": "exp-b", "value": 3.0},
        {"experiment_id": "exp-a", "value": value},
    ]))
    result = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert result["kind"] == "interval_prediction_score"
    assert result["checkpoint_sha256"] == checkpoint.sha256
    assert result["predictions_sha256"] == prediction.sha256
    assert result["prediction_kind"] == "feasible_envelope"
    assert result["prediction_scope"] == case.config["prediction_scope"]
    assert result["assessment"] == "feasible_envelope_consistency_not_point_forecast"
    assert result["complete"] is False
    assert (result["n_total"], result["n_successful"], result["n_scorable"]) == (3, 2, 2)
    assert result["status_counts"] == {"ok": 2, "failed": 0, "infeasible": 0, "unsupported": 1}
    assert result["successful_fraction_total"] == pytest.approx(2 / 3)
    assert result["scorable_fraction_total"] == pytest.approx(2 / 3)
    assert result["n_contained"] == 1
    assert result["n_violated"] == 1
    assert result["containment_on_scorable"] == 0.5
    assert result["containment_fraction_total"] == pytest.approx(1 / 3)
    assert result["mean_width_on_scorable"] == 2.5
    assert result["mean_violation_on_scorable"] == 0.5
    assert result["max_violation_on_scorable"] == 1.0
    rows = result["rows"]
    assert [row["experiment_id"] for row in rows] == ["exp-a", "exp-b", "exp-c"]
    assert [row["observed"] for row in rows] == [value, 3.0, 5.0]
    assert [row["contained"] for row in rows] == [False, True, None]
    assert [row["width"] for row in rows] == [2.0, 3.0, None]
    assert [row["violation"] for row in rows] == [1.0, 0.0, None]
    assert [row["violation_side"] for row in rows] == [side, None, None]
    assert rows[-1]["reason"] == "required native constraint absent"
    assert rows[-1]["lower"] is None and rows[-1]["upper"] is None
    json.dumps(result, allow_nan=False)


def test_interval_scoring_is_invariant_to_label_row_order(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    first = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    payload = json.loads(case.label_path.read_text())
    payload["observations"].reverse()
    case.label_path.write_text(json.dumps(payload))
    second = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert first["rows"] == second["rows"]
    for key in ("n_total", "n_successful", "n_scorable", "containment_on_scorable",
                "mean_width_on_scorable", "mean_violation_on_scorable"):
        assert first[key] == second[key]


def test_all_failed_infeasible_or_unsupported_intervals_keep_total_denominator(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    case.intervals = [
        {"experiment_id": row["experiment_id"], "status": status, "lower": None, "upper": None,
         "reason": "no finite feasible envelope available"}
        for row, status in zip(case.intervals, ("failed", "infeasible", "unsupported"))
    ]
    prediction = seal_intervals(case, checkpoint)
    result = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert (result["n_total"], result["n_successful"], result["n_scorable"]) == (3, 0, 0)
    assert result["status_counts"] == {"ok": 0, "failed": 1, "infeasible": 1, "unsupported": 1}
    assert result["successful_fraction_total"] == result["scorable_fraction_total"] == 0.0
    assert result["n_contained"] == result["n_violated"] == 0
    assert result["containment_fraction_total"] == 0.0
    for key in ("containment_on_scorable", "mean_width_on_scorable", "mean_violation_on_scorable",
                "max_violation_on_scorable"):
        assert result[key] is None
    assert result["complete"] is False
    assert len(result["rows"]) == 3
    assert all(row["contained"] is None for row in result["rows"])


def test_closed_interval_endpoints_and_degenerate_intervals_remain_intervals(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    case.intervals = [
        {"experiment_id": identifier, "status": "ok", "lower": lower, "upper": upper,
         "reason": None}
        for identifier, lower, upper in (("exp-a", 1.0, 1.0), ("exp-b", 1.0, 3.0),
                                        ("exp-c", 5.0, 5.0))
    ]
    prediction = seal_intervals(case, checkpoint)
    result = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert result["complete"] is True
    assert result["containment_on_scorable"] == result["containment_fraction_total"] == 1.0
    assert [row["width"] for row in result["rows"]] == [0.0, 2.0, 0.0]
    assert all(row["contained"] for row in result["rows"])


def test_broad_interval_report_has_no_point_accuracy_fields(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    case.intervals = [
        {"experiment_id": row["experiment_id"], "status": "ok", "lower": 0.0,
         "upper": 1e12, "reason": None}
        for row in case.intervals
    ]
    prediction = seal_intervals(case, checkpoint)
    result = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert result["containment_on_scorable"] == 1.0
    assert result["mean_width_on_scorable"] == 1e12
    assert result["assessment"] == "feasible_envelope_consistency_not_point_forecast"
    assert set(result) == {
        "schema_version", "kind", "sha256", "checkpoint_sha256", "predictions_sha256",
        "label_manifest", "quantity", "unit", "prediction_kind", "prediction_scope", "assessment",
        "rows", "complete", "status_counts", "n_total", "n_successful", "n_scorable",
        "successful_fraction_total", "scorable_fraction_total", "n_contained", "n_violated",
        "containment_on_scorable", "containment_fraction_total", "mean_width_on_scorable",
        "mean_violation_on_scorable", "max_violation_on_scorable",
    }
    for row in result["rows"]:
        assert set(row) == {
            "experiment_id", "status", "lower", "upper", "reason", "observed", "label_source_path",
            "contained", "width", "violation", "violation_side",
        }


@pytest.mark.parametrize("problem", ["duplicate", "missing", "extra", "nonnumeric", "unit", "quantity"])
def test_interval_labels_keep_exact_id_value_and_quantity_guards(interval_case, problem):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    payload = json.loads(case.label_path.read_text())
    if problem == "duplicate":
        payload["observations"].append(payload["observations"][0])
    elif problem == "missing":
        payload["observations"].pop()
    elif problem == "extra":
        payload["observations"].append({"experiment_id": "extra", "value": 1.0})
    elif problem == "nonnumeric":
        payload["observations"][0]["value"] = None
    else:
        payload[problem] = "different"
    case.label_path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="duplicate|missing|extra|numeric|unit|quantity"):
        score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))


@pytest.mark.parametrize("phase", ["native", "task"])
def test_interval_scoring_preserves_label_input_disjointness(interval_case, phase):
    case = interval_case
    path = case.contract.inputs[0].path if phase == "native" else "tasks/conditions.json"
    (case.root / path).write_bytes(case.label_path.read_bytes())
    if phase == "native":
        case.contract = replace(case.contract, inputs=(input_spec(
            case.root, path, "fit", "native_physiology"
        ),))
    else:
        case.task_inputs = (input_spec(
            case.root, path, "task_conditions", "evaluation_task"
        ), *case.task_inputs[1:])
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    with pytest.raises(ValueError, match="overlap|reus"):
        score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))


def test_interval_scoring_requires_seal_and_verified_labels(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    missing = ArtifactRef(case.root / "artifacts/not-sealed.json", "0" * 64)
    expected_labels = labels(case)
    case.label_path.write_text("POISON")
    with pytest.raises(FileNotFoundError, match="not-sealed.json"):
        score_interval_predictions(missing, **live(case, checkpoint), labels=expected_labels)
    prediction = seal_intervals(case, checkpoint)
    with pytest.raises(ValueError, match="SHA-256"):
        score_interval_predictions(prediction, **live(case, checkpoint), labels=expected_labels)


def test_changed_labels_change_interval_score_without_refitting_or_resealing(interval_case):
    case = interval_case
    checkpoint = freeze(case)
    prediction = seal_intervals(case, checkpoint)
    model = deepcopy(case.parameters), deepcopy(case.config)
    before = checkpoint.path.read_bytes(), prediction.path.read_bytes()
    first = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    case.label_path.write_text(label_bytes([
        {"experiment_id": "exp-a", "value": -3.0},
        {"experiment_id": "exp-b", "value": 17.0},
        {"experiment_id": "exp-c", "value": 1.0},
    ]))
    second = score_interval_predictions(prediction, **live(case, checkpoint), labels=labels(case))
    assert first["label_manifest"]["sha256"] != second["label_manifest"]["sha256"]
    assert first["containment_on_scorable"] != second["containment_on_scorable"]
    assert first["checkpoint_sha256"] == second["checkpoint_sha256"] == checkpoint.sha256
    assert first["predictions_sha256"] == second["predictions_sha256"] == prediction.sha256
    assert model == (case.parameters, case.config)
    assert before == (checkpoint.path.read_bytes(), prediction.path.read_bytes())


@pytest.fixture(params=["scalar", "interval"])
def prediction_case(request):
    interval = request.param == "interval"
    case = request.getfixturevalue("interval_case" if interval else "case")
    case.sealer = seal_intervals if interval else seal
    case.validator = validate_interval_predictions if interval else validate_predictions
    case.scorer = score_interval_predictions if interval else score_predictions
    return case


def rehashed_manifest(manifest):
    payload = {key: value for key, value in manifest.items() if key != "sha256"}
    return {**payload, "sha256": digest(payload)}


def forged_artifact(path, payload):
    data = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode() + b"\n"
    path.write_bytes(data)
    return ArtifactRef(path, hashlib.sha256(data).hexdigest())


def test_evaluator_code_role_may_differ_but_frozen_identity_must_match(prediction_case):
    case = prediction_case
    checkpoint = freeze(case)
    task = prepare(case, checkpoint)
    saved = validate_checkpoint(checkpoint, contract=case.contract)
    native = next(entry for entry in saved["training_manifest"]["code"]
                  if entry["path"] == "tasks/evaluate.py")
    evaluator = next(entry for entry in task["inputs"] if entry["role"] == "evaluation_code")
    assert native["role"] == "training_code"
    assert evaluator["role"] == "evaluation_code"
    assert native["source"]["scope"] != evaluator["source"]["scope"]
    for field in ("path", "resolved_path", "sha256"):
        assert native[field] == evaluator[field]
    prediction = case.sealer(case, checkpoint, manifest=task)
    case.validator(prediction, **live(case, checkpoint))
    case.scorer(prediction, **live(case, checkpoint), labels=labels(case))


@pytest.mark.parametrize("repin", [False, True])
def test_uncommitted_evaluator_cannot_be_added_or_repinned_after_freeze(prediction_case, repin):
    case = prediction_case
    case.contract = replace(case.contract, code=tuple(
        entry for entry in case.contract.code if entry.path != "tasks/evaluate.py"
    ))
    checkpoint = freeze(case)
    if repin:
        (case.root / "tasks/evaluate.py").write_text("def evaluate(model, inputs):\n    return 17\n")
        case.task_inputs = (*case.task_inputs[:2], input_spec(
            case.root, "tasks/evaluate.py", "evaluation_code", "evaluation_code"
        ))
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        prepare(case, checkpoint)


def test_each_evaluator_dependency_must_be_committed_not_only_the_entrypoint(prediction_case):
    case = prediction_case
    checkpoint = freeze(case)
    (case.root / "tasks/new_helper.py").write_text("def helper(x):\n    return x\n")
    case.task_inputs = (*case.task_inputs, input_spec(
        case.root, "tasks/new_helper.py", "evaluation_code", "evaluation_code"
    ))
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        prepare(case, checkpoint)


@pytest.mark.parametrize("alias_type", ["copy", "symlink", "hardlink"])
def test_same_evaluator_bytes_do_not_authorize_an_uncommitted_path(prediction_case, alias_type):
    case = prediction_case
    checkpoint = freeze(case)
    original = case.root / "tasks/evaluate.py"
    alias = case.root / "tasks/alias.py"
    if alias_type == "copy":
        alias.write_bytes(original.read_bytes())
    elif alias_type == "symlink":
        alias.symlink_to(original)
    else:
        alias.hardlink_to(original)
    evaluator = next(spec for spec in case.task_inputs if spec.role == "evaluation_code")
    case.task_inputs = (*case.task_inputs[:2], replace(evaluator, path="tasks/alias.py"))
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        prepare(case, checkpoint)


@pytest.mark.parametrize("field", ["authority", "source_id", "reference"])
def test_evaluator_provenance_cannot_be_relabelled_after_freeze(prediction_case, field):
    case = prediction_case
    checkpoint = freeze(case)
    evaluator = case.task_inputs[-1]
    case.task_inputs = (*case.task_inputs[:-1], replace(
        evaluator, source={**evaluator.source, field: "replacement"}
    ))
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        prepare(case, checkpoint)


@pytest.mark.parametrize("stage", ["seal", "validate", "score"])
def test_manual_task_manifest_cannot_introduce_unfrozen_evaluator(prediction_case, stage):
    case = prediction_case
    checkpoint = freeze(case)
    task = prepare(case, checkpoint)
    prediction = case.sealer(case, checkpoint, manifest=task)
    payload = case.validator(prediction, **live(case, checkpoint))
    rogue_path = "tasks/late_evaluator.py"
    (case.root / rogue_path).write_text("def evaluate(model, inputs):\n    return 17\n")
    rogue = input_spec(case.root, rogue_path, "evaluation_code", "evaluation_code")
    case.task_inputs = (*case.task_inputs[:2], rogue)
    task["inputs"] = [entry for entry in task["inputs"] if entry["role"] != "evaluation_code"] + [{
        "path": rogue.path, "resolved_path": rogue.path, "sha256": rogue.sha256,
        "size_bytes": (case.root / rogue.path).stat().st_size,
        "role": rogue.role, "source": rogue.source,
    }]
    task["inputs"].sort(key=lambda entry: entry["path"])
    task = rehashed_manifest(task)
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        if stage == "seal":
            case.sealer(case, checkpoint, name="manual-seal.json", manifest=task)
        else:
            payload["task_manifest"] = task
            manual = forged_artifact(case.root / "artifacts/manual-predictions.json", payload)
            if stage == "validate":
                case.validator(manual, **live(case, checkpoint))
            else:
                case.scorer(manual, **live(case, checkpoint), labels=labels(case))


def test_uncommitted_evaluator_is_rejected_before_task_bytes_are_opened(prediction_case):
    case = prediction_case
    case.contract = replace(case.contract, code=tuple(
        entry for entry in case.contract.code if entry.path != "tasks/evaluate.py"
    ))
    checkpoint = freeze(case)
    (case.root / "tasks/conditions.json").write_text("POISON: invalid task bytes")
    with pytest.raises(ValueError, match="evaluation_code.*checkpoint"):
        prepare(case, checkpoint)


@pytest.mark.parametrize("metadata", [
    {"quantity": "response"}, {"unit": "a.u."},
    {"quantity": "response", "unit": None}, {"quantity": None, "unit": "a.u."},
    {"quantity": "", "unit": "a.u."}, {"quantity": "response", "unit": 7},
])
def test_partial_or_invalid_checkpoint_readout_metadata_is_rejected(prediction_case, metadata):
    case = prediction_case
    case.config.update(metadata)
    with pytest.raises(ValueError, match="quantity|unit"):
        freeze(case)


def test_partial_readout_in_a_manually_built_checkpoint_is_rejected(prediction_case):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "a.u."})
    checkpoint = freeze(case)
    payload = validate_checkpoint(checkpoint, contract=case.contract)
    del payload["model"]["config"]["unit"]
    manual = forged_artifact(case.root / "artifacts/partial-checkpoint.json", payload)
    with pytest.raises(ValueError, match="quantity|unit"):
        validate_checkpoint(manual, contract=case.contract)


def test_declared_checkpoint_quantity_and_unit_are_preserved_through_scoring(prediction_case):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "a.u."})
    checkpoint = freeze(case)
    task = prepare(case, checkpoint)
    prediction = case.sealer(case, checkpoint, manifest=task)
    result = case.scorer(prediction, **live(case, checkpoint), labels=labels(case))
    assert (task["quantity"], task["unit"]) == ("response", "a.u.")
    assert (result["quantity"], result["unit"]) == ("response", "a.u.")


def test_legacy_config_without_readout_pair_keeps_explicit_task_metadata(prediction_case):
    case = prediction_case
    assert "quantity" not in case.config and "unit" not in case.config
    checkpoint = freeze(case)
    task = prepare_evaluation(
        **live(case, checkpoint), experiment_inputs=case.experiment_inputs,
        quantity="legacy_response", unit="legacy_unit",
    )
    prediction = case.sealer(case, checkpoint, manifest=task)
    payload = json.loads(case.label_path.read_text())
    payload.update({"quantity": "legacy_response", "unit": "legacy_unit"})
    case.label_path.write_text(json.dumps(payload))
    result = case.scorer(prediction, **live(case, checkpoint), labels=labels(case))
    assert (result["quantity"], result["unit"]) == ("legacy_response", "legacy_unit")


@pytest.mark.parametrize("field,replacement", [("quantity", "different_response"), ("unit", "mg/L")])
@pytest.mark.parametrize("stage", ["prepare", "seal", "validate", "score"])
def test_frozen_readout_cannot_drift_at_any_task_seam(prediction_case, field, replacement, stage):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "a.u."})
    checkpoint = freeze(case)
    task = prepare(case, checkpoint)
    prediction = case.sealer(case, checkpoint, manifest=task)
    payload = case.validator(prediction, **live(case, checkpoint))
    task[field] = replacement
    task = rehashed_manifest(task)
    with pytest.raises(ValueError, match="quantity|unit"):
        if stage == "prepare":
            prepare_evaluation(
                **live(case, checkpoint), experiment_inputs=case.experiment_inputs,
                quantity=task["quantity"], unit=task["unit"],
            )
        elif stage == "seal":
            case.sealer(case, checkpoint, name="readout-drift.json", manifest=task)
        else:
            payload["task_manifest"] = task
            manual = forged_artifact(case.root / "artifacts/readout-drift.json", payload)
            if stage == "validate":
                case.validator(manual, **live(case, checkpoint))
            else:
                outcomes = json.loads(case.label_path.read_text())
                outcomes[field] = replacement
                case.label_path.write_text(json.dumps(outcomes))
                case.scorer(manual, **live(case, checkpoint), labels=labels(case))


def test_frozen_dry_weight_quantity_cannot_be_relabelled_per_litre(prediction_case):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "mg/gDW"})
    checkpoint = freeze(case)
    with pytest.raises(ValueError, match="quantity|unit"):
        prepare_evaluation(
            **live(case, checkpoint), experiment_inputs=case.experiment_inputs,
            quantity="response", unit="mg/L",
        )


def test_readout_mismatch_is_rejected_before_task_bytes_are_opened(prediction_case):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "a.u."})
    checkpoint = freeze(case)
    (case.root / "tasks/conditions.json").write_text("POISON: invalid task bytes")
    with pytest.raises(ValueError, match="quantity|unit"):
        prepare_evaluation(
            **live(case, checkpoint), experiment_inputs=case.experiment_inputs,
            quantity="response", unit="different_unit",
        )


@pytest.mark.parametrize("field", ["quantity", "unit"])
def test_labels_must_match_the_frozen_readout_not_only_the_id_set(prediction_case, field):
    case = prediction_case
    case.config.update({"quantity": "response", "unit": "a.u."})
    checkpoint = freeze(case)
    prediction = case.sealer(case, checkpoint)
    payload = json.loads(case.label_path.read_text())
    payload[field] = "different"
    case.label_path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="quantity|unit"):
        case.scorer(prediction, **live(case, checkpoint), labels=labels(case))
