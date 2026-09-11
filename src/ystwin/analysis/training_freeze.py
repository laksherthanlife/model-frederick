from __future__ import annotations

import hashlib
import json
import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path


@dataclass(frozen=True)
class InputSpec:
    path: str
    sha256: str
    role: str
    source: Mapping[str, object]


@dataclass(frozen=True)
class TrainingContract:
    root: Path
    inputs: Sequence[InputSpec]
    code: Sequence[InputSpec]


@dataclass(frozen=True)
class ArtifactRef:
    path: Path
    sha256: str


_NATIVE_ROLES = frozenset({"fit", "selection", "prior", "provenance"})
_TASK_ROLES = frozenset({"task_structure", "task_conditions", "evaluation_code"})
_SCOPES = {
    **dict.fromkeys(_NATIVE_ROLES, "native_physiology"),
    "training_code": "native_training_code",
    "task_structure": "evaluation_task",
    "task_conditions": "evaluation_task",
    "evaluation_code": "evaluation_code",
    "outcome_labels": "heldout_outcomes",
}
_SOURCE_FIELDS = frozenset({"authority", "source_id", "reference", "scope"})
_INTERVAL_FAILURES = ("failed", "infeasible", "unsupported")


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be a nonempty, unpadded string")
    return value


def _sha256(value: object) -> str:
    if (not isinstance(value, str) or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)):
        raise ValueError("SHA-256 must be 64 lowercase hexadecimal characters")
    return value


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be numeric and finite")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric and finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be numeric and finite")
    return result


def _json_value(value: object) -> object:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        return _number(value, "JSON number")
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("JSON mapping keys must be strings")
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    raise ValueError("payload must contain only finite JSON values")


def _json_bytes(value: object) -> bytes:
    return json.dumps(
        _json_value(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _digest(value: object) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _manifest(kind: str, **fields: object) -> dict:
    result = _json_value({"schema_version": 1, "kind": kind, **fields})
    result["sha256"] = _digest(result)
    return result


def _keys(value: object, expected: set[str], name: str) -> None:
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ValueError(f"{name} has invalid fields")


def _root(contract: TrainingContract) -> Path:
    root = Path(contract.root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("contract root must be a directory")
    return root


def _relative(root: Path, path: str | Path) -> Path:
    if not isinstance(path, (str, Path)):
        raise ValueError("allowlist path must name an exact file")
    candidate = Path(path)
    if ".." in candidate.parts:
        raise ValueError("path escape is not allowed")
    if candidate.is_absolute():
        if not candidate.is_relative_to(root):
            raise ValueError("path escape from contract root")
        candidate = candidate.relative_to(root)
    if not candidate.parts:
        raise ValueError("allowlist path must name an exact file")
    return candidate


def _resolve_file(root: Path, path: str | Path) -> Path:
    resolved = (root / _relative(root, path)).resolve(strict=True)
    if not resolved.is_relative_to(root):
        raise ValueError("path escape through a symlink")
    if not resolved.is_file():
        raise ValueError("allowlist path must name an exact regular file")
    return resolved


def _source(spec: InputSpec, allowed_roles: frozenset[str]) -> dict:
    if spec.role not in allowed_roles:
        raise ValueError("input role is not allowed in this phase")
    source = _json_value(spec.source)
    if not isinstance(source, dict) or not _SOURCE_FIELDS <= source.keys():
        raise ValueError("source requires authority, source_id, reference, and scope")
    for field in _SOURCE_FIELDS:
        _text(source[field], f"source {field}")
    if source["scope"] != _SCOPES[spec.role]:
        raise ValueError("source scope does not authorize this input role")
    return source


def _verified_bytes(root: Path, entry: Mapping) -> bytes:
    resolved = _resolve_file(root, entry["path"])
    if resolved.relative_to(root).as_posix() != entry["resolved_path"]:
        raise ValueError("input manifest resolved path changed")
    data = resolved.read_bytes()
    if hashlib.sha256(data).hexdigest() != entry["sha256"]:
        raise ValueError(f"SHA-256 mismatch for input {entry['path']}")
    return data


def _verify_entries(root: Path, entries: Sequence[Mapping]) -> None:
    for entry in entries:
        resolved = _resolve_file(root, entry["path"])
        if resolved.relative_to(root).as_posix() != entry["resolved_path"]:
            raise ValueError("input manifest resolved path changed")
        with resolved.open("rb") as handle:
            actual = hashlib.file_digest(handle, "sha256").hexdigest()
        if actual != entry["sha256"]:
            raise ValueError(f"SHA-256 mismatch for input {entry['path']}")


def _entries(
    root: Path, specs: Sequence[InputSpec], allowed_roles: frozenset[str], *, verify: bool = True
) -> list[dict]:
    if not specs or isinstance(specs, (str, bytes)):
        raise ValueError("an explicit nonempty file allowlist is required")
    entries = []
    identities = set()
    for spec in specs:
        if not isinstance(spec, InputSpec):
            raise ValueError("allowlist entries must be InputSpec records")
        source = _source(spec, allowed_roles)
        expected_sha256 = _sha256(spec.sha256)
        relative = _relative(root, spec.path)
        resolved = _resolve_file(root, relative)
        stat = resolved.stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity in identities:
            raise ValueError("duplicate physical file in allowlist")
        identities.add(identity)
        entries.append({
            "path": relative.as_posix(),
            "resolved_path": resolved.relative_to(root).as_posix(),
            "sha256": expected_sha256,
            "size_bytes": stat.st_size,
            "role": spec.role,
            "source": source,
        })
    entries.sort(key=lambda entry: entry["path"])
    if verify:
        _verify_entries(root, entries)
    return entries


def build_training_manifest(contract: TrainingContract) -> dict:
    root = _root(contract)
    for specs, roles in ((contract.inputs, _NATIVE_ROLES),
                         (contract.code, frozenset({"training_code"}))):
        if not specs:
            raise ValueError("native inputs and explicit code allowlists must be nonempty")
        for spec in specs:
            if not isinstance(spec, InputSpec):
                raise ValueError("allowlist entries must be InputSpec records")
            _source(spec, roles)
    entries = _entries(root, (*contract.inputs, *contract.code),
                       _NATIVE_ROLES | {"training_code"})
    return _manifest(
        "training_input_manifest", root=str(root), dependency_mode="explicit",
        inputs=[entry for entry in entries if entry["role"] != "training_code"],
        code=[entry for entry in entries if entry["role"] == "training_code"],
    )


def validate_training_manifest(manifest: Mapping, contract: TrainingContract) -> dict:
    expected = build_training_manifest(contract)
    if _json_bytes(manifest) != _json_bytes(expected):
        raise ValueError("training manifest differs from the explicit source/code allowlists")
    return expected


def _declared_readout(config: Mapping) -> tuple[str, str] | None:
    fields = {"quantity", "unit"} & config.keys()
    if not fields:
        return None
    if fields != {"quantity", "unit"}:
        raise ValueError("config quantity and unit must be declared together")
    return _text(config["quantity"], "config quantity"), _text(config["unit"], "config unit")


def _model(parameters: Mapping[str, Real], config: Mapping) -> dict:
    if not isinstance(parameters, Mapping) or not parameters:
        raise ValueError("parameters must be a nonempty native parameter mapping")
    values = {_text(key, "parameter name"): _number(value, "parameter")
              for key, value in parameters.items()}
    if not isinstance(config, Mapping):
        raise ValueError("config must be a JSON mapping")
    config = _json_value(config)
    _declared_readout(config)
    return {"parameters": values, "config": config}


def _write_once(path: str | Path, root: Path, payload: Mapping) -> ArtifactRef:
    relative = _relative(root, path)
    parent = (root / relative.parent).resolve(strict=True)
    if not parent.is_relative_to(root):
        raise ValueError("artifact path escape through a symlink")
    destination = parent / relative.name
    data = _json_bytes(payload) + b"\n"
    with destination.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    return ArtifactRef(path=destination, sha256=hashlib.sha256(data).hexdigest())


def _json_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON field")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise ValueError("JSON numbers must be finite")


def _parse_json(data: bytes) -> object:
    return json.loads(data, object_pairs_hook=_json_object, parse_constant=_invalid_constant)


def _read_artifact(ref: ArtifactRef, root: Path, kind: str, fields: set[str]) -> dict:
    expected = _sha256(ref.sha256)
    data = _resolve_file(root, ref.path).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"SHA-256 mismatch for {kind} artifact")
    payload = _parse_json(data)
    _keys(payload, {"schema_version", "kind", *fields}, kind)
    if payload["schema_version"] != 1 or payload["kind"] != kind:
        raise ValueError("artifact kind or schema version does not match")
    return payload


def freeze_checkpoint(
    path: str | Path, *, contract: TrainingContract, training_manifest: Mapping,
    parameters: Mapping[str, Real], config: Mapping,
) -> ArtifactRef:
    manifest = validate_training_manifest(training_manifest, contract)
    model = _model(parameters, config)
    return _write_once(path, _root(contract), {
        "schema_version": 1, "kind": "native_checkpoint",
        "training_manifest": manifest, "model": model,
    })


def validate_checkpoint(
    checkpoint: ArtifactRef, *, contract: TrainingContract,
    parameters: Mapping[str, Real] | None = None, config: Mapping | None = None,
) -> dict:
    payload = _read_artifact(
        checkpoint, _root(contract), "native_checkpoint", {"training_manifest", "model"}
    )
    validate_training_manifest(payload["training_manifest"], contract)
    _keys(payload["model"], {"parameters", "config"}, "native model")
    model = _model(**payload["model"])
    if _json_bytes(model) != _json_bytes(payload["model"]):
        raise ValueError("checkpoint model is not a canonical native parameter payload")
    live_model = _model(
        model["parameters"] if parameters is None else parameters,
        model["config"] if config is None else config,
    )
    if _json_bytes(live_model) != _json_bytes(model):
        raise ValueError("live model parameters or config differ from frozen checkpoint")
    return payload


def _unique_rows(rows: Sequence[Mapping], fields: set[str], name: str) -> dict[str, dict]:
    if not isinstance(rows, (list, tuple)) or not rows:
        raise ValueError(f"{name} must contain nonempty explicit experiment rows")
    result = {}
    for row in rows:
        _keys(row, fields, name)
        identifier = _text(row["experiment_id"], "experiment_id")
        if identifier in result:
            raise ValueError(f"duplicate experiment_id in {name}: {identifier}")
        result[identifier] = _json_value(row)
    return result


def _same_ids(actual: set[str], expected: set[str], name: str) -> None:
    missing, extra = sorted(expected - actual), sorted(actual - expected)
    if missing or extra:
        raise ValueError(f"{name} experiment IDs differ: missing={missing}, extra={extra}")


def _code_identity(entry: Mapping) -> tuple:
    return (
        entry["path"], entry["resolved_path"], entry["sha256"],
        entry["source"]["authority"], entry["source"]["source_id"], entry["source"]["reference"],
    )


def _require_frozen_evaluation_code(entries: Sequence[Mapping], code: Sequence[Mapping]) -> None:
    approved = {_code_identity(entry) for entry in code}
    for entry in entries:
        if entry["role"] == "evaluation_code" and _code_identity(entry) not in approved:
            raise ValueError(
                f"evaluation_code {entry['path']} identity is not committed in checkpoint training code"
            )


def _task_manifest(
    root: Path, checkpoint_sha256: str, task_inputs: Sequence[InputSpec],
    experiment_inputs: Sequence[Mapping], quantity: str, unit: str, *, frozen: Mapping,
) -> dict:
    quantity, unit = _text(quantity, "quantity"), _text(unit, "unit")
    readout = _declared_readout(frozen["model"]["config"])
    if readout is not None and readout != (quantity, unit):
        raise ValueError("task quantity or unit differs from frozen checkpoint config")
    entries = _entries(root, task_inputs, _TASK_ROLES, verify=False)
    roles = {entry["role"] for entry in entries}
    if "evaluation_code" not in roles:
        raise ValueError("task inputs require an explicit evaluation_code dependency list")
    _require_frozen_evaluation_code(entries, frozen["training_manifest"]["code"])
    paths = {entry["path"] for entry in entries if entry["role"] != "evaluation_code"}
    if not paths:
        raise ValueError("task inputs require explicit chemistry or condition files")
    rows = _unique_rows(experiment_inputs, {"experiment_id", "input_paths"}, "experiment inputs")
    used = set()
    for row in rows.values():
        values = row["input_paths"]
        if not isinstance(values, list) or not values:
            raise ValueError("each experiment requires explicit input_paths")
        values = [_relative(root, path).as_posix() for path in values]
        if len(set(values)) != len(values) or not set(values) <= paths:
            raise ValueError("duplicate or non-allowlisted experiment input path")
        row["input_paths"] = sorted(values)
        used.update(values)
    if used != paths:
        raise ValueError("unused task input paths are not bound to an experiment")
    _verify_entries(root, entries)
    return _manifest(
        "evaluation_task_manifest", root=str(root), checkpoint_sha256=checkpoint_sha256,
        inputs=entries, experiment_inputs=[rows[key] for key in sorted(rows)],
        quantity=quantity, unit=unit,
    )


def prepare_evaluation(
    *, checkpoint: ArtifactRef, contract: TrainingContract, task_inputs: Sequence[InputSpec],
    experiment_inputs: Sequence[Mapping], quantity: str, unit: str,
    parameters: Mapping[str, Real], config: Mapping,
) -> dict:
    saved = validate_checkpoint(checkpoint, contract=contract, **_model(parameters, config))
    return _task_manifest(
        _root(contract), checkpoint.sha256, task_inputs, experiment_inputs, quantity, unit,
        frozen=saved,
    )


def _validate_task_manifest(
    manifest: Mapping, root: Path, checkpoint: ArtifactRef, task_inputs: Sequence[InputSpec],
    *, frozen: Mapping,
) -> dict:
    _keys(manifest, {
        "schema_version", "kind", "root", "checkpoint_sha256", "inputs",
        "experiment_inputs", "quantity", "unit", "sha256",
    }, "task manifest")
    expected = _task_manifest(
        root, checkpoint.sha256, task_inputs,
        manifest["experiment_inputs"], manifest["quantity"], manifest["unit"], frozen=frozen,
    )
    if _json_bytes(manifest) != _json_bytes(expected):
        raise ValueError("task manifest differs from checkpoint or explicit input allowlist")
    return expected


def _prediction_rows(predictions: Sequence[Mapping], expected_ids: set[str]) -> list[dict]:
    rows = _unique_rows(
        predictions, {"experiment_id", "status", "value", "reason"}, "predictions"
    )
    _same_ids(set(rows), expected_ids, "predictions")
    for row in rows.values():
        if row["status"] == "ok":
            row["value"] = _number(row["value"], "prediction")
            if row["reason"] is not None:
                raise ValueError("successful prediction reason must be null")
        elif row["status"] in ("failed", "infeasible"):
            if row["value"] is not None:
                raise ValueError("failed or infeasible prediction value must be null, not zero")
            _text(row["reason"], "prediction failure reason")
        else:
            raise ValueError("prediction status must be ok, failed, or infeasible")
    return [rows[key] for key in sorted(rows)]


def _interval_rows(predictions: Sequence[Mapping], expected_ids: set[str]) -> list[dict]:
    rows = _unique_rows(
        predictions, {"experiment_id", "status", "lower", "upper", "reason"}, "interval predictions"
    )
    _same_ids(set(rows), expected_ids, "interval predictions")
    for row in rows.values():
        if row["status"] == "ok":
            row["lower"] = _number(row["lower"], "interval lower bound")
            row["upper"] = _number(row["upper"], "interval upper bound")
            if row["lower"] > row["upper"]:
                raise ValueError("interval lower bound must not exceed upper bound")
            _number(row["upper"] - row["lower"], "interval width")
            if row["reason"] is not None:
                raise ValueError("successful interval prediction reason must be null")
        elif row["status"] in _INTERVAL_FAILURES:
            if row["lower"] is not None or row["upper"] is not None:
                raise ValueError("unsuccessful interval bounds must both be null, not zero")
            _text(row["reason"], "interval failure reason")
        else:
            raise ValueError("interval status must be ok, failed, infeasible, or unsupported")
    return [rows[key] for key in sorted(rows)]


def _prediction_metadata(config: Mapping, interval: bool) -> dict:
    if not interval:
        if config.get("prediction_kind", "point") != "point":
            raise ValueError("scalar predictions require prediction_kind point, not an envelope")
        return {}
    if config.get("prediction_kind") != "feasible_envelope":
        raise ValueError("interval prediction_kind must be feasible_envelope in frozen config")
    return {
        "prediction_kind": "feasible_envelope",
        "prediction_scope": _text(config.get("prediction_scope"), "prediction_scope"),
    }


def _seal_prediction_artifact(
    path: str | Path, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_manifest: Mapping, task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real],
    config: Mapping, predictions: Sequence[Mapping], interval: bool,
) -> ArtifactRef:
    saved = validate_checkpoint(checkpoint, contract=contract, **_model(parameters, config))
    metadata = _prediction_metadata(saved["model"]["config"], interval)
    root = _root(contract)
    task = _validate_task_manifest(task_manifest, root, checkpoint, task_inputs, frozen=saved)
    rows = (_interval_rows if interval else _prediction_rows)(
        predictions, {row["experiment_id"] for row in task["experiment_inputs"]}
    )
    return _write_once(path, root, {
        "schema_version": 1, "kind": "sealed_interval_predictions" if interval else "sealed_predictions",
        "checkpoint_sha256": checkpoint.sha256, "model_sha256": _digest(saved["model"]),
        "task_manifest": task, "predictions": rows, **metadata,
    })


def _validate_prediction_artifact(
    predictions: ArtifactRef, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real], config: Mapping, interval: bool,
) -> dict:
    saved = validate_checkpoint(checkpoint, contract=contract, **_model(parameters, config))
    metadata = _prediction_metadata(saved["model"]["config"], interval)
    root = _root(contract)
    kind = "sealed_interval_predictions" if interval else "sealed_predictions"
    payload = _read_artifact(predictions, root, kind, {
        "checkpoint_sha256", "model_sha256", "task_manifest", "predictions", *metadata,
    })
    if payload["checkpoint_sha256"] != checkpoint.sha256:
        raise ValueError("predictions are bound to a different checkpoint")
    if payload["model_sha256"] != _digest(saved["model"]):
        raise ValueError("prediction model differs from checkpoint")
    if any(payload[key] != value for key, value in metadata.items()):
        raise ValueError("prediction kind or scope differs from frozen config")
    task = _validate_task_manifest(
        payload["task_manifest"], root, checkpoint, task_inputs, frozen=saved
    )
    rows = (_interval_rows if interval else _prediction_rows)(
        payload["predictions"], {row["experiment_id"] for row in task["experiment_inputs"]}
    )
    if _json_bytes(rows) != _json_bytes(payload["predictions"]):
        raise ValueError("sealed predictions are not canonical experiment rows")
    return payload


def seal_predictions(
    path: str | Path, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_manifest: Mapping, task_inputs: Sequence[InputSpec],
    parameters: Mapping[str, Real], config: Mapping, predictions: Sequence[Mapping],
) -> ArtifactRef:
    return _seal_prediction_artifact(
        path, checkpoint=checkpoint, contract=contract, task_manifest=task_manifest,
        task_inputs=task_inputs, parameters=parameters, config=config, predictions=predictions,
        interval=False,
    )


def validate_predictions(
    predictions: ArtifactRef, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real], config: Mapping,
) -> dict:
    return _validate_prediction_artifact(
        predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
        parameters=parameters, config=config, interval=False,
    )


def seal_interval_predictions(
    path: str | Path, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_manifest: Mapping, task_inputs: Sequence[InputSpec],
    parameters: Mapping[str, Real], config: Mapping, predictions: Sequence[Mapping],
) -> ArtifactRef:
    return _seal_prediction_artifact(
        path, checkpoint=checkpoint, contract=contract, task_manifest=task_manifest,
        task_inputs=task_inputs, parameters=parameters, config=config, predictions=predictions,
        interval=True,
    )


def validate_interval_predictions(
    predictions: ArtifactRef, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real], config: Mapping,
) -> dict:
    return _validate_prediction_artifact(
        predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
        parameters=parameters, config=config, interval=True,
    )


def _reject_label_overlap(labels: Sequence[Mapping], inputs: Sequence[Mapping]) -> None:
    for label in labels:
        for entry in inputs:
            same_path = bool(
                {label["path"], label["resolved_path"]}
                & {entry["path"], entry["resolved_path"]}
            )
            same_source = (
                label["source"]["authority"], label["source"]["source_id"]
            ) == (entry["source"]["authority"], entry["source"]["source_id"])
            if same_path or same_source or label["sha256"] == entry["sha256"]:
                raise ValueError("outcome labels overlap an input; labels cannot be reused as inputs")


def _scoring_labels(
    sealed: Mapping, contract: TrainingContract, labels: Sequence[InputSpec]
) -> tuple[dict, dict]:
    root = _root(contract)
    task = sealed["task_manifest"]
    label_entries = _entries(root, labels, frozenset({"outcome_labels"}), verify=False)
    training = build_training_manifest(contract)
    _reject_label_overlap(label_entries, [*training["inputs"], *training["code"], *task["inputs"]])
    observed = {}
    for entry in label_entries:
        payload = _parse_json(_verified_bytes(root, entry))
        _keys(payload, {"quantity", "unit", "observations"}, "outcome labels")
        if payload["quantity"] != task["quantity"] or payload["unit"] != task["unit"]:
            raise ValueError("label quantity or unit differs from frozen prediction task")
        rows = _unique_rows(payload["observations"], {"experiment_id", "value"}, "outcome labels")
        for identifier, row in rows.items():
            if identifier in observed:
                raise ValueError(f"duplicate experiment_id across label files: {identifier}")
            observed[identifier] = {
                "value": _number(row["value"], "outcome label"), "source_path": entry["path"],
            }
    _same_ids(set(observed), {row["experiment_id"] for row in sealed["predictions"]}, "labels")
    return _manifest("outcome_labels_manifest", root=str(root), inputs=label_entries), observed


def _mean(values: Sequence[float]) -> float | None:
    return math.fsum(value / len(values) for value in values) if values else None


def score_predictions(
    predictions: ArtifactRef, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real], config: Mapping,
    labels: Sequence[InputSpec],
) -> dict:
    sealed = validate_predictions(
        predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
        parameters=parameters, config=config,
    )
    label_manifest, observed = _scoring_labels(sealed, contract, labels)
    task = sealed["task_manifest"]
    rows, errors = [], []
    for prediction in sealed["predictions"]:
        label = observed[prediction["experiment_id"]]
        residual = None
        if prediction["status"] == "ok":
            residual = _number(prediction["value"] - label["value"], "scoring residual")
            errors.append(abs(residual))
        rows.append({
            "experiment_id": prediction["experiment_id"], "status": prediction["status"],
            "predicted": prediction["value"], "observed": label["value"], "residual": residual,
            "reason": prediction["reason"], "label_source_path": label["source_path"],
        })
    maximum = max(errors, default=0.0)
    rmse = (maximum * math.sqrt(math.fsum((error / maximum) ** 2 for error in errors) / len(errors))
            if maximum else 0.0)
    return _manifest(
        "prediction_score", checkpoint_sha256=checkpoint.sha256, predictions_sha256=predictions.sha256,
        label_manifest=label_manifest, quantity=task["quantity"], unit=task["unit"], rows=rows,
        n_total=len(rows), n_scored=len(errors), n_failed=len(rows) - len(errors),
        complete=len(errors) == len(rows), mae_on_successes=_mean(errors),
        rmse_on_successes=rmse if errors else None,
    )


def score_interval_predictions(
    predictions: ArtifactRef, *, checkpoint: ArtifactRef, contract: TrainingContract,
    task_inputs: Sequence[InputSpec], parameters: Mapping[str, Real], config: Mapping,
    labels: Sequence[InputSpec],
) -> dict:
    sealed = validate_interval_predictions(
        predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
        parameters=parameters, config=config,
    )
    label_manifest, observed = _scoring_labels(sealed, contract, labels)
    task = sealed["task_manifest"]
    rows, widths, violations = [], [], []
    status_counts = dict.fromkeys(("ok", *_INTERVAL_FAILURES), 0)
    for prediction in sealed["predictions"]:
        label = observed[prediction["experiment_id"]]
        status_counts[prediction["status"]] += 1
        width = contained = violation = side = None
        if prediction["status"] == "ok":
            lower, upper, value = prediction["lower"], prediction["upper"], label["value"]
            width = upper - lower
            if value < lower:
                violation, side = _number(lower - value, "interval violation"), "below"
            elif value > upper:
                violation, side = _number(value - upper, "interval violation"), "above"
            else:
                violation = 0.0
            contained = side is None
            widths.append(width)
            violations.append(violation)
        rows.append({
            **prediction, "observed": label["value"], "label_source_path": label["source_path"],
            "contained": contained, "width": width, "violation": violation, "violation_side": side,
        })
    total, scorable = len(rows), len(widths)
    contained = sum(row["contained"] is True for row in rows)
    return _manifest(
        "interval_prediction_score", checkpoint_sha256=checkpoint.sha256,
        predictions_sha256=predictions.sha256, label_manifest=label_manifest,
        quantity=task["quantity"], unit=task["unit"], prediction_kind=sealed["prediction_kind"],
        prediction_scope=sealed["prediction_scope"],
        assessment="feasible_envelope_consistency_not_point_forecast", rows=rows,
        complete=scorable == total, status_counts=status_counts,
        n_total=total, n_successful=status_counts["ok"], n_scorable=scorable,
        successful_fraction_total=status_counts["ok"] / total, scorable_fraction_total=scorable / total,
        n_contained=contained, n_violated=scorable - contained,
        containment_on_scorable=contained / scorable if scorable else None,
        containment_fraction_total=contained / total, mean_width_on_scorable=_mean(widths),
        mean_violation_on_scorable=_mean(violations),
        max_violation_on_scorable=max(violations, default=None),
    )
