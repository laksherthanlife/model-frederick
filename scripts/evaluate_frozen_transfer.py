from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
import multiprocessing
from pathlib import Path
import sys

from ystwin.analysis.blind_transfer import BalancedCondition, HostFluxRoles, predict_balanced_envelope
from ystwin.analysis.native_physiology import NativeExchangeModel
from ystwin.analysis.training_freeze import (
    ArtifactRef, InputSpec, TrainingContract, prepare_evaluation, score_interval_predictions,
    seal_interval_predictions, validate_checkpoint, validate_interval_predictions,
)
from ystwin.fba.chemical_task import ChemicalTask
from ystwin.fba.native_obligations import NativeCarbonObligation, NativeObligationRoles
from ystwin.fba.solver import load_model

_WORKER = None


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON field")
        result[key] = value
    return result


def _read(path):
    return json.loads(path.read_text(), object_pairs_hook=_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError("nonfinite JSON number")))


def _write(path, payload):
    serialized = json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n"
    with path.open("x") as handle:
        handle.write(serialized)


def _contract(path):
    payload = _read(path)
    return TrainingContract(Path(payload["root"]), tuple(InputSpec(**item) for item in payload["inputs"]),
                            tuple(InputSpec(**item) for item in payload["code"]))


def _assert_frozen_dependencies(contract):
    root = Path(contract.root).resolve()
    allowed = {(root / item.path).resolve() for item in contract.code}
    actual = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name = getattr(module, "__name__", "")
        filename = getattr(module, "__file__", None)
        if filename is not None and (name == "ystwin" or name.startswith("ystwin.")):
            actual.add(Path(filename).resolve())
    if not actual <= allowed:
        raise ValueError("unfrozen evaluation dependency: the coordinator and all loaded model modules must be committed before task introduction")


def _inside(root, path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    resolved = candidate.resolve(strict=True)
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise ValueError("evaluation inputs must be explicit files inside the contract root")
    return resolved


def _spec(root, path, role, source, expected=None):
    path = _inside(root, path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if expected is not None and expected != digest:
        raise ValueError("declared evaluation input checksum differs from its actual bytes")
    return InputSpec(str(path.relative_to(root)), digest, role, source)


def _prepare_bundle(root, path, code_specs):
    path = _inside(root, path)
    bundle = _read(path)
    if set(bundle) != {"schema_version", "source", "experiments"} or bundle["schema_version"] != 1:
        raise ValueError("invalid condition-bundle schema")
    specs = [_spec(root, path, "task_conditions", bundle["source"])]
    experiment_inputs, jobs, chemistry = [], [], {}
    row_fields = {"experiment_id", "chemistry_path", "chemistry_sha256", "condition", "added_nacl_molar"}
    ids = set()
    for row in bundle["experiments"]:
        if set(row) != row_fields or row["experiment_id"] in ids:
            raise ValueError("invalid or duplicate evaluation experiment")
        ids.add(row["experiment_id"])
        condition = BalancedCondition(**row["condition"])
        if condition.experiment_id != row["experiment_id"]:
            raise ValueError("condition identity differs from the declared experiment")
        chemical_path = _inside(root, row["chemistry_path"])
        if chemical_path not in chemistry:
            payload = _read(chemical_path)
            if set(payload) != {"schema_version", "source", "task"} or payload["schema_version"] != 1:
                raise ValueError("invalid chemical-input schema")
            chemical_spec = _spec(root, chemical_path, "task_structure", payload["source"], row["chemistry_sha256"])
            task = ChemicalTask.from_dict(payload["task"])
            specs.append(chemical_spec)
            chemistry[chemical_path] = (task.to_dict(), chemical_spec)
        task, chemical_spec = chemistry[chemical_path]
        if row["chemistry_sha256"] != chemical_spec.sha256:
            raise ValueError("shared chemical input has conflicting declared digests")
        experiment_inputs.append({"experiment_id": row["experiment_id"],
                                  "input_paths": [str(path.relative_to(root)), chemical_spec.path]})
        jobs.append({"condition": asdict(condition), "task": task, "added_nacl_molar": row["added_nacl_molar"]})
    if not code_specs:
        raise ValueError("evaluation requires the pre-frozen executable dependency roster")
    specs.extend(InputSpec(item.path, item.sha256, "evaluation_code",
                           {**dict(item.source), "scope": "evaluation_code"}) for item in code_specs)
    return tuple(specs), experiment_inputs, jobs


def _initialize_worker(root, config):
    global _WORKER
    host, _ = load_model(Path(root) / config["host_model_path"])
    role_values = dict(config["host_roles"])
    role_values["native"] = NativeObligationRoles(**role_values["native"])
    _WORKER = (host, HostFluxRoles(**role_values), NativeExchangeModel.from_dict(config["native_exchange_model"]), config)


def _predict(job):
    host, roles, physiology, config = _WORKER
    condition = BalancedCondition(**job["condition"])
    dose = float(job["added_nacl_molar"])
    summary = config["native_obligation_scenarios"].get(format(dose, ".17g"))
    if summary is None:
        return ({"experiment_id": condition.experiment_id, "status": "unsupported", "lower": None,
                 "upper": None, "reason": "osmotic perturbation is outside the frozen native scenario support"}, {})
    obligations = tuple(NativeCarbonObligation(summary[key], summary["source"], summary["source_medium"],
                                              {"native_training_scope": "glucose_only",
                                               "summary": summary["scope"], "nacl_molar": dose})
                        for key in ("minimum", "maximum"))
    return predict_balanced_envelope(host, task=ChemicalTask.from_dict(job["task"]), condition=condition,
                                     physiology=physiology, obligations=obligations, roles=roles)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Predict from a frozen native model, seal outputs, then separately score held-out labels")
    parser.add_argument("mode", choices=("predict", "score"))
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--checkpoint-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--conditions", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--predictions-sha256")
    parser.add_argument("--label-manifest", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args(argv)
    contract = _contract(args.contract)
    root = Path(contract.root).resolve()
    checkpoint = ArtifactRef(args.checkpoint, args.checkpoint_sha256)
    frozen = validate_checkpoint(checkpoint, contract=contract)
    _assert_frozen_dependencies(contract)
    parameters, config = frozen["model"]["parameters"], frozen["model"]["config"]
    output = args.output_dir.resolve()
    if not output.is_relative_to(root) or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise ValueError("use a new or empty output directory inside the contract root")
    if args.mode == "predict":
        if args.conditions is None or args.workers < 1:
            raise ValueError("prediction requires explicit conditions and a positive worker count")
        task_inputs, experiment_inputs, jobs = _prepare_bundle(root, args.conditions, contract.code)
        task_manifest = prepare_evaluation(checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
                                           experiment_inputs=experiment_inputs, quantity=config["quantity"],
                                           unit=config["unit"], parameters=parameters, config=config)
        output.mkdir(parents=True, exist_ok=True)
        _write(output / "task_manifest.json", task_manifest)
        if args.workers == 1:
            _initialize_worker(root, config)
            results = [_predict(job) for job in jobs]
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs)),
                                     mp_context=multiprocessing.get_context("spawn"),
                                     initializer=_initialize_worker, initargs=(str(root), config)) as executor:
                results = list(executor.map(_predict, jobs))
        predictions = seal_interval_predictions(output / "predictions.json", checkpoint=checkpoint, contract=contract,
                                                task_manifest=task_manifest, task_inputs=task_inputs,
                                                parameters=parameters, config=config,
                                                predictions=[result[0] for result in results])
        _write(output / "prediction_details.json", {result[0]["experiment_id"]: result[1] for result in results})
        _write(output / "evaluation_receipt.json", {"checkpoint_sha256": checkpoint.sha256,
                                                    "predictions": {"path": str(predictions.path.relative_to(root)),
                                                                    "sha256": predictions.sha256},
                                                    "task_inputs": [asdict(item) for item in task_inputs]})
        print(json.dumps({"predictions": str(predictions.path.relative_to(root)), "sha256": predictions.sha256,
                          "labels_opened": False, "experiments": len(results)}, indent=2))
    else:
        if args.receipt is None or args.predictions_sha256 is None or args.label_manifest is None:
            raise ValueError("scoring requires a sealed-prediction receipt, its retained digest, and a label manifest")
        receipt = _read(_inside(root, args.receipt))
        if receipt["checkpoint_sha256"] != checkpoint.sha256 or receipt["predictions"]["sha256"] != args.predictions_sha256:
            raise ValueError("retained checkpoint/prediction digests do not match the receipt")
        task_inputs = tuple(InputSpec(**item) for item in receipt["task_inputs"])
        predictions = ArtifactRef(Path(receipt["predictions"]["path"]), args.predictions_sha256)
        validate_interval_predictions(predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
                                      parameters=parameters, config=config)
        labels = tuple(InputSpec(**item) for item in _read(_inside(root, args.label_manifest))["labels"])
        report = score_interval_predictions(predictions, checkpoint=checkpoint, contract=contract, task_inputs=task_inputs,
                                             parameters=parameters, config=config, labels=labels)
        output.mkdir(parents=True, exist_ok=True)
        _write(output / "score.json", report)
        print(json.dumps({key: report[key] for key in ("assessment", "n_total", "n_scorable",
                                                      "containment_fraction_total", "mean_width_on_scorable")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
