from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import pytest

from scripts import evaluate_frozen_transfer as runner
from ystwin.analysis.training_freeze import InputSpec, TrainingContract, build_training_manifest, freeze_checkpoint


def source(scope):
    return {"authority": "fixture source", "source_id": scope, "reference": "fixture", "scope": scope}


def bundle(root, monkeypatch):
    script = root / "coordinator.py"
    script.write_text("fixture coordinator")
    monkeypatch.setattr(runner, "__file__", str(script))
    chemical = root / "chemistry.json"
    chemical.write_text(json.dumps({"schema_version": 1, "source": source("evaluation_task"), "task": {
        "metabolites": [{"id": "opaque_output", "formula": "C6H12O6", "charge": 0, "compartment": "c"}],
        "reactions": [{"id": "opaque_conversion", "stoichiometry": [["native_feed", -1], ["opaque_output", 1]],
                       "source": "independent chemistry"}], "output_metabolite_id": "opaque_output"}}))
    conditions = {"schema_version": 1, "source": source("evaluation_task"), "experiments": [{
        "experiment_id": "a", "chemistry_path": "chemistry.json",
        "chemistry_sha256": hashlib.sha256(chemical.read_bytes()).hexdigest(), "added_nacl_molar": 0,
        "condition": {"experiment_id": "a", "growth_rate_per_h": 0.1, "carbon_input": "glucose",
                      "readout_basis": "balanced_growth_intracellular_content",
                      "native_context_assumption": "explicit native transfer assumption"}}]}
    path = root / "conditions.json"
    path.write_text(json.dumps(conditions))
    code_specs = (InputSpec("coordinator.py", hashlib.sha256(script.read_bytes()).hexdigest(),
                            "training_code", source("native_training_code")),)
    return path, conditions, code_specs


def test_condition_reader_passes_only_explicit_chemical_and_condition_fields(tmp_path, monkeypatch):
    path, _, code_specs = bundle(tmp_path, monkeypatch)
    specs, experiments, jobs = runner._prepare_bundle(tmp_path, path, code_specs)
    assert {item.role for item in specs} == {"task_structure", "task_conditions", "evaluation_code"}
    assert experiments[0]["experiment_id"] == jobs[0]["condition"]["experiment_id"] == "a"
    assert set(jobs[0]) == {"condition", "task", "added_nacl_molar"}


def test_outcome_fields_cannot_be_hidden_in_condition_inputs(tmp_path, monkeypatch):
    path, payload, code_specs = bundle(tmp_path, monkeypatch)
    payload["experiments"][0]["condition"]["observed_content"] = 5
    path.write_text(json.dumps(payload))
    with pytest.raises(TypeError):
        runner._prepare_bundle(tmp_path, path, code_specs)


def test_invalid_checkpoint_is_rejected_before_task_inputs_are_opened(tmp_path, monkeypatch):
    native = tmp_path / "native.txt"
    code = tmp_path / "code.py"
    native.write_text("native input")
    code.write_text("native code")
    contract = TrainingContract(tmp_path,
        (InputSpec("native.txt", hashlib.sha256(native.read_bytes()).hexdigest(), "fit", source("native_physiology")),),
        (InputSpec("code.py", hashlib.sha256(code.read_bytes()).hexdigest(), "training_code", source("native_training_code")),))
    checkpoint = freeze_checkpoint(tmp_path / "checkpoint.json", contract=contract,
                                   training_manifest=build_training_manifest(contract), parameters={"native_rate": 1},
                                   config={"prediction_kind": "feasible_envelope", "prediction_scope": "native capacity"})
    contract_file = tmp_path / "contract.json"
    contract_file.write_text(json.dumps({"root": str(tmp_path), "inputs": [asdict(item) for item in contract.inputs],
                                         "code": [asdict(item) for item in contract.code]}))
    monkeypatch.setattr(runner, "_prepare_bundle", lambda *args: pytest.fail("task inputs opened before checkpoint validation"))
    with pytest.raises(ValueError, match="SHA-256"):
        runner.main(["predict", "--contract", str(contract_file), "--checkpoint", str(checkpoint.path),
                     "--checkpoint-sha256", "0" * 64, "--output-dir", str(tmp_path / "prediction"),
                     "--conditions", str(tmp_path / "unopened.json")])


def test_duplicate_json_fields_are_not_silently_accepted(tmp_path):
    path = Path(tmp_path / "duplicate.json")
    path.write_text('{"experiment_id":"a","experiment_id":"b"}')
    with pytest.raises(ValueError, match="duplicate"):
        runner._read(path)


@pytest.mark.parametrize("value", [float("inf"), -float("inf"), float("nan"), object()])
def test_invalid_diagnostics_never_leave_a_partial_json_file(tmp_path, value):
    destination = tmp_path / "prediction_details.json"
    with pytest.raises((ValueError, TypeError)):
        runner._write(destination, {"sealed_predictions_modified": False, "bounds": [0.0, value]})
    assert not destination.exists()


def test_diagnostic_writer_remains_exclusive_and_readable(tmp_path):
    destination = tmp_path / "prediction_details.json"
    record = {"status": "infeasible", "prediction": None}
    runner._write(destination, record)
    assert runner._read(destination) == record
    previous = destination.read_bytes()
    with pytest.raises(FileExistsError):
        runner._write(destination, {"status": "different"})
    assert destination.read_bytes() == previous
