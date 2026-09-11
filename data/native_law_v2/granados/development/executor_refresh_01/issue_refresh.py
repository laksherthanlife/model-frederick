#!/usr/bin/env python3
"""Metadata-only executor refresh. Never fit, score, or regenerate data."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import stat
import sys
import types

import jsonschema

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
DEV = OUT.parent
INPUT = ROOT / "outputs/native_population_development/executor_refresh_01"
PROTOCOL = ROOT / "data/native_law_v2/development_protocol.json"
MODULE = ROOT / "src/ystwin/analysis/native_population_development.py"
PINS = {
    "protocol": "0f53f40699ebd0db6df62f00e0a0e6bfc5b56b374ec84a6e36846a8488e6dad6",
    "refresh_request": "0b9539cfc3c53dd88186541613c1e6cdb338ecb79d9406db592bdbb215236d91",
    "request": "3bc3f5d0441288a13c67bc3a17380002ed09916340f78102a5d39fb16c782caf",
    "configuration": "7e1aa35c8f6c248fc921f2348d362b8fc76e23bd8cfff070b19d21b2d9195cdb",
    "module": "654dd5729d204c5963f5e51576a4830b61f736296c877ffc419efa77d009e2a0",
    "phase_b_manifest": "3b8c67c79191bf4dec0a1e7e0514b1a145359c30061868bbfdd136f0ee1cd3dd",
    "previous_approval": "87536dc42cccde2733daf4cb52f54024661a7d58386acd180623adf266126533",
    "training": "fc8000c7720e554756c90e0da79f2d0ae4a7928d5923f0efa1a2ed5823d37c9e",
    "development_inputs": "41ee2a8ce6b09cb664859873fd59fd005e31e6eaa054c3703098c91e3f790613",
}
PROTOCOL_ID = "granados-sfp1-published-cohort-development-01"
REFRESH_FIELDS = {
    "document_type", "schema_version", "protocol_id", "scientific_scope", "recorded_at_utc",
    "phase_b_manifest", "previous_request", "previous_approval", "previous_executor_artifacts",
    "new_request", "new_configuration", "fixes", "reuse_existing_projections", "scientific_recipe_changed",
    "native_fits_performed_by_this_request", "development_responses_accessed", "data_reexport_requested",
    "required_custodian_admission", "strong_claim_authorized",
}
ADMISSION_FIELDS = [
    "document_type", "schema_version", "protocol_id", "protocol_sha256", "phase_b_manifest",
    "previous_approval_sha256", "request", "approval", "configuration", "custodian_identity",
    "operational_event_reference", "reuse_existing_projections", "scientific_recipe_changed",
    "data_reexport_permitted", "development_responses_released", "strong_claim_authorized",
]
VERIFIED = {}


class RefreshError(Exception):
    pass


def require(value, reason):
    if not value:
        raise RefreshError(reason)


def sha(content):
    return hashlib.sha256(content).hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def strict_json(content):
    def pairs(items):
        output = {}
        for key, value in items:
            require(key not in output, "Duplicate JSON key")
            output[key] = value
        return output

    def number(token):
        value = float(token)
        require(math.isfinite(value), "Nonfinite JSON numeric token")
        return value

    def invalid(_):
        raise RefreshError("Nonfinite JSON constant")

    return json.loads(content, object_pairs_hook=pairs, parse_float=number, parse_constant=invalid)


def regular(path):
    path = Path(path)
    require(path.is_absolute() and path.is_relative_to(ROOT), "Artifact is outside the repository")
    require(not path.is_symlink() and not any(p.is_symlink() for p in path.parents), "Symlink artifact refused")
    require(path.resolve() == path and path.is_file() and stat.S_ISREG(path.stat().st_mode), "Artifact is not a canonical regular file")
    return path


def checked(path, expected_sha, expected_bytes=None):
    path = regular(path)
    content = path.read_bytes()
    require(sha(content) == expected_sha, "Fresh artifact SHA256 mismatch: " + relative(path))
    require(expected_bytes is None or len(content) == expected_bytes, "Fresh artifact length mismatch: " + relative(path))
    VERIFIED[relative(path)] = {"path": relative(path), "sha256": expected_sha, "bytes": len(content)}
    return content


def check_ref(reference, protocol, role=None, parse=True):
    validate_type(protocol, reference, "ArtifactRef")
    require(role is None or reference["role"] == role, "Artifact role mismatch")
    content = checked(ROOT / reference["path"], reference["sha256"], reference["bytes"])
    return strict_json(content) if parse else content


def validate_type(protocol, value, name):
    schemas = protocol["export_contract"]["schemas"]
    schema = {"$schema": schemas["$schema"], "$ref": "#/$defs/" + name, "$defs": schemas["$defs"]}
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.Draft202012Validator(schema).validate(value)


def write_new(path, value):
    path = Path(path)
    require(path.is_relative_to(OUT) and path.resolve().is_relative_to(OUT), "Refresh writer scope violation")
    require(not path.is_symlink() and not any(p.is_symlink() for p in path.parents), "Refresh output symlink refused")
    content = (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()
    if path.exists():
        require(regular(path).read_bytes() == content, "Existing refresh record differs; never overwrite history")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
        path.chmod(0o600)
    return {"path": relative(path), "sha256": sha(content), "bytes": len(content)}


def fresh_executor():
    content = checked(MODULE, PINS["module"])
    module = types.ModuleType("custody_checked_executor_refresh_01")
    module.__file__ = str(MODULE)
    sys.modules[module.__name__] = module
    # Execute verified source bytes, not a possibly stale .pyc. Module loading
    # defines functions/imports libraries only; no fit/prediction function runs.
    exec(compile(content, str(MODULE), "exec"), module.__dict__)
    return module


def verify_inputs():
    protocol = strict_json(checked(PROTOCOL, PINS["protocol"]))
    context = strict_json(checked(INPUT / "executor_refresh_request.json", PINS["refresh_request"]))
    require(set(context) == REFRESH_FIELDS, "Unexpected executor refresh-request fields")
    require(context["document_type"] == "population_development_executor_refresh_request_v1" and context["schema_version"] == 1 and context["protocol_id"] == PROTOCOL_ID, "Wrong refresh request version/protocol")
    require(context["scientific_scope"] == "retrospective_source_cohort_development", "Wrong scientific scope")
    require(all(context[key] is False for key in ("scientific_recipe_changed", "native_fits_performed_by_this_request", "development_responses_accessed", "data_reexport_requested", "strong_claim_authorized")), "Refresh request widens scope or reports prior fitting/access")
    require(context["new_request"]["sha256"] == PINS["request"] and context["new_configuration"]["sha256"] == PINS["configuration"], "New request/configuration pins differ")
    require(context["phase_b_manifest"]["sha256"] == PINS["phase_b_manifest"] and context["previous_approval"]["sha256"] == PINS["previous_approval"], "Original phase-B pins differ")
    request = check_ref(context["new_request"], protocol, "operational_event_log")
    configuration = check_ref(context["new_configuration"], protocol, "execution_code")
    manifest = check_ref(context["phase_b_manifest"], protocol, "operational_event_log")
    previous_approval = check_ref(context["previous_approval"], protocol, "operational_approval")
    previous_request = check_ref(context["previous_request"], protocol, "operational_event_log")
    for value, definition in ((request, "DevelopmentReleaseRequest"), (previous_request, "DevelopmentReleaseRequest"), (previous_approval, "DevelopmentCustodyDecision")):
        validate_type(protocol, value, definition)
    require(request["protocol_sha256"] == configuration["protocol_sha256"] == manifest["protocol_sha256"] == PINS["protocol"], "Protocol binding mismatch")
    require(context["previous_approval"]["sha256"] == manifest["custody_decision"]["sha256"] and context["previous_request"]["sha256"] == manifest["release_request"]["sha256"], "Refresh rewrites original custody history")
    require(context["previous_executor_artifacts"] == manifest["frozen_executor_artifacts"], "Historical executor inventory was altered")
    for ref in context["previous_executor_artifacts"]:
        check_ref(ref, protocol, "execution_code", parse=False)
    check_ref(manifest["protocol_snapshot"], protocol, "protocol_snapshot")
    require(context["reuse_existing_projections"] == manifest["learner_allowlist"], "Projection references changed")
    require([ref["sha256"] for ref in context["reuse_existing_projections"]] == [PINS["training"], PINS["development_inputs"]], "Phase-B data digests changed")
    for ref in context["reuse_existing_projections"]:
        check_ref(ref, protocol, parse=False)  # Hash only; never parse/recreate measurement arrays.
    for ref in manifest["source_binding_artifacts"]:
        check_ref(ref, protocol, "source_binding_audit")
    for path, expected in manifest["original_strong_artifact_sha256_unchanged"].items():
        checked(ROOT / path, expected)
    source = protocol["source"]
    checked(ROOT / source["custodian_only_path"], source["source_sha256"], source["source_bytes"])
    require(configuration["document_type"] == "population_development_executor_configuration" and configuration["requested_phase"] == "after_phase_b_export_before_native_fitting", "Wrong configuration phase or type")
    require(request["executor_code"] == configuration["code_artifacts"] + [context["new_configuration"]], "Request does not bind the exact current code/configuration")
    require(request["exporter_code"] == manifest["exporter_code_artifacts"], "Exporter freeze changed")
    for ref in request["executor_code"] + request["exporter_code"]:
        check_ref(ref, protocol, "execution_code", parse=False)
    require(len(configuration["code_artifacts"]) == len(configuration["code_snapshots"]) == 3, "Incomplete code or snapshot inventory")
    for current, snapshot in zip(configuration["code_artifacts"], configuration["code_snapshots"]):
        require(current["sha256"] == snapshot["sha256"] and current["bytes"] == snapshot["bytes"], "Snapshot differs from actual current code")
        check_ref(snapshot, protocol, "execution_code", parse=False)
    executor = fresh_executor()
    require(executor.code_artifacts(ROOT) == configuration["code_artifacts"], "Actual current source/runner/test inventory differs")
    require(executor.environment_record() == configuration["environment"], "Actual runtime/PRNG/optimizer environment differs")
    contract = executor.executor_refresh_contract()
    require(context["required_custodian_admission"] == contract and contract["required_fields"] == ADMISSION_FIELDS, "Closed consumer admission contract differs from the request")
    executor.validate_protocol(protocol)
    protocol_from_gate, _, provenance = executor._admission(ROOT, PROTOCOL, ROOT / context["phase_b_manifest"]["path"], "training_and_conditioning_export")
    require(protocol_from_gate == protocol and provenance["adapter_version"] == "phase_b_release_manifest_v1" and provenance["data_ready"] is True and provenance["executor_ready"] is False, "Original phase-B adapter did not validate or is already refreshed")
    require(provenance["approved_data_artifacts"] == context["reuse_existing_projections"], "Consumer sees different phase-B projections")
    for path, ref in list(VERIFIED.items()):
        checked(ROOT / path, ref["sha256"], ref["bytes"])
    return protocol, context, request, configuration, manifest, previous_approval, executor


def issue():
    protocol, context, request, configuration, manifest, previous_approval, executor = verify_inputs()
    event_path = OUT / "refresh_event.json"
    if event_path.exists():
        event = strict_json(regular(event_path).read_bytes())
        require(event["refresh_request_sha256"] == PINS["refresh_request"] and event["request_sha256"] == PINS["request"] and event["configuration_sha256"] == PINS["configuration"], "Existing refresh event has different bindings")
    else:
        # Only filenames are inspected; no fitting results or model scores are
        # read. This is an operational marker check, not proof against unlogged
        # activity outside the accountable workflow.
        markers = []
        output_root = ROOT / "outputs/native_population_development"
        for name in ("training_access_intent.json", "development_prediction_freeze.json"):
            markers.extend(relative(path) for path in output_root.rglob(name))
        require(not markers, "Native fit/access marker exists; before-first-fit refresh requires renewed custody review")
        event = {
            "document_type": "executor_refresh_custody_event_not_external_attestation",
            "event_id": "granados-executor-refresh-01", "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "phase": "after_phase_b_export_before_first_native_fit",
            "chronology": "Phase-B data were already exported. This is a NEW executor/configuration decision after POP001/POP003/POP004; neither the original approval nor the historical executor SHA is rewritten or backdated.",
            "refresh_request_path": relative(INPUT / "executor_refresh_request.json"),
            "refresh_request_sha256": PINS["refresh_request"], "request_sha256": PINS["request"],
            "configuration_sha256": PINS["configuration"], "original_phase_b_manifest_sha256": PINS["phase_b_manifest"],
            "original_data_approval_sha256": PINS["previous_approval"],
            "previous_executor_artifacts": context["previous_executor_artifacts"],
            "new_executor_artifacts": configuration["code_artifacts"],
            "observed_native_fit_access_markers_before_issue": markers,
            "fit_history_qualification": "The parent/request state before-first-native-fit; no designated native fit/access marker was found. No external immutable chronology or protection against unlogged activity is claimed.",
            "new_data_exported": False, "native_fits_or_scoring_performed_by_custodian": False,
            "development_responses_released": False, "strong_claim_authorized": False,
        }
    write_new(event_path, event)
    identity = previous_approval["custodian_identity"]
    require(identity != request["requester"], "Custodian and executor requester must be distinct roles")
    decision = {
        "document_type": "development_custody_decision_not_external_attestation", "protocol_id": PROTOCOL_ID,
        "request_sha256": PINS["request"], "decision": "approve_development_exports",
        "custodian_identity": identity, "operational_event_reference": relative(event_path) + "#granados-executor-refresh-01",
        "reviewed_source_asset_id": protocol["source"]["asset_id"], "reviewed_source_sha256": protocol["source"]["source_sha256"],
        "permitted_group_ids": request["train_ids"] + request["development_ids"],
        "conditions": [
            "This new decision admits the exact POP001/POP003/POP004 executor/request/configuration refresh for reuse of the already approved phase-B data. It is after data export and before the first native fit, not a backdated pre-value decision.",
            "All bound current source, runner, independent-test, snapshot, exporter, configuration and protocol artifacts must retain the freshly verified SHA256 values. The historical executor artifacts and old request/approval remain unchanged records, not asserted equal to the new code.",
            "Reuse ONLY original training projection fc8000c7720e554756c90e0da79f2d0ae4a7928d5923f0efa1a2ed5823d37c9e and development prefix/timing projection 41ee2a8ce6b09cb664859873fd59fd005e31e6eaa054c3703098c91e3f790613 under original phase-B manifest 3b8c67c79191bf4dec0a1e7e0514b1a145359c30061868bbfdd136f0ee1cd3dd and data approval 87536dc42cccde2733daf4cb52f54024661a7d58386acd180623adf266126533. No data regeneration or additional export is authorized.",
            "Do not change the scientific recipe, groups, raw ratio operator, windows, prefix cohort rule, budget or model-family set. Source-verifier review refresh is separate and is not represented here as completed.",
            "Development response values remain withheld until all six canonical fits, all attempts and all 288 development predictions are frozen and the matching digests receive a separate custody review. Rep6, aliases and every other unapproved group/field/window remain unexported.",
            "This is limited retrospective published-processed-cohort development only. No original strong-protocol status, biological ready grade, certified prospective forecast or external authority is granted.",
        ],
        "limitations": previous_approval["limitations"], "strong_claim_authorized": False,
    }
    validate_type(protocol, decision, "DevelopmentCustodyDecision")
    executor.validate_document(protocol, decision, "DevelopmentCustodyDecision")
    approval_ref = write_new(OUT / "custody_decision.json", decision)
    approval_ref["role"] = "operational_approval"
    admission = {
        "document_type": "population_development_executor_refresh_admission_v1", "schema_version": 1,
        "protocol_id": PROTOCOL_ID, "protocol_sha256": PINS["protocol"],
        "phase_b_manifest": context["phase_b_manifest"], "previous_approval_sha256": PINS["previous_approval"],
        "request": context["new_request"], "approval": approval_ref, "configuration": context["new_configuration"],
        "custodian_identity": identity, "operational_event_reference": decision["operational_event_reference"],
        "reuse_existing_projections": context["reuse_existing_projections"], "scientific_recipe_changed": False,
        "data_reexport_permitted": False, "development_responses_released": False, "strong_claim_authorized": False,
    }
    require(list(admission) == ADMISSION_FIELDS, "Admission fields/order differ from the explicit v1 contract")
    admission_path = OUT / "executor_refresh_admission.json"
    admission_ref = write_new(admission_path, admission)
    # Exercise the exact consumer admission path only. Do NOT call
    # execute_training, fit_training, predict_development, or scoring routines.
    gate = executor.check_custody_gate(ROOT, PROTOCOL, ROOT / context["phase_b_manifest"]["path"], executor_admission_path=admission_path)
    require(gate["data_ready"] is True and gate["executor_ready"] is True and gate["strong_claim_authorized"] is False, "Fresh consumer did not accept the refreshed operational admission")
    require(gate["measurement_arrays_opened"] is False, "Admission check unexpectedly opened measurements")
    for path, ref in list(VERIFIED.items()):
        checked(ROOT / path, ref["sha256"], ref["bytes"])
    report = {
        "document_type": "executor_refresh_custody_verification_not_biological_validation", "protocol_id": PROTOCOL_ID,
        "decision": approval_ref, "admission": admission_ref, "refresh_request_sha256": PINS["refresh_request"],
        "verified_artifacts": list(VERIFIED.values()),
        "consumer_function": "native_population_development.check_custody_gate with explicit executor_admission_path",
        "consumer_contract": "population_development_executor_refresh_admission_v1",
        "consumer_data_admission_accepted": gate["data_ready"], "consumer_executor_admission_accepted": gate["executor_ready"],
        "consumer_measurement_arrays_opened": gate["measurement_arrays_opened"],
        "model_fitting_or_scoring_performed": False, "scientific_recipe_changed": False,
        "data_regenerated_or_reexported": False, "development_responses_released": False,
        "source_verifier_review_refreshed_here": False, "strong_claim_authorized": False,
    }
    verification_ref = write_new(OUT / "verification.json", report)
    print(json.dumps({"decision": approval_ref, "admission": admission_ref, "verification": verification_ref, "operational_consumer_admission_accepted": True, "development_responses_released": False, "strong_claim_authorized": False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("inspect", "issue"))
    args = parser.parse_args()
    if args.command == "inspect":
        protocol, context, request, configuration, manifest, previous_approval, executor = verify_inputs()
        write_new(OUT / "views/refresh_request.json", context)
        write_new(OUT / "views/release_request.json", request)
        write_new(OUT / "views/execution_configuration.json", configuration)
        print(json.dumps({"all_input_bindings_verified": True, "artifact_count": len(VERIFIED), "admission_contract": context["required_custodian_admission"], "no_fitting_or_data_regeneration": True}))
    else:
        issue()


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        if isinstance(error, RefreshError):
            sys.exit(str(error))
        sys.exit("Executor refresh stopped without data release: " + type(error).__name__)
