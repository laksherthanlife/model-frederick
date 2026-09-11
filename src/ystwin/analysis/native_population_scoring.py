from __future__ import annotations

import hashlib
import math
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ystwin.analysis import native_population_development as frozen


ADAPTER_VERSION = "native_population_phase_c_scoring_v2"
MANIFEST_SHA256 = "0e2d18b1d239f33fd21165219c512b52546635b733ec47d67c85f726a919ec82"
FREEZE_SHA256 = "e03f6b9407b29e7e76da7a0be566ce9036c5cf27ac65bfaa9a4d3d7e009b5753"
MODULE_SHA256 = "654dd5729d204c5963f5e51576a4830b61f736296c877ffc419efa77d009e2a0"
PROTOCOL_SHA256 = "0f53f40699ebd0db6df62f00e0a0e6bfc5b56b374ec84a6e36846a8488e6dad6"
MANIFEST_PATH = "data/native_law_v2/granados/development/phase_c_01/release_manifest.json"
PROTOCOL_PATH = "data/native_law_v2/development_protocol.json"
SCIENTIFIC_SCOPE = "retrospective_source_cohort_development"
REFERENCE_ROLES = {
    "release_request": "operational_event_log",
    "custody_decision": "operational_approval",
    "freeze_verification": "source_binding_audit",
    "prediction_freeze": "development_predictions",
    "attempt_log": "attempt_log",
    "conditioning_checkpoint": "model_checkpoint",
    "training_run_summary": "operational_event_log",
    "phase_b_manifest": "operational_event_log",
    "original_training_data_approval": "operational_approval",
    "executor_refresh_admission": "operational_event_log",
    "executor_refresh_approval": "operational_approval",
    "configuration": "execution_code",
    "response_projection": "development_responses",
    "previously_released_development_inputs": "development_inputs",
    "source_binding_audit": "source_binding_audit",
    "source_replay_verification": "source_binding_audit",
    "coverage_audit": "source_binding_audit",
}
MANIFEST_FIELDS = set(REFERENCE_ROLES) | {
    "document_type", "protocol_id", "protocol_sha256", "phase", "frozen_executor_artifacts",
    "response_data_allowlist", "operational_event_reference", "development_responses_released",
    "previous_data_regenerated", "frozen_cohorts_or_coordinates_changed",
    "all_registered_observation_means_and_coordinates_available", "current_module_sha256",
    "legacy_execute_scoring_approval_equality_compatible", "scoring_integration_requirement",
    "selected_development_model", "development_scores_computed", "strong_claim_authorized",
    "final_test_defined", "limitations",
}


class ScoringError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ScoringError(message)


@dataclass(frozen=True)
class TrainingAdmission:
    phase_b_manifest: dict
    data_approval: dict
    executor_approval: dict
    executor_admission: dict
    provenance: dict


@dataclass(frozen=True)
class PredictionFreeze:
    reference: dict
    data: dict
    conditioning: dict
    attempt_log: dict
    summary: dict
    map_audit: dict
    checkpoints: dict


@dataclass(frozen=True)
class PhaseCResponseRelease:
    manifest_reference: dict
    manifest: dict
    request: dict
    decision: dict
    source_binding: dict
    source_replay: dict
    freeze_verification: dict
    event: dict


@dataclass(frozen=True, order=True)
class ArtifactPin:
    path: str
    sha256: str
    bytes: int
    role: str

    @classmethod
    def from_reference(cls, reference):
        require(all(type(reference[key]) is str for key in ("path", "sha256", "role")) and type(reference["bytes"]) is int, "artifact pins must contain immutable scalar identities")
        return cls(reference["path"], reference["sha256"], reference["bytes"], reference["role"])


@dataclass(frozen=True)
class AdmissionPins:
    authority_root: str
    manifest: ArtifactPin
    artifacts: tuple[ArtifactPin, ...]

    def __post_init__(self):
        require(type(self.authority_root) is str and isinstance(self.manifest, ArtifactPin) and type(self.artifacts) is tuple and all(isinstance(pin, ArtifactPin) for pin in self.artifacts), "admission pins must be immutable, typed identities")


@dataclass(frozen=True)
class ScoringBundle:
    protocol: dict
    training: TrainingAdmission
    prediction: PredictionFreeze
    release: PhaseCResponseRelease
    verified_artifacts: tuple
    admission_pins: AdmissionPins


class MetadataReader:
    def __init__(self, root, protocol, deferred=()):
        self.root = Path(root).resolve()
        self.protocol = protocol
        self.deferred = set(deferred)
        self.verified = {}

    def read(self, reference, role, parse=True):
        require(role not in {"training_projection", "development_inputs", "development_responses"}, "measurement projections require the separately admitted value-reading step")
        require(reference["path"] not in self.deferred, "outcome-bearing audit cannot be read during metadata admission")
        value = frozen.read_artifact(self.root, reference, self.protocol, role, parse=parse)
        self.verified[reference["path"]] = dict(reference)
        return value


def _manifest_contract(protocol, manifest):
    require(isinstance(manifest, dict) and set(manifest) == MANIFEST_FIELDS, "phase-C release-manifest-v1 fields mismatch")
    require(manifest["document_type"] == "approved_development_phase_c_response_release_v1" and manifest["phase"] == "development_response_export", "wrong release phase or version")
    require(manifest["protocol_id"] == frozen.PROTOCOL_ID and manifest["protocol_sha256"] == PROTOCOL_SHA256, "phase-C scientific protocol mismatch")
    require(manifest["development_responses_released"] is True, "phase-C responses have not been released")
    for field in ("previous_data_regenerated", "frozen_cohorts_or_coordinates_changed", "development_scores_computed", "strong_claim_authorized", "final_test_defined"):
        require(manifest[field] is False, f"phase-C invariant failed: {field}")
    require(manifest["selected_development_model"] is None, "a release cannot preselect a development winner")
    require(manifest["response_data_allowlist"] == [manifest["response_projection"]], "response allowlist must contain only the exact released projection")
    for key, role in REFERENCE_ROLES.items():
        frozen.validate_document(protocol, manifest[key], "ArtifactRef")
        require(manifest[key]["role"] == role, f"wrong artifact role: {key}")
    require(len(manifest["frozen_executor_artifacts"]) == 3 and manifest["current_module_sha256"] == MODULE_SHA256, "frozen executor inventory mismatch")


def inspect_phase_c(root, manifest_path=MANIFEST_PATH, expected_manifest_sha=MANIFEST_SHA256):
    root = Path(root).resolve()
    path = Path(manifest_path)
    relative = path.relative_to(root) if path.is_absolute() else path
    path = frozen.safe_path(root, relative)
    frozen._artifact_path_policy(root, path, "operational_event_log")
    require(path.is_file() and path.stat().st_size < 1024 * 1024, "phase-C metadata file is missing or exceeds the metadata byte budget")
    content = path.read_bytes()
    require(hashlib.sha256(content).hexdigest() == expected_manifest_sha, "phase-C manifest digest mismatch")
    manifest = frozen.strict_json(content)
    manifest_ref = {"path": relative.as_posix(), "sha256": expected_manifest_sha, "bytes": len(content), "role": "operational_event_log"}
    module_ref = frozen.artifact_reference(root, "src/ystwin/analysis/native_population_development.py", "execution_code")
    require(module_ref["sha256"] == MODULE_SHA256, "frozen scientific helper code digest changed")
    protocol = frozen.scoped_protocol(root, PROTOCOL_PATH)
    protocol_ref = frozen.artifact_reference(root, PROTOCOL_PATH, "protocol_snapshot")
    require(protocol_ref["sha256"] == PROTOCOL_SHA256, "scientific protocol bytes changed")
    _manifest_contract(protocol, manifest)
    require(manifest["prediction_freeze"]["sha256"] == FREEZE_SHA256, "phase-C release is not bound to the registered prediction freeze")
    deferred = [manifest[key]["path"] for key in ("response_projection", "previously_released_development_inputs", "coverage_audit")]
    reader = MetadataReader(root, protocol, deferred)
    documents = {key: reader.read(reference, REFERENCE_ROLES[key]) for key, reference in manifest.items() if key in REFERENCE_ROLES and key not in {"response_projection", "previously_released_development_inputs", "coverage_audit"}}
    _, _, training_provenance = frozen._admission(root, root / PROTOCOL_PATH, root / manifest["phase_b_manifest"]["path"], "training_and_conditioning_export", executor_admission_path=root / manifest["executor_refresh_admission"]["path"])
    require(training_provenance.get("executor_ready") is True, "training executor admission is not valid")
    require(training_provenance["phase_b_manifest"] == manifest["phase_b_manifest"] and training_provenance["approval_artifact"] == manifest["original_training_data_approval"] and training_provenance["executor_approval_artifact"] == manifest["executor_refresh_approval"] and training_provenance["executor_admission_sha256"] == manifest["executor_refresh_admission"]["sha256"] and training_provenance["configuration_artifact"] == manifest["configuration"], "training/executor lineage is not the referenced phase-B lineage")
    require(training_provenance["frozen_code_artifacts"] == manifest["frozen_executor_artifacts"], "admitted code copies differ from the prediction-seal code")
    for reference in manifest["frozen_executor_artifacts"]:
        reader.read(reference, "execution_code", parse=False)
    for reference in documents["release_request"]["executor_code"] + documents["release_request"]["exporter_code"]:
        reader.read(reference, "execution_code", parse=False)
    event_path, separator, event_id = manifest["operational_event_reference"].partition("#")
    require(separator and event_id == "granados-phase-c-01", "phase-C event reference mismatch")
    event_ref = frozen.artifact_reference(root, event_path, "operational_event_log")
    event = reader.read(event_ref, "operational_event_log")
    sealed = documents["prediction_freeze"]
    frozen.validate_freeze(protocol, sealed)
    require(documents["training_run_summary"]["nonnegative_map_audit"] == sealed["environment"]["nonnegative_map_audit"], "domain-map audit lineage mismatch")
    map_audit = reader.read(sealed["environment"]["nonnegative_map_audit"], "model_checkpoint")
    checkpoints = {}
    for attempt in sealed["attempts"]:
        if attempt["status"] != "ok":
            continue
        checkpoint = reader.read(attempt["checkpoint"], "model_checkpoint")
        require(checkpoint["model_id"] == attempt["model_id"] and checkpoint["parameters"] == attempt["parameters"] and checkpoint["start_index"] == attempt["start_index"] and checkpoint["seed"] == attempt["seed"] and checkpoint["training_mean_mse"] == attempt["training_mean_mse"] and checkpoint["training_projection_sha256"] == sealed["training_projection_sha256"], "checkpoint is not the logged frozen training attempt")
        require(checkpoint["bounds"] == frozen.parameter_bounds(attempt["model_id"], sealed["H_train"]) and checkpoint["H_train"] == sealed["H_train"] and checkpoint["S_train"] == sealed["S_train"], "checkpoint bounds/scales differ from the training freeze")
        checkpoints[attempt["checkpoint"]["path"]] = checkpoint
    references = [manifest_ref, module_ref, protocol_ref, *reader.verified.values(), *(manifest[key] for key in REFERENCE_ROLES)]
    pins = AdmissionPins(str(root), ArtifactPin.from_reference(manifest_ref), tuple(sorted({ArtifactPin.from_reference(reference) for reference in references})))
    bundle = ScoringBundle(protocol, TrainingAdmission(manifest["phase_b_manifest"], manifest["original_training_data_approval"], manifest["executor_refresh_approval"], manifest["executor_refresh_admission"], training_provenance), PredictionFreeze(manifest["prediction_freeze"], sealed, documents["conditioning_checkpoint"], documents["attempt_log"], documents["training_run_summary"], map_audit, checkpoints), PhaseCResponseRelease(manifest_ref, manifest, documents["release_request"], documents["custody_decision"], documents["source_binding_audit"], documents["source_replay_verification"], documents["freeze_verification"], event), tuple(reader.verified.values()), pins)
    validate_lineage(bundle)
    return bundle


def validate_lineage(bundle):
    protocol, manifest = bundle.protocol, bundle.release.manifest
    sealed, training = bundle.prediction.data, bundle.training
    _manifest_contract(protocol, manifest)
    approvals = (training.data_approval["sha256"], training.executor_approval["sha256"], manifest["custody_decision"]["sha256"])
    require(len(set(approvals)) == 3, "training, executor and response approvals must be distinct linked records")
    require(manifest["original_training_data_approval"] == training.data_approval and manifest["executor_refresh_approval"] == training.executor_approval and manifest["executor_refresh_admission"] == training.executor_admission and manifest["phase_b_manifest"] == training.phase_b_manifest, "approval-role references were reused or substituted")
    require(sealed["approval_sha256"] == approvals[0], "prediction freeze must retain its original training data approval")
    provenance = sealed["environment"]["provenance"]
    require(provenance["approval_artifact"] == training.data_approval and provenance["executor_approval_artifact"] == training.executor_approval and provenance["executor_admission_sha256"] == training.executor_admission["sha256"], "frozen training admission does not link the correct executor refresh")
    require(provenance["phase_b_manifest"] == training.phase_b_manifest and provenance["configuration_artifact"] == manifest["configuration"], "freeze configuration/training manifest lineage mismatch")
    require(sealed["code_artifacts"] == manifest["frozen_executor_artifacts"] == provenance["frozen_code_artifacts"], "frozen code lineage differs across release and checkpoint")
    require(provenance["development_inputs"] == manifest["previously_released_development_inputs"] and sealed["development_inputs_sha256"] == manifest["previously_released_development_inputs"]["sha256"] and sealed["training_projection_sha256"] == provenance["training_projection"]["sha256"], "frozen training/conditioning data lineage mismatch")
    frozen.validate_freeze(protocol, sealed)
    require(bundle.prediction.reference == manifest["prediction_freeze"], "substituted prediction-freeze reference")
    require(bundle.prediction.summary["prediction_freeze"] == manifest["prediction_freeze"] and bundle.prediction.summary["attempt_log"] == manifest["attempt_log"] and bundle.prediction.summary["conditioning_checkpoint"] == manifest["conditioning_checkpoint"], "run summary references differ from the response release")
    require(sealed["environment"]["attempt_log"] == manifest["attempt_log"] and sealed["environment"]["conditioning_checkpoint"] == manifest["conditioning_checkpoint"], "frozen attempt/conditioning references differ from the release")
    require(bundle.prediction.conditioning["prefix_summaries"] == sealed["prefix_summaries"] and bundle.prediction.conditioning["H_train"] == sealed["H_train"] and bundle.prediction.conditioning["S_train"] == sealed["S_train"], "frozen cohort or training-only scales changed")
    require(set(bundle.prediction.attempt_log["families"]) == set(frozen.MODEL_IDS) and set(bundle.prediction.summary["selected_checkpoints"]) == set(frozen.MODEL_IDS), "canonical family omission in the attempt log or summary")
    for model_id in frozen.MODEL_IDS:
        official = [row for row in sealed["attempts"] if row["model_id"] == model_id]
        detailed = bundle.prediction.attempt_log["families"][model_id]["attempts"]
        require(len(official) == len(detailed), "attempt log count differs from the prediction freeze")
        for row, full in zip(official, detailed):
            require(all(row[key] == full[key] for key in ("model_id", "start_index", "seed", "parameters", "loss_evaluations", "status", "training_mean_mse", "reason")), "attempt log changed after training")
        selected = frozen._select_start(official, frozen.PARAMETERS[model_id])
        expected = selected["checkpoint"] if selected else None
        require(bundle.prediction.summary["selected_checkpoints"][model_id] == expected and sealed["environment"]["selected_checkpoints"][model_id] == expected, "training-selected start/checkpoint changed")
    decision, request = bundle.release.decision, bundle.release.request
    frozen.validate_document(protocol, decision, "DevelopmentCustodyDecision")
    frozen.validate_document(protocol, request, "DevelopmentReleaseRequest")
    require(decision["decision"] == "approve_development_exports" and decision["request_sha256"] == manifest["release_request"]["sha256"] and decision["custodian_identity"] != request["requester"], "response release lacks its own matching independent decision")
    require(request["protocol_sha256"] == manifest["protocol_sha256"] and request["executor_code"] == provenance["code_artifacts"] + [manifest["configuration"]], "response request does not bind the executed scientific code/configuration")
    require(set(decision["permitted_group_ids"]) == {row["group_id"] for row in frozen.group_specs(protocol)}, "response decision group inventory mismatch")
    require(decision["operational_event_reference"] == manifest["operational_event_reference"], "response decision is not linked to the release event")
    event = bundle.release.event
    event_links = {"prediction_freeze_sha256": manifest["prediction_freeze"]["sha256"], "freeze_verification_sha256": manifest["freeze_verification"]["sha256"], "original_phase_b_manifest_sha256": training.phase_b_manifest["sha256"], "original_data_approval_sha256": approvals[0], "executor_admission_sha256": training.executor_admission["sha256"], "request_sha256": manifest["release_request"]["sha256"]}
    require(event["document_type"] == "phase_c_custody_approval_event_not_external_attestation" and event["event_id"] == "granados-phase-c-01" and event["phase"] == "after_complete_prediction_freeze_before_response_extraction" and all(event[key] == value for key, value in event_links.items()), "response event does not link training admission, seal and separate response release")
    require(all(event[key] is False for key in ("fitting_scoring_or_family_selection_performed", "scientific_recipe_changed", "prior_releases_or_events_rewritten", "strong_claim_authorized")), "response event changed the frozen scientific procedure")
    verification = bundle.release.freeze_verification
    require(verification["document_type"] == "phase_c_sealed_freeze_verification_not_scoring" and verification["protocol_id"] == frozen.PROTOCOL_ID and verification["protocol_sha256"] == manifest["protocol_sha256"] and verification["freeze"] == manifest["prediction_freeze"] and verification["attempt_log"] == manifest["attempt_log"] and verification["conditioning_checkpoint"] == manifest["conditioning_checkpoint"] and verification["training_run_summary"] == manifest["training_run_summary"] and verification["original_data_approval_sha256"] == approvals[0], "pre-response freeze verification binds different artifacts")
    require(verification["canonical_families"] == list(frozen.MODEL_IDS) and verification["attempt_count"] == len(sealed["attempts"]) and verification["attempt_counts"] == dict(Counter(row["model_id"] for row in sealed["attempts"])) and verification["prediction_count"] == len(sealed["predictions"]), "custodial freeze inventory differs from the actual sealed inventory")
    actual_status = {model: dict(Counter(row["status"] for row in sealed["predictions"] if row["model_id"] == model)) for model in frozen.MODEL_IDS}
    require(verification["prediction_status_counts"] == actual_status and all(verification[key] is False for key in ("model_fitting_performed", "family_selection_performed", "development_scoring_performed", "development_responses_opened")), "freeze verification is not a pre-response verification of these predictions")
    binding, replay = bundle.release.source_binding, bundle.release.source_replay
    require(binding["document_type"] == "phase_c_source_scalar_binding_not_biological_validation" and binding["protocol_id"] == frozen.PROTOCOL_ID and binding["protocol_sha256"] == manifest["protocol_sha256"] and binding["source_asset_id"] == frozen.ASSET_ID and binding["source_sha256"] == frozen.SOURCE_SHA256, "response source binding has a different source or protocol")
    require(binding["prediction_freeze_sha256"] == manifest["prediction_freeze"]["sha256"] and binding["response_release_approval_sha256"] == approvals[2] and binding["original_data_approval_sha256"] == approvals[0] and binding["freeze_verification"] == manifest["freeze_verification"], "response binding conflates or substitutes phase-specific approvals")
    dev_ids = [row["group_id"] for row in frozen.group_specs(protocol)[3:]]
    require(binding["group_ids"] == dev_ids and binding["frames_1based"] == list(frozen.RESPONSE_FRAMES) and binding["expected_source_cells"] == 8448 and binding["expected_source_scalars"] == 25344 and binding["strong_claim_authorized"] is False, "response source window/inventory mismatch")
    require(replay["document_type"] == "phase_c_serialized_source_replay_not_scoring" and replay["protocol_id"] == frozen.PROTOCOL_ID and replay["prediction_freeze_sha256"] == manifest["prediction_freeze"]["sha256"] and replay["response_projection"] == manifest["response_projection"] and replay["source_binding_audit"] == manifest["source_binding_audit"], "source replay binds a different response artifact or seal")
    require(replay["frozen_prefix_cohorts_used"] == {row["group_id"]: row["eligible_cell_indices"] for row in sealed["prefix_summaries"] if row["group_id"] in dev_ids}, "response replay changed the frozen cohorts")
    checks = replay["checks"]
    require(checks["source_cell_count"] == 8448 and checks["source_scalar_count"] == 25344 and all(checks[key] is True for key in ("all_source_scalars_equal", "nulls_preserved", "source_rows_and_pointers_verified", "group_window_inventory_verified", "all_times_equal_previously_released_source_times")), "custodial replay did not preserve the complete source scalar/time/row inventory")
    require(all(replay[key] is False for key in ("development_scores_computed", "model_family_selected", "strong_claim_authorized")), "response replay already selected or scored model outcomes")
    return True


def reload_admission(pins):
    require(isinstance(pins, AdmissionPins), "reloading requires the immutable pins established during actual admission")
    fresh = inspect_phase_c(pins.authority_root, pins.manifest.path, pins.manifest.sha256)
    require(fresh.admission_pins == pins, "pinned artifact identity, inventory or digest changed since admission")
    return fresh


def revalidate_admission(bundle):
    require(isinstance(bundle, ScoringBundle), "a genuinely inspected scoring bundle is required")
    return reload_admission(bundle.admission_pins)


def load_released_observations(root, bundle):
    bundle = revalidate_admission(bundle)
    manifest, protocol = bundle.release.manifest, bundle.protocol
    try:
        response = frozen.read_artifact(root, manifest["response_projection"], protocol, "development_responses")
    except (frozen.DevelopmentError, OSError) as error:
        raise ScoringError("response projection must bind its separate phase-C approval and the pinned source row/frame inventory; authoritative digest/byte verification failed") from error
    frozen.validate_projection(protocol, response, "development_responses")
    require(response["approval_sha256"] == manifest["custody_decision"]["sha256"], "response projection must bind its separate phase-C approval, not a training or executor approval")
    require(response["source_binding_audit"] == manifest["source_binding_audit"], "response projection source-binding reference differs from the authorized release")
    payload = {key: response[key] for key in ("group_ids", "source_cells", "source_times")}
    require(frozen.digest(payload) == bundle.release.source_binding["ordered_payload_sha256"], "response payload digest differs from direct source-scalar replay")
    nulls = {key: sum(row[key] is None for row in response["source_cells"]) for key in ("max5", "median", "relative_time_min")}
    require(nulls == bundle.release.source_replay["checks"]["null_counts"], "response null preservation differs from the source replay")
    conditioning = frozen.read_artifact(root, manifest["previously_released_development_inputs"], protocol, "development_inputs")
    frozen.validate_projection(protocol, conditioning, "development_inputs")
    require(conditioning["approval_sha256"] == bundle.training.data_approval["sha256"], "previous conditioning projection lost its original training-data approval")
    prefix = {row["group_id"]: row for row in bundle.prediction.data["prefix_summaries"]}
    coordinates = {(row["group_id"], row["frame_1based"]): row["time_min"] for row in bundle.prediction.conditioning["development_coordinates"]}
    frames = []
    for group in frozen.group_specs(protocol)[3:]:
        cells = [row for row in response["source_cells"] if row["group_id"] == group["group_id"]]
        times = [row for row in conditioning["source_times"] if row["group_id"] == group["group_id"]]
        time_index = {(row["cell_index"], row["frame_1based"]): row["relative_time_min"] for row in times}
        require(all(row["relative_time_min"] == time_index[row["cell_index"], row["frame_1based"]] for row in cells), "response times changed from the pre-response timing projection")
        derived = frozen.derive_frames(group, cells, times, prefix[group["group_id"]])
        require(all(row["time_min"] == coordinates[row["group_id"], row["frame_1based"]] for row in derived), "aggregation recomputed different frozen response coordinates")
        frames.extend(derived)
    require(len(frames) == 48, "all 48 response frame summaries must be retained")
    coverage = frozen.read_artifact(root, manifest["coverage_audit"], protocol, "source_binding_audit")
    require(coverage["document_type"] == "phase_c_frozen_cohort_response_coverage_not_model_scoring" and coverage["protocol_id"] == frozen.PROTOCOL_ID and coverage["prediction_freeze_sha256"] == bundle.prediction.reference["sha256"] and coverage["cohorts_recomputed_from_responses"] is False and coverage["no_family_selected_or_scored"] is True and coverage["strong_claim_authorized"] is False, "coverage audit is not linked to the frozen cohorts and predictions")
    derived_index = {(row["group_id"], row["frame_1based"]): row for row in frames}
    seen = set()
    for group in coverage["groups"]:
        require(group["frozen_prefix_cell_indices"] == prefix[group["group_id"]]["eligible_cell_indices"], "coverage audit substituted a response-based cohort")
        for audited in group["frames"]:
            key = group["group_id"], audited["frame_1based"]
            require(key in derived_index and key not in seen, "coverage audit has duplicate or unregistered frames")
            seen.add(key)
            derived = derived_index[key]
            require(audited["frozen_time_min"] == derived["time_min"] and all(audited[field] == derived[field] for field in ("prefix_cells", "valid_cells", "source_missing_cells", "invalid_cells")), "derived frame validity/time differs from independent source coverage")
            require(audited["mean_is_defined"] == (derived["observed_mean"] is not None) and audited["within_cell_sd_is_defined"] == (derived["cell_sample_sd"] is not None), "source coverage and computed aggregate availability disagree")
            require([row["cell_index"] for row in audited["source_missing_details"]] == [row["cell_index"] for row in derived["source_missing_details"]] and audited["invalid_details"] == derived["invalid_details"], "response exclusions differ from the source-bound audit")
    require(seen == set(derived_index), "coverage audit omitted a registered development frame")
    complete = all(row["observed_mean"] is not None and row["time_min"] is not None for row in frames)
    require(coverage["all_registered_means_and_coordinates_available"] == complete == manifest["all_registered_observation_means_and_coordinates_available"], "declared completeness differs from the actual response aggregates")
    return frames


def score_frozen_means(bundle, observations):
    bundle = revalidate_admission(bundle)
    report = frozen.build_development_report(bundle.protocol, bundle.prediction.data, observations, bundle.prediction.reference["sha256"], bundle.release.manifest["response_projection"]["sha256"])
    models, bound_warnings, map_bindings = {}, {}, {}
    scale = bundle.prediction.data["S_train"]
    for model_id in frozen.MODEL_IDS:
        groups = [row for row in report["group_evaluations"] if row["model_id"] == model_id]
        complete = all(row["complete"] for row in groups)
        mse = float(np.mean([row["mean_mse"] for row in groups])) if complete else None
        rmse = math.sqrt(mse) if complete else None
        mae = float(np.mean([row["mae"] for row in groups])) if complete else None
        models[model_id] = {"complete": complete, "mean_mse": mse, "rmse": rmse, "mae": mae, "scaled_rmse": rmse / scale if complete and scale > 0 else None, "scaled_mae": mae / scale if complete and scale > 0 else None, "groups": {row["group_id"]: {key: row[key] for key in ("complete", "mean_mse", "rmse", "mae")} for row in groups}}
        reference = bundle.prediction.summary["selected_checkpoints"][model_id]
        if reference is not None:
            checkpoint = bundle.prediction.checkpoints[reference["path"]]
            hits = {name: "lower" if math.isclose(value, checkpoint["bounds"][name][0], rel_tol=1e-9, abs_tol=1e-9) else "upper" for name, value in checkpoint["parameters"].items() if any(math.isclose(value, bound, rel_tol=1e-9, abs_tol=1e-9) for bound in checkpoint["bounds"][name])}
            if hits:
                bound_warnings[model_id] = hits
        map_bindings[model_id] = {part: sum(bool(row["nonnegative_map_binds"]) for row in bundle.prediction.map_audit[part] if row["model_id"] == model_id) for part in ("training", "development")}
    summary = {"adapter_version": ADAPTER_VERSION, "scientific_scope": SCIENTIFIC_SCOPE, "status": report["status"], "models": models, "selected_development_model": report["selected_development_model"], "H_train": bundle.prediction.data["H_train"], "S_train": scale, "parameter_bound_warnings": bound_warnings, "nonnegative_output_map_bindings": map_bindings, "prediction_count": len(bundle.prediction.data["predictions"]), "evaluated_frame_rows": sum(len(row["frames"]) for row in report["group_evaluations"]), "limitations": [*frozen.LIMITATIONS, "This scorer/admission adapter was created after prediction sealing and response availability as an operational integration fix; it is not a preregistered new biological analysis."], "strong_claim_authorized": False, "final_test_defined": False}
    return report, summary


def new_score_store(root, bundle, directory):
    relative = Path(directory)
    require(relative.parts[:2] == ("outputs", "native_population_development"), "scoring must use a separate development output directory")
    frozen_directory = Path(bundle.prediction.reference["path"]).parent
    require(relative != frozen_directory and not relative.is_relative_to(frozen_directory), "the sealed training directory is immutable")
    path = frozen.safe_path(root, relative)
    require(not path.exists(), "score directory already exists; no repeat scoring or overwrite is permitted")
    return frozen.ArtifactStore(root, relative)


def score_once(root, manifest_path, output_directory, expected_manifest_sha=MANIFEST_SHA256):
    root = Path(root).resolve()
    bundle = inspect_phase_c(root, manifest_path, expected_manifest_sha)
    store = new_score_store(root, bundle, output_directory)
    code_names = ("src/ystwin/analysis/native_population_scoring.py", "scripts/score_native_population_development.py", "tests/test_native_population_scoring.py")
    adapter_refs = [frozen.artifact_reference(root, name, "execution_code") for name in code_names]
    snapshots = [store.write_bytes(f"adapter_code/{Path(reference['path']).name}", frozen.read_artifact(root, reference, bundle.protocol, "execution_code", parse=False), "execution_code") for reference in adapter_refs]
    intent = {"document_type": "phase_c_scoring_operational_adapter_v1", "adapter_version": ADAPTER_VERSION, "scientific_scope": SCIENTIFIC_SCOPE, "created_at_utc": datetime.now(timezone.utc).isoformat(), "created_after_prediction_sealing": True, "created_after_response_availability": True, "preregistered_new_biology": False, "scientific_metric_or_prediction_change": False, "phase_c_manifest": bundle.release.manifest_reference, "training_admission": {"data_approval": bundle.training.data_approval, "executor_approval": bundle.training.executor_approval, "executor_admission": bundle.training.executor_admission}, "prediction_freeze": bundle.prediction.reference, "response_release_approval": bundle.release.manifest["custody_decision"], "response_projection": bundle.release.manifest["response_projection"], "adapter_code": snapshots, "unchanged_frozen_scientific_code": bundle.release.manifest["frozen_executor_artifacts"], "verified_lineage_artifacts": list(bundle.verified_artifacts), "strong_claim_authorized": False}
    intent_ref = store.write("scoring_access_intent.json", intent, "operational_event_log")
    observations = load_released_observations(root, bundle)
    report, summary = score_frozen_means(bundle, observations)
    report_ref = store.write("development_report.json", report, "development_report")
    observed_ref = store.write("response_aggregation_audit.json", {"scientific_scope": SCIENTIFIC_SCOPE, "response_projection": bundle.release.manifest["response_projection"], "fixed_prefix_cohorts": bundle.release.source_replay["frozen_prefix_cohorts_used"], "frames": observations, "strong_claim_authorized": False}, "development_report")
    warnings_ref = store.write("retained_training_domain_warnings.json", {"scientific_scope": SCIENTIFIC_SCOPE, "parameter_bound_warnings": summary["parameter_bound_warnings"], "nonnegative_output_map_bindings": summary["nonnegative_output_map_bindings"], "source_audit": bundle.prediction.data["environment"]["nonnegative_map_audit"], "bounds_or_predictions_changed": False, "strong_claim_authorized": False}, "development_report")
    summary.update(phase_c_manifest=bundle.release.manifest_reference, prediction_freeze=bundle.prediction.reference, response_projection=bundle.release.manifest["response_projection"], report=report_ref, aggregation_audit=observed_ref, warnings=warnings_ref, scoring_access_intent=intent_ref, adapter_code=snapshots)
    summary_ref = store.write("scoring_summary.json", summary, "development_report")
    return {**summary, "summary_artifact": summary_ref}


def verify_published_results(root, output_directory):
    root = Path(root).resolve()
    bundle = inspect_phase_c(root)
    published_ref = {"path": "outputs/native_population_development/native_training_01_phase_c_01_scoring/scoring_summary.json", "sha256": "3b2093a057f5a6819c93138fc22c7b973a4235410f99b916e74acfaa434dd7ca", "bytes": 7647, "role": "development_report"}
    published = frozen.read_artifact(root, published_ref, bundle.protocol, "development_report")
    report = frozen.read_artifact(root, published["report"], bundle.protocol, "development_report")
    frozen.validate_document(bundle.protocol, report, "DevelopmentReport")
    require(published["prediction_freeze"] == bundle.prediction.reference and published["phase_c_manifest"] == bundle.release.manifest_reference and published["response_projection"] == bundle.release.manifest["response_projection"], "published results have different sealed inputs")
    intent = frozen.read_artifact(root, published["scoring_access_intent"], bundle.protocol, "operational_event_log")
    for reference in intent["adapter_code"]:
        frozen.read_artifact(root, reference, bundle.protocol, "execution_code", parse=False)
    observed = load_released_observations(root, bundle)
    observations = {(row["group_id"], row["frame_1based"]): row for row in observed}
    predictions = {(row["model_id"], row["group_id"], row["frame_1based"]): row for row in bundle.prediction.data["predictions"]}
    expected_groups = {(model_id, group["group_id"]) for model_id in frozen.MODEL_IDS for group in frozen.group_specs(bundle.protocol)[3:]}
    require({(row["model_id"], row["group_id"]) for row in report["group_evaluations"]} == expected_groups and len(report["group_evaluations"]) == 12, "published group inventory changed")
    for group in report["group_evaluations"]:
        errors = []
        require([row["frame_1based"] for row in group["frames"]] == list(frozen.RESPONSE_FRAMES), "published frame inventory changed")
        for frame in group["frames"]:
            actual = observations[frame["group_id"], frame["frame_1based"]]
            prediction = predictions[frame["model_id"], frame["group_id"], frame["frame_1based"]]
            require(frame["status"] == "quantified" and frame["predicted_mean"] == prediction["mean"] and frame["observed_mean"] == actual["observed_mean"] and frame["time_min"] == prediction["time_min"] == actual["time_min"], "published prediction/observation values changed")
            require(all(frame[key] == actual[key] for key in ("prefix_cells", "valid_cells", "source_missing_cells", "invalid_cells", "cell_sample_sd")), "published source diagnostics changed")
            error = prediction["mean"] - actual["observed_mean"]
            require(frame["squared_error"] == error * error and frame["absolute_error"] == abs(error), "published frame residuals changed")
            errors.append(error)
        mse, mae = float(np.mean([error * error for error in errors])), float(np.mean([abs(error) for error in errors]))
        require(group["complete"] is True and group["mean_mse"] == mse and group["rmse"] == math.sqrt(mse) and group["mae"] == mae, "published group metrics changed")
    for model_id in frozen.MODEL_IDS:
        groups = [row for row in report["group_evaluations"] if row["model_id"] == model_id]
        mse, mae = float(np.mean([row["mean_mse"] for row in groups])), float(np.mean([row["mae"] for row in groups]))
        previous = published["models"][model_id]
        require(previous["mean_mse"] == mse and previous["rmse"] == math.sqrt(mse) and previous["mae"] == mae, "published family metrics changed")
        require(previous["groups"] == {row["group_id"]: {key: row[key] for key in ("complete", "mean_mse", "rmse", "mae")} for row in groups}, "published per-group summary changed")
        require(previous["scaled_rmse"] == math.sqrt(mse) / published["S_train"] and previous["scaled_mae"] == mae / published["S_train"], "published scaled diagnostics changed")
    require(report["selected_development_model"] == published["selected_development_model"], "published selection labels disagree")
    store = new_score_store(root, bundle, output_directory)
    code_names = ("src/ystwin/analysis/native_population_scoring.py", "scripts/score_native_population_development.py", "tests/test_native_population_scoring.py")
    code = [frozen.artifact_reference(root, name, "execution_code") for name in code_names]
    snapshots = [store.write_bytes(f"adapter_code/{Path(reference['path']).name}", frozen.read_artifact(root, reference, bundle.protocol, "execution_code", parse=False), "execution_code") for reference in code]
    verification = {"document_type": "phase_c_reusable_interface_fix_verification_v2", "adapter_version": ADAPTER_VERSION, "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "issue": "SCORING001", "original_adapter_version": published["adapter_version"], "original_bug": "Mutable cached metadata allowed unauthorized response copies and caller-recomputed hashes to replace the admitted projection.", "original_adapter_was_bug_free": False, "published_fresh_score_once_bypass_demonstrated": False, "fix": "Immutable pins captured during actual admission; authoritative manifest/receipt bytes and typed phase lineage are revalidated on reuse; cached metadata is not authority.", "authority_pins": asdict(bundle.admission_pins), "original_adapter_snapshots": intent["adapter_code"], "fixed_adapter_snapshots": snapshots, "published_summary": published_ref, "published_report": published["report"], "prediction_freeze": bundle.prediction.reference, "response_projection": bundle.release.manifest["response_projection"], "verified_frame_rows": 288, "verified_group_rows": 12, "verified_families": list(frozen.MODEL_IDS), "published_metrics_exactly_unchanged": True, "published_selected_development_model": published["selected_development_model"], "published_model_metrics": published["models"], "new_fit_performed": False, "new_predictions_generated": False, "new_selection_performed": False, "historical_reports_overwritten": False, "scientific_scope": SCIENTIFIC_SCOPE, "strong_claim_authorized": False}
    reference = store.write("reusable_interface_fix_verification.json", verification, "operational_event_log")
    return {"verification_artifact": reference, "adapter_version": ADAPTER_VERSION, "published_metrics_exactly_unchanged": True, "verified_frame_rows": 288, "verified_group_rows": 12, "new_fit_performed": False, "new_predictions_generated": False, "new_selection_performed": False, "published_selected_development_model": published["selected_development_model"], "scientific_scope": SCIENTIFIC_SCOPE, "strong_claim_authorized": False, "fixed_adapter_snapshots": snapshots}
