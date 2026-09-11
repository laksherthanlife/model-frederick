from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re

from .training_freeze import (
    InputSpec, _NATIVE_ROLES, _TASK_ROLES, _reject_label_overlap,
    _require_frozen_evaluation_code, _source,
)


__all__ = [
    "PortableArtifact", "PortableEvidenceBundle", "PortableEvidenceError",
    "export_portable_evidence", "load_portable_evidence",
]

CHECKPOINT_PATH = "outputs/native_training_run_02/checkpoint.json"
CHECKPOINT_SHA256 = "384175e19b2451f542ab08e347ee9d9021e22772fbc205e0c404bbc28416ec94"
_POPULATION = "outputs/native_population_development/native_training_01"
POPULATION_FREEZE_PATH = f"{_POPULATION}/development_prediction_freeze.json"
POPULATION_FREEZE_SHA256 = "e03f6b9407b29e7e76da7a0be566ce9036c5cf27ac65bfaa9a4d3d7e009b5753"
_FIX = f"{_POPULATION}_phase_c_01_scoring_fix_v2/reusable_interface_fix_verification.json"
_REPORT = f"{_POPULATION}_phase_c_01_scoring/development_report.json"
_SUMMARY = f"{_POPULATION}_phase_c_01_scoring/scoring_summary.json"
_NATIVE = "outputs/native_training_run_02"
_HOLDOUT = "data/holdout_transfer"
_EC = "outputs/ecmodel_prior_port_01"
_VARIANTS = ("nad", "oxygen_peroxide", "oxygen_water")

HISTORICAL_SHA256 = {
    f"{_HOLDOUT}/source_index.json": "6ad84e0aacea5b861d982b930f8f0b1eb6a5278c2899e36c416ad76d3826220b",
    CHECKPOINT_PATH: CHECKPOINT_SHA256,
    f"{_NATIVE}/contract.json": "e51bd185a8082634dc6b2926041ca659936591356d0141223dbc4cff3bde9c05",
    f"{_NATIVE}/training_manifest.json": "19f7e62901817e0173516c33bd124e22feff4a5ea63ffcff043283e8b8f7a589",
    "outputs/heldout_transfer_nad/predictions.json": "e14b024cdeb43681b9ff9f2c0867c17581eae5659975094b30a5895debabafe5",
    "outputs/heldout_transfer_nad/task_manifest.json": "adb481cda9d438a5902523f4aad3485a7cf8ead5668bb0389c64b7c2cd56ead5",
    "outputs/heldout_transfer_oxygen_peroxide/predictions.json": "c63078f25008b304f3979079da5f565ecb80e0b9f615100a7d6e2e0fc6a2fbdd",
    "outputs/heldout_transfer_oxygen_peroxide/task_manifest.json": "b988afb1df6bc422ce32bc42acfd70b5ac1ca0c001553afc37aa7897f6b99dab",
    "outputs/heldout_transfer_oxygen_water/predictions.json": "aae9284d276cd8ecf0f7b728f300adf7232e060c23b415ab350f4f01221418e0",
    "outputs/heldout_transfer_oxygen_water/task_manifest.json": "1aedca7103197e9ca6601155bf3dc42fab528990c21835db84c24592169f0842",
    "outputs/heldout_score_nad/score.json": "8e2d3ab40c8135f2358ebf0f5d080864d6a9caa7e89b2790d6f36c19d82461fc",
    "outputs/heldout_score_oxygen_peroxide/score.json": "aad8443c965c20c930024291182b2e3fdf8125154b841131dbec4b94c2cf0838",
    "outputs/heldout_score_oxygen_water/score.json": "1938fe11861365c45bd312b0ee8e2ce3e99f4c1dede428b4f714341b59ec3ae0",
    "outputs/biology_readiness_01/readiness.json": "66d23ad3e6acde015d81fb59fc5f74758690d6696f1b5274e68c52c2cd19c442",
    f"{_EC}/ecmodel_isoprenoid_prior.json": "522612c5806ee439c181dcf379e99df874905ebdd1a87550b8638aae49fb671e",
    _FIX: "34f98090270d78f309d1d4bddc1f28d2f5cf780e5f1f844133512bc8a49398a0",
}
_RETAINED = {
    **HISTORICAL_SHA256,
    POPULATION_FREEZE_PATH: POPULATION_FREEZE_SHA256,
    f"{_NATIVE}/checkpoint_receipt.json": "26acdc29bacca86a34e900fe9997631e6f76927819f4695418b119290bf7aa4b",
    f"{_NATIVE}/native_hog_fit.json": "59b8d21c772adcd90a9f0d1ea34aa7c6fa2acbd0727e35f42ce1610df028cf4b",
    f"{_NATIVE}/native_exchange_fit.json": "e6d14896ab52c267b22f7de74b1e704f9a8194d5fe5c9786c7cf5ec650cc9274",
    f"{_NATIVE}/native_exchange_validation.json": "8e787ff7ae6d4b2eeff7a10ad2dba8d9a0698faf06b92d61da9b06e99e58a2a2",
    f"{_NATIVE}/native_obligations.json": "9e5748a0d91b8ab3ae4a5f83405b15601fa08cfbed0c778446a664e260f9b7d0",
    f"{_NATIVE}/native_training_observations.csv": "e90064161ea3c516101c0616bef576cd457604c8241b014145e659834fa16d1d",
    f"{_NATIVE}/read_audit.json": "7f83a82d46d229c4795db2f08be7efb8a0b32a954f4444cec407083770ef3174",
    f"{_HOLDOUT}/chemistry_evidence.json": "bb9ff21723f73bd94d6ec8e602d8ba8ef0c396a2627b35886e796461ac88a518",
    f"{_HOLDOUT}/condition_records.json": "198c1cbc7dd6e89bacc6ed4763f891bf9bbcb7f0592582da0c533e56fae3cd82",
    f"{_HOLDOUT}/condition_verification.json": "36f4771951205d6ae21e9f1356e38a950cf5ede7d35073a953f845f37045c379",
    f"{_HOLDOUT}/input_resolution.json": "96a6ddf90297d5ce7508d1efbb905591480a6e2422de78a4fa52f5c5bbda7195",
    f"{_HOLDOUT}/label_provenance.json": "9415e0e0d7c27dbe47778dc2088e9c95d74f366647dfee981a1b163314fc24dd",
    f"{_HOLDOUT}/labels_manifest.json": "80054cb3a96d1dc6ad89d9c62e4ef6b882751a67031e9cad18094773c0ac4757",
    f"{_HOLDOUT}/labels.json": "1a8e6cd6171408594e4f99389aceb9dbedbd73214fa707da36dac06480c8dc2c",
    f"{_HOLDOUT}/hypotheses/manifest.json": "df946c96374f7a1c6974564603bbc17ef3835e82b646fded73329a2b49867354",
    f"{_EC}/ecmodel_isoprenoid_kcats.csv": "6cff27ba91f1e9a4489bb587ffaa4a438eeb39710f91c765e5633d0b47118928",
    f"{_EC}/enzyme_capacity_sweep.csv": "ae8bb9d4caeea64fee769c33a146a83beef6fcaeef4e1cfc73199338c98b524c",
}
_PUBLIC_COMPANIONS = frozenset({
    "data/gem/ecYeastGEM_batch.xml.gz", "data/hog2013/model_wt.xml",
    "data/hog2013/observations.xls", "data/hog2013/methods.pdf",
    "data/hog2013/sources.json", "data/hog2013/native_training_protocol.json",
    "data/physiology/chemostatData_VanHoek1998.tsv",
    "data/native_law_v2/development_protocol.json",
    "outputs/native_population_development/executor_refresh_01/execution_configuration.json",
    f"{_POPULATION}/attempt_log.json", f"{_POPULATION}/conditioning_checkpoint.json",
    f"{_POPULATION}/nonnegative_map_audit.json", f"{_POPULATION}/training_run_summary.json",
    _REPORT, _SUMMARY,
    *(f"{_HOLDOUT}/hypotheses/{kind}_{variant}.json"
      for kind in ("conditions", "chemistry") for variant in _VARIANTS),
    *(f"{_POPULATION}/checkpoints/{family}_start_{index}.json"
      for family in ("prefix_persistence", "training_constant", "prefix_offset",
                     "prefix_affine_time", "first_order_response", "transient_response")
      for index in range(1 if family == "prefix_persistence" else 8)),
})
_STORAGE = {
    CHECKPOINT_PATH: {"/training_manifest/root": "."},
    f"{_NATIVE}/contract.json": {"/root": "."},
    f"{_NATIVE}/training_manifest.json": {"/root": "."},
    f"{_HOLDOUT}/source_index.json": {
        "/original_assets/0/local_absolute_path": "data/carotenoid/elizondo2025_steady_states.tsv",
        "/original_assets/1/local_absolute_path": "data/carotenoid/SOURCE.md",
    },
    "outputs/biology_readiness_01/readiness.json": {"/claims/0/assessment/verification_root": "."},
    f"{_EC}/ecmodel_isoprenoid_prior.json": {
        "/source_code/0/path": "scripts/ecmodel_isoprenoid_prior.py",
        "/source_code/1/path": "src/ystwin/pathway/proteome.py",
        "/source_code/2/path": "src/ystwin/pathway/enzyme_capacity.py",
        "/source_model/path": "data/gem/ecYeastGEM_batch.xml.gz",
    },
    _FIX: {"/authority_pins/authority_root": "."},
    **{f"outputs/heldout_transfer_{variant}/{name}.json": {pointer: "."}
       for variant in _VARIANTS
       for name, pointer in (("predictions", "/task_manifest/root"), ("task_manifest", "/root"))},
    **{f"outputs/heldout_score_{variant}/score.json": {"/label_manifest/root": "."}
       for variant in _VARIANTS},
}
_REFERENCE_PAIRS = (
    ("path", "sha256"), ("path", "file_sha256"), ("path", "sha256_before"),
    ("path", "download_sha256"), ("original_file_path", "sha256"),
    ("chemistry_path", "chemistry_sha256"), ("condition_packet_path", "condition_packet_sha256"),
    ("hypothesis_manifest_path", "hypothesis_manifest_sha256"),
    ("strict_chemistry_evidence_path", "strict_chemistry_evidence_sha256"),
    ("native_host_path", "native_host_sha256"), ("interface_path", "interface_sha256"),
    ("custodian_only_path", "source_sha256"),
)
_SCOPE = {
    "verification": "Externally pinned portable manifest, public content and complete scientific payload digests",
    "original_verification": "Export checked retained original byte digests; loading does not open or authenticate private originals",
    "historical_digests": "Payload digests and seals identify historical originals, never rewritten portable bytes",
    "scientific_projection": "Only the enumerated storage-location JSON pointers are replaced; every other value and key is preserved",
    "transform_value_digest_basis": "UTF-8 canonical JSON value; sorted keys, compact separators, finite numbers",
    "reference_closure": "Exported scientific artifacts plus explicit lineage-only boundary nodes; boundary bytes are not read",
    "not_verified": [
        "Private original storage values and byte-for-byte reconstruction of rewritten originals",
        "Private Granados raw data, private coverage, response exports and custody/source replays",
        "Historical source-code execution, operator blindness, chronology or independent custody attestation",
        "New fitting, selection, predictions, independent validation or biological readiness",
    ],
    "new_scientific_claim_authorized": False,
    "new_export_was_sealed_before_data": False,
}
_PRIVATE_LOCATION = re.compile(r"[/](?:Users[/]|home[/][a-z])")


class PortableEvidenceError(ValueError):
    pass


def _require(condition, message):
    if not condition:
        raise PortableEvidenceError(message)


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _digest_text(value):
    _require(isinstance(value, str) and len(value) == 64
             and set(value) <= set("0123456789abcdef"), "invalid SHA-256 identity")
    return value


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _serialized(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, "duplicate JSON field")
        result[key] = value
    return result


def _invalid_constant(value):
    raise PortableEvidenceError("nonfinite JSON number")


def _parse(data):
    value = json.loads(data, object_pairs_hook=_pairs, parse_constant=_invalid_constant)
    _canonical(value)
    return value


def _relative(value):
    _require(isinstance(value, str) and value and "\\" not in value and ":" not in value,
             "expected repository-relative identity")
    path = PurePosixPath(value)
    _require(not path.is_absolute() and ".." not in path.parts and path.as_posix() == value
             and value != ".", "unsafe repository-relative identity")
    return value


def _file(root, relative):
    path = root / _relative(relative)
    _require(path.resolve() == path and path.is_file(),
             f"missing or symlinked evidence file: {relative}")
    return path


def _read_original(root, relative):
    _require(relative in _RETAINED or relative in _PUBLIC_COMPANIONS,
             "source is outside the explicit export allowlist")
    return _file(root, relative).read_bytes()


def _pointer(value, pointer):
    try:
        for part in pointer.split("/")[1:]:
            key = part.replace("~1", "/").replace("~0", "~")
            value = value[int(key)] if isinstance(value, list) else value[key]
        return value
    except (KeyError, IndexError, TypeError, ValueError):
        raise PortableEvidenceError(f"missing or mistyped declared JSON pointer: {pointer}") from None


def _set_pointer(value, pointer, replacement):
    parent, _, key = pointer.rpartition("/")
    target = _pointer(value, parent)
    key = key.replace("~1", "/").replace("~0", "~")
    target[int(key) if isinstance(target, list) else key] = replacement


def _walk(value, pointer=""):
    yield pointer, value
    if isinstance(value, dict):
        for key, item in sorted(value.items()):
            yield from _walk(item, pointer + "/" + key.replace("~", "~0").replace("/", "~1"))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk(item, pointer + "/" + str(index))


def _public(data, origin):
    encodings = ("utf-8", "utf-16-le", "utf-16-be") if origin.endswith(".xls") else ("utf-8",)
    for encoding in encodings:
        _require(not _PRIVATE_LOCATION.search(data.decode(encoding, errors="replace")),
                 f"undeclared private storage text in {origin}")
    if origin.endswith(".gz"):
        _public(gzip.decompress(data), origin[:-3])


def _transform(origin, original, historical_root):
    payload, transforms = deepcopy(original), []
    for pointer, replacement in sorted(_STORAGE.get(origin, {}).items()):
        value = _pointer(original, pointer)
        expected = historical_root if replacement == "." else historical_root + "/" + replacement
        _require(isinstance(value, str) and value == expected,
                 f"storage schema or original location mismatch: {origin} {pointer}")
        _set_pointer(payload, pointer, replacement)
        transforms.append({
            "json_pointer": pointer, "operation": "replace_storage_location",
            "original_value_sha256": _sha(_canonical(value)), "replacement": replacement,
        })
    _public(_canonical(payload), origin)
    return payload, transforms


_PHASE_MANIFESTS = {
    CHECKPOINT_PATH: ("/training_manifest", "training_input_manifest"),
    f"{_NATIVE}/training_manifest.json": ("", "training_input_manifest"),
    **{f"outputs/heldout_transfer_{variant}/{name}.json": (pointer, "evaluation_task_manifest")
       for variant in _VARIANTS
       for name, pointer in (("predictions", "/task_manifest"), ("task_manifest", ""))},
    **{f"outputs/heldout_score_{variant}/score.json": ("/label_manifest", "outcome_labels_manifest")
       for variant in _VARIANTS},
}
_ACTIVE_ROLES = _NATIVE_ROLES | _TASK_ROLES | {"training_code", "outcome_labels"}
_SPEC_FIELDS = {"path", "sha256", "role", "source"}
_ARCHIVED_EXECUTION = {
    "population_development_executor_configuration": ("/code_artifacts/", "/code_snapshots/"),
    "development_prediction_freeze_not_final_validation": (
        "/code_artifacts/", "/environment/provenance/",
    ),
    "phase_c_reusable_interface_fix_verification_v2": (
        "/authority_pins/artifacts/", "/fixed_adapter_snapshots/", "/original_adapter_snapshots/",
    ),
}


def _phase_entries(origin, payload):
    if origin == f"{_NATIVE}/contract.json":
        _require(set(payload) == {"root", "inputs", "code"}, "invalid native training contract")
        prefix, manifest, expanded = "", payload, False
        sections = {"inputs": _NATIVE_ROLES, "code": {"training_code"}}
    elif origin == f"{_HOLDOUT}/labels_manifest.json":
        _require(set(payload) == {"labels"}, "invalid outcome label contract")
        prefix, manifest, expanded = "", payload, False
        sections = {"labels": {"outcome_labels"}}
    elif origin in _PHASE_MANIFESTS:
        prefix, kind = _PHASE_MANIFESTS[origin]
        manifest, expanded = _pointer(payload, prefix), True
        _require(isinstance(manifest, dict) and manifest.get("kind") == kind
                 and type(manifest.get("schema_version")) is int and manifest["schema_version"] == 1,
                 "invalid owning dependency phase schema")
        if kind == "training_input_manifest":
            _require(manifest.get("dependency_mode") == "explicit", "native dependencies must be explicit")
            sections = {"inputs": _NATIVE_ROLES, "code": {"training_code"}}
        else:
            sections = {"inputs": _TASK_ROLES if kind == "evaluation_task_manifest" else {"outcome_labels"}}
    else:
        return {}
    result, paths = {}, set()
    for field, roles in sections.items():
        entries = manifest.get(field)
        _require(isinstance(entries, list) and bool(entries), "missing explicit phase dependencies")
        for index, entry in enumerate(entries):
            fields = _SPEC_FIELDS | ({"resolved_path", "size_bytes"} if expanded else set())
            _require(isinstance(entry, dict) and set(entry) == fields, "invalid active dependency fields")
            path = _relative(entry["path"])
            _require(path not in paths, "duplicate active dependency identity")
            paths.add(path)
            _digest_text(entry["sha256"])
            _require(isinstance(entry["role"], str) and entry["role"] in roles,
                     "input role is not allowed by its owning phase")
            source = _source(InputSpec(**{key: entry[key] for key in _SPEC_FIELDS}), frozenset(roles))
            if expanded:
                _require(entry["resolved_path"] == path and type(entry["size_bytes"]) is int
                         and entry["size_bytes"] > 0, "invalid resolved active dependency identity")
            if entry["role"] in {"training_code", "evaluation_code"}:
                _require(source["source_id"] == path, "active code path differs from its declared source identity")
                _relative(source["reference"])
            result[f"{prefix}/{field}/{index}"] = entry
    return result


def _spec_identity(entry):
    return _canonical({key: entry[key] for key in _SPEC_FIELDS})


def _dependency_contracts(payloads, records):
    phases = {origin: _phase_entries(origin, payload) for origin, payload in payloads.items()
              if isinstance(payload, dict)}
    contract = payloads[f"{_NATIVE}/contract.json"]
    declarations = {_spec_identity(entry) for entry in phases[f"{_NATIVE}/contract.json"].values()}
    for origin in (CHECKPOINT_PATH, f"{_NATIVE}/training_manifest.json"):
        _require({_spec_identity(entry) for entry in phases[origin].values()} == declarations,
                 "active training dependencies differ from the owning contract")
    selection = {entry["path"] for entry in contract["inputs"] if entry["role"] == "selection"}
    frozen = payloads[CHECKPOINT_PATH]["training_manifest"]
    for entries in phases.values():
        for entry in entries.values():
            role, path = entry["role"], entry["path"]
            if role in {"training_code", "evaluation_code"}:
                _require(entry["source"]["reference"] in selection,
                         "active code lacks its owning protocol reference")
                if role == "evaluation_code":
                    _require_frozen_evaluation_code([entry], frozen["code"])
                else:
                    _require(_spec_identity(entry) in declarations,
                             "active training code lacks an approved contract identity")
            else:
                _require(path in records and records[path]["original_sha256"] == entry["sha256"],
                         "active data dependency lacks matching public content")
                if "size_bytes" in entry:
                    _require(entry["size_bytes"] == records[path]["original_size_bytes"],
                             "active data dependency byte count differs from its content identity")
                if role == "outcome_labels":
                    _reject_label_overlap([{**entry, "resolved_path": path}], frozen["inputs"] + frozen["code"])
    return phases


def _archival_execution_reference(origin, payload, pointer):
    kind = payload.get("document_type")
    prefixes = _ARCHIVED_EXECUTION.get(kind, ()) if isinstance(kind, str) else ()
    if origin == _SUMMARY:
        _require(payload.get("adapter_version") == "native_population_phase_c_scoring_v1"
                 and payload.get("strong_claim_authorized") is False,
                 "invalid historical population score schema")
        prefixes = ("/adapter_code/",)
    return any(pointer.startswith(prefix) for prefix in prefixes)


def _references(origin, payload):
    refs, active = [], _phase_entries(origin, payload)
    for pointer, node in _walk(payload):
        if (origin == "data/native_law_v2/development_protocol.json"
                and (pointer == "/export_contract/schemas" or pointer.startswith("/export_contract/schemas/"))):
            continue
        if not isinstance(node, dict):
            continue
        role = node.get("role")
        if isinstance(role, str) and role in _ACTIVE_ROLES:
            _require(pointer in active, "active dependency is outside an owning phase contract")
        elif role == "execution_code":
            _require(_archival_execution_reference(origin, payload, pointer),
                     "execution code is not a declared archival provenance reference")
            _relative(node.get("path"))
        for path_key, hash_key in _REFERENCE_PAIRS:
            if path_key in node and hash_key in node:
                path, digest = node[path_key], node[hash_key]
                _require(isinstance(path, str) and not PurePosixPath(path).is_absolute(),
                         f"nonportable reference in {origin} {pointer}")
                refs.append((pointer + "/" + path_key, pointer + "/" + hash_key,
                             path, _digest_text(digest)))
        if "checkpoint_sha256" in node:
            refs.append((pointer + "/checkpoint_sha256", pointer + "/checkpoint_sha256",
                         CHECKPOINT_PATH, _digest_text(node["checkpoint_sha256"])))
        if "prediction_freeze_sha256" in node:
            refs.append((pointer + "/prediction_freeze_sha256", pointer + "/prediction_freeze_sha256",
                         POPULATION_FREEZE_PATH, _digest_text(node["prediction_freeze_sha256"])))
    if origin == f"{_EC}/ecmodel_isoprenoid_prior.json":
        for name, value in sorted(payload["outputs"].items()):
            refs.append(("/outputs/" + name, "/outputs/" + name + "/sha256",
                         f"{_EC}/{name}", _digest_text(value["sha256"])))
    for variant in _VARIANTS:
        if origin == f"outputs/heldout_score_{variant}/score.json":
            refs.append(("/predictions_sha256", "/predictions_sha256",
                         f"outputs/heldout_transfer_{variant}/predictions.json",
                         _digest_text(payload["predictions_sha256"])))
    return refs


def _graph(payloads, records):
    phases = _dependency_contracts(payloads, records)
    links, external, required = [], {}, set(_RETAINED)
    for origin, payload in sorted(payloads.items()):
        if not isinstance(payload, dict):
            continue
        for pointer, digest_pointer, target, digest in _references(origin, payload):
            active = phases[origin].get(pointer.rpartition("/")[0])
            if active is not None:
                exported = active["role"] not in {"training_code", "evaluation_code"} or target in records
            else:
                exported = target in _RETAINED or target in _PUBLIC_COMPANIONS
            if exported:
                required.add(target)
                _require(target in records and records[target]["original_sha256"] == digest,
                         f"reference identity mismatch or missing companion: {origin} {pointer}")
                public_digest = records[target]["public_sha256"]
            else:
                public_digest = None
                external[target, digest] = {
                    "origin_path": target, "original_sha256": digest,
                    "integrity_scope": "lineage_only_not_exported_or_independently_verified",
                }
            links.append({
                "source_origin_path": origin, "json_pointer": pointer,
                "digest_pointer": digest_pointer, "target_origin_path": target,
                "target_original_sha256": digest, "target_public_sha256": public_digest,
            })
    _require(set(records) == required, "portable record inventory is not the declared reference closure")
    return sorted(links, key=lambda row: _canonical(row)), [external[key] for key in sorted(external)]


def _invariants(payloads):
    checkpoint = payloads[CHECKPOINT_PATH]
    _require(set(checkpoint) == {"kind", "schema_version", "model", "training_manifest"}
             and checkpoint["kind"] == "native_checkpoint" and checkpoint["schema_version"] == 1,
             "invalid native checkpoint schema")
    model = checkpoint["model"]
    _require(set(model) == {"parameters", "config"} and isinstance(model["parameters"], dict)
             and bool(model["parameters"]), "missing scientific model parameters")
    _require(all(type(value) in (int, float) for value in model["parameters"].values()),
             "invalid scientific model parameters")
    _require(checkpoint["training_manifest"] == payloads[f"{_NATIVE}/training_manifest.json"],
             "checkpoint and training manifest differ")
    labels = payloads[f"{_HOLDOUT}/labels.json"]
    observations = {row["experiment_id"]: row["value"] for row in labels["observations"]}
    for variant in _VARIANTS:
        prediction = payloads[f"outputs/heldout_transfer_{variant}/predictions.json"]
        task = payloads[f"outputs/heldout_transfer_{variant}/task_manifest.json"]
        score = payloads[f"outputs/heldout_score_{variant}/score.json"]
        _require(prediction["task_manifest"] == task, "embedded task manifest mismatch")
        _require(prediction["model_sha256"] == _sha(_canonical(model)),
                 "prediction model identity mismatch")
        for key in ("quantity", "unit"):
            _require(task[key] == score[key] == labels[key] == model["config"][key],
                     f"scientific {key} mismatch")
        predictions = {row["experiment_id"]: row for row in prediction["predictions"]}
        ids = [row["experiment_id"] for row in task["experiment_inputs"]]
        _require(set(ids) == set(predictions) == set(observations)
                 and len(ids) == len(prediction["predictions"]) == len(observations),
                 "prediction and observation inventory mismatch")
        _require(len(score["rows"]) == len(ids), "score inventory mismatch")
        for row in score["rows"]:
            expected = predictions[row["experiment_id"]]
            _require(all(row[key] == expected[key] for key in ("status", "lower", "upper", "reason"))
                     and row["observed"] == observations[row["experiment_id"]],
                     "prediction or observation differs from the historical score")
    population = payloads[POPULATION_FREEZE_PATH]
    _require(population["strong_claim_authorized"] is False, "population freeze claim changed")
    for attempt in population["attempts"]:
        checkpoint = payloads[attempt["checkpoint"]["path"]]
        _require(all(checkpoint[key] == attempt[key]
                     for key in ("model_id", "parameters", "seed", "start_index", "training_mean_mse")),
                 "population checkpoint differs from its frozen attempt")
    _require(payloads[_FIX]["prediction_freeze"]["sha256"] == POPULATION_FREEZE_SHA256,
             "population seal lineage mismatch")
    _require(payloads["outputs/biology_readiness_01/readiness.json"]["ready_for_biological_law_claim"] is False,
             "historical readiness claim changed")


def _build(root):
    initial = _read_original(root, CHECKPOINT_PATH)
    _require(_sha(initial) == CHECKPOINT_SHA256, "original checkpoint SHA-256 mismatch")
    historical_root = _parse(initial)["training_manifest"]["root"]
    _require(isinstance(historical_root, str) and PurePosixPath(historical_root).is_absolute()
             and ".." not in PurePosixPath(historical_root).parts,
             "invalid historical repository context")
    expected, originals, payloads, records, files = dict(_RETAINED), {}, {}, {}, {}
    while set(expected) - set(records):
        origin = min(set(expected) - set(records))
        data = _read_original(root, origin)
        _require(_sha(data) == expected[origin], f"original SHA-256 mismatch: {origin}")
        originals[origin] = data
        is_json = origin.endswith(".json")
        if is_json:
            payload, transforms = _transform(origin, _parse(data), historical_root)
            _require(isinstance(payload, dict), f"expected JSON evidence object: {origin}")
            content = _serialized({
                "schema_version": 1, "kind": "portable_evidence_artifact",
                "origin": {"path": origin, "sha256": expected[origin]}, "payload": payload,
            })
            science = _sha(_canonical(payload))
            for _, _, target, digest in _references(origin, payload):
                if target in _RETAINED or target in _PUBLIC_COMPANIONS:
                    _require(target not in expected or expected[target] == digest,
                             f"conflicting original identity: {target}")
                    expected[target] = digest
        else:
            payload, content, transforms, science = data, data, [], _sha(data)
        _public(content, origin)
        public_path = "artifacts/" + origin
        records[origin] = {
            "origin_path": origin, "original_sha256": expected[origin],
            "original_size_bytes": len(data), "public_path": public_path,
            "public_sha256": _sha(content), "public_size_bytes": len(content),
            "format": "json_envelope" if is_json else "original_bytes",
            "transforms": transforms,
            "scientific_payload": {
                "original_sha256": science, "public_sha256": science,
                "basis": "canonical_json_after_declared_storage_transforms" if is_json else "exact_bytes",
            },
        }
        payloads[origin], files[public_path] = payload, content
    _invariants(payloads)
    links, external = _graph(payloads, records)
    manifest = {
        "schema_version": 1, "kind": "portable_evidence_manifest", "profile": "native_v1",
        "integrity_scope": deepcopy(_SCOPE), "records": [records[key] for key in sorted(records)],
        "links": links, "external_references": external,
    }
    _public(_serialized(manifest), "manifest.json")
    return manifest, files, originals


def export_portable_evidence(root, destination) -> dict:
    root = Path(root).resolve(strict=True)
    destination = Path(destination)
    destination = destination if destination.is_absolute() else root / destination
    _require(destination.resolve() == destination, "export destination must not traverse symlinks")
    if destination.exists():
        raise FileExistsError("portable evidence destination is write-once and must not exist")
    manifest, files, originals = _build(root)
    for origin, data in originals.items():
        _require(_read_original(root, origin) == data, "original changed while preparing export")
    destination.mkdir(parents=True, exist_ok=False)
    for relative, data in sorted(files.items()):
        path = destination / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as handle:
            handle.write(data)
    for origin, data in originals.items():
        _require(_read_original(root, origin) == data, "original changed while writing export")
    with (destination / "manifest.json").open("xb") as handle:
        handle.write(_serialized(manifest))
    return manifest


@dataclass(frozen=True)
class PortableArtifact:
    payload: dict | bytes
    integrity_scope: dict


class PortableEvidenceBundle:
    def __init__(self, manifest, manifest_sha256, payloads):
        self._manifest = deepcopy(manifest)
        self._manifest_sha256 = manifest_sha256
        self._payloads = deepcopy(payloads)
        self._records = {record["origin_path"]: record for record in self._manifest["records"]}

    @property
    def manifest(self):
        return deepcopy(self._manifest)

    @property
    def manifest_sha256(self):
        return self._manifest_sha256

    def fetch(self, original_relative_path, *, original_sha256=None) -> PortableArtifact:
        origin = _relative(original_relative_path)
        _require(origin in self._records, "origin is not exported; lineage-only references are not payloads")
        record = self._records[origin]
        if original_sha256 is not None:
            _require(_digest_text(original_sha256) == record["original_sha256"],
                     "expected original identity does not match portable lineage")
        scope = {
            **deepcopy(_SCOPE), "manifest_sha256": self.manifest_sha256,
            "origin_path": origin, "original_sha256": record["original_sha256"],
            "public_sha256": record["public_sha256"], "public_content_verified": True,
            "original_bytes_verified_on_load": record["format"] == "original_bytes",
            "scientific_payload": deepcopy(record["scientific_payload"]),
            "transforms": deepcopy(record["transforms"]),
        }
        return PortableArtifact(deepcopy(self._payloads[origin]), scope)

    def fetch_public(self, public_sha256) -> PortableArtifact:
        digest = _digest_text(public_sha256)
        matches = [path for path, record in self._records.items() if record["public_sha256"] == digest]
        _require(len(matches) == 1, "public content identity is absent or ambiguous")
        return self.fetch(matches[0])


def load_portable_evidence(root, manifest_path, expected_sha256) -> PortableEvidenceBundle:
    root = Path(root).resolve(strict=True)
    path = Path(manifest_path)
    if path.is_absolute():
        _require(path.is_relative_to(root), "manifest escapes the supplied root")
        path = path.relative_to(root)
    path = _file(root, path.as_posix())
    data = path.read_bytes()
    _require(_sha(data) == _digest_text(expected_sha256), "portable manifest SHA-256 mismatch")
    _public(data, "manifest.json")
    manifest = _parse(data)
    _public(_canonical(manifest), "manifest.json")
    _require(isinstance(manifest, dict) and set(manifest) == {"schema_version", "kind", "profile", "integrity_scope", "records",
                              "links", "external_references"}
             and manifest["schema_version"] == 1 and manifest["kind"] == "portable_evidence_manifest"
             and manifest["profile"] == "native_v1" and manifest["integrity_scope"] == _SCOPE,
             "unsupported portable evidence manifest or integrity scope")
    _require(isinstance(manifest["records"], list), "invalid portable records")
    records, payloads = {}, {}
    for record in manifest["records"]:
        _require(isinstance(record, dict) and set(record) == {
            "origin_path", "original_sha256", "original_size_bytes", "public_path", "public_sha256",
            "public_size_bytes", "format", "transforms", "scientific_payload",
        }, "invalid portable record schema")
        origin = _relative(record["origin_path"])
        _require(origin not in records and (origin in _RETAINED or origin in _PUBLIC_COMPANIONS),
                 "duplicate or unapproved portable origin")
        _digest_text(record["original_sha256"])
        if origin in _RETAINED:
            _require(record["original_sha256"] == _RETAINED[origin], "retained original identity changed")
        _require(record["public_path"] == "artifacts/" + origin, "public storage link changed")
        content = _file(path.parent, record["public_path"]).read_bytes()
        _require(_sha(content) == _digest_text(record["public_sha256"])
                 and len(content) == record["public_size_bytes"], "public artifact SHA-256 mismatch")
        _public(content, origin)
        _require(type(record["original_size_bytes"]) is int and record["original_size_bytes"] > 0,
                 "invalid original byte count")
        if origin.endswith(".json"):
            _require(record["format"] == "json_envelope", "JSON must use a distinct portable envelope")
            envelope = _parse(content)
            _public(_canonical(envelope), origin)
            _require(isinstance(envelope, dict) and set(envelope) == {"schema_version", "kind", "origin", "payload"}
                     and envelope["schema_version"] == 1 and envelope["kind"] == "portable_evidence_artifact"
                     and envelope["origin"] == {"path": origin, "sha256": record["original_sha256"]},
                     "portable envelope identity mismatch")
            payload = envelope["payload"]
            _require(isinstance(payload, dict), "missing scientific payload")
            expected = _STORAGE.get(origin, {})
            _require(isinstance(record["transforms"], list) and len(record["transforms"]) == len(expected),
                     "undeclared or missing storage transforms")
            seen = set()
            for transform in record["transforms"]:
                _require(set(transform) == {"json_pointer", "operation", "original_value_sha256", "replacement"},
                         "invalid transform record")
                pointer = transform["json_pointer"]
                _require(pointer in expected and pointer not in seen
                         and transform["operation"] == "replace_storage_location"
                         and transform["replacement"] == expected[pointer]
                         and _pointer(payload, pointer) == expected[pointer],
                         "undeclared scientific redaction or storage transform")
                _digest_text(transform["original_value_sha256"])
                seen.add(pointer)
            science, basis = _sha(_canonical(payload)), "canonical_json_after_declared_storage_transforms"
        else:
            _require(record["format"] == "original_bytes" and record["transforms"] == []
                     and record["original_sha256"] == _sha(content)
                     and record["original_size_bytes"] == len(content), "original-byte companion changed")
            payload, science, basis = content, _sha(content), "exact_bytes"
        _require(record["scientific_payload"] == {
            "original_sha256": science, "public_sha256": science, "basis": basis,
        }, "scientific payload digest mismatch")
        records[origin], payloads[origin] = record, payload
    _require(set(_RETAINED) <= set(records), "missing retained evidence")
    links, external = _graph(payloads, records)
    _require(manifest["links"] == links and manifest["external_references"] == external,
             "portable reference graph changed")
    _invariants(payloads)
    return PortableEvidenceBundle(manifest, expected_sha256, payloads)
