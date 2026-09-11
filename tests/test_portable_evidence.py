from copy import deepcopy
from graphlib import TopologicalSorter
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import shutil
from types import SimpleNamespace

import pytest

from ystwin.analysis import portable_evidence as portable


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "data/frozen_evidence/native_v1"
MANIFEST_SHA256 = "eff317f689fe096f6b3df35e050f557b51cd7e51cdb192bb95775eef8ef529ec"
CHECKPOINT = portable.CHECKPOINT_PATH
PREDICTIONS = "outputs/heldout_transfer_nad/predictions.json"
SCORE = "outputs/heldout_score_nad/score.json"
FIXTURE_PINSET_SHA256 = "ced333923d4fc65253743409d9147c4138b2ab703bfd42317c1d29d6a9e16ae8"
FIXTURE_MANIFEST_SHA256 = "a525bb53d957a18474d5354fafde903705d32be6a358e88323330dc253ca2440"


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _at(value, pointer):
    for part in pointer.split("/")[1:]:
        part = part.replace("~1", "/").replace("~0", "~")
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value


def _replace(value, pointer, replacement):
    parent, _, key = pointer.rpartition("/")
    container = _at(value, parent)
    key = key.replace("~1", "/").replace("~0", "~")
    container[int(key) if isinstance(container, list) else key] = replacement


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    destination = tmp_path_factory.mktemp("portable_evidence") / "export"
    shutil.copytree(FROZEN, destination)
    data = (destination / "manifest.json").read_bytes()
    assert _digest(data) == MANIFEST_SHA256
    return destination, json.loads(data), _digest(data)


@pytest.fixture(scope="module", autouse=True)
def preserve_available_historical_originals(exported, tmp_path_factory):
    originals = {}
    for record in exported[1]["records"]:
        path = ROOT / record["origin_path"]
        if path.is_file():
            assert path.resolve() == path
            data = path.read_bytes()
            assert _digest(data) == record["original_sha256"]
            originals[path] = data
    if len(originals) == len(exported[1]["records"]):
        destination = tmp_path_factory.mktemp("historical_export") / "bundle"
        assert portable.export_portable_evidence(ROOT, destination) == exported[1]
        portable.load_portable_evidence(destination, "manifest.json", MANIFEST_SHA256)
    yield
    for path, data in originals.items():
        assert path.read_bytes() == data


def _synthetic_sources(bundle):
    manifest = bundle.manifest
    records = {record["origin_path"]: record for record in manifest["records"]}
    dependencies = {origin: set() for origin in records}
    links = {origin: [] for origin in records}
    for link in manifest["links"]:
        if link["target_origin_path"] in records:
            dependencies[link["source_origin_path"]].add(link["target_origin_path"])
            links[link["source_origin_path"]].append(link)
    embedded = {CHECKPOINT: ("training_manifest", "outputs/native_training_run_02/training_manifest.json")}
    for variant in ("nad", "oxygen_peroxide", "oxygen_water"):
        base = f"outputs/heldout_transfer_{variant}"
        embedded[f"{base}/predictions.json"] = ("task_manifest", f"{base}/task_manifest.json")
    for origin, (_, target) in embedded.items():
        dependencies[origin].add(target)
    historical_root = str(PurePosixPath("/") / "Users" / "portable-fixture" / "repository")
    data, pins, payloads = {}, {}, {}
    for origin in TopologicalSorter(dependencies).static_order():
        record = records[origin]
        payload = bundle.fetch(origin).payload
        if record["format"] == "original_bytes":
            data[origin] = payload
        else:
            for transform in record["transforms"]:
                replacement = transform["replacement"]
                value = historical_root if replacement == "." else historical_root + "/" + replacement
                _replace(payload, transform["json_pointer"], value)
            for link in links[origin]:
                target = link["target_origin_path"]
                _replace(payload, link["digest_pointer"], pins[target])
                reference = _at(payload, link["digest_pointer"].rpartition("/")[0])
                if isinstance(reference, dict) and reference.get("path") == target:
                    for field in ("size_bytes", "bytes"):
                        if field in reference:
                            reference[field] = len(data[target])
            if origin in embedded:
                field, target = embedded[origin]
                payload[field] = deepcopy(payloads[target])
            for _, node in reversed(list(portable._walk(payload))):
                if (isinstance(node, dict) and node.get("kind") in (
                        "training_input_manifest", "evaluation_task_manifest",
                        "outcome_labels_manifest", "interval_prediction_score")):
                    node["sha256"] = _digest(_canonical({key: value for key, value in node.items() if key != "sha256"}))
            payloads[origin] = payload
            data[origin] = _json_bytes(payload)
        pins[origin] = _digest(data[origin])
    assert payloads[CHECKPOINT]["model"] == bundle.fetch(CHECKPOINT).payload["model"]
    assert payloads[portable.POPULATION_FREEZE_PATH]["predictions"] == bundle.fetch(portable.POPULATION_FREEZE_PATH).payload["predictions"]
    return data, pins


@pytest.fixture(scope="module")
def exporter_material(exported):
    data, pins = _synthetic_sources(_load(exported))
    assert _digest(_canonical(pins)) == FIXTURE_PINSET_SHA256
    return data, pins


@pytest.fixture
def exporter_source(exporter_material, tmp_path, monkeypatch):
    data, pins = exporter_material
    root = tmp_path / "synthetic_originals"
    for origin, raw in data.items():
        path = root / origin
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    monkeypatch.setattr(portable, "_RETAINED", {origin: pins[origin] for origin in portable._RETAINED})
    monkeypatch.setattr(portable, "CHECKPOINT_SHA256", pins[CHECKPOINT])
    monkeypatch.setattr(portable, "POPULATION_FREEZE_SHA256", pins[portable.POPULATION_FREEZE_PATH])
    return SimpleNamespace(root=root, pins=pins, originals=data)


@pytest.fixture
def mutable_export(exported, tmp_path):
    source, manifest, digest = exported
    destination = tmp_path / "portable"
    shutil.copytree(source, destination)
    return destination, deepcopy(manifest), digest


def _load(exported):
    destination, _, digest = exported
    return portable.load_portable_evidence(destination, "manifest.json", digest)


def _record(manifest, origin):
    return next(record for record in manifest["records"] if record["origin_path"] == origin)


def _tamper(exported, origin, pointer, replacement, *, refresh_content=False, refresh_science=False):
    destination, manifest, digest = exported
    record = _record(manifest, origin)
    path = destination / record["public_path"]
    envelope = json.loads(path.read_bytes())
    _replace(envelope["payload"], pointer, replacement)
    data = _json_bytes(envelope)
    path.write_bytes(data)
    if refresh_content:
        record["public_sha256"] = _digest(data)
        record["public_size_bytes"] = len(data)
        for link in manifest["links"]:
            if link["target_origin_path"] == origin:
                link["target_public_sha256"] = record["public_sha256"]
        if refresh_science:
            science = _digest(_canonical(envelope["payload"]))
            record["scientific_payload"]["original_sha256"] = science
            record["scientific_payload"]["public_sha256"] = science
        data = _json_bytes(manifest)
        (destination / "manifest.json").write_bytes(data)
        digest = _digest(data)
    return destination, manifest, digest


def test_shipped_native_v1_loads_with_retained_manifest_pin(monkeypatch):
    def forbidden(*args):
        raise AssertionError("shipped public evidence must not open originals")

    monkeypatch.setattr(portable, "_read_original", forbidden)
    bundle = portable.load_portable_evidence(ROOT, "data/frozen_evidence/native_v1/manifest.json", MANIFEST_SHA256)
    assert bundle.manifest_sha256 == MANIFEST_SHA256
    assert bundle.fetch(CHECKPOINT).integrity_scope["original_sha256"] == portable.CHECKPOINT_SHA256
    assert bundle.fetch(portable.POPULATION_FREEZE_PATH).integrity_scope["original_sha256"] == portable.POPULATION_FREEZE_SHA256


def test_public_bundle_has_distinct_lineage_and_content_identities(exported):
    bundle = _load(exported)
    artifact = bundle.fetch(CHECKPOINT, original_sha256=portable.CHECKPOINT_SHA256)
    record = _record(bundle.manifest, CHECKPOINT)
    assert len(bundle.manifest["records"]) == 96
    assert sum(len(item["transforms"]) for item in bundle.manifest["records"]) == 20
    assert len(portable.HISTORICAL_SHA256) == 16
    assert record["public_sha256"] != portable.CHECKPOINT_SHA256
    assert artifact.payload["training_manifest"]["root"] == "."
    assert artifact.integrity_scope["original_bytes_verified_on_load"] is False
    assert artifact.integrity_scope["new_export_was_sealed_before_data"] is False
    assert artifact.integrity_scope["new_scientific_claim_authorized"] is False
    assert artifact.integrity_scope["public_content_verified"] is True
    assert bundle.fetch_public(record["public_sha256"]).payload == artifact.payload
    with pytest.raises(portable.PortableEvidenceError, match="original identity"):
        bundle.fetch(CHECKPOINT, original_sha256=record["public_sha256"])
    with pytest.raises(portable.PortableEvidenceError, match="absent or ambiguous"):
        bundle.fetch_public(portable.CHECKPOINT_SHA256)
    wrapper = json.loads((exported[0] / record["public_path"]).read_bytes())
    assert wrapper["kind"] == "portable_evidence_artifact"
    assert wrapper["kind"] != artifact.payload["kind"]
    artifact.payload["model"]["parameters"].clear()
    assert bundle.fetch(CHECKPOINT).payload["model"]["parameters"]
    manifest = bundle.manifest
    manifest["records"].clear()
    assert bundle.manifest["records"]


def test_all_original_bytes_and_scientific_payloads_are_unchanged(exporter_source, tmp_path):
    source = exporter_source
    originals = {origin: (source.root / origin).read_bytes() for origin in source.pins}
    destination = tmp_path / "unchanged"
    manifest = portable.export_portable_evidence(source.root, destination)
    assert _digest((destination / "manifest.json").read_bytes()) == FIXTURE_MANIFEST_SHA256
    bundle = portable.load_portable_evidence(destination, "manifest.json", FIXTURE_MANIFEST_SHA256)
    assert len(manifest["records"]) == 96
    assert sum(len(record["transforms"]) for record in manifest["records"]) == 20
    assert _digest(originals[CHECKPOINT]) == source.pins[CHECKPOINT]
    assert _digest(originals[portable.POPULATION_FREEZE_PATH]) == source.pins[portable.POPULATION_FREEZE_PATH]
    for record in manifest["records"]:
        origin = record["origin_path"]
        assert (source.root / origin).read_bytes() == originals[origin]
        assert _digest(originals[origin]) == record["original_sha256"]
        artifact = bundle.fetch(origin)
        if record["format"] == "original_bytes":
            assert artifact.payload == originals[origin]
            assert artifact.integrity_scope["original_bytes_verified_on_load"] is True
            continue
        expected = json.loads(originals[origin])
        for transform in record["transforms"]:
            pointer = transform["json_pointer"]
            original_value = _at(expected, pointer)
            assert isinstance(original_value, str)
            assert _digest(_canonical(original_value)) == transform["original_value_sha256"]
            assert original_value != transform["replacement"]
            _replace(expected, pointer, transform["replacement"])
        assert artifact.payload == expected
        assert record["scientific_payload"]["original_sha256"] == _digest(_canonical(expected))
        assert record["scientific_payload"]["public_sha256"] == _digest(_canonical(artifact.payload))


def test_export_is_deterministic_after_input_and_output_relocation(exporter_source, tmp_path):
    source = exporter_source
    first = tmp_path / "first_export"
    expected = portable.export_portable_evidence(source.root, first)
    assert _digest((first / "manifest.json").read_bytes()) == FIXTURE_MANIFEST_SHA256
    relocated = tmp_path / "relocated_inputs"
    relocated.mkdir()
    for origin in source.pins:
        path = relocated / origin
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((source.root / origin).read_bytes())
    destination = tmp_path / "second_export"
    manifest = portable.export_portable_evidence(relocated, destination)
    assert manifest == expected
    assert (destination / "manifest.json").read_bytes() == (first / "manifest.json").read_bytes()
    for record in manifest["records"]:
        relative = record["public_path"]
        assert (destination / relative).read_bytes() == (first / relative).read_bytes()
    portable.load_portable_evidence(destination, "manifest.json", FIXTURE_MANIFEST_SHA256)


def test_portable_only_load_never_opens_originals_or_private_lineage(exported, tmp_path, monkeypatch):
    destination = tmp_path / "isolated"
    shutil.copytree(exported[0], destination)
    original_read = Path.read_bytes
    accessed = []

    def guarded_read(path):
        assert path.is_relative_to(destination)
        accessed.append(path.relative_to(destination).as_posix())
        return original_read(path)

    def forbidden_original(*args):
        raise AssertionError("portable loading must not consult original files")

    monkeypatch.setattr(Path, "read_bytes", guarded_read)
    monkeypatch.setattr(portable, "_read_original", forbidden_original)
    bundle = portable.load_portable_evidence(destination, "manifest.json", exported[2])
    assert len(accessed) == len(bundle.manifest["records"]) + 1
    assert not (destination / CHECKPOINT).exists()
    assert bundle.fetch(PREDICTIONS).payload["predictions"][1]["upper"] == 16.574032086242326
    assert bundle.fetch(SCORE).payload["rows"][0]["observed"] == 0.5025304583095278
    assert bundle.fetch(portable.POPULATION_FREEZE_PATH).integrity_scope["original_sha256"] == portable.POPULATION_FREEZE_SHA256
    assert bundle.fetch("outputs/biology_readiness_01/readiness.json").payload["ready_for_biological_law_claim"] is False


def test_reference_closure_resolves_origin_and_new_content_without_substitution(exported):
    bundle = _load(exported)
    records = {record["origin_path"]: record for record in bundle.manifest["records"]}
    external = {(ref["origin_path"], ref["original_sha256"])
                for ref in bundle.manifest["external_references"]}
    assert external
    for link in bundle.manifest["links"]:
        origin, old_digest = link["target_origin_path"], link["target_original_sha256"]
        if link["target_public_sha256"] is None:
            assert (origin, old_digest) in external
            assert origin not in records
        else:
            assert records[origin]["original_sha256"] == old_digest
            old = bundle.fetch(origin, original_sha256=old_digest)
            new = bundle.fetch_public(link["target_public_sha256"])
            assert old.payload == new.payload
    with pytest.raises(portable.PortableEvidenceError, match="lineage-only"):
        bundle.fetch("src/ystwin/paths.py")
    assert not any(path.startswith("data/native_law_v2/granados/") for path in records)
    payloads = {origin: bundle.fetch(origin).payload for origin in records}
    training_path = "outputs/native_training_run_02/training_manifest.json"
    changed = payloads[training_path]["code"][0]
    changed["path"] = changed["resolved_path"] = changed["source"]["source_id"] = "unexported/active_input.py"
    payloads[CHECKPOINT]["training_manifest"] = deepcopy(payloads[training_path])
    with pytest.raises(portable.PortableEvidenceError, match="owning contract"):
        portable._dependency_contracts(payloads, records)
    with pytest.raises(portable.PortableEvidenceError, match="owning phase"):
        portable._references("citation.json", {"dependency": changed})


def test_export_reads_only_explicit_nonprivate_closure(exporter_source, tmp_path, monkeypatch):
    source = exporter_source
    allowed, actual = set(source.pins), set()
    read = portable._read_original

    def guarded(root, relative):
        assert root == source.root
        assert relative in allowed
        actual.add(relative)
        return read(root, relative)

    monkeypatch.setattr(portable, "_read_original", guarded)
    destination = tmp_path / "allowlisted"
    portable.export_portable_evidence(source.root, destination)
    assert _digest((destination / "manifest.json").read_bytes()) == FIXTURE_MANIFEST_SHA256
    assert actual == allowed


def test_write_once_does_not_reopen_sources(exported, monkeypatch):
    def forbidden(*args):
        raise AssertionError("an existing destination must be rejected before reading originals")

    monkeypatch.setattr(portable, "_read_original", forbidden)
    with pytest.raises(FileExistsError, match="write-once"):
        portable.export_portable_evidence(ROOT, exported[0])


def test_original_pin_failure_precedes_all_output_writes(exporter_source, tmp_path, monkeypatch):
    source = exporter_source
    read = portable._read_original

    def altered(root, relative):
        data = read(root, relative)
        if relative == CHECKPOINT:
            payload = json.loads(data)
            parameter = next(iter(payload["model"]["parameters"]))
            payload["model"]["parameters"][parameter] += 1.0
            return _json_bytes(payload)
        return data

    monkeypatch.setattr(portable, "_read_original", altered)
    destination = tmp_path / "rejected"
    with pytest.raises(portable.PortableEvidenceError, match="original checkpoint SHA-256"):
        portable.export_portable_evidence(source.root, destination)
    assert not destination.exists()
    assert (source.root / CHECKPOINT).read_bytes() == source.originals[CHECKPOINT]


@pytest.mark.parametrize("origin,pointer,value", [
    (CHECKPOINT, "/model/parameters", {"altered_parameter": 1.0}),
    (PREDICTIONS, "/predictions/1/upper", 999.0),
    (SCORE, "/unit", "different_unit"),
    (PREDICTIONS, "/checkpoint_sha256", "0" * 64),
    (SCORE, "/rows/0/observed", 0.0),
])
def test_tampered_scientific_bytes_are_rejected(mutable_export, origin, pointer, value):
    changed = _tamper(mutable_export, origin, pointer, value)
    with pytest.raises(portable.PortableEvidenceError, match="public artifact SHA-256"):
        _load(changed)


@pytest.mark.parametrize("origin,pointer,value", [
    (CHECKPOINT, "/model/parameters", {"altered_parameter": 1.0}),
    (PREDICTIONS, "/predictions/1/upper", 999.0),
    (SCORE, "/unit", "different_unit"),
])
def test_rehashed_public_content_still_requires_scientific_digest(mutable_export, origin, pointer, value):
    changed = _tamper(mutable_export, origin, pointer, value, refresh_content=True)
    with pytest.raises(portable.PortableEvidenceError, match="scientific payload digest"):
        _load(changed)


@pytest.mark.parametrize("origin,pointer,value,reason", [
    (CHECKPOINT, "/model/parameters", {"altered_parameter": 1.0}, "prediction model identity"),
    (PREDICTIONS, "/predictions/1/upper", 999.0, "historical score"),
    (SCORE, "/unit", "different_unit", "scientific unit"),
    (PREDICTIONS, "/checkpoint_sha256", "0" * 64, "reference identity"),
])
def test_cross_artifact_links_reject_rehashed_inconsistent_payloads(mutable_export, origin, pointer, value, reason):
    changed = _tamper(mutable_export, origin, pointer, value, refresh_content=True, refresh_science=True)
    with pytest.raises(portable.PortableEvidenceError, match=reason):
        _load(changed)


def test_manifest_requires_an_externally_retained_digest(mutable_export):
    destination, manifest, digest = mutable_export
    manifest["links"] = []
    (destination / "manifest.json").write_bytes(_json_bytes(manifest))
    with pytest.raises(portable.PortableEvidenceError, match="manifest SHA-256"):
        portable.load_portable_evidence(destination, "manifest.json", digest)
    changed_digest = _digest((destination / "manifest.json").read_bytes())
    with pytest.raises(portable.PortableEvidenceError, match="reference graph changed"):
        portable.load_portable_evidence(destination, "manifest.json", changed_digest)


def test_scientific_fields_cannot_be_declared_storage_redactions(mutable_export):
    destination, manifest, _ = mutable_export
    record = _record(manifest, CHECKPOINT)
    record["transforms"][0]["json_pointer"] = "/model/config/unit"
    data = _json_bytes(manifest)
    (destination / "manifest.json").write_bytes(data)
    with pytest.raises(portable.PortableEvidenceError, match="scientific redaction"):
        portable.load_portable_evidence(destination, "manifest.json", _digest(data))


def test_missing_scientific_companion_fails_closed(mutable_export):
    destination, manifest, digest = mutable_export
    record = _record(manifest, "data/holdout_transfer/labels.json")
    (destination / record["public_path"]).unlink()
    with pytest.raises(portable.PortableEvidenceError, match="missing or symlinked"):
        portable.load_portable_evidence(destination, "manifest.json", digest)


def test_public_paths_cannot_escape_or_follow_symlinks(mutable_export, tmp_path):
    destination, manifest, _ = mutable_export
    record = _record(manifest, CHECKPOINT)
    original = destination / record["public_path"]
    outside = tmp_path / "outside.json"
    outside.write_bytes(original.read_bytes())
    original.unlink()
    original.symlink_to(outside)
    with pytest.raises(portable.PortableEvidenceError, match="symlinked"):
        _load(mutable_export)
    record["public_path"] = "../outside.json"
    data = _json_bytes(manifest)
    (destination / "manifest.json").write_bytes(data)
    with pytest.raises(portable.PortableEvidenceError, match="public storage link"):
        portable.load_portable_evidence(destination, "manifest.json", _digest(data))


def test_undeclared_private_text_is_refused_not_rewritten(exported):
    checkpoint = _load(exported).fetch(CHECKPOINT).payload
    synthetic_root = str(PurePosixPath("/") / "Users" / "synthetic" / "repository")
    checkpoint["training_manifest"]["root"] = synthetic_root
    checkpoint["model"]["config"]["unit"] = synthetic_root + "/measurement"
    with pytest.raises(portable.PortableEvidenceError, match="undeclared private storage"):
        portable._transform(CHECKPOINT, checkpoint, synthetic_root)
    checkpoint["model"]["config"]["unit"] = "mg/gDW"
    del checkpoint["training_manifest"]["root"]
    with pytest.raises(portable.PortableEvidenceError, match="missing or mistyped"):
        portable._transform(CHECKPOINT, checkpoint, synthetic_root)


def test_escaped_private_json_values_are_also_rejected(mutable_export):
    destination, manifest, _ = mutable_export
    private = str(PurePosixPath("/") / "Users" / "synthetic" / "repository")
    record = _record(manifest, SCORE)
    path = destination / record["public_path"]
    envelope = json.loads(path.read_bytes())
    envelope["payload"]["unit"] = private
    data = _json_bytes(envelope).replace(private.encode(), private.replace("/", "\\u002f").encode())
    path.write_bytes(data)
    record["public_sha256"], record["public_size_bytes"] = _digest(data), len(data)
    data = _json_bytes(manifest)
    (destination / "manifest.json").write_bytes(data)
    with pytest.raises(portable.PortableEvidenceError, match="undeclared private storage"):
        portable.load_portable_evidence(destination, "manifest.json", _digest(data))


@pytest.mark.parametrize("data", [b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":1e999}'])
def test_json_integrity_rejects_duplicate_or_nonfinite_values(data):
    with pytest.raises(ValueError):
        portable._parse(data)


def test_new_tracked_candidates_contain_no_literal_home_paths(exported):
    pattern = re.compile(r"[/](?:Users[/]|home[/][a-z])")
    files = [ROOT / "src/ystwin/analysis/portable_evidence.py",
             ROOT / "scripts/export_portable_evidence.py", Path(__file__),
             exported[0] / "manifest.json"]
    files.extend(exported[0] / record["public_path"] for record in exported[1]["records"])
    for path in files:
        data = path.read_bytes()
        assert pattern.search(data.decode("utf-8", errors="replace")) is None
        portable._public(data, path.name)


def test_cli_reports_new_identity_without_machine_location(exporter_source, tmp_path, capsys):
    source = exporter_source
    spec = importlib.util.spec_from_file_location("portable_export_cli", ROOT / "scripts/export_portable_evidence.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    destination = tmp_path / "cli_export"
    assert module.main(["--root", str(source.root), "--destination", str(destination)]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["sha256"] == FIXTURE_MANIFEST_SHA256
    assert result["records"] == 96
    assert result["storage_transforms"] == 20
    assert result["new_scientific_claim_authorized"] is False
    assert str(ROOT) not in json.dumps(result)
    assert str(source.root) not in json.dumps(result)
    portable.load_portable_evidence(destination, "manifest.json", FIXTURE_MANIFEST_SHA256)
