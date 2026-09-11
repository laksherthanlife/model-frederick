from __future__ import annotations

import copy
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from ystwin import artifacts


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("_provenance_audit", ROOT / "scripts/audit_reproducibility.py")
repro = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = repro
SPEC.loader.exec_module(repro)


def _write_registry(root, document):
    (root / "data/artifact_registry.json").write_text(json.dumps(document), encoding="utf-8")


def _snapshot(root):
    return {path.relative_to(root).as_posix(): (path.read_bytes(), path.stat().st_mtime_ns)
            for path in root.rglob("*") if path.is_file()}


@pytest.fixture
def specimen(tmp_path):
    root = tmp_path / "repository"
    for directory in ("scripts", "src/ystwin", "data/a", "data/b", "outputs/a", "outputs/b"):
        (root / directory).mkdir(parents=True)
    (root / "src/ystwin/__init__.py").write_text("")
    (root / "src/ystwin/writer.py").write_text(
        "from pathlib import Path\n"
        "def write(destination):\n"
        "    Path(destination).write_text(Path('data/a/source.csv').read_text())\n")
    (root / "scripts/run.py").write_text(
        "import sys\nfrom ystwin.writer import write\nwrite(sys.argv[1])\n")
    payload = "id,value,unit\na,1,mg\nb,,mg\n"
    for path in ("data/a/source.csv", "data/b/source.csv", "outputs/a/result.csv"):
        (root / path).write_text(payload)
    run = {"id": "a", "lifecycle": "current", "status": "investigation_pending",
           "reason": "Synthetic current output with an explicit execution recipe",
           "producer": {"path": "scripts/run.py", "callable": "write", "libraries": ["src/ystwin/writer.py"]},
           "dependencies": [{"path": "data/a/source.csv", "role": "input", "sha256": None}],
           "parameters": {"selection": "all rows"},
           "model_identity": {"kind": "none", "reason": "identity transformation"},
           "runtime": {"kind": "python_standard_library"},
           "reproduction": {"argv": ["{python}", "-B", "scripts/run.py", "outputs/a/result.csv"]},
           "artifacts": [{"path": "outputs/a/result.csv", "comparison": {
               "format": "csv", "row_keys": ["id"], "units": {"value": "mg"}, "unit_columns": ["unit"]}}]}
    document = {"schema_version": 1, "runs": [run], "frozen_bundles": []}
    _write_registry(root, document)
    return root, document


def _verify(specimen, destination):
    root, document = specimen
    receipt = artifacts.reproduce_run(root, "a", destination)
    assert receipt["success"], receipt
    evidence = root / "data/receipt.json"
    evidence.write_text(json.dumps(receipt))
    document["runs"][0]["status"] = "verified"
    document["runs"][0]["verification"] = {"path": "data/receipt.json", "sha256": artifacts.sha256(evidence)}
    _write_registry(root, document)
    return receipt


@pytest.mark.parametrize("identity", ["result.json", "/tmp/result.json", "outputs/../result.json",
                                      "outputs//result.json", "outputs/./result.json", "outputs\\result.json",
                                      "outputs/*.json", "outputs/run/", "C:/outputs/result.json"])
def test_noncanonical_or_basename_identity_is_rejected(identity):
    with pytest.raises(artifacts.ArtifactError):
        artifacts.artifact_id(identity)


def test_nested_duplicate_basenames_are_different_identities(specimen):
    root, document = specimen
    second = copy.deepcopy(document["runs"][0])
    second["id"] = "b"
    second["producer"]["path"] = "scripts/other.py"
    second["artifacts"][0]["path"] = "outputs/b/result.csv"
    document["runs"].append(second)
    _write_registry(root, document)
    registry = artifacts.load_registry(root)
    assert registry["artifacts"]["outputs/a/result.csv"]["run"] == "a"
    assert registry["artifacts"]["outputs/b/result.csv"]["run"] == "b"


def test_same_full_identity_cannot_be_registered_twice(specimen):
    root, document = specimen
    second = copy.deepcopy(document["runs"][0])
    second["id"] = "b"
    document["runs"].append(second)
    _write_registry(root, document)
    with pytest.raises(artifacts.ArtifactError, match="registered twice"):
        artifacts.load_registry(root)


def test_duplicate_json_identity_is_not_silently_overwritten(specimen):
    root, _ = specimen
    (root / "data/artifact_registry.json").write_text('{"schema_version":1,"runs":[],"runs":[]}')
    rows = artifacts.audit_registry(root)
    assert rows[0]["status"] == "FAIL"
    assert "duplicate JSON key" in rows[0]["detail"]


def test_shipped_registry_names_every_tracked_retained_output():
    tracked = subprocess.run(["git", "ls-files", "-z", "--", "outputs"], cwd=ROOT, capture_output=True, text=True)
    if tracked.returncode:
        pytest.skip("a Git index is required for retained-output coverage")
    registry = artifacts.load_registry(ROOT)
    identities = {path for path in tracked.stdout.split("\0") if path}
    assert identities <= set(registry["artifacts"])


def test_dependency_digest_includes_full_paths(specimen):
    root, _ = specimen
    first = artifacts.producer_digest(root, [root / "data/a/source.csv"])
    second = artifacts.producer_digest(root, [root / "data/b/source.csv"])
    assert first != second
    assert len(first) == 64


def test_registry_not_filename_specificity_controls_the_producer(specimen):
    root, _ = specimen
    wrong_map = {"bogus.py": ["outputs/a/result.csv"]}
    rows = repro.check_output_provenance(root, wrong_map)
    assert "scripts/run.py::write" in rows[0]["detail"]
    assert "bogus.py" not in rows[0]["detail"]


def test_unknown_artifact_is_a_failure_not_automatically_local_or_historical(specimen):
    root, _ = specimen
    (root / "outputs/b/result.csv").write_text("id,value\na,2\n")
    row = next(row for row in artifacts.audit_registry(root) if row["check"].endswith("outputs/b/result.csv"))
    assert row["status"] == "FAIL"
    assert row["actual"] == "investigation_pending"


def test_absent_registry_and_legacy_stamp_do_not_certify_outputs(specimen):
    root, _ = specimen
    (root / "data/artifact_registry.json").rename(root / "data/not-the-registry.json")
    (root / "outputs/provenance.json").write_text('{"result.csv":"0123456789abcdef"}')
    before = _snapshot(root)
    rows = repro.check_output_provenance(root, {"run.py": ["result.csv"]})
    assert rows[0]["status"] == "FAIL"
    assert _snapshot(root) == before


def test_pending_then_reproduced_artifact_goes_from_red_to_green(specimen, tmp_path):
    root, _ = specimen
    assert artifacts.audit_registry(root)[0]["status"] == "FAIL"
    _verify(specimen, tmp_path / "reproduction")
    assert [row["status"] for row in artifacts.audit_registry(root)] == ["PASS"]


def test_library_and_input_attribution_are_both_bound(specimen, tmp_path):
    root, _ = specimen
    receipt = _verify(specimen, tmp_path / "reproduction")
    assert "src/ystwin/writer.py" in {entry["path"] for entry in receipt["identity"]["code"]}
    assert receipt["identity"]["dependencies"][0]["path"] == "data/a/source.csv"
    (root / "data/b/source.csv").write_text("unrelated input\n")
    assert artifacts.audit_registry(root)[0]["status"] == "PASS"
    (root / "data/a/source.csv").write_text("id,value,unit\na,3,mg\nb,,mg\n")
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "numerical effect is untested" in row["detail"]


def test_source_change_does_not_automatically_mean_wrong_numbers(specimen, tmp_path):
    root, _ = specimen
    _verify(specimen, tmp_path / "reproduction")
    (root / "src/ystwin/writer.py").write_text((root / "src/ystwin/writer.py").read_text() + "\n")
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "numerical effect is untested" in row["detail"]
    assert "wrong" not in row["detail"]


def test_a_stale_receipt_cannot_certify_changed_output(specimen, tmp_path):
    root, _ = specimen
    _verify(specimen, tmp_path / "reproduction")
    (root / "outputs/a/result.csv").write_text("id,value,unit\na,4,mg\nb,,mg\n")
    assert "retained artifact changed" in artifacts.audit_registry(root)[0]["detail"]


def test_declaring_verified_without_reproducing_does_not_pass(specimen):
    root, document = specimen
    document["runs"][0]["status"] = "verified"
    _write_registry(root, document)
    assert "receipt absent" in artifacts.audit_registry(root)[0]["detail"]


def test_a_failed_execution_cannot_be_verified(specimen, tmp_path):
    root, _ = specimen
    (root / "scripts/run.py").write_text((root / "scripts/run.py").read_text() + "raise SystemExit(9)\n")
    before = _snapshot(root)
    receipt = artifacts.reproduce_run(root, "a", tmp_path / "reproduction")
    assert receipt["exit_code"] == 9
    assert not receipt["success"]
    assert receipt["comparisons"]["outputs/a/result.csv"]["matches"]
    assert _snapshot(root) == before


def test_success_exit_without_output_is_not_successful_reproduction(specimen, tmp_path):
    root, _ = specimen
    (root / "scripts/run.py").write_text("pass\n")
    receipt = artifacts.reproduce_run(root, "a", tmp_path / "reproduction")
    assert receipt["exit_code"] == 0
    assert not receipt["success"]
    assert not receipt["comparisons"]["outputs/a/result.csv"]["dimensions"]["missingness"]


def test_reproduction_refuses_to_write_inside_repository(specimen):
    root, _ = specimen
    with pytest.raises(artifacts.ArtifactError, match="outside"):
        artifacts.reproduce_run(root, "a", root / "outputs/new-run")


@pytest.mark.parametrize("change", ["other-producer", "original-output", "environment-output"])
def test_reproduction_cannot_bypass_registered_producer_or_workspace(specimen, tmp_path, change):
    root, document = specimen
    recipe = document["runs"][0]["reproduction"]
    if change == "other-producer":
        recipe["argv"] = ["{python}", "-c", "print('not the registered producer')"]
    elif change == "original-output":
        recipe["argv"][-1] = str(root / "outputs/a/result.csv")
    else:
        recipe["environment"] = {"YSTWIN_OUTPUTS": str(root / "outputs")}
    _write_registry(root, document)
    before = _snapshot(root)
    with pytest.raises(artifacts.ArtifactError):
        artifacts.reproduce_run(root, "a", tmp_path / "reproduction")
    assert _snapshot(root) == before


def test_blind_restamp_is_refused_without_writing(specimen):
    root, _ = specimen
    before = _snapshot(root)
    with pytest.raises(artifacts.ArtifactError, match="blind restamping"):
        repro.check_output_provenance(root, {}, restamp=True)
    assert _snapshot(root) == before


@pytest.mark.parametrize("change,dimension", [
    ("id,value,unit\na,2,mg\nb,,mg\n", "values"),
    ("id,value,unit\na,1,mg\nb,0,mg\n", "missingness"),
    ("id,value,unit\na,1,g\nb,,mg\n", "units"),
    ("id,value,unit\nz,1,mg\nb,,mg\n", "row_keys"),
    ("id,value,unit\na,1,mg\na,,mg\n", "row_keys"),
])
def test_csv_semantics_cover_all_four_dimensions(specimen, tmp_path, change, dimension):
    root, document = specimen
    comparison = document["runs"][0]["artifacts"][0]["comparison"]
    candidate = tmp_path / "candidate.csv"
    candidate.write_text(change)
    report = artifacts.semantic_comparison(root / "outputs/a/result.csv", candidate, comparison)
    assert not report["matches"]
    assert not report["dimensions"][dimension]


def test_csv_row_order_is_not_row_identity(specimen, tmp_path):
    root, document = specimen
    candidate = tmp_path / "candidate.csv"
    candidate.write_text("unit,value,id\nmg,,b\nmg,1.0,a\n")
    comparison = document["runs"][0]["artifacts"][0]["comparison"]
    assert artifacts.semantic_comparison(root / "outputs/a/result.csv", candidate, comparison)["matches"]


def test_json_row_identity_missingness_and_units(tmp_path):
    original, candidate = tmp_path / "one.json", tmp_path / "two.json"
    contract = {"format": "json", "row_keys": {"/rows": ["id"]}, "units": {"/unit": "mg"}}
    data = {"rows": [{"id": "a", "value": 1}, {"id": "b", "value": None}], "unit": "mg"}
    original.write_text(json.dumps(data))
    data["rows"].reverse()
    candidate.write_text(json.dumps(data))
    assert artifacts.semantic_comparison(original, candidate, contract)["matches"]
    data["rows"][0]["value"] = 0
    data["unit"] = "g"
    candidate.write_text(json.dumps(data))
    report = artifacts.semantic_comparison(original, candidate, contract)
    assert not report["dimensions"]["missingness"]
    assert not report["dimensions"]["units"]


def test_frozen_manifest_compares_original_bytes_not_current_code(specimen):
    root, document = specimen
    frozen_dir = root / "data/frozen"
    frozen_dir.mkdir()
    public = frozen_dir / "result.json"
    public.write_text('{"frozen":true}')
    manifest = frozen_dir / "manifest.json"
    manifest.write_text(json.dumps({"records": [{"origin_path": "outputs/b/frozen.json", "public_path": "result.json",
                                                "original_sha256": artifacts.sha256(public), "public_sha256": artifacts.sha256(public)}]}))
    document["frozen_bundles"] = [{"id": "historical", "manifest": "data/frozen/manifest.json", "sha256": artifacts.sha256(manifest),
                                   "reproduction_command": "use the retained original runtime, not current source"}]
    _write_registry(root, document)
    (root / "scripts/run.py").write_text("changed = True\n")
    rows = artifacts.audit_registry(root)
    frozen = next(row for row in rows if row["check"].endswith("data/frozen/result.json"))
    assert frozen["status"] == "PASS"
    assert "not a new fit or a runtime replay" in frozen["detail"]
    public.write_text('{"frozen":false}')
    assert next(row for row in artifacts.audit_registry(root) if row["check"].endswith("data/frozen/result.json"))["status"] == "FAIL"


def test_default_audit_cli_is_read_only_and_keeps_failure_exit_code(specimen, monkeypatch):
    root, _ = specimen
    monkeypatch.setattr(repro, "REPO", root)
    monkeypatch.setattr(repro, "collect", lambda *args, **kwargs: artifacts.audit_registry(root))
    before = _snapshot(root)
    assert repro.main(["--quiet"]) == 1
    assert _snapshot(root) == before


def test_explicit_audit_output_is_json_and_csv_outside_repository(specimen, tmp_path, monkeypatch):
    root, _ = specimen
    monkeypatch.setattr(repro, "REPO", root)
    monkeypatch.setattr(repro, "collect", lambda *args, **kwargs: artifacts.audit_registry(root))
    before = _snapshot(root)
    target = tmp_path / "audit"
    assert repro.main(["--quiet", "--output-dir", str(target)]) == 1
    assert _snapshot(root) == before
    report = json.loads((target / "audit_reproducibility.json").read_text())
    assert report["summary"]["failed"] == 1
    assert report["investigation_queue"][0]["check"].endswith("outputs/a/result.csv")
    assert (target / "audit_reproducibility.csv").is_file()


def test_audit_cli_rejects_in_repository_output_and_restamp(specimen, monkeypatch):
    root, _ = specimen
    monkeypatch.setattr(repro, "REPO", root)
    before = _snapshot(root)
    for arguments in (["--output-dir", str(root / "outputs")], ["--restamp"]):
        with pytest.raises(SystemExit) as error:
            repro.main(arguments)
        assert error.value.code == 2
    assert _snapshot(root) == before


def test_resolver_introspection_does_not_create_outputs(tmp_path, monkeypatch):
    from ystwin import paths

    target = tmp_path / "does-not-exist"
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(target))
    repro.describe_resolvers(paths)
    assert not target.exists()


def test_nested_json_parse_failure_is_a_real_audit_failure(specimen):
    root, _ = specimen
    (root / "outputs/b/result.json").write_text('{"x":1,"x":2}')
    row = next(row for row in repro.check_output_tables(root) if row["check"].endswith("outputs/b/result.json"))
    assert row["status"] == "FAIL"


def test_path_discovery_preserves_nested_reads_and_writes(specimen):
    root, _ = specimen
    source = root / "scripts/run.py"
    source.write_text(
        "from pathlib import Path\n"
        "out = Path('outputs/a')\n"
        "data = Path('data/a')\n"
        "source = out / 'input.json'\n"
        "source.read_text()\n"
        "(data / 'result.json').read_text()\n"
        "(out / 'result.json').write_text('{}')\n")
    candidates = repro.writers_from_scripts(root)
    assert candidates["run.py"] == ["outputs/a/result.json"]
    assert repro.outputs_read_by(source) == {"outputs/a/input.json"}
    assert repro._writers_of("outputs/b/result.json", candidates) == []
    assert repro._writers_of("data/a/result.json", candidates) == []


def test_library_owned_exclusive_json_writer_is_discovered(specimen):
    root, _ = specimen
    (root / "src/ystwin/writer.py").write_text(
        "from pathlib import Path\n"
        "def emit():\n"
        "    with Path('outputs/a/result.json').open('x') as stream:\n"
        "        stream.write('{}')\n")
    (root / "scripts/run.py").write_text("from ystwin.writer import emit\nemit()\n")
    assert repro.writers_from_scripts(root)["run.py"] == ["outputs/a/result.json"]


def test_filename_allowlist_cannot_excuse_an_unregistered_writer(specimen):
    root, _ = specimen
    (root / "data/artifact_registry.json").rename(root / "data/unregistered.json")
    (root / "outputs/consolidation_checks").mkdir()
    (root / "outputs/consolidation_checks/new-result.json").write_text('{"value":42}')
    rows = repro.check_output_freshness(root, {})
    row = next(row for row in rows if row["check"].endswith("outputs/consolidation_checks/new-result.json"))
    assert row["status"] == "FAIL"
    assert row["actual"] == "investigation_pending"


def test_registered_writer_map_never_attributes_an_input_as_output(specimen):
    root, _ = specimen
    assert artifacts.registered_writers(root) == {"run.py": ["outputs/a/result.csv"]}


def test_full_audit_subprocess_does_not_create_outputs_or_bytecode(tmp_path):
    root = tmp_path / "minimal"
    (root / "scripts").mkdir(parents=True)
    (root / "src/ystwin").mkdir(parents=True)
    (root / "scripts/audit_reproducibility.py").write_bytes((ROOT / "scripts/audit_reproducibility.py").read_bytes())
    (root / "src/ystwin/artifacts.py").write_bytes((ROOT / "src/ystwin/artifacts.py").read_bytes())
    (root / "src/ystwin/__init__.py").write_text("")
    (root / "src/ystwin/paths.py").write_text(
        "from pathlib import Path\nREPO_ROOT = Path(__file__).resolve().parents[2]\n"
        "__all__ = ['outputs_dir', 'data_dir']\n"
        "def _resolve(*args):\n    return None\n"
        "def outputs_dir():\n    out = REPO_ROOT / 'outputs'\n    out.mkdir()\n    return out\n"
        "def data_dir():\n    return REPO_ROOT / 'data'\n")
    (root / "pyproject.toml").write_text('[project]\nrequires-python = ">=3.11"\ndependencies = []\n')
    before = _snapshot(root)
    environment = {key: value for key, value in os.environ.items() if key != "PYTHONDONTWRITEBYTECODE"}
    result = subprocess.run([sys.executable, "scripts/audit_reproducibility.py", "--quiet"],
                            cwd=root, env=environment, capture_output=True, text=True, timeout=30)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "read-only; no report files written" in result.stdout
    assert not (root / "outputs").exists()
    assert _snapshot(root) == before


def test_indexed_artifacts_are_scoped_to_their_full_run_directory(specimen):
    root, document = specimen
    index = root / "outputs/a/index.json"
    index.write_text('[{"report":"cases/result.json"}]')
    run = document["runs"][0]
    run["artifact_indexes"] = [{"path": "outputs/a/index.json", "membership_sha256": artifacts.sha256(index),
                                "path_field": "report", "base_directory": "outputs/a"}]
    run["artifacts"].append({"path": "outputs/a/index.json"})
    _write_registry(root, document)
    registry = artifacts.load_registry(root)
    assert "outputs/a/cases/result.json" in registry["artifacts"]
    assert "outputs/b/cases/result.json" not in registry["artifacts"]
    index.write_text('[{"report":"cases/other.json"}]')
    with pytest.raises(artifacts.ArtifactError, match="membership index changed"):
        artifacts.load_registry(root)


def test_input_cannot_be_reclassified_as_the_runs_output(specimen):
    root, document = specimen
    document["runs"][0]["dependencies"] = [{"path": "outputs/a/result.csv", "role": "input"}]
    _write_registry(root, document)
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "input is also an output" in row["detail"]


def test_arbitrary_historical_label_does_not_turn_red_green(specimen):
    root, document = specimen
    run = document["runs"][0]
    run["status"] = "verified"
    run["lifecycle"] = "historical_frozen"
    run["evidence"] = ["data/a/source.csv"]
    run["artifacts"][0]["sha256"] = artifacts.sha256(root / "outputs/a/result.csv")
    run["frozen_runtime"] = {"claimed": "original"}
    _write_registry(root, document)
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "original code snapshot" in row["detail"]


def test_original_runtime_and_code_snapshot_not_current_source_bind_frozen_evidence(specimen):
    root, document = specimen
    snapshot = root / "data/original.py"
    snapshot.write_text((root / "scripts/run.py").read_text())
    run = document["runs"][0]
    run["status"] = "verified"
    run["lifecycle"] = "historical_frozen"
    run["evidence"] = [{"path": "data/a/source.csv", "sha256": artifacts.sha256(root / "data/a/source.csv")}]
    run["artifacts"][0]["sha256"] = artifacts.sha256(root / "outputs/a/result.csv")
    run["frozen_runtime"] = {"python": "original-runtime-version", "packages": {},
                             "code": [{"path": "data/original.py", "sha256": artifacts.sha256(snapshot)}]}
    _write_registry(root, document)
    (root / "scripts/run.py").write_text("current_implementation_changed = True\n")
    assert artifacts.audit_registry(root)[0]["status"] == "PASS"
    snapshot.write_text("original_implementation_changed = True\n")
    assert artifacts.audit_registry(root)[0]["status"] == "FAIL"


def test_numeric_looking_text_is_not_retyped_without_a_numeric_contract(tmp_path):
    left, right = tmp_path / "a.csv", tmp_path / "b.csv"
    left.write_text("id,label,value\na,001,9007199254740992\n")
    right.write_text("id,label,value\na,1,9007199254740993\n")
    report = artifacts.semantic_comparison(left, right, {"format": "csv", "row_keys": ["id"], "units": {"value": "count"}})
    assert not report["matches"]
    assert report["differences"] == ["values: ('a',): label", "values: ('a',): value"]


def test_dynamic_import_library_writer_and_builtin_open_helpers_are_preserved(specimen):
    root, _ = specimen
    (root / "src/ystwin/writer.py").write_text(
        "def save(path, payload):\n    with open(path, 'w') as stream:\n        stream.write(payload)\n"
        "def emit():\n    save('outputs/a/result.json', '{}')\n")
    (root / "scripts/run.py").write_text("import importlib\nwriter = importlib.import_module('ystwin.writer')\nwriter.emit()\n")
    assert repro.writers_from_scripts(root)["run.py"] == ["outputs/a/result.json"]


def _staged(specimen, tmp_path):
    root, _ = specimen
    directory = tmp_path / "stage"
    receipt = artifacts.reproduce_run(root, "a", directory)
    return directory / "reproduction.json", receipt


def _review(root, path, destination, *, negative=False):
    review = artifacts.review_template(root, path)
    review.update(approved=True, reviewer="fixture reviewer", reason="Explicit test-only scientific change")
    for approval in review["artifacts"].values():
        approval["reason"] = "Checked this exact baseline/candidate pair"
        approval["accepted_dimensions"] = approval["observed_changed_dimensions"]
        if negative:
            approval["approved_negative_case_changes_sha256"] = approval["required_negative_case_changes_sha256"]
    destination.write_text(json.dumps(review))
    return destination


def test_durable_adoption_changes_registry_not_retained_output(specimen, tmp_path):
    root, _ = specimen
    path, receipt = _staged(specimen, tmp_path)
    assert receipt["execution_complete"] and receipt["success"]
    before = (root / "outputs/a/result.csv").read_bytes()
    result = artifacts.adopt_reproduction(root, path)
    assert result["outputs_replaced"] == []
    assert (root / "outputs/a/result.csv").read_bytes() == before
    registry = artifacts.load_registry(root)
    assert result["receipt_id"] in registry["receipts"]
    assert artifacts.audit_registry(root)[0]["status"] == "PASS"
    assert not (root / "outputs/provenance.json").exists()


@pytest.mark.parametrize("change", ["input", "code", "parameters", "runtime", "candidate", "baseline", "baseline_policy", "child_record", "retained"])
def test_adoption_rejects_stale_or_tampered_receipts(specimen, tmp_path, monkeypatch, change):
    root, document = specimen
    path, receipt = _staged(specimen, tmp_path)
    assert receipt["success"]
    if change == "input":
        (root / "data/a/source.csv").write_text("id,value,unit\na,9,mg\nb,,mg\n")
    elif change == "code":
        source = root / "src/ystwin/writer.py"
        source.write_text(source.read_text() + "\n")
    elif change == "parameters":
        document["runs"][0]["parameters"]["selection"] = "different selection"
        _write_registry(root, document)
    elif change == "baseline_policy":
        document["baseline_ref"] = "a" * 40
        _write_registry(root, document)
    elif change == "runtime":
        original = artifacts.runtime_identity
        monkeypatch.setattr(artifacts, "runtime_identity", lambda r: {**original(r), "different_runtime": True})
    elif change in {"candidate", "baseline"}:
        parent = "workspace" if change == "candidate" else "baseline"
        (path.parent / parent / "outputs/a/result.csv").write_text("id,value,unit\na,9,mg\nb,,mg\n")
    elif change == "child_record":
        execution = path.parent / "child_execution.json"
        value = json.loads(execution.read_text())
        value["exit_code"] = 8
        execution.write_text(json.dumps(value))
    else:
        (root / "outputs/a/result.csv").write_text("id,value,unit\na,9,mg\nb,,mg\n")
    before = _snapshot(root)
    with pytest.raises(artifacts.ArtifactError):
        artifacts.adopt_reproduction(root, path, install=True)
    assert _snapshot(root) == before


def test_exit_zero_without_child_execution_record_cannot_be_adopted(specimen, tmp_path):
    root, _ = specimen
    script = root / "scripts/run.py"
    script.write_text(script.read_text() + "import os\nos._exit(0)\n")
    path, receipt = _staged(specimen, tmp_path)
    assert receipt["exit_code"] == 0
    assert not receipt["execution_complete"]
    with pytest.raises(artifacts.ArtifactError, match="child execution"):
        artifacts.adopt_reproduction(root, path)


def test_scientific_changes_need_flag_and_exact_review_before_installation(specimen, tmp_path):
    root, _ = specimen
    (root / "data/a/source.csv").write_text("id,value,unit\na,2,mg\nb,,mg\n")
    path, receipt = _staged(specimen, tmp_path)
    assert receipt["execution_complete"] and not receipt["success"]
    review = _review(root, path, tmp_path / "review.json")
    before = _snapshot(root)
    for options in ({}, {"accept_scientific_changes": True}, {"review_path": review}):
        with pytest.raises(artifacts.ArtifactError, match="reviewed acceptance"):
            artifacts.adopt_reproduction(root, path, install=True, **options)
        assert _snapshot(root) == before
    decision = artifacts.adopt_reproduction(root, path, accept_scientific_changes=True, review_path=review)
    assert not decision["decisions"]["outputs/a/result.csv"]["installed"]
    assert artifacts.audit_registry(root)[0]["status"] == "FAIL"
    result = artifacts.adopt_reproduction(root, path, accept_scientific_changes=True, review_path=review, install=True)
    assert result["outputs_replaced"] == ["outputs/a/result.csv"]
    assert artifacts.audit_registry(root)[0]["status"] == "PASS"


def test_negative_case_removal_requires_explicit_negative_review(specimen, tmp_path):
    root, document = specimen
    original = "id,value,unit,status\na,1,mg,ok\nb,,mg,failed\n"
    (root / "outputs/a/result.csv").write_text(original)
    (root / "data/a/source.csv").write_text("id,value,unit,status\na,1,mg,ok\n")
    document["runs"][0]["artifacts"][0]["comparison"]["negative_cases"] = {"status": {"values": ["failed"]}}
    _write_registry(root, document)
    path, receipt = _staged(specimen, tmp_path)
    comparison = receipt["comparisons"]["outputs/a/result.csv"]
    assert comparison["negative_case_changes"][0]["before"] == "failed"
    review = _review(root, path, tmp_path / "review.json")
    with pytest.raises(artifacts.ArtifactError, match="negative cases"):
        artifacts.adopt_reproduction(root, path, accept_scientific_changes=True, review_path=review, install=True)
    assert (root / "outputs/a/result.csv").read_text() == original
    review = _review(root, path, tmp_path / "review.json", negative=True)
    result = artifacts.adopt_reproduction(root, path, accept_scientific_changes=True, review_path=review, install=True)
    stored = artifacts.load_registry(root)["receipts"][result["receipt_id"]]
    assert stored["comparisons"]["outputs/a/result.csv"]["negative_case_changes"] == comparison["negative_case_changes"]
    assert stored["retained_sha256"]["outputs/a/result.csv"] != stored["candidate_sha256"]["outputs/a/result.csv"]


def test_changed_input_after_review_cannot_be_approved_by_an_old_review(specimen, tmp_path):
    root, _ = specimen
    (root / "data/a/source.csv").write_text("id,value,unit\na,2,mg\nb,,mg\n")
    path, _ = _staged(specimen, tmp_path)
    review = _review(root, path, tmp_path / "review.json")
    (root / "data/a/source.csv").write_text("id,value,unit\na,3,mg\nb,,mg\n")
    before = _snapshot(root)
    with pytest.raises(artifacts.ArtifactError, match="stale"):
        artifacts.adopt_reproduction(root, path, accept_scientific_changes=True, review_path=review, install=True)
    assert _snapshot(root) == before


def test_unreviewed_worktree_replacement_does_not_replace_the_git_baseline(specimen, tmp_path, monkeypatch):
    root, document = specimen
    baseline = (root / "outputs/a/result.csv").read_bytes()
    document["baseline_ref"] = "a" * 40
    _write_registry(root, document)
    for name in ("data/a/source.csv", "outputs/a/result.csv"):
        (root / name).write_text("id,value,unit\na,2,mg\nb,,mg\n")
    original_run = subprocess.run

    def run(command, **kwargs):
        if command[:2] == ["git", "show"]:
            return subprocess.CompletedProcess(command, 0, baseline, b"")
        return original_run(command, **kwargs)

    monkeypatch.setattr(artifacts.subprocess, "run", run)
    path, receipt = _staged(specimen, tmp_path)
    assert receipt["execution_complete"] and not receipt["success"]
    assert receipt["working_sha256"]["outputs/a/result.csv"] == receipt["candidate_sha256"]["outputs/a/result.csv"]
    with pytest.raises(artifacts.ArtifactError, match="scientific changes"):
        artifacts.adopt_reproduction(root, path)


def test_multiset_contract_preserves_duplicate_negative_rows_and_their_count(tmp_path):
    original, candidate = tmp_path / "old.csv", tmp_path / "new.csv"
    original.write_text("kind,seed,status\nrefusal,,failed\nrefusal,,failed\n")
    candidate.write_text("kind,seed,status\nrefusal,,failed\n")
    contract = {"format": "csv", "row_keys": ["kind", "seed"], "nullable_row_keys": ["seed"],
                "duplicate_keys": "multiset", "units": {}, "negative_cases": {"status": {"values": ["failed"]}}}
    assert artifacts.semantic_comparison(original, original, contract)["matches"]
    result = artifacts.semantic_comparison(original, candidate, contract)
    assert not result["matches"]
    assert len(result["negative_case_changes"]) == 1


def test_unclassified_local_is_visible_but_tracked_unowned_is_blocking(specimen):
    root, _ = specimen
    subprocess.run(["git", "init", "-q"], cwd=root, check=True, capture_output=True)
    extra = root / "outputs/b/result.csv"
    extra.write_text("id,value\na,3\n")
    row = next(row for row in artifacts.audit_registry(root) if row["check"].endswith("outputs/b/result.csv"))
    assert row["status"] == "SKIP" and row["actual"] == "local_decision_pending"
    subprocess.run(["git", "add", "outputs/b/result.csv"], cwd=root, check=True, capture_output=True)
    row = next(row for row in artifacts.audit_registry(root) if row["check"].endswith("outputs/b/result.csv"))
    assert row["status"] == "FAIL"


def test_declared_local_run_never_exempts_a_later_tracked_member(specimen):
    root, document = specimen
    subprocess.run(["git", "init", "-q"], cwd=root, check=True, capture_output=True)
    metadata = root / "outputs/b/protocol.json"
    metadata.write_text(json.dumps({"mode": "synthetic", "code_hashes": {"scripts/run.py": "b" * 64}}))
    decision = artifacts.describe_local_run(root, "outputs/b", "outputs/b/protocol.json", "scripts/run.py",
                                            {"/mode": "synthetic"}, "Explicit local synthetic fixture")
    document["local_runs"] = [decision]
    _write_registry(root, document)
    row = next(row for row in artifacts.audit_registry(root) if row["check"].endswith("outputs/b/protocol.json"))
    assert row["status"] == "SKIP" and row["actual"] == "local_exploratory"
    subprocess.run(["git", "add", "outputs/b/protocol.json"], cwd=root, check=True, capture_output=True)
    row = next(row for row in artifacts.audit_registry(root) if row["check"].endswith("outputs/b/protocol.json"))
    assert row["status"] == "FAIL"


def test_private_input_outside_repository_is_not_read(specimen, tmp_path, monkeypatch):
    root, document = specimen
    private = tmp_path / "private.csv"
    private.write_text("not an approved input")
    document["runs"][0]["dependencies"] = [{"path": "data/missing.csv", "role": "input",
                                             "source": {"environment": "TEST_PRIVATE_INPUT"}}]
    monkeypatch.setenv("TEST_PRIVATE_INPUT", str(private))
    _write_registry(root, document)
    with pytest.raises(artifacts.ArtifactError, match="outside-repository input access is disabled"):
        artifacts.run_identity(root, artifacts.load_registry(root)["runs"]["a"])


def test_changing_an_accepted_hash_cannot_resign_an_unwitnessed_output(specimen, tmp_path):
    root, _ = specimen
    path, _ = _staged(specimen, tmp_path)
    artifacts.adopt_reproduction(root, path)
    output = root / "outputs/a/result.csv"
    output.write_text("id,value,unit\na,999,mg\nb,,mg\n")
    document = artifacts.read_json(root / artifacts.REGISTRY_PATH)
    document["runs"][0]["verification"]["artifacts"]["outputs/a/result.csv"]["accepted_sha256"] = artifacts.sha256(output)
    _write_registry(root, document)
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "not witnessed" in row["detail"]


def test_nested_json_keys_preserve_negative_cases_and_synthetic_scope(tmp_path):
    left, right = tmp_path / "old.json", tmp_path / "new.json"
    source = {"biological_validation": False, "results": [
        {"sample": "a", "allocation": {"policy": "mixture", "status": "infeasible", "unit": "1/h"}},
        {"sample": "a", "allocation": {"policy": "fixed", "status": "optimal", "unit": "1/h"}},
    ]}
    left.write_text(json.dumps(source))
    contract = {"format": "json", "row_keys": {"/results": ["sample", "/allocation/policy"]},
                "units": {}, "unit_field_names": ["unit"], "invariants": {"/biological_validation": False},
                "negative_json_fields": {"/results/*/allocation/status": {"values": ["infeasible"]}}}
    right.write_text(json.dumps({**source, "results": source["results"][::-1]}))
    assert artifacts.semantic_comparison(left, right, contract)["matches"]
    right.write_text(json.dumps({**source, "results": source["results"][1:]}))
    changed = artifacts.semantic_comparison(left, right, contract)
    assert not changed["matches"] and len(changed["negative_case_changes"]) == 1
    for value in (True, 0):
        right.write_text(json.dumps({**source, "biological_validation": value}))
        assert not artifacts.semantic_comparison(right, right, contract)["matches"]


def test_registered_failed_export_preserves_error_and_original_runtime(tmp_path):
    root = tmp_path / "repo"
    (root / "outputs").mkdir(parents=True)
    (root / "data").mkdir()
    partial = root / "outputs/partial.json"
    partial.write_text('{"incomplete":')
    snapshot = root / "data/original.py"
    snapshot.write_text("original_exporter = True\n")
    evidence = root / "data/export_receipt.json"
    evidence.write_text(json.dumps({"diagnostic_export": {"status": "failed_after_prediction_seal"}}))
    run = {"id": "a", "producer": {"path": "data/original.py", "callable": "main", "libraries": []},
           "reason": "Recorded synthetic fixture export failure, not a current prediction"}
    run.update(lifecycle="historical_frozen", status="verified",
               frozen_runtime={"python": "original", "packages": {}, "code": [{"path": "data/original.py", "sha256": artifacts.sha256(snapshot)}]})
    run["artifacts"] = [{"path": "outputs/partial.json", "scientific_status": "failed_after_prediction_seal",
                          "origin_sha256": artifacts.sha256(partial), "error_evidence": {
                              "path": "data/export_receipt.json", "sha256": artifacts.sha256(evidence),
                              "pointer": "/diagnostic_export/status"}}]
    _write_registry(root, {"schema_version": 1, "runs": [run]})
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "SKIP" and row["actual"] == "failed_after_prediction_seal"
    assert repro.check_output_tables(root)[0]["status"] == "SKIP"
    partial.write_text('{}')
    assert artifacts.audit_registry(root)[0]["status"] == "FAIL"


def test_durable_receipt_cannot_change_child_payload_without_its_execution_digest(specimen, tmp_path):
    root, _ = specimen
    path, _ = _staged(specimen, tmp_path)
    adopted = artifacts.adopt_reproduction(root, path)
    document = artifacts.read_json(root / artifacts.REGISTRY_PATH)
    receipt = document["receipts"][adopted["receipt_id"]]
    receipt["execution"]["read_paths"].append("data/not-actually-read.csv")
    changed_id = artifacts.canonical_digest(receipt)
    document["receipts"] = {changed_id: receipt}
    document["runs"][0]["verification"]["artifacts"]["outputs/a/result.csv"]["receipt_id"] = changed_id
    _write_registry(root, document)
    row = artifacts.audit_registry(root)[0]
    assert row["status"] == "FAIL"
    assert "child binding failed" in row["detail"]


@pytest.mark.parametrize("run_id,relative,role", [
    ("env-to-product", "data/physiology/chemostatData_VanHoek1998.tsv", "input"),
    ("product-environment-sweep", "data/physiology/chemostatData_VanHoek1998.tsv", "input"),
    ("registered-prediction-d018", "data/physiology/chemostatData_VanHoek1998.tsv", "input"),
    ("scenarios", "data/physiology/chemostatData_VanHoek1998.tsv", "input"),
    ("flux-null-baseline", "data/phb/kocharin2013_chemostat_states.tsv", "input"),
    ("thiolase-threshold", "data/gem/yeast-GEM.xml.gz", "model"),
])
def test_current_recipes_declare_the_public_inputs_their_executed_paths_require(run_id, relative, role):
    run = artifacts.load_registry(ROOT)["runs"][run_id]
    assert any(entry["path"] == relative and entry["role"] == role for entry in run["dependencies"])
    assert relative in artifacts.dependency_sources(ROOT, run)


def test_flux_null_recipe_keeps_descriptive_candidates_separate_from_random_designs():
    run = artifacts.load_registry(ROOT)["runs"]["flux-null-baseline"]
    assert {item["path"] for item in run["artifacts"]} == {
        "outputs/flux_null_baseline.csv", "outputs/flux_candidate_scores.csv"}
    candidate = next(item for item in run["artifacts"] if item["path"] == "outputs/flux_candidate_scores.csv")
    assert candidate["comparison"]["negative_cases"]["inference_status"] == {
        "values": ["pending_exchangeability"]}


def test_public_source_alias_does_not_follow_an_in_repository_symlink(specimen, tmp_path, monkeypatch):
    root, document = specimen
    outside = tmp_path / "outside.csv"
    outside.write_text("must not read this external file")
    alias = root / "data/alias.csv"
    alias.symlink_to(outside)
    document["runs"][0]["dependencies"][0]["source"] = {"environment": "PUBLIC_SOURCE_ALIAS"}
    monkeypatch.setenv("PUBLIC_SOURCE_ALIAS", str(alias))
    _write_registry(root, document)
    with pytest.raises(artifacts.ArtifactError, match="outside-repository input access is disabled|symlink"):
        artifacts.run_identity(root, artifacts.load_registry(root)["runs"]["a"])


@pytest.fixture(scope="module")
def linked_registry():
    return artifacts.load_registry(ROOT)


def _contract_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def test_linked_split_recipe_requires_the_full_public_population(linked_registry):
    run = linked_registry["runs"]["splits"]
    assert run["parameters"]["panel_only"] is False
    assert run["parameters"]["original_invocation"] == "unrecorded"
    assert run["parameters"]["real_table"] == "outputs/sensor_characterisation.csv"
    assert run["reproduction"]["argv"] == [
        "{python}", "-B", "scripts/make_splits.py", "--real-table", "outputs/sensor_characterisation.csv"]
    assert {"outputs/sensor_characterisation.csv", "data/plates/manifest.csv"} <= {
        entry["path"] for entry in run["dependencies"]}
    assert run["lifecycle"] == "current"


def test_linked_gate_recipe_keeps_august3_and_adds_current_xpt_and_manifest(linked_registry):
    run = linked_registry["runs"]["gates"]
    assert run["reproduction"]["argv"] == ["{python}", "-B", "scripts/run_gates.py", "--committed"]
    assert run["parameters"]["source_sets"] == ["july_2026_07", "newprotocol"]
    assert run["parameters"]["committed_only"] is True
    stems = ["20260701_ER_preliminary_(RAW)", "20260709_ER_stress_1st_(RAW)",
             "20260722_ERandOxidativeStress_NewProtocol_ANALYSED",
             "20260728_ERandOxidativeStress_NewProtocol_Replicate2",
             "20260803_ERandoxidativestress_Replicate3", "20260804_ER_Oxidative_Replicate4"]
    expected = {f"outputs/{prefix}_{stem}.csv" for prefix in ("g1", "d2_g1passed") for stem in stems}
    assert {item["path"] for item in run["artifacts"]} == expected | {"outputs/gates_manifest.json"}
    assert run["lifecycle"] == "current"


def test_linked_gate_csv_schema_and_units_match_the_emitted_quantities(linked_registry):
    g1 = linked_registry["artifacts"]["outputs/g1_20260804_ER_Oxidative_Replicate4.csv"]["comparison"]
    d2 = linked_registry["artifacts"]["outputs/d2_g1passed_20260804_ER_Oxidative_Replicate4.csv"]["comparison"]
    assert g1["columns"] == ["well", "passed", "failures", "max_raw_od", "min_od_above_blank",
                             "decline_fraction", "od_fold_change", "max_reporter_above_background"]
    assert d2["columns"] == ["well", "dilution_r2", "dilution_slope", "qss_fraction", "naive_fold",
                             "activity_fold", "negative_activity_fraction", "implied_min_k_deg",
                             "excess_loss_rate", "mu_min", "mu_max", "verdict"]
    assert g1["units"]["max_raw_od"] == g1["units"]["min_od_above_blank"] == "OD600"
    assert g1["units"]["max_reporter_above_background"] == "RFU"
    assert all(d2["units"][name] == "1/h" for name in ("mu_min", "mu_max", "excess_loss_rate", "implied_min_k_deg"))
    assert all(d2["units"][name] == "dimensionless" for name in (
        "dilution_r2", "dilution_slope", "qss_fraction", "naive_fold", "activity_fold", "negative_activity_fraction"))


@pytest.mark.parametrize("verdict", ["model-inadequate", "not-estimable", "extra-loss"])
@pytest.mark.parametrize("change", ["drop", "promote"])
def test_linked_d2_contract_retains_every_negative_label(linked_registry, tmp_path, verdict, change):
    contract = linked_registry["artifacts"]["outputs/d2_g1passed_20260701_ER_preliminary_(RAW).csv"]["comparison"]
    columns = ["well", *contract["units"], "verdict"]
    negative = {**dict.fromkeys(contract["units"], "0.25"), "well": "A1", "verdict": verdict}
    other = {**negative, "well": "B1", "verdict": "mixed"}
    left, right = tmp_path / "old.csv", tmp_path / "new.csv"
    _contract_csv(left, columns, [negative, other])
    replacement = [other] if change == "drop" else [{**negative, "verdict": "mixed"}, other]
    _contract_csv(right, columns, replacement)
    report = artifacts.semantic_comparison(left, right, contract)
    assert not report["matches"]
    assert {"row_key": ["A1"], "column": "verdict", "before": verdict,
            "after": None if change == "drop" else "mixed", "row_retained": change != "drop"} in report["negative_case_changes"]


def _gate_manifest_fixture():
    return {"schema_version": 1, "complete": False, "committed_only": True,
            "independent_unit": "biological plate; file representations are not independent replicates",
            "n_sources_processed": 1, "n_sources_refused": 1,
            "exports": [
                {"file": "unresolved.xlsx", "status": "refused", "source_set": "newprotocol",
                 "source_sha256": "a" * 64, "note": "correction-state contradiction"},
                {"file": "current.xpt", "status": "processed", "source_set": "newprotocol",
                 "source_sha256": "b" * 64, "n_cultures": 3, "n_pass": 2,
                 "n_d2_attempted": 2, "n_d2_scored": 0, "blank_wells": ["H1"], "unrecorded_wells": [],
                 "d2_exclusions": [{"well": "A1", "reason": "too few readings"},
                                   {"well": "A2", "reason": "nonfinite trace"}]}]}


def test_linked_gate_manifest_keys_are_source_and_well_not_array_order(linked_registry, tmp_path):
    contract = linked_registry["artifacts"]["outputs/gates_manifest.json"]["comparison"]
    left, right = tmp_path / "old.json", tmp_path / "new.json"
    payload = _gate_manifest_fixture()
    left.write_text(json.dumps(payload))
    payload["exports"][1]["d2_exclusions"].reverse()
    payload["exports"].reverse()
    right.write_text(json.dumps(payload))
    assert artifacts.semantic_comparison(left, right, contract)["matches"]
    payload["exports"].append(copy.deepcopy(payload["exports"][0]))
    right.write_text(json.dumps(payload))
    assert not artifacts.semantic_comparison(right, right, contract)["matches"]


@pytest.mark.parametrize("change", ["drop_refusal", "promote_refusal", "hide_count", "complete", "drop_exclusion", "unit", "not_committed"])
def test_linked_gate_manifest_cannot_hide_refused_sources_or_units(linked_registry, tmp_path, change):
    contract = linked_registry["artifacts"]["outputs/gates_manifest.json"]["comparison"]
    left, right = tmp_path / "old.json", tmp_path / "new.json"
    payload = _gate_manifest_fixture()
    left.write_text(json.dumps(payload))
    if change == "drop_refusal":
        payload["exports"].pop(0)
    elif change == "promote_refusal":
        payload["exports"][0]["status"] = "processed"
    elif change == "hide_count":
        payload["n_sources_refused"] = 0
    elif change == "complete":
        payload["complete"] = True
    elif change == "drop_exclusion":
        payload["exports"][1]["d2_exclusions"].pop()
    elif change == "unit":
        payload["independent_unit"] = "well"
    else:
        payload["committed_only"] = False
    right.write_text(json.dumps(payload))
    report = artifacts.semantic_comparison(left, right, contract)
    assert not report["matches"]
    if change in {"drop_refusal", "promote_refusal", "hide_count", "complete"}:
        assert report["negative_case_changes"]
    if change == "unit":
        assert not report["dimensions"]["units"]


_CROSS_TRANSFER_LEGACY = ["family", "perturbs", "transfer_r2", "oracle_r2", "alignment", "null_median",
                          "p_value", "beats_null", "skill_over_null", "delta_vs_baseline"]
_CROSS_TRANSFER_METADATA = ["model", "configuration", "target", "evaluation_unit", "selection_scope", "oracle_role"]
_CROSS_OPTIMISM_LEGACY = ["ucbog", "mean_gap", "alpha", "n_bootstrap", "n_reference", "n_candidate_domains",
                         "n_reference_domains", "clipped_fraction", "candidate_states", "candidate_design",
                         "baseline_transfer_r2", "chosen_on_mean", "held_out_mean"]
_CROSS_OPTIMISM_METADATA = ["model", "transfer_configuration", "selected_configuration", "candidate_configuration",
                           "score_target", "configuration_grid", "split_target", "split_evaluation_unit",
                           "split_selection_scope", "selected_pipeline_optimism_status", "bound_target",
                           "bound_resampling_unit", "bound_selection_scope", "bound_method", "bound_coverage_status"]


@pytest.mark.parametrize("name,legacy,metadata,json_columns", [
    ("cross_family_transfer", _CROSS_TRANSFER_LEGACY, _CROSS_TRANSFER_METADATA, ["configuration"]),
    ("cross_family_optimism", _CROSS_OPTIMISM_LEGACY, _CROSS_OPTIMISM_METADATA,
     ["transfer_configuration", "selected_configuration", "candidate_configuration", "configuration_grid",
      "split_selection_scope", "bound_selection_scope"]),
])
def test_linked_cross_family_schema_migration_and_json_semantics(linked_registry, tmp_path, name, legacy, metadata, json_columns):
    contract = linked_registry["artifacts"][f"outputs/{name}.csv"]["comparison"]
    assert contract["column_variants"] == [legacy, legacy + metadata]
    assert contract["json_columns"] == json_columns
    row = {column: "descriptive fixture" for column in legacy + metadata}
    row.update(dict.fromkeys(contract["units"], "0.25"))
    row.update(dict.fromkeys(contract.get("boolean_columns", []), "False"))
    for column in json_columns:
        row[column] = '{"n_states":3,"replicates":3}'
    left, right = tmp_path / "old.csv", tmp_path / "new.csv"
    _contract_csv(left, legacy + metadata, [row])
    reordered = {**row, **dict.fromkeys(json_columns, '{ "replicates": 3, "n_states": 3 }')}
    _contract_csv(right, legacy + metadata, [reordered])
    assert artifacts.semantic_comparison(left, right, contract)["matches"]
    for column in json_columns:
        _contract_csv(right, legacy + metadata, [{**row, column: '{"n_states":2,"replicates":3}'}])
        report = artifacts.semantic_comparison(left, right, contract)
        assert not report["matches"] and not report["dimensions"]["values"]
    _contract_csv(right, legacy + metadata[:-1], [{key: value for key, value in row.items() if key != metadata[-1]}])
    assert not artifacts.semantic_comparison(right, right, contract)["matches"]
    _contract_csv(right, legacy, [{key: row[key] for key in legacy}])
    assert artifacts.semantic_comparison(right, right, contract)["matches"]
    assert not artifacts.semantic_comparison(left, right, contract)["matches"]


def test_linked_cross_family_recipe_identifies_the_spota_master_seed(linked_registry):
    run = linked_registry["runs"]["cross-family"]
    assert "candidate_seed" not in run["parameters"]
    assert run["parameters"]["spota_seed"] == 0
    assert run["parameters"]["family_seed"] == run["parameters"]["data_seed"] == 0
    assert run["parameters"]["scale"] is False
    assert run["parameters"]["transfer_configuration"]["n_states"] == 3
    assert run["parameters"]["n_reference"] == 8
    assert run["parameters"]["n_candidate_domains"] == 4
    assert run["parameters"]["n_reference_domains"] == 2
    assert run["parameters"]["original_invocation"] == "unrecorded"
    assert run["lifecycle"] == "current"


def test_linked_thiolase_inputs_are_exact_loader_dependencies_not_the_download_directory(linked_registry):
    run = linked_registry["runs"]["thiolase-threshold"]
    required = {"data/thermo/aliases.tsv", "data/thermo/thermo_data.thermodb",
                "data/thermo/heterologous_metabolites.tsv", "data/thermo/refuted_energies.tsv"}
    thermo = {entry["path"] for entry in run["dependencies"] if entry["path"].startswith("data/thermo")}
    assert thermo == required
    assert required <= artifacts.dependency_sources(ROOT, run).keys()
    public = artifacts.read_json(ROOT / "data/public_inputs.json")
    for name in ("data/thermo/aliases.tsv", "data/thermo/thermo_data.thermodb"):
        pin = next(entry for entry in public["files"] if entry["path"] == name)
        assert artifacts.sha256(ROOT / name) == pin["sha256"]
    assert run["reproduction"]["environment"]["YSTWIN_THERMO"] == "{workspace}/data/thermo"


def test_linked_public_plan_is_read_only_and_preserves_upstream_order():
    watched = [ROOT / name for name in ("data/artifact_registry.json", "data/current_claims.json", "src/ystwin/artifacts.py")]
    before = {path: (artifacts.sha256(path), path.stat().st_mtime_ns) for path in watched}
    plans = artifacts.plan_current_runs(ROOT)
    indexed = {row["run_id"]: row for row in plans}
    assert indexed["splits"]["upstream_runs"] == ["sensor-characterisation"]
    assert [row["run_id"] for row in plans].index("sensor-characterisation") < [row["run_id"] for row in plans].index("splits")
    for run_id in ("splits", "gates", "cross-family", "thiolase-threshold"):
        assert not indexed[run_id]["input_issues"]
    assert before == {path: (artifacts.sha256(path), path.stat().st_mtime_ns) for path in watched}


@pytest.mark.parametrize("binding", ["argv", "metadata_contract", "input_roster"])
def test_linked_recipe_changes_invalidate_identity_without_rewriting_a_receipt(specimen, binding):
    root, document = specimen
    run = artifacts.load_registry(root)["runs"]["a"]
    original_identity = artifacts.run_identity(root, run)
    unchanged_receipt = json.dumps({"identity": original_identity}, sort_keys=True)
    if binding == "argv":
        document["runs"][0]["reproduction"]["argv"].extend(["--real-table", "outputs/source.csv"])
    elif binding == "metadata_contract":
        document["runs"][0]["artifacts"][0]["comparison"]["json_columns"] = ["configuration"]
    else:
        document["runs"][0]["dependencies"].append({"path": "data/b/source.csv", "role": "input"})
    _write_registry(root, document)
    current_identity = artifacts.run_identity(root, artifacts.load_registry(root)["runs"]["a"])
    assert current_identity != original_identity
    assert artifacts.canonical_digest(current_identity) != artifacts.canonical_digest(original_identity)
    assert json.dumps({"identity": original_identity}, sort_keys=True) == unchanged_receipt


@pytest.mark.parametrize("run_id,lifecycle", [("public-replay-record", "historical_frozen")])
def test_linked_noncurrent_dispositions_pin_original_git_bytes_and_role_evidence(linked_registry, run_id, lifecycle):
    import hashlib

    run = linked_registry["runs"][run_id]
    assert run["lifecycle"] == lifecycle and run["status"] == "verified"
    assert run["retention"] == "retained"
    assert not run.get("verification")
    assert not artifacts._pinned_references(ROOT, run["evidence"], "lifecycle evidence")
    for item in run["artifacts"]:
        original = subprocess.run(["git", "show", f"{linked_registry['document']['baseline_ref']}:{item['path']}"],
                                  cwd=ROOT, capture_output=True, check=True)
        assert item["sha256"] == hashlib.sha256(original.stdout).hexdigest() == artifacts.sha256(ROOT / item["path"])
    if lifecycle == "historical_frozen":
        assert not artifacts._historical_runtime_issues(ROOT, run)
        report = artifacts.read_json(ROOT / run["artifacts"][0]["path"])
        assert report["biological_validation"] is False and report["independent_test"] is False
        assert report["new_fitting"] is False and report["new_selection"] is False
        evidence = artifacts.read_json(ROOT / run["evidence"][0]["path"])
        assert evidence["evidence"]["public_replay_sha256"] == run["artifacts"][0]["sha256"]
    else:
        context = artifacts.read_json(ROOT / run["evidence"][0]["path"])
        assert context["historical_demo_overwritten"] is False and context["biological_evidence"] is False
        assert run["parameters"]["original_invocation"] == "unrecorded"
        assert not run.get("frozen_runtime")
        assert all((ROOT / relative).is_file() for relative in run["superseded_by"])


_DISPOSITION_BASELINE = "84561f51502b5b4bf66c052dc88c853038917934"
_DIAGNOSTIC_DISPOSITIONS = (
    "outputs/audit_claims.csv", "outputs/audit_determinism.csv", "outputs/audit_output_tables.csv",
    "outputs/audit_reproducibility.csv", "outputs/provenance.json",
    "outputs/biology_learning_audit_01/diagnostic_plan.json",
    "outputs/biology_learning_audit_01/ledger.json",
    "outputs/biology_learning_audit_01/native_diagnostics.json",
    "outputs/biology_learning_audit_02_gate_hardening/diagnostic_plan.json",
    "outputs/biology_learning_audit_02_gate_hardening/ledger.json",
    "outputs/biology_learning_audit_02_gate_hardening/native_diagnostics.json",
    "outputs/biology_learning_review.json", "outputs/biology_observability_gauge.json",
    "outputs/heldout_transfer_summary.json", "outputs/partial_order_v2/run_context.json",
    "outputs/consolidation_checks/portable_review_resolved.json",
    "outputs/consolidation_checks/reconciliation.json", "outputs/consolidation_checks/verification.json",
    "outputs/native_reconciliation_20260907_01/findings.json",
    "outputs/native_reconciliation_20260907_01/followup_checks.json",
    "outputs/native_reconciliation_20260907_01/source_lineage_audit.json",
    "outputs/native_reconciliation_20260907_01/verification.json",
    "outputs/native_population_development/executor_refresh_01/development_protocol.json",
    "outputs/native_population_development/executor_refresh_01/executor_refresh_request.json",
    "outputs/native_population_development/native_training_01_phase_c_01_scoring/independent_review.json",
)
_RECONCILIATION_DISPOSITION_BLOCKERS = tuple(
    f"outputs/native_reconciliation_20260907_01/final/{name}" for name in (
        "case_summary.csv", "case_summary.json", "completion.json", "native_conditions.json",
        "native_structure.json", "protocol.json"))


@pytest.fixture(scope="module")
def disposition_git_bytes():
    return {name: subprocess.run(["git", "show", f"{_DISPOSITION_BASELINE}:{name}"], cwd=ROOT,
                                 capture_output=True, check=True).stdout
            for name in (*_DIAGNOSTIC_DISPOSITIONS, *_RECONCILIATION_DISPOSITION_BLOCKERS)}


def _assert_disposition_evidence(document, originals):
    import hashlib

    assert document["baseline_ref"] == _DISPOSITION_BASELINE
    assert document["records"] == []
    assert document["current_scientific_support"] is False
    assert document["historical_execution_verified"] is False
    archived = document["external_references"]
    blocked = document["blocked_records"]
    assert len(archived) == 25 and {item["origin_path"] for item in archived} == set(_DIAGNOSTIC_DISPOSITIONS)
    assert len(blocked) == 6 and {item["origin_path"] for item in blocked} == set(_RECONCILIATION_DISPOSITION_BLOCKERS)
    roles = {"historical_repository_diagnostic", "historical_diagnostic_protocol", "historical_source_audit",
             "historical_retrospective_diagnostic", "historical_negative_assessment",
             "historical_operational_record", "historical_engineering_protocol"}
    for item in (*archived, *blocked):
        original = originals[item["origin_path"]]
        assert item["original_sha256"] == hashlib.sha256(original).hexdigest()
        assert item["original_bytes"] == len(original)
        assert item["basis"]
        if item in archived:
            assert item["role"] in roles
        for pointer, expected in item.get("source_assertions", {}).items():
            actual = artifacts.json_pointer(json.loads(original), pointer)
            assert type(actual) is type(expected) and actual == expected
    aliases = {item["origin_path"]: item["recorded_alias_of"] for item in archived if "recorded_alias_of" in item}
    assert aliases == {
        "outputs/biology_learning_audit_02_gate_hardening/native_diagnostics.json":
            "outputs/biology_learning_audit_01/native_diagnostics.json",
        "outputs/native_population_development/executor_refresh_01/development_protocol.json":
            "data/native_law_v2/development_protocol.json",
    }
    for path, alias in aliases.items():
        source = originals.get(alias)
        if source is None:
            source = subprocess.run(["git", "show", f"{_DISPOSITION_BASELINE}:{alias}"], cwd=ROOT,
                                    capture_output=True, check=True).stdout
        assert originals[path] == source
    assert not any(name in json.dumps(document) for name in ("data/artifact_registry.json", "data/repository_inventory.json"))


def test_noncurrent_disposition_ledger_is_original_git_integrity_not_a_reproduction_receipt(disposition_git_bytes):
    ledger = artifacts.read_json(ROOT / "data/artifact_dispositions.json")
    _assert_disposition_evidence(ledger, disposition_git_bytes)
    registry = artifacts.load_registry(ROOT)
    bundle = next(item for item in registry["frozen_bundles"] if item["id"] == "historical-diagnostic-records")
    assert bundle["manifest"] == "data/artifact_dispositions.json"
    assert bundle["sha256"] == artifacts.sha256(ROOT / bundle["manifest"])
    assert bundle["include_external_outputs"] is True
    for name in _DIAGNOSTIC_DISPOSITIONS:
        item = registry["artifacts"][name]
        assert item["lifecycle"] == "historical_frozen" and item["bundle"] == bundle
        assert item["run"] not in registry["runs"]
    restoration, = ledger["restoration_requests"]
    assert restoration["path"] == "outputs/audit_output_tables.csv"
    assert restoration["required_sha256"] == "7fb36b639fb8a0797e784a8e1a1e19c7ba300c2d059850eb7ee64b6858305bb9"


def test_noncurrent_dispositions_do_not_downgrade_current_claims_or_recipe_inputs():
    from ystwin.analysis.claims import load_registry as load_claims

    registry = artifacts.load_registry(ROOT)
    claims = load_claims(ROOT)
    selected = set(_DIAGNOSTIC_DISPOSITIONS)
    assert not any(source["path"] in selected for claim in claims["claims"] if claim["role"] != "historical"
                   for source in claim["evidence"])
    assert not selected.intersection(claims["marked_values"]["tables"])
    assert not any(item["path"] in selected for run in registry["runs"].values() if run["lifecycle"] == "current"
                   for item in run["artifacts"])
    consumers = {(run["id"], entry["path"], entry["role"]) for run in registry["runs"].values()
                 if run["lifecycle"] == "current" for entry in run["dependencies"] if entry["path"] in selected}
    assert consumers == {("partial-orders-v2", "outputs/partial_order_v2/run_context.json", "protocol")}


@pytest.mark.parametrize("forgery", ["hash", "baseline", "scope", "role", "assertion", "membership", "alias"])
def test_noncurrent_resigned_ledger_cannot_replace_independent_git_evidence(disposition_git_bytes, forgery):
    ledger = artifacts.read_json(ROOT / "data/artifact_dispositions.json")
    if forgery == "hash":
        ledger["external_references"][0]["original_sha256"] = "0" * 64
    elif forgery == "baseline":
        ledger["baseline_ref"] = "0" * 40
    elif forgery == "scope":
        ledger["current_scientific_support"] = True
    elif forgery == "role":
        ledger["external_references"][0]["role"] = "current_scientific_result"
    elif forgery == "assertion":
        entry = next(item for item in ledger["external_references"] if item["origin_path"] ==
                     "outputs/native_population_development/native_training_01_phase_c_01_scoring/independent_review.json")
        entry["source_assertions"]["/overall_verdict"] = "PASS"
    elif forgery == "alias":
        entry = next(item for item in ledger["external_references"] if "recorded_alias_of" in item)
        entry["recorded_alias_of"] = "outputs/not-the-original/native_diagnostics.json"
    else:
        ledger["external_references"].pop()
    assert len(artifacts.canonical_digest(ledger)) == 64
    with pytest.raises(AssertionError):
        _assert_disposition_evidence(ledger, disposition_git_bytes)


def test_noncurrent_final_reconciliation_preserves_recovered_source_without_promoting_replay(disposition_git_bytes):
    registry = artifacts.load_registry(ROOT)
    run = registry["runs"]["native-reconciliation-final"]
    assert run["lifecycle"] == run["status"] == "investigation_pending"
    assert not run.get("frozen_runtime") and not run.get("verification")
    pins = {item["origin_path"]: item["original_sha256"]
            for item in artifacts.read_json(ROOT / "data/artifact_dispositions.json")["blocked_records"]}
    for item in run["artifacts"]:
        if item["path"] in _RECONCILIATION_DISPOSITION_BLOCKERS:
            assert item["sha256"] == pins[item["path"]]
    final = json.loads(disposition_git_bytes["outputs/native_reconciliation_20260907_01/final/protocol.json"])
    initial = artifacts.read_json(ROOT / "outputs/native_reconciliation_20260907_01/protocol.json")
    assert final["implementation_sha256"]["scripts/run_native_reconciliation.py"] == "73a4aafb649bd0e65ea8368b75f0b5f32fc8e651c37efdb4c5b9769e2a550621"
    assert initial["implementation_sha256"]["scripts/run_native_reconciliation.py"] == "9579fee36be9af4a1e6efd0a8f98c081956e8422a932da2ece27a2ed5382c739"
    assert any(final["implementation_sha256"]["scripts/run_native_reconciliation.py"] in issue for issue in run["issues"])
    assert not run.get("retention")
    assert run["artifact_indexes"][0]["membership_sha256"] == "5a324ab6eca8626b7571849dc72b345c7e81517fd8a6d03808c2b86f8751e8ab"
    reference = run["source_recovery"]
    assert reference["path"] == "data/native_reconciliation_recovery.json"
    assert artifacts.sha256(ROOT / reference["path"]) == reference["sha256"]
    recovery = artifacts.read_json(ROOT / reference["path"])
    assert recovery["run_id"] == run["id"]
    assert {key: recovery["original_runtime"][key] for key in ("python", "packages")} == final["environment"]
    assert recovery["original_runtime"]["matched_before_execution"] is True
    assert artifacts.sha256(ROOT / recovery["authority"]["protocol_path"]) == recovery["authority"]["protocol_sha256"]
    assert {entry["logical_path"] for entry in recovery["recovered_sources"]} == {
        "scripts/run_native_reconciliation.py", "tests/test_native_reconciliation.py"}
    for entry in recovery["recovered_sources"]:
        stored = ROOT / entry["storage_path"]
        assert entry["storage_path"].endswith(".py.source") and entry["storage_path"] != entry["logical_path"]
        assert artifacts.sha256(stored) == entry["sha256"] == final["implementation_sha256"][entry["logical_path"]]
        assert stored.stat().st_size == entry["size_bytes"]
        blob = subprocess.run(["git", "hash-object", "--stdin"], cwd=ROOT, input=stored.read_bytes(),
                              capture_output=True, check=True).stdout.decode().strip()
        assert blob == entry["git_blob"]
    attempt = recovery["isolated_attempt"]
    assert attempt["exit_code"] == 1 and attempt["status"] == "refused_before_numerical_execution"
    assert attempt["differing_training_manifest_fields"] == ["root", "sha256"]
    assert all(attempt[key] is True for key in ("input_records_equal", "code_records_equal",
                                               "source_and_inputs_unchanged", "repository_originals_unchanged",
                                               "frozen_validator_unchanged"))
    assert attempt["numerical_comparison"] == "not_run"
    assert all(recovery[key] is False for key in ("numerical_reproduction_verified",
                                                 "historical_execution_authenticated", "calibration_adopted",
                                                 "biological_validation", "outputs_replaced"))
    assert "--artifact-source" not in run["reproduction"]["argv"]


@pytest.fixture
def diagnostic_disposition_fixture(tmp_path):
    root = tmp_path / "historical"
    (root / "data").mkdir(parents=True)
    (root / "outputs").mkdir()
    output = root / "outputs/diagnostic.json"
    output.write_text('{"status":"refused","biological_validation":false,"value":null}')
    ledger = {"schema_version": 1, "records": [], "external_references": [
        {"origin_path": "outputs/diagnostic.json", "original_sha256": artifacts.sha256(output)}]}
    manifest = root / "data/dispositions.json"
    manifest.write_text(json.dumps(ledger))
    document = {"schema_version": 1, "runs": [], "frozen_bundles": [{
        "id": "diagnostic", "manifest": "data/dispositions.json", "sha256": artifacts.sha256(manifest),
        "include_external_outputs": True, "reason": "Recorded refusal fixture; no current scientific support",
        "reproduction_command": "Original record integrity only; historical execution is not asserted"}]}
    _write_registry(root, document)
    return root, document


def test_noncurrent_integrity_audit_does_not_require_an_invented_executable(diagnostic_disposition_fixture):
    root, _ = diagnostic_disposition_fixture
    before = _snapshot(root)
    rows = artifacts.audit_registry(root)
    assert all(row["status"] == "PASS" for row in rows)
    assert all("not a new fit or a runtime replay" in row["detail"] for row in rows)
    assert _snapshot(root) == before


@pytest.mark.parametrize("change", ["historical_bytes", "refusal_promoted", "evidence_bytes"])
def test_noncurrent_integrity_refuses_altered_bytes_or_forged_evidence(diagnostic_disposition_fixture, change):
    root, _ = diagnostic_disposition_fixture
    if change == "evidence_bytes":
        path = root / "data/dispositions.json"
        payload = json.loads(path.read_text())
        payload["external_references"][0]["original_sha256"] = "0" * 64
        path.write_text(json.dumps(payload))
    else:
        path = root / "outputs/diagnostic.json"
        text = path.read_text()
        path.write_text(text + "\n" if change == "historical_bytes" else text.replace('"refused"', '"passed"'))
    before = _snapshot(root)
    rows = artifacts.audit_registry(root)
    assert any(row["status"] == "FAIL" for row in rows)
    assert _snapshot(root) == before


def test_noncurrent_integrity_cannot_be_promoted_to_current_provenance(diagnostic_disposition_fixture):
    root, document = diagnostic_disposition_fixture
    (root / "scripts").mkdir()
    (root / "scripts/run.py").write_text("raise RuntimeError('not executed by this test')\n")
    document["frozen_bundles"] = []
    document["runs"] = [{"id": "promoted", "lifecycle": "current", "status": "verified",
                         "producer": {"path": "scripts/run.py", "callable": "main", "libraries": []},
                         "reason": "Claimed promotion without execution", "dependencies": [],
                         "parameters": {"original_invocation": "unrecorded"},
                         "model_identity": {"kind": "fixture"}, "runtime": {"kind": "fixture"},
                         "reproduction": {"argv": ["{python}", "-B", "scripts/run.py"]},
                         "artifacts": [{"path": "outputs/diagnostic.json",
                                        "sha256": artifacts.sha256(root / "outputs/diagnostic.json"),
                                        "comparison": {"format": "json", "row_keys": {}, "units": {}}}]}]
    _write_registry(root, document)
    row, = artifacts.audit_registry(root)
    assert row["status"] == "FAIL" and "receipt absent" in row["detail"]


# Declared artifacts that cannot ship, with the reason. Each records the absolute root of
# the machine that produced it, which scripts/audit_claims.py rightly refuses outside the
# path resolver. The defect is in the producer, not the decision to withhold them: a run
# that writes its own checkout path into its output cannot be published from any other
# machine. Re-running native training with a repository-relative root is what empties this.
_UNSHIPPABLE_MACHINE_LOCAL_ARTIFACTS = {
    "outputs/native_training_run_01/checkpoint.json",
    "outputs/native_training_run_01/contract.json",
    "outputs/native_training_run_01/training_manifest.json",
}


def test_every_registered_run_artifact_is_tracked_in_git():
    """A declared artifact that is not committed exists only on the author's machine.

    Until 2026-09-11 nineteen of them were in that state across five runs, and nothing
    caught it: the checks that existed asked whether the path was on disk, which is true
    for whoever produced it and false for every clone. Ask git instead. This reads the
    registry document rather than load_registry, because the loader also expands
    artifact_indexes into their per-case members -- bulk .json.gz that .gitignore excludes
    on purpose as regenerable. Those are exempt by structure; a run's own artifacts list
    is not.
    """
    tracked = set(subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True,
                                 check=True, text=True).stdout.split())
    document = artifacts.read_json(ROOT / "data/artifact_registry.json")
    runs = document["runs"]
    declared = {(run["id"], item["path"])
                for run in (runs.values() if isinstance(runs, dict) else runs)
                for item in (run.get("artifacts") or [])}
    assert declared, "registry declared no artifacts at all"
    untracked = sorted(f"{run_id}: {path}" for run_id, path in declared
                       if path not in tracked and path not in _UNSHIPPABLE_MACHINE_LOCAL_ARTIFACTS)
    assert not untracked, "registered artifacts absent from a fresh clone:\n" + "\n".join(untracked)


def test_the_unshippable_artifacts_are_exactly_the_ones_carrying_a_machine_path():
    """The exemption above must stay earned: shrink when a producer is fixed, never grow quietly."""
    for name in sorted(_UNSHIPPABLE_MACHINE_LOCAL_ARTIFACTS):
        path = ROOT / name
        if not path.exists():
            continue
        assert str(ROOT) in path.read_text(encoding="utf-8"), (
            f"{name} no longer carries a machine-local path and must now be tracked")
