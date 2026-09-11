from copy import deepcopy
from email.parser import BytesParser
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tomllib

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("repository_inventory", ROOT / "scripts/repository_inventory.py")
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


def _write(root, name, content):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _json(root, name, value):
    _write(root, name, json.dumps(value, indent=2) + "\n")


def _blob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


@pytest.fixture
def fixture_repo(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    files = {
        "src/ystwin/__init__.py": "",
        "src/ystwin/runtime.py": "from .core import value\nanswer = value\n",
        "src/ystwin/core.py": "value = 1\n",
        "src/ystwin/extension.py": "__all__ = ['extension']\ndef extension():\n    return 1\n",
        "src/ystwin/test_only.py": "def test_only():\n    return 2\n",
        "scripts/run_checks.py": "from ystwin.core import value\nprint(value)\n",
        "tests/test_validation.py": "from ystwin import extension, test_only\ndef test_values():\n    assert extension.extension() == test_only.test_only() - 1\n",
        "scripts/repository_inventory.py": "def generate():\n    return 1\n",
        "tests/test_repository_inventory.py": "def test_inventory():\n    assert True\n",
        "outputs/negative.json": '{"ready": false}\n',
        "docs/CLAIM.md": "Recorded refusal, not biological readiness.\n",
        "S.csv": "DO_NOT_READ_THIS_STATUS_ONLY_FILE\n",
        "jws.html": "DO_NOT_READ_THIS_STATUS_ONLY_FILE\n",
        ".ystwin.env": "DO_NOT_READ_THIS_SECRET\n",
    }
    for name, content in files.items():
        _write(root, name, content)
    seed = {
        "schema_version": 1,
        "kind": "repository_inventory",
        "classification_policy": {
            "roles": list(inventory.ROLES),
            "active_roots": ["src/ystwin/runtime.py"],
            "supported_public_apis": [{"path": "src/ystwin/extension.py", "work_item": "W14", "purpose": "Public extension needed by the finished architecture, not justified by its test."}],
            "explicit_path_roles": [
                {"path": "S.csv", "role": "unowned_requires_decision", "work_item": "W16", "purpose": "Unknown root CSV, status only"},
                {"path": "jws.html", "role": "unowned_requires_decision", "work_item": "W16", "purpose": "Unknown local HTML, status only"},
            ],
            "reviewed_data_domains": [],
            "do_not_read": list(inventory.STATUS_ONLY),
        },
        "retained_negative_registry": [{
            "id": "N01", "claim": "Readiness refused", "scope": "Fixture software status only",
            "artifact_paths": ["outputs/negative.json"], "artifact_prefixes": [],
            "supporting_code": ["src/ystwin/runtime.py"],
            "supporting_tests": ["tests/test_validation.py"],
            "claim_documents": ["docs/CLAIM.md"], "preservation": "Keep every failed result",
        }],
        "decisions": [{"id": "D01", "status": "open", "owner": "maintainer", "blocks": ["W16"], "paths": ["S.csv"], "question": "Resolve provenance without deletion"}],
    }
    nodes = []
    for identity in ("W00", "W01", "W09", "W13", "W14", "W15", "W16"):
        nodes.append({
            "id": identity, "deliverable": f"Deliverable {identity}", "owner": "fixture-owner",
            "status": "not_started", "depends_on": [] if identity == "W00" else ["W00"],
            "existing_modules": [], "existing_paths": [],
            "tests": ["tests/test_validation.py"],
            "completion_conditions": ["Actual acceptance evidence required"], "acceptance_evidence": [],
        })
    nodes[0]["existing_paths"] = list(inventory.OWNED_PATHS) + ["src/ystwin/runtime.py", "src/ystwin/core.py", "scripts/run_checks.py"]
    workflow = {"schema_version": 1, "scope_complete": False, "baseline": {"commit": "0" * 40}, "nodes": nodes}
    deletions = {"schema_version": 1, "approval_required": True, "execution_allowed": False, "proposals": [{
        "path": "S.csv", "approval_required": True, "approved": False, "owner": "maintainer",
        "proposed_action": "retain_pending_decision", "reason": "Unknown source", "replacement_paths": [],
        "replacement_contract": "None identified", "approval_conditions": ["Identify provenance first"],
    }]}
    for name, value in zip(inventory.DOCUMENTS, (seed, workflow, deletions)):
        _json(root, f"data/{name}", value)
    tracked = {name: {"mode": "100644", "object_id": _blob((root / name).read_bytes())} for name in set(files) | {"data/" + name for name in inventory.DOCUMENTS} if name not in {"jws.html", ".ystwin.env"}}
    original = {name: {**metadata, "kind": "blob"} for name, metadata in tracked.items()}

    def catalog(_root, baseline):
        return deepcopy(tracked), deepcopy(original), {
            "baseline_commit": baseline, "baseline_tree_object_id": "1" * 40,
            "tracked_pathset_sha256": inventory.digest(sorted(tracked)),
        }

    def status(_root, name, names):
        return {"path": name, "tracked": name in names, "exists": (_root / name).exists(), "symlink": False,
                "ignored": name == "jws.html", "status": "tracked" if name in names else "untracked_ignored", "contents_read": False}

    monkeypatch.setattr(inventory, "git_catalog", catalog)
    monkeypatch.setattr(inventory, "path_status", status)
    return root, seed, workflow, deletions, tracked


def test_actual_git_catalog_covers_every_tracked_path_without_opening_contents():
    baseline = "84561f51502b5b4bf66c052dc88c853038917934"
    current, original, metadata = inventory.git_catalog(ROOT, baseline)
    authoritative = inventory.git(ROOT, "ls-files", "-z").decode().split("\0")
    assert set(current) == set(authoritative) - {""}
    assert "S.csv" in original
    assert "jws.html" not in current
    assert metadata["baseline_commit"] == baseline
    assert all(len(item["object_id"]) in (40, 64) for item in original.values())


def test_every_path_has_one_role_and_test_only_is_not_runtime(fixture_repo):
    root, seed, workflow, deletions, tracked = fixture_repo
    result, _, _ = inventory.build_inventory(root, seed, workflow, deletions)
    entries = {entry["path"]: entry for entry in result["entries"]}
    assert set(tracked) <= set(entries)
    assert len(entries) == len(result["entries"])
    assert all(entry["role"] in inventory.ROLES and entry["purpose"] and entry["work_items"] for entry in entries.values())
    assert entries["src/ystwin/core.py"]["role"] == "active_runtime"
    assert entries["src/ystwin/extension.py"]["role"] == "supported_public_library"
    assert entries["src/ystwin/test_only.py"]["role"] == "unowned_requires_decision"
    assert entries["src/ystwin/test_only.py"]["imported_by"]["test"] == ["tests/test_validation.py"]
    assert entries["outputs/negative.json"]["role"] == "recorded_negative_result"
    assert sum(result["summary"]["tracked_by_role"].values()) == len(tracked)


def test_status_only_artifacts_and_secret_env_are_never_opened(fixture_repo, monkeypatch):
    root, seed, workflow, deletions, _ = fixture_repo
    original = Path.read_bytes

    def guard(path):
        if path.name in inventory.STATUS_ONLY:
            raise AssertionError("Status-only or secret content was opened")
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guard)
    result, _, _ = inventory.build_inventory(root, seed, workflow, deletions)
    status = {item["path"]: item for item in result["special_root_artifacts"]}
    assert status["S.csv"]["tracked"]
    assert status["jws.html"]["status"] == "untracked_ignored"
    assert not any(item["contents_read"] for item in status.values())
    assert result["scan"]["secrets_read"] is False


def test_relative_deferred_and_literal_dynamic_imports_report_actual_name_uses():
    sources = {
        "src/ystwin/__init__.py": "",
        "src/ystwin/core.py": "value = 1\n",
        "src/ystwin/driver.py": "from .core import value as local\nanswer = local\ndef load():\n    from . import core\n    return core.value\nimport importlib\nimportlib.import_module('ystwin.core')\n",
    }
    graph = inventory.analyze_python(sources, set(sources))
    edges = [edge for edge in graph["imports"] if edge["target"] == "src/ystwin/core.py"]
    assert any(edge["bound_names"] == ["local"] and edge["name_use_lines"] == [2] for edge in edges)
    assert any(edge["deferred"] and edge["name_use_lines"] == [5] for edge in edges)
    assert any(edge["kind"] == "literal_dynamic_import" for edge in edges)
    assert not graph["dynamic_resolution_gaps"]


def test_computed_dynamic_import_is_a_gap_not_a_deleted_module():
    sources = {"scripts/check.py": "import importlib\nname = input()\nimportlib.import_module(name)\n"}
    graph = inventory.analyze_python(sources, set(sources))
    assert graph["dynamic_resolution_gaps"][0]["kind"] == "dynamic_module_import"
    assert graph["imports"] == []


def test_type_checking_edges_do_not_establish_runtime_reachability():
    sources = {"src/ystwin/a.py": "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from .b import b\n",
               "src/ystwin/b.py": "b = 1\n"}
    graph = inventory.analyze_python(sources, set(sources))
    assert inventory.reachable(["src/ystwin/a.py"], graph["imports"]) == {"src/ystwin/a.py"}


def test_import_graph_may_cycle_but_work_dag_cannot(fixture_repo):
    sources = {"src/ystwin/a.py": "from . import b\n", "src/ystwin/b.py": "def use():\n    from . import a\n    return a\n"}
    graph = inventory.analyze_python(sources, set(sources))
    assert inventory.strongly_connected(graph["imports"]) == [["src/ystwin/a.py", "src/ystwin/b.py"]]
    root, seed, workflow, deletions, _ = fixture_repo
    workflow["nodes"][0]["depends_on"] = ["W15"]
    with pytest.raises(inventory.InventoryError, match="cycle"):
        inventory.build_inventory(root, seed, workflow, deletions)


def test_implementation_cannot_be_complete_from_existing_files(fixture_repo):
    root, seed, workflow, deletions, _ = fixture_repo
    workflow["nodes"][0]["status"] = "completed"
    with pytest.raises(inventory.InventoryError, match="acceptance evidence"):
        inventory.build_inventory(root, seed, workflow, deletions)


def test_negative_registry_requires_real_support_and_does_not_delete(fixture_repo):
    root, seed, workflow, deletions, _ = fixture_repo
    result, _, manifest = inventory.build_inventory(root, seed, workflow, deletions)
    negative = result["retained_negative_registry"][0]
    assert negative["retained_artifact_paths"] == ["outputs/negative.json"]
    assert negative["tests_are_scoped_software_checks_not_biological_authorization"]
    assert "outputs/negative.json" in manifest["protected_exact_paths"]
    seed["retained_negative_registry"][0]["supporting_tests"] = ["tests/not_real.py"]
    with pytest.raises(inventory.InventoryError, match="missing path"):
        inventory.build_inventory(root, seed, workflow, deletions)


@pytest.mark.parametrize("mutation", ["approval", "execution", "protected", "wildcard"])
def test_deletion_manifest_cannot_authorize_or_target_evidence(fixture_repo, mutation):
    root, seed, workflow, deletions, _ = fixture_repo
    if mutation == "approval":
        deletions["proposals"][0]["approved"] = True
    elif mutation == "execution":
        deletions["execution_allowed"] = True
    elif mutation == "protected":
        deletions["proposals"][0]["path"] = "outputs/negative.json"
    else:
        deletions["proposals"][0]["path"] = "outputs/*.json"
    with pytest.raises(inventory.InventoryError):
        inventory.build_inventory(root, seed, workflow, deletions)
    assert (root / "outputs/negative.json").is_file()


def test_static_writer_uses_ast_destination_not_filename_resemblance():
    sources = {"scripts/writer.py": "from ystwin import paths\nOUT = paths.outputs_dir()\nframe.to_csv(OUT / 'exact.csv')\nmentioned = 'outputs/not_written.csv'\n"}
    graph = inventory.analyze_python(sources, set(sources) | {"outputs/exact.csv", "outputs/not_written.csv"})
    assert {edge["target"] for edge in graph["static_writers"]} == {"outputs/exact.csv"}
    assert any(edge["target"] == "outputs/not_written.csv" for edge in graph["path_references"])


def test_check_is_read_only_and_does_not_invalidate_on_ordinary_body_edits(fixture_repo, monkeypatch):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    paths = [root / "data" / name for name in inventory.DOCUMENTS]
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}
    _write(root, "src/ystwin/core.py", "value = 2\n")

    def refuse_write(*args, **kwargs):
        raise AssertionError("--check attempted a write")

    monkeypatch.setattr(Path, "write_text", refuse_write)
    monkeypatch.setattr(Path, "write_bytes", refuse_write)
    monkeypatch.setattr(Path, "mkdir", refuse_write)
    assert inventory.run(["--root", str(root), "--check"]) == 0
    assert before == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


def test_check_detects_changed_import_topology(fixture_repo):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    _write(root, "src/ystwin/core.py", "from .test_only import test_only\nvalue = test_only()\n")
    assert inventory.run(["--root", str(root), "--check"]) == 1


def test_output_directory_isolated_and_exact_snapshot_hashes_are_opt_in(fixture_repo, tmp_path, capsys):
    root, _, _, _, _ = fixture_repo
    original = {name: (root / "data" / name).read_bytes() for name in inventory.DOCUMENTS}
    output = tmp_path / "final"
    assert inventory.run(["--root", str(root), "--output-dir", str(output), "--snapshot-hashes"]) == 0
    assert original == {name: (root / "data" / name).read_bytes() for name in inventory.DOCUMENTS}
    assert inventory.run(["--root", str(root), "--output-dir", str(output), "--check", "--verify-source-hashes"]) == 0
    _write(root, "src/ystwin/core.py", "value = 8\n")
    assert inventory.run(["--root", str(root), "--output-dir", str(output), "--check"]) == 0
    capsys.readouterr()
    assert inventory.run(["--root", str(root), "--output-dir", str(output), "--check", "--verify-source-hashes"]) == 1
    assert "src/ystwin/core.py" in json.loads(capsys.readouterr().out)["errors"][0]


def test_unknown_new_tracked_path_is_explicitly_unresolved(fixture_repo):
    root, seed, workflow, deletions, tracked = fixture_repo
    _write(root, "new_artifact.bin", "unclassified\n")
    tracked["new_artifact.bin"] = {"mode": "100644", "object_id": "2" * 40}
    result, _, _ = inventory.build_inventory(root, seed, workflow, deletions)
    assert "new_artifact.bin" in result["summary"]["unresolved_role_paths"]


def test_untracked_source_is_not_read_without_explicit_inclusion(fixture_repo):
    root, seed, workflow, deletions, _ = fixture_repo
    _write(root, "src/ystwin/new_engine.py", "this is not valid python !!!\n")
    result, _, _ = inventory.build_inventory(root, seed, workflow, deletions)
    assert "src/ystwin/new_engine.py" not in {entry["path"] for entry in result["entries"]}
    with pytest.raises(inventory.InventoryError, match="invalid Python"):
        inventory.build_inventory(root, seed, workflow, deletions, ["src/ystwin/new_engine.py"])


def test_symlink_source_cannot_escape_read_allowlist(fixture_repo, tmp_path):
    root, seed, workflow, deletions, tracked = fixture_repo
    outside = tmp_path / "private.py"
    outside.write_text("secret_value = 'do not read'\n")
    (root / "src/ystwin/link.py").symlink_to(outside)
    tracked["src/ystwin/link.py"] = {"mode": "120000", "object_id": "3" * 40}
    with pytest.raises(inventory.InventoryError, match="symlinked"):
        inventory.build_inventory(root, seed, workflow, deletions)


@pytest.mark.parametrize("name", ["../secret", "/absolute", "data/../secret", "a\\b", "a//b", "./x"])
def test_unsafe_relative_paths_are_rejected(name):
    with pytest.raises(inventory.InventoryError):
        inventory.relative_path(name)


def test_duplicate_json_keys_and_nonfinite_values_are_rejected():
    with pytest.raises(inventory.InventoryError):
        inventory.parse_json(b'{"a": 1, "a": 2}')
    with pytest.raises(inventory.InventoryError):
        inventory.parse_json(b'{"a": NaN}')


def test_real_retained_bundle_is_verified_without_importing_runtime_or_opening_originals():
    baseline = "84561f51502b5b4bf66c052dc88c853038917934"
    current, original, _ = inventory.git_catalog(ROOT, baseline)
    document = inventory.parse_json((ROOT / "data/deletion_manifest.json").read_bytes())
    cleanup = inventory.validate_cleanup(ROOT, document, current, original)
    evidence = inventory.frozen_evidence(ROOT, set(current), original, cleanup["relocations"])
    assert evidence["manifest"] == inventory.FROZEN_MANIFEST
    assert evidence["records"]
    assert all(check["verified"] for check in evidence["integrity_checks"])
    assert evidence["references"]["outputs/native_training_run_02/checkpoint.json"]


def test_inventory_cli_uses_no_package_or_local_environment_loader():
    result = subprocess.run([sys.executable, "-B", str(ROOT / "scripts/repository_inventory.py"), "--help"], capture_output=True, text=True, check=False)
    assert result.returncode == 0
    assert "--output-dir" in result.stdout and "--check" in result.stdout
    assert "--verify-source-hashes" in result.stdout


def test_serialization_is_readable_lossless_finite_json():
    value = {"small": [1, 2, 3], "records": [{"name": "a", "passed": False}], "nested": {"x": {"long": "x" * 200}}}
    encoded = inventory.serialized(value)
    assert json.loads(encoded) == value
    assert "[1, 2, 3]" in encoded
    with pytest.raises(ValueError):
        inventory.serialized({"not_finite": float("nan")})


def test_type_checking_else_branch_retains_runtime_import():
    sources = {
        "src/ystwin/a.py": "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from .b import value\nelse:\n    from .c import value\n",
        "src/ystwin/b.py": "value = 1\n", "src/ystwin/c.py": "value = 2\n",
    }
    graph = inventory.analyze_python(sources, set(sources))
    assert inventory.reachable(["src/ystwin/a.py"], graph["imports"]) == {"src/ystwin/a.py", "src/ystwin/c.py"}


def test_literal_external_import_is_not_an_unresolved_repository_dependency():
    sources = {"scripts/check.py": "import importlib\nimportlib.import_module('numpy')\n"}
    graph = inventory.analyze_python(sources, set(sources))
    assert "numpy" in graph["python_files"]["scripts/check.py"]["external_imports"]
    assert graph["dynamic_resolution_gaps"] == []


def test_structured_module_use_comes_from_verified_column_and_consumer(tmp_path):
    table = "data/modules.tsv"
    _write(tmp_path, table, "code_module\tverdict\ncore\tIN_MODULE\n-\tNOT_METABOLIC\n")
    sources = {
        "src/ystwin/consumer.py": "import importlib\ndef check(record):\n    importlib.import_module(f'ystwin.{record.code_module}')\n",
        "src/ystwin/core.py": "value = 1\n",
    }
    contract = {"table": table, "consumer": "src/ystwin/consumer.py", "function": "check", "column": "code_module", "module_prefix": "ystwin."}
    hashes = {}
    edges = inventory.structured_imports(tmp_path, [contract], sources, set(sources) | {table}, hashes)
    assert len(edges) == 1
    assert edges[0]["target"] == "src/ystwin/core.py" and edges[0]["registry_rows"] == [2]
    assert inventory.reachable(["src/ystwin/consumer.py"], edges, include_validation_imports=False) == {"src/ystwin/consumer.py"}
    assert "src/ystwin/core.py" in inventory.reachable(["src/ystwin/consumer.py"], edges)
    assert table in hashes
    sources["src/ystwin/consumer.py"] = "def check(record):\n    return record\n"
    with pytest.raises(inventory.InventoryError, match="consumer"):
        inventory.structured_imports(tmp_path, [contract], sources, set(sources) | {table}, {})


def test_historical_executor_hashes_do_not_follow_current_runtime(tmp_path):
    name = "outputs/frozen/code/executor.py"
    _write(tmp_path, name, "value = 1\n")
    baseline = {name: {"kind": "blob", "object_id": _blob((tmp_path / name).read_bytes())}}
    assert inventory.frozen_evidence(tmp_path, {name}, baseline)["integrity_checks"][0]["verified"]
    _write(tmp_path, name, "value = 2\n")
    with pytest.raises(inventory.InventoryError, match="changed from baseline"):
        inventory.frozen_evidence(tmp_path, {name}, baseline)
    with pytest.raises(inventory.InventoryError, match="no longer inventoried"):
        inventory.frozen_evidence(tmp_path, set(), baseline)


def test_frozen_manifest_tampering_is_refused_without_opening_originals(monkeypatch):
    monkeypatch.setattr(inventory, "read_public", lambda root, name: b"{}")
    with pytest.raises(inventory.InventoryError, match="independently pinned"):
        inventory.frozen_evidence(ROOT, {inventory.FROZEN_MANIFEST}, {})


@pytest.mark.parametrize("expected", [True, 0, "false"])
def test_frozen_claim_assertions_refuse_false_readiness_and_type_substitution(tmp_path, expected):
    name = "data/frozen_evidence/fixture.json"
    _json(tmp_path, name, {"payload": {"ready": False}})
    frozen = {"records": {name: {"public_sha256": hashlib.sha256((tmp_path / name).read_bytes()).hexdigest()}}}
    correct = {"path": name, "pointer": "/payload/ready", "equals": False}
    assert inventory.check_frozen_claims(tmp_path, [correct], frozen)[0]["verified"]
    with pytest.raises(inventory.InventoryError, match="contradicts"):
        inventory.check_frozen_claims(tmp_path, [{**correct, "equals": expected}], frozen)


def test_real_negative_assertions_are_bound_to_pinned_artifacts():
    seed = inventory.parse_json((ROOT / "data/repository_inventory.json").read_bytes())
    current, original, _ = inventory.git_catalog(ROOT, "84561f51502b5b4bf66c052dc88c853038917934")
    document = inventory.parse_json((ROOT / "data/deletion_manifest.json").read_bytes())
    cleanup = inventory.validate_cleanup(ROOT, document, current, original)
    frozen = inventory.frozen_evidence(ROOT, set(current), original, cleanup["relocations"])
    assertions = [assertion for record in seed["retained_negative_registry"] for assertion in record.get("frozen_assertions", [])]
    assert len(assertions) >= 14
    checks = inventory.check_frozen_claims(ROOT, assertions, frozen)
    assert all(check["verified"] for check in checks)
    assert any(check["pointer"] == "/payload/final_test_defined" and check["equals"] is False for check in checks)
    assert any(check["pointer"] == "/payload/ready_for_biological_law_claim" and check["equals"] is False for check in checks)


def test_summary_tampering_does_not_pass_check(fixture_repo):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    path = root / "data/repository_inventory.json"
    document = inventory.parse_json(path.read_bytes())
    document["summary"]["architecture_complete"] = True
    _json(root, "data/repository_inventory.json", document)
    assert inventory.run(["--root", str(root), "--check"]) == 1


def test_deletion_caller_changes_require_review_even_without_import_changes(fixture_repo):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    _write(root, "scripts/run_checks.py", "from ystwin.core import value\nprint(value)\npath = 'S.csv'\n")
    assert inventory.run(["--root", str(root), "--check"]) == 1


def test_reviewed_session_artifact_inclusion_is_metadata_only(fixture_repo, monkeypatch):
    root, seed, workflow, deletions, _ = fixture_repo
    name = "data/new_public_registry.json"
    _write(root, name, "Not needed for role/path inventory; do not open.\n")
    workflow["reviewed_session_additions"] = [{"path": name, "work_item": "W13"}]
    original = Path.read_bytes

    def guard(path):
        if path == root / name:
            raise AssertionError("Artifact metadata inclusion read its content")
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guard)
    result, _, _ = inventory.build_inventory(root, seed, workflow, deletions)
    entry = next(item for item in result["entries"] if item["path"] == name)
    assert not entry["tracked"] and entry["role"] == "validation_reproduction"
    assert entry["work_items"] == ["W13"]


def test_clean_clone_need_not_have_the_local_ignored_html(fixture_repo, monkeypatch):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    original = inventory.path_status

    def absent_html(repo, name, tracked):
        result = original(repo, name, tracked)
        if name == "jws.html":
            result["exists"] = False
        return result

    monkeypatch.setattr(inventory, "path_status", absent_html)
    assert inventory.run(["--root", str(root), "--check"]) == 0


def test_require_resolved_never_mistakes_consistency_for_completion(fixture_repo):
    root, _, _, _, _ = fixture_repo
    assert inventory.run(["--root", str(root)]) == 0
    assert inventory.run(["--root", str(root), "--check", "--require-resolved"]) == 2


def _approved_cleanup_manifest(completed):
    # One request per approval batch: a contract is only valid inside the batch that authorised
    # it, so a fixture that pooled every action into one request would not exercise the contract.
    proposals, by_batch = [], {}
    for name, contract in inventory.CLEANUP_CONTRACTS.items():
        batch = contract["batch"]
        group = by_batch.setdefault(batch, {"actions": [], "receipts": [], "completed": []})
        action = {"path": name, "action": contract["action"], "replacement_paths": contract["replacement_paths"]}
        if name == inventory.RELEASE_SOURCE:
            action.update(replacement=inventory.RELEASE_ASSET, sha256=contract["original"]["sha256"])
        group["actions"].append(action)
        proposals.append({
            "path": name, "proposed_action": contract["action"], "approval_required": True,
            "approved": True, "action_executed": name in completed, "owner": "fixture-reviewer",
            "approval_request_id": batch, "reason": "Exact reviewed fixture",
            "replacement_paths": contract["replacement_paths"], "replacement_contract": "Reviewed replacement role",
            "approval_conditions": ["Actual user approval and verified replacements"],
        })
        if name in completed:
            group["completed"].append(name)
            group["receipts"].append({
                "path": name, "action": contract["action"], "approval_request_id": batch,
                "original": contract["original"], "replacement_paths": contract["replacement_paths"],
                "command": inventory.cleanup_command(name), "returncode": 0,
                "operator": "fixture-reviewer", "date": inventory.CLEANUP_APPROVALS[batch]["date"],
            })

    def request(batch):
        group = by_batch[batch]
        return {"id": batch, "status": "approved", "approval": inventory.CLEANUP_APPROVALS[batch],
                "actions": group["actions"], "executed_paths": sorted(group["completed"]),
                "execution_receipts": group["receipts"]}

    extra = [batch for batch in by_batch if batch != inventory.CLEANUP_APPROVAL_ID]
    return deepcopy({
        "schema_version": 1, "approval_required": True, "execution_allowed": False, "proposals": proposals,
        "approval_request": request(inventory.CLEANUP_APPROVAL_ID),
        "additional_approval_requests": [request(batch) for batch in sorted(extra)],
    })


@pytest.fixture
def completed_cleanup_repo(tmp_path, monkeypatch):
    root = tmp_path.resolve() / "reviewed_cleanup"
    root.mkdir()
    completed = set(inventory.CLEANUP_CONTRACTS) - {inventory.TASK_TABLE}
    document = _approved_cleanup_manifest(completed)
    baseline, tracked, blobs = {}, {}, {}
    for name, contract in inventory.CLEANUP_CONTRACTS.items():
        identity = contract["original"]
        if identity["tracked"]:
            baseline[name] = {"mode": "100644", "kind": "blob", "object_id": identity["git_object_id"]}
            blobs[f"{identity['git_commit']}:{name}"] = inventory.git(ROOT, "show", f"{identity['git_commit']}:{name}")
        if name not in completed:
            _write(root, name, blobs[f"{identity['git_commit']}:{name}"].decode())
            tracked[name] = {key: baseline[name][key] for key in ("mode", "object_id")}
        else:
            for replacement in contract["replacement_paths"]:
                _write(root, replacement, "{}\n")
    identity = inventory.CLEANUP_CONTRACTS[inventory.RELEASE_SOURCE]["original"]
    (root / inventory.RELEASE_ASSET).write_bytes(blobs[f"{identity['git_commit']}:{inventory.RELEASE_SOURCE}"])
    tracked[inventory.RELEASE_ASSET] = {"mode": "100644", "object_id": identity["git_object_id"]}

    def read_git(_root, *arguments, **kwargs):
        assert arguments[0] == "show", "Cleanup validation attempted a non-read Git operation"
        return blobs[arguments[1]]

    monkeypatch.setattr(inventory, "git", read_git)
    return root, document, tracked, baseline


def test_completed_approved_actions_are_verified_without_execution(completed_cleanup_repo, monkeypatch):
    root, document, tracked, baseline = completed_cleanup_repo

    def no_mutation(*args, **kwargs):
        raise AssertionError("Inventory attempted to execute cleanup")

    for method in ("unlink", "rename", "replace", "write_text", "write_bytes", "mkdir"):
        monkeypatch.setattr(Path, method, no_mutation)
    result = inventory.validate_cleanup(root, document, tracked, baseline)
    assert len(result["completed"]) == len(inventory.CLEANUP_CONTRACTS) - 1
    assert result["relocations"] == {inventory.RELEASE_SOURCE: inventory.RELEASE_ASSET}
    assert result["replacements"] == {}
    assert set(result["approved"]) - set(result["completed"]) == {inventory.TASK_TABLE}
    assert (root / inventory.TASK_TABLE).is_file()


@pytest.mark.parametrize("mutation", [
    "request_id", "approval_status", "actor", "selection", "missing_action", "extra_action",
    "path", "action", "replacement", "source_sha", "executed_paths", "duplicate_receipt",
    "receipt_path", "receipt_approval", "receipt_hash", "receipt_size", "receipt_original_git",
    "receipt_replacement", "receipt_command", "receipt_exit", "boolean_exit",
])
def test_forged_or_mismatched_approvals_and_receipts_are_rejected(completed_cleanup_repo, mutation):
    root, document, tracked, baseline = completed_cleanup_repo
    request = document["approval_request"]
    action = request["actions"][0]
    receipt = request["execution_receipts"][0]
    if mutation == "request_id":
        request["id"] = "not-approved"
    elif mutation == "approval_status":
        request["status"] = "pending"
    elif mutation in {"actor", "selection"}:
        request["approval"][mutation] = "forged"
    elif mutation == "missing_action":
        request["actions"].pop()
    elif mutation == "extra_action":
        request["actions"].append({**action, "path": "outputs/negative.json"})
    elif mutation == "path":
        action["path"] = "src/ystwin.egg-info/*"
    elif mutation == "action":
        action["action"] = "remove_directory"
    elif mutation == "replacement":
        action["replacement_paths"] = ["README.md"]
    elif mutation == "source_sha":
        request["actions"][-1]["sha256"] = "0" * 64
    elif mutation == "executed_paths":
        request["executed_paths"].append(inventory.TASK_TABLE)
    elif mutation == "duplicate_receipt":
        request["execution_receipts"].append(deepcopy(receipt))
    elif mutation == "receipt_path":
        receipt["path"] = "S-other.csv"
    elif mutation == "receipt_approval":
        receipt["approval_request_id"] = "not-approved"
    elif mutation == "receipt_hash":
        receipt["original"]["sha256"] = "0" * 64
    elif mutation == "receipt_size":
        receipt["original"]["size_bytes"] += 1
    elif mutation == "receipt_original_git":
        receipt["original"]["git_commit"] = "0" * 40
    elif mutation == "receipt_replacement":
        receipt["replacement_paths"] = ["README.md"]
    elif mutation == "receipt_command":
        receipt["command"] = ["git", "clean", "-fdx"]
    elif mutation == "receipt_exit":
        receipt["returncode"] = 1
    else:
        receipt["returncode"] = False
    with pytest.raises(inventory.InventoryError):
        inventory.validate_cleanup(root, document, tracked, baseline)


@pytest.mark.parametrize("mutation", ["present", "tracked", "directory", "symlink", "missing_replacement", "replacement_symlink", "baseline"])
def test_completed_cleanup_checks_actual_paths_not_receipt_booleans(completed_cleanup_repo, mutation):
    root, document, tracked, baseline = completed_cleanup_repo
    original = "src/ystwin.egg-info/PKG-INFO"
    path = root / original
    if mutation == "present":
        _write(root, original, "reappeared\n")
    elif mutation == "tracked":
        tracked[original] = baseline[original]
    elif mutation in {"directory", "symlink"}:
        path.parent.mkdir(parents=True)
        if mutation == "directory":
            path.mkdir()
        else:
            path.symlink_to(root / "pyproject.toml")
    elif mutation == "missing_replacement":
        (root / "pyproject.toml").unlink()
    elif mutation == "replacement_symlink":
        (root / "pyproject.toml").unlink()
        (root / "pyproject.toml").symlink_to(root / "requirements-quality.txt")
    else:
        baseline[original]["object_id"] = "0" * 40
    with pytest.raises(inventory.InventoryError):
        inventory.validate_cleanup(root, document, tracked, baseline)


@pytest.mark.parametrize("mutation", ["bytes", "same_size", "missing", "symlink", "untracked", "git_blob", "git_mode", "both_paths"])
def test_relocated_source_cannot_be_changed_aliased_or_rebound(completed_cleanup_repo, mutation):
    root, document, tracked, baseline = completed_cleanup_repo
    asset = root / inventory.RELEASE_ASSET
    if mutation in {"bytes", "same_size"}:
        data = asset.read_bytes()
        asset.write_bytes(b"!" + data[1:] if mutation == "same_size" else data + b"\n")
    elif mutation in {"missing", "symlink"}:
        asset.unlink()
        if mutation == "symlink":
            asset.symlink_to(root / inventory.TASK_TABLE)
    elif mutation == "untracked":
        tracked.pop(inventory.RELEASE_ASSET)
    elif mutation == "git_blob":
        tracked[inventory.RELEASE_ASSET]["object_id"] = "0" * 40
    elif mutation == "git_mode":
        tracked[inventory.RELEASE_ASSET]["mode"] = "120000"
    else:
        (root / inventory.RELEASE_SOURCE).write_bytes(asset.read_bytes())
    with pytest.raises(inventory.InventoryError):
        inventory.validate_cleanup(root, document, tracked, baseline)


def test_relocation_preserves_frozen_and_negative_logical_identities(completed_cleanup_repo):
    root, document, tracked, baseline = completed_cleanup_repo
    cleanup = inventory.validate_cleanup(root, document, tracked, baseline)
    paths = set(tracked)
    frozen = inventory.frozen_evidence(root, paths, baseline, cleanup["relocations"])
    check = frozen["integrity_checks"][0]
    assert check["path"] == inventory.RELEASE_ASSET
    assert check["logical_source_path"] == inventory.RELEASE_SOURCE
    assert check["sha256"] == inventory.CLEANUP_CONTRACTS[inventory.RELEASE_SOURCE]["original"]["sha256"]
    registry = [{
        "id": "N01", "claim": "Historical refusal", "scope": "Retained original source",
        "artifact_paths": [inventory.RELEASE_SOURCE], "artifact_prefixes": [],
        "supporting_code": [], "supporting_tests": [inventory.TASK_TABLE],
        "claim_documents": [inventory.TASK_TABLE], "preservation": "Original identity remains frozen",
    }]
    result, membership = inventory.expand_negative_registry(registry, paths, set(baseline), frozen, root, cleanup["relocations"])
    assert result[0]["artifact_paths"] == [inventory.RELEASE_SOURCE]
    assert result[0]["retained_logical_artifact_paths"] == [inventory.RELEASE_SOURCE]
    assert result[0]["retained_artifact_paths"] == [inventory.RELEASE_ASSET]
    assert membership[inventory.RELEASE_ASSET] == {"N01"}
    with pytest.raises(inventory.InventoryError, match="no longer inventoried"):
        inventory.frozen_evidence(root, paths, baseline)
    with pytest.raises(inventory.InventoryError, match="Unreviewed"):
        inventory.frozen_evidence(root, paths, baseline, {"outputs/other.py": inventory.RELEASE_ASSET})


@pytest.mark.parametrize("name", ["S.csv", "src/ystwin/core.py", "outputs/negative.json"])
def test_unapproved_baseline_untracking_is_never_exempt(fixture_repo, name):
    root, seed, workflow, deletions, tracked = fixture_repo
    tracked.pop(name)
    with pytest.raises(inventory.InventoryError, match="Baseline paths lost"):
        inventory.build_inventory(root, seed, workflow, deletions)
    assert (root / name).exists()


@pytest.mark.parametrize("command", ["rm", "mv", "clean", "reset", "add", "commit"])
def test_inventory_git_wrapper_has_no_mutation_capability(command, monkeypatch):
    def refuse_subprocess(*args, **kwargs):
        raise AssertionError("Mutation reached subprocess")

    monkeypatch.setattr(subprocess, "run", refuse_subprocess)
    with pytest.raises(inventory.InventoryError, match="read-only"):
        inventory.git(ROOT, command, "--", "S.csv")


@pytest.mark.parametrize("mutation", [None, "approved_flag", "executed_flag", "approval_id", "action", "replacement", "missing_proposal"])
def test_proposals_must_match_verified_completed_actions(completed_cleanup_repo, mutation):
    root, document, tracked, baseline = completed_cleanup_repo
    cleanup = inventory.validate_cleanup(root, document, tracked, baseline)
    proposal = document["proposals"][0]
    if mutation == "approved_flag":
        proposal["approved"] = False
    elif mutation == "executed_flag":
        proposal["action_executed"] = False
    elif mutation == "approval_id":
        proposal["approval_request_id"] = "forged"
    elif mutation == "action":
        proposal["proposed_action"] = "remove_directory"
    elif mutation == "replacement":
        proposal["replacement_paths"] = ["README.md"]
    elif mutation == "missing_proposal":
        document["proposals"].pop()
    graph = {"imports": [], "path_references": [], "document_references": []}
    protected = {inventory.RELEASE_SOURCE, inventory.RELEASE_ASSET, inventory.TASK_TABLE}
    arguments = (document, set(tracked), protected, graph, {"references": {}}, {}, cleanup)
    if mutation is not None:
        with pytest.raises(inventory.InventoryError):
            inventory.validate_deletions(*arguments)
    else:
        result = inventory.validate_deletions(*arguments)
        assert sum(proposal["action_executed"] for proposal in result["proposals"]) == len(inventory.CLEANUP_CONTRACTS) - 1
        assert result["approval_request"]["pending_paths"] == [inventory.TASK_TABLE]
        assert len(result["execution_verification"]) == len(inventory.CLEANUP_CONTRACTS) - 1


@pytest.mark.parametrize("mutation", [None, "missing_biomass", "reduced_model", "success", "biology", "unbound_hash", "old_claim", "unregistered", "missing_report"])
def test_task_replacement_requires_the_full_refusal_and_migrated_registries(tmp_path, mutation):
    root = tmp_path.resolve()
    tasks = ["biomass", "ethanol", "CO2"]
    model = {"tasks": tasks, "status": "infeasible", "refusal": "full task set infeasible",
             "order": [], "constrained_growth": None, "biological_validation": False}
    report = {"tasks": tasks, "status": "completed", "biological_validation": False,
              "models": {"plain": deepcopy(model), "ec": deepcopy(model)}}
    if mutation == "missing_biomass":
        report["tasks"] = tasks[1:]
    elif mutation == "reduced_model":
        report["models"]["ec"]["tasks"] = tasks[1:]
    elif mutation == "success":
        report["models"]["plain"]["status"] = "optimal"
    elif mutation == "biology":
        report["biological_validation"] = True
    _json(root, inventory.TASK_REPORT, report)
    _json(root, "data/current_claims.json", {"profiles": {"cosmic": {"evidence": [
        {"path": inventory.TASK_TABLE if mutation == "old_claim" else inventory.TASK_REPORT},
    ]}}})
    _json(root, "data/artifact_registry.json", {"runs": [{"artifacts": [] if mutation == "unregistered" else [{"path": inventory.TASK_REPORT}]}]})
    receipt = {"replacement_sha256": hashlib.sha256((root / inventory.TASK_REPORT).read_bytes()).hexdigest()}
    if mutation == "unbound_hash":
        receipt["replacement_sha256"] = "0" * 64
    elif mutation == "missing_report":
        (root / inventory.TASK_REPORT).unlink()
    if mutation is None:
        inventory.validate_task_replacement(root, receipt)
    else:
        with pytest.raises(inventory.InventoryError):
            inventory.validate_task_replacement(root, receipt)


@pytest.mark.parametrize("mutation", [None, "missing_report", "missing_history", "unknown_history", "dropped_negative", "restored_ranking", "unapproved_history"])
def test_n06_migration_retains_adopted_report_and_exact_original_history(mutation):
    registry = [{
        "id": "N06", "claim": "Upper bounds do not require product flux", "scope": "Retained scoped refusal",
        "artifact_paths": [inventory.TASK_REPORT], "historical_artifact_paths": [inventory.TASK_TABLE],
        "supporting_code": [], "supporting_tests": ["tests/test_run_regulation.py"],
        "claim_documents": ["data/current_claims.json"], "preservation": "Keep original Git identity and full-task refusal",
    }]
    paths = {inventory.TASK_REPORT, "tests/test_run_regulation.py", "data/current_claims.json"}
    replacements = {inventory.TASK_TABLE: inventory.TASK_REPORT}
    if mutation == "missing_report":
        paths.remove(inventory.TASK_REPORT)
    elif mutation == "missing_history":
        registry[0]["historical_artifact_paths"] = []
    elif mutation == "unknown_history":
        registry[0]["historical_artifact_paths"] = ["outputs/unreviewed.csv"]
    elif mutation == "dropped_negative":
        registry = []
    elif mutation == "restored_ranking":
        registry[0]["artifact_paths"].append(inventory.TASK_TABLE)
    elif mutation == "unapproved_history":
        replacements = {}
    arguments = (registry, paths, {inventory.TASK_TABLE}, {"records": {}}, ROOT)
    if mutation is not None:
        with pytest.raises(inventory.InventoryError, match="N06|historical artifact"):
            inventory.expand_negative_registry(*arguments, replacements=replacements)
        return
    expanded, membership = inventory.expand_negative_registry(*arguments, replacements=replacements)
    record = expanded[0]
    assert record["artifact_paths"] == record["retained_artifact_paths"] == [inventory.TASK_REPORT]
    assert set(record["retained_logical_artifact_paths"]) == {inventory.TASK_TABLE, inventory.TASK_REPORT}
    history = record["approved_result_replacements"][0]
    assert history["original_path"] == inventory.TASK_TABLE
    assert history["replacement_path"] == inventory.TASK_REPORT
    assert history["git_object_id"] == "764ef3aa9eb9c40f0b642243098d15002874c41d"
    assert history["git_commit"] == "84561f51502b5b4bf66c052dc88c853038917934"
    assert history["sha256"] == "6ed033a695c8d723a67034aa20ec35899c34ec9655772dceadbe9a4b27041ded"
    assert history["approval_request_id"] == inventory.CLEANUP_APPROVAL_ID
    assert membership == {inventory.TASK_REPORT: {"N06"}}


def test_reviewed_test_helper_has_explicit_use_and_owner_without_becoming_a_test_case(fixture_repo):
    root, seed, workflow, deletions, _ = fixture_repo
    helper = "tests/_frozen_helpers.py"
    _write(root, helper, "def run_original():\n    return 1\n")
    _write(root, "tests/test_validation.py", "from _frozen_helpers import run_original\ndef test_probe():\n    assert run_original() == 1\n")
    purpose = "Test-only original-source subprocess probe; not a production receipt generator"
    workflow["reviewed_session_additions"] = [{"path": helper, "work_item": "W13", "purpose": purpose}]
    result, refreshed, _ = inventory.build_inventory(root, seed, workflow, deletions)
    entry = next(item for item in result["entries"] if item["path"] == helper)
    assert entry["role"] == "validation_reproduction" and entry["purpose"] == purpose
    assert entry["work_items"] == ["W13"]
    assert entry["imported_by"] == {"test": ["tests/test_validation.py"]}
    assert all(helper not in node["tests"] for node in refreshed["nodes"])
    assert helper in refreshed["inventory_coverage"]["W13"]


@pytest.mark.parametrize("mutation", ["duplicate", "empty_purpose", "nontext_purpose", "unknown_owner"])
def test_reviewed_session_use_contracts_cannot_be_ambiguous(fixture_repo, mutation):
    root, seed, workflow, deletions, _ = fixture_repo
    addition = {"path": "tests/test_validation.py", "work_item": "W13", "purpose": "Original-source regression"}
    workflow["reviewed_session_additions"] = [addition]
    if mutation == "duplicate":
        workflow["reviewed_session_additions"].append(deepcopy(addition))
    elif mutation == "empty_purpose":
        addition["purpose"] = " "
    elif mutation == "nontext_purpose":
        addition["purpose"] = True
    else:
        addition["work_item"] = "W99"
    with pytest.raises(inventory.InventoryError):
        inventory.build_inventory(root, seed, workflow, deletions)


def test_real_refresh_covers_reviewed_untracked_and_changed_sources_and_every_receipt(monkeypatch):
    seed, workflow, deletions = [inventory.parse_json((ROOT / "data" / name).read_bytes()) for name in inventory.DOCUMENTS]
    approval_before = deepcopy(deletions["approval_request"])

    def no_mutation(*args, **kwargs):
        raise AssertionError("Inventory integration verification must be read-only")

    for method in ("unlink", "rename", "replace", "write_text", "write_bytes", "mkdir"):
        monkeypatch.setattr(Path, method, no_mutation)
    result, refreshed, manifest = inventory.build_inventory(ROOT, seed, workflow, deletions)
    entries = {item["path"]: item for item in result["entries"]}
    reviewed = {item["path"]: item for item in workflow["reviewed_session_additions"]}
    untracked = set(inventory.git(ROOT, "ls-files", "--others", "--exclude-standard", "-z").decode().split("\0")) - {""}
    assert untracked <= reviewed.keys(), sorted(untracked - reviewed.keys())
    assert untracked <= entries.keys()
    assert all(item["purpose"] and item["work_item"] in entries[name]["work_items"] for name, item in reviewed.items())
    declared_sources = set(reviewed)
    for node in workflow["nodes"]:
        declared_sources.update(node.get("existing_paths", []))
        declared_sources.update("src/" + name.replace(".", "/") + ".py" for name in node.get("existing_modules", []))
    changed = set(inventory.git(ROOT, "ls-files", "--modified", "-z").decode().split("\0")) - {""}
    changed_sources = {name for name in changed if name.endswith(".py") and name.startswith(("src/ystwin/", "scripts/"))}
    assert changed_sources <= declared_sources, sorted(changed_sources - declared_sources)
    assert not result["summary"]["unresolved_role_paths"]
    assert result["summary"]["completed_approved_cleanup_actions"] == len(inventory.CLEANUP_CONTRACTS)
    assert result["summary"]["pending_approved_cleanup_paths"] == []
    assert manifest["approval_request"] == approval_before
    assert manifest["execution_allowed"] is False
    assert len(manifest["execution_verification"]) == len(inventory.CLEANUP_CONTRACTS)
    assert {proposal["path"] for proposal in manifest["proposals"]} == set(inventory.CLEANUP_CONTRACTS)
    review = manifest["additional_cleanup_review"]
    assert review["approved"] is review["action_executed"] is True
    assert review["execution_allowed"] is False
    assert review["approval"]["actor"] == "user"
    assert review["approval"]["selection"] == "Approve both changes"
    assert review["exact_whole_file_deletion_proposals"] == []
    assert {proposal["id"] for proposal in review["partial_cleanup_proposals"]} == {"P01", "P02"}
    for proposal in review["partial_cleanup_proposals"]:
        assert proposal["path"] in entries and proposal["approval_required"] is True
        assert proposal["approved"] is proposal["action_executed"] is True
        assert "approval_request_id" not in proposal
        assert all(proposal[field] for field in ("scope", "role", "reason", "replacement_contract", "callers",
                                                 "unique_coverage", "historical_impact", "expected_verification"))
    assert result["summary"]["architecture_complete"] is False and refreshed["scope_complete"] is False
    assert [node["status"] for node in refreshed["nodes"]] == [node["status"] for node in workflow["nodes"]]
    n06 = next(record for record in result["retained_negative_registry"] if record["id"] == "N06")
    assert inventory.TASK_REPORT in n06["retained_artifact_paths"]
    assert inventory.TASK_TABLE not in entries
    assert entries[inventory.TASK_REPORT]["retained_negative_ids"] == ["N06"]
    assert entries[inventory.TASK_REPORT]["retirement_protected"] is True
    receipt = next(item for item in result["completed_cleanup_actions"] if item["path"] == inventory.TASK_TABLE)
    assert receipt["replacement_sha256"] == "5fff23d5c8e1d228783164546a6366291d2a77d1ef868b5853ded319b535eca4"
    assert n06["approved_result_replacements"][0]["git_object_id"] == receipt["original"]["git_object_id"]
    assert entries["src/ystwin/analysis/cosmic_data.py"]["imported_by"]["script"] == ["scripts/inspect_cosmic_data.py"]
    helper = entries["tests/_frozen_runtime_helpers.py"]
    assert helper["role"] == "validation_reproduction" and "W13" in helper["work_items"]
    assert {"tests/test_biology_learning_audit.py", "tests/test_native_reconciliation.py",
            "tests/test_portable_evidence_adversarial.py", "tests/test_portable_replay.py"} <= set(helper["imported_by"]["test"])


def test_clean_sdist_recreates_metadata_and_preserves_all_package_sources(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    generated = [ROOT / name for name, contract in inventory.CLEANUP_CONTRACTS.items()
                 if contract["action"] == "remove_generated_file"]
    before = {path: path.read_bytes() if path.exists() else None for path in generated}
    shutil.copy2(ROOT / "pyproject.toml", source / "pyproject.toml")
    shutil.copytree(ROOT / "src/ystwin", source / "src/ystwin", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    project = tomllib.loads((source / "pyproject.toml").read_text())["project"]
    assert not (source / "src/ystwin.egg-info").exists()
    result = subprocess.run(
        [sys.executable, "-B", "-c", "from setuptools.build_meta import build_sdist; build_sdist('dist')"],
        cwd=source, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    distributions = list((source / "dist").glob("*.tar.gz"))
    assert len(distributions) == 1
    prefix = f"{project['name']}-{project['version']}/"
    expected = {path.relative_to(source).as_posix(): path.read_bytes() for path in (source / "src/ystwin").rglob("*.py")}
    with tarfile.open(distributions[0], "r:gz") as archive:
        actual = {member.name.removeprefix(prefix): archive.extractfile(member).read()
                  for member in archive.getmembers() if member.isfile() and member.name.endswith(".py")}
        assert actual == expected
        metadata = BytesParser().parsebytes(archive.extractfile(prefix + "PKG-INFO").read())
    assert metadata["Name"] == project["name"] and metadata["Version"] == project["version"]
    requirements = metadata.get_all("Requires-Dist", [])
    assert set(project["dependencies"]) <= set(requirements)
    for extra, dependencies in project["optional-dependencies"].items():
        assert {dependency + f'; extra == "{extra}"' for dependency in dependencies} <= set(requirements)
    assert before == {path: path.read_bytes() if path.exists() else None for path in generated}
