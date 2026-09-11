from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from ystwin.analysis import frozen_runtime as runtime
from ystwin.analysis import portable_replay as replay
from ystwin.analysis.portable_evidence import load_portable_evidence


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def public_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("historical_public_only").resolve(strict=True)
    directory = Path(replay.DEFAULT_MANIFEST_PATH).parent
    shutil.copytree(ROOT / directory, root / directory)
    return root


@pytest.fixture(scope="module")
def snapshot(public_root):
    with runtime.materialize_frozen_runtime(public_root, source_repository=ROOT) as value:
        yield value


@pytest.fixture(scope="module")
def historical_report(public_root):
    return runtime.replay_historical_evidence(public_root, source_repository=ROOT)


@pytest.fixture(scope="module")
def audit_snapshot(public_root):
    with runtime.materialize_frozen_runtime(public_root, source_repository=ROOT, operation="audit") as value:
        yield value


@pytest.fixture(scope="module")
def historical_audit(public_root):
    return runtime.audit_historical_evidence(public_root, source_repository=ROOT,
                                             required_claims=["portable_scientific_content"],
                                             reuse_native_diagnostics=True)


def _copy_public(tmp_path, public_root):
    root = tmp_path.resolve(strict=True) / "public_copy"
    shutil.copytree(public_root, root)
    return root


def _forbidden_child(*args, **kwargs):
    raise AssertionError("unverified historical inputs must not reach the child interpreter")


def _runner():
    spec = importlib.util.spec_from_file_location("historical_replay_cli", ROOT / "scripts/replay_portable_evidence.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def _file_digests(root):
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file() and ".git" not in path.relative_to(root).parts}


def test_registry_binds_the_original_public_checkpoint_and_runtime(public_root):
    registry = runtime.load_frozen_runtime(ROOT)
    bundle = load_portable_evidence(public_root, replay.DEFAULT_MANIFEST_PATH, replay.DEFAULT_MANIFEST_SHA256)
    checkpoint = bundle.fetch(replay.CHECKPOINT_PATH, original_sha256=replay.CHECKPOINT_SHA256).payload
    assert registry["code_ref"] == "84561f51502b5b4bf66c052dc88c853038917934"
    assert registry["environment"] == checkpoint["model"]["config"]["environment"]
    assert registry["environment"]["python"] == "3.14.2"
    assert registry["checkpoint"]["model_sha256"] == replay.MODEL_SHA256
    assert registry["manifest"]["sha256"] == replay.DEFAULT_MANIFEST_SHA256
    assert len(checkpoint["training_manifest"]["code"]) == 20


def test_snapshot_provisions_missing_sources_without_reconstructing_private_originals(public_root, snapshot):
    before = _file_digests(public_root)
    assert not (public_root / "data/hog2013/model_wt.xml").exists()
    assert not (snapshot.root / replay.CHECKPOINT_PATH).exists()
    assert not (snapshot.root / ".git").exists()
    assert snapshot.root == snapshot.root.resolve(strict=True)
    assert len(snapshot.supplied_data) == 14
    assert all(row["status"] == "absent_provisioned_from_verified_git" for row in snapshot.supplied_data)
    for path, identity in snapshot.identities.items():
        content = (snapshot.root / path).read_bytes()
        assert hashlib.sha256(content).hexdigest() == identity["sha256"]
        assert len(content) == identity["size_bytes"]
    assert snapshot.identities["pyproject.toml"]["sha256"] == "0bf5440e9a8107c9d5e0e7ad1a596466be43f1f8af958129750fece811c2fc90"
    assert _file_digests(public_root) == before


def test_historical_replay_records_real_isolated_execution(historical_report):
    execution = historical_report["historical_execution"]
    assert historical_report["verified"] and execution["child_exit_code"] == 0
    assert execution["isolated_interpreter"] and execution["bytecode_disabled"]
    assert execution["pythonpath"] == ["src"] and len(execution["loaded_modules"]) == 19
    assert execution["git_bytes_verified_before_execution"]
    assert execution["snapshot_and_supplied_inputs_unchanged"]
    assert not execution["current_worktree_code_used"]
    assert execution["runtime_versions"] == historical_report["integrity"]["runtime"]["recorded_versions"]
    assert not historical_report["new_fitting"] and not historical_report["biological_validation"]
    assert historical_report["integrity"]["original_checkpoint_bytes"]["status"] == "not_checked"
    assert str(ROOT) not in json.dumps(historical_report)


def test_changed_current_code_cannot_change_historical_numerics_or_mutate_tracked_files(
    tmp_path, snapshot, historical_report,
):
    root = tmp_path.resolve(strict=True) / "changed_current_worktree"
    shutil.copytree(snapshot.root, root)
    for relative in ("scripts/replay_portable_evidence.py", "src/ystwin/analysis/portable_replay.py",
                     "src/ystwin/analysis/frozen_runtime.py"):
        shutil.copy2(ROOT / relative, root / relative)
    (root / "pyproject.toml").write_text("[project]\nname = 'changed-current-project'\nversion = '9.9.9'\n")
    (root / "src/ystwin/mech/hog.py").write_text("raise AssertionError('the current HOG kernel must never execute')\n")
    (root / "scripts/export_granados_training.py").write_text("raise AssertionError('lineage-only exporter must never execute')\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        replay.load_portable_native_checkpoint(root)
    subprocess.run(["git", "init", "--quiet", str(root)], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True, capture_output=True)
    before = _file_digests(root)
    environment = dict(os.environ, PYTHONPATH=str(root / "src"), PYTHONDONTWRITEBYTECODE="1",
                       YSTWIN_HOG_REFERENCE="/not/a/historical/source", YSTWIN_OUTPUTS=str(root / "forbidden_outputs"))
    process = subprocess.run(
        [sys.executable, str(root / "scripts/replay_portable_evidence.py"), "--root", str(root),
         "--frozen-code-ref", runtime.FROZEN_CODE_REF, "--source-repository", str(ROOT)],
        cwd=root, env=environment, capture_output=True, text=True, timeout=120,
    )
    assert process.returncode == 0, process.stderr
    report = json.loads(process.stdout)
    for field in ("native_hog", "native_obligations", "native_exchange", "transfer_scores", "population_metrics"):
        assert report[field] == historical_report[field]
    assert report["historical_execution"]["child_exit_code"] == 0
    assert report["historical_execution"]["loaded_modules"] == historical_report["historical_execution"]["loaded_modules"]
    assert _file_digests(root) == before
    status = subprocess.run(["git", "-C", str(root), "diff", "--exit-code"], capture_output=True, text=True)
    assert status.returncode == 0, status.stdout + status.stderr
    untracked = subprocess.run(["git", "-C", str(root), "ls-files", "--others", "--exclude-standard"],
                               check=True, capture_output=True, text=True)
    assert not untracked.stdout


@pytest.mark.parametrize("operation", ["replay_historical_evidence", "audit_historical_evidence"])
@pytest.mark.parametrize("relative", ["data/hog2013/model_wt.xml", "data/holdout_transfer/labels.json"])
def test_altered_supplied_source_data_is_not_hidden_by_git_provisioning(tmp_path, public_root, snapshot, monkeypatch, relative, operation):
    root = _copy_public(tmp_path, public_root)
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((snapshot.root / relative).read_bytes() + b"\n")
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="supplied source data SHA-256 mismatch"):
        getattr(runtime, operation)(root, source_repository=ROOT)


@pytest.mark.parametrize("relative", ["pyproject.toml", "src/ystwin/analysis/portable_replay.py",
                                      "data/hog2013/model_wt.xml", replay.DEFAULT_MANIFEST_PATH])
def test_git_bytes_are_checked_against_frozen_hashes_before_execution(public_root, monkeypatch, relative):
    original = runtime._git_blobs

    def corrupt(*args, **kwargs):
        blobs = original(*args, **kwargs)
        blobs[relative] += b"\n"
        return blobs

    monkeypatch.setattr(runtime, "_git_blobs", corrupt)
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="Git blob SHA-256 mismatch"):
        runtime.replay_historical_evidence(public_root, source_repository=ROOT)


@pytest.mark.parametrize("field", ["code_ref", "manifest", "checkpoint", "environment", "replay_support"])
def test_changed_runtime_registry_cannot_restamp_historical_identities(tmp_path, field):
    registry = runtime.load_frozen_runtime(ROOT)
    if field == "code_ref":
        registry[field] = "0" * 40
    elif field == "environment":
        registry[field]["python"] = "0.0.0"
    elif field == "replay_support":
        registry[field][0]["sha256"] = "0" * 64
    else:
        registry[field]["sha256"] = "0" * 64
    path = tmp_path / runtime.RUNTIME_REGISTRY_PATH
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(registry))
    with pytest.raises(ValueError, match="runtime registry SHA-256 mismatch"):
        runtime.load_frozen_runtime(tmp_path)


@pytest.mark.parametrize("repin", [False, True])
def test_changed_public_manifest_is_not_replaced_by_its_git_copy(tmp_path, public_root, monkeypatch, repin):
    root = _copy_public(tmp_path, public_root)
    path = root / replay.DEFAULT_MANIFEST_PATH
    data = path.read_bytes() + b"\n"
    path.write_bytes(data)
    expected = hashlib.sha256(data).hexdigest() if repin else None
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="manifest SHA-256|registered frozen manifest identity"):
        runtime.replay_historical_evidence(root, source_repository=ROOT, expected_sha256=expected)


def test_changed_public_checkpoint_is_not_replaced_by_its_git_copy(tmp_path, public_root, monkeypatch):
    root = _copy_public(tmp_path, public_root)
    path = root / Path(replay.DEFAULT_MANIFEST_PATH).parent / "artifacts" / replay.CHECKPOINT_PATH
    path.write_bytes(path.read_bytes() + b"\n")
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="public artifact SHA-256"):
        runtime.replay_historical_evidence(root, source_repository=ROOT)


@pytest.mark.parametrize("reference", ["HEAD", "84561f5", "0" * 40, "--help", "../84561f5"])
def test_only_the_explicit_registered_commit_can_execute(public_root, monkeypatch, reference):
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="registered full frozen Git commit"):
        runtime.replay_historical_evidence(public_root, source_repository=ROOT, frozen_code_ref=reference)


def test_missing_commit_blocks_instead_of_using_current_source(tmp_path, public_root, monkeypatch):
    repository = tmp_path.resolve(strict=True) / "empty_repository"
    repository.mkdir()
    subprocess.run(["git", "init", "--quiet", str(repository)], check=True, capture_output=True)
    (repository / "data").mkdir()
    shutil.copy2(ROOT / runtime.RUNTIME_REGISTRY_PATH, repository / runtime.RUNTIME_REGISTRY_PATH)
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="frozen Git snapshot unavailable"):
        runtime.replay_historical_evidence(public_root, source_repository=repository)


@pytest.mark.parametrize("path", ["../outside.py", "/outside.py", "src/./kernel.py", "src\\kernel.py", "https://example.invalid/kernel.py"])
def test_snapshot_dependency_paths_cannot_alias_or_escape(path):
    with pytest.raises(ValueError, match="repository-relative"):
        runtime._identity({"path": path, "sha256": "0" * 64, "size_bytes": 1})


def test_snapshot_resolved_paths_cannot_alias():
    with pytest.raises(ValueError, match="aliases"):
        runtime._identity({"path": "src/kernel.py", "resolved_path": "src/other.py", "sha256": "0" * 64, "size_bytes": 1})


def test_supplied_data_symlinks_are_rejected_even_if_target_bytes_match(tmp_path, public_root, snapshot, monkeypatch):
    root = _copy_public(tmp_path, public_root)
    (root / "data/hog2013").symlink_to(snapshot.root / "data/hog2013", target_is_directory=True)
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="nonsymlink canonical path"):
        runtime.replay_historical_evidence(root, source_repository=ROOT)


def test_git_tree_symlinks_are_not_materialized(public_root, monkeypatch):
    original = runtime._git

    def alias(repository, *arguments, **kwargs):
        result = original(repository, *arguments, **kwargs)
        if arguments[0] == "ls-tree":
            result = result.replace(b"100644 blob ", b"120000 blob ", 1)
        return result

    monkeypatch.setattr(runtime, "_git", alias)
    with pytest.raises(ValueError, match="regular files, never symlinks"):
        runtime._git_blobs(ROOT, runtime.FROZEN_CODE_REF, ["pyproject.toml"])


def test_git_environment_cannot_redirect_the_registered_repository(tmp_path, monkeypatch, snapshot):
    monkeypatch.setenv("GIT_DIR", str(tmp_path / "wrong_git_directory"))
    monkeypatch.setenv("GIT_WORK_TREE", str(tmp_path))
    content = runtime._git_blobs(ROOT, runtime.FROZEN_CODE_REF, ["pyproject.toml"])["pyproject.toml"]
    assert hashlib.sha256(content).hexdigest() == snapshot.identities["pyproject.toml"]["sha256"]


@pytest.mark.parametrize("field", ["python", "packages"])
def test_child_measures_runtime_versions_instead_of_trusting_requested_metadata(snapshot, field):
    registry = deepcopy(snapshot.registry)
    if field == "python":
        registry["environment"]["python"] = "0.0.0"
    else:
        registry["environment"]["packages"]["numpy"] = "0.0.0"
    altered = replace(snapshot, registry=registry)
    process = runtime._run_child(altered, sys.executable, None, 60)
    assert process.returncode != 0 and "frozen runtime versions differ" in process.stderr
    assert not process.stdout


def test_real_subprocess_failure_is_blocking(public_root):
    executable = shutil.which("false")
    assert executable is not None
    with pytest.raises(runtime.HistoricalReplayError) as error:
        runtime.replay_historical_evidence(public_root, source_repository=ROOT, python_executable=executable)
    assert error.value.returncode == 1
    assert error.value.report["child_exit_code"] == 1 and not error.value.report["verified"]


def test_nonzero_child_exit_cannot_be_overruled_by_success_shaped_stdout(public_root, monkeypatch):
    monkeypatch.setattr(runtime, "_run_child", lambda *args: subprocess.CompletedProcess(
        args=[], returncode=23, stdout='{"verified": true}', stderr="child stopped"))
    with pytest.raises(runtime.HistoricalReplayError) as error:
        runtime.replay_historical_evidence(public_root, source_repository=ROOT)
    assert error.value.returncode == 23 and error.value.report["stderr"] == "child stopped"
    assert not error.value.report["verified"]


@pytest.mark.parametrize("stdout", ["", "{}", '{"verified": true}'])
def test_zero_exit_requires_a_complete_verified_child_report(snapshot, stdout):
    process = subprocess.CompletedProcess(args=[], returncode=0, stdout=stdout, stderr="")
    with pytest.raises(runtime.HistoricalReplayError, match="invalid historical child report") as error:
        runtime._child_report(process, snapshot)
    assert error.value.returncode == 0 and not error.value.report["verified"]


def test_child_report_cannot_substitute_current_module_hashes(snapshot, historical_report):
    report = deepcopy(historical_report)
    report["historical_execution"]["loaded_modules"][0]["sha256"] = "0" * 64
    process = subprocess.CompletedProcess(args=[], returncode=0, stdout=json.dumps(report), stderr="")
    with pytest.raises(runtime.HistoricalReplayError, match="child loaded source identities"):
        runtime._child_report(process, snapshot)


def test_timeout_does_not_report_success(public_root, monkeypatch):
    def timeout(*args):
        raise subprocess.TimeoutExpired("historical child", 0.01, output=b"partial", stderr=b"interrupted")

    monkeypatch.setattr(runtime, "_run_child", timeout)
    with pytest.raises(runtime.HistoricalReplayError, match="timed out") as error:
        runtime.replay_historical_evidence(public_root, source_repository=ROOT)
    assert error.value.report["child_exit_code"] is None
    assert error.value.report["stdout"] == "partial" and error.value.report["stderr"] == "interrupted"


def test_cli_propagates_child_exit_and_never_writes_a_success_record(tmp_path, public_root, monkeypatch, capsys):
    monkeypatch.setattr(runtime, "_run_child", lambda *args: subprocess.CompletedProcess(
        args=[], returncode=29, stdout="", stderr="blocked child"))
    output = tmp_path.resolve(strict=True) / "must_not_exist.json"
    result = _runner().main(["--root", str(public_root), "--historical", "--source-repository", str(ROOT),
                             "--output", str(output)])
    assert result == 29 and not output.exists()
    captured = capsys.readouterr()
    assert not captured.out
    report = json.loads(captured.err)
    assert report["child_exit_code"] == 29 and not report["verified"]


def test_cli_historical_options_are_never_implicit(public_root):
    with pytest.raises(SystemExit) as error:
        _runner().main(["--root", str(public_root), "--source-repository", str(ROOT)])
    assert error.value.code == 2


def test_cli_does_not_resolve_away_an_aliased_output_parent(tmp_path, public_root, monkeypatch):
    target = tmp_path.resolve(strict=True) / "actual_parent"
    target.mkdir()
    alias = tmp_path.resolve(strict=True) / "aliased_parent"
    alias.symlink_to(target, target_is_directory=True)
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="nonsymlink parent"):
        _runner().main(["--root", str(public_root), "--historical", "--source-repository", str(ROOT),
                        "--output", str(alias / "replay.json")])
    assert not (target / "replay.json").exists()


def test_historical_audit_uses_original_collector_and_preserves_diagnostic_history(historical_audit, audit_snapshot):
    report = historical_audit
    assert report["kind"] == runtime.HISTORICAL_AUDIT_KIND
    assert report["verified"] and report["required_claims_satisfied"]
    assert report["historical_execution"]["child_exit_code"] == 0
    assert report["historical_execution"]["operation"] == "audit"
    assert len(report["historical_execution"]["loaded_modules"]) == 20
    assert report["ledger"]["evidence_predicates"]["checkpoint_integrity"]["value"] is None
    assert report["integrity"]["original_checkpoint_bytes"]["status"] == "not_checked"
    assert [row["claim"] for row in report["authorized_claims"]] == ["portable_scientific_content"]
    assert not report["new_fitting"] and not report["biological_validation"]
    assert not report["historical_diagnostic_code_reexecuted"]
    origin = report["diagnostics_origin"]
    assert origin["mode"] == "reused_unchanged_retrospective_evidence"
    assert not origin["numerical_diagnostics_rerun"] and not origin["biological_validation_credit"]
    assert origin["source"]["sha256"] == "b88d338879fc617b1ebfd82fe88ee9fe6fd0fd9cb927f9db973939b134f19ec8"
    plan = json.loads((audit_snapshot.root / origin["original_plan"]["path"]).read_bytes())
    assert plan["audit_code"][0]["sha256"] == "bac49aa1a6a8421886e78e006a709a81845a98a7d4fbe75b639d6cd4d97a43d8"
    loaded = {row["path"]: row["sha256"] for row in report["historical_execution"]["loaded_modules"]}
    assert loaded["src/ystwin/analysis/biology_learning_audit.py"] == "aefb6953fa2cd9ac50d541b6ed57e5ff6795062ae3da309e1d38fc3b0895e2e5"
    assert str(ROOT) not in json.dumps(report) and str(audit_snapshot.root) not in json.dumps(report)


def test_historical_claim_gate_rechecks_a_serialized_ledger_without_trusting_it(public_root, historical_audit):
    ledger = json.loads(json.dumps(historical_audit["ledger"]))
    ledger["untrusted_extra_metadata"] = "x" * 200000
    result = runtime.audit_historical_evidence(public_root, source_repository=ROOT, ledger=ledger,
                                               required_claims=["portable_scientific_content"])
    assert result["required_claims_satisfied"] and result["historical_execution"]["child_exit_code"] == 0
    assert "untrusted_extra_metadata" not in result["ledger"]
    assert result["diagnostics_origin"]["mode"] == "not_run_historical_source_audit"


@pytest.mark.parametrize("mutation", ["reference", "parameter", "original_integrity"])
def test_historical_claim_gate_rejects_changed_ledger_fields(public_root, historical_audit, mutation):
    ledger = deepcopy(historical_audit["ledger"])
    if mutation == "reference":
        ledger["evidence_predicates"]["native_cells_bound"]["evidence"][0]["path"] = "../unapproved.json"
    elif mutation == "parameter":
        hog = next(row for row in ledger["components"] if row["id"] == "native_hog")
        hog["absolute_parameter_overrides"]["kv16f_1"] *= 2
    else:
        ledger["evidence_predicates"]["checkpoint_integrity"]["value"] = True
    with pytest.raises(runtime.HistoricalReplayError) as error:
        runtime.audit_historical_evidence(public_root, source_repository=ROOT, ledger=ledger,
                                          required_claims=["portable_scientific_content"])
    assert error.value.returncode == 1
    assert "authoritative collector" in error.value.report["stderr"]


@pytest.mark.parametrize("relative", ["data/gem/MANIFEST.md", "outputs/biology_learning_audit_01/native_diagnostics.json"])
def test_historical_audit_does_not_hide_tampered_provenance_or_diagnostics(
    tmp_path, public_root, audit_snapshot, monkeypatch, relative,
):
    root = _copy_public(tmp_path, public_root)
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((audit_snapshot.root / relative).read_bytes() + b"\n")
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="supplied source data SHA-256 mismatch"):
        runtime.audit_historical_evidence(root, source_repository=ROOT, reuse_native_diagnostics=True)


def test_historical_audit_checks_git_collector_bytes_before_import(public_root, monkeypatch):
    original = runtime._git_blobs

    def corrupt(*args, **kwargs):
        blobs = original(*args, **kwargs)
        blobs["src/ystwin/analysis/biology_learning_audit.py"] += b"\n"
        return blobs

    monkeypatch.setattr(runtime, "_git_blobs", corrupt)
    monkeypatch.setattr(runtime, "_run_child", _forbidden_child)
    with pytest.raises(ValueError, match="Git blob SHA-256 mismatch"):
        runtime.audit_historical_evidence(public_root, source_repository=ROOT)


def test_lineage_exporter_pin_is_retained_without_requiring_current_exporter_bytes(public_root, snapshot, audit_snapshot):
    bundle = load_portable_evidence(public_root, replay.DEFAULT_MANIFEST_PATH, replay.DEFAULT_MANIFEST_SHA256)
    path = "scripts/export_granados_training.py"
    reference = next(row for row in bundle.manifest["external_references"] if row["origin_path"] == path)
    assert reference["original_sha256"] == "6447608e622c6bdd47e16f4f930f91f9297077aaccfbc5a0745b359b1bbd7159"
    assert reference["integrity_scope"] == "lineage_only_not_exported_or_independently_verified"
    original = runtime._git_blobs(ROOT, runtime.FROZEN_CODE_REF, [path])[path]
    assert hashlib.sha256(original).hexdigest() == reference["original_sha256"]
    assert path not in snapshot.identities and path not in audit_snapshot.identities


def _audit_cli(root, output, *arguments):
    environment = {key: value for key, value in os.environ.items() if not key.startswith(("PYTHON", "YSTWIN_"))}
    environment.update(PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts/audit_biology_learning.py"), "--root", str(root),
         "--output-dir", str(output), *arguments], cwd=root, env=environment,
        capture_output=True, text=True, timeout=120,
    )


def test_historical_audit_cli_reuses_original_diagnostic_bytes_and_is_write_once(tmp_path, public_root, audit_snapshot):
    output = tmp_path.resolve(strict=True) / "new_historical_audit"
    arguments = ("--historical", "--source-repository", str(ROOT), "--reuse-native-diagnostics",
                 "--require-claim", "portable_scientific_content")
    process = _audit_cli(public_root, output, *arguments)
    assert process.returncode == 0, process.stderr
    report = json.loads(process.stdout)
    assert report == json.loads((output / "audit_report.json").read_bytes())
    assert report["ledger"] == json.loads((output / "ledger.json").read_bytes())
    original = audit_snapshot.root / report["diagnostics_origin"]["source"]["path"]
    assert (output / "native_diagnostics.json").read_bytes() == original.read_bytes()
    before = _file_digests(output)
    again = _audit_cli(public_root, output, *arguments)
    assert again.returncode != 0 and "write-once" in again.stderr
    assert _file_digests(output) == before


@pytest.mark.parametrize("claim", ["integrity_established", "independent_mechanistic_evidence"])
def test_historical_audit_cli_blocks_unsupported_required_claims(tmp_path, public_root, claim):
    output = tmp_path.resolve(strict=True) / "rejected_claim_audit"
    process = _audit_cli(public_root, output, "--historical", "--source-repository", str(ROOT),
                         "--require-claim", claim)
    assert process.returncode == 2, process.stderr
    report = json.loads(process.stdout)
    assert report["verified"] and not report["required_claims_satisfied"]
    assert report["historical_execution"]["child_exit_code"] == 0
    assert not report["authorized_claims"]
    assert report["rejected_required_claims"][0]["claim"] == claim


def test_default_audit_cli_does_not_fall_back_to_public_historical_evidence(tmp_path, audit_snapshot):
    root = tmp_path.resolve(strict=True) / "strict_public_only"
    shutil.copytree(audit_snapshot.root, root)
    output = root / "outputs/biology_learning_audit_strict_refusal"
    process = _audit_cli(root, output, "--reuse-native-diagnostics")
    assert process.returncode == 1 and "nonsymlink file" in process.stderr
    assert not (output / "ledger.json").exists()
    assert not (output / "audit_report.json").exists()


def test_historical_audit_cli_never_reaches_changed_current_audit_or_native_code(
    tmp_path, audit_snapshot, historical_audit,
):
    root = tmp_path.resolve(strict=True) / "changed_current_audit"
    shutil.copytree(audit_snapshot.root, root)
    for relative in ("scripts/audit_biology_learning.py", "src/ystwin/analysis/frozen_runtime.py",
                     "src/ystwin/analysis/portable_replay.py"):
        shutil.copy2(ROOT / relative, root / relative)
    for relative in ("src/ystwin/analysis/biology_learning_audit.py", "src/ystwin/mech/hog.py",
                     "scripts/export_granados_training.py"):
        (root / relative).write_text("raise AssertionError('current code must not execute in historical auditing')\n")
    (root / "pyproject.toml").write_text("[project]\nname = 'changed-current-audit'\n")
    before = _file_digests(root)
    output = tmp_path.resolve(strict=True) / "historical_audit_from_changed_code"
    environment = dict(os.environ, PYTHONPATH=str(root / "src"), PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.run(
        [sys.executable, str(root / "scripts/audit_biology_learning.py"), "--root", str(root),
         "--historical", "--source-repository", str(ROOT), "--output-dir", str(output),
         "--require-claim", "portable_scientific_content", "--reuse-native-diagnostics"],
        cwd=root, env=environment, capture_output=True, text=True, timeout=120,
    )
    assert process.returncode == 0, process.stderr
    report = json.loads(process.stdout)
    assert report["ledger"] == historical_audit["ledger"]
    assert report["findings"] == historical_audit["findings"]
    assert _file_digests(root) == before


def test_historical_audit_propagates_subprocess_failure_without_creating_output(tmp_path, public_root, monkeypatch, capsys):
    monkeypatch.setattr(runtime, "_run_child", lambda *args, **kwargs: subprocess.CompletedProcess(
        args=[], returncode=31, stdout="", stderr="audit child stopped"))
    spec = importlib.util.spec_from_file_location("frozen_audit_cli", ROOT / "scripts/audit_biology_learning.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    output = tmp_path.resolve(strict=True) / "must_not_be_an_audit"
    result = runner.main(["--root", str(public_root), "--historical", "--source-repository", str(ROOT),
                          "--output-dir", str(output)])
    assert result == 31 and not output.exists()
    captured = capsys.readouterr()
    assert not captured.out and json.loads(captured.err)["child_exit_code"] == 31
