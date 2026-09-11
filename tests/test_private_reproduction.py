from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import stat

import pytest

from ystwin import artifacts


ROOT = Path(__file__).resolve().parents[1]
POLICY = "data/private_reproduction_policy.json"
PRIVATE = "data/gen5/declared.xpt"
PUBLIC = "data/public.csv"
OUTPUT = "outputs/result.csv"
CANARY = "private-author-metadata-canary"


def _json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")


def _save(root, document):
    _json(root / artifacts.REGISTRY_PATH, document)


def _stat_bytes(path):
    info = path.stat()
    return path.read_bytes(), info.st_mtime_ns, stat.S_IMODE(info.st_mode)


@pytest.fixture
def private_specimen(tmp_path):
    root = tmp_path / "repository"
    corpus = tmp_path / "approved-corpus"
    for directory in (root / "scripts", root / "src/ystwin", root / "data", root / "outputs", corpus):
        directory.mkdir(parents=True)
    (root / "src/ystwin/__init__.py").write_text("")
    payload = "id,value,unit\na,1,mg\nb,,mg\n"
    (corpus / "declared.xpt").write_text(CANARY + "\n" + payload)
    (root / PUBLIC).write_text("required public calibration\n")
    (root / OUTPUT).write_text(payload)
    (root / "scripts/run.py").write_text(
        "from pathlib import Path\n"
        f"source = Path({PRIVATE!r})\n"
        "metadata, table = source.read_text().split('\\n', 1)\n"
        "print(metadata, source.resolve())\n"
        f"assert Path({PUBLIC!r}).read_text() == 'required public calibration\\n'\n"
        f"Path({OUTPUT!r}).write_text(table)\n")
    run = {
        "id": "private-fixture", "lifecycle": "current", "status": "investigation_pending",
        "reason": "Explicit synthetic private fixture, never a public replay",
        "producer": {"path": "scripts/run.py", "callable": "main", "libraries": []},
        "dependencies": [{"path": PRIVATE, "role": "input", "sha256": None},
                         {"path": PUBLIC, "role": "input", "sha256": None}],
        "parameters": {"selection": "all rows"},
        "model_identity": {"kind": "none", "reason": "identity transformation"},
        "runtime": {"kind": "python_standard_library"},
        "reproduction": {"argv": ["{python}", "-B", "scripts/run.py"],
                         "environment": {"YSTWIN_GEN5_XPT": "{workspace}/data/gen5"}},
        "artifacts": [{"path": OUTPUT, "comparison": {
            "format": "csv", "row_keys": ["id"], "units": {"value": "mg"}, "unit_columns": ["unit"]}}],
    }
    document = {"schema_version": 1, "runs": [run], "frozen_bundles": []}
    _save(root, document)
    policy = {
        "schema_version": 1, "permission": "explicit_env_file_read_only",
        "ci_private_bytes": "skip_unless_explicitly_reverified",
        "sources": {"data/gen5": {"environment": "YSTWIN_GEN5_XPT", "suffix": ".xpt", "resolver": "gen5_xpt_dir"}},
        "runs": {run["id"]: {PRIVATE: {"kind": "file"}}},
    }
    _json(root / POLICY, policy)
    env_file = root / ".ystwin.env"
    env_file.write_text(f"YSTWIN_GEN5_XPT={corpus}\nUNRELATED_SECRET=do-not-evaluate-this\n")
    return root, document, corpus, env_file


def _stage(specimen, tmp_path):
    root, document, _, env_file = specimen
    directory = tmp_path / "private-stage"
    receipt = artifacts.reproduce_run(root, document["runs"][0]["id"], directory,
                                      private_inputs_from_env=env_file)
    return directory, receipt


def _adopt(specimen, tmp_path):
    directory, receipt = _stage(specimen, tmp_path)
    assert receipt["success"], receipt
    result = artifacts.adopt_reproduction(specimen[0], directory / "reproduction.json",
                                          private_inputs_from_env=specimen[3])
    return directory, receipt, result


def _artifact_row(rows):
    return next(row for row in rows if row["check"] == f"artifact provenance: {OUTPUT}")


def test_default_denial_even_with_ambient_env_and_local_private_copy(private_specimen, tmp_path, monkeypatch):
    root, document, corpus, _ = private_specimen
    monkeypatch.setenv("YSTWIN_GEN5_XPT", str(corpus))
    (root / PRIVATE).parent.mkdir()
    (root / PRIVATE).write_bytes((corpus / "declared.xpt").read_bytes())
    before = _stat_bytes(corpus / "declared.xpt")
    with pytest.raises(artifacts.ArtifactError, match="explicit.*private-inputs-from-env"):
        artifacts.reproduce_run(root, document["runs"][0]["id"], tmp_path / "denied")
    assert not (tmp_path / "denied").exists()
    assert _stat_bytes(corpus / "declared.xpt") == before


def test_explicit_private_staging_binds_logical_sha_without_path_or_metadata_leaks(private_specimen, tmp_path):
    root, _, corpus, env_file = private_specimen
    before = {path: _stat_bytes(path) for path in (corpus / "declared.xpt", root / OUTPUT, env_file)}
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["success"] and receipt["execution_complete"] and receipt["exit_code"] == 0
    private = next(entry for entry in receipt["identity"]["dependencies"] if entry["path"] == PRIVATE)
    assert private == {"path": PRIVATE, "role": "input", "sha256": artifacts.sha256(corpus / "declared.xpt"), "private": True}
    assert receipt["identity"]["private_inputs"]["policy"] == {"path": POLICY, "sha256": artifacts.sha256(root / POLICY)}
    assert receipt["execution"]["read_paths"] == [PRIVATE, PUBLIC, "scripts/run.py"]
    assert before == {path: _stat_bytes(path) for path in before}
    for path in [directory, *directory.rglob("*")]:
        assert stat.S_IMODE(path.stat().st_mode) == (0o700 if path.is_dir() else 0o600), path.relative_to(directory)
    for name in ("execution_plan.json", "child_execution.json", "reproduction.json", "stdout.txt", "stderr.txt"):
        text = (directory / name).read_text()
        assert str(corpus) not in text and str(root) not in text and CANARY not in text
    assert "suppressed" in (directory / "stdout.txt").read_text().lower()
    assert not (root / PRIVATE).exists()


def test_missing_private_is_not_partial_execution_and_missing_public_never_becomes_private(private_specimen, tmp_path):
    root, document, corpus, env_file = private_specimen
    private = corpus / "declared.xpt"
    old = private.read_bytes()
    private.unlink()
    with pytest.raises(artifacts.ArtifactError, match="private input.*unavailable") as failure:
        _stage(private_specimen, tmp_path)
    assert str(corpus) not in str(failure.value)
    private.write_bytes(old)
    (root / PUBLIC).unlink()
    with pytest.raises(artifacts.ArtifactError, match="required input unavailable: data/public.csv"):
        artifacts.reproduce_run(root, document["runs"][0]["id"], tmp_path / "public-missing",
                                private_inputs_from_env=env_file)
    assert not (tmp_path / "private-stage").exists()


@pytest.mark.parametrize("bad_source", [
    {"environment": "UNRELATED_SECRET"},
    {"environment": "YSTWIN_QPCR_RAW"},
    {"environment": "YSTWIN_GEN5_XPT", "relative_path": "../outside.xpt"},
    {"environment": "YSTWIN_GEN5_XPT", "relative_path": "/elsewhere/file.xpt"},
    {"environment": "YSTWIN_GEN5_XPT", "relative_path": "other.xpt"},
    {"environment": "YSTWIN_GEN5_XPT", "resolver": "other_resolver"},
])
def test_registry_source_cannot_redirect_policy_permission(private_specimen, tmp_path, bad_source):
    root, document, _, _ = private_specimen
    document["runs"][0]["dependencies"][0]["source"] = bad_source
    _save(root, document)
    with pytest.raises(artifacts.ArtifactError, match="private source.*policy"):
        _stage(private_specimen, tmp_path)


@pytest.mark.parametrize("mode", ["file", "parent", "root", "env"])
def test_private_source_symlinks_are_rejected_before_ingestion(private_specimen, tmp_path, mode):
    root, _, corpus, env_file = private_specimen
    if mode == "file":
        original = corpus / "declared.xpt"
        saved = tmp_path / "other.xpt"
        original.rename(saved)
        original.symlink_to(saved)
    elif mode in {"parent", "root"}:
        link = tmp_path / "alias"
        link.symlink_to(corpus if mode == "root" else tmp_path, target_is_directory=True)
        location = link if mode == "root" else link / corpus.name
        env_file.write_text(f"YSTWIN_GEN5_XPT={location}\n")
    else:
        saved = tmp_path / "settings.env"
        env_file.rename(saved)
        env_file.symlink_to(saved)
    with pytest.raises(artifacts.ArtifactError, match="symlink") as failure:
        _stage(private_specimen, tmp_path)
    assert str(tmp_path) not in str(failure.value)
    assert not (tmp_path / "private-stage").exists()


@pytest.mark.parametrize("setting", [
    "YSTWIN_GEN5_XPT=\n", "YSTWIN_GEN5_XPT=relative/path\n",
    "YSTWIN_GEN5_XPT=$(touch forbidden)\n", "YSTWIN_GEN5_XPT=\"unterminated\n",
    "YSTWIN_GEN5_XPT=/missing\nYSTWIN_GEN5_XPT=/also-missing\n",
    "OTHER_SOURCE=/missing\n",
])
def test_env_parser_is_data_not_shell_and_has_no_fallback(private_specimen, tmp_path, setting):
    private_specimen[3].write_text(setting)
    with pytest.raises(artifacts.ArtifactError):
        _stage(private_specimen, tmp_path)
    assert not (tmp_path / "private-stage").exists()


def test_env_root_path_escape_is_rejected(private_specimen, tmp_path):
    _, _, corpus, env_file = private_specimen
    env_file.write_text(f"YSTWIN_GEN5_XPT={corpus}/../{corpus.name}\n")
    with pytest.raises(artifacts.ArtifactError, match="canonical"):
        _stage(private_specimen, tmp_path)


def test_quoted_env_root_is_supported_without_evaluating_other_names(private_specimen, tmp_path):
    _, _, corpus, env_file = private_specimen
    env_file.write_text(f"export YSTWIN_GEN5_XPT='{corpus}'\nOTHER=$(never-run)\n")
    assert _stage(private_specimen, tmp_path)[1]["success"]


def test_unknown_private_file_and_unrelated_extension_are_blocking(private_specimen, tmp_path):
    root, document, _, _ = private_specimen
    document["runs"][0]["dependencies"].append({"path": "data/gen5/unapproved.xpt", "role": "input"})
    _save(root, document)
    with pytest.raises(artifacts.ArtifactError, match="unapproved private input"):
        _stage(private_specimen, tmp_path)
    document["runs"][0]["dependencies"][-1]["path"] = "data/gen5/unrelated.pem"
    policy = artifacts.read_json(root / POLICY)
    policy["runs"]["private-fixture"]["data/gen5/unrelated.pem"] = {"kind": "file"}
    _json(root / POLICY, policy)
    _save(root, document)
    with pytest.raises(artifacts.ArtifactError, match="extension"):
        _stage(private_specimen, tmp_path)


def test_declared_private_hash_is_mandatory_when_present(private_specimen, tmp_path):
    root, document, _, _ = private_specimen
    document["runs"][0]["dependencies"][0]["sha256"] = "0" * 64
    _save(root, document)
    with pytest.raises(artifacts.ArtifactError, match="identity changed"):
        _stage(private_specimen, tmp_path)


def test_corpus_is_only_top_level_authorized_extension(private_specimen, tmp_path, monkeypatch):
    root, document, corpus, _ = private_specimen
    policy = artifacts.read_json(root / POLICY)
    policy["runs"]["private-fixture"] = {"data/gen5": {"kind": "directory", "glob": "*.xpt"}}
    _json(root / POLICY, policy)
    document["runs"][0]["dependencies"][0] = {"path": "data/gen5", "kind": "directory", "role": "input"}
    _save(root, document)
    (corpus / "second.xpt").write_text("declared corpus member")
    (corpus / "unrelated.txt").write_text("MUST NOT READ")
    (corpus / "nested").mkdir()
    (corpus / "nested/hidden.xpt").write_text("MUST NOT READ")
    original_open = Path.open

    def guarded(path, *args, **kwargs):
        assert path not in {corpus / "unrelated.txt", corpus / "nested/hidden.xpt"}
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded)
    directory, receipt = _stage(private_specimen, tmp_path)
    private = {entry["path"] for entry in receipt["identity"]["dependencies"] if entry.get("private")}
    assert private == {PRIVATE, "data/gen5/second.xpt"}
    assert sorted(path.name for path in (directory / "workspace/data/gen5").iterdir()) == ["declared.xpt", "second.xpt"]
    assert receipt["success"]


def test_matching_corpus_symlink_is_not_followed(private_specimen, tmp_path):
    root, document, corpus, _ = private_specimen
    policy = artifacts.read_json(root / POLICY)
    policy["runs"]["private-fixture"] = {"data/gen5": {"kind": "directory", "glob": "*.xpt"}}
    _json(root / POLICY, policy)
    document["runs"][0]["dependencies"][0] = {"path": "data/gen5", "kind": "directory", "role": "input"}
    _save(root, document)
    (corpus / "escape.xpt").symlink_to(tmp_path / "outside.xpt")
    with pytest.raises(artifacts.ArtifactError, match="symlink"):
        _stage(private_specimen, tmp_path)


def test_repaired_qpcr_requires_opt_in_but_not_an_ambient_env_override(private_specimen, tmp_path, monkeypatch):
    root, document, corpus, env_file = private_specimen
    repaired = "data/qpcr_repaired/replicate.xlsx"
    (root / repaired).parent.mkdir()
    (root / repaired).write_bytes((corpus / "declared.xpt").read_bytes())
    script = root / "scripts/run.py"
    script.write_text(script.read_text().replace(PRIVATE, repaired))
    document["runs"][0]["dependencies"][0]["path"] = repaired
    document["runs"][0]["reproduction"]["environment"] = {"YSTWIN_QPCR": "{workspace}/data/qpcr_repaired"}
    _save(root, document)
    policy = artifacts.read_json(root / POLICY)
    policy["sources"] = {"data/qpcr_repaired": {"environment": "YSTWIN_QPCR", "suffix": ".xlsx", "resolver": "qpcr_dir", "repository_root": "data/qpcr_repaired"}}
    policy["runs"]["private-fixture"] = {repaired: {"kind": "file"}}
    _json(root / POLICY, policy)
    monkeypatch.setenv("YSTWIN_QPCR", str(tmp_path / "unapproved"))
    before = _stat_bytes(root / repaired)
    with pytest.raises(artifacts.ArtifactError, match="private-inputs-from-env"):
        artifacts.reproduce_run(root, "private-fixture", tmp_path / "denied")
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["success"]
    assert (directory / "workspace" / repaired).is_file()
    assert _stat_bytes(root / repaired) == before


def test_producer_cannot_write_original_or_staged_inputs(private_specimen, tmp_path):
    root, _, corpus, _ = private_specimen
    script = root / "scripts/run.py"
    script.write_text(script.read_text() + f"Path({PRIVATE!r}).write_text('mutated')\n")
    before = _stat_bytes(corpus / "declared.xpt")
    _, receipt = _stage(private_specimen, tmp_path)
    assert receipt["exit_code"] != 0 and not receipt["execution_complete"]
    assert _stat_bytes(corpus / "declared.xpt") == before


def test_receipt_ci_scope_skips_only_private_bytes_and_never_reads_originals(private_specimen, tmp_path, monkeypatch):
    root, _, corpus, _ = private_specimen
    directory, receipt, result = _adopt(private_specimen, tmp_path)
    (corpus / "declared.xpt").unlink()
    original_open = Path.open

    def guarded(path, *args, **kwargs):
        assert not path.is_relative_to(corpus)
        assert not path.is_relative_to(directory)
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded)
    rows = artifacts.audit_registry(root)
    output = _artifact_row(rows)
    assert output["status"] == "PASS"
    assert output["actual"] == "verified_prior_private_execution"
    assert "not a public replay" in output["detail"]
    skipped = [row for row in rows if row["status"] == "SKIP"]
    assert len(skipped) == 1
    assert skipped[0]["actual"] == "private_bytes_not_reverified"
    assert PRIVATE in skipped[0]["check"]
    private_sha = next(entry["sha256"] for entry in receipt["identity"]["dependencies"] if entry["path"] == PRIVATE)
    assert private_sha in skipped[0]["detail"]
    assert result["receipt_id"] in artifacts.load_registry(root)["receipts"]


@pytest.mark.parametrize("change", ["public_missing", "public_hash", "output", "code", "policy", "receipt", "unknown_private"])
def test_ci_receipt_scope_never_excuses_tampering_or_public_input_gaps(private_specimen, tmp_path, change):
    root, _, corpus, _ = private_specimen
    _, _, result = _adopt(private_specimen, tmp_path)
    (corpus / "declared.xpt").unlink()
    if change == "public_missing":
        (root / PUBLIC).unlink()
    elif change == "public_hash":
        (root / PUBLIC).write_text("changed public calibration")
    elif change == "output":
        (root / OUTPUT).write_text("id,value,unit\na,9,mg\n")
    elif change == "code":
        source = root / "scripts/run.py"
        source.write_text(source.read_text() + "\n")
    elif change == "policy":
        (root / POLICY).write_text((root / POLICY).read_text() + "\n")
    else:
        document = artifacts.read_json(root / artifacts.REGISTRY_PATH)
        receipt = document["receipts"][result["receipt_id"]]
        if change == "receipt":
            receipt["execution"]["exit_code"] = 9
        else:
            entry = next(entry for entry in receipt["identity"]["dependencies"] if entry["path"] == PRIVATE)
            entry["path"] = "data/gen5/unknown.xpt"
            digest = artifacts.canonical_digest(receipt["identity"])
            for key in ("before_binding_sha256", "after_binding_sha256", "expected_binding_sha256"):
                receipt["execution"][key] = digest
            receipt["execution_sha256"] = artifacts.canonical_digest(receipt["execution"])
            new_id = artifacts.canonical_digest(receipt)
            document["receipts"] = {new_id: receipt}
            document["runs"][0]["verification"]["artifacts"][OUTPUT]["receipt_id"] = new_id
        _save(root, document)
    rows = artifacts.audit_registry(root)
    assert _artifact_row(rows)["status"] == "FAIL", rows


def test_private_unavailability_without_prior_execution_remains_blocking(private_specimen):
    root, document, corpus, _ = private_specimen
    (corpus / "declared.xpt").unlink()
    document["runs"][0]["status"] = "verified"
    _save(root, document)
    assert _artifact_row(artifacts.audit_registry(root))["status"] == "FAIL"


@pytest.mark.parametrize("change", ["private_bytes", "staged_bytes", "child", "private_policy", "private_metadata"])
def test_private_review_and_adoption_revalidate_sources_and_staged_binding(private_specimen, tmp_path, change):
    root, _, corpus, env_file = private_specimen
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["success"]
    if change == "private_bytes":
        (corpus / "declared.xpt").write_text("changed")
    elif change == "staged_bytes":
        (directory / "workspace" / PRIVATE).write_text("changed")
    elif change == "child":
        ledger = artifacts.read_json(directory / "child_execution.json")
        ledger["exit_code"] = 9
        _json(directory / "child_execution.json", ledger)
    elif change == "private_policy":
        (root / POLICY).write_text((root / POLICY).read_text() + "\n")
    else:
        del receipt["identity"]["private_inputs"]
        _json(directory / "reproduction.json", receipt)
    with pytest.raises(artifacts.ArtifactError):
        artifacts.review_template(root, directory / "reproduction.json", private_inputs_from_env=env_file)
    with pytest.raises(artifacts.ArtifactError):
        artifacts.adopt_reproduction(root, directory / "reproduction.json", private_inputs_from_env=env_file)


def test_private_validation_is_opt_in_and_changed_comparisons_are_safe(private_specimen, tmp_path):
    root, _, corpus, env_file = private_specimen
    (corpus / "declared.xpt").write_text(CANARY + "\nid,value,unit\na,3,mg\nb,,mg\n")
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["execution_complete"] and not receipt["success"]
    with pytest.raises(artifacts.ArtifactError, match="private-inputs-from-env"):
        artifacts.review_template(root, directory / "reproduction.json")
    review = artifacts.review_template(root, directory / "reproduction.json", private_inputs_from_env=env_file)
    assert review["artifacts"][OUTPUT]["observed_changed_dimensions"] == ["values"]
    assert CANARY not in json.dumps(receipt)
    assert "changed_cell_samples" not in receipt["comparisons"][OUTPUT]


def test_private_cli_has_explicit_flag_and_does_not_source_shell(private_specimen, tmp_path, monkeypatch, capsys):
    root, _, _, env_file = private_specimen
    spec = importlib.util.spec_from_file_location("_private_regenerate_cli", ROOT / "scripts/regenerate_current_artifacts.py")
    cli = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    monkeypatch.setattr(cli, "ROOT", root)
    code = cli.main(["--run", "private-fixture", "--output-dir", str(tmp_path / "cli-stage"),
                     "--private-inputs-from-env", str(env_file)])
    assert code == 0
    result = json.loads(capsys.readouterr().out)
    assert result["execution_complete"]
    assert CANARY not in json.dumps(result)


def test_shipped_policy_contains_only_scoped_logical_inputs():
    policy = artifacts.read_json(ROOT / POLICY)
    assert set(policy["runs"]) == {"afl-circuit", "autofluorescence-measured", "gain-linearity", "g4"}
    # A leak check, not a key check: the key-set assertion above inspects no VALUE, so it
    # cannot catch a developer's absolute path serialised into the shipped policy.
    text = json.dumps(policy)
    # Segments, not literals: audit_claims.py forbids a literal home path in a tracked file,
    # and spelling one here to test for it would trip that gate instead of this one.
    leaked = [root for root in (f"/{'Users'}/", f"/{'home'}/", "/private/tmp/") if root in text]
    assert not leaked, leaked
    registry = artifacts.load_registry(ROOT)
    for run_id, entries in policy["runs"].items():
        run = registry["runs"][run_id]
        declared = {entry["path"] for entry in run["dependencies"]}
        assert set(entries) <= declared
        assert set(artifacts.private_input_declarations(ROOT, run)) == set(entries)
    assert policy["runs"]["gain-linearity"] == {"data/gen5": {"kind": "directory", "glob": "*.xpt"}}


def test_normalized_source_declarations_can_be_added_without_resigning_a_receipt(private_specimen, tmp_path):
    root, document, _, env_file = private_specimen
    directory, receipt = _stage(private_specimen, tmp_path)
    run = artifacts.load_registry(root)["runs"]["private-fixture"]
    declarations = artifacts.private_input_declarations(root, run)
    assert declarations == {PRIVATE: {"environment": "YSTWIN_GEN5_XPT", "resolver": "gen5_xpt_dir", "relative_path": "declared.xpt"}}
    document["runs"][0]["dependencies"][0]["source"] = declarations[PRIVATE]
    _save(root, document)
    _, _, verified = artifacts.validate_reproduction(root, directory / "reproduction.json", private_inputs_from_env=env_file)
    assert artifacts.canonical_digest(receipt) == artifacts.canonical_digest(verified)


@pytest.mark.parametrize("mode", ["missing_policy", "unapproved_source_type", "staged_bypass"])
def test_private_permission_cannot_be_inferred_from_a_path_alone(private_specimen, tmp_path, mode):
    root, document, _, _ = private_specimen
    if mode == "missing_policy":
        (root / POLICY).unlink()
    elif mode == "unapproved_source_type":
        policy = artifacts.read_json(root / POLICY)
        policy["sources"]["data/gen5"]["environment"] = "UNRELATED_SECRET"
        _json(root / POLICY, policy)
    else:
        with pytest.raises(artifacts.ArtifactError):
            artifacts.run_identity(root, document["runs"][0], staged=True)
        return
    with pytest.raises(artifacts.ArtifactError):
        _stage(private_specimen, tmp_path)


@pytest.mark.parametrize("mode", ["permissive_directory", "symlink_directory", "repository_directory"])
def test_private_staging_rejects_unprotected_locations(private_specimen, tmp_path, mode):
    root, _, _, env_file = private_specimen
    target = tmp_path / "target"
    if mode == "permissive_directory":
        target.mkdir(mode=0o755)
        target.chmod(0o755)
    elif mode == "symlink_directory":
        target.symlink_to(tmp_path, target_is_directory=True)
    else:
        target = root / "private-stage"
    with pytest.raises(artifacts.ArtifactError):
        artifacts.reproduce_run(root, "private-fixture", target, private_inputs_from_env=env_file)


@pytest.mark.parametrize("operation", ["unlink", "rename", "directory_read", "network", "subprocess"])
def test_private_child_cannot_mutate_originals_or_read_other_roots_or_upload(private_specimen, tmp_path, operation):
    root, _, corpus, _ = private_specimen
    original = corpus / "declared.xpt"
    script = root / "scripts/run.py"
    actions = {
        "unlink": f"Path({str(original)!r}).unlink()\n",
        "rename": f"Path({str(original)!r}).rename({str(corpus / 'renamed.xpt')!r})\n",
        "directory_read": f"import os\nos.open({str(corpus)!r}, os.O_RDONLY | os.O_DIRECTORY)\n",
        "network": "import socket\nsocket.socket().connect(('127.0.0.1', 9))\n",
        "subprocess": "import subprocess\nsubprocess.run(['false'])\n",
    }
    script.write_text(script.read_text() + actions[operation])
    before = _stat_bytes(original)
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["exit_code"] != 0 and not receipt["execution_complete"]
    assert _stat_bytes(original) == before
    assert str(corpus) not in json.dumps(receipt)
    assert CANARY not in (directory / "stderr.txt").read_text()


def test_private_malformed_output_metadata_is_not_copied_into_receipts(private_specimen, tmp_path):
    _, _, corpus, _ = private_specimen
    (corpus / "declared.xpt").write_text(CANARY + f"\nid,value,unit\na,{CANARY},mg\n")
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["exit_code"] == 0 and not receipt["execution_complete"]
    assert CANARY not in json.dumps(receipt)
    assert CANARY not in (directory / "reproduction.json").read_text()


def test_private_negative_case_changes_keep_exact_review_obligations(private_specimen, tmp_path):
    root, document, corpus, env_file = private_specimen
    (root / OUTPUT).write_text("id,value,unit,status\na,1,mg,ok\nb,,mg,failed\n")
    (corpus / "declared.xpt").write_text(CANARY + "\nid,value,unit,status\na,1,mg,ok\n")
    document["runs"][0]["artifacts"][0]["comparison"]["negative_cases"] = {"status": {"values": ["failed"]}}
    _save(root, document)
    directory, receipt = _stage(private_specimen, tmp_path)
    assert receipt["execution_complete"] and not receipt["success"]
    changes = receipt["comparisons"][OUTPUT]["negative_case_changes"]
    assert len(changes) == 1 and len(changes[0]["change_sha256"]) == 64
    review = artifacts.review_template(root, directory / "reproduction.json", private_inputs_from_env=env_file)
    assert review["artifacts"][OUTPUT]["required_negative_case_changes_sha256"] == artifacts.canonical_digest(changes)
    assert review["artifacts"][OUTPUT]["approved_negative_case_changes_sha256"] is None
    with pytest.raises(artifacts.ArtifactError, match="reviewed acceptance"):
        artifacts.adopt_reproduction(root, directory / "reproduction.json", accept_scientific_changes=True,
                                     private_inputs_from_env=env_file)


def test_private_corpus_roster_change_invalidates_review_and_nonmatching_receipt_member_blocks_ci(private_specimen, tmp_path):
    root, document, corpus, env_file = private_specimen
    policy = artifacts.read_json(root / POLICY)
    policy["runs"]["private-fixture"] = {"data/gen5": {"kind": "directory", "glob": "*.xpt"}}
    _json(root / POLICY, policy)
    document["runs"][0]["dependencies"][0] = {"path": "data/gen5", "kind": "directory", "role": "input"}
    _save(root, document)
    directory, _, result = _adopt(private_specimen, tmp_path)
    (corpus / "extra.xpt").write_text("new corpus member")
    with pytest.raises(artifacts.ArtifactError, match="stale"):
        artifacts.review_template(root, directory / "reproduction.json", private_inputs_from_env=env_file)
    document = artifacts.read_json(root / artifacts.REGISTRY_PATH)
    receipt = document["receipts"][result["receipt_id"]]
    member = next(entry for entry in receipt["identity"]["dependencies"] if entry.get("private"))
    member["path"] = "data/gen5/unrelated.pem"
    digest = artifacts.canonical_digest(receipt["identity"])
    for key in ("before_binding_sha256", "after_binding_sha256", "expected_binding_sha256"):
        receipt["execution"][key] = digest
    receipt["execution_sha256"] = artifacts.canonical_digest(receipt["execution"])
    new_id = artifacts.canonical_digest(receipt)
    document["receipts"] = {new_id: receipt}
    document["runs"][0]["verification"]["artifacts"][OUTPUT]["receipt_id"] = new_id
    _save(root, document)
    assert _artifact_row(artifacts.audit_registry(root))["status"] == "FAIL"
