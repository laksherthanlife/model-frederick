from __future__ import annotations

from copy import deepcopy
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tomllib
from types import SimpleNamespace
import urllib.error
import xml.etree.ElementTree as ET

import pytest
import yaml


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import conftest as shared
import provision_public_inputs as inputs
import verify_quality as quality


ROOT = Path(__file__).resolve().parents[1]
pytest_plugins = ["pytester"]


def _git(root, *arguments):
    return subprocess.run(["git", *arguments], cwd=root, capture_output=True, text=True, check=True)


@pytest.fixture
def repository(tmp_path):
    root = tmp_path / "checkout"
    root.mkdir()
    _git(root, "init", "--quiet")
    (root / ".gitignore").write_text("data/external/\n", encoding="utf-8")
    _git(root, "add", ".gitignore")
    return root.resolve()


def _manifest(payload=b"verified public input\n", path="data/external/input.tsv"):
    commit = "a" * 40
    return {
        "schema_version": 1,
        "sources": {
            "checkout": {"kind": "checkout"},
            "public": {"kind": "git", "repository": "example/scientific-data", "commit": commit,
                       "base_url": f"https://raw.githubusercontent.com/example/scientific-data/{commit}/"},
        },
        "files": [{"path": path, "source": "public", "remote_path": "input.tsv",
                   "sha256": hashlib.sha256(payload).hexdigest(), "size_bytes": len(payload)}],
        "collections": [],
        "environment": {},
        "optional_private_inputs": [],
    }


def _opener(monkeypatch, payload):
    calls = []

    def open_request(request, timeout):
        calls.append((request.full_url, timeout))
        return io.BytesIO(payload)

    monkeypatch.setattr(inputs.urllib.request, "build_opener",
                        lambda *args: SimpleNamespace(open=open_request))
    return calls


def _no_network(*args, **kwargs):
    raise AssertionError("this check must not access the network")


def test_workflow_is_one_unfiltered_blocking_job_with_pinned_actions_and_runtime():
    workflow = yaml.load((ROOT / ".github/workflows/quality.yml").read_text(), Loader=yaml.BaseLoader)
    assert set(workflow["on"]) == {"push", "pull_request", "workflow_dispatch"}
    assert all(value == "" for value in workflow["on"].values())
    assert set(workflow["jobs"]) == {"quality"}
    job = workflow["jobs"]["quality"]
    assert job["runs-on"] == "ubuntu-24.04"
    assert "continue-on-error" not in job
    assert workflow["permissions"] == {"contents": "read"}
    steps = job["steps"]
    assert all("continue-on-error" not in step for step in steps)
    actions = [step for step in steps if "uses" in step]
    for step in actions:
        action, revision = step["uses"].split("@")
        assert action.startswith("actions/")
        assert len(revision) == 40 and all(c in "0123456789abcdef" for c in revision)
    checkout = next(step for step in steps if step.get("uses", "").startswith("actions/checkout@"))
    assert checkout["with"]["fetch-depth"] == "0"
    setup = next(step for step in steps if step.get("uses", "").startswith("actions/setup-python@"))
    config = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert setup["with"]["python-version"] == config["tool"]["ystwin"]["quality"]["python-version"] == "3.14.2"
    runtime = next(step for step in steps if step.get("id") == "source_runtimes")
    pins = config["tool"]["ystwin"]["quality"]["source-runtimes"]
    assert pins["octave"] == {"linux-version": "8.4.0", "ubuntu-package": "octave",
                              "ubuntu-package-version": "8.4.0-1build5"}
    assert pins["pdftotext"] == {"linux-version": "24.02.0", "ubuntu-package": "poppler-utils",
                                 "ubuntu-package-version": "24.02.0-1ubuntu9.9"}
    assert "sudo apt-get update" in runtime["run"]
    assert "sudo apt-get install --yes --no-install-recommends" in runtime["run"]
    for pin in pins.values():
        assert f"{pin['ubuntu-package']}={pin['ubuntu-package-version']}" in runtime["run"]
    assert "||" not in runtime["run"] and "--allow-unauthenticated" not in runtime["run"]
    check_step = next(step for step in steps if step.get("id") == "checks")
    assert steps.index(runtime) < steps.index(check_step)
    assert check_step["if"] == "${{ !cancelled() }}"
    run = check_step["run"]
    assert "verify_quality.py --install" in run
    assert '--output-dir "$RUNNER_TEMP/ystwin-quality"' in run
    assert "||" not in run and "--restamp" not in run
    upload = next(step for step in steps if step.get("uses", "").startswith("actions/upload-artifact@"))
    assert upload["if"] == "failure()"
    assert upload["with"]["path"] == "${{ runner.temp }}/ystwin-quality"
    restore = next(step for step in steps if step.get("uses", "").startswith("actions/cache/restore@"))
    save = next(step for step in steps if step.get("uses", "").startswith("actions/cache/save@"))
    assert "data/public_inputs.json" in restore["with"]["key"]
    assert "!cancelled()" in save["if"]
    assert save["with"]["path"] == restore["with"]["path"]


def test_numerical_solver_and_all_quality_extras_are_exactly_pinned():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text())
    dependencies = config["project"]["dependencies"]
    assert {"numpy==2.4.1", "scipy==1.17.0", "pandas==2.3.3", "cobra==0.32.1",
            "optlang==1.9.1", "swiglpk==5.0.13", "sympy==1.14.0"} <= set(dependencies)
    extras = config["tool"]["ystwin"]["quality"]["extras"]
    assert extras == ["dev", "thermo", "sbml-audit"]
    optional = config["project"]["optional-dependencies"]
    assert "libroadrunner==2.10.0" in optional["sbml-audit"]
    assert {"component-contribution==0.7.0", "equilibrator-cache==0.7.1"} <= set(optional["thermo"])
    for requirement in dependencies + [item for extra in extras for item in optional[extra]]:
        assert requirement.count("==") == 1
        assert not any(symbol in requirement for symbol in "<>;,*")


def test_check_plan_runs_full_pytest_audits_replay_source_clis_and_final_integrity(tmp_path):
    output, cache = tmp_path / "reports", tmp_path / "cache"
    checks = quality.build_checks(ROOT, output, cache, install=True)
    commands = {check.name: check.argv for check in checks}
    assert commands["install-dependencies"][-3:] == ("--require-hashes", "-r", "requirements-quality.txt")
    assert commands["install"][1:3] == ("scripts/verify_quality.py", "--install-project")
    assert commands["install"][-2:] == ("--output-dir", str(output / "installation"))
    assert "-e" not in commands["install"]
    assert commands["pytest"][1:4] == ("-m", "pytest", "-ra")
    assert "--require-public-inputs" in commands["pytest"]
    assert "--junitxml" in commands["pytest"]
    assert not any(arg.startswith("tests/") or arg in {"-k", "-x", "--lf", "--ignore"}
                   for arg in commands["pytest"])
    assert commands["ruff"][1:] == ("-m", "ruff", "check", ".")
    for name in ("claims", "reproducibility"):
        assert commands[name][1] == f"scripts/audit_{name}.py"
        assert commands[name][2:] == ("--output-dir", str(output / name))
    assert commands["portable-native-replay"][1] == "scripts/replay_portable_evidence.py"
    assert "--historical" in commands["portable-native-replay"]
    assert commands["portable-native-replay"][-2:] == ("--output", str(output / "portable-native-replay.json"))
    assert commands["source-runtimes"][1:] == ("scripts/verify_quality.py", "--check-source-runtimes")
    assert set(quality.SOURCE_TESTS) <= set(commands["source-model-tests"])
    assert "tests/test_native_reference_models.py" in commands["source-model-tests"]
    assert "--require-public-inputs" in commands["source-model-tests"]
    assert "--gem-check" in commands["hog-source"]
    assert commands["ecmodel-source"][1] == "scripts/ecmodel_isoprenoid_prior.py"
    assert commands["stress-map-source"][1] == "scripts/audit_stress_map.py"
    assert commands["public-input-integrity"][2] == "--verify-only"
    assert [check.name for check in checks[-2:]] == ["worktree-after", "index-after"]
    assert commands["worktree-after"] == ("git", "diff", "--exit-code", "--")
    assert commands["index-after"] == ("git", "diff", "--cached", "--exit-code", "--")
    assert all("--restamp" not in check.argv for check in checks)


def test_quality_environment_keeps_committed_outputs_readable_and_reports_external(tmp_path, monkeypatch):
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path / "wrong-inputs"))
    monkeypatch.setenv("PYTEST_ADDOPTS", "-x -k narrowed")
    environment = quality.quality_environment(ROOT, tmp_path / "reports", tmp_path / "cache")
    assert "YSTWIN_OUTPUTS" not in environment
    assert environment["PYTEST_ADDOPTS"] == ""
    assert environment["CI"] == "true"
    for name in ("MPLCONFIGDIR", "RUFF_CACHE_DIR", "XDG_CACHE_HOME", "TMPDIR"):
        assert not Path(environment[name]).is_relative_to(ROOT)
    manifest = inputs.load_manifest()
    assert quality.PUBLIC_ENVIRONMENT == manifest["environment"]


def test_runner_aggregates_real_failures_and_still_runs_later_children(repository, tmp_path):
    checks = [
        quality.Check("fails", (sys.executable, "-c", "import sys; print('failure detail', file=sys.stderr); sys.exit(7)")),
        quality.Check("passes", (sys.executable, "-c", "print('later child ran')")),
    ]
    output = tmp_path / "reports"
    result = quality.run_checks(repository, output, checks, dict(os.environ))
    assert result["complete"] is True and result["passed"] is False
    assert [row["returncode"] for row in result["checks"]] == [7, 0]
    assert "failure detail" in (output / "fails.log").read_text()
    assert "later child ran" in (output / "passes.log").read_text()
    assert json.loads((output / "status.json").read_text()) == result


def test_runner_does_not_turn_a_launch_failure_or_signal_into_success(repository, tmp_path):
    checks = [
        quality.Check("missing", (str(tmp_path / "absent-executable"),)),
        quality.Check("signal", (sys.executable, "-c", "import os, signal; os.kill(os.getpid(), signal.SIGTERM)")),
        quality.Check("later", (sys.executable, "-c", "pass")),
    ]
    result = quality.run_checks(repository, tmp_path / "reports", checks, dict(os.environ))
    assert result["passed"] is False
    assert result["checks"][0]["returncode"] is None
    assert result["checks"][0]["launch_error"]
    assert result["checks"][1]["returncode"] == -signal.SIGTERM
    assert result["checks"][2]["returncode"] == 0


def test_final_tracked_diff_fails_after_a_child_mutates_a_tracked_file(repository, tmp_path):
    tracked = repository / "tracked.txt"
    tracked.write_text("original\n")
    _git(repository, "add", "tracked.txt")
    checks = [
        quality.Check("before", ("git", "diff", "--exit-code", "--")),
        quality.Check("mutation", (sys.executable, "-c", "from pathlib import Path; Path('tracked.txt').write_text('changed\\n')")),
        quality.Check("after", ("git", "diff", "--exit-code", "--")),
    ]
    result = quality.run_checks(repository, tmp_path / "reports", checks, dict(os.environ))
    assert [row["returncode"] for row in result["checks"]] == [0, 0, 1]
    assert result["passed"] is False
    assert "changed" in (tmp_path / "reports/after.log").read_text()


def test_runner_refuses_checkout_reports_and_reusing_existing_diagnostics(repository, tmp_path):
    checks = [quality.Check("one", (sys.executable, "-c", "pass"))]
    with pytest.raises(ValueError, match="outside"):
        quality.run_checks(repository, repository / "reports", checks, dict(os.environ))
    output = tmp_path / "reports"
    quality.run_checks(repository, output, checks, dict(os.environ))
    with pytest.raises(FileExistsError, match="fresh"):
        quality.run_checks(repository, output, checks, dict(os.environ))


@pytest.fixture
def package_repository(repository):
    (repository / "pyproject.toml").write_bytes((ROOT / "pyproject.toml").read_bytes())
    package = repository / "src/ystwin"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("VERSION = '0.1.0'\n")
    metadata = repository / "src/ystwin.egg-info"
    metadata.mkdir()
    for filename in ("PKG-INFO", "SOURCES.txt", "dependency_links.txt", "requires.txt", "top_level.txt"):
        (metadata / filename).write_text("stale generated metadata\n")
    _git(repository, "add", "src/ystwin.egg-info")
    return repository


@pytest.mark.parametrize("returncode", [0, 7])
def test_installation_uses_a_fresh_source_copy_without_touching_tracked_metadata(package_repository, tmp_path, monkeypatch, returncode):
    root = package_repository
    before = quality._metadata_hashes(root)
    real_run = subprocess.run
    invoked = []

    def installer(argv, **kwargs):
        if argv[:3] != [sys.executable, "-m", "pip"]:
            return real_run(argv, **kwargs)
        invoked.append(argv)
        stage = kwargs["cwd"]
        assert not stage.is_relative_to(root)
        assert (stage / "src/ystwin/__init__.py").read_bytes() == (root / "src/ystwin/__init__.py").read_bytes()
        assert not list(stage.rglob("*.egg-info"))
        assert "PYTHONPATH" not in kwargs["env"]
        assert {"--no-deps", "--no-build-isolation", "--no-index"} <= set(argv)
        assert "-e" not in argv
        return subprocess.CompletedProcess(argv, returncode)

    monkeypatch.setattr(quality.subprocess, "run", installer)
    output = tmp_path / "installation"
    assert quality.install_project(root, output) == returncode
    report = json.loads((output / "installation.json").read_text())
    assert report["pip_returncode"] == returncode
    assert report["passed"] is (returncode == 0)
    assert report["tracked_metadata_unchanged"] is True
    assert report["metadata_before"] == report["metadata_after"] == before
    assert len(report["omitted_generated_metadata"]) == 5
    assert len(invoked) == 1


def test_installation_detects_an_unexpected_tracked_metadata_write(package_repository, tmp_path, monkeypatch):
    root = package_repository
    real_run = subprocess.run

    def unsafe_installer(argv, **kwargs):
        if argv[:3] != [sys.executable, "-m", "pip"]:
            return real_run(argv, **kwargs)
        (root / "src/ystwin.egg-info/PKG-INFO").write_text("unexpected rewrite\n")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(quality.subprocess, "run", unsafe_installer)
    output = tmp_path / "installation"
    assert quality.install_project(root, output) == 1
    report = json.loads((output / "installation.json").read_text())
    assert report["pip_returncode"] == 0
    assert report["tracked_metadata_unchanged"] is False
    assert report["passed"] is False


def test_runtime_probe_fails_for_a_missing_numerical_dependency(monkeypatch):
    actual_version = quality.importlib.metadata.version

    def version(name):
        if name == "swiglpk":
            raise quality.importlib.metadata.PackageNotFoundError(name)
        return actual_version(name)

    monkeypatch.setattr(quality.importlib.metadata, "version", version)
    report = quality.check_environment(ROOT)
    assert report["passed"] is False
    assert report["dependencies"]["swiglpk"] == {"expected": "5.0.13", "actual": None}
    assert any("swiglpk" in error for error in report["errors"])


@pytest.fixture
def source_runtime_probe(monkeypatch):
    monkeypatch.setattr(quality.platform, "system", lambda: "Linux")
    monkeypatch.setattr(quality.shutil, "which", lambda name: f"/usr/bin/{name}")
    responses = {
        "octave": SimpleNamespace(returncode=0, stdout='{"version":"8.4.0"}\n', stderr=""),
        "pdftotext": SimpleNamespace(returncode=0, stdout="", stderr="pdftotext version 24.02.0\nCopyright\n"),
    }
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        response = responses[Path(argv[0]).name]
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr(quality.subprocess, "run", run)
    return responses, calls


def test_public_source_runtime_probe_records_versions_and_disables_local_octave_startup(source_runtime_probe):
    _, calls = source_runtime_probe
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is True
    assert report["platform"] == "Linux"
    assert set(report["runtimes"]) == {"octave", "pdftotext"}
    for name, version in (("octave", "8.4.0"), ("pdftotext", "24.02.0")):
        runtime = report["runtimes"][name]
        assert runtime["actual_version"] == runtime["expected_version"] == version
        assert runtime["returncode"] == 0
    octave = next(argv for argv, _ in calls if Path(argv[0]).name == "octave")
    assert {"--no-gui", "--quiet", "--no-history", "--no-init-file", "--no-site-file", "--eval"} <= set(octave)
    assert all(kwargs["timeout"] == 30 and kwargs["check"] is False for _, kwargs in calls)


@pytest.mark.parametrize("missing", ["octave", "pdftotext"])
def test_missing_public_source_runtime_is_a_failure_not_a_source_only_success(source_runtime_probe, monkeypatch, missing):
    _, calls = source_runtime_probe
    monkeypatch.setattr(quality.shutil, "which", lambda name: None if name == missing else f"/usr/bin/{name}")
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is False
    assert report["runtimes"][missing]["executable"] is None
    assert report["runtimes"][missing]["returncode"] is None
    assert any(missing in error and "required" in error for error in report["errors"])
    assert len(calls) == 1


@pytest.mark.parametrize("runtime", ["octave", "pdftotext"])
def test_linux_source_runtime_version_drift_fails(source_runtime_probe, runtime):
    responses, _ = source_runtime_probe
    responses[runtime] = SimpleNamespace(returncode=0, stdout='{"version":"99.0.0"}\n',
                                         stderr="pdftotext version 99.0.0\n")
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is False
    assert report["runtimes"][runtime]["actual_version"] == "99.0.0"
    assert any(runtime in error and "expected" in error for error in report["errors"])


@pytest.mark.parametrize("returncode", [7, -signal.SIGTERM])
def test_public_source_runtime_probe_keeps_failure_and_signal_statuses(source_runtime_probe, returncode):
    responses, calls = source_runtime_probe
    responses["octave"] = SimpleNamespace(returncode=returncode, stdout="", stderr="runtime failure")
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is False
    assert report["runtimes"]["octave"]["returncode"] == returncode
    assert report["runtimes"]["pdftotext"]["returncode"] == 0
    assert len(calls) == 2


@pytest.mark.parametrize("error", [OSError("runtime disappeared"), subprocess.TimeoutExpired("octave", 30)])
def test_public_source_runtime_probe_records_launch_and_timeout_failures(source_runtime_probe, error):
    responses, calls = source_runtime_probe
    responses["octave"] = error
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is False
    assert report["runtimes"]["octave"]["returncode"] is None
    assert report["runtimes"]["octave"]["error"]
    assert report["runtimes"]["pdftotext"]["returncode"] == 0
    assert len(calls) == 2


@pytest.mark.parametrize("runtime", ["octave", "pdftotext"])
def test_public_source_runtime_probe_refuses_unparseable_version_output(source_runtime_probe, runtime):
    responses, _ = source_runtime_probe
    responses[runtime] = SimpleNamespace(returncode=0, stdout="{}", stderr="unrecognized version banner")
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is False
    assert report["runtimes"][runtime]["actual_version"] is None


def test_local_source_runtime_versions_are_recorded_without_claiming_the_linux_pins(source_runtime_probe, monkeypatch):
    responses, _ = source_runtime_probe
    monkeypatch.setattr(quality.platform, "system", lambda: "Darwin")
    responses["octave"].stdout = '{"version":"11.3.0"}\n'
    responses["pdftotext"].stderr = "pdftotext version 26.04.0\n"
    report = quality.check_source_runtimes(ROOT)
    assert report["passed"] is True
    assert report["platform"] == "Darwin"
    assert all(row["expected_version"] is None for row in report["runtimes"].values())
    assert report["runtimes"]["octave"]["actual_version"] == "11.3.0"
    assert report["runtimes"]["pdftotext"]["actual_version"] == "26.04.0"


def test_source_runtime_cli_fails_and_reports_the_actual_child_status(source_runtime_probe, capsys):
    responses, _ = source_runtime_probe
    responses["octave"].returncode = 7
    assert quality.main(["--check-source-runtimes"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["passed"] is False
    assert report["runtimes"]["octave"]["returncode"] == 7


def test_download_is_pinned_verified_cached_and_reused_without_network(repository, tmp_path, monkeypatch):
    payload = b"verified public input\n"
    manifest = _manifest(payload)
    calls = _opener(monkeypatch, payload)
    cache = tmp_path / "cache"
    result = inputs.provision(repository, manifest, cache)
    assert result["passed"] is True
    assert len(calls) == 1
    assert calls[0][0] == manifest["sources"]["public"]["base_url"] + "input.tsv"
    assert (repository / manifest["files"][0]["path"]).read_bytes() == payload
    assert (cache / "sha256" / hashlib.sha256(payload).hexdigest()).read_bytes() == payload
    other = tmp_path / "other-checkout"
    other.mkdir()
    _git(other, "init", "--quiet")
    (other / ".gitignore").write_text("data/external/\n")
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    replayed = inputs.provision(other, manifest, cache)
    assert replayed["passed"] is True
    assert (other / manifest["files"][0]["path"]).read_bytes() == payload


@pytest.mark.parametrize("payload", [b"wrong data of same len", b"too much data" * 100])
def test_bad_download_never_creates_an_input_or_a_cache_entry(repository, tmp_path, monkeypatch, payload):
    manifest = _manifest()
    _opener(monkeypatch, payload)
    cache = tmp_path / "cache"
    result = inputs.provision(repository, manifest, cache)
    assert result["passed"] is False
    assert not (repository / manifest["files"][0]["path"]).exists()
    assert not (cache / "sha256" / manifest["files"][0]["sha256"]).exists()


def test_existing_bad_input_is_not_overwritten_by_a_correct_cache(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    entry = manifest["files"][0]
    destination = repository / entry["path"]
    destination.parent.mkdir(parents=True)
    destination.write_bytes(b"changed")
    cached = tmp_path / "cache/sha256" / entry["sha256"]
    cached.parent.mkdir(parents=True)
    cached.write_bytes(b"verified public input\n")
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    report = inputs.provision(repository, manifest, tmp_path / "cache")
    assert report["passed"] is False
    assert destination.read_bytes() == b"changed"
    assert "mismatch" in report["inputs"][0]["error"]


def test_corrupted_cache_is_reverified_not_trusted_or_silently_replaced(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    cached = tmp_path / "cache/sha256" / manifest["files"][0]["sha256"]
    cached.parent.mkdir(parents=True)
    cached.write_bytes(b"corrupted cache")
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    result = inputs.provision(repository, manifest, tmp_path / "cache")
    assert result["passed"] is False
    assert not (repository / manifest["files"][0]["path"]).exists()
    assert cached.read_bytes() == b"corrupted cache"


def test_checksum_pinned_mutable_source_downloads_matching_bytes(repository, tmp_path, monkeypatch):
    payload = b"verified public input\n"
    manifest = _manifest(payload)
    entry = manifest["files"][0]
    url = "https://example.invalid/current/input.tsv"
    manifest["sources"]["public"] = {"kind": "checksum", "urls": {entry["path"]: url}}
    inputs.validate_entry(entry, manifest["sources"])
    original = deepcopy(manifest)
    calls = _opener(monkeypatch, payload)
    result = inputs.provision(repository, manifest, tmp_path / "cache")
    assert result["passed"] is True
    assert calls == [(url, 120)]
    assert result["inputs"][0]["url"] == url
    assert (repository / entry["path"]).read_bytes() == payload
    assert manifest == original


def test_checksum_drift_fails_without_restamping_or_skipping_other_inputs(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    first = manifest["files"][0]
    second = {**first, "path": "data/external/present.tsv"}
    manifest["files"].append(second)
    manifest["sources"]["public"] = {
        "kind": "checksum", "urls": {entry["path"]: f"https://example.invalid/current/{Path(entry['path']).name}"
                                      for entry in manifest["files"]}}
    original = deepcopy(manifest)
    target = repository / second["path"]
    target.parent.mkdir(parents=True)
    target.write_bytes(b"verified public input\n")
    calls = _opener(monkeypatch, b"wrong data of same len")
    result = inputs.provision(repository, manifest, tmp_path / "cache")
    assert result["passed"] is False
    assert [row["status"] for row in result["inputs"]] == ["failed", "verified"]
    assert "sha256 mismatch" in result["inputs"][0]["error"]
    assert len(calls) == 1
    assert not (repository / first["path"]).exists()
    assert manifest == original


def test_authentication_failure_stops_network_operations_but_records_remaining_inputs(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    manifest["files"].append({**manifest["files"][0], "path": "data/external/second.tsv"})
    calls = []

    def denied(request, timeout):
        calls.append(request.full_url)
        raise urllib.error.HTTPError(request.full_url, 403, "Forbidden", {}, None)

    monkeypatch.setattr(inputs.urllib.request, "build_opener", lambda *args: SimpleNamespace(open=denied))
    report = inputs.provision(repository, manifest, tmp_path / "cache")
    assert report["passed"] is False
    assert len(calls) == 1
    assert len(report["inputs"]) == 2
    assert "already stopped" in report["inputs"][1]["error"]


def test_missing_tracked_input_is_not_repaired_by_provisioning(repository, tmp_path, monkeypatch):
    manifest = _manifest(path="data/tracked.tsv")
    manifest["files"][0]["source"] = "checkout"
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    report = inputs.provision(repository, manifest, tmp_path / "cache")
    assert report["passed"] is False
    assert "mandatory public input is absent" in report["inputs"][0]["error"]
    assert not (repository / "data/tracked.tsv").exists()


def test_downloader_cannot_write_unignored_or_tracked_paths(repository, tmp_path, monkeypatch):
    manifest = _manifest(path="data/unignored.tsv")
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    result = inputs.provision(repository, manifest, tmp_path / "cache")
    assert result["passed"] is False
    assert "gitignored" in result["inputs"][0]["error"]
    with pytest.raises(ValueError, match="tracked"):
        inputs.authorize_destination(repository, repository / "data/tracked.tsv", {"data/tracked.tsv"})


def test_verify_only_checks_required_absence_without_creating_cache_or_downloading(repository, tmp_path, monkeypatch):
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    report = inputs.provision(repository, _manifest(), verify_only=True)
    assert report["passed"] is False
    assert "absent" in report["inputs"][0]["error"]
    assert not (repository / "data").exists()


@pytest.mark.parametrize("value", ["../outside", "/absolute", "data/../outside", "data/./file", "data\\file", ""])
def test_input_paths_reject_escapes_and_noncanonical_aliases(value):
    with pytest.raises(ValueError):
        inputs.relative_path(value)


@pytest.mark.parametrize("revision", ["main", "latest", "v0.9.1", "a" * 39])
def test_download_sources_reject_mutable_git_refs(revision):
    manifest = _manifest()
    manifest["sources"]["public"]["commit"] = revision
    with pytest.raises(ValueError, match="digest"):
        inputs.source_url(manifest["sources"]["public"], manifest["files"][0])


def test_download_source_url_cannot_escape_its_declared_repository():
    manifest = _manifest()
    source = manifest["sources"]["public"]
    source["base_url"] = "https://elsewhere.invalid/latest/"
    with pytest.raises(ValueError, match="immutable source identity"):
        inputs.source_url(source, manifest["files"][0])


@pytest.mark.parametrize("url", ["http://example.invalid/input", "https://user@example.invalid/input",
                                 "https://example.invalid/input#fragment", "file:///tmp/input"])
def test_checksum_source_urls_still_require_credential_free_https(url):
    manifest = _manifest()
    entry = manifest["files"][0]
    source = {"kind": "checksum", "urls": {entry["path"]: url}}
    with pytest.raises(ValueError, match="HTTPS"):
        inputs.source_url(source, entry)


def test_checksum_source_never_replaces_a_missing_expected_digest():
    manifest = _manifest()
    entry = manifest["files"][0]
    manifest["sources"]["public"] = {"kind": "checksum", "urls": {entry["path"]: "https://example.invalid/input"}}
    entry["sha256"] = ""
    with pytest.raises(ValueError, match="digest"):
        inputs.validate_entry(entry, manifest["sources"])


def test_redirects_cannot_downgrade_https():
    with pytest.raises(inputs.NetworkFailure, match="HTTPS"):
        inputs.HTTPSOnlyRedirect().redirect_request(None, None, 302, "", {}, "http://example.invalid/input")


def test_symlinks_are_rejected_before_reading_or_writing_inputs(repository, tmp_path, monkeypatch):
    target = tmp_path / "outside"
    target.mkdir()
    (repository / "data").symlink_to(target, target_is_directory=True)
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    report = inputs.provision(repository, _manifest(), tmp_path / "cache")
    assert report["passed"] is False
    assert "symlink" in report["inputs"][0]["error"]
    assert not list(target.iterdir())


def test_compressed_model_and_its_uncompressed_source_are_both_verified(tmp_path):
    payload = b"<sbml>source model</sbml>"
    compressed = gzip.compress(payload)
    path = tmp_path / "model.xml.gz"
    path.write_bytes(compressed)
    entry = {"sha256": hashlib.sha256(compressed).hexdigest(), "size_bytes": len(compressed),
             "uncompressed_sha256": hashlib.sha256(payload).hexdigest()}
    assert inputs.verify_file(path, entry) == len(compressed)
    entry["uncompressed_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="uncompressed SBML"):
        inputs.verify_file(path, entry)


def test_manifest_covers_existing_source_bundles_and_all_external_validation_platforms():
    manifest = inputs.load_manifest()
    entries, errors = inputs.inventory(ROOT, manifest)
    assert errors == []
    paths = {entry["path"] for entry in entries}
    assert len(paths) == len(entries)
    source = json.loads((ROOT / "data/hog2013/sources.json").read_text())
    assert {f"data/hog2013/{entry['filename']}" for entry in source["assets"].values()} <= paths
    portable = json.loads((ROOT / quality.PORTABLE_MANIFEST).read_text())
    base = str(Path(quality.PORTABLE_MANIFEST).parent)
    assert {f"{base}/{entry['public_path']}" for entry in portable["records"]} <= paths
    assert next(item["sha256"] for item in manifest["collections"] if item["format"] == "portable_native") == quality.PORTABLE_SHA256
    from ystwin.analysis.frozen_runtime import RUNTIME_REGISTRY_PATH, RUNTIME_REGISTRY_SHA256

    assert next(item["sha256"] for item in manifest["files"] if item["path"] == RUNTIME_REGISTRY_PATH) == RUNTIME_REGISTRY_SHA256
    import fetch_external_validation_data as fetch

    assert {f"data/external/gasch2000/GSE18-GPL{gpl}_series_matrix.txt.gz" for gpl in fetch.PLATFORMS} <= paths
    assert {f"data/external/gasch2000/GPL{gpl}.annot.gz" for gpl in fetch.PLATFORMS if gpl not in fetch.NO_ANNOTATION} <= paths
    assert {f"data/external/sgd/regulation/{gene}.json" for gene in fetch.FACTORS} <= paths
    assert "data/external/sgd/gene_association.sgd.gaf.gz" in paths
    assert "data/kaggle/BioNumbers_Nov2024.csv" in paths
    assert {entry["path"] for entry in entries if entry.get("storage") == "equilibrator_cache"} == {"compounds.sqlite", "cc_params.npz"}
    expected_geo = {f"data/external/gasch2000/GSE18-GPL{gpl}_series_matrix.txt.gz": fetch.GEO_MATRIX.format(gpl=gpl)
                    for gpl in fetch.PLATFORMS}
    expected_geo.update({f"data/external/gasch2000/GPL{gpl}.annot.gz": fetch.GEO_ANNOT.format(gpl=gpl)
                         for gpl in fetch.PLATFORMS if gpl not in fetch.NO_ANNOTATION})
    assert manifest["sources"]["geo"]["urls"] == expected_geo
    assert manifest["sources"]["sgd_regulation"]["urls"] == {
        f"data/external/sgd/regulation/{gene}.json": fetch.SGD_REGULATION.format(gene=gene)
        for gene in fetch.FACTORS}
    assert manifest["sources"]["sgd_gaf"]["urls"] == {
        "data/external/sgd/gene_association.sgd.gaf.gz":
            "https://downloads.yeastgenome.org/curation/literature/archive/gene_association.sgd.20260825.gaf.gz"}
    checksum_inputs = [entry for entry in manifest["files"] if manifest["sources"][entry["source"]]["kind"] == "checksum"]
    assert len(checksum_inputs) == 50
    assert all(inputs.source_url(manifest["sources"][entry["source"]], entry) for entry in checksum_inputs)
    for entry in manifest["files"]:
        inputs.validate_entry(entry, manifest["sources"])


def test_manifest_hash_is_checked_before_expanding_a_collection(repository):
    manifest = _manifest()
    collection = repository / "data/source.json"
    collection.parent.mkdir()
    collection.write_text('{"assets": {"bad": {"filename": "../outside"}}}')
    manifest["collections"] = [{"path": "data/source.json", "sha256": "0" * 64, "format": "hog_sources"}]
    entries, errors = inputs.inventory(repository, manifest)
    assert len(entries) == 2
    assert len(errors) == 1 and "sha256 mismatch" in errors[0]


def test_public_resolver_override_cannot_hide_a_missing_or_different_input(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    entry = manifest["files"][0]
    path = repository / entry["path"]
    path.parent.mkdir(parents=True)
    path.write_bytes(b"verified public input\n")
    manifest["environment"] = {"YSTWIN_YEAST_GEM": entry["path"]}
    monkeypatch.setenv("YSTWIN_YEAST_GEM", str(tmp_path / "missing-model"))
    report = inputs.provision(repository, manifest, verify_only=True)
    assert report["passed"] is False
    assert any("YSTWIN_YEAST_GEM" in error for error in report["errors"])


def test_cli_returns_failure_and_writes_only_external_diagnostics(repository, tmp_path, monkeypatch):
    manifest = _manifest()
    path = repository / "manifest.json"
    path.write_text(json.dumps(manifest))
    monkeypatch.setattr(inputs, "download_verified", _no_network)
    output = tmp_path / "reports"
    assert inputs.main(["--root", str(repository), "--manifest", str(path), "--verify-only", "--output-dir", str(output)]) == 1
    report = json.loads((output / "public-inputs.json").read_text())
    assert report["passed"] is False
    assert not (repository / "data").exists()
    assert inputs.main(["--root", str(repository), "--manifest", str(path), "--verify-only",
                        "--output-dir", str(repository / "reports")]) == 1
    assert not (repository / "reports").exists()


@pytest.mark.parametrize("filename,fixtures,ancestry,private", [
    ("tests/test_synergy_parser_real.py", ["real_er_prelim"], [], True),
    ("tests/test_gain_linearity.py", ["xpt_dir"], [], True),
    ("tests/test_gain_linearity.py", [], [], False),
    ("tests/test_gen5_xpt.py", [], ["tests/test_gen5_xpt.py::TestTheWholeCorpus"], True),
    ("tests/test_gen5_xpt.py", [], ["tests/test_gen5_xpt.py::TestDecoding"], False),
    ("tests/test_product_panel.py", ["yeast_gem"], [], False),
    ("tests/test_pipeline_real.py", ["prepared", "public_er_prelim"], [], False),
    ("tests/test_pipeline_real.py", ["public_er_stress"], [], False),
])
def test_private_scope_is_explicit_and_does_not_exempt_public_or_synthetic_tests(filename, fixtures, ancestry, private):
    item = SimpleNamespace(path=ROOT / filename, fixturenames=fixtures,
                           listchain=lambda: [SimpleNamespace(nodeid=node) for node in ancestry])
    catalogue = inputs.load_manifest()["optional_private_inputs"]
    assert (shared._private_input_for(item, catalogue) is not None) is private


@pytest.mark.parametrize("ci_name", ["CI", "GITHUB_ACTIONS"])
@pytest.mark.parametrize("ci_value", ["1", "true", "TRUE"])
def test_private_local_environment_is_not_loaded_in_ci(tmp_path, monkeypatch, ci_name, ci_value):
    environment = tmp_path / ".ystwin.env"
    environment.write_text("YSTWIN_YEAST_GEM=private-machine-path\n")
    monkeypatch.setattr(shared, "_ENV_FILE", environment)
    monkeypatch.delenv("YSTWIN_YEAST_GEM", raising=False)
    for name in ("CI", "GITHUB_ACTIONS"):
        monkeypatch.setenv(name, ci_value if name == ci_name else "false")
    shared._load_local_env()
    assert "YSTWIN_YEAST_GEM" not in os.environ


def test_public_fixture_absence_fails_in_ci_but_remains_optional_locally(monkeypatch):
    monkeypatch.setenv("CI", "true")
    with pytest.raises(pytest.fail.Exception, match="mandatory public input"):
        shared._public_unavailable("missing model")
    monkeypatch.delenv("CI")
    monkeypatch.delenv("GITHUB_ACTIONS", raising=False)
    with pytest.raises(pytest.skip.Exception, match="missing model"):
        shared._public_unavailable("missing model")


@pytest.fixture
def policy_runner(pytester, monkeypatch):
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join((str(ROOT), str(ROOT / "src"))))
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    pytester.makeini("[pytest]\naddopts = --strict-markers\nmarkers =\n    private_data: optional original wet-lab input\n")
    pytester.makeconftest('''
        from types import SimpleNamespace
        import tests.conftest as shared
        pytest_plugins = ["tests.conftest"]
        def pytest_configure(config):
            shared._verify_public_inputs = lambda: SimpleNamespace(returncode=0, stdout="verified")
            shared._verify_public_source_runtimes = lambda: SimpleNamespace(returncode=0, stdout="verified runtimes")
    ''')
    return pytester


def test_ci_turns_an_unclassified_skip_into_a_failure_and_junit_failure(policy_runner):
    policy_runner.makepyfile("import pytest\ndef test_public():\n    pytest.skip('mandatory input absent')\n")
    result = policy_runner.runpytest_subprocess("--require-public-inputs", "-q", "--junitxml=junit.xml")
    assert result.ret == 1
    result.assert_outcomes(failed=1)
    xml = ET.parse(policy_runner.path / "junit.xml")
    assert len(xml.findall(".//failure")) == 1
    assert not xml.findall(".//skipped")


def test_ci_retains_only_explicit_optional_private_skips(policy_runner):
    policy_runner.makepyfile('''
        import pytest
        @pytest.mark.private_data("YSTWIN_GEN5_XPT")
        def test_private_original():
            pytest.skip("original instrument file absent")
    ''')
    result = policy_runner.runpytest_subprocess("--require-public-inputs", "-q", "--junitxml=junit.xml")
    assert result.ret == 0
    result.assert_outcomes(skipped=1)
    skip = ET.parse(policy_runner.path / "junit.xml").find(".//skipped")
    assert skip is not None
    assert "optional private input YSTWIN_GEN5_XPT" in skip.attrib["message"]


def test_ci_module_importorskip_is_a_collection_error_not_a_green_missing_dependency(policy_runner):
    policy_runner.makepyfile("import pytest\npytest.importorskip('ystwin_nonexistent_required_module')\ndef test_public(): pass\n")
    result = policy_runner.runpytest_subprocess("--require-public-inputs", "-q")
    assert result.ret == 2
    result.assert_outcomes(errors=1)
    assert "mandatory test module skipped" in result.stdout.str()


def test_preflight_public_input_failure_produces_nonzero_exit_and_junit(policy_runner):
    policy_runner.makeconftest('''
        from types import SimpleNamespace
        import tests.conftest as shared
        pytest_plugins = ["tests.conftest"]
        def pytest_configure(config):
            shared._verify_public_inputs = lambda: SimpleNamespace(returncode=7, stdout="missing public model")
    ''')
    policy_runner.makepyfile("def test_never_runs():\n    assert False\n")
    result = policy_runner.runpytest_subprocess("--require-public-inputs", "-q", "--junitxml=junit.xml")
    assert result.ret == 1
    assert "exit 7" in result.stdout.str() + result.stderr.str()
    assert (policy_runner.path / "junit.xml").is_file()


def test_preflight_public_source_runtime_failure_cannot_pass_as_source_only(policy_runner):
    policy_runner.makeconftest('''
        from types import SimpleNamespace
        import tests.conftest as shared
        pytest_plugins = ["tests.conftest"]
        def pytest_configure(config):
            shared._verify_public_inputs = lambda: SimpleNamespace(returncode=0, stdout="verified")
            shared._verify_public_source_runtimes = lambda: SimpleNamespace(returncode=7, stdout="required octave absent")
    ''')
    policy_runner.makepyfile("def test_source_only_fallback():\n    pass\n")
    result = policy_runner.runpytest_subprocess("--require-public-inputs", "-q", "--junitxml=junit.xml")
    assert result.ret == 1
    assert "Mandatory public source runtime verification failed (exit 7)" in result.stdout.str() + result.stderr.str()
    assert (policy_runner.path / "junit.xml").is_file()
    assert not ET.parse(policy_runner.path / "junit.xml").findall(".//skipped")


def test_local_optional_skip_contract_is_unchanged_without_strict_ci(policy_runner):
    policy_runner.makepyfile("import pytest\ndef test_local():\n    pytest.skip('local input absent')\n")
    result = policy_runner.runpytest_subprocess("-q")
    assert result.ret == 0
    result.assert_outcomes(skipped=1)
