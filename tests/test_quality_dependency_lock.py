from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import tomllib

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
import pytest

from ystwin.artifacts import load_registry


ROOT = Path(__file__).resolve().parents[1]


def test_quality_lock_binds_direct_and_optional_runtime_dependencies():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    records = (ROOT / "requirements-quality.txt").read_text().replace("\\\n", " ").splitlines()
    locked = {}
    for record in records:
        if not record.strip():
            continue
        requirement, *hashes = record.split("--hash=sha256:")
        parsed = Requirement(requirement.strip())
        assert hashes, f"unhashed distribution: {parsed.name}"
        assert all(len(digest.strip()) == 64 for digest in hashes)
        assert len(parsed.specifier) == 1
        assert next(iter(parsed.specifier)).operator == "=="
        locked.setdefault(canonicalize_name(parsed.name), set()).add(str(parsed.specifier))
    direct = list(project["project"]["dependencies"])
    for extra in project["tool"]["ystwin"]["quality"]["extras"]:
        direct.extend(project["project"]["optional-dependencies"][extra])
    for item in direct:
        required = Requirement(item)
        assert str(required.specifier) in locked[canonicalize_name(required.name)]


def test_pytfa_runtime_explicitly_declares_its_distutils_provider():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    thermo = {canonicalize_name(Requirement(item).name): Requirement(item)
              for item in project["project"]["optional-dependencies"]["thermo"]}
    assert "setuptools" in thermo
    assert str(thermo["setuptools"].specifier) == "==82.0.0"


def test_frozen_parser_schema_validator_is_a_pinned_development_dependency():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    dev = {canonicalize_name(Requirement(item).name): Requirement(item)
           for item in project["project"]["optional-dependencies"]["dev"]}
    assert "jsonschema" in dev
    assert len(dev["jsonschema"].specifier) == 1
    assert next(iter(dev["jsonschema"].specifier)).operator == "=="


@pytest.fixture
def gitignore_repo(tmp_path):
    subprocess.run(["git", "init", "--quiet", str(tmp_path)], check=True, capture_output=True)
    shutil.copyfile(ROOT / ".gitignore", tmp_path / ".gitignore")
    return tmp_path


def test_quality_lock_is_not_hidden_by_the_scratch_file_ignore_rule(gitignore_repo):
    result = subprocess.run(["git", "check-ignore", "--no-index", "--quiet", "requirements-quality.txt"],
                            cwd=gitignore_repo, capture_output=True)
    assert result.returncode == 1, "The quality dependency lock must be eligible for version control"


def test_registered_current_artifacts_are_not_hidden_by_output_ignore_rules(gitignore_repo):
    registry = load_registry(ROOT)
    paths = [item["path"] for run in registry["runs"].values() if run["lifecycle"] == "current"
             for item in run["artifacts"]]
    result = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"],
                            input="\n".join(paths) + "\n", cwd=gitignore_repo, capture_output=True, text=True)
    assert result.returncode == 1, f"Registered current artifacts would be missing from a clone:\n{result.stdout}"


def test_output_exceptions_do_not_allow_unregistered_private_or_scratch_files(gitignore_repo):
    paths = ["scratch.txt", "outputs/private.json", "outputs/unregistered.npz",
             "outputs/partial_order_v2/raw.csv", "outputs/state_allocation/private.json",
             "outputs/state_allocation/synthetic_demo/raw.csv"]
    result = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"],
                            input="\n".join(paths) + "\n", cwd=gitignore_repo, capture_output=True, text=True)
    assert result.returncode == 0
    assert set(result.stdout.splitlines()) == set(paths)
