"""Explicit original-source fixtures for the historical population releases.

The listed active source paths are selected from the registered Git commit.
Public custody, saved code, projections and outputs must independently match
their supplied bytes; Git never replaces changed evidence or repairs a tampered
admitted root. Production APIs retain strict current-root admission behavior.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
from types import ModuleType

import pytest

from ystwin.analysis import frozen_runtime as runtime
from ystwin.analysis import native_population_development as frozen
from ystwin.analysis import native_population_scoring as scoring


ROOT = Path(__file__).resolve().parents[1]
CODE_REF = "84561f51502b5b4bf66c052dc88c853038917934"
EXPORTER_PATH = "scripts/export_granados_training.py"
EXPORTER_SHA = "6447608e622c6bdd47e16f4f930f91f9297077aaccfbc5a0745b359b1bbd7159"
RELEASE_PATH = "data/native_law_v2/granados/development/phase_c_01/release_phase_c.py"
RELEASE_ASSET_PATH = RELEASE_PATH + ".source"
RELEASE_SHA = "68c31da88770e4718bc5668791edd38de6a054fa821aaa074a2ef48d6c676d8a"
SOURCE_ASSET_RELOCATIONS = {
    RELEASE_PATH: {"path": RELEASE_ASSET_PATH, "sha256": RELEASE_SHA, "size_bytes": 42204},
}
HELPER_PATH = "src/ystwin/analysis/native_population_development.py"
ACTIVE_SOURCE_PATHS = (
    HELPER_PATH,
    "scripts/run_native_population_development.py",
    "tests/test_native_population_development.py",
    EXPORTER_PATH,
    "src/ystwin/analysis/native_population_scoring.py",
    "scripts/score_native_population_development.py",
    "tests/test_native_population_scoring.py",
)
PUBLIC_PATHS = (
    scoring.PROTOCOL_PATH,
    "data/native_law_v2/granados/development",
    frozen.admission_contract()["projection_root"],
    "outputs/native_population_development",
)


def _check_bytes(content, sha256, size, path):
    assert len(content) == size, f"historical byte inventory mismatch: {path}"
    assert hashlib.sha256(content).hexdigest() == sha256, f"historical digest mismatch: {path}"
    return content


def _supplied_public_source(source_root, relative):
    path = frozen.safe_path(source_root, relative)
    relocation = SOURCE_ASSET_RELOCATIONS.get(relative)
    if relocation is None:
        return path.read_bytes()
    asset = frozen.safe_path(source_root, relocation["path"])
    supplied = [candidate for candidate in (path, asset) if candidate.exists()]
    assert len(supplied) == 1, f"historical source requires exactly one original or registered asset: {relative}"
    return _check_bytes(supplied[0].read_bytes(), relocation["sha256"], relocation["size_bytes"], relative)


def _original_population_files(source_root, source_repository=ROOT):
    source_root = Path(source_root).resolve(strict=True)
    repository = Path(source_repository).resolve(strict=True)
    tree = runtime._git(repository, "ls-tree", "-r", "-z", "--name-only", "--full-tree",
                        CODE_REF, "--", *PUBLIC_PATHS, *ACTIVE_SOURCE_PATHS)
    paths = [path.decode("utf-8") for path in tree.split(b"\0") if path]
    assert set(ACTIVE_SOURCE_PATHS) <= set(paths)
    blobs = runtime._git_blobs(repository, CODE_REF, paths)
    _check_bytes(blobs[scoring.MANIFEST_PATH], scoring.MANIFEST_SHA256, 8574, scoring.MANIFEST_PATH)
    manifest = frozen.strict_json(blobs[scoring.MANIFEST_PATH])
    request_ref = manifest["release_request"]
    request = frozen.strict_json(_check_bytes(blobs[request_ref["path"]], request_ref["sha256"],
                                              request_ref["bytes"], request_ref["path"]))
    for reference in request["executor_code"] + request["exporter_code"]:
        _check_bytes(blobs[reference["path"]], reference["sha256"], reference["bytes"], reference["path"])
    _check_bytes(blobs[EXPORTER_PATH], EXPORTER_SHA, 104189, EXPORTER_PATH)
    _check_bytes(blobs[RELEASE_PATH], RELEASE_SHA, 42204, RELEASE_PATH)
    for relative in ACTIVE_SOURCE_PATHS[-3:]:
        snapshot = ("outputs/native_population_development/native_training_01_phase_c_01_scoring_fix_v2/"
                    f"adapter_code/{Path(relative).name}")
        assert blobs[relative] == blobs[snapshot], f"historical adapter snapshot mismatch: {relative}"
    for relative, original in blobs.items():
        frozen.safe_path(source_root, relative)
        if relative not in ACTIVE_SOURCE_PATHS:
            supplied = _supplied_public_source(source_root, relative)
            _check_bytes(supplied, hashlib.sha256(original).hexdigest(), len(original), relative)
    return blobs


def _materialize_original_population(directory, source_root=ROOT, source_repository=ROOT):
    blobs = _original_population_files(source_root, source_repository)
    directory = Path(directory).resolve(strict=True)
    for relative, content in blobs.items():
        destination = frozen.safe_path(directory, relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as stream:
            stream.write(content)
    return blobs


@pytest.fixture(scope="module")
def historical_population_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("original_population_sources").resolve(strict=True)
    blobs = _materialize_original_population(root)
    yield root
    for relative, content in blobs.items():
        path = frozen.safe_path(root, relative)
        _check_bytes(path.read_bytes(), hashlib.sha256(content).hexdigest(), len(content), relative)


def _copy_population(tmp_path, historical_population_root):
    root = tmp_path.resolve(strict=True) / "population_copy"
    shutil.copytree(historical_population_root, root)
    return root


def _load_original_helper(root):
    path = frozen.safe_path(root, HELPER_PATH)
    content = _check_bytes(path.read_bytes(), scoring.MODULE_SHA256, 92464, HELPER_PATH)
    module = ModuleType("verified_original_population_helper")
    module.__file__ = str(path)
    exec(compile(content, str(path), "exec"), module.__dict__)
    return module


def _forbid_git(*args, **kwargs):
    raise AssertionError("strict admission must not replace changed artifacts from Git")


def test_fixture_uses_original_source_bytes_and_actual_public_release(historical_population_root):
    root = historical_population_root
    _check_bytes((root / EXPORTER_PATH).read_bytes(), EXPORTER_SHA, 104189, EXPORTER_PATH)
    _check_bytes((root / RELEASE_PATH).read_bytes(), RELEASE_SHA, 42204, RELEASE_PATH)
    assert not (root / RELEASE_ASSET_PATH).exists()
    bundle = scoring.inspect_phase_c(root)
    assert bundle.admission_pins.authority_root == str(root)
    assert bundle.release.manifest_reference["sha256"] == scoring.MANIFEST_SHA256
    assert bundle.prediction.reference["sha256"] == scoring.FREEZE_SHA256
    request = bundle.release.request
    for reference in request["executor_code"] + request["exporter_code"]:
        content = frozen.read_artifact(root, reference, bundle.protocol, "execution_code", parse=False)
        assert hashlib.sha256(content).hexdigest() == reference["sha256"]
    assert not (root / ".git").exists()
    assert not (root / "data/native_law_v2/granados/sealed").exists()
    assert not (root / "data/native_law_v2/granados/development/private").exists()


def test_original_helper_executes_only_after_source_verification(historical_population_root):
    root = historical_population_root
    original = _load_original_helper(root)
    phase_b = "data/native_law_v2/granados/development/release_manifest.json"
    refresh = "data/native_law_v2/granados/development/executor_refresh_01/executor_refresh_admission.json"
    before = original.check_custody_gate(root, root / scoring.PROTOCOL_PATH, root / phase_b)
    assert before["data_ready"] is True and before["executor_ready"] is False
    admitted = original.check_custody_gate(root, root / scoring.PROTOCOL_PATH, root / phase_b,
                                          executor_admission_path=root / refresh)
    assert admitted["ready"] is True and admitted["measurement_arrays_opened"] is False
    assert original.__file__ == str(root / HELPER_PATH)


@pytest.mark.parametrize("relative", [EXPORTER_PATH, HELPER_PATH, "tests/test_native_population_development.py"])
def test_strict_current_admission_rejects_changed_source_without_git_fallback(
    tmp_path, historical_population_root, monkeypatch, relative,
):
    root = _copy_population(tmp_path, historical_population_root)
    bundle = scoring.inspect_phase_c(root)
    path = root / relative
    path.write_bytes(path.read_bytes() + b"\n")
    monkeypatch.setattr(runtime, "_git_blobs", _forbid_git)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError)):
        scoring.inspect_phase_c(root)
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError)):
        scoring.revalidate_admission(bundle)
    gate = frozen.check_custody_gate(
        root, root / scoring.PROTOCOL_PATH,
        root / bundle.release.manifest["phase_b_manifest"]["path"],
        executor_admission_path=root / bundle.release.manifest["executor_refresh_admission"]["path"],
    )
    assert gate["ready"] is False and gate["measurement_arrays_opened"] is False


@pytest.mark.parametrize("relative", [scoring.MANIFEST_PATH, RELEASE_PATH,
                                      f"{frozen.admission_contract()['projection_root']}/development_responses.json",
                                      "data/native_law_v2/granados/development/phase_b/approval_event.json",
                                      "outputs/native_population_development/executor_refresh_01/code/native_population_development.py"])
def test_changed_supplied_historical_evidence_is_never_replaced_from_git(
    tmp_path, historical_population_root, relative,
):
    root = _copy_population(tmp_path, historical_population_root)
    path = root / relative
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(AssertionError, match="historical (byte inventory|digest) mismatch"):
        _original_population_files(root)


@pytest.mark.parametrize("relative", [EXPORTER_PATH, HELPER_PATH, RELEASE_PATH])
@pytest.mark.parametrize("same_size", [False, True])
def test_corrupt_git_source_is_rejected_before_fixture_execution(monkeypatch, relative, same_size):
    original = runtime._git_blobs

    def corrupt(*args, **kwargs):
        blobs = original(*args, **kwargs)
        content = blobs[relative]
        blobs[relative] = b"!" + content[1:] if same_size else content + b"\n"
        return blobs

    monkeypatch.setattr(runtime, "_git_blobs", corrupt)
    with pytest.raises(AssertionError, match="historical (byte inventory|digest) mismatch"):
        _original_population_files(ROOT)


def test_changed_materialized_helper_cannot_execute(tmp_path, historical_population_root):
    root = _copy_population(tmp_path, historical_population_root)
    path = root / HELPER_PATH
    path.write_bytes(path.read_bytes() + b"\nraise RuntimeError('unverified code executed')\n")
    with pytest.raises(AssertionError, match="historical byte inventory mismatch"):
        _load_original_helper(root)


@pytest.mark.parametrize("relative", [EXPORTER_PATH, scoring.MANIFEST_PATH, RELEASE_PATH])
def test_historical_sources_and_metadata_cannot_escape_through_symlinks(
    tmp_path, historical_population_root, relative,
):
    root = tmp_path.resolve(strict=True) / "aliased_population"
    root.mkdir()
    link = root / relative
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(historical_population_root / relative)
    with pytest.raises(frozen.DevelopmentError, match="symlink"):
        frozen.safe_path(root, relative)


def test_original_fixture_rejects_supplied_metadata_symlink(tmp_path, historical_population_root):
    root = tmp_path.resolve(strict=True) / "aliased_metadata"
    root.mkdir()
    (root / "data").symlink_to(historical_population_root / "data", target_is_directory=True)
    with pytest.raises(frozen.DevelopmentError, match="symlink"):
        _original_population_files(root)


def test_registered_source_asset_materializes_only_original_logical_path(tmp_path, historical_population_root):
    root = _copy_population(tmp_path, historical_population_root)
    (root / RELEASE_PATH).rename(root / RELEASE_ASSET_PATH)
    destination = tmp_path / "materialized"
    destination.mkdir()
    blobs = _materialize_original_population(destination, root)
    assert RELEASE_PATH in blobs and RELEASE_ASSET_PATH not in blobs
    _check_bytes((destination / RELEASE_PATH).read_bytes(), RELEASE_SHA, 42204, RELEASE_PATH)
    assert not (destination / RELEASE_ASSET_PATH).exists()
    assert scoring.inspect_phase_c(destination).release.manifest_reference["sha256"] == scoring.MANIFEST_SHA256
    with pytest.raises((scoring.ScoringError, frozen.DevelopmentError)):
        scoring.inspect_phase_c(root)


@pytest.mark.parametrize("same_size", [False, True])
def test_changed_source_asset_is_not_repaired_by_pristine_git(tmp_path, historical_population_root, same_size):
    root = _copy_population(tmp_path, historical_population_root)
    asset = root / RELEASE_ASSET_PATH
    (root / RELEASE_PATH).rename(asset)
    original = asset.read_bytes()
    asset.write_bytes(b"!" + original[1:] if same_size else original + b"\n")
    with pytest.raises(AssertionError, match="historical (byte inventory|digest) mismatch"):
        _original_population_files(root)
    assert not (root / RELEASE_PATH).exists()


@pytest.mark.parametrize("changed", [None, RELEASE_PATH, RELEASE_ASSET_PATH])
def test_original_and_relocated_source_are_ambiguous_even_if_one_is_valid(
    tmp_path, historical_population_root, changed,
):
    root = _copy_population(tmp_path, historical_population_root)
    (root / RELEASE_ASSET_PATH).write_bytes((root / RELEASE_PATH).read_bytes())
    if changed is not None:
        path = root / changed
        path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(AssertionError, match="exactly one original or registered asset"):
        _original_population_files(root)


@pytest.mark.parametrize("unregistered_alias", [False, True])
def test_missing_original_and_registered_asset_are_not_sourced_from_git(
    tmp_path, historical_population_root, unregistered_alias,
):
    root = _copy_population(tmp_path, historical_population_root)
    path = root / RELEASE_PATH
    if unregistered_alias:
        path.rename(path.with_suffix(".backup"))
    else:
        path.unlink()
    with pytest.raises(AssertionError, match="exactly one original or registered asset"):
        _original_population_files(root)


def test_registered_source_asset_rejects_symlink_even_with_valid_bytes(tmp_path, historical_population_root):
    root = _copy_population(tmp_path, historical_population_root)
    (root / RELEASE_PATH).unlink()
    (root / RELEASE_ASSET_PATH).symlink_to(historical_population_root / RELEASE_PATH)
    with pytest.raises(frozen.DevelopmentError, match="symlink"):
        _original_population_files(root)
