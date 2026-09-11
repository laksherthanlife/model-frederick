from __future__ import annotations

import base64
from copy import deepcopy
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
from types import SimpleNamespace

import cobra
import pandas as pd
import pytest

from ystwin import paths
from ystwin.fba.solver import configure, load_model


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("gem_version_diff", ROOT / "scripts/gem_version_diff.py")
gem_diff = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gem_diff)


@pytest.fixture
def toy():
    model = cobra.Model("growth_probe")
    carbon = cobra.Metabolite("carbon", compartment="c")
    oxygen = cobra.Metabolite("oxygen", compartment="c")
    glucose = cobra.Reaction("r_1714", lower_bound=-10.0, upper_bound=1000.0)
    glucose.add_metabolites({carbon: -1.0})
    o2 = cobra.Reaction("r_1992", lower_bound=-20.0, upper_bound=1000.0)
    o2.add_metabolites({oxygen: -1.0})
    biomass = cobra.Reaction("r_2111", lower_bound=0.0, upper_bound=1000.0)
    biomass.add_metabolites({carbon: -1.0})
    model.add_reactions([glucose, o2, biomass])
    model.objective = biomass
    configure(model)
    return model


def test_optimal_zero_is_real_growth_not_nan(toy):
    toy.reactions.r_2111.upper_bound = 0.0
    solution = toy.optimize()
    assert solution.status == "optimal"
    assert solution.objective_value == 0.0
    assert gem_diff._growth(toy, -1.0) == 0.0


def test_infeasible_is_not_zero_growth(toy):
    toy.reactions.r_1714.lower_bound = -1.0
    toy.reactions.r_2111.lower_bound = 2.0
    with pytest.warns(UserWarning, match="infeasible"):
        solution = toy.optimize()
    assert solution.status == "infeasible"
    with pytest.warns(UserWarning, match="infeasible"):
        value = gem_diff._growth(toy, -1.0)
    assert math.isnan(value)


def test_growth_restores_both_exchange_bounds(toy):
    before = {reaction.id: reaction.bounds for reaction in toy.reactions}
    assert gem_diff._growth(toy, -1.5) == pytest.approx(1.5)
    assert {reaction.id: reaction.bounds for reaction in toy.reactions} == before


@pytest.mark.parametrize("status", ["infeasible", "unbounded", "undefined"])
def test_structured_nonoptimal_status_is_checked_before_objective(toy, monkeypatch, status):
    class FailedSolution:
        @property
        def objective_value(self):
            raise AssertionError("a nonoptimal objective must not be inspected")

    solution = FailedSolution()
    solution.status = status
    monkeypatch.setattr(toy, "optimize", lambda: solution)
    assert math.isnan(gem_diff._growth(toy, -1.0))


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -float("inf")])
def test_optimal_missing_or_nonfinite_objective_is_not_growth(toy, monkeypatch, value):
    monkeypatch.setattr(toy, "optimize", lambda: SimpleNamespace(status="optimal", objective_value=value))
    assert math.isnan(gem_diff._growth(toy, -1.0))


@pytest.mark.parametrize("candidate", [0.0, 1.0, -1.0])
def test_zero_reference_percentage_is_explicitly_undefined(candidate):
    assert gem_diff._percent_change(0.0, candidate) == "undefined (zero reference)"


@pytest.mark.parametrize("reference,candidate", [(float("nan"), 1.0), (1.0, float("nan")),
                                               (float("inf"), 1.0), (1.0, -float("inf"))])
def test_nonfinite_percentage_is_explicitly_undefined(reference, candidate):
    assert gem_diff._percent_change(reference, candidate) == "undefined (non-finite growth)"


@pytest.mark.parametrize("reference,candidate,expected", [(2.0, 1.0, "-50.00%"),
                                                         (1.0, 2.0, "+100.00%"),
                                                         (1.0, 0.0, "-100.00%"),
                                                         (1.0, 1.0, "+0.00%")])
def test_positive_reference_percentage_is_unchanged(reference, candidate, expected):
    assert gem_diff._percent_change(reference, candidate) == expected


@pytest.mark.parametrize("reference,expected", [(0.0, "undefined (zero reference)"),
                                               (float("nan"), "undefined (non-finite growth)")])
def test_comparison_writes_undefined_instead_of_nan_or_zero_percent(toy, tmp_path, monkeypatch, reference, expected):
    candidate = tmp_path / "other-candidate.xml"
    candidate.write_text("test loader supplies a model")
    old, new = toy, toy.copy()
    monkeypatch.setattr(gem_diff.paths, "yeast_gem", lambda: tmp_path / "reference.xml")
    monkeypatch.setattr(gem_diff, "load_model", lambda path: (new if path == candidate else old, "test GLPK"))
    monkeypatch.setattr(gem_diff, "_growth", lambda model, bound: reference if model is old else 1.0)
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    assert gem_diff.main(["gem_version_diff.py", str(candidate)]) == 0
    frame = pd.read_csv(tmp_path / "gem_version_diff.csv", dtype=str, keep_default_na=False)
    growth = frame[frame["check"].str.startswith("growth /h")]
    assert len(growth) == 3
    assert growth["delta"].tolist() == [expected] * 3
    assert growth["vendored"].tolist() == (["0.0"] if reference == 0 else [""]) * 3


def _no_network(*args, **kwargs):
    raise AssertionError("ordinary verification/comparison must not invoke network acquisition")


def _snapshot(files):
    return {str(path): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
            for path in files}


def _write_entry(root, entry):
    (root / "data/public_inputs.json").write_text(json.dumps({"files": [entry]}), encoding="utf-8")


@pytest.fixture
def candidate_checkout(tmp_path):
    root = tmp_path.resolve() / "checkout"
    (root / "data/gem").mkdir(parents=True)
    entry = deepcopy(gem_diff.candidate_entry(ROOT))
    payload = b"<sbml>synthetic comparison asset</sbml>\n"
    compressed = gem_diff._compress_candidate(payload)
    entry.update(sha256=hashlib.sha256(compressed).hexdigest(), size_bytes=len(compressed),
                 uncompressed_sha256=hashlib.sha256(payload).hexdigest(), uncompressed_size_bytes=len(payload))
    entry["upstream"]["git_blob"] = hashlib.sha1(
        f"blob {len(payload)}\0".encode() + payload, usedforsecurity=False).hexdigest()
    _write_entry(root, entry)
    for name in ("yeast-GEM.xml.gz", "ecYeastGEM_batch.xml.gz", "ecYeastGEM_yeast902.xml.gz"):
        (root / "data/gem" / name).write_bytes(b"protected model fixture")
    return root, entry, payload, compressed


def _github_responses(monkeypatch, entry, payload, *, metadata_changes=None, blob_changes=None):
    upstream = entry["upstream"]
    metadata = {"type": "file", "path": upstream["path"], "sha": upstream["git_blob"],
                "size": entry["uncompressed_size_bytes"]}
    metadata.update(metadata_changes or {})
    blob = {"sha": upstream["git_blob"], "size": entry["uncompressed_size_bytes"],
            "encoding": "base64", "content": base64.b64encode(payload).decode()}
    blob.update(blob_changes or {})
    endpoints = [f"repos/{upstream['repository']}/contents/{upstream['path']}?ref={upstream['commit']}",
                 f"repos/{upstream['repository']}/git/blobs/{upstream['git_blob']}"]
    calls = []

    def run(argv, *, capture_output, check, timeout):
        assert capture_output is True and check is True and timeout == 120
        assert argv == ["gh", "api", "--method", "GET", endpoints[len(calls)]]
        calls.append(argv)
        return SimpleNamespace(stdout=json.dumps(metadata if len(calls) == 1 else blob).encode())

    monkeypatch.setattr(gem_diff.subprocess, "run", run)
    return calls


def test_explicit_fetch_creates_only_verified_deterministic_candidate(candidate_checkout, monkeypatch):
    root, entry, payload, compressed = candidate_checkout
    protected = [path for path in root.rglob("*") if path.is_file()]
    before = _snapshot(protected)
    calls = _github_responses(monkeypatch, entry, payload)
    path, actual_entry = gem_diff.provision_candidate(root=root, fetch=True)
    assert path == root / gem_diff.CANDIDATE
    assert actual_entry == entry
    assert path.read_bytes() == compressed
    assert len(calls) == 2
    assert _snapshot(protected) == before
    assert set(path for path in root.rglob("*") if path.is_file()) == set(protected) | {path}
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    created = _snapshot([path, *protected])
    gem_diff.provision_candidate(root=root)
    gem_diff.provision_candidate(root=root, fetch=True)
    assert _snapshot([path, *protected]) == created


def test_missing_candidate_verification_does_not_fetch_or_write(candidate_checkout, monkeypatch):
    root, _, _, _ = candidate_checkout
    before = _snapshot(path for path in root.rglob("*") if path.is_file())
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    with pytest.raises(FileNotFoundError):
        gem_diff.provision_candidate(root=root)
    assert _snapshot(path for path in root.rglob("*") if path.is_file()) == before


@pytest.mark.parametrize("fetch", [False, True])
def test_existing_bad_candidate_is_never_overwritten_or_refetched(candidate_checkout, monkeypatch, fetch):
    root, _, _, _ = candidate_checkout
    path = root / gem_diff.CANDIDATE
    path.write_bytes(b"not the pinned candidate")
    before = _snapshot(path for path in root.rglob("*") if path.is_file())
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    with pytest.raises(ValueError, match="mismatch"):
        gem_diff.provision_candidate(root=root, fetch=fetch)
    assert _snapshot(path for path in root.rglob("*") if path.is_file()) == before


@pytest.mark.parametrize("metadata_changes,blob_changes", [
    ({"type": "symlink"}, {}), ({"path": "model/other.xml"}, {}),
    ({"sha": "0" * 40}, {}), ({"size": 1}, {}),
    ({}, {"sha": "0" * 40}), ({}, {"size": 1}), ({}, {"encoding": "none"}),
    ({}, {"content": "not valid base64!"}),
])
def test_fetch_rejects_wrong_commit_path_blob_metadata_before_creating_file(
        candidate_checkout, monkeypatch, metadata_changes, blob_changes):
    root, entry, payload, _ = candidate_checkout
    _github_responses(monkeypatch, entry, payload, metadata_changes=metadata_changes, blob_changes=blob_changes)
    with pytest.raises(ValueError):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert not (root / gem_diff.CANDIDATE).exists()


def test_fetch_rejects_changed_sbml_even_with_matching_api_metadata(candidate_checkout, monkeypatch):
    root, entry, payload, _ = candidate_checkout
    _github_responses(monkeypatch, entry, b"X" + payload[1:])
    with pytest.raises(ValueError, match="SBML SHA-256 mismatch"):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert not (root / gem_diff.CANDIDATE).exists()


def test_fetch_does_not_publish_wrong_compression_bytes(candidate_checkout, monkeypatch):
    root, entry, payload, compressed = candidate_checkout
    _github_responses(monkeypatch, entry, payload)
    monkeypatch.setattr(gem_diff, "_compress_candidate", lambda payload: b"X" + compressed[1:])
    with pytest.raises(ValueError, match="compressed candidate SHA-256 mismatch"):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert not (root / gem_diff.CANDIDATE).exists()


def test_failed_github_request_does_not_retry_or_publish(candidate_checkout, monkeypatch):
    root, _, _, _ = candidate_checkout
    calls = []

    def unavailable(argv, **kwargs):
        calls.append(argv)
        raise subprocess.CalledProcessError(1, argv, stderr=b"unavailable")

    monkeypatch.setattr(gem_diff.subprocess, "run", unavailable)
    with pytest.raises(subprocess.CalledProcessError):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert len(calls) == 1
    assert not (root / gem_diff.CANDIDATE).exists()


@pytest.mark.parametrize("key,value,match", [
    ("sha256", "0" * 64, "compressed candidate SHA-256"),
    ("size_bytes", 1, "compressed candidate byte count"),
    ("uncompressed_sha256", "0" * 64, "candidate SBML SHA-256"),
    ("uncompressed_size_bytes", 1, "candidate SBML byte count"),
    ("git_blob", "0" * 40, "candidate SBML Git blob"),
])
def test_verification_rechecks_each_recorded_byte_identity(candidate_checkout, key, value, match):
    root, entry, _, compressed = candidate_checkout
    path = root / gem_diff.CANDIDATE
    path.write_bytes(compressed)
    (entry["upstream"] if key == "git_blob" else entry)[key] = value
    with pytest.raises(ValueError, match=match):
        gem_diff.verify_candidate(path, entry)


@pytest.mark.parametrize("mutation", ["duplicate", "source", "compression", "commit", "url", "size"])
def test_invalid_candidate_manifest_cannot_authorize_acquisition(candidate_checkout, monkeypatch, mutation):
    root, entry, _, _ = candidate_checkout
    if mutation == "source":
        entry["source"] = "download"
    elif mutation == "compression":
        entry["compression"]["mtime"] = 1
    elif mutation == "commit":
        entry["upstream"]["commit"] = "v9.1.1"
    elif mutation == "url":
        entry["upstream"]["url"] = "https://example.invalid/model.xml"
    elif mutation == "size":
        entry["size_bytes"] = 0
    _write_entry(root, entry)
    if mutation == "duplicate":
        (root / "data/public_inputs.json").write_text(json.dumps({"files": [entry, entry]}))
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    with pytest.raises(ValueError):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert not (root / gem_diff.CANDIDATE).exists()


def test_fetch_refuses_symlink_destination(candidate_checkout, monkeypatch, tmp_path):
    root, _, _, _ = candidate_checkout
    outside = tmp_path / "outside.xml.gz"
    (root / gem_diff.CANDIDATE).symlink_to(outside)
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    with pytest.raises(ValueError, match="symlink"):
        gem_diff.provision_candidate(root=root, fetch=True)
    assert not outside.exists()


@pytest.mark.parametrize("mode", ["--verify-candidate", "--fetch-candidate"])
def test_candidate_cli_modes_do_not_load_models_or_write_tables(candidate_checkout, monkeypatch, mode):
    root, _, _, compressed = candidate_checkout
    (root / gem_diff.CANDIDATE).write_bytes(compressed)
    monkeypatch.setattr(gem_diff.paths, "REPO_ROOT", root)
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    monkeypatch.setattr(gem_diff, "load_model", _no_network)
    monkeypatch.setattr(gem_diff.paths, "outputs_dir", _no_network)
    before = _snapshot(path for path in root.rglob("*") if path.is_file())
    assert gem_diff.main(["gem_version_diff.py", mode]) == 0
    assert _snapshot(path for path in root.rglob("*") if path.is_file()) == before


@pytest.mark.parametrize("mode", ["--verify-candidate", "--fetch-candidate"])
def test_candidate_cli_modes_cannot_target_existing_reference(mode, monkeypatch):
    monkeypatch.setattr(gem_diff, "provision_candidate", _no_network)
    with pytest.raises(SystemExit) as error:
        gem_diff.main(["gem_version_diff.py", mode, str(ROOT / "data/gem/yeast-GEM.xml.gz")])
    assert error.value.code == 2


def test_comparison_rejects_changed_registered_candidate_before_loading(candidate_checkout, monkeypatch):
    root, _, _, _ = candidate_checkout
    path = root / gem_diff.CANDIDATE
    path.write_bytes(b"corrupt candidate")
    monkeypatch.setattr(gem_diff.paths, "REPO_ROOT", root)
    monkeypatch.setattr(gem_diff, "load_model", _no_network)
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    monkeypatch.setattr(gem_diff.paths, "outputs_dir", _no_network)
    assert gem_diff.main(["gem_version_diff.py", str(path)]) == 1


def test_missing_registered_candidate_cli_fails_offline(candidate_checkout, monkeypatch):
    root, _, _, _ = candidate_checkout
    monkeypatch.setattr(gem_diff.paths, "REPO_ROOT", root)
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    assert gem_diff.main(["gem_version_diff.py", "--verify-candidate"]) == 1
    assert gem_diff.main(["gem_version_diff.py", str(root / gem_diff.CANDIDATE)]) == 1
    assert not (root / gem_diff.CANDIDATE).exists()


def test_real_candidate_matches_upstream_identity_and_deterministic_gzip(monkeypatch):
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    path, entry = gem_diff.provision_candidate()
    compressed = path.read_bytes()
    payload = gzip.decompress(compressed)
    assert entry["upstream"]["commit"] == "2d594ae1c4a2d550ccef120d96a58c7bbf586255"
    assert entry["upstream"]["git_blob"] == "04531c1e56dacaa2874c2a62bb48e67a1e3c2fd3"
    assert len(payload) == 11568534
    assert hashlib.sha256(payload).hexdigest() == "30842b15eefb0ef7e36cbdea86a9efddfacf69a871c8b054165faa9af9f6c8eb"
    assert len(compressed) == 332223
    assert hashlib.sha256(compressed).hexdigest() == "d1c00f217c9f654d884d8f9740e050eeb5d2193cb6e3eed71210deb65fb89e14"
    assert compressed[:10].hex() == "1f8b08000000000002ff"
    assert gem_diff._compress_candidate(payload) == gem_diff._compress_candidate(payload) == compressed
    assert entry["license"]["spdx"] == "CC-BY-4.0"
    assert entry["license"]["git_blob"] == "d7244c0d623e30a760d098e4ed829f4fc37ac262"


def test_public_provisioner_only_verifies_candidate_and_leaves_default_resolvers_unchanged(monkeypatch):
    spec = importlib.util.spec_from_file_location("_gem_public_inputs", ROOT / "scripts/provision_public_inputs.py")
    inputs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inputs)
    document = inputs.load_manifest(ROOT)
    entry = next(entry for entry in document["files"] if entry["path"] == gem_diff.CANDIDATE.as_posix())
    assert inputs.source_url(document["sources"][entry["source"]], entry) is None
    before = _snapshot([ROOT / gem_diff.CANDIDATE])
    for _ in range(2):
        assert inputs.verify_file(ROOT / gem_diff.CANDIDATE, entry) == 332223
    assert _snapshot([ROOT / gem_diff.CANDIDATE]) == before
    monkeypatch.delenv("YSTWIN_YEAST_GEM", raising=False)
    monkeypatch.delenv("YSTWIN_EC_YEAST_GEM", raising=False)
    assert paths.yeast_gem() == ROOT / "data/gem/yeast-GEM.xml.gz"
    assert paths.ec_yeast_gem() == ROOT / "data/gem/ecYeastGEM_batch.xml.gz"
    assert document["environment"]["YSTWIN_YEAST_GEM"] == "data/gem/yeast-GEM.xml.gz"
    assert document["environment"]["YSTWIN_EC_YEAST_GEM"] == "data/gem/ecYeastGEM_batch.xml.gz"


@pytest.mark.integration
def test_real_glpk_comparison_is_byte_identical_without_changing_models_or_retained_csv(tmp_path, monkeypatch):
    protected = [ROOT / gem_diff.CANDIDATE, ROOT / "data/gem/yeast-GEM.xml.gz",
                 ROOT / "data/gem/ecYeastGEM_batch.xml.gz", ROOT / "data/gem/ecYeastGEM_yeast902.xml.gz",
                 ROOT / "data/public_inputs.json", ROOT / "outputs/gem_version_diff.csv"]
    before = _snapshot(protected)
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    monkeypatch.setenv("YSTWIN_YEAST_GEM", str(ROOT / "data/gem/yeast-GEM.xml.gz"))
    monkeypatch.setattr(cobra.Configuration(), "solver", "glpk_exact")
    monkeypatch.setattr(gem_diff.subprocess, "run", _no_network)
    loaded = []

    def pinned_load(path):
        model, settings = load_model(path)
        assert type(model.solver).__module__ == "optlang.glpk_interface"
        assert settings.solver == "glpk" and settings.tolerance == 1e-7
        assert model.tolerance == pytest.approx(1e-7)
        assert model.objective.direction == "max"
        assert cobra.util.solver.linear_reaction_coefficients(model) == {model.reactions.r_2111: 1.0}
        loaded.append(path)
        return model, settings

    monkeypatch.setattr(gem_diff, "load_model", pinned_load)
    assert gem_diff.main(["gem_version_diff.py", str(ROOT / gem_diff.CANDIDATE)]) == 0
    assert loaded == [ROOT / "data/gem/yeast-GEM.xml.gz", ROOT / gem_diff.CANDIDATE]
    assert (tmp_path / "gem_version_diff.csv").read_bytes() == (ROOT / "outputs/gem_version_diff.csv").read_bytes()
    assert _snapshot(protected) == before
