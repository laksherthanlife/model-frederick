"""Run the frozen phase C parser test against the registered release source.

``data/native_law_v2/granados/development/phase_c_01/test_phase_c_parser.py`` cannot be
edited and cannot be collected where it sits. It is byte-frozen custody evidence pinned to
``CODE_REF`` by :mod:`tests.test_population_historical_sources`, and it loads its subject
with ``runpy.run_path`` from ``release_phase_c.py`` -- a path an approved storage
relocation moved to ``release_phase_c.py.source``, so collecting it in place raises
``FileNotFoundError`` before any assertion runs. ``pyproject.toml`` sets
``testpaths = ["tests"]``, which hides that rather than fixing it.

The relocation record answers this itself: the registered asset is verified and then
materialized at the *original* logical path. Doing that in ``tmp_path`` gives the frozen
test the tree it was written against, so it runs unmodified, from bytes checked against
the digest the relocation was approved with. Nothing is written into the repository and
the frozen file keeps its pinned bytes.

The byte identity of everything under ``.../granados/development`` is
:mod:`tests.test_population_historical_sources`'s job, not this module's. Here the only
verified artifact is the relocated release source, because that is the one this test has
to reconstruct a name for.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import unittest

import pytest


ROOT = Path(__file__).resolve().parents[1]
DELETION_MANIFEST = "data/deletion_manifest.json"
PHASE_C_DIR = "data/native_law_v2/granados/development/phase_c_01"
PARSER_TEST_PATH = f"{PHASE_C_DIR}/test_phase_c_parser.py"
HELPER_PATH = "data/native_law_v2/granados/development/executor_refresh_01/issue_refresh.py"
RELEASE_PATH = f"{PHASE_C_DIR}/release_phase_c.py"


def registered_relocation():
    """The approved storage relocation for the phase C release script.

    Read from the deletion manifest rather than retyped, so the digest checked here
    is the one the relocation was reviewed and approved against.
    """
    manifest = json.loads((ROOT / DELETION_MANIFEST).read_text(encoding="utf-8"))
    record = manifest["review_evidence"]["historical_source"]
    assert record["logical_source_path"] == RELEASE_PATH, (
        f"registered relocation does not describe {RELEASE_PATH}")
    return record


def supplied_release_source(record):
    """The release bytes, from whichever of the two registered names is present.

    Both present is ambiguity and neither is loss; either way there is nothing to
    verify against, and Git must not be asked to supply the difference.
    """
    candidates = [ROOT / record["logical_source_path"], ROOT / record["storage_path"]]
    supplied = [path for path in candidates if path.exists()]
    assert len(supplied) == 1, (
        f"release source requires exactly one original or registered asset: {RELEASE_PATH}")
    content = supplied[0].read_bytes()
    assert len(content) == record["size_bytes"], (
        f"registered release byte inventory mismatch: {supplied[0]}")
    assert hashlib.sha256(content).hexdigest() == record["sha256"], (
        f"registered release digest mismatch: {supplied[0]}")
    return content


def materialize_original_phase_c(directory):
    """Rebuild the tree the frozen parser test was written against, under ``directory``.

    Only what the test reaches at import: itself, the release script at its original
    logical name, and the refresh helper that ``release_phase_c.py`` runs on load.
    """
    directory = Path(directory).resolve(strict=True)
    for relative in (PARSER_TEST_PATH, HELPER_PATH):
        destination = directory / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, destination)
    release = directory / RELEASE_PATH
    release.parent.mkdir(parents=True, exist_ok=True)
    release.write_bytes(supplied_release_source(registered_relocation()))
    return directory / PARSER_TEST_PATH


def load_parser_test_module(path):
    spec = importlib.util.spec_from_file_location("materialized_phase_c_parser_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def materialized_parser_test(tmp_path_factory):
    root = tmp_path_factory.mktemp("original_phase_c_release").resolve(strict=True)
    return load_parser_test_module(materialize_original_phase_c(root))


def test_relocated_release_source_is_the_registered_bytes():
    record = registered_relocation()
    assert record["storage_path"] == RELEASE_PATH + ".source"
    assert not (ROOT / RELEASE_PATH).exists(), (
        "the original logical path must stay absent; both names present is ambiguity")
    assert hashlib.sha256(supplied_release_source(record)).hexdigest() == record["sha256"]


def test_frozen_parser_test_cannot_load_its_subject_in_place():
    """The defect this module exists for, asserted rather than described.

    If the release script ever returns to its logical path this fails, and the
    frozen test can then be collected where it sits instead.
    """
    assert (ROOT / PARSER_TEST_PATH).is_file()
    assert not (ROOT / RELEASE_PATH).exists()


def test_materialized_parser_test_loads_the_verified_release(materialized_parser_test):
    module = materialized_parser_test
    loaded = Path(module.MODULE["__file__"]).resolve()
    assert loaded.name == Path(RELEASE_PATH).name
    assert not loaded.is_relative_to(ROOT), "the release must be loaded from the sandbox"
    assert module.Reader.__name__ == "SelectedResponseReader"


def test_frozen_parser_assertions_pass_unmodified(materialized_parser_test):
    suite = unittest.defaultTestLoader.loadTestsFromModule(materialized_parser_test)
    assert suite.countTestCases() == 5
    result = unittest.TextTestRunner(verbosity=0).run(suite)
    assert not result.failures and not result.errors, (
        [str(entry) for entry in result.failures + result.errors])
