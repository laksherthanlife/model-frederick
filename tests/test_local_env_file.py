"""The untracked file that tells a machine where its own data is.

`ystwin.paths` deliberately has no machine-local defaults: a path found by convention makes
a result reproducible on one machine only, and both audit gates refuse one. The cost was
that a machine holding every wet-lab file still ran the suite blind, because eleven separate
variables had to be exported by hand and nothing said so -- fourteen tests skipped on the
machine that had all the data they needed.

`.ystwin.env` fixes that without reintroducing the defect: it is untracked, it is read only
by the test suite, and an exported variable always beats it.
"""

from __future__ import annotations

import os
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture
def local_loader(tmp_path, monkeypatch):
    """Exercise local loading with an isolated environment, never the checkout's file.

    Only these unit tests simulate a non-CI machine. The real suite keeps both its
    CI guard and its public-input preflight, and variables loaded here cannot leak.
    """
    from tests import conftest

    monkeypatch.setattr(conftest, "_ENV_FILE", tmp_path / ".ystwin.env")
    monkeypatch.setattr(os, "environ", os.environ.copy())
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.delenv("GITHUB_ACTIONS", raising=False)
    return conftest


class TestTheTemplateShipsAndTheRealFileDoesNot:
    def test_the_example_is_tracked(self):
        """A clone has to be able to see what the variables are called."""
        listed = subprocess.run(
            ["git", "ls-files", ".ystwin.env.example"],
            cwd=REPO, capture_output=True, text=True).stdout.strip()

        assert listed == ".ystwin.env.example"

    def test_the_real_file_is_ignored(self):
        """It names directories outside the checkout and, on the lab machine, a workbook
        whose document properties carry a private individual's name."""
        ignored = subprocess.run(
            ["git", "check-ignore", ".ystwin.env"],
            cwd=REPO, capture_output=True, text=True)

        assert ignored.returncode == 0

    def test_the_real_file_is_not_tracked(self):
        listed = subprocess.run(
            ["git", "ls-files", ".ystwin.env"],
            cwd=REPO, capture_output=True, text=True).stdout.strip()

        assert listed == ""

    def test_the_example_names_every_variable_the_loader_accepts(self):
        """A template missing a variable is how the next person loses a day."""
        text = (REPO / ".ystwin.env.example").read_text()

        for name in ("YSTWIN_GEN5_XPT", "YSTWIN_PLATES", "YSTWIN_IGEM_RESULTS",
                     "YSTWIN_QPCR_RAW", "YSTWIN_CROSSTALK_WORKBOOK"):
            assert name in text, name

    def test_the_example_carries_no_actual_path(self):
        """It is a template. A filled-in one would be the machine-local default again."""
        for line in (REPO / ".ystwin.env.example").read_text().splitlines():
            if line.startswith("YSTWIN_"):
                assert line.rstrip().endswith("="), line


class TestTheLoaderIsSafeWithoutTheFile:
    """A clone has no `.ystwin.env`, and the suite has to behave exactly as before."""

    def test_it_returns_quietly_when_the_file_is_absent(self, local_loader):
        assert not local_loader._ENV_FILE.exists()

        local_loader._load_local_env()  # must not raise

    def test_an_exported_variable_always_wins(self, local_loader, monkeypatch):
        """A one-off `YSTWIN_... = ... pytest` is unaffected by a stale local file."""
        local_loader._ENV_FILE.write_text("YSTWIN_TEST_ONLY=from_the_file\n")
        monkeypatch.setenv("YSTWIN_TEST_ONLY", "from_the_environment")

        local_loader._load_local_env()

        assert os.environ["YSTWIN_TEST_ONLY"] == "from_the_environment"

    def test_it_ignores_anything_not_a_ystwin_variable(self, local_loader, monkeypatch):
        """It is not a general dotenv loader and must not become a way to set PATH."""
        local_loader._ENV_FILE.write_text("PATH=/nowhere\nYSTWIN_TEST_OTHER=fine\n")
        monkeypatch.delenv("YSTWIN_TEST_OTHER", raising=False)
        original_path = os.environ.get("PATH")

        local_loader._load_local_env()

        assert os.environ.get("PATH") == original_path
        assert os.environ["YSTWIN_TEST_OTHER"] == "fine"

    def test_comments_and_blank_lines_are_skipped(self, local_loader, monkeypatch):
        local_loader._ENV_FILE.write_text("# a comment\n\nYSTWIN_TEST_THIRD=value\n")
        monkeypatch.delenv("YSTWIN_TEST_THIRD", raising=False)

        local_loader._load_local_env()

        assert os.environ["YSTWIN_TEST_THIRD"] == "value"
