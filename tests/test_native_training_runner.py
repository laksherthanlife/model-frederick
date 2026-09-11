from __future__ import annotations

import json
import os
import subprocess
import sys

from scripts import run_native_training as runner
from ystwin import paths


def test_training_contract_is_an_explicit_native_file_allowlist():
    protocol = json.loads((paths.REPO_ROOT / runner.PROTOCOL_PATH).read_text())
    contract = runner.build_contract(paths.REPO_ROOT, protocol)
    actual = {item.path for item in contract.inputs}
    expected = {entry["path"] for entry in protocol["training_inputs"]} | {runner.PROTOCOL_PATH}
    assert actual == expected
    assert all((contract.root / item.path).is_file() for item in (*contract.inputs, *contract.code))
    assert all(item.source["scope"] == "native_physiology" for item in contract.inputs)


def run_guard(root, allowed, output, operation):
    code = (
        "from pathlib import Path; from scripts.run_native_training import install_project_read_guard; "
        f"install_project_read_guard(Path({str(root)!r}), [Path({str(allowed)!r})], Path({str(output)!r})); "
        f"{operation}"
    )
    return subprocess.run([sys.executable, "-c", code], cwd=paths.REPO_ROOT, capture_output=True,
                          text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, timeout=30)


def test_training_read_guard_blocks_unapproved_outcomes_in_a_child_process(tmp_path):
    native = tmp_path / "native.json"
    withheld = tmp_path / "outcome.json"
    output = tmp_path / "output"
    native.write_text("native")
    withheld.write_text("not available to training")
    output.mkdir()
    accepted = run_guard(tmp_path, native, output, f"print(Path({str(native)!r}).read_text())")
    assert accepted.returncode == 0, accepted.stderr
    denied = run_guard(tmp_path, native, output, f"Path({str(withheld)!r}).read_text()")
    assert denied.returncode != 0
    assert "outside the explicit native allowlist" in denied.stderr


def test_training_guard_does_not_follow_a_project_symlink_to_an_unapproved_input(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    native = root / "native.json"
    output = root / "output"
    output.mkdir()
    native.write_text("native")
    outside = tmp_path / "outside.json"
    outside.write_text("not training data")
    alias = root / "alias.json"
    alias.symlink_to(outside)
    result = run_guard(root, native, output, f"Path({str(alias)!r}).read_text()")
    assert result.returncode != 0
    assert "outside the explicit native allowlist" in result.stderr


def test_training_guard_preserves_native_inputs(tmp_path):
    native = tmp_path / "native.json"
    native.write_text("unchanged")
    output = tmp_path / "output"
    output.mkdir()
    result = run_guard(tmp_path, native, output, f"Path({str(native)!r}).write_text('changed')")
    assert result.returncode != 0
    assert native.read_text() == "unchanged"
