from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(params=["audit_determinism", "audit_output_tables"])
def audit(request, monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location(request.param, ROOT / "scripts" / f"{request.param}.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, request.param, module)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "REPO", tmp_path)
    if request.param == "audit_determinism":
        monkeypatch.setattr(module, "collect", lambda **kwargs: [
            {"check": "example", "status": "PASS", "expected": "same", "actual": "same", "detail": ""}])
    else:
        monkeypatch.setattr(module, "REGISTRY", {"example": ("example.csv",)})
        monkeypatch.setattr(module, "EXEMPT", {})
        # DRIFT rows fail unconditionally, so an unpatched one leaks in and flips main() to 1,
        # stopping these tests at the exit-code precondition before their real assertions run.
        monkeypatch.setattr(module, "DRIFT", {})
        monkeypatch.setattr(module, "check_pair", lambda *args: [
            {"script": "example", "table": "example.csv", "ok": True, "detail": "reproduces"}])
        monkeypatch.setattr(module, "check_coverage", lambda: [])
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path / "outputs"))
    return request.param, module


def test_verification_does_not_mutate_retained_output(audit, monkeypatch, tmp_path):
    name, module = audit
    destination = tmp_path / "outputs"
    destination.mkdir()
    report = destination / f"{name}.csv"
    report.write_bytes(b"retained evidence\n")
    monkeypatch.setattr(sys, "argv", [name, "--quiet"])
    assert module.main() == 0
    assert report.read_bytes() == b"retained evidence\n"


def test_explicit_report_destination_is_honoured(audit, monkeypatch, tmp_path):
    name, module = audit
    destination = tmp_path / "diagnostics"
    monkeypatch.setattr(sys, "argv", [name, "--quiet", "--output-dir", str(destination)])
    assert module.main() == 0
    assert (destination / f"{name}.csv").is_file()
    assert not (tmp_path / "outputs").exists()
