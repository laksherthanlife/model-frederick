"""Exercise the main claim gate against deliberately wrong synthetic repositories.

Repository-reference checks remain separate from scientific support. The driver must
consume the declared inventory once, preserve historical skips and verified refusals,
and block current claims with missing evidence or provenance. Unmarked prose is not an
inventory, regardless of how many numbers it contains.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
import pathlib
import subprocess
import sys

import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "audit_claims.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("audit_claims", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _registry(repo):
    return json.loads((repo / "data/current_claims.json").read_text())


def _save(repo, registry):
    (repo / "data/current_claims.json").write_text(json.dumps(registry))


@pytest.fixture
def repo(tmp_path, script, monkeypatch):
    """A minimal repository satisfying the reference and declared-evidence checks."""
    for relative in ("src/ystwin", "docs", "scripts", "tests", "outputs", "data"):
        (tmp_path / relative).mkdir(parents=True)
    (tmp_path / "src/ystwin/reporter.py").write_text("VALUE = 1\n")
    (tmp_path / "scripts/run_thing.py").write_text("VALUE = 1\n")
    (tmp_path / "tests/test_thing.py").write_text("def test_it():\n    assert True\n")
    (tmp_path / "outputs/table.csv").write_text("dose_mM,skill\n0.2,0.368\n")
    (tmp_path / "outputs/empty.csv").write_text("")
    (tmp_path / "pyproject.toml").write_text('dependencies = [\n  "numpy==2.4.1",\n]\n')
    (tmp_path / "README.md").write_text("The reporter lives in `reporter.py`.\n")
    surfaces = [
        {"id": "readme", "glob": "README.md", "kind": "markdown", "role": "current"},
        {"id": "docs", "glob": "docs/**/*.md", "exclude": ["docs/superseded/**"],
         "kind": "markdown", "role": "current"},
        {"id": "historical", "glob": "docs/superseded/**/*.md", "kind": "markdown", "role": "historical"},
        {"id": "library", "glob": "src/**/*.py", "kind": "python_docstrings", "role": "current"},
        {"id": "fields", "glob": "src/**/*.py", "kind": "python_field_comments", "role": "current"},
        {"id": "scripts", "glob": "scripts/**/*.py", "kind": "python_docstrings", "role": "current"},
        {"id": "tests", "glob": "tests/**/*.py", "kind": "python_docstrings", "role": "forbidden"},
    ]
    for surface in surfaces:
        surface.update(unmarked_policy="registered_only", reason="Explicit fixture surface")
    registry = {
        "schema_version": 1, "surfaces": surfaces,
        "families": {"fixture": {"boundary": "A synthetic table, not biological evidence",
                                   "locations": [{"scope": "readme", "path": "README.md", "anchor": "reporter"}]}},
        "profiles": {"fixture": {
            "family": "fixture", "method": {"description": "Single fixture cell", "implementation": ["scripts/run_thing.py"]},
            "conditions": {"dose_mM": 0.2}, "model": "fixture model",
            "configuration": {"id": "fixture", "binding": {"status": "verified", "checks": [
                {"field": "dose_mM", "value": 0.2}]}},
            "provenance": {"status": "verified", "scope": "reproduced_result", "artifacts": [
                {"path": path, "sha256": hashlib.sha256((tmp_path / path).read_bytes()).hexdigest()}
                for path in ("outputs/table.csv", "scripts/run_thing.py")]},
            "uncertainty": {"status": "not_applicable", "reason": "Fixture identity comparison"},
            "coverage": {"expected_rows": 1}, "independent_unit": {"name": "fixture row"},
            "evidence": [{"id": "main", "path": "outputs/table.csv", "format": "csv",
                          "row_keys": ["dose_mM"], "unit": "dimensionless"}],
        }},
        "claims": [{"id": "fixture.skill", "profile": "fixture", "statement": "The fixture value is 0.368",
                    "role": "current", "evaluation": {"kind": "numeric", "aggregation": "identity", "column": "skill"},
                    "expected": {"value": 0.368}}],
        "marked_values": {"baseline_occurrences": 0, "tables": {}},
    }
    _save(tmp_path, registry)
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.email=t@example.com", "-c", "user.name=t",
                    "commit", "-q", "--allow-empty", "-m", "initial"], cwd=tmp_path, check=True)
    monkeypatch.setattr(script, "REPO", tmp_path)
    monkeypatch.setattr(script, "SRC", tmp_path / "src/ystwin")
    monkeypatch.setattr(script, "README", tmp_path / "README.md")
    return tmp_path


def _statuses(rows) -> list[str]:
    return [row["status"] for row in rows]


def _marker(literal="0.368", extra=""):
    return f'<!-- audit:value table=outputs/table.csv column=skill row="dose_mM=0.2"{extra} -->{literal}'


def _declare_marker(repo):
    registry = _registry(repo)
    registry["marked_values"]["tables"] = {"outputs/table.csv": {"profile": "fixture", "selectors": [
        {"where": {"dose_mM": 0.2}, "columns": {"skill": "dimensionless"}}]}}
    _save(repo, registry)


def _reported(script, repo):
    code = script.main(["--quiet", "--output-dir", str(repo / "diagnostics")])
    return code, pd.read_csv(repo / "diagnostics/audit_claims.csv")


class TestTheGateFailsOnAReferenceThatDoesNotResolve:
    def test_a_readme_naming_a_module_that_does_not_exist_fails(self, script, repo):
        rows = script.check_module_references("the entry point is `src/ystwin/gone.py`")
        assert _statuses(rows) == ["FAIL"]

    def test_a_module_that_does_exist_passes(self, script, repo):
        rows = script.check_module_references("see `reporter.py` and `scripts/run_thing.py`")
        assert _statuses(rows) == ["PASS", "PASS"]

    def test_a_documented_command_naming_a_missing_script_fails(self, script, repo):
        rows = script.check_documented_commands(
            "run `python3 scripts/gone.py` then `python3 scripts/run_thing.py`")
        assert _statuses(rows) == ["FAIL", "PASS"]


class TestTheGateFailsOnATableNobodyCanOpen:
    def test_a_cited_table_that_is_absent_fails(self, script, repo):
        assert _statuses(script.check_output_tables("the numbers are in `outputs/absent.csv`")) == ["FAIL"]

    def test_a_cited_table_that_is_empty_fails(self, script, repo):
        rows = script.check_output_tables("the numbers are in `outputs/empty.csv`")
        assert (rows[0]["status"], rows[0]["actual"]) == ("FAIL", "empty")

    def test_a_cited_table_with_content_passes(self, script, repo):
        assert _statuses(script.check_output_tables("the numbers are in `outputs/table.csv`")) == ["PASS"]

    def test_a_table_cited_as_a_markdown_link_is_not_reported_missing(self, script, repo):
        rows = script.check_output_tables("see [the table](outputs/table.csv) for it")
        assert _statuses(rows) == ["PASS"]

    def test_a_filename_whose_own_parens_balance_keeps_them(self, script, repo):
        (repo / "outputs/export_(RAW).csv").write_text("a,b\n1,2\n")
        assert _statuses(script.check_output_tables("raw is in `outputs/export_(RAW).csv`")) == ["PASS"]

    def test_a_balanced_name_inside_a_markdown_link_survives_both_rules(self, script, repo):
        (repo / "outputs/export_(RAW).csv").write_text("a,b\n1,2\n")
        rows = script.check_output_tables("see [raw](outputs/export_(RAW).csv) for it")
        assert _statuses(rows) == ["PASS"]


class TestExistenceChecksStopAtTheRepositoryBoundary:
    @staticmethod
    def _doc(repo, text):
        (repo / "docs/research").mkdir(parents=True, exist_ok=True)
        (repo / "docs/research/survey.md").write_text(text)

    def test_a_path_belonging_to_another_codebase_is_not_flagged(self, script, repo):
        self._doc(repo, "we read `pytfa/thermo/metabolite.py` and `flapjack_api/analysis/analysis.py`.\n")
        assert _statuses(script.check_docs_module_references()) == ["PASS"]

    def test_a_missing_path_under_an_owned_prefix_is_still_flagged(self, script, repo):
        self._doc(repo, "the fit happens in `src/ystwin/gone.py`.\n")
        rows = script.check_docs_module_references()
        assert (rows[0]["status"], rows[0]["detail"]) == ("FAIL", "src/ystwin/gone.py")

    def test_a_bare_module_name_is_still_checked(self, script, repo):
        self._doc(repo, "described in `nowhere.py`.\n")
        assert _statuses(script.check_docs_module_references()) == ["FAIL"]


class TestANumberClaimedInTheReadmeIsCheckedNotAssumed:
    def test_a_marked_count_the_repository_contradicts_is_not_accepted(self, script, repo):
        rows = script.check_test_count("<!-- audit:test_count --> The suite is 300 tests.\n")
        assert _statuses(rows) == ["FAIL"]

    def test_a_readme_making_no_count_claim_is_recorded_as_making_none(self, script, repo):
        rows = script.check_test_count("no numbers here at all\n")
        assert (rows[0]["status"], rows[0]["actual"]) == ("PASS", "none found")

    def test_an_unmarked_number_in_the_prose_is_not_read_as_a_claim(self, script):
        rows = script.check_test_count("an earlier README claimed 149 tests when the suite had 1,204.\n")
        assert (rows[0]["status"], rows[0]["actual"]) == ("PASS", "none found")

    def test_one_valid_count_cannot_hide_an_empty_marker(self, script, repo):
        rows = script.check_test_count("<!-- audit:test_count --> 1 test\n<!-- audit:test_count --> empty\n")
        assert _statuses(rows) == ["FAIL"]

    def test_architecture_collection_failure_is_not_reported_as_success(self, script, repo, monkeypatch):
        (repo / "docs/ARCHITECTURE.md").write_text(
            "Scope: every module under `src/ystwin/` (1 code modules + 0 package `__init__.py`, ~1 LOC), "
            "every script under `scripts/` (1), and the contracts the 1-test suite\n")
        monkeypatch.setattr(script, "_collected_test_count", lambda: None)
        rows = script.check_architecture_scope()
        row = next(row for row in rows if row["check"] == "architecture scope: tests collected")
        assert row["blocking"] and row["status"] == "FAIL" and row["actual"] == "could not count"


class TestTheGateGates:
    def test_a_repository_its_readme_describes_correctly_exits_zero(self, script, repo):
        assert script.main(["--quiet"]) == 0

    def test_a_readme_naming_a_module_that_does_not_exist_exits_non_zero(self, script, repo):
        (repo / "README.md").write_text("the entry point is `src/ystwin/gone.py`\n")
        assert script.main(["--quiet"]) == 1

    def test_the_failing_check_is_named_in_the_table_it_writes(self, script, repo):
        (repo / "README.md").write_text("the entry point is `src/ystwin/gone.py`\n")
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame[frame.blocking].check.tolist() == ["module referenced: src/ystwin/gone.py"]

    def test_a_missing_declared_registry_is_blocking_not_optional(self, script, repo):
        (repo / "data/current_claims.json").unlink()
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame[frame.blocking].id.tolist() == ["registry.schema"]

    def test_historical_records_are_skips_not_failures(self, script, repo, capsys):
        registry = _registry(repo)
        old = deepcopy(registry["claims"][0])
        old.update(id="fixture.old", role="historical", expected={"value": 999})
        registry["claims"].append(old)
        _save(repo, registry)
        assert script.main([]) == 0
        output = capsys.readouterr().out
        assert "[skip] structured claim: fixture.old" in output
        assert "[FAIL]" not in output
        assert "1 SKIP" in output

    def test_correct_scientific_refusal_passes_and_is_not_a_positive_result(self, script, repo):
        registry = _registry(repo)
        claim = registry["claims"][0]
        claim.update(role="refused", evaluation={"kind": "refusal"}, acceptance={
            "predicates": [{"field": "skill", "operator": "lt", "value": 0.5}]})
        claim.pop("expected")
        _save(repo, registry)
        code, frame = _reported(script, repo)
        assert code == 0
        row = frame[frame.id == "fixture.skill"].iloc[0]
        assert row.status == "PASS" and row.verdict == "REFUSAL_CONFIRMED" and not row.blocking

    def test_matching_numbers_do_not_rescue_unknown_provenance(self, script, repo):
        registry = _registry(repo)
        registry["profiles"]["fixture"]["provenance"] = {"status": "pending", "reason": "No verified receipt"}
        _save(repo, registry)
        code, frame = _reported(script, repo)
        assert code == 1
        row = frame[frame.id == "fixture.skill"].iloc[0]
        assert row.comparison_status == "PASS" and row.verdict == "BLOCKED"
        assert pd.isna(row.actual) and row.observed == pytest.approx(0.368)

    def test_exit_code_uses_blocking_contract_not_status_spelling(self, script, repo, monkeypatch):
        monkeypatch.setattr(script, "audit_claims", lambda root: [{
            "id": "fixture.blocking", "check": "fixture blocker", "record_type": "claim", "status": "SKIP",
            "blocking": True, "expected": "known evidence", "actual": None, "detail": "unresolved"}])
        assert script.main(["--quiet"]) == 1


class TestDeclaredMarkersReachTheDriverOnce:
    def test_inventory_and_marker_evaluation_are_called_once(self, script, repo, monkeypatch):
        _declare_marker(repo)
        (repo / "README.md").write_text(_marker() + " and " + _marker() + "\n")
        calls = []
        delegate = script.audit_claims

        def spy(root):
            calls.append(root)
            return delegate(root)

        monkeypatch.setattr(script, "audit_claims", spy)
        code, frame = _reported(script, repo)
        assert code == 0 and calls == [repo]
        assert len(frame[frame.record_type == "marked_value"]) == 2
        assert frame.id.dropna().is_unique

    @pytest.mark.parametrize("literal,code", [("0.368", 0), ("0.37", 0), ("0.38", 1), ("999.0", 1)])
    def test_explicit_cell_and_its_rounding_control_the_result(self, script, repo, literal, code):
        _declare_marker(repo)
        (repo / "README.md").write_text(_marker(literal) + "\n")
        assert script.main(["--quiet"]) == code

    def test_a_new_undeclared_selector_fails(self, script, repo):
        (repo / "README.md").write_text(_marker())
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame.detail.str.contains("unregistered authoritative selector", na=False).any()

    def test_removing_a_declared_marker_cannot_shrink_coverage_silently(self, script, repo):
        _declare_marker(repo)
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame.detail.str.contains("no remaining surface marker", na=False).any()

    def test_historical_marker_disagreement_is_nonblocking(self, script, repo, capsys):
        (repo / "docs/superseded").mkdir()
        (repo / "docs/superseded/old.md").write_text(_marker("999"))
        assert script.main([]) == 0
        output = capsys.readouterr().out
        assert "[skip]" in output and "[FAIL]" not in output

    def test_a_historical_alias_cannot_replace_current_marker_coverage(self, script, repo):
        _declare_marker(repo)
        (repo / "docs/superseded").mkdir()
        (repo / "docs/superseded/old.md").write_text(_marker())
        assert script.main(["--quiet"]) == 1

    @pytest.mark.parametrize("relative,body", [
        ("docs/current.md", "{marker}"),
        ("src/ystwin/claims.py", '"""{marker}"""\n'),
        ("src/ystwin/claims.py", "#: {marker}\nVALUE = 1\n"),
        ("scripts/claims.py", '"""{marker}"""\n'),
    ])
    def test_declared_documentation_surfaces_all_check_their_markers(self, script, repo, relative, body):
        _declare_marker(repo)
        path = repo / relative
        path.write_text(body.format(marker=_marker("999")))
        assert script.main(["--quiet"]) == 1
        path.write_text(body.format(marker=_marker()))
        assert script.main(["--quiet"]) == 0

    def test_executable_fixture_literals_are_not_docstring_claims(self, script, repo):
        (repo / "tests/fixture.py").write_text(f"FIXTURE = {_marker('999')!r}\n")
        assert script.main(["--quiet"]) == 0
        (repo / "tests/fixture.py").write_text(f"'''{_marker('999')}'''\n")
        assert script.main(["--quiet"]) == 1

    def test_unmarked_decimals_do_not_create_or_exempt_current_claims(self, script, repo):
        prose = "In 2025 there were 1,204 tests; see section 6.5, 2.718x and 95%. " * 1000
        (repo / "README.md").write_text(prose)
        (repo / "docs/current.md").write_text(prose)
        (repo / "src/ystwin/claims.py").write_text(f'"""{prose}"""\n#: {prose}\nVALUE = 1\n')
        (repo / "scripts/claims.py").write_text(f'"""{prose}"""\n')
        code, frame = _reported(script, repo)
        assert code == 0
        assert frame[frame.record_type == "claim"].id.tolist() == ["fixture.skill"]
        _declare_marker(repo)
        (repo / "README.md").write_text(prose + "\n" + _marker("999"))
        assert script.main(["--quiet"]) == 1


class TestReadOnlyVerification:
    def test_default_audit_does_not_create_or_modify_reports(self, script, repo, monkeypatch):
        report = repo / "outputs/audit_claims.csv"
        report.write_text("retained,report\noriginal,evidence\n")
        before = {path.relative_to(repo): path.read_bytes() for path in repo.rglob("*") if path.is_file()}
        monkeypatch.setattr(sys, "argv", ["audit_claims.py", "--quiet"])
        assert script.main() == 0
        after = {path.relative_to(repo): path.read_bytes() for path in repo.rglob("*") if path.is_file()}
        assert after == before

    def test_explicit_output_directory_gets_the_report(self, script, repo, monkeypatch):
        destination = repo / "diagnostics"
        monkeypatch.setattr(sys, "argv", ["audit_claims.py", "--quiet", "--output-dir", str(destination)])
        assert script.main() == 0
        assert (destination / "audit_claims.csv").is_file()
        assert not (repo / "outputs/audit_claims.csv").exists()

    def test_failed_collection_cannot_supply_a_successful_count(self, script, monkeypatch):
        result = subprocess.CompletedProcess(["pytest"], 2, "5 tests collected\n", "error")
        monkeypatch.setattr(script.subprocess, "run", lambda *args, **kwargs: result)
        assert script._collected_test_count() is None
        assert _statuses(script.check_test_count("<!-- audit:test_count --> 5 tests\n")) == ["FAIL"]

    @pytest.mark.parametrize("scale", ["nan", "inf", "-inf"])
    def test_nonfinite_marker_scales_are_refused(self, script, repo, scale):
        _declare_marker(repo)
        (repo / "README.md").write_text(_marker(extra=f" scale={scale}"))
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame.detail.str.contains("nonfinite", na=False).any()

    def test_duplicate_marker_keys_are_refused(self, script, repo):
        _declare_marker(repo)
        (repo / "README.md").write_text(_marker(extra=" table=outputs/other.csv"))
        code, frame = _reported(script, repo)
        assert code == 1
        assert frame.detail.str.contains("duplicate marker field", na=False).any()
