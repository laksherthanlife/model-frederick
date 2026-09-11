from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import shutil

import pytest

from ystwin import artifacts
from ystwin.analysis.claims import ClaimRegistryError, audit_claims, load_registry, selector_identity


_ROOT = Path(__file__).resolve().parents[1]
_TABLE = (
    "key,value,reference,low,high,config,model,plate,weight\n"
    "a,1,1,0.9,1.1,v1,m1,p1,1\n"
    "b,2,1,1.7,2.3,v1,m1,p2,2\n"
    "c,99,1,95,103,v1,m1,p3,3\n"
)


def _pin(root, registry):
    for profile in registry["profiles"].values():
        if profile["provenance"].get("status") != "verified":
            continue
        paths = [source["path"] for source in profile["evidence"]]
        paths += profile["method"].get("implementation", [])
        profile["provenance"]["artifacts"] = [
            {"path": path, "sha256": hashlib.sha256((root / path).read_bytes()).hexdigest()}
            for path in dict.fromkeys(paths)
        ]


def _save(root, registry, *, pin=False, relative="data/current_claims.json"):
    if pin:
        _pin(root, registry)
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(registry), encoding="utf-8")
    return path


def _record(root, identity="fixture.median", registry_path=None):
    rows = audit_claims(root, registry_path)
    return next(row for row in rows if row["id"] == identity)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    for relative in ("data", "outputs", "src", "docs", "tests"):
        (root / relative).mkdir(parents=True)
    (root / "README.md").write_text("# Current\nA scientifically scoped statement.\n", encoding="utf-8")
    (root / "outputs/table.csv").write_text(_TABLE, encoding="utf-8")
    (root / "src/model.py").write_text("MODEL = 'm1'\n", encoding="utf-8")
    registry = {
        "schema_version": 1,
        "surfaces": [{"id": "readme", "glob": "README.md", "kind": "markdown", "role": "current",
                      "unmarked_policy": "registered_only", "required": True,
                      "reason": "Fixture current scientific surface"}],
        "families": {"fixture": {"boundary": "Only fixture measurements, not archival commentary",
                                  "locations": [{"scope": "readme", "path": "README.md", "anchor": "Current"}]}},
        "profiles": {"fixture": {
            "family": "fixture", "method": {"description": "Median over all independently identified rows",
                                               "implementation": ["src/model.py"]},
            "conditions": {"population": "all three fixture observations"}, "model": "m1",
            "configuration": {"id": "v1", "binding": {"status": "verified", "checks": [
                {"field": "config", "value": "v1"}, {"field": "model", "value": "m1"}]}},
            "provenance": {"status": "verified", "scope": "reproduced_result", "artifacts": []},
            "uncertainty": {"status": "not_applicable", "reason": "Fixture descriptive median"},
            "coverage": {"expected_rows": 3},
            "independent_unit": {"name": "plate", "columns": ["plate"], "expected_count": 3},
            "evidence": [{"id": "main", "path": "outputs/table.csv", "format": "csv", "row_keys": ["key"],
                          "unit": "dimensionless"}],
        }},
        "claims": [{"id": "fixture.median", "profile": "fixture", "statement": "The complete median is two",
                    "role": "current", "evaluation": {"kind": "numeric", "column": "value", "aggregation": "median"},
                    "expected": {"value": 2}}],
        "marked_values": {"baseline_occurrences": 0, "tables": {}},
    }
    _save(root, registry, pin=True)
    return root, registry


@pytest.fixture(scope="module")
def renderer():
    spec = importlib.util.spec_from_file_location("render_current_claims", _ROOT / "scripts/render_current_claims.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _add_marker(registry, where=None, column="value", scale=1):
    table = registry["marked_values"]["tables"].setdefault("outputs/table.csv", {"profile": "fixture", "selectors": []})
    table["selectors"].append({"where": where or {"key": "a"},
                              "columns": {column: {"unit": "dimensionless", "scale": scale}}})


def _marker(value="1.00", row="key=a", column="value", extra=""):
    return f'<!-- audit:value table=outputs/table.csv column={column} row="{row}"{extra} -->{value}'


def test_complete_registry_schema_is_expanded_and_audited_read_only(repo):
    root, registry = repo
    before = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
    expanded = load_registry(root)
    assert expanded["claims"][0]["method"] == registry["profiles"]["fixture"]["method"]
    row = _record(root)
    assert row["status"] == "PASS"
    assert row["actual"] == 2
    assert row["supported"]
    assert row["coverage_observed"] == {"selected_rows": 3, "independent_unit": "plate", "independent_units": 3}
    assert [result["value"] for result in row["result_rows"]] == [1, 2, 99]
    after = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
    assert before == after


def test_unmarked_decimals_years_and_historical_test_counts_do_not_create_claims(repo):
    root, _ = repo
    (root / "README.md").write_text("In 2025 the historical report said 1,204 tests and 2.718 error, see section 6.5.\n")
    rows = audit_claims(root)
    assert [row["id"] for row in rows if row["record_type"] == "claim"] == ["fixture.median"]
    assert not [row for row in rows if row["record_type"] == "marked_value"]
    assert not any(row["blocking"] for row in rows)


@pytest.mark.parametrize("field", ["method", "conditions", "model", "configuration", "provenance", "uncertainty", "coverage", "independent_unit"])
def test_required_scientific_metadata_cannot_be_dropped(repo, field):
    root, registry = repo
    del registry["profiles"]["fixture"][field]
    _save(root, registry)
    rows = audit_claims(root)
    assert rows[0]["id"] == "registry.schema"
    assert rows[0]["blocking"]


@pytest.mark.parametrize("change", ["no_surfaces", "duplicate_surface", "no_locations", "unknown_scope", "unscoped_family", "infer_decimals"])
def test_scope_is_an_explicit_validated_contract(repo, change):
    root, registry = repo
    if change == "no_surfaces":
        registry.pop("surfaces")
    elif change == "duplicate_surface":
        registry["surfaces"].append(deepcopy(registry["surfaces"][0]))
    elif change == "no_locations":
        registry["families"]["fixture"]["locations"] = []
    elif change == "unknown_scope":
        registry["families"]["fixture"]["locations"][0]["scope"] = "absent"
    elif change == "unscoped_family":
        registry["claims"][0]["family"] = "undeclared"
    else:
        registry["surfaces"][0]["unmarked_policy"] = "scan_decimals"
    _save(root, registry)
    assert audit_claims(root)[0]["id"] == "registry.schema"


def test_location_outside_declared_surface_blocks_audit(repo):
    root, registry = repo
    (root / "docs/elsewhere.md").write_text("outside the surface")
    registry["families"]["fixture"]["locations"][0]["path"] = "docs/elsewhere.md"
    _save(root, registry)
    assert any(row["blocking"] and row["id"].startswith("scope.family") for row in audit_claims(root))


def test_duplicate_stable_claim_ids_fail(repo):
    root, registry = repo
    registry["claims"].append(deepcopy(registry["claims"][0]))
    _save(root, registry)
    assert "duplicate claim id" in audit_claims(root)[0]["detail"]


def test_duplicate_json_keys_fail_instead_of_last_value_winning(repo):
    root, _ = repo
    path = root / "data/current_claims.json"
    path.write_text('{"schema_version": 1, "schema_version": 1}')
    assert "duplicate JSON key" in audit_claims(root)[0]["detail"]


def test_single_cell_selector_must_be_unique(repo):
    root, registry = repo
    registry["claims"][0]["evaluation"]["aggregation"] = "identity"
    _save(root, registry)
    row = _record(root)
    assert row["blocking"] and row["actual"] is None
    assert "matched 3 rows" in row["detail"]


@pytest.mark.parametrize("change, message", [
    ("duplicate_key", "duplicate row key"), ("blank_key", "missing row key"),
    ("duplicate_header", "duplicate columns"), ("no_match", "matched no rows"),
    ("missing_column", "missing column"), ("missing_independent_unit", "missing independent unit"),
    ("missing_row", "coverage mismatch"), ("bad_row", "malformed row"),
])
def test_tabular_missingness_and_identity_never_pass(repo, change, message):
    root, registry = repo
    text = _TABLE
    if change == "duplicate_key":
        text = text.replace("b,2,", "a,2,")
    elif change == "blank_key":
        text = text.replace("b,2,", ",2,")
    elif change == "duplicate_header":
        text = text.replace("reference,", "value,", 1)
    elif change == "no_match":
        registry["claims"][0]["evaluation"]["where"] = {"key": "missing"}
    elif change == "missing_column":
        registry["claims"][0]["evaluation"]["column"] = "not_a_column"
    elif change == "missing_independent_unit":
        text = text.replace(",p2,", ",,")
    elif change == "missing_row":
        text = text[:text.rfind("c,99,")]
    else:
        text += "malformed,row\n"
    (root / "outputs/table.csv").write_text(text)
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"] and not row["supported"]
    assert message in row["detail"]


@pytest.mark.parametrize("value", ["", "NaN", "inf", "-inf", "not-a-number"])
def test_numeric_aggregation_does_not_drop_missing_or_nonfinite_values(repo, value):
    root, registry = repo
    (root / "outputs/table.csv").write_text(_TABLE.replace("b,2,", f"b,{value},"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"]
    assert row["actual"] is None and row["observed"] is None
    assert row["comparison_status"] == "NOT_EVALUATED"


@pytest.mark.parametrize("aggregation, expected", [("mean", 34), ("median", 2), ("max", 99), ("min", 1), ("sum", 102), ("count", 3), ("nunique", 3), ("weighted_mean", 302 / 6)])
def test_declared_aggregation_is_used_not_a_convenient_row(repo, aggregation, expected):
    root, registry = repo
    registry["claims"][0]["evaluation"].update(aggregation=aggregation, weight_column="weight")
    registry["claims"][0]["expected"] = {"value": expected, "atol": 1e-10}
    _save(root, registry)
    row = _record(root)
    assert row["status"] == "PASS", row["detail"]
    assert row["actual"] == pytest.approx(expected)


def test_relative_error_and_symmetric_fold_keep_all_targets(repo):
    root, registry = repo
    evaluation = registry["claims"][0]["evaluation"]
    evaluation.update(transform={"kind": "absolute_relative_error", "reference_column": "reference"}, scale=100)
    registry["claims"][0]["expected"] = {"value": 100}
    _save(root, registry)
    assert _record(root)["actual"] == 100
    evaluation.update(aggregation="max", scale=1,
                      transform={"kind": "symmetric_fold_error", "reference_column": "reference"})
    registry["claims"][0]["expected"] = {"value": 99}
    _save(root, registry)
    assert _record(root)["actual"] == 99


def test_nonpositive_reference_cannot_be_fabricated_into_a_fold(repo):
    root, registry = repo
    (root / "outputs/table.csv").write_text(_TABLE.replace("b,2,1,", "b,2,0,"))
    registry["claims"][0]["evaluation"]["transform"] = {"kind": "ratio", "reference_column": "reference"}
    _save(root, registry, pin=True)
    assert "positive reference" in _record(root)["detail"]


@pytest.mark.parametrize("status", ["pending", "failed"])
def test_provenance_failure_blocks_even_a_matching_numeric_value(repo, status):
    root, registry = repo
    registry["profiles"]["fixture"]["provenance"] = {"status": status, "reason": "regeneration was not verified"}
    _save(root, registry)
    row = _record(root)
    assert row["comparison_status"] == "PASS"
    assert row["observed"] == 2
    assert row["actual"] is None and not row["supported"] and row["blocking"]
    assert all(value["evidence_status"] == "UNVERIFIED" for value in row["result_rows"])
    assert f"provenance {status}" in row["detail"]


def test_unchanged_values_with_changed_producer_fail_provenance(repo):
    root, _ = repo
    (root / "src/model.py").write_text("MODEL = 'm1'; CHANGED = True\n")
    row = _record(root)
    assert row["comparison_status"] == "PASS" and row["actual"] is None
    assert "provenance mismatch" in row["detail"] and "do not restamp" in row["detail"]


def test_bare_verified_assertion_is_not_provenance(repo):
    root, registry = repo
    registry["profiles"]["fixture"]["provenance"]["artifacts"] = []
    _save(root, registry)
    assert "pre-existing artifact identities" in _record(root)["detail"]


def test_evidence_and_producers_both_need_provenance_bindings(repo):
    root, registry = repo
    profile = registry["profiles"]["fixture"]
    profile["provenance"]["artifacts"] = profile["provenance"]["artifacts"][:1]
    _save(root, registry)
    assert "unbound producer" in _record(root)["detail"]
    profile["provenance"]["artifacts"] = [{"path": "src/model.py", "sha256": hashlib.sha256((root / "src/model.py").read_bytes()).hexdigest()}]
    _save(root, registry)
    assert "unbound evidence" in _record(root)["detail"]


def test_model_and_configuration_are_checked_even_if_values_match(repo):
    root, registry = repo
    (root / "outputs/table.csv").write_text(_TABLE.replace("v1,m1", "v2,m2"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["comparison_status"] == "PASS"
    assert row["blocking"] and row["actual"] is None
    assert "configuration mismatch" in row["detail"]


def test_pending_investigation_cannot_be_promoted_by_numeric_agreement(repo):
    root, registry = repo
    registry["claims"][0]["role"] = "investigation_pending"
    _save(root, registry)
    row = _record(root)
    assert row["observed"] == 2 and row["actual"] is None
    assert row["verdict"] == "BLOCKED"


def test_historical_claim_never_becomes_current_support(repo):
    root, registry = repo
    registry["claims"][0]["role"] = "historical"
    registry["claims"][0]["expected"]["value"] = 500
    _save(root, registry)
    row = _record(root)
    assert row["status"] == "SKIP" and row["verdict"] == "HISTORICAL"
    assert not row["supported"] and not row["blocking"] and row["actual"] is None


def _refusal(root, registry):
    payload = {"curation": {"config": "v1", "model": "m1", "status": "refused", "prediction": None,
                            "reason": "No independently established chemistry"}}
    (root / "data/refusal.json").write_text(json.dumps(payload))
    profile = registry["profiles"]["fixture"]
    profile.update(evidence=[{"id": "main", "path": "data/refusal.json", "format": "json", "pointer": "/curation",
                              "row_keys": [], "unit": "qualitative refusal"}],
                   coverage={"expected_rows": 1}, independent_unit={"name": "curation record"})
    claim = registry["claims"][0]
    claim.update(role="refused", evaluation={"kind": "refusal"}, acceptance={"predicates": [
        {"field": "status", "value": "refused"}, {"field": "prediction", "operator": "missing"},
        {"field": "reason", "operator": "nonempty"}]})
    claim.pop("expected")
    _save(root, registry, pin=True)
    return payload


def test_refusal_is_evaluated_from_fields_not_fabricated_as_zero(repo):
    root, registry = repo
    _refusal(root, registry)
    row = _record(root)
    assert row["status"] == "PASS" and row["verdict"] == "REFUSAL_CONFIRMED"
    assert row["actual"] == {"predicates_satisfied": True, "checked": 3}
    assert row["result_rows"][1]["value"] is None


def test_refusal_detects_an_invented_numeric_prediction(repo):
    root, registry = repo
    payload = _refusal(root, registry)
    payload["curation"]["prediction"] = 0
    (root / "data/refusal.json").write_text(json.dumps(payload))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["verdict"] == "NOT_SUPPORTED" and not row["supported"]


def test_failed_provenance_blocks_qualitative_refusals_too(repo):
    root, registry = repo
    _refusal(root, registry)
    (root / "data/refusal.json").write_text((root / "data/refusal.json").read_text() + "\n")
    row = _record(root)
    assert row["comparison_status"] == "PASS" and row["verdict"] == "BLOCKED"


def test_refused_role_cannot_use_numeric_evaluator(repo):
    root, registry = repo
    registry["claims"][0]["role"] = "refused"
    _save(root, registry)
    assert "refusal evaluator" in audit_claims(root)[0]["detail"]


def _equivalence_claim(registry):
    profile = registry["profiles"]["fixture"]
    profile["independent_unit"] = {"name": "one interval over independent plates"}
    claim = registry["claims"][0]
    claim["evaluation"] = {"kind": "equivalence", "where": {"key": "a"}, "expected_rows": 1,
                           "value_column": "value", "lower_column": "low", "upper_column": "high"}
    claim["acceptance"] = {"margins": [0.8, 1.25], "null": 1, "confidence": 0.95,
                           "margin_source": "Prespecified fixture SESOI, not chosen from the observed interval"}
    claim.pop("expected")
    return claim


def test_equivalence_requires_the_full_interval_inside_declared_margins(repo):
    root, registry = repo
    _equivalence_claim(registry)
    _save(root, registry)
    row = _record(root)
    assert row["supported"] and row["actual"]["equivalent"]


@pytest.mark.parametrize("low,high", [(0.7, 1.1), (0.9, 1.5), (0.1, 4), (0.8, 1.25)])
def test_interval_containing_one_is_not_equivalence(repo, low, high):
    root, registry = repo
    _equivalence_claim(registry)
    (root / "outputs/table.csv").write_text(_TABLE.replace("0.9,1.1", f"{low},{high}"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["verdict"] == "NOT_SUPPORTED"
    assert not row["actual"]["equivalent"]


@pytest.mark.parametrize("change", ["no_margin", "no_justification", "reversed_interval", "nonfinite_interval"])
def test_equivalence_refuses_missing_margin_or_uncertainty(repo, change):
    root, registry = repo
    claim = _equivalence_claim(registry)
    if change == "no_margin":
        claim["acceptance"]["margins"] = None
    elif change == "no_justification":
        claim["acceptance"].pop("margin_source")
    elif change == "reversed_interval":
        (root / "outputs/table.csv").write_text(_TABLE.replace("0.9,1.1", "1.1,0.9"))
    else:
        (root / "outputs/table.csv").write_text(_TABLE.replace("0.9,1.1", "0.9,NaN"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["verdict"] == "BLOCKED" and row["actual"] is None


def _nis_claim(root, registry):
    table = (
        "key,plate,channel,k_steps,mean_nis,lag1_autocorr,min_ess,tested,config,model\n"
        "a,p1,od,5,0.3,0.8,200,True,v1,m1\n"
        "b,p2,od,25,2,0.8,200,True,v1,m1\n"
        "c,p3,od,100,20,0.8,200,True,v1,m1\n"
        "d,p3,od,25,,,5,False,v1,m1\n"
    )
    (root / "outputs/table.csv").write_text(table)
    profile = registry["profiles"]["fixture"]
    profile["coverage"] = {"expected_rows": 4}
    claim = registry["claims"][0]
    claim.update(evaluation={"kind": "nis", "where": {"channel": "od"}, "alpha": 0.05,
                             "ess_floor": 100, "white_fraction_threshold": 0.5},
                 acceptance={"nis_verdict": "over_confident", "whiteness_verdict": "not_white"})
    claim.pop("expected")
    _save(root, registry, pin=True)
    return table


def test_nis_aggregates_eligible_wells_and_uses_each_actual_step_count(repo):
    root, registry = repo
    _nis_claim(root, registry)
    row = _record(root)
    assert row["status"] == "PASS", row["detail"]
    assert row["actual"]["median_mean_nis"] == 2
    assert row["actual"]["tested_wells"] == 3
    assert row["actual"]["inconclusive_wells"] == 1
    assert row["actual"]["fraction_white"] == pytest.approx(1 / 3)
    assert len({result["band_high"] for result in row["result_rows"] if result["tested"]}) == 3
    assert row["coverage_observed"]["independent_units"] == 3
    assert "serial dependence" in row["actual"]["interpretation"]


def test_nis_no_ess_eligible_wells_is_inconclusive_not_zero_or_pass(repo):
    root, registry = repo
    table = _nis_claim(root, registry)
    (root / "outputs/table.csv").write_text(table.replace("200,True", "5,False"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"] and row["actual"] is None
    assert "inconclusive" in row["detail"]
    assert row["observed"]["median_mean_nis"] is None
    assert row["observed"]["inconclusive_wells"] == 4
    assert len(row["result_rows"]) == 4


def test_nis_flags_inconsistent_eligibility_and_requires_whiteness(repo):
    root, registry = repo
    table = _nis_claim(root, registry)
    (root / "outputs/table.csv").write_text(table.replace("200,True", "5,True"))
    _save(root, registry, pin=True)
    assert "tested flag" in _record(root)["detail"]
    (root / "outputs/table.csv").write_text(table)
    registry["claims"][0]["acceptance"].pop("whiteness_verdict")
    _save(root, registry, pin=True)
    assert "assessed together" in _record(root)["detail"]


def test_marker_selectors_are_authoritative_and_ids_survive_line_changes(repo):
    root, registry = repo
    _add_marker(registry)
    (root / "README.md").write_text(_marker() + "\n")
    _save(root, registry)
    first = next(row for row in audit_claims(root) if row["record_type"] == "marked_value")
    (root / "README.md").write_text("Historical 123 tests in 2020, section 2.718.\n\n" + _marker() + "\n")
    second = next(row for row in audit_claims(root) if row["record_type"] == "marked_value")
    assert first["id"] == second["id"]
    assert first["location"]["line"] == 1 and second["location"]["line"] == 3
    assert second["status"] == "PASS" and second["actual"] == 1


def test_numeric_selector_identity_normalizes_order_and_decimal_serialization():
    assert selector_identity("outputs/a.csv", "x", {"a": 0.2, "b": "case"}) == selector_identity(
        "outputs/a.csv", "x", {"b": "case", "a": "0.20"})
    assert selector_identity("outputs/a.csv", "x", {"a": "case"}) != selector_identity(
        "outputs/a.csv", "x", {"a": "other"})


def test_multiple_markers_on_one_line_keep_separate_literals_and_scales(repo):
    root, registry = repo
    _add_marker(registry, scale=100)
    _add_marker(registry, {"key": "b"})
    (root / "README.md").write_text(_marker("1.00e2", extra=" scale=100") + " and " + _marker("2.00", "key=b"))
    _save(root, registry)
    rows = [row for row in audit_claims(root) if row["record_type"] == "marked_value"]
    assert len(rows) == 2 and all(row["status"] == "PASS" for row in rows)
    assert [row["actual"] for row in rows] == [100, 2]


@pytest.mark.parametrize("marker", [
    _marker("not a value; year 2025"),
    _marker(row="key=a;key=a"),
    _marker(extra=" scale=nan"),
    _marker(extra=" column=value"),
    '<!-- audit:value table=outputs/table.csv column=value row="key=a"',
])
def test_malformed_explicit_markers_block_instead_of_becoming_exemptions(repo, marker):
    root, registry = repo
    _add_marker(registry)
    (root / "README.md").write_text(marker)
    _save(root, registry)
    assert any(row["blocking"] and row["record_type"] == "scope" for row in audit_claims(root))


def test_new_selector_and_removed_selector_both_require_inventory_review(repo):
    root, registry = repo
    _add_marker(registry)
    (root / "README.md").write_text(_marker("2", "key=b"))
    _save(root, registry)
    rows = audit_claims(root)
    assert any("unregistered authoritative selector" in row["detail"] for row in rows)
    assert any("no remaining surface marker" in row["detail"] for row in rows)


def test_duplicate_selector_declarations_fail(repo):
    root, registry = repo
    _add_marker(registry)
    _add_marker(registry)
    _save(root, registry)
    assert "duplicate marker selector" in audit_claims(root)[0]["detail"]


def test_archival_marked_prose_is_not_current_support(repo):
    root, registry = repo
    _add_marker(registry)
    (root / "docs/history.md").write_text(_marker("999.0"))
    registry["surfaces"].append({"id": "archive", "glob": "docs/history.md", "kind": "markdown",
                                 "role": "historical", "unmarked_policy": "registered_only", "reason": "Retained old claim"})
    _save(root, registry)
    row = next(row for row in audit_claims(root) if row["record_type"] == "marked_value")
    assert row["verdict"] == "HISTORICAL" and not row["supported"] and row["actual"] is None


def test_ast_scope_distinguishes_fixture_literals_from_docstring_claims(repo):
    root, registry = repo
    _add_marker(registry)
    fixture = _marker()
    (root / "tests/example.py").write_text(f"FIXTURE = {fixture!r}\n")
    registry["surfaces"].append({"id": "test", "glob": "tests/**/*.py", "kind": "python_docstrings",
                                 "role": "forbidden", "unmarked_policy": "registered_only", "reason": "No scientific test-docstring claims"})
    (root / "README.md").write_text(fixture)
    _save(root, registry)
    rows = audit_claims(root)
    assert len([row for row in rows if row["record_type"] == "marked_value"]) == 1
    (root / "tests/example.py").write_text(f"'''{fixture}'''\n")
    assert any("forbidden surface" in row["detail"] for row in audit_claims(root))


def test_docstring_and_field_surfaces_are_classified_by_python_syntax(repo):
    root, registry = repo
    _add_marker(registry)
    fixture = _marker()
    (root / "src/claims.py").write_text(f"'''{fixture}'''\nFIELD_LITERAL = '#: not a field comment'\n#: {fixture}\n")
    for kind in ("python_docstrings", "python_field_comments"):
        registry["surfaces"].append({"id": kind, "glob": "src/**/*.py", "kind": kind, "role": "current",
                                     "unmarked_policy": "registered_only", "reason": "Explicit Python documentation surface"})
    _save(root, registry)
    rows = [row for row in audit_claims(root) if row["record_type"] == "marked_value"]
    assert len(rows) == 2 and all(row["status"] == "PASS" for row in rows)


def test_unparseable_declared_python_surface_is_blocking(repo):
    root, registry = repo
    (root / "src/broken.py").write_text("def broken(:\n")
    registry["surfaces"].append({"id": "python", "glob": "src/**/*.py", "kind": "python_docstrings", "role": "current",
                                 "unmarked_policy": "registered_only", "reason": "A broken AST is not a verified empty surface"})
    _save(root, registry)
    assert any(row["id"] == "scope.python" and row["blocking"] for row in audit_claims(root))


def test_relative_registry_and_evidence_paths_are_relative_to_root_not_cwd(repo, monkeypatch, tmp_path):
    root, registry = repo
    _save(root, registry, relative="config/claims.json")
    monkeypatch.chdir(tmp_path)
    assert _record(root, registry_path="config/claims.json")["actual"] == 2
    assert _record(root, registry_path=root / "config/claims.json")["actual"] == 2


@pytest.mark.parametrize("path", ["../outside.csv", "/outside.csv", "C:\\outside.csv", "outputs/../outside.csv"])
def test_evidence_paths_cannot_escape_declared_repository(repo, path):
    root, registry = repo
    registry["profiles"]["fixture"]["evidence"][0]["path"] = path
    _save(root, registry)
    assert audit_claims(root)[0]["blocking"]
    with pytest.raises(ClaimRegistryError):
        load_registry(root)


def test_symlinked_evidence_cannot_escape_root(repo, tmp_path):
    root, registry = repo
    outside = tmp_path / "outside.csv"
    outside.write_text(_TABLE)
    (root / "outputs/link.csv").symlink_to(outside)
    registry["profiles"]["fixture"]["evidence"][0]["path"] = "outputs/link.csv"
    _save(root, registry)
    row = _record(root)
    assert row["blocking"] and "escapes repository" in row["detail"]


def test_render_requires_explicit_output_dir_and_makes_no_implicit_files(repo, renderer):
    root, _ = repo
    before = set(root.rglob("*"))
    with pytest.raises(SystemExit) as error:
        renderer.main(["--root", str(root)])
    assert error.value.code == 2
    assert set(root.rglob("*")) == before


def test_renderer_writes_only_an_isolated_dir_and_never_restamps(repo, renderer):
    root, registry = repo
    stamp = root / "outputs/provenance.json"
    stamp.write_text('{"unchanged": "historical"}')
    before = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
    assert renderer.main(["--root", str(root), "--output-dir", "results/claims"]) == 0
    destination = root / "results/claims"
    assert {path.name for path in destination.iterdir()} == {"claim_inventory.csv", "claim_results.csv", "claim_report.md", "claim_verdicts.json"}
    assert all((root / path).read_bytes() == value for path, value in before.items())
    assert "Supported current result" in (destination / "claim_report.md").read_text()
    assert json.loads((destination / "claim_verdicts.json").read_text())[0]["actual"] == 2
    assert renderer.main(["--root", str(root), "--output-dir", "results/claims"]) == 2
    assert renderer.main(["--root", str(root), "--output-dir", "outputs"]) == 2
    assert all((root / path).read_bytes() == value for path, value in before.items())


def test_renderer_preserves_blockers_and_labels_cached_results_as_unverified(repo, renderer):
    root, registry = repo
    registry["profiles"]["fixture"]["provenance"] = {"status": "pending", "reason": "regenerate first"}
    _save(root, registry)
    assert renderer.main(["--root", str(root), "--output-dir", "pending-report"]) == 1
    content = (root / "pending-report/claim_report.md").read_text()
    assert "NOT AUTHORIZED" in content and "regenerate first" in content
    assert "UNVERIFIED" in (root / "pending-report/claim_results.csv").read_text()


@pytest.mark.parametrize("field,value", [("artifacts", None), ("artifacts", [None]), ("checks", None), ("checks", [None])])
def test_malformed_provenance_or_configuration_returns_a_blocking_record(repo, field, value):
    root, registry = repo
    profile = registry["profiles"]["fixture"]
    target = profile["provenance"] if field == "artifacts" else profile["configuration"]["binding"]
    target[field] = value
    _save(root, registry)
    rows = audit_claims(root)
    assert rows[0]["id"] == "registry.schema" and rows[0]["blocking"]


def test_equal_row_count_is_not_equal_population_coverage(repo):
    root, registry = repo
    registry["profiles"]["fixture"]["coverage"]["expected_keys"] = [{"key": key} for key in ("a", "b", "c")]
    (root / "outputs/table.csv").write_text(_TABLE.replace("c,99,", "other,99,"))
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"] and "declared row identities" in row["detail"]


def test_marker_scale_cannot_change_without_changing_its_declaration(repo):
    root, registry = repo
    _add_marker(registry, scale=100)
    (root / "README.md").write_text(_marker("0", extra=" scale=0"))
    _save(root, registry)
    rows = audit_claims(root)
    assert any("unregistered authoritative selector" in row["detail"] for row in rows)
    assert any("no remaining surface marker" in row["detail"] for row in rows)


def test_archive_alias_cannot_replace_a_declared_current_marker(repo):
    root, registry = repo
    _add_marker(registry)
    (root / "docs/history.md").write_text(_marker())
    registry["surfaces"].append({"id": "archive", "glob": "docs/history.md", "kind": "markdown",
                                 "role": "historical", "unmarked_policy": "registered_only", "reason": "Retained history"})
    _save(root, registry)
    assert any("no remaining surface marker" in row["detail"] for row in audit_claims(root))


def test_archived_unregistered_markers_do_not_become_current_claims(repo):
    root, registry = repo
    (root / "docs/history.md").write_text(_marker())
    registry["surfaces"].append({"id": "archive", "glob": "docs/history.md", "kind": "markdown",
                                 "role": "historical", "unmarked_policy": "registered_only", "reason": "Retained history"})
    _save(root, registry)
    rows = audit_claims(root)
    assert not any(row["blocking"] for row in rows)
    assert next(row for row in rows if row["record_type"] == "marked_value")["verdict"] == "HISTORICAL"


@pytest.mark.parametrize("destination", ["data/generated", "docs/generated", "src/generated", "scripts/generated", "tests/generated", ".git/generated"])
def test_renderer_cannot_use_repository_source_or_documentation_trees(repo, renderer, destination):
    root, _ = repo
    assert renderer.main(["--root", str(root), "--output-dir", destination]) == 2
    assert not (root / destination).exists()


@pytest.fixture
def receipt_repo(repo):
    root, registry = repo
    (root / "scripts").mkdir()
    (root / "data/source.csv").write_text(_TABLE)
    (root / "scripts/run.py").write_text(
        "from pathlib import Path\n"
        "Path('outputs/table.csv').write_text(Path('data/source.csv').read_text())\n")
    run = {
        "id": "fixture-run", "lifecycle": "current", "status": "investigation_pending",
        "reason": "Synthetic execution fixture, not retained scientific output",
        "producer": {"path": "scripts/run.py", "callable": "main", "libraries": ["src/model.py"]},
        "dependencies": [{"path": "data/source.csv", "role": "input", "sha256": None}],
        "parameters": {"config": "v1", "selection": "all"}, "model_identity": {"id": "m1"},
        "runtime": {"kind": "python_standard_library"},
        "reproduction": {"argv": ["{python}", "-B", "scripts/run.py"]},
        "artifacts": [{"path": "outputs/table.csv", "comparison": {
            "format": "csv", "row_keys": ["key"],
            "units": {name: "dimensionless" for name in ("value", "reference", "low", "high", "weight")}}}],
    }
    document = {"schema_version": 1, "runs": [run], "frozen_bundles": []}
    (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    profile = registry["profiles"]["fixture"]
    profile["provenance"] = {"source": "artifact_registry", "run_id": "fixture-run"}
    profile["configuration"]["binding"] = {"source": "artifact_registry", "checks": [
        {"field": "/parameters/config", "value": "v1"}, {"field": "/model_identity/id", "value": "m1"},
        {"source": "main", "field": "config", "value": "v1"}]}
    _save(root, registry)
    return root, registry, document


def _adopt_fixture(receipt_repo, destination):
    root, _, _ = receipt_repo
    receipt = artifacts.reproduce_run(root, "fixture-run", destination)
    assert receipt["success"], receipt
    result = artifacts.adopt_reproduction(root, destination / "reproduction.json")
    return receipt, result


def test_current_claim_reads_adopted_durable_receipt_without_manual_pins(receipt_repo, tmp_path, monkeypatch):
    root, registry, _ = receipt_repo
    before = _record(root)
    assert before["blocking"] and before["comparison_status"] == "PASS"
    assert "durable reproduction receipt absent" in before["detail"]
    _, adoption = _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    shutil.rmtree(tmp_path / "reproduction")
    snapshot = {path.relative_to(root): (path.read_bytes(), path.stat().st_mtime_ns)
                for path in root.rglob("*") if path.is_file()}

    def forbidden(*args, **kwargs):
        pytest.fail("a read-only claim audit must never reproduce or adopt")

    monkeypatch.setattr(artifacts, "reproduce_run", forbidden)
    monkeypatch.setattr(artifacts, "adopt_reproduction", forbidden)
    row = _record(root)
    assert row["status"] == "PASS" and row["actual"] == 2
    assert row["provenance"] == {"source": "artifact_registry", "run_id": "fixture-run"}
    observed = row["provenance_observed"]
    assert observed["status"] == "verified"
    assert observed["artifacts"][0]["receipt_id"] == adoption["receipt_id"]
    assert observed["identity"]["parameters"]["config"] == "v1"
    assert snapshot == {path.relative_to(root): (path.read_bytes(), path.stat().st_mtime_ns)
                        for path in root.rglob("*") if path.is_file()}
    assert json.loads((root / "data/current_claims.json").read_text()) == registry


def test_scientific_change_needs_reviewed_and_installed_receipt_before_claim_support(receipt_repo, tmp_path):
    root, registry, _ = receipt_repo
    (root / "data/source.csv").write_text(_TABLE.replace("b,2,", "b,3,"))
    stage = tmp_path / "reproduction"
    receipt = artifacts.reproduce_run(root, "fixture-run", stage)
    assert receipt["execution_complete"] and not receipt["success"]
    receipt_path = stage / "reproduction.json"
    review = artifacts.review_template(root, receipt_path)
    review.update(approved=True, reviewer="synthetic fixture reviewer", reason="Test-only changed median")
    for approval in review["artifacts"].values():
        approval["reason"] = "Reviewed the exact synthetic baseline and candidate"
        approval["accepted_dimensions"] = approval["observed_changed_dimensions"]
    review_path = tmp_path / "review.json"
    review_path.write_text(json.dumps(review))
    registry["claims"][0]["expected"]["value"] = 3
    _save(root, registry)
    artifacts.adopt_reproduction(root, receipt_path, accept_scientific_changes=True, review_path=review_path)
    row = _record(root)
    assert row["blocking"] and row["actual"] is None
    assert "reviewed candidate not installed" in row["detail"]
    artifacts.adopt_reproduction(root, receipt_path, accept_scientific_changes=True,
                                review_path=review_path, install=True)
    row = _record(root)
    assert row["status"] == "PASS" and row["actual"] == 3
    assert row["provenance_observed"]["artifacts"][0]["mode"] == "reviewed_change"


def test_shared_receipt_validator_is_called_once_per_artifact_per_audit(receipt_repo, tmp_path, monkeypatch):
    root, registry, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    _add_marker(registry)
    (root / "README.md").write_text(_marker() + " and " + _marker())
    _save(root, registry)
    original = artifacts._receipt_issues
    calls = []

    def spy(root, registry, run, item, current):
        calls.append(item["path"])
        return original(root, registry, run, item, current)

    monkeypatch.setattr(artifacts, "_receipt_issues", spy)
    rows = audit_claims(root)
    assert not any(row["blocking"] for row in rows)
    assert len([row for row in rows if row["record_type"] == "marked_value"]) == 2
    assert calls == ["outputs/table.csv"]


def test_verified_receipt_can_confirm_a_scientific_refusal(receipt_repo, tmp_path):
    root, registry, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    claim = registry["claims"][0]
    claim.update(role="refused", evaluation={"kind": "refusal"}, acceptance={
        "predicates": [{"field": "reference", "value": 1}]})
    claim.pop("expected")
    _save(root, registry)
    row = _record(root)
    assert row["status"] == "PASS" and row["verdict"] == "REFUSAL_CONFIRMED"
    assert row["provenance_observed"]["status"] == "verified"


def test_receipt_configuration_checks_cannot_be_removed_after_successful_adoption(receipt_repo, tmp_path):
    root, registry, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    registry["profiles"]["fixture"]["configuration"]["binding"]["checks"] = []
    _save(root, registry)
    row = _record(root)
    assert row["blocking"] and row["comparison_status"] == "PASS"
    assert "configuration needs authoritative field checks" in row["detail"]


def test_partial_run_adoption_can_authorize_only_its_verified_evidence(receipt_repo, tmp_path):
    root, registry, document = receipt_repo
    run = document["runs"][0]
    other = deepcopy(run["artifacts"][0])
    other["path"] = "outputs/other.csv"
    run["artifacts"].append(other)
    (root / "outputs/other.csv").write_text(_TABLE)
    script = root / "scripts/run.py"
    script.write_text(script.read_text() + "Path('outputs/other.csv').write_text(Path('data/source.csv').read_text())\n")
    (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    stage = tmp_path / "reproduction"
    receipt = artifacts.reproduce_run(root, "fixture-run", stage)
    assert receipt["success"], receipt
    artifacts.adopt_reproduction(root, stage / "reproduction.json", artifact_ids=["outputs/table.csv"])
    assert artifacts.load_registry(root)["runs"]["fixture-run"]["status"] == "investigation_pending"
    assert _record(root)["status"] == "PASS"
    registry["profiles"]["fixture"]["evidence"].append({
        "id": "other", "path": "outputs/other.csv", "format": "csv", "row_keys": ["key"], "unit": "dimensionless"})
    _save(root, registry)
    assert "durable reproduction receipt absent: outputs/other.csv" in _record(root)["detail"]


@pytest.mark.parametrize("change", ["producer", "input", "output", "parameters", "model", "runtime", "baseline"])
def test_receipt_freshness_is_validated_by_artifact_api(receipt_repo, tmp_path, monkeypatch, change):
    root, _, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    if change in {"producer", "input", "output"}:
        path = root / {"producer": "src/model.py", "input": "data/source.csv", "output": "outputs/table.csv"}[change]
        path.write_text(path.read_text() + "\n")
    elif change == "runtime":
        original = artifacts.runtime_identity
        monkeypatch.setattr(artifacts, "runtime_identity", lambda root: {**original(root), "python": "different"})
    else:
        document = json.loads((root / artifacts.REGISTRY_PATH).read_text())
        if change == "parameters":
            document["runs"][0]["parameters"]["config"] = "v2"
        elif change == "model":
            document["runs"][0]["model_identity"]["id"] = "m2"
        else:
            document["baseline_ref"] = "different-baseline-policy"
        (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    row = _record(root)
    assert row["blocking"] and row["actual"] is None
    assert row["observed"] == 2 and row["comparison_status"] == "PASS"
    assert "provenance unverified" in row["detail"]


def test_run_verified_status_without_a_receipt_cannot_authorize_a_claim(receipt_repo):
    root, _, document = receipt_repo
    document["runs"][0]["status"] = "verified"
    (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    assert "durable reproduction receipt absent" in _record(root)["detail"]


@pytest.mark.parametrize("change", ["receipt_hash", "failed_execution", "unwitnessed_output", "semantic_failure"])
def test_tampered_durable_proof_is_not_a_manual_freshness_assertion(receipt_repo, tmp_path, change):
    root, _, _ = receipt_repo
    _, adoption = _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    document = json.loads((root / artifacts.REGISTRY_PATH).read_text())
    reference = document["runs"][0]["verification"]["artifacts"]["outputs/table.csv"]
    receipt = document["receipts"][adoption["receipt_id"]]
    if change == "receipt_hash":
        receipt["success"] = False
    elif change == "unwitnessed_output":
        reference["accepted_sha256"] = "0" * 64
    else:
        if change == "failed_execution":
            receipt["exit_code"] = 9
        else:
            receipt["candidate_validation"]["outputs/table.csv"]["matches"] = False
        replacement = artifacts.canonical_digest(receipt)
        document["receipts"][replacement] = receipt
        reference["receipt_id"] = replacement
    (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    row = _record(root)
    assert row["blocking"] and row["actual"] is None and row["comparison_status"] == "PASS"
    assert "artifact receipt invalid" in row["detail"]


@pytest.mark.parametrize("change", ["unknown_run", "other_member", "unbound_producer", "noncurrent_run"])
def test_receipts_cannot_authorize_another_run_or_unbound_method(receipt_repo, tmp_path, change):
    root, registry, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    profile = registry["profiles"]["fixture"]
    if change == "unknown_run":
        profile["provenance"]["run_id"] = "other"
    elif change == "other_member":
        (root / "outputs/other.csv").write_text(_TABLE)
        profile["evidence"][0]["path"] = "outputs/other.csv"
    elif change == "unbound_producer":
        (root / "src/other.py").write_text("MODEL = 'm1'\n")
        profile["method"]["implementation"].append("src/other.py")
    else:
        document = json.loads((root / artifacts.REGISTRY_PATH).read_text())
        document["runs"][0]["lifecycle"] = "historical_frozen"
        (root / artifacts.REGISTRY_PATH).write_text(json.dumps(document))
    _save(root, registry)
    row = _record(root)
    assert row["blocking"] and row["actual"] is None and row["comparison_status"] == "PASS"


def test_configuration_is_checked_against_verified_receipt_not_only_copied_cells(receipt_repo, tmp_path):
    root, registry, _ = receipt_repo
    _adopt_fixture(receipt_repo, tmp_path / "reproduction")
    registry["profiles"]["fixture"]["configuration"]["binding"]["checks"][0]["value"] = "v2"
    _save(root, registry)
    row = _record(root)
    assert row["blocking"] and row["actual"] is None and row["observed"] == 2
    assert row["provenance_observed"]["status"] == "verified"
    assert "configuration mismatch: /parameters/config" in row["detail"]


@pytest.mark.parametrize("change", ["manual_status", "no_run", "unknown_source", "bad_pointer", "no_checks"])
def test_receipt_binding_schema_cannot_be_a_bare_status_or_empty_contract(receipt_repo, change):
    root, registry, _ = receipt_repo
    profile = registry["profiles"]["fixture"]
    if change == "manual_status":
        profile["provenance"]["status"] = "verified"
    elif change == "no_run":
        profile["provenance"].pop("run_id")
    elif change == "unknown_source":
        profile["provenance"]["source"] = "current_file_hash"
    elif change == "bad_pointer":
        profile["configuration"]["binding"]["checks"][0]["field"] = "config"
    else:
        profile["configuration"]["binding"]["checks"] = []
    _save(root, registry)
    assert any(row["blocking"] for row in audit_claims(root))
    if change != "no_checks":
        assert audit_claims(root)[0]["id"] == "registry.schema"


def _retire_marker(registry):
    entry = {"path": "outputs/table.csv", "where": {"key": "a"}, "columns": {"value": "dimensionless"},
             "reason": "Single-cell headline replaced with complete-population aggregation",
             "replacement_family": "fixture", "replacement_claims": ["fixture.median"]}
    registry["marked_values"]["retired_selectors"] = [entry]
    registry["marked_values"]["baseline_occurrences"] = 246
    registry["marked_values"]["baseline_distinct_selectors"] = 169
    return entry


def test_explicit_retirement_preserves_history_without_a_numeric_count_target(repo):
    root, registry = repo
    _retire_marker(registry)
    _save(root, registry)
    rows = audit_claims(root)
    assert not any(row["blocking"] for row in rows)
    retired = next(row for row in rows if row["record_type"] == "marker_retirement")
    assert retired["status"] == "SKIP" and retired["verdict"] == "RETIRED" and not retired["supported"]
    assert "fixture.median" in retired["detail"]
    summary = next(row for row in rows if row["id"] == "registry.marker_inventory")
    assert summary["actual"]["distinct_selectors"] == summary["expected"]["declared_selectors"] == 0
    assert summary["expected"]["baseline_distinct_selectors"] == 169
    assert summary["expected"]["baseline_occurrences"] == 246
    assert summary["expected"]["retired_selectors"] == 1


@pytest.mark.parametrize("change", ["reason", "family", "replacement", "historical_replacement", "duplicate", "still_active"])
def test_retirement_requires_reason_and_current_replacement_coverage(repo, change):
    root, registry = repo
    entry = _retire_marker(registry)
    if change == "reason":
        entry["reason"] = ""
    elif change == "family":
        entry["replacement_family"] = "absent"
    elif change == "replacement":
        entry["replacement_claims"] = []
    elif change == "historical_replacement":
        registry["claims"][0]["role"] = "historical"
    elif change == "duplicate":
        registry["marked_values"]["retired_selectors"].append(deepcopy(entry))
    else:
        _add_marker(registry)
    _save(root, registry)
    assert audit_claims(root)[0]["id"] == "registry.schema"


def test_retirement_never_exempts_a_reintroduced_current_marker(repo):
    root, registry = repo
    _retire_marker(registry)
    (root / "README.md").write_text(_marker())
    _save(root, registry)
    assert any("retired selector reintroduced as current" in row["detail"] and row["blocking"]
               for row in audit_claims(root))


def test_retired_selector_can_remain_on_explicitly_historical_surface(repo):
    root, registry = repo
    _retire_marker(registry)
    (root / "docs/history.md").write_text(_marker("999"))
    registry["surfaces"].append({"id": "archive", "glob": "docs/history.md", "kind": "markdown",
                                 "role": "historical", "unmarked_policy": "registered_only", "reason": "Retained history"})
    _save(root, registry)
    rows = audit_claims(root)
    assert not any(row["blocking"] for row in rows)
    assert next(row for row in rows if row["record_type"] == "marked_value")["status"] == "SKIP"


@pytest.mark.parametrize("channel", ["od", "rfu"])
@pytest.mark.parametrize("change", [None, "tested", "reason", "missing_well"])
def test_current_nis_contract_accepts_only_complete_recorded_refusal(repo, channel, change):
    root, registry = repo
    current = next(claim for claim in load_registry(_ROOT)["claims"] if claim["id"] == f"nis.current.{channel}")
    assert current["role"] == "refused" and current["evaluation"]["kind"] == "refusal"
    rows = [f"{name}-{index},p{index % 4},{name},False,ess_below_floor,v1,m1"
            for name in ("od", "rfu") for index in range(319)]
    selected = 0 if channel == "od" else 319
    if change == "tested":
        rows[selected] = rows[selected].replace("False", "True")
    elif change == "reason":
        rows[selected] = rows[selected].replace("ess_below_floor", "")
    elif change == "missing_well":
        rows.pop(selected)
    (root / "outputs/table.csv").write_text("key,plate,channel,tested,exclusion_reason,config,model\n" + "\n".join(rows) + "\n")
    profile = registry["profiles"]["fixture"]
    profile["coverage"] = deepcopy(current["coverage"])
    profile["independent_unit"] = deepcopy(current["independent_unit"])
    claim = registry["claims"][0]
    claim.update(role="refused", evaluation=deepcopy(current["evaluation"]), acceptance=deepcopy(current["acceptance"]))
    claim.pop("expected")
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"] == (change is not None)
    if change is None:
        assert row["verdict"] == "REFUSAL_CONFIRMED"
        assert row["coverage_observed"]["selected_rows"] == 319
        assert row["coverage_observed"]["independent_units"] == 4


@pytest.mark.parametrize("change", [None, "drop_biomass", "promote_feasible", "fabricate_growth", "invent_order", "no_reason"])
def test_current_cosmic_contract_keeps_both_complete_infeasible_task_sets(repo, change):
    root, registry = repo
    current = next(claim for claim in load_registry(_ROOT)["claims"] if claim["id"] == "cosmic.full_task_set.refused")
    payload = _refusal(root, registry)
    task = {"status": "infeasible", "tasks": ["biomass", "ethanol", "CO2"], "order": [],
            "constrained_growth": None, "refusal": "Full task set is infeasible"}
    payload["curation"]["models"] = {"plain": deepcopy(task), "ec": deepcopy(task)}
    ec = payload["curation"]["models"]["ec"]
    if change == "drop_biomass":
        ec["tasks"].remove("biomass")
    elif change == "promote_feasible":
        ec["status"] = "optimal"
    elif change == "fabricate_growth":
        ec["constrained_growth"] = 0
    elif change == "invent_order":
        ec["order"] = [{"task": "ethanol"}]
    elif change == "no_reason":
        ec["refusal"] = ""
    (root / "data/refusal.json").write_text(json.dumps(payload))
    registry["claims"][0]["acceptance"] = deepcopy(current["acceptance"])
    _save(root, registry, pin=True)
    row = _record(root)
    assert row["blocking"] == (change is not None)
    if change is None:
        assert row["verdict"] == "REFUSAL_CONFIRMED"


@pytest.mark.parametrize("change", [None, "estimable", "low", "high", "missing_fold", "missing_row"])
def test_invalid_bootstrap_interval_is_a_refusal_with_the_point_retained(repo, change):
    root, registry = repo
    current = next(claim for claim in load_registry(_ROOT)["claims"]
                   if claim["id"] == "biosensor.invalid_tail.refusal")
    values = {"fold": "1.0", "low": "", "high": "", "estimable": "False"}
    if change == "estimable":
        values[change] = "True"
    elif change in ("low", "high"):
        values[change] = "0.5"
    elif change == "missing_fold":
        values["fold"] = ""
    header = "key,construct,stressor,dose_mM,late_fraction,fold,low,high,estimable,config,model,plate\n"
    row = "a,AlteredYap1,H2O2,1,0.75," + ",".join(values.values()) + ",v1,m1,p1\n"
    (root / "outputs/table.csv").write_text(header + ("" if change == "missing_row" else row))
    profile = registry["profiles"]["fixture"]
    profile["coverage"] = {"expected_rows": 1}
    profile["independent_unit"]["expected_count"] = 1
    claim = registry["claims"][0]
    claim.update(role=current["role"], evaluation=deepcopy(current["evaluation"]),
                 acceptance=deepcopy(current["acceptance"]))
    claim.pop("expected")
    _save(root, registry, pin=True)
    result = _record(root)
    assert result["blocking"] == (change is not None)
    if change is None:
        assert result["verdict"] == "REFUSAL_CONFIRMED"


def test_invalid_bootstrap_endpoints_are_explicitly_retired_not_silently_dropped():
    registry = load_registry(_ROOT)
    where = {"construct": "AlteredYap1", "stressor": "H2O2", "dose_mM": 1, "late_fraction": 0.75}
    entry = next(entry for entry in registry["marked_values"]["retired_selectors"]
                 if entry["path"] == "outputs/late_window_sensitivity.csv" and entry["where"] == where)
    assert set(entry["columns"]) == {"low", "high"}
    assert entry["replacement_claims"] == ["biosensor.invalid_tail.refusal"]
    active = next(selector for selector in registry["marked_values"]["tables"][entry["path"]]["selectors"]
                  if selector["where"] == where)
    assert set(active["columns"]) == {"naive_fold", "fold"}


def test_unsupported_no_response_claims_are_withdrawn_without_inventing_margins():
    registry = load_registry(_ROOT)
    withdrawn = [claim for claim in registry["claims"] if claim["evaluation"]["kind"] == "equivalence"]
    assert len(withdrawn) == 5
    for claim in withdrawn:
        assert claim["role"] == "historical"
        assert claim["withdrawal_reason"]
        assert claim["acceptance"]["margins"] is None and claim["acceptance"]["margin_source"] is None


def test_current_generated_profiles_name_authoritative_artifact_runs_and_members():
    registry = load_registry(_ROOT)
    artifact_registry = artifacts.load_registry(_ROOT)
    for profile in registry["profiles"].values():
        provenance = profile["provenance"]
        if provenance.get("source") != "artifact_registry":
            continue
        run_id = provenance["run_id"]
        assert artifact_registry["runs"][run_id]["lifecycle"] == "current"
        assert "status" not in provenance and "artifacts" not in provenance
        assert profile["configuration"]["binding"]["source"] == "artifact_registry"
        for source in profile["evidence"]:
            assert artifact_registry["artifacts"][source["path"]]["run"] == run_id


def test_repository_inventory_has_all_declared_marker_selectors_and_both_targets():
    registry = load_registry(_ROOT)
    rows = audit_claims(_ROOT)
    ids = {claim["id"] for claim in registry["claims"]}
    assert {"product.loso.predictions", "product.loso.strains", "product.beta.median_error", "product.beta.worst_fold",
            "product.lycopene.median_error", "product.lycopene.worst_fold", "nis.current.od", "nis.current.rfu",
            "native.strict_chemistry.refused", "validation.synthetic.rotated_null", "validation.real.rotated_null"} <= ids
    markers = [row for row in rows if row["record_type"] == "marked_value"]
    summary = next(row for row in rows if row["id"] == "registry.marker_inventory")
    assert summary["actual"]["distinct_selectors"] == summary["expected"]["declared_selectors"]
    assert len(markers) == summary["actual"]["occurrences"]
    assert registry["marked_values"]["baseline_occurrences"] == 246
    assert registry["marked_values"]["baseline_distinct_selectors"] == 169
    declared = {selector_identity(path, column, selector["where"], spec.get("scale", 1) if isinstance(spec, dict) else 1)
                for path, table in registry["marked_values"]["tables"].items()
                for selector in table["selectors"] for column, spec in selector["columns"].items()}
    current = {row["selector_id"] for row in markers if row.get("role") != "historical"}
    assert declared == current
    assert not [row for row in rows if row["id"].startswith(("marker_inventory.", "scope.", "registry.schema"))]
    assert len({row["id"] for row in rows}) == len(rows)
    assert all(row["evaluation"]["aggregation"] == "identity" for row in markers)
    retired_widths = {row["selector_id"] for row in rows if row["record_type"] == "marker_retirement"
                      and row["selector"]["path"] == "outputs/content_band_coverage.csv"}
    widths = ("share_alpha", "share_crtyb_capacity", "flux_only_too_narrow_by",
              "with_kinetics_over_truth", "independent_over_joint")
    expected_widths = {selector_identity("outputs/content_band_coverage.csv", column, {"growth_rate_per_h": rate})
                       for rate, columns in ((0.101, widths), (0.15, ("with_kinetics_over_truth",)),
                                             (0.2, ("with_kinetics_over_truth",)), (0.254, widths))
                       for column in columns}
    assert retired_widths == expected_widths and len(retired_widths) == 12
    assert not retired_widths & declared
    claims = {claim["id"]: claim for claim in registry["claims"]}
    for entry in registry["marked_values"]["retired_selectors"]:
        assert entry["reason"] and entry["replacement_claims"]
        for identity in entry["replacement_claims"]:
            assert claims[identity]["role"] in {"current", "refused"}
            assert claims[identity]["family"] == entry["replacement_family"]
    for row in rows:
        if row.get("provenance", {}).get("status") == "pending" and row.get("role") != "historical":
            assert row["blocking"] and not row["supported"] and row["actual"] is None
        if (row.get("provenance", {}).get("source") == "artifact_registry"
                and row.get("supported")):
            assert row["provenance_observed"]["status"] == "verified"


@pytest.fixture
def linked_cross_family_script(tmp_path, monkeypatch):
    from ystwin import paths

    monkeypatch.setattr(paths, "outputs_dir", lambda: tmp_path)
    spec = importlib.util.spec_from_file_location("_cross_family_metadata_contract", _ROOT / "scripts/run_cross_family.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_linked_cross_family_profile_binds_actual_fixed_split_and_spota_models(linked_cross_family_script):
    script = linked_cross_family_script
    current = next(claim for claim in load_registry(_ROOT)["claims"] if claim["id"] == "validation.synthetic.rotated_null")
    run = artifacts.load_registry(_ROOT)["runs"]["cross-family"]
    checks = {(check.get("source", "receipt"), check["field"]): check["value"]
              for check in current["configuration"]["binding"]["checks"]}
    fixed = script.configuration((3, "recommended pairs"), (0.25, 0.5, 1.0, 2.0), 3)
    candidate = script.configuration((2, "recommended pairs"), (0.25, 0.5, 1.0, 2.0), 3)
    assert checks["receipt", "/parameters/spota_seed"] == 0
    assert ("receipt", "/parameters/candidate_seed") not in checks
    assert checks["receipt", "/parameters/transfer_configuration"] == json.loads(fixed)
    assert run["parameters"]["transfer_configuration"] == json.loads(fixed)
    assert checks["main", "model"] == checks["optimism", "model"] == script.MODEL
    assert checks["main", "configuration"] == checks["optimism", "transfer_configuration"] == fixed
    assert checks["optimism", "selected_configuration"] == fixed
    assert checks["optimism", "candidate_configuration"] == candidate
    assert checks["optimism", "candidate_states"] == 2
    assert checks["optimism", "configuration_grid"] == script.compact_json({"n_states": script.K_GRID, "designs": script.DESIGNS})
    split = json.loads(checks["optimism", "split_selection_scope"])
    assert split == {"training_families": list(script.TRAINING), "held_out_families": list(script.HELD_OUT),
                     "family_seed": 0, "data_seed": 0}
    bound = json.loads(checks["optimism", "bound_selection_scope"])
    assert bound == {"sampler": "ystwin.generator.families.sample_family", "family_registry": list(script.family_names()),
                     "family_weights": "uniform", "seed": 0, "candidate_fixed_for_references": True}
    assert checks["optimism", "selected_pipeline_optimism_status"] == script.PIPELINE_OPTIMISM_STATUS
    assert current["coverage"]["expected_keys"] == [{"family": name} for name in script.family_names()]
    assert {"scripts/run_cross_family.py", "src/ystwin/analysis/transfer.py", "src/ystwin/analysis/optimism.py"} <= set(current["method"]["implementation"])
    assert current["role"] == "current" and current["evaluation"] == {"kind": "qualitative"}
    assert current["acceptance"] == {"predicates": [{"field": "beats_null", "value": False}]}


def _cross_family_claim_fixture(root, registry, script):
    from ystwin.analysis.optimism import OptimismBound

    current = next(claim for claim in load_registry(_ROOT)["claims"] if claim["id"] == "validation.synthetic.rotated_null")
    fixed = script.configuration((3, "recommended pairs"), (0.25, 0.5, 1.0, 2.0), 3)
    transfer = [{"family": name, "beats_null": False, "model": script.MODEL,
                 "configuration": fixed, "target": script.SCORE_TARGET,
                 "evaluation_unit": "one synthetic family dataset; stressor folds share training data",
                 "selection_scope": "fixed dimension/co-doses; latent weights refit per outer training fold",
                 "oracle_role": "in-sample diagnostic, not a held-out prediction"} for name in script.family_names()]
    optimism = [{"candidate_design": "recommended pairs", "candidate_states": 2,
                 "alpha": 0.05, "n_bootstrap": 1000, "n_reference": 8,
                 "n_candidate_domains": 4, "n_reference_domains": 2,
                 "model": script.MODEL, "transfer_configuration": fixed, "selected_configuration": fixed,
                 "candidate_configuration": script.configuration((2, "recommended pairs"), (0.25, 0.5, 1.0, 2.0), 3),
                 "configuration_grid": script.compact_json({"n_states": script.K_GRID, "designs": script.DESIGNS}),
                 "score_target": script.SCORE_TARGET,
                 "split_target": "mean family transfer on named training/held-out sets; descriptive difference",
                 "split_evaluation_unit": "named synthetic family; fixed split, not iid family draws",
                 "split_selection_scope": script.compact_json({"training_families": script.TRAINING,
                                                              "held_out_families": script.HELD_OUT,
                                                              "family_seed": 0, "data_seed": 0}),
                 "selected_pipeline_optimism_status": script.PIPELINE_OPTIMISM_STATUS,
                 "bound_target": OptimismBound.target,
                 "bound_resampling_unit": "independent reference family set; one gap per reference solution",
                 "bound_selection_scope": script.compact_json({"sampler": "ystwin.generator.families.sample_family",
                                                              "family_registry": script.family_names(),
                                                              "family_weights": "uniform", "seed": 0,
                                                              "candidate_fixed_for_references": True}),
                 "bound_method": OptimismBound.method, "bound_coverage_status": OptimismBound.coverage_status}]
    profile = registry["profiles"]["fixture"]
    for key in ("conditions", "model", "coverage", "independent_unit", "evidence"):
        profile[key] = deepcopy(current[key])
    profile["configuration"] = {"id": current["configuration"]["id"], "binding": {"status": "verified", "checks": [
        deepcopy(check) for check in current["configuration"]["binding"]["checks"] if check.get("source") in {"main", "optimism"}]}}
    claim = registry["claims"][0]
    claim.update(role=current["role"], evaluation=deepcopy(current["evaluation"]), acceptance=deepcopy(current["acceptance"]))
    claim.pop("expected")
    return transfer, optimism


@pytest.mark.parametrize("change", [None, "transfer_width", "spota_width", "selected_width", "split_population", "bound_population",
                                    "pipeline_promoted", "oracle_promoted", "missing_configuration", "missing_family", "wrong_family", "rejects_null"])
def test_linked_cross_family_claim_rejects_scope_drift_not_merely_changed_numbers(repo, linked_cross_family_script, change):
    root, registry = repo
    script = linked_cross_family_script
    transfer, optimism = _cross_family_claim_fixture(root, registry, script)
    if change in {"transfer_width", "spota_width", "selected_width"}:
        row, field = ((transfer[0], "configuration") if change == "transfer_width" else
                      (optimism[0], "candidate_configuration" if change == "spota_width" else "selected_configuration"))
        configuration = json.loads(row[field])
        configuration["n_states"] = 3 if change == "spota_width" else 2
        row[field] = script.compact_json(configuration)
    elif change == "split_population":
        scope = json.loads(optimism[0]["split_selection_scope"])
        scope["training_families"].append(scope["held_out_families"].pop())
        optimism[0]["split_selection_scope"] = script.compact_json(scope)
    elif change == "bound_population":
        scope = json.loads(optimism[0]["bound_selection_scope"])
        scope["candidate_fixed_for_references"] = False
        optimism[0]["bound_selection_scope"] = script.compact_json(scope)
    elif change == "pipeline_promoted":
        optimism[0]["selected_pipeline_optimism_status"] = "validated"
    elif change == "oracle_promoted":
        transfer[0]["oracle_role"] = "independent held-out prediction"
    elif change == "missing_configuration":
        transfer[0]["configuration"] = ""
    elif change == "missing_family":
        transfer.pop()
    elif change == "wrong_family":
        transfer[-1]["family"] = "unregistered"
    elif change == "rejects_null":
        transfer[-1]["beats_null"] = True
    for relative, rows in (("outputs/cross_family_transfer.csv", transfer), ("outputs/cross_family_optimism.csv", optimism)):
        with (root / relative).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    _save(root, registry, pin=True)
    result = _record(root)
    assert result["blocking"] == (change is not None), result["detail"]
    if change is None:
        assert result["verdict"] == "SUPPORTED"
        assert "equivalent" not in result["actual"]
        assert result["coverage_observed"]["selected_rows"] == 7
    elif change == "rejects_null":
        assert result["verdict"] == "NOT_SUPPORTED"
    else:
        assert result["actual"] is None and not result["supported"]


_BOUNDED_SCOPE_CASES = [
    ("power.synthetic_scope", "power", 28, 7,
     "/artifacts/0/comparison/units/n_simulations"),
    ("maintenance.regime_specific", "maintenance-scale", 16, 16,
     "/artifacts/0/comparison/units/ngam_slope_per_h_per_mmol_atp"),
    ("environment.gem.scope_and_coverage", "gem-environment", 16, 4,
     "/artifacts/0/comparison/units/points_requested"),
    ("kinetics.parsimony.not_equivalence", "carotenoid-kinetics", 3, 3,
     "/artifacts/2/comparison/units/value"),
    ("biosensor.plate_identity.recorded_population", "gates", 4, 4,
     "/artifacts/12/comparison/quantity_units/~1n_sources_processed"),
]


@pytest.fixture(scope="module")
def bounded_evidence_audit():
    return {row["id"]: row for row in audit_claims(_ROOT) if row["record_type"] == "claim"}


@pytest.mark.parametrize("identity,run_id,n_rows,n_groups,unit_pointer", _BOUNDED_SCOPE_CASES)
def test_bounded_evidence_scopes_use_existing_receipts_and_complete_groups(
    bounded_evidence_audit, identity, run_id, n_rows, n_groups, unit_pointer,
):
    row = bounded_evidence_audit[identity]
    assert row["status"] == "PASS" and row["supported"], row["detail"]
    assert row["role"] == "current"
    assert row["provenance"] == {"source": "artifact_registry", "run_id": run_id}
    assert row["provenance_observed"]["status"] == "verified"
    assert row["coverage_observed"]["selected_rows"] == n_rows
    assert row["coverage_observed"]["independent_units"] == n_groups
    assert all(result["unit"] and result["evidence_status"] == "VERIFIED"
               for result in row["result_rows"])
    checks = row["configuration"]["binding"]["checks"]
    assert any(check["field"] == unit_pointer for check in checks)
    assert row["uncertainty"]["status"] not in {"biological_validation", "equivalence"}


def _recorded_scope_fixture(repo, identity):
    # Copy existing observations, never run a scientific producer. Re-pinned fixture
    # identities deliberately isolate the field/coverage guards from byte-integrity guards.
    root, registry = repo
    current = next(claim for claim in load_registry(_ROOT)["claims"] if claim["id"] == identity)
    profile = registry["profiles"]["fixture"]
    for field in ("conditions", "model", "configuration", "uncertainty", "coverage", "independent_unit", "evidence"):
        profile[field] = deepcopy(current[field])
    document = artifacts.load_registry(_ROOT)
    run = document["runs"][current["provenance"]["run_id"]]
    reference = run["verification"]["artifacts"][current["evidence"][0]["path"]]
    documents = {"data/recorded_identity.json": deepcopy(document["receipts"][reference["receipt_id"]]["identity"])}
    for source in current["evidence"]:
        path = _ROOT / source["path"]
        if source["format"] == "json":
            documents[source["path"]] = json.loads(path.read_text())
        else:
            with path.open(newline="", encoding="utf-8") as stream:
                documents[source["path"]] = list(csv.DictReader(stream))
    profile["evidence"].append({"id": "recorded_identity", "path": "data/recorded_identity.json",
                                "format": "json", "pointer": "", "row_keys": [],
                                "unit": "fixture copy of the previously verified run identity"})
    binding = profile["configuration"]["binding"]
    binding.pop("source")
    binding["status"] = "verified"
    for check in binding["checks"]:
        if check.get("source", "receipt") == "receipt":
            check["source"] = "recorded_identity"
    claim = {"id": "fixture.median", "profile": "fixture", "statement": current["statement"],
             "role": current["role"], "evaluation": deepcopy(current["evaluation"])}
    for field in ("acceptance", "expected"):
        if field in current:
            claim[field] = deepcopy(current[field])
    registry["claims"] = [claim]
    return root, registry, documents


def _save_recorded_scope(root, registry, documents):
    for relative, document in documents.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".json":
            path.write_text(json.dumps(document), encoding="utf-8")
        else:
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(document[0]))
                writer.writeheader()
                writer.writerows(document)
    _save(root, registry, pin=True)


def _scope_rows(registry, documents, source_id=None):
    claim = registry["claims"][0]
    source_id = source_id or claim["evaluation"].get("source", "main")
    source = next(source for source in registry["profiles"]["fixture"]["evidence"] if source["id"] == source_id)
    document = documents[source["path"]]
    rows = artifacts.json_pointer(document, source["pointer"]) if source["format"] == "json" else document
    return source, rows if isinstance(rows, list) else [rows]


@pytest.mark.parametrize("identity,run_id,n_rows,n_groups,unit_pointer", _BOUNDED_SCOPE_CASES)
@pytest.mark.parametrize("change", [None, "missing_row", "wrong_identity", "unverified", "wrong_argv", "wrong_units", "stale_output"])
def test_bounded_evidence_rejects_wrong_population_configuration_units_and_provenance(
    repo, identity, run_id, n_rows, n_groups, unit_pointer, change,
):
    root, registry, documents = _recorded_scope_fixture(repo, identity)
    source, rows = _scope_rows(registry, documents)
    where = registry["claims"][0]["evaluation"].get("where", {})
    selected = next(row for row in rows if all(row[field] == value for field, value in where.items()))
    if change == "missing_row":
        rows.remove(selected)
    elif change == "wrong_identity":
        selected[source["row_keys"][0]] = "unregistered-population"
    elif change == "wrong_argv":
        documents["data/recorded_identity.json"]["reproduction"]["argv"][-1] = "--wrong-configuration"
    elif change == "wrong_units":
        parent, _, key = unit_pointer.rpartition("/")
        artifacts.json_pointer(documents["data/recorded_identity.json"], parent)[key.replace("~1", "/")] = "wrong units"
    _save_recorded_scope(root, registry, documents)
    if change == "unverified":
        registry["profiles"]["fixture"]["provenance"] = {"status": "pending", "reason": "No verified execution"}
        _save(root, registry)
    elif change == "stale_output":
        path = root / source["path"]
        path.write_text(path.read_text() + "\n")
    result = _record(root)
    assert result["blocking"] == (change is not None), result["detail"]
    if change is None:
        assert result["coverage_observed"]["selected_rows"] == n_rows
        assert result["coverage_observed"]["independent_units"] == n_groups
    else:
        assert not result["supported"]


@pytest.mark.parametrize("change", ["simulations", "noise", "seeds", "design_count"])
def test_bounded_evidence_power_keeps_simulation_budget_noise_and_seed_identity(repo, change):
    root, registry, documents = _recorded_scope_fixture(repo, "power.synthetic_scope")
    _, rows = _scope_rows(registry, documents)
    attribute = next(row for row in rows if row["question"] == "attribute")
    if change == "simulations":
        attribute["n_simulations"] = "12"
    elif change == "noise":
        attribute["well_cv"] = "0.12"
    elif change == "seeds":
        documents["data/recorded_identity.json"]["parameters"]["candidate_seeds"] = [1, 2]
    else:
        for row in rows:
            if row["design"] == "dose x 3 nutrients":
                row["design"] = "dose x 4 nutrients"
    _save_recorded_scope(root, registry, documents)
    assert _record(root)["blocking"]


@pytest.mark.parametrize("change", ["regime", "reference", "no_growth", "unreachable_penalty", "probe_count"])
def test_bounded_evidence_maintenance_does_not_promote_a_regime_or_infeasible_probe(repo, change):
    root, registry, documents = _recorded_scope_fixture(repo, "maintenance.regime_specific")
    _, rows = _scope_rows(registry, documents)
    transition = next(row for row in rows if row["glucose_lower_bound_mmol_per_gdcw_h"] == "-1.0"
                      and row["oxygen_lower_bound_mmol_per_gdcw_h"] == "-3.7")
    unreachable = next(row for row in rows if row["glucose_lower_bound_mmol_per_gdcw_h"] == "-1.0"
                       and row["oxygen_lower_bound_mmol_per_gdcw_h"] == "-1.0")
    if change == "regime":
        transition["regime_at_top_probe"] = "respiratory"
    elif change == "reference":
        transition["measured_total_ngam_mmol_atp_per_gdw_h"] = "7.2"
    elif change == "no_growth":
        unreachable["no_growth_ngam_mmol_atp_per_gdw_h"] = "0"
    elif change == "unreachable_penalty":
        unreachable["growth_penalty_asserted_percent"] = "100"
    else:
        unreachable["slope_probes_used"] = "5"
    _save_recorded_scope(root, registry, documents)
    assert _record(root)["blocking"]


@pytest.mark.parametrize("change", ["feasible_count", "prediction_status", "headroom", "scoring_supply", "fabricated_infeasible"])
def test_bounded_evidence_gem_cannot_promote_partial_or_infeasible_coverage(repo, change):
    root, registry, documents = _recorded_scope_fixture(repo, "environment.gem.scope_and_coverage")
    _, rows = _scope_rows(registry, documents)
    partial = next(row for row in rows if row["product"] == "phb" and row["growth_rate_per_h"] == "0.1")
    if change == "feasible_count":
        partial["points_feasible"] = "6"
    elif change == "prediction_status":
        partial["is_a_prediction"] = "True"
    elif change == "headroom":
        partial["min_growth_headroom_required"] = "0"
    elif change == "scoring_supply":
        _, scores = _scope_rows(registry, documents, "scores")
        scores[0]["scoring_supply_mmol_c_per_gdcw_h"] = "11"
    else:
        infeasible = next(row for row in rows if row["product"] == "phb" and row["growth_rate_per_h"] == "0.2")
        infeasible["ethanol_over_glucose"] = "1"
    _save_recorded_scope(root, registry, documents)
    assert _record(root)["blocking"]


@pytest.mark.parametrize("change", ["full_params", "strain_count", "p_value", "rss", "degrees_of_freedom"])
def test_bounded_evidence_parsimony_keeps_nested_models_and_nonrejection_not_equivalence(repo, change):
    root, registry, documents = _recorded_scope_fixture(repo, "kinetics.parsimony.not_equivalence")
    _, quantities = _scope_rows(registry, documents, "main")
    by_quantity = {row["quantity"]: row for row in quantities}
    if change == "full_params":
        _, fits = _scope_rows(registry, documents, "fits")
        next(row for row in fits if row["model"] == "capacity per strain x mu")["free params"] = "3"
    elif change == "strain_count":
        documents["data/recorded_identity.json"]["parameters"]["strains"] = ["b-car2", "b-car3"]
    elif change == "p_value":
        by_quantity["nested F p"]["value"] = "0.001"
    elif change == "rss":
        by_quantity["rss full"]["value"] = ""
    else:
        by_quantity["nested F p"]["note"] = "F(2, 20)"
    _save_recorded_scope(root, registry, documents)
    assert _record(root)["blocking"]


@pytest.mark.parametrize("change", ["incomplete", "export_count", "culture_count", "superseded_as_replicate"])
def test_bounded_evidence_plate_population_does_not_count_exports_or_blanks_as_plates(repo, change):
    root, registry, documents = _recorded_scope_fixture(repo, "biosensor.plate_identity.recorded_population")
    manifest = documents["outputs/gates_manifest.json"]
    if change == "incomplete":
        manifest["complete"] = False
    elif change == "export_count":
        manifest["n_sources_processed"] = 4
    elif change == "culture_count":
        next(row for row in manifest["exports"] if row.get("source_set") == "newprotocol")["n_cultures"] = 96
    else:
        next(row for row in manifest["exports"] if row["status"] == "superseded")["source_set"] = "newprotocol"
    _save_recorded_scope(root, registry, documents)
    assert _record(root)["blocking"]


def _read_only_producer(name):
    spec = importlib.util.spec_from_file_location(f"_bounded_evidence_{name}", _ROOT / f"scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bounded_evidence_power_metadata_and_wilson_bounds_match_the_registered_code(repo):
    from ystwin.analysis import power

    script = _read_only_producer("run_power")
    _, registry, documents = _recorded_scope_fixture(repo, "power.synthetic_scope")
    conditions = registry["profiles"]["fixture"]["conditions"]
    assert conditions["doses_mM"] == list(script.DOSES)
    assert conditions["folds"] == list(script.FOLDS)
    assert conditions["simulations"] == script.SIMS == 120
    defaults = inspect.signature(power.discrimination_power).parameters
    assert conditions["reader_cv"] == defaults["reader_cv"].default
    assert conditions["between_plate_log_sd"] == defaults["biological_cv"].default
    assert conditions["between_well_log_sd"] == defaults["well_cv"].default == power.DEFAULT_WELL_CV
    assert power._ALPHA == 0.05
    _, rows = _scope_rows(registry, documents)
    checked = 0
    for row in rows:
        if row["question"] == "detect":
            assert row["n_simulations"] == row["power_n6_low"] == row["power_n6_high"] == ""
            continue
        n = int(float(row["n_simulations"]))
        for prefix in (("power_n6",) if row["question"] == "attribute" else ("power_plate_n4", "power_chemostat_n4")):
            point = float(row[prefix])
            estimate = power.PowerEstimate(round(point * n), n)
            assert estimate.confidence == 0.95
            assert point == pytest.approx(float(estimate))
            assert float(row[f"{prefix}_low"]) == pytest.approx(estimate.low)
            assert float(row[f"{prefix}_high"]) == pytest.approx(estimate.high)
            checked += 1
    assert checked == 20


def test_bounded_evidence_parsimony_degrees_of_freedom_are_not_inferred_from_a_p_value(repo):
    from scipy.stats import f as f_dist
    from ystwin.kinetic.carotenoid import calibration_states

    script = _read_only_producer("fit_carotenoid_kinetics")
    _, registry, documents = _recorded_scope_fixture(repo, "kinetics.parsimony.not_equivalence")
    conditions = registry["profiles"]["fixture"]["conditions"]
    frame = calibration_states()
    assert conditions["states"] == len(frame) == 6
    assert conditions["strains"] == sorted(frame.strain.unique()) == list(script.STRAINS)
    assert frame.groupby("strain").size().tolist() == [2, 2, 2]
    restricted = len(script.CANDIDATES[conditions["restricted_model"]][1])
    full = len(script.CANDIDATES[conditions["full_model"]][1])
    assert conditions["restricted_parameters"] == restricted == 2
    assert conditions["full_parameters"] == full == 4
    df_num, df_den = full - restricted, len(frame) - full
    assert conditions["numerator_degrees_of_freedom"] == df_num == 2
    assert conditions["denominator_degrees_of_freedom"] == df_den == 2
    _, quantities = _scope_rows(registry, documents, "main")
    values = {row["quantity"]: float(row["value"]) for row in quantities}
    statistic = ((values["rss restricted"] - values["rss full"]) / df_num) / (values["rss full"] / df_den)
    assert statistic == pytest.approx(values["nested F statistic"], abs=5e-5)
    assert f_dist.sf(statistic, df_num, df_den) == pytest.approx(values["nested F p"], abs=5e-5)
    assert values["nested F p"] > 0.05
    assert registry["claims"][0]["evaluation"]["kind"] != "equivalence"


def test_bounded_evidence_keeps_unavailable_claims_blocked_without_inferred_refusals(bounded_evidence_audit):
    pending = {
        "biosensor.plate_identity_and_qpcr", "native.biological_law.refused",
        "native.population.retrospective_scope",
        "validation.teacher_vs_changed_equation", "validation.hog.holdout_lineage",
        "validation.prospective.not_yet_scored",
    }
    # This pass owns these claim contracts, not the changing status of unrelated runs.
    # The full read-only audit still reports every other blocker independently.
    assert pending <= bounded_evidence_audit.keys()
    for identity in pending:
        row = bounded_evidence_audit[identity]
        assert row["verdict"] == "BLOCKED" and row["actual"] is None and not row["supported"]
        assert row["role"] != "historical" and row["detail"]
    native = bounded_evidence_audit["native.biological_law.refused"]
    assert native["comparison_status"] == "PASS"
    assert native["expected"]["predicates"] == [
        {"field": "ready_for_biological_law_claim", "value": False},
        {"field": "calibration_adopted", "value": False},
    ]
    real = bounded_evidence_audit["validation.real.rotated_null"]
    assert real["status"] == "PASS" and real["supported"], real["detail"]
    assert real["provenance_observed"]["status"] == "verified"
    assert real["coverage_observed"]["selected_rows"] == 1
    assert real["expected"]["predicates"] == [{"field": "beats_null", "value": False}]
    assert real["evaluation"]["where"] == {"test": "transfer (median own-target score)", "null": "rotated_subspace"}
