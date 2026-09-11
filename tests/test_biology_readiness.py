from __future__ import annotations

import pytest

from scripts import check_biology_readiness as readiness
from ystwin import paths
from ystwin.fba.native_reconciliation import ec_native_observables, native_conditions
from ystwin.fba.solver import load_model


@pytest.fixture(scope="module")
def native_host():
    model, _ = load_model(paths.ec_yeast_gem())
    return model, ec_native_observables(model)


def test_native_readiness_keeps_reported_growth_fixed_while_testing_printed_precision(native_host):
    model, observables = native_host
    condition = next(item for item in native_conditions() if item["growth_rate_per_h"] == 0.1)
    exact = readiness._condition_probe(model, condition, observables, "quantified_point", 15)
    rounded = readiness._condition_probe(model, condition, observables, "printed_rounding", 15)
    assert exact["consistent"] is False
    assert rounded["consistent"] is True
    assert rounded["growth_rate_per_h"] == exact["growth_rate_per_h"] == 0.1
    assert rounded["uncertainty_scope"] == "printed rounding, not a confidence interval"
    assert rounded["physical_audit"]["max_abs_mass_balance"] <= 1e-7


def test_infeasible_conditions_are_not_removed_from_readiness(native_host):
    model, observables = native_host
    condition = next(item for item in native_conditions() if item["growth_rate_per_h"] == 0.4)
    result = readiness._condition_probe(model, condition, observables, "printed_rounding", 15)
    assert result["condition_id"] == condition["condition_id"]
    assert result["consistent"] is False
    assert result["status"] == "infeasible"
    assert result["physical_audit"] is None


def test_gate_refuses_unknown_native_constraint_mode_before_solving(native_host):
    model, observables = native_host
    with pytest.raises(ValueError, match="constraint mode"):
        readiness._condition_probe(model, native_conditions()[0], observables, "auto_relax", 15)


def test_readiness_revalidates_collector_claims_and_keeps_failed_conditions(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(readiness, "collect_learning_evidence", lambda root: {"evidence_predicates": {}})

    def require(ledger, claim, *, root):
        calls.append((claim, root))
        if claim == "independent_mechanistic_evidence":
            raise ValueError("no authoritative evidence assessor")
        return {"authorized": True}

    monkeypatch.setattr(readiness, "require_claim", require)
    monkeypatch.setattr(readiness, "load_model", lambda path: (object(), "fixture"))
    monkeypatch.setattr(readiness, "ec_native_observables", lambda model: {})
    monkeypatch.setattr(readiness, "native_conditions", lambda: [True, False])
    monkeypatch.setattr(readiness, "_condition_probe", lambda model, condition, *args: {"consistent": condition})
    report = readiness.run_readiness(tmp_path, tmp_path / "report")
    assert report["n_native_conditions"] == 2
    assert report["n_native_consistent"] == 1
    assert report["ready_for_biological_law_claim"] is False
    assert report["all_frozen_files_unchanged"] is True
    assert [claim for claim, _ in calls].count("integrity_established") == 2
    assert all(root == tmp_path for _, root in calls)
    assert (tmp_path / "report" / "readiness.json").is_file()
