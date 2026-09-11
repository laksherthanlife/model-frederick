"""Real development data checks; mutated rows below are engineering fixtures only."""
from __future__ import annotations

from copy import deepcopy
import csv
from decimal import Decimal
import importlib.util
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ystwin.analysis import partial_orders as po
from ystwin.fba.native_reconciliation import native_conditions


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/run_native_order_robustness.py"


@pytest.fixture(scope="module")
def demo():
    spec = importlib.util.spec_from_file_location("native_order_demo", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def conditions():
    return native_conditions()


def by_growth(conditions, growth):
    return next(row for row in conditions if row["growth_rate_per_h"] == growth)


def query(domain, coefficients):
    return po.affine_extrema(
        domain, coefficients, estimand="test exchange contrast", unit=domain.order.scope.unit,
    )


def test_authoritative_tokens_units_and_censoring_are_retained(demo, conditions):
    original = deepcopy(conditions)
    built = demo.build_uptake_domains(conditions, "GlucoseUptake")
    assert len(built["rows"]) == 10
    row = next(row for row in built["rows"] if row["growth_rate_per_h"] == 0.1)
    assert row["condition_id"] == "vanHoek1998:D=0.1"
    assert row["source_column"] == "GlucoseUptake"
    assert row["source_line"] == 4
    assert row["evidence_scope"] == "native_development_training"
    assert row["unit"] == "mmol/gDW/h"
    assert row["primary_printed_token"] == "1.1"
    assert row["rounding_half_width"] == 0.05
    assert row["rounding_is_confidence_interval"] is False
    assert row["rounding_band"] == {"lower": 1.05, "upper": 1.15,
                                     "decimal_lower": "1.05", "decimal_upper": "1.15"}
    assert conditions == original
    scope = built["order"].scope
    assert scope.relation_kind == "magnitude"
    assert scope.quantity == "GlucoseUptake"
    assert "DS28911" in scope.context
    assert conditions[0]["source"]["doi"] in scope.calibration.family
    assert all(bound.source.kind == "assumed" for bound in built["bounded"].bounds)


def test_bands_use_verified_primary_precision_not_float_formatting_or_ci(demo, conditions):
    observations = by_growth(conditions, 0.25)["observations"]
    assert observations["O2uptake"]["primary_printed_token"] == "7.0"
    assert demo.nearest_rounding_band(observations["O2uptake"]) == (
        Decimal("6.95"), Decimal("7.05"),
    )
    record = deepcopy(observations["O2uptake"])
    record["rounding_is_confidence_interval"] = True
    with pytest.raises(ValueError, match="confidence"):
        demo.nearest_rounding_band(record)
    record["rounding_is_confidence_interval"] = False
    record["rounding_half_width"] = 0.5
    with pytest.raises(ValueError, match="precision|rounding"):
        demo.nearest_rounding_band(record)


def test_source_zeros_are_censored_not_zero_anchors_or_numeric_limits(demo, conditions):
    zeros = [record for condition in conditions for record in condition["observations"].values()
             if record["observation_status"] == "below_detection_limit"]
    assert len(zeros) == 26
    for record in zeros:
        assert record["reported_value"] == 0
        assert record["value"] is record["rounding_half_width"] is record["detection_limit"] is None
        assert demo.nearest_rounding_band(record) is None
    # Mutation fixture: a censored source record must stay incomparable in a readout scope.
    changed = deepcopy(conditions)
    target = by_growth(changed, 0.1)
    target["observations"]["GlucoseUptake"] = deepcopy(target["observations"]["Glycerol"])
    built = demo.build_uptake_domains(changed, "GlucoseUptake")
    item = target["condition_id"]
    assert not any(item in pair for pair in built["order"].closure)
    assert not any(bound.item == item for bound in built["bounded"].bounds)
    assert len([pair for pair in built["order"].incomparable_pairs if item in pair]) == 9


@pytest.mark.parametrize("token", ["0.3", "0.4"])
def test_overlapping_or_touching_decimal_bands_remain_incomparable(demo, conditions, token):
    # 0.3 + 0.05 and 0.4 - 0.05 differ in binary floating point, but both equal 0.35.
    changed = deepcopy(conditions)
    record = by_growth(changed, 0.05)["observations"]["GlucoseUptake"]
    record.update(primary_printed_token=token, value=float(token), reported_value=float(token))
    built = demo.build_uptake_domains(changed, "GlucoseUptake")
    low = by_growth(changed, 0.025)["condition_id"]
    high = by_growth(changed, 0.05)["condition_id"]
    assert built["order"].compare(low, high) == "incomparable"
    assert built["order"].ties == ()
    ranks = demo.rank_ranges(built["order"])
    assert ranks[low] == ranks[high] == {"minimum_rank": 1, "maximum_rank": 2}


def test_missing_precision_does_not_reconstruct_bounds_from_the_token(demo, conditions):
    changed = deepcopy(conditions)
    target = by_growth(changed, 0.1)
    target["observations"]["GlucoseUptake"]["rounding_half_width"] = None
    built = demo.build_uptake_domains(changed, "GlucoseUptake")
    item = target["condition_id"]
    assert not any(item in pair for pair in built["order"].closure)
    assert not any(bound.item == item for bound in built["bounded"].bounds)
    assert demo.rank_ranges(built["order"])[item] == {"minimum_rank": 1, "maximum_rank": 10}
    batch = po.sample_uniform_box(
        built["bounded"], n=2, seed=1, budget=4, source=demo.PROPOSAL_SOURCE,
    ).report()
    assert batch["proposal_status"] == "undefined_missing_bounds"
    assert batch["proposal_distribution"] is batch["distribution"] is None
    assert batch["attempted"] == 0


@pytest.mark.parametrize("field,value", [
    ("unit", "mol/gDW/h"), ("observation_status", "unknown"),
    ("value", float("nan")), ("reported_value", 99.0),
])
def test_incompatible_or_malformed_measurements_fail_closed(demo, conditions, field, value):
    conditions[0]["observations"]["GlucoseUptake"][field] = value
    with pytest.raises(ValueError):
        demo.build_uptake_domains(conditions, "GlucoseUptake")


def test_no_cross_quantity_context_or_reserved_scope_pooling(demo, conditions):
    glc = demo.build_uptake_domains(conditions, "GlucoseUptake")
    oxygen = demo.build_uptake_domains(conditions, "O2uptake")
    assert glc["order"].scope != oxygen["order"].scope
    with pytest.raises(ValueError, match="scope|conversion"):
        po.compile_order(glc["order"].items, glc["order"].scope, oxygen["order"].claims)
    with pytest.raises(ValueError, match="uptake"):
        demo.build_uptake_domains(conditions, "BiomassYield_g_per_g")
    changed = deepcopy(conditions)
    # The adapter shares one source dict across conditions; break that alias for a mixed-source fixture.
    changed[0]["source"] = deepcopy(changed[0]["source"])
    changed[0]["source"]["conditions"]["strain"] = "engineering-fixture-other-strain"
    with pytest.raises(ValueError, match="source|context"):
        demo.build_uptake_domains(changed, "GlucoseUptake")
    changed = deepcopy(conditions)
    changed[0]["evidence_scope"] = "reserved_validation"
    with pytest.raises(ValueError, match="development"):
        demo.build_uptake_domains(changed, "GlucoseUptake")


def test_actual_hierarchies_and_rank_computation_do_not_substitute_growth_order(demo, conditions):
    expected_growth_orders = {
        "GlucoseUptake": [0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.35, 0.4],
        "O2uptake": [0.025, 0.05, 0.1, 0.4, 0.15, 0.35, 0.2, 0.3, 0.25, 0.28],
    }
    for key, growth_order in expected_growth_orders.items():
        built = demo.build_uptake_domains(list(reversed(conditions)), key)
        order = built["order"]
        assert len(order.closure) == len(order.claims) == 45
        assert order.incomparable_pairs == ()  # Honest: neither actual readout needs invented unknowns.
        expected_items = [by_growth(conditions, growth)["condition_id"] for growth in growth_order]
        ranks = demo.rank_ranges(order)
        for rank, item in enumerate(expected_items, 1):
            assert ranks[item] == {"minimum_rank": rank, "maximum_rank": rank}
        assert set(map(tuple, demo.cover_edges(order))) == set(zip(expected_items, expected_items[1:]))


def test_order_only_has_no_absolute_scale_and_keeps_semibounded_contrasts(demo, conditions):
    assert demo.PREDECLARED_CONTRASTS == {"GlucoseUptake": (0.4, 0.1), "O2uptake": (0.4, 0.28)}
    for key, (a, b) in demo.PREDECLARED_CONTRASTS.items():
        domain = demo.build_uptake_domains(conditions, key)["order_only"]
        assert domain.bounds == domain.anchors == ()
        assert domain.gap == 0.0
        assert domain.report()["embedding"] == "weak_order_preserving"
        high, low = (by_growth(conditions, growth)["condition_id"] for growth in (a, b))
        absolute = query(domain, {high: 1.0})
        assert absolute["status"] == "unbounded"
        assert absolute["minimum"]["value"] is absolute["maximum"]["value"] is None
        contrast = query(domain, {high: 1.0, low: -1.0})
        assert contrast["status"] == "unbounded"
        finite_side, infinite_side = ("minimum", "maximum") if key == "GlucoseUptake" else ("maximum", "minimum")
        assert contrast[finite_side]["value"] == pytest.approx(0.0)
        assert contrast[infinite_side]["status"] == "unbounded"
        assert contrast[infinite_side]["value"] is contrast[infinite_side]["point"] is None
        assert contrast["sign"] == ("nonnegative_within_tolerance" if key == "GlucoseUptake"
                                     else "nonpositive_within_tolerance")


@pytest.mark.parametrize("key,expected,sign", [
    ("GlucoseUptake", (9.9, 10.1), "positive"), ("O2uptake", (-3.8, -3.6), "negative"),
])
def test_native_bounded_affine_certificates_match_independent_interval_arithmetic(
    demo, conditions, key, expected, sign,
):
    domain = demo.build_uptake_domains(conditions, key)["bounded"]
    a, b = demo.PREDECLARED_CONTRASTS[key]
    high, low = (by_growth(conditions, growth)["condition_id"] for growth in (a, b))
    result = query(domain, {high: 1, low: -1})
    assert result["status"] == "bounded" and result["sign"] == sign
    assert (result["minimum"]["value"], result["maximum"]["value"]) == pytest.approx(expected)
    assert "numerical" in result["certificate"] and "not a rational" in result["certificate"]
    for endpoint in (result["minimum"], result["maximum"], result["feasibility"]):
        assert endpoint["validated"]
        assert domain.contains(endpoint["point"], tolerance=domain.tolerance)
        assert all(0 <= error <= domain.tolerance for error in endpoint["residuals"].values())


@pytest.mark.parametrize("gain", [0.1, 7.0, 1000.0])
def test_positive_shared_scale_covariance_is_explicit_not_identification(demo, conditions, gain):
    for key, growth_pair in demo.PREDECLARED_CONTRASTS.items():
        original = demo.build_uptake_domains(conditions, key)["bounded"]
        scaled = demo.rescale_domain(original, gain)
        assert scaled.order.closure == original.order.closure
        assert scaled.order.incomparable_pairs == original.order.incomparable_pairs
        assert scaled.order.scope.quantity == original.order.scope.quantity
        assert scaled.order.scope != original.order.scope
        assert scaled.order.conversions[0].factor == gain
        a, b = (by_growth(conditions, growth)["condition_id"] for growth in growth_pair)
        base_result, scaled_result = (query(domain, {a: 1, b: -1}) for domain in (original, scaled))
        for side in ("minimum", "maximum"):
            assert scaled_result[side]["value"] == pytest.approx(gain * base_result[side]["value"])
        assert scaled_result["sign"] == base_result["sign"]


@pytest.mark.parametrize("gain", [0, -1, True, float("nan"), float("inf")])
def test_invalid_shared_scales_are_refused(demo, conditions, gain):
    with pytest.raises(ValueError):
        demo.rescale_domain(demo.build_uptake_domains(conditions, "GlucoseUptake")["bounded"], gain)


def test_apparent_yield_has_correct_mmol_to_mass_units_and_inverse_scale(demo):
    assert demo.YIELD_GROWTH_RATE_PER_H == 0.1
    assert demo.GLUCOSE_MOLAR_MASS_G_PER_MOL == pytest.approx(6 * 12.011 + 12 * 1.008 + 6 * 15.999)
    value = demo.apparent_specific_yield(0.1, 1.1)
    assert value == pytest.approx(0.1 / (1.1 * 180.156 / 1000))
    assert 0.5 < value < 0.51
    assert demo.apparent_specific_yield(0.1, 7 * 1.1) == pytest.approx(value / 7)
    for bad in (0, -1, None, float("nan"), float("inf"), True):
        with pytest.raises(ValueError):
            demo.apparent_specific_yield(0.1, bad)


def test_real_report_separates_conditional_lp_and_seeded_sample_envelope(demo):
    report = demo.build_demo(samples=12, seed=1729, budget=12)
    repeated = demo.build_demo(samples=12, seed=1729, budget=12)
    assert report == repeated
    assert report["real_native_data"] and report["law_ready"] is False
    assert report["biological_probability"] is None
    assert report["source_audit"]["n_conditions"] == 10
    assert report["source_audit"]["n_censored"] == 26
    assert report["source_audit"]["evidence_scope_counts"] == {
        "native_development_training": 7, "native_development_interpolation": 3,
    }
    assert report["gem"]["loaded"] is report["gem"]["bounds_modified"] is False
    assumptions = " ".join(report["assumptions"])
    assert "not independent" in assumptions
    assert "not statistical" in assumptions
    for scoped in report["scopes"].values():
        assert scoped["hierarchy"]["n_incomparable_pairs"] == 0
        assert scoped["sampling"]["order_only"]["distribution"] is None
        bounded = scoped["sampling"]["bounded"]
        assert bounded["proposal_status"] == "defined"
        assert bounded["conditioning_status"] == "positive_probability"
        assert bounded["accepted"] == bounded["attempted"] == 12
        assert bounded["biological_posterior"] is False
        assert scoped["scaling"]["verified"] is True
    forward = report["nonlinear"]
    assert forward["condition_id"] == "vanHoek1998:D=0.1"
    assert forward["unit"] == "gDW/g_glucose"
    assert forward["globally_certified"] is forward["is_confidence_interval"] is False
    assert "sampled envelope" in forward["envelope_label"]
    assert forward["accepted"] == forward["successful"] == 12
    assert forward["failed"] == 0
    assert forward["sampled_minimum"] >= demo.apparent_specific_yield(0.1, 1.15)
    assert forward["sampled_maximum"] <= demo.apparent_specific_yield(0.1, 1.05)
    for record in forward["records"]:
        assert record["error"] is record["error_type"] is None
        assert record["value"] == pytest.approx(
            demo.apparent_specific_yield(0.1, record["point"][forward["condition_id"]]),
        )
    json.dumps(report, allow_nan=False)


def test_zero_budget_keeps_definedness_separate_from_sample_success(demo):
    report = demo.build_demo(samples=8, seed=3, budget=0)
    assert report["run_status"] == "incomplete_or_failed"
    forward = report["nonlinear"]
    assert forward["sampled_minimum"] is forward["sampled_maximum"] is None
    assert forward["accepted"] == forward["successful"] == forward["failed"] == 0
    assert forward["sign_stability"] == "inconclusive"
    for scoped in report["scopes"].values():
        batch = scoped["sampling"]["bounded"]
        assert batch["status"] == "budget_exhausted"
        assert batch["proposal_status"] == "defined"
        assert batch["conditioning_status"] == "positive_probability"
        assert batch["distribution"] is not None
    json.dumps(report, allow_nan=False)


def test_callback_failures_remain_nullable_with_coordinates_and_denominators(demo, monkeypatch):
    def fail(*args):
        raise RuntimeError("engineering-only injected callback failure")

    monkeypatch.setattr(demo, "apparent_specific_yield", fail)
    report = demo.build_demo(samples=4, seed=3, budget=4)
    forward = report["nonlinear"]
    assert report["run_status"] == "incomplete_or_failed"
    assert forward["failed"] == forward["accepted"] == 4
    assert forward["successful"] == 0
    assert forward["sampled_minimum"] is forward["sampled_maximum"] is None
    assert forward["sign_stability"] == "inconclusive"
    assert all(row["value"] is None and row["error_type"] == "RuntimeError"
               and row["point"] for row in forward["records"])
    summary = list(csv.DictReader(io.StringIO(demo.summary_csv(report))))
    row = next(row for row in summary if row["query"] == "apparent_specific_yield")
    assert row["failed_samples"] == "4" and row["successful_samples"] == "0"
    assert row["minimum"] == row["maximum"] == ""
    json.dumps(report, allow_nan=False)


def test_solver_failure_keeps_null_endpoints_and_scale_errors_not_fake_zeros(demo, monkeypatch):
    monkeypatch.setattr(po, "linprog", lambda *a, **k: SimpleNamespace(
        status=1, message="engineering-only injected solve limit",
    ))
    report = demo.build_demo(samples=4, seed=3, budget=4)
    assert report["run_status"] == "incomplete_or_failed"
    for scoped in report["scopes"].values():
        for result in scoped["affine"].values():
            assert result["status"] == "numerical_failure"
            assert result["sign"] == "unidentified"
            assert result["minimum"]["value"] is result["maximum"]["value"] is None
        assert scoped["scaling"]["verified"] is None
        assert all(error is None for error in scoped["scaling"]["endpoint_errors"].values())
        assert scoped["sampling"]["bounded"]["distribution"] is not None
    json.dumps(report, allow_nan=False)


def test_orders_and_bands_do_not_count_as_independent_range_evidence(demo, conditions):
    for key, growth_pair in demo.PREDECLARED_CONTRASTS.items():
        bounded = demo.build_uptake_domains(conditions, key)["bounded"]
        no_claims = po.compile_order(bounded.order.items, bounded.order.scope)
        bounds_alone = po.MagnitudeDomain(no_claims, bounded.bounds, gap=0, gap_source=demo.GAP_SOURCE)
        a, b = (by_growth(conditions, growth)["condition_id"] for growth in growth_pair)
        with_order, without_order = (query(domain, {a: 1, b: -1}) for domain in (bounded, bounds_alone))
        for side in ("minimum", "maximum"):
            assert with_order[side]["value"] == pytest.approx(without_order[side]["value"])


def test_runner_opens_only_the_authoritative_native_table(demo, conditions, monkeypatch):
    expected = (ROOT / conditions[0]["source"]["table_relative_path"]).resolve()
    original_open = Path.open
    opened = []

    def checked_open(path, *args, **kwargs):
        assert path.resolve() == expected, f"unexpected data access: {path}"
        opened.append(path.resolve())
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", checked_open)
    report = demo.build_demo(samples=2, seed=3, budget=2)
    assert opened == [expected]
    assert report["gem"]["loaded"] is report["gem"]["bounds_modified"] is False


def test_cli_budget_exhaustion_writes_nullable_results_and_returns_failure(demo, tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "__file__", str(tmp_path / "scripts" / SCRIPT.name))
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    assert demo.main(["--output-dir", str(outputs), "--budget", "0"]) == 1
    report = json.loads((outputs / demo.REPORT_NAME).read_text())
    assert report["run_status"] == "incomplete_or_failed"
    assert report["nonlinear"]["sampled_minimum"] is report["nonlinear"]["sampled_maximum"] is None
    with (outputs / demo.SUMMARY_NAME).open(newline="") as stream:
        row = next(row for row in csv.DictReader(stream) if row["query"] == "apparent_specific_yield")
    assert row["minimum_status"] == row["maximum_status"] == "unavailable_no_successful_samples"
    assert row["proposal_status"] == "defined"
    assert row["conditioning_status"] == "positive_probability"
    assert row["sampling_status"] == "budget_exhausted"


def test_fresh_outputs_include_audit_columns_and_refuse_overwrite_before_loading(demo, tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "__file__", str(tmp_path / "scripts" / SCRIPT.name))
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    assert demo.main(["--output-dir", str(outputs), "--samples", "6", "--budget", "6"]) == 0
    paths = [outputs / name for name in (demo.REPORT_NAME, demo.SUMMARY_NAME, demo.OBSERVATIONS_NAME)]
    assert all(path.name.startswith("native_order_") for path in paths)
    before = {path: path.read_bytes() for path in paths}
    report = json.loads(paths[0].read_text())
    assert report["run_status"] == "complete"
    with paths[1].open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 9
    for row in rows:
        assert row["source_doi"] == report["source"]["doi"]
        assert row["n_incomparable_pairs"] == "0"
        assert row["assumptions"] and row["guarantee"] and row["proposal_status"]
        assert row["failed_samples"] != ""
        assert row["details_file"] == demo.REPORT_NAME
        assert row["estimand"]
        if row["query"] == "bounded_contrast":
            key = row["observable_id"]
            assert json.loads(row["coefficients"]) == report["scopes"][key]["predeclared_contrast"]["coefficients"]
            assert row["minimum_validated"] == row["maximum_validated"] == "True"
            assert float(row["maximum_certificate_residual"]) <= 1e-8
        if row["query"] == "apparent_specific_yield":
            assert row["fixed_condition_id"] == "vanHoek1998:D=0.1"
            assert row["minimum_validated"] == row["maximum_validated"] == "False"
            assert row["maximum_certificate_residual"] == ""
        if row["sampling_status"] != "not_requested":
            assert row["requested_proposal"] and row["proposal_source"]
            assert row["rejected_samples"] == "0"
    with paths[2].open(newline="") as stream:
        observations = list(csv.DictReader(stream))
    assert len(observations) == 20
    assert all(row["rounding_is_confidence_interval"] == "False" for row in observations)

    def must_not_load():
        raise AssertionError("must reject output collision before reading any data")

    monkeypatch.setattr(demo, "native_conditions", must_not_load)
    with pytest.raises(SystemExit):
        demo.main(["--output-dir", str(outputs)])
    assert {path: path.read_bytes() for path in paths} == before
    with pytest.raises(SystemExit):
        demo.main(["--output-dir", str(tmp_path)])
