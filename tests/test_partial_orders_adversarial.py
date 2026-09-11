"""Independent, synthetic engineering checks; these fixtures authorize no biology.

All sources, contexts, calibrations, ranges, gaps and sampling proposals below are
caller-supplied assumptions. Even a fixture labelled ``measured`` tests a label's
round trip, not independently verified evidence. No external data are opened.
"""
from __future__ import annotations

import json
from dataclasses import asdict, replace
from itertools import permutations, product
from math import fsum, isfinite
from types import SimpleNamespace

import numpy as np
import pytest

from ystwin.analysis import partial_orders as po


ASSUMPTION = po.Source(
    "assumed", "engineering-only:adversarial-fixture", "constructed algebra, not an assay"
)
SYNTHETIC = po.Source(
    "synthetic", "engineering-only:synthetic-fixture", "deterministic test construction"
)
SCOPE = po.Scope(
    "magnitude", "synthetic-context", "abstract-capacity", "fixture-unit",
    po.Calibration("synthetic-family", "known"),
)


def _domain(items, edges=(), *, limits=None, anchors=None, ties=(), gap=0.0, scope=SCOPE):
    order = po.compile_order(
        items, scope,
        [po.OrderClaim(left, right, scope, ASSUMPTION) for left, right in edges],
        ties=[po.DeclaredTie(tuple(group), scope, ASSUMPTION) for group in ties],
    )
    return po.MagnitudeDomain(
        order,
        bounds=tuple(po.Bound(name, lower, upper, scope, ASSUMPTION)
                     for name, (lower, upper) in (limits or {}).items()),
        anchors=tuple(po.Anchor(name, value, scope, ASSUMPTION)
                      for name, value in (anchors or {}).items()),
        gap=gap, gap_source=ASSUMPTION,
    )


def _query(domain, coefficients, offset=0.0):
    return po.affine_extrema(
        domain, coefficients, estimand="engineering affine query", unit="fixture-output",
        offset=offset,
    )


def _explicit_residual(domain, point):
    """Check public input rows directly, without the implementation's LP/membership code."""
    assert set(point) == set(domain.order.items)
    assert all(isfinite(value) for value in point.values())
    errors = [0.0]
    conversions = (*domain.order.conversions, *domain.conversions)

    def convert(value, scope):
        if scope == domain.order.scope:
            return value
        rows = [row for row in conversions if row.origin == scope]
        assert len(rows) == 1
        return rows[0].factor * value + rows[0].offset

    for row in domain.bounds:
        if row.lower is not None:
            errors.append(convert(row.lower, row.scope) - point[row.item])
        if row.upper is not None:
            errors.append(point[row.item] - convert(row.upper, row.scope))
    for row in domain.anchors:
        errors.append(abs(point[row.item] - convert(row.value, row.scope)))
    for row in domain.order.ties:
        errors.extend(abs(point[name] - point[row.members[0]]) for name in row.members)
    for row in domain.order.claims:
        errors.append(domain.gap - (point[row.upper] - point[row.lower]))
    return max(errors)


def _assert_endpoint(domain, endpoint, coefficients, offset=0.0):
    assert endpoint["status"] == "optimal", endpoint
    assert endpoint["validated"] is True
    point = endpoint["point"]
    residual = _explicit_residual(domain, point)
    assert residual <= domain.tolerance
    value = fsum(coefficients.get(name, 0.0) * point[name] for name in point) + offset
    assert endpoint["value"] == pytest.approx(value, abs=domain.tolerance, rel=1e-10)
    assert endpoint["residuals"]["primal_absolute"] == pytest.approx(residual, abs=1e-12)
    assert set(endpoint["residuals"]) == {
        "primal_absolute", "dual_sign_absolute", "absent_bound_dual_absolute",
        "stationarity_relative", "duality_gap_relative", "objective_relative",
        "complementarity_relative",
    }
    assert all(isfinite(value) and 0 <= value <= domain.tolerance
               for value in endpoint["residuals"].values())


def _assert_extrema_at_corners(domain, coefficients, corners, offset=0.0):
    expected = [fsum(coefficients.get(name, 0.0) * value for name, value in point.items())
                + offset for point in corners]
    result = _query(domain, coefficients, offset)
    assert result["status"] == "bounded", result
    assert result["minimum"]["value"] == pytest.approx(min(expected), abs=1e-8)
    assert result["maximum"]["value"] == pytest.approx(max(expected), abs=1e-8)
    for endpoint in (result["minimum"], result["maximum"]):
        _assert_endpoint(domain, endpoint, coefficients, offset)
    _assert_endpoint(domain, result["feasibility"], {})
    json.dumps(result, allow_nan=False)
    return result


def _scope_from_dict(row):
    return po.Scope(**{**row, "calibration": po.Calibration(**row["calibration"])})


def _conversion_from_dict(row):
    return po.Conversion(
        _scope_from_dict(row["origin"]), _scope_from_dict(row["target"]),
        row["factor"], row["offset"], po.Source(**row["source"]),
    )


def _domain_from_report(report):
    """Reconstruct public constructors from their JSON rows, not a new module API."""
    order_report = report["order"]
    scope = _scope_from_dict(order_report["scope"])
    claims = [po.OrderClaim(row["lower"], row["upper"], _scope_from_dict(row["scope"]),
                            po.Source(**row["source"])) for row in order_report["claims"]]
    ties = [po.DeclaredTie(tuple(row["members"]), _scope_from_dict(row["scope"]),
                           po.Source(**row["source"])) for row in order_report["declared_ties"]]
    order = po.compile_order(
        order_report["items"], scope, claims, ties=ties,
        conversions=[_conversion_from_dict(row) for row in order_report["conversions"]],
    )
    return po.MagnitudeDomain(
        order,
        bounds=tuple(po.Bound(row["item"], row["lower"], row["upper"],
                              _scope_from_dict(row["scope"]), po.Source(**row["source"]))
                     for row in report["bounds"]),
        anchors=tuple(po.Anchor(row["item"], row["value"], _scope_from_dict(row["scope"]),
                                po.Source(**row["source"])) for row in report["anchors"]),
        conversions=tuple(_conversion_from_dict(row)
                          for row in report["additional_conversions"]),
        gap=report["gap"]["value"], gap_source=po.Source(**report["gap"]["source"]),
        tolerance=report["lp_tolerance"],
    )


@pytest.mark.parametrize("kind", ["temporal_precedence", "causal_dependency", "effect", "importance"])
def test_magnitude_claims_cannot_be_relabelled_as_other_relations(kind):
    target = replace(SCOPE, relation_kind=kind,
                     estimand="explicit synthetic estimand" if kind in {"effect", "importance"} else None)
    claim = po.OrderClaim("a", "b", SCOPE, ASSUMPTION)
    with pytest.raises(ValueError):
        po.compile_order(("a", "b"), target, [claim])
    with pytest.raises(ValueError):
        po.Conversion(SCOPE, target, 1.0, 0.0, ASSUMPTION)


@pytest.mark.parametrize("kind", ["temporal_precedence", "causal_dependency"])
def test_nonnumeric_precedence_cannot_become_a_magnitude_domain(kind):
    scope = replace(SCOPE, relation_kind=kind)
    order = po.compile_order(("a", "b"), scope, [po.OrderClaim("a", "b", scope, ASSUMPTION)])
    with pytest.raises(ValueError):
        po.MagnitudeDomain(order, gap=0.0, gap_source=ASSUMPTION)


@pytest.mark.parametrize("kind", ["effect", "importance"])
def test_importance_and_effect_require_a_named_estimand(kind):
    with pytest.raises(ValueError):
        replace(SCOPE, relation_kind=kind)
    with pytest.raises(ValueError):
        replace(SCOPE, relation_kind=kind, estimand=" ")


@pytest.mark.parametrize("kind", ["magnitude", "temporal_precedence", "causal_dependency"])
def test_feedback_cycle_is_rejected_with_actual_cycle_and_sources(kind):
    scope = replace(SCOPE, relation_kind=kind)
    claims = tuple(po.OrderClaim(a, b, scope, po.Source("assumed", f"fixture-edge-{i}", "constructed"))
                   for i, (a, b) in enumerate((("a", "b"), ("b", "c"), ("c", "a"))))
    with pytest.raises(po.CycleError) as caught:
        po.compile_order(("a", "b", "c", "unconnected"), scope, claims)
    error = caught.value
    assert error.cycle[0] == error.cycle[-1]
    assert set(error.claims) == set(claims)
    assert {group[0] for group in error.cycle} == {"a", "b", "c"}
    for claim in claims:
        assert claim.source.reference in str(error)
    if kind == "causal_dependency":
        assert "feedback" in str(error)


def test_strict_edge_inside_declared_tie_is_a_contradiction_not_a_sort():
    with pytest.raises(po.CycleError) as caught:
        po.compile_order(
            ("a", "b", "c"), SCOPE, [po.OrderClaim("a", "c", SCOPE, ASSUMPTION)],
            ties=[po.DeclaredTie(("a", "b"), SCOPE, ASSUMPTION),
                  po.DeclaredTie(("b", "c"), SCOPE, ASSUMPTION)],
        )
    assert caught.value.cycle == (("a", "b", "c"), ("a", "b", "c"))


def test_unknown_and_disconnected_items_do_not_become_ties_or_total_ranks():
    order = po.compile_order(
        ("isolated", "c", "b", "a"), SCOPE,
        [po.OrderClaim("a", "b", SCOPE, ASSUMPTION)],
        ties=[po.DeclaredTie(("b", "c"), SCOPE, ASSUMPTION)],
    )
    assert order.compare("a", "c") == "less"
    assert order.compare("c", "a") == "greater"
    assert order.compare("b", "c") == "tied"
    assert order.compare("a", "a") == "same"
    assert order.compare("a", "isolated") == "incomparable"
    assert order.compare("b", "isolated") == "incomparable"
    assert order.incomparable_pairs == (("a", "isolated"), ("b", "isolated"), ("c", "isolated"))
    with pytest.raises(ValueError):
        order.compare("a", "absent")


def test_absent_evidence_never_invents_an_order_even_with_equal_numeric_anchors():
    domain = _domain(("z", "a", "m"), anchors={"z": 2.0, "a": 2.0})
    assert domain.order.claims == ()
    assert domain.order.closure == ()
    assert len(domain.order.incomparable_pairs) == 3
    assert domain.order.compare("z", "a") == "incomparable"
    assert domain.order.report()["source_kinds"] == []
    assert domain.feasibility()["status"] == "optimal"
    assert domain.order.compare("z", "a") == "incomparable"


@pytest.mark.parametrize("change", [
    {"unit": "incompatible-unit"},
    {"context": "other-synthetic-context"},
    {"quantity": "other-quantity"},
    {"calibration": po.Calibration("other-family", "known")},
    {"calibration": po.Calibration("synthetic-family", "shared_unknown_positive")},
])
def test_unsupported_comparisons_and_range_transfers_are_refused(change):
    other = replace(SCOPE, **change)
    with pytest.raises(ValueError):
        po.compile_order(("a", "b"), SCOPE, [po.OrderClaim("a", "b", other, ASSUMPTION)])
    with pytest.raises(ValueError):
        po.compile_order(("a", "b"), SCOPE, ties=[po.DeclaredTie(("a", "b"), other, ASSUMPTION)])
    order = po.compile_order(("a",), SCOPE)
    with pytest.raises(ValueError):
        po.MagnitudeDomain(order, bounds=(po.Bound("a", 0, 1, other, ASSUMPTION),),
                           gap=0, gap_source=ASSUMPTION)
    with pytest.raises(ValueError):
        po.MagnitudeDomain(order, anchors=(po.Anchor("a", 0, other, ASSUMPTION),),
                           gap=0, gap_source=ASSUMPTION)


@pytest.mark.parametrize("factor", [0.0, -1.0, float("inf"), float("nan"), True])
def test_nonpositive_or_nonfinite_conversions_cannot_preserve_order(factor):
    with pytest.raises(ValueError):
        po.Conversion(replace(SCOPE, unit="input-unit"), SCOPE, factor, 0.0, ASSUMPTION)


def test_explicit_affine_conversion_transforms_bounds_and_anchors_consistently():
    origin = replace(SCOPE, unit="input-unit", context="input-fixture-context")
    conversion = po.Conversion(origin, SCOPE, 2.0, 3.0, ASSUMPTION)
    order = po.compile_order(
        ("a", "b"), SCOPE, [po.OrderClaim("a", "b", origin, ASSUMPTION)],
        conversions=[conversion],
    )
    domain = po.MagnitudeDomain(
        order, bounds=(po.Bound("b", 2.0, 4.0, origin, ASSUMPTION),),
        anchors=(po.Anchor("a", 1.0, origin, ASSUMPTION),),
        gap=1.0, gap_source=ASSUMPTION,
    )
    result = _assert_extrema_at_corners(
        domain, {"a": 3.0, "b": -2.0}, [{"a": 5.0, "b": 7.0}, {"a": 5.0, "b": 11.0}],
    )
    assert result["domain"]["order"]["conversions"][0]["source"] == asdict(ASSUMPTION)


@pytest.mark.parametrize("gain", [1e-6, 0.25, 7.0, 1e6])
def test_shared_positive_gain_preserves_only_compatible_observed_order(gain):
    scope = replace(SCOPE, calibration=po.Calibration("shared-fixture-gain", "shared_unknown_positive"))
    domain = _domain(("a", "b", "c"), (("a", "b"), ("b", "c")), scope=scope)
    point = {"a": 2.0 * gain, "b": 3.0 * gain, "c": 5.0 * gain}
    assert domain.contains(point)
    assert domain.order.compare("a", "c") == "less"
    result = _query(domain, {"b": 1.0})
    assert result["status"] == "unbounded"
    assert result["minimum"]["value"] is None
    assert result["maximum"]["value"] is None
    assert result["domain"]["bounds"] == []
    assert result["domain"]["anchors"] == []


@pytest.mark.parametrize("factor_quantity", ["E", "kcat"])
def test_product_order_does_not_identify_e_or_kcat_separately(factor_quantity):
    product_scope = replace(SCOPE, quantity="E_times_kcat")
    factor_scope = replace(product_scope, quantity=factor_quantity)
    order = po.compile_order(
        ("a", "b"), product_scope, [po.OrderClaim("a", "b", product_scope, ASSUMPTION)],
    )
    with pytest.raises(ValueError):
        po.Conversion(product_scope, factor_scope, 1.0, 0.0, ASSUMPTION)
    # Same observed products, opposite orders of either factor. This is algebra,
    # not an assertion that these synthetic values describe biological parameters.
    factors_one = {"a": (1.0, 6.0), "b": (4.0, 3.0)}
    factors_two = {"a": (6.0, 1.0), "b": (3.0, 4.0)}
    domain = po.MagnitudeDomain(order, gap=0, gap_source=ASSUMPTION)
    for factors in (factors_one, factors_two):
        assert domain.contains({name: e * kcat for name, (e, kcat) in factors.items()})
    column = 0 if factor_quantity == "E" else 1
    assert ((factors_one["a"][column] < factors_one["b"][column])
            != (factors_two["a"][column] < factors_two["b"][column]))


def test_bounds_anchors_and_gap_require_explicit_provenance():
    order = po.compile_order(("a",), SCOPE)
    with pytest.raises(TypeError):
        po.Bound("a", 0, 1, SCOPE)
    with pytest.raises(TypeError):
        po.Anchor("a", 0, SCOPE)
    with pytest.raises(TypeError):
        po.MagnitudeDomain(order)
    with pytest.raises(TypeError):
        po.MagnitudeDomain(order, gap=0)
    with pytest.raises(TypeError):
        po.MagnitudeDomain(order, gap=0, gap_source=None)
    with pytest.raises(ValueError):
        po.Bound("a", None, None, SCOPE, ASSUMPTION)


def test_empty_item_set_is_explicitly_refused_instead_of_fabricated():
    with pytest.raises(ValueError):
        po.compile_order((), SCOPE)


def test_order_only_domain_has_no_default_scale_prior_or_finite_interval():
    domain = _domain(("a", "b"), (("a", "b"),))
    assert domain.contains({"a": -1e9, "b": -1e8})
    assert domain.contains({"a": 1e8, "b": 1e9})
    assert all(row["lower"] is None and row["upper"] is None
               for row in domain.report()["effective_box"])
    result = _query(domain, {"a": 2.0, "b": 1.0})
    assert result["status"] == "unbounded"
    for name in ("minimum", "maximum"):
        assert result[name]["status"] == "unbounded"
        assert result[name]["value"] is None
        assert result[name]["point"] is None
        assert result[name]["validated"] is False
    assert domain.report()["probability_distribution"] == "none unless an explicit sampling proposal is requested"
    json.dumps(result, allow_nan=False)


def test_requesting_uniform_samples_cannot_define_a_distribution_on_an_unbounded_box():
    domain = _domain(("a", "b"), (("a", "b"),))
    batch = po.sample_uniform_box(domain, n=5, seed=4, budget=50, source=ASSUMPTION)
    report = batch.report()
    assert report["status"] == "unidentified"
    assert report["accepted"] == report["attempted"] == report["rejected"] == 0
    assert report["points"] == []
    assert report["acceptance_rate"] is None
    assert report["biological_posterior"] is False
    # There is no uniform probability distribution on an unbounded box. A
    # requested proposal may be documented separately, but not reported as one.
    assert report["distribution"] is None


def test_closed_weak_embedding_does_not_claim_order_reflection_or_a_strict_gap():
    domain = _domain(("a", "b", "c"), (("a", "b"),), anchors={"a": 2, "b": 2, "c": 2})
    assert domain.contains({"a": 2.0, "b": 2.0, "c": 2.0})
    assert domain.order.compare("a", "b") == "less"
    assert domain.order.compare("a", "c") == "incomparable"
    assert domain.report()["embedding"] == "weak_order_preserving"
    assert domain.report()["order_reflecting"] is False
    assert domain.report()["positive_gap_exceeds_tolerance"] is False


@pytest.mark.parametrize("coefficients", [
    {"a": 2.0, "b": -1.0}, {"a": -2.0, "b": 1.0},
    {"a": 1.0, "b": 1.0}, {"a": -1.0, "b": -1.0}, {},
])
def test_triangle_lp_minimum_and_maximum_match_three_known_corners(coefficients):
    domain = _domain(("a", "b"), (("a", "b"),), limits={"a": (0, 2), "b": (0, 2)})
    corners = [{"a": 0.0, "b": 0.0}, {"a": 0.0, "b": 2.0}, {"a": 2.0, "b": 2.0}]
    _assert_extrema_at_corners(domain, coefficients, corners, offset=1.25)


@pytest.mark.parametrize("coefficients", [
    {"a": -3.0, "b": 2.0, "c": -4.0}, {"a": 3.0, "b": -2.0, "c": 4.0},
    {"a": 0.0, "b": 1.0, "c": -1.0}, {"a": 1.0}, {"b": -1.0},
])
def test_branched_order_polytope_matches_independently_enumerated_binary_corners(coefficients):
    domain = _domain(
        ("a", "b", "c"), (("a", "b"), ("a", "c")),
        limits={name: (0, 1) for name in ("a", "b", "c")},
    )
    corners = [dict(zip(("a", "b", "c"), values)) for values in product((0.0, 1.0), repeat=3)
               if values[0] <= values[1] and values[0] <= values[2]]
    assert len(corners) == 5
    _assert_extrema_at_corners(domain, coefficients, corners, offset=-0.75)


@pytest.mark.parametrize("coefficients", [
    {"a": 3.0, "b": -2.0, "c": 1.0}, {"a": -3.0, "b": 2.0, "c": -1.0}, {},
])
def test_fixed_and_disconnected_nodes_are_not_dropped_from_affine_objectives(coefficients):
    domain = _domain(
        ("a", "b", "c"), (("a", "b"),), limits={"b": (0, 4), "c": (-1, 3)}, anchors={"a": 2},
    )
    corners = [{"a": 2.0, "b": b, "c": c} for b, c in product((2.0, 4.0), (-1.0, 3.0))]
    _assert_extrema_at_corners(domain, coefficients, corners)


def test_positive_gap_accumulates_along_paths_and_extrema_are_not_corner_only_guesses():
    domain = _domain(
        ("a", "b", "c"), (("a", "b"), ("b", "c")), gap=0.25,
        limits={name: (0, 1) for name in ("a", "b", "c")},
    )
    # Subtract (0, .25, .5): a weak chain in [0, .5], whose four vertices are known.
    corners = [dict(zip(("a", "b", "c"), values)) for values in
               ((0.0, 0.25, 0.5), (0.0, 0.25, 1.0), (0.0, 0.75, 1.0), (0.5, 0.75, 1.0))]
    result = _assert_extrema_at_corners(domain, {"a": -2.0, "b": 3.0, "c": -1.0}, corners)
    assert result["minimum"]["value"] == pytest.approx(0.75 - 1.0)
    assert _query(domain, {"c": 1, "a": -1})["minimum"]["value"] == pytest.approx(0.5)


def test_tied_coefficients_are_summed_and_intersected_bounds_are_respected():
    domain = _domain(
        ("a", "b"), ties=(("a", "b"),), limits={"a": (-1, 2), "b": (1, 3)},
    )
    _assert_extrema_at_corners(domain, {"a": 4, "b": -1}, [{"a": 1, "b": 1}, {"a": 2, "b": 2}])
    result = _query(domain, {"a": 1, "b": -1}, offset=7)
    assert result["minimum"]["value"] == result["maximum"]["value"] == 7


@pytest.mark.parametrize("offset", [-3.0, 0.0, 5.0])
def test_constant_objective_is_bounded_even_when_parameter_domain_is_unbounded(offset):
    domain = _domain(("a", "unconnected"))
    result = _query(domain, {}, offset=offset)
    assert result["status"] == "bounded"
    assert result["minimum"]["value"] == result["maximum"]["value"] == offset
    for name in ("minimum", "maximum"):
        _assert_endpoint(domain, result[name], {}, offset)
    assert all(row["lower"] is None and row["upper"] is None
               for row in result["domain"]["effective_box"])


@pytest.mark.parametrize("case", ["tie-conflict", "anchor-conflict", "order-conflict", "gap-conflict"])
def test_empty_feasible_domains_are_infeasible_even_for_constant_objectives(case):
    if case == "tie-conflict":
        domain = _domain(("a", "b"), ties=(("a", "b"),), limits={"a": (0, 1), "b": (2, 3)})
    elif case == "anchor-conflict":
        domain = _domain(("a",), limits={"a": (0, 1)}, anchors={"a": 2})
    elif case == "order-conflict":
        domain = _domain(("a", "b"), (("a", "b"),), anchors={"a": 2, "b": 1})
    else:
        domain = _domain(("a", "b"), (("a", "b"),), anchors={"a": 1, "b": 1}, gap=0.25)
    result = _query(domain, {}, offset=13)
    assert result["status"] == "infeasible"
    assert result["sign"] == "unidentified"
    for name in ("minimum", "maximum", "feasibility"):
        assert result[name]["status"] == "infeasible"
        assert result[name]["value"] is None
        assert result[name]["point"] is None
        assert result[name]["validated"] is False
    batch = po.sample_uniform_box(domain, n=3, seed=4, budget=20, source=ASSUMPTION)
    assert batch.status == "infeasible"
    assert batch.attempted == 0
    assert batch.points == ()


def test_disconnected_unbounded_nuisance_does_not_spoil_an_identified_query():
    domain = _domain(("a", "nuisance"), limits={"a": (-2, 3)})
    result = _query(domain, {"a": 2.0}, offset=1)
    assert result["status"] == "bounded"
    assert result["minimum"]["value"] == -3
    assert result["maximum"]["value"] == 7
    for endpoint in (result["minimum"], result["maximum"]):
        _assert_endpoint(domain, endpoint, {"a": 2}, 1)
    assert _query(domain, {"nuisance": 1})["status"] == "unbounded"


@pytest.mark.parametrize("bounds,coefficients,known_endpoint,expected", [
    ((2.0, None), {"a": 1.0}, "minimum", 2.0),
    ((None, -2.0), {"a": 1.0}, "maximum", -2.0),
    ((2.0, None), {"a": -1.0}, "maximum", -2.0),
])
def test_one_sided_unbounded_queries_preserve_the_finite_endpoint(bounds, coefficients, known_endpoint, expected):
    domain = _domain(("a",), limits={"a": bounds})
    result = _query(domain, coefficients)
    assert result["status"] == "unbounded"
    assert result[known_endpoint]["value"] == expected
    _assert_endpoint(domain, result[known_endpoint], coefficients)
    other = "maximum" if known_endpoint == "minimum" else "minimum"
    assert result[other]["status"] == "unbounded"
    assert result[other]["value"] is None


@pytest.mark.parametrize("point,expected", [
    ({"a": 0.5, "b": 0.5, "c": 1.0}, 0.0),
    ({"a": 0.5, "b": 0.6, "c": 1.0}, 0.1),
    ({"a": 0.5, "b": 0.5, "c": 0.6}, 0.15),
    ({"a": -0.2, "b": -0.2, "c": 1.0}, 0.2),
])
def test_membership_residual_is_checked_independently_for_ties_bounds_and_gaps(point, expected):
    domain = _domain(
        ("a", "b", "c"), (("b", "c"),), ties=(("a", "b"),), gap=0.25,
        limits={"a": (0, 1), "b": (0, 1), "c": (0, 2)},
    )
    assert domain.violation(point) == pytest.approx(expected)
    assert _explicit_residual(domain, point) == pytest.approx(expected)
    assert domain.contains(point) is (expected == 0)


@pytest.mark.parametrize("corruption", ["primal", "objective", "dual-sign", "stationarity", "nan"])
def test_a_success_flag_from_solver_is_not_itself_a_valid_witness(monkeypatch, corruption):
    domain = _domain(("a",), limits={"a": (0, 1)})
    result = SimpleNamespace(
        status=0, message="synthetic success flag; not trusted", x=np.array([0.0]), fun=0.0,
        ineqlin=SimpleNamespace(marginals=np.array([])),
        lower=SimpleNamespace(marginals=np.array([0.0])),
        upper=SimpleNamespace(marginals=np.array([0.0])),
    )
    if corruption == "primal":
        result.x[0] = -1.0
    elif corruption == "objective":
        result.fun = 2.0
    elif corruption == "dual-sign":
        result.lower.marginals[0] = -1.0
        result.upper.marginals[0] = 1.0
    elif corruption == "stationarity":
        result.lower.marginals[0] = 1.0
    else:
        result.x[0] = np.nan
    monkeypatch.setattr(po, "linprog", lambda *args, **kwargs: result)
    feasibility = domain.feasibility()
    assert feasibility["status"] == "numerical_failure"
    assert feasibility["validated"] is False
    assert feasibility["value"] is None
    assert feasibility["point"] is None


def test_sampling_budget_records_rejections_and_every_accepted_point_is_strictly_admissible():
    domain = _domain(("a", "b"), (("a", "b"),), gap=0.25, limits={"a": (0, 1), "b": (0, 1)})
    batch = po.sample_uniform_box(domain, n=100, seed=27, budget=50, source=ASSUMPTION)
    report = batch.report()
    assert report["status"] == "budget_exhausted"
    assert report["attempted"] == report["budget"] == 50
    assert 0 < report["accepted"] < report["attempted"]
    assert report["rejected"] == report["attempted"] - report["accepted"]
    assert report["acceptance_rate"] == report["accepted"] / report["attempted"]
    assert report["membership_tolerance"] == 0.0
    assert report["proposal_source"] == asdict(ASSUMPTION)
    assert report["biological_posterior"] is False
    for point in report["points"]:
        assert 0 <= point["a"] <= 1
        assert 0 <= point["b"] <= 1
        assert point["b"] - point["a"] >= 0.25
        assert _explicit_residual(domain, point) == 0.0
    assert batch.points == po.sample_uniform_box(domain, n=100, seed=27, budget=50, source=ASSUMPTION).points


def test_zero_budget_and_zero_measure_slices_never_get_a_biased_fallback_sample():
    domain = _domain(("a", "b"), (("a", "b"),), limits={"a": (0, 1)}, anchors={"b": 0})
    for budget in (0, 29):
        batch = po.sample_uniform_box(domain, n=2, seed=101, budget=budget, source=ASSUMPTION)
        assert batch.status == "budget_exhausted"
        assert batch.points == ()
        assert batch.attempted == budget
        assert batch.report()["rejected"] == budget
        assert batch.report()["feasibility"]["status"] == "optimal"


def test_inferred_finite_bounds_are_not_silently_promoted_to_a_sampling_prior():
    domain = _domain(("a", "b", "c"), (("a", "b"), ("b", "c")), anchors={"a": 0, "c": 1})
    assert _query(domain, {"b": 1})["status"] == "bounded"
    batch = po.sample_uniform_box(domain, n=2, seed=9, budget=30, source=ASSUMPTION)
    assert batch.status == "unidentified"
    assert batch.points == ()
    assert batch.attempted == 0
    assert domain.bounds == ()


def test_declared_ties_and_fixed_coordinates_are_sampled_in_the_quotient_not_by_near_equal_rejection():
    domain = _domain(
        ("a", "b", "c"), (("a", "c"),), ties=(("a", "b"),),
        limits={"a": (0, 0.5), "b": (0.1, 0.6)}, anchors={"c": 1}, gap=0.25,
    )
    batch = po.sample_uniform_box(domain, n=23, seed=2, budget=23, source=ASSUMPTION)
    assert batch.status == "complete"
    assert batch.attempted == 23
    assert batch.report()["rejected"] == 0
    for point in batch.report()["points"]:
        assert point["a"] == point["b"]
        assert 0.1 <= point["a"] <= 0.5
        assert point["c"] == 1.0
        assert _explicit_residual(domain, point) == 0.0


@pytest.mark.parametrize("n,seed,budget", [
    (0, 0, 1), (True, 0, 1), (1.5, 0, 1), (1, -1, 1), (1, True, 1), (1, 0, -1), (1, 0, 2.5),
])
def test_finite_sampling_budget_is_explicit_and_integer(n, seed, budget):
    domain = _domain(("a",), limits={"a": (0, 1)})
    with pytest.raises(ValueError):
        po.sample_uniform_box(domain, n=n, seed=seed, budget=budget, source=ASSUMPTION)


def test_failed_forward_calls_are_counted_and_cannot_support_sign_stability():
    domain = _domain(("a",), limits={"a": (0, 1)})
    batch = po.sample_uniform_box(domain, n=40, seed=8, budget=40, source=ASSUMPTION)
    seen = []

    def forward(point):
        seen.append(point["a"])
        if point["a"] < 0.5:
            raise RuntimeError("deliberate engineering forward failure")
        return point["a"] + 1.0

    result = po.evaluate_forward(batch, forward, estimand="synthetic response", unit="fixture-output")
    expected_failures = sum(point[0] < 0.5 for point in batch.points)
    assert 0 < expected_failures < len(batch.points)
    assert len(seen) == len(batch.points) == result["accepted"] == 40
    assert result["failed"] == expected_failures
    assert result["successful"] == 40 - expected_failures
    assert result["sign_stability"] == "inconclusive"
    assert result["globally_certified"] is False
    assert [row["sample_index"] for row in result["records"]] == list(range(40))
    for record, point in zip(result["records"], batch.report()["points"]):
        assert record["point"] == point
        if point["a"] < 0.5:
            assert record["status"] == "failed"
            assert record["value"] is None
            assert record["error_type"] == "RuntimeError"
        else:
            assert record["status"] == "ok"
            assert record["value"] == point["a"] + 1
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf"), -float("inf"), True, "1", None])
def test_nonfinite_or_non_numeric_forward_outputs_are_failures_not_dropped(bad_value):
    batch = po.sample_uniform_box(_domain(("a",), limits={"a": (0, 1)}), n=3, seed=7, budget=3, source=ASSUMPTION)
    result = po.evaluate_forward(batch, lambda point: bad_value, estimand="bad-output test", unit="fixture-output")
    assert result["accepted"] == result["failed"] == 3
    assert result["successful"] == 0
    assert result["sampled_minimum"] is result["sampled_maximum"] is None
    assert result["sign_stability"] == "inconclusive"
    assert result["successful_sample_sign"] == "unidentified"
    assert result["globally_certified"] is False
    assert all(row["status"] == "failed" and row["value"] is None for row in result["records"])
    json.dumps(result, allow_nan=False)


def test_forward_cannot_mutate_accepted_samples():
    domain = _domain(("a",), limits={"a": (0, 1)})
    batch = po.sample_uniform_box(domain, n=3, seed=1, budget=3, source=ASSUMPTION)
    before = batch.points

    def forward(point):
        point["a"] = 99
        return 1.0

    result = po.evaluate_forward(batch, forward, estimand="mutation test", unit="fixture-output")
    assert result["failed"] == 3
    assert all(row["error_type"] == "TypeError" for row in result["records"])
    assert batch.points == before


@pytest.mark.parametrize("bad_point", [(-1.0,), (float("nan"),), (), (0.5, 0.5)])
def test_manually_supplied_invalid_sample_never_reaches_forward_callback(bad_point):
    domain = _domain(("a",), limits={"a": (0, 1)})
    batch = po.sample_uniform_box(domain, n=1, seed=1, budget=1, source=ASSUMPTION)
    batch = replace(batch, points=(bad_point,))
    called = []
    result = po.evaluate_forward(
        batch, lambda point: called.append(point) or 1.0,
        estimand="invalid sample test", unit="fixture-output",
    )
    assert called == []
    assert result["failed"] == 1
    assert result["successful"] == 0
    assert result["globally_certified"] is False


def test_adversarial_nonlinear_interior_peak_is_not_a_global_sample_certificate():
    domain = _domain(("x",), limits={"x": (0, 1)})
    batch = po.sample_uniform_box(domain, n=11, seed=314, budget=11, source=ASSUMPTION)
    coordinates = sorted([0.0, 1.0, *(point[0] for point in batch.points)])
    left, right = max(zip(coordinates, coordinates[1:]), key=lambda pair: pair[1] - pair[0])
    center, radius = (left + right) / 2.0, (right - left) / 4.0
    assert 0 < center < 1 and radius > 0

    def interior_peak(point):
        u = (point["x"] - center) / radius
        return -1.0 + 2.0 * max(0.0, 1.0 - u * u) ** 2

    assert domain.contains({"x": center})
    assert interior_peak({"x": center}) == 1.0
    assert interior_peak({"x": 0.0}) == interior_peak({"x": 1.0}) == -1.0
    result = po.evaluate_forward(batch, interior_peak, estimand="constructed interior peak", unit="fixture-output")
    assert result["successful"] == 11
    assert result["failed"] == 0
    assert result["sampled_minimum"] == result["sampled_maximum"] == -1.0
    assert result["sampled_maximum"] < interior_peak({"x": center})
    assert result["successful_sample_sign"] == "negative"
    assert result["sign_stability"] == "negative_on_samples"
    assert result["globally_certified"] is False


@pytest.mark.parametrize("names", list(permutations(("z", "middle", "a"))))
def test_parameter_renaming_and_coefficient_insertion_order_do_not_change_affine_extrema(names):
    original = ("left", "right", "isolated")
    renamed = dict(zip(original, names))
    coefficients = {"left": 3.0, "right": -2.0, "isolated": 4.0}
    limits = {"left": (-1, 2), "right": (0, 3), "isolated": (-2, 1)}
    baseline = _query(_domain(original, (("left", "right"),), limits=limits, gap=0.5), coefficients, 2.0)
    domain = _domain(
        reversed(names), ((renamed["left"], renamed["right"]),), gap=0.5,
        limits={renamed[name]: limits[name] for name in reversed(original)},
    )
    transformed = {renamed[name]: coefficients[name] for name in reversed(original)}
    result = _query(domain, transformed, 2.0)
    assert result["status"] == baseline["status"] == "bounded"
    assert result["sign"] == baseline["sign"]
    for endpoint in ("minimum", "maximum"):
        assert result[endpoint]["value"] == pytest.approx(baseline[endpoint]["value"], abs=1e-8)
        _assert_endpoint(domain, result[endpoint], transformed, 2.0)
    assert domain.order.compare(renamed["left"], renamed["isolated"]) == "incomparable"


def test_serialization_round_trip_preserves_relations_limits_sources_and_seeded_samples():
    origin = replace(SCOPE, unit="input-unit")
    conversion = po.Conversion(origin, SCOPE, 2.0, 1.0, ASSUMPTION)
    order = po.compile_order(
        ("c", "b", "a"), SCOPE, [po.OrderClaim("a", "c", origin, SYNTHETIC)],
        ties=[po.DeclaredTie(("a", "b"), SCOPE, ASSUMPTION)], conversions=[conversion],
    )
    domain = po.MagnitudeDomain(
        order, bounds=(po.Bound("a", 0, 1, origin, ASSUMPTION),
                       po.Bound("b", 1, 4, SCOPE, ASSUMPTION)),
        anchors=(po.Anchor("c", 5, SCOPE, ASSUMPTION),), gap=0.25, gap_source=ASSUMPTION,
    )
    serialized = json.loads(json.dumps(domain.report(), allow_nan=False, sort_keys=True))
    restored = _domain_from_report(serialized)
    assert restored.report() == domain.report()
    coefficients = {"a": 2, "b": -1, "c": -1}
    first, second = _query(domain, coefficients), _query(restored, coefficients)
    for endpoint in ("minimum", "maximum"):
        assert first[endpoint]["value"] == second[endpoint]["value"]
        _assert_endpoint(restored, second[endpoint], coefficients)
    left = po.sample_uniform_box(domain, n=10, seed=991, budget=30, source=SYNTHETIC)
    right = po.sample_uniform_box(restored, n=10, seed=991, budget=30, source=SYNTHETIC)
    assert left.report() == right.report()


@pytest.mark.parametrize("label", ["assumed", "synthetic", "simulated", "measured"])
def test_source_classifications_remain_caller_supplied_labels_not_biological_authorization(label):
    source = po.Source(label, "unverified-fixture-label", "not independently verified")
    order = po.compile_order(("a", "b"), SCOPE, [po.OrderClaim("a", "b", SCOPE, source)])
    assert order.report()["source_kinds"] == [label]
    assert order.report()["claims"][0]["source"] == asdict(source)
    domain = po.MagnitudeDomain(
        order, bounds=tuple(po.Bound(name, 0, 1, SCOPE, source) for name in ("a", "b")),
        gap=0, gap_source=source,
    )
    batch = po.sample_uniform_box(domain, n=2, seed=19, budget=20, source=source)
    assert batch.report()["biological_posterior"] is False
    result = po.evaluate_forward(batch, lambda point: point["a"] ** 2,
                                 estimand="engineering label test", unit="fixture-output")
    assert result["globally_certified"] is False
    assert result["sampling"]["proposal_source"] == asdict(source)
