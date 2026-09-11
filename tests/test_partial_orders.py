from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from dataclasses import replace
from itertools import combinations, permutations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ystwin.analysis import partial_orders as po


def test_structured_forward_solver_status_does_not_turn_infeasible_into_zero():
    admissible = domain()
    batch = sample(admissible, n=2)
    responses = iter((po.ForwardEvaluation("optimal", 0.0),
                      po.ForwardEvaluation("infeasible", None, {"solver": "toy"})))
    result = po.evaluate_forward(batch, lambda _: next(responses),
                                 estimand="growth", unit="1/h")
    assert result["successful"] == 1
    assert result["failed"] == 1
    assert result["records"][0]["value"] == 0.0
    assert result["records"][1]["status"] == "infeasible"
    assert result["records"][1]["value"] is None
    assert result["sign_stability"] == "inconclusive"


SOURCE = po.Source("synthetic", "neutral-test-fixture", "constructed, not biological data")
SCOPE = po.Scope(
    "magnitude", "neutral-context", "quantity", "arbitrary_unit",
    po.Calibration("neutral-family", "known"),
)


def claim(lower, upper, scope=SCOPE, source=SOURCE):
    return po.OrderClaim(lower, upper, scope, source)


def bound(item, lower, upper, scope=SCOPE, source=SOURCE):
    return po.Bound(item, lower, upper, scope, source)


def order(items=("a", "b", "c"), edges=(("a", "b"),), scope=SCOPE, **kwargs):
    return po.compile_order(items, scope, [claim(a, b, scope) for a, b in edges], **kwargs)


def domain(poset=None, lower=0.0, upper=1.0, gap=0.0, **kwargs):
    poset = order() if poset is None else poset
    return po.MagnitudeDomain(
        poset, [bound(item, lower, upper, poset.scope) for item in poset.items],
        gap=gap, gap_source=SOURCE, **kwargs,
    )


def extrema(admissible, coefficients, offset=0.0):
    return po.affine_extrema(
        admissible, coefficients, estimand="neutral affine query",
        unit="query_unit", offset=offset,
    )


def sample(admissible, n=64, seed=3, budget=10000):
    return po.sample_uniform_box(admissible, n=n, seed=seed, budget=budget, source=SOURCE)


def test_closure_does_not_totalize_unknown_pairs():
    poset = order(("d", "c", "b", "a"), (("a", "b"), ("b", "c")))
    assert poset.items == ("a", "b", "c", "d")
    assert set(poset.closure) == {("a", "b"), ("a", "c"), ("b", "c")}
    assert set(poset.incomparable_pairs) == {("a", "d"), ("b", "d"), ("c", "d")}
    assert poset.compare("a", "a") == "same"
    assert poset.compare("a", "c") == "less"
    assert poset.compare("c", "a") == "greater"
    assert poset.compare("a", "d") == "incomparable"
    with pytest.raises(ValueError, match="unknown"):
        poset.compare("a", "missing")


@pytest.mark.parametrize("edges", [(("a", "a"),), (("a", "b"), ("b", "a")),
                                  (("a", "b"), ("b", "c"), ("c", "a"))])
def test_strict_cycles_are_diagnostic_and_never_ties(edges):
    with pytest.raises(po.CycleError) as caught:
        order(edges=edges)
    error = caught.value
    assert error.cycle[0] == error.cycle[-1]
    assert len(error.claims) == len(edges)
    assert "neutral-test-fixture" in str(error)
    assert "cycle" in str(error).lower()


def test_overlapping_declared_ties_form_a_quotient_only_when_explicit():
    ties = [po.DeclaredTie(("a", "b"), SCOPE, SOURCE),
            po.DeclaredTie(("b", "c"), SCOPE, SOURCE)]
    poset = order(("a", "b", "c", "d"), (("a", "d"),), ties=ties)
    assert poset.classes == (("a", "b", "c"), ("d",))
    assert poset.compare("a", "c") == "tied"
    assert poset.compare("b", "d") == "less"
    assert set(poset.closure) == {("a", "d"), ("b", "d"), ("c", "d")}
    with pytest.raises(po.CycleError):
        order(edges=(("a", "b"),), ties=ties)


def test_noisy_opposite_claims_are_not_merged_even_if_sources_differ():
    other = po.Source("measured", "neutral-assay-replicate", "neutral measurement fixture")
    with pytest.raises(po.CycleError) as caught:
        po.compile_order(("a", "b"), SCOPE, [claim("a", "b"), claim("b", "a", source=other)])
    assert {row.source.kind for row in caught.value.claims} == {"synthetic", "measured"}


@pytest.mark.parametrize("kind", ["temporal_precedence", "causal_dependency"])
def test_non_magnitude_posets_cannot_be_used_as_magnitude_domains(kind):
    scope = replace(SCOPE, relation_kind=kind)
    poset = order(scope=scope)
    assert poset.compare("a", "b") == "less"
    with pytest.raises(ValueError, match="magnitude"):
        domain(poset)
    with pytest.raises(ValueError, match="relation kind"):
        po.Conversion(scope, SCOPE, 1.0, 0.0, SOURCE)


def test_causal_feedback_is_rejected_as_non_poset():
    scope = replace(SCOPE, relation_kind="causal_dependency")
    with pytest.raises(po.CycleError, match="feedback"):
        order(scope=scope, edges=(("a", "b"), ("b", "a")))


@pytest.mark.parametrize("kind", ["importance", "effect"])
def test_importance_and_effect_require_a_named_estimand(kind):
    with pytest.raises(ValueError, match="estimand"):
        replace(SCOPE, relation_kind=kind)
    scoped = replace(SCOPE, relation_kind=kind, estimand="change in output per unit intervention")
    assert order(scope=scoped).scope.estimand == scoped.estimand


@pytest.mark.parametrize("change", [
    {"context": "other-context"}, {"unit": "other-unit"},
    {"quantity": "other-quantity"},
    {"calibration": po.Calibration("other-family", "known")},
    {"calibration": po.Calibration("neutral-family", "shared_unknown_positive")},
])
def test_mixed_scope_claims_bounds_and_ties_need_conversion_evidence(change):
    other = replace(SCOPE, **change)
    with pytest.raises(ValueError, match="conversion"):
        po.compile_order(("a", "b", "c"), SCOPE,
                         [claim("a", "b"), claim("b", "c", other)])
    with pytest.raises(ValueError, match="conversion"):
        po.compile_order(("a", "b"), SCOPE, ties=[po.DeclaredTie(("a", "b"), other, SOURCE)])
    with pytest.raises(ValueError, match="conversion"):
        po.MagnitudeDomain(order(), [bound("a", 0, 1, other)], gap=0, gap_source=SOURCE)
    with pytest.raises(ValueError, match="conversion"):
        po.MagnitudeDomain(order(), anchors=[po.Anchor("a", 1, other, SOURCE)],
                           gap=0, gap_source=SOURCE)


def test_explicit_positive_affine_conversion_keeps_evidence_and_converts_bounds():
    origin = replace(SCOPE, context="other-assay-context", unit="other_unit",
                     calibration=po.Calibration("other-family", "known"))
    conversion = po.Conversion(origin, SCOPE, 2.0, 1.0, SOURCE)
    poset = po.compile_order(("a", "b"), SCOPE, [claim("a", "b", origin)],
                             conversions=[conversion])
    admissible = po.MagnitudeDomain(
        poset, [bound("a", 0, 1, origin), bound("b", 0, 1, origin)],
        gap=0, gap_source=SOURCE,
    )
    report = extrema(admissible, {"a": 1})
    assert report["status"] == "bounded"
    assert report["minimum"]["value"] == pytest.approx(1)
    assert report["maximum"]["value"] == pytest.approx(3)
    assert report["domain"]["order"]["conversions"][0]["source"]["reference"] == SOURCE.reference


@pytest.mark.parametrize("factor", [0, -1, math.inf, math.nan])
def test_conversion_must_be_finite_and_increasing(factor):
    with pytest.raises(ValueError):
        po.Conversion(replace(SCOPE, unit="other"), SCOPE, factor, 0.0, SOURCE)


def test_conversions_cannot_relabel_quantities_or_effect_estimands():
    with pytest.raises(ValueError, match="quantity"):
        po.Conversion(replace(SCOPE, quantity="unrelated"), SCOPE, 1.0, 0.0, SOURCE)
    scoped = replace(SCOPE, relation_kind="effect", estimand="effect on output X")
    with pytest.raises(ValueError, match="estimand"):
        po.Conversion(replace(scoped, estimand="effect on output Y"), scoped, 1, 0, SOURCE)


@pytest.mark.parametrize("factor,offset", [(2, 0), (1, 2)])
def test_a_nonidentity_conversion_cannot_silently_redefine_its_own_scope(factor, offset):
    with pytest.raises(ValueError, match="distinct"):
        po.Conversion(SCOPE, SCOPE, factor, offset, SOURCE)


def test_ambiguous_conversion_origins_are_rejected_not_chosen_by_row_order():
    origin = replace(SCOPE, unit="other_unit")
    first = po.Conversion(origin, SCOPE, 2, 0, SOURCE)
    second = po.Conversion(origin, SCOPE, 3, 0, SOURCE)
    with pytest.raises(ValueError, match="ambiguous"):
        po.compile_order(("a", "b"), SCOPE, [claim("a", "b", origin)],
                         conversions=[first, second])


def test_shared_unknown_positive_gain_preserves_absolute_quantity_order_without_bounds():
    scope = replace(SCOPE, unit="molecules_per_item",
                    calibration=po.Calibration("shared-assay-family", "shared_unknown_positive"))
    poset = order(scope=scope)
    assert poset.compare("a", "b") == "less"
    admissible = po.MagnitudeDomain(poset, gap=0, gap_source=SOURCE)
    absolute = extrema(admissible, {"a": 1})
    assert absolute["status"] == "unbounded"
    assert absolute["minimum"]["status"] == absolute["maximum"]["status"] == "unbounded"
    assert absolute["minimum"]["value"] is None
    contrast = extrema(admissible, {"b": 1, "a": -1})
    assert contrast["minimum"]["value"] == pytest.approx(0)
    assert contrast["maximum"]["status"] == "unbounded"
    assert contrast["sign"] == "nonnegative_within_tolerance"
    assert sample(admissible).status == "unidentified"


def test_simulation_derived_orders_keep_their_origin_not_measured_status():
    simulated = po.Source("simulated", "synthetic-sensitivity-run", "finite difference on toy map")
    poset = po.compile_order(("a", "b"), SCOPE, [claim("a", "b", source=simulated)])
    report = poset.report()
    assert report["claims"][0]["source"]["kind"] == "simulated"
    assert report["source_kinds"] == ["simulated"]


@pytest.mark.parametrize("factory", [
    lambda: po.Source("", "id", "assay"),
    lambda: po.Source("posterior", "id", "assay"),
    lambda: po.Source("measured", "", "assay"),
    lambda: po.Source("measured", "id", ""),
    lambda: po.Calibration("", "known"),
    lambda: po.Calibration("family", "unknown"),
    lambda: replace(SCOPE, unit=""),
    lambda: replace(SCOPE, context=""),
    lambda: replace(SCOPE, relation_kind="unknown"),
    lambda: po.compile_order(("a", "a"), SCOPE),
    lambda: po.compile_order((), SCOPE),
    lambda: po.compile_order(("a",), SCOPE, [claim("a", "b")]),
    lambda: po.DeclaredTie(("a",), SCOPE, SOURCE),
    lambda: bound("a", None, None),
    lambda: bound("a", math.nan, 1),
    lambda: bound("a", 0, math.inf),
    lambda: bound("a", 2, 1),
    lambda: po.Anchor("a", math.inf, SCOPE, SOURCE),
    lambda: domain(gap=-1),
    lambda: domain(gap=math.nan),
    lambda: domain(tolerance=0),
    lambda: domain(tolerance=1e-12),
])
def test_invalid_or_unsourced_inputs_fail_early(factory):
    with pytest.raises((ValueError, TypeError)):
        factory()


def test_missing_bounds_do_not_trigger_scipy_default_nonnegativity():
    admissible = po.MagnitudeDomain(order(("a",), ()), gap=0, gap_source=SOURCE)
    assert admissible.feasibility()["status"] == "optimal"
    assert admissible.contains({"a": -123.0})
    report = extrema(admissible, {"a": 1})
    assert report["minimum"]["status"] == "unbounded"
    assert report["maximum"]["status"] == "unbounded"
    assert extrema(admissible, {}, offset=2)["minimum"]["value"] == pytest.approx(2)


def test_weak_embedding_does_not_claim_strict_numeric_separation():
    poset = order(("a", "b"), (("a", "b"),))
    weak = domain(poset)
    assert weak.contains({"a": 0.5, "b": 0.5})
    assert weak.report()["embedding"] == "weak_order_preserving"
    assert poset.compare("a", "b") == "less"
    separated = domain(poset, gap=0.25)
    assert not separated.contains({"a": 0.5, "b": 0.5})
    assert separated.report()["embedding"] == "positive_gap_order_preserving"
    assert extrema(separated, {"b": 1, "a": -1})["minimum"]["value"] == pytest.approx(0.25)


def test_numeric_equalities_never_rewrite_qualitative_incomparability():
    poset = order(("a", "b"), ())
    admissible = domain(poset, lower=2, upper=2)
    assert admissible.contains({"a": 2, "b": 2})
    assert poset.compare("a", "b") == "incomparable"
    assert poset.classes == (("a",), ("b",))


def test_anchors_and_orders_can_bound_queries_but_do_not_invent_a_proposal_box():
    poset = order(("a", "b", "c"), (("a", "b"), ("b", "c")))
    admissible = po.MagnitudeDomain(
        poset, anchors=[po.Anchor("a", 0, SCOPE, SOURCE), po.Anchor("c", 2, SCOPE, SOURCE)],
        gap=0.25, gap_source=SOURCE,
    )
    report = extrema(admissible, {"b": 1})
    assert report["minimum"]["value"] == pytest.approx(0.25)
    assert report["maximum"]["value"] == pytest.approx(1.75)
    assert sample(admissible).status == "unidentified"
    assert len(report["domain"]["anchors"]) == 2


def test_affine_extrema_are_sharp_with_primal_dual_residuals():
    admissible = domain()
    robust = extrema(admissible, {"b": 1, "a": -1}, offset=0.125)
    variable = extrema(admissible, {"c": 1}, offset=-0.5)
    assert robust["minimum"]["value"] == pytest.approx(0.125)
    assert robust["maximum"]["value"] == pytest.approx(1.125)
    assert robust["sign"] == "positive"
    assert variable["minimum"]["value"] == pytest.approx(-0.5)
    assert variable["maximum"]["value"] == pytest.approx(0.5)
    assert variable["sign"] == "crosses_zero"
    assert "numerical" in robust["certificate"] and "rational" in robust["certificate"]
    assert robust["tolerance"] == admissible.tolerance
    for endpoint in (robust["minimum"], robust["maximum"]):
        assert admissible.contains(endpoint["point"], tolerance=admissible.tolerance)
        assert all(value <= admissible.tolerance for value in endpoint["residuals"].values())
        assert endpoint["validated"]
    json.dumps(robust, allow_nan=False)


def test_ties_intersect_bounds_and_sum_affine_coefficients():
    poset = order(ties=[po.DeclaredTie(("a", "c"), SCOPE, SOURCE)])
    admissible = po.MagnitudeDomain(
        poset, [bound("a", 0, 2), bound("c", 1, 3), bound("b", 0, 4)],
        gap=0.5, gap_source=SOURCE,
    )
    tied_difference = extrema(admissible, {"a": 1, "c": -1})
    assert tied_difference["minimum"]["value"] == pytest.approx(0)
    assert tied_difference["maximum"]["value"] == pytest.approx(0)
    assert extrema(admissible, {"a": 1})["minimum"]["value"] == pytest.approx(1)
    assert not admissible.contains({"a": 1, "c": 2, "b": 3})


@pytest.mark.parametrize("admissible", [
    lambda: domain(gap=2),
    lambda: po.MagnitudeDomain(order(), [bound("a", 2, 3), bound("b", 0, 1)],
                              gap=0, gap_source=SOURCE),
    lambda: po.MagnitudeDomain(order(), [bound("a", 0, 1)],
                              anchors=[po.Anchor("a", 2, SCOPE, SOURCE)],
                              gap=0, gap_source=SOURCE),
    lambda: po.MagnitudeDomain(order(ties=[po.DeclaredTie(("a", "c"), SCOPE, SOURCE)]),
                              [bound("a", 0, 1), bound("c", 2, 3)],
                              gap=0, gap_source=SOURCE),
])
def test_infeasibility_is_not_reported_as_robust_or_sampled(admissible):
    admissible = admissible()
    assert admissible.feasibility()["status"] == "infeasible"
    report = extrema(admissible, {"a": 1})
    assert report["status"] == "infeasible"
    assert report["sign"] == "unidentified"
    assert report["minimum"]["value"] is None
    batch = sample(admissible)
    assert batch.status == "infeasible"
    assert not batch.points


@pytest.mark.parametrize("target", ["x", "fun", "dual"])
def test_solver_success_is_not_enough_without_residual_validation(monkeypatch, target):
    real = po.linprog

    def corrupt(*args, **kwargs):
        result = real(*args, **kwargs)
        if result.status == 0:
            if target == "x":
                result.x = np.full_like(result.x, 99.0)
            elif target == "fun":
                result.fun = 123.0
            else:
                result.lower.marginals = np.full_like(result.x, -20.0)
        return result

    monkeypatch.setattr(po, "linprog", corrupt)
    result = extrema(domain(), {"a": 1})
    assert result["status"] == "numerical_failure"
    assert result["sign"] == "unidentified"
    assert result["minimum"]["value"] is None


@pytest.mark.parametrize("solver_status", [1, 3, 4])
def test_structured_solver_failure_is_not_misread_as_feasibility(monkeypatch, solver_status):
    monkeypatch.setattr(po, "linprog", lambda *a, **k: SimpleNamespace(
        status=solver_status, message="failed zero-objective feasibility solve",
    ))
    report = extrema(domain(), {"a": 1})
    assert report["status"] == "numerical_failure"
    assert report["sign"] == "unidentified"
    assert sample(domain()).status == "numerical_failure"


def test_standalone_feasibility_carries_tolerance_and_numerical_certificate_label():
    report = domain().feasibility()
    assert report["tolerance"] == 1e-8
    assert "numerical" in report["certificate"]
    assert "rational" in report["certificate"]


def test_generator_conversion_inputs_do_not_disappear_from_range_provenance():
    origin = replace(SCOPE, unit="other_unit")
    conversion = po.Conversion(origin, SCOPE, 2.0, 0.0, SOURCE)
    admissible = po.MagnitudeDomain(
        order(), [bound("a", 0, 1, origin)],
        conversions=iter([conversion]), gap=0, gap_source=SOURCE,
    )
    assert len(admissible.report()["additional_conversions"]) == 1
    assert admissible.report()["effective_box"][0]["upper"] == 2


def test_extreme_or_overflowed_inputs_never_become_numeric_certificates():
    with pytest.raises(ValueError, match="finite"):
        bound("a", 0, 10 ** 1000)
    tied = order(ties=[po.DeclaredTie(("a", "c"), SCOPE, SOURCE)])
    with pytest.raises(ValueError, match="coefficient"):
        extrema(domain(tied), {"a": 1e308, "c": 1e308})


def test_forward_revalidates_supplied_sample_coordinates_and_keeps_invalid_rows():
    batch = sample(domain(), n=2)
    corrupt = replace(batch, points=((100.0, 0.0, 0.0), batch.points[1]))
    calls = []

    def forward(point):
        calls.append(point)
        return 1.0

    report = po.evaluate_forward(corrupt, forward, estimand="validation fixture", unit="toy")
    assert report["accepted"] == 2
    assert report["failed"] == 1
    assert len(calls) == 1
    assert report["records"][0]["point"]["a"] == 100
    assert "outside" in report["records"][0]["error"]


@pytest.mark.parametrize("coefficients", [{"missing": 1}, {"a": math.inf}, {"b": math.nan}])
def test_invalid_affine_queries_are_rejected(coefficients):
    with pytest.raises(ValueError):
        extrema(domain(), coefficients)


def test_finite_positive_gain_and_shift_transform_sharp_extrema_equivariantly():
    rng = np.random.default_rng(152)
    base = domain(order(("a", "b"), (("a", "b"),)), lower=-2, upper=3, gap=0.1)
    base_result = extrema(base, {"a": 2, "b": -1})
    for gain in (0.01, 0.1, 1.0, 13.0, 100.0):
        shift = float(rng.uniform(-4, 4))
        transformed = domain(base.order, lower=-2 * gain + shift, upper=3 * gain + shift,
                             gap=base.gap * gain)
        result = extrema(transformed, {"a": 2, "b": -1})
        for endpoint in ("minimum", "maximum"):
            assert result[endpoint]["value"] == pytest.approx(
                gain * base_result[endpoint]["value"] + shift,
            )


def test_all_small_dags_have_transitive_antisymmetric_closure_and_keep_unknowns():
    names = ("a", "b", "c", "d")
    possible = tuple(combinations(names, 2))
    for mask in range(1 << len(possible)):
        edges = [edge for i, edge in enumerate(possible) if mask & (1 << i)]
        poset = order(names, edges)
        closure = set(poset.closure)
        expected = set(edges)
        for middle in names:
            expected |= {(a, b) for a in names for b in names
                         if (a, middle) in expected and (middle, b) in expected}
        assert closure == expected
        assert all((b, a) not in closure for a, b in closure)
        for a, b in combinations(names, 2):
            assert ((a, b) in poset.incomparable_pairs) == (
                (a, b) not in closure and (b, a) not in closure
            )


def test_order_input_permutations_do_not_choose_an_arbitrary_total_order():
    edges = [("a", "c"), ("b", "c")]
    baseline = order(("a", "b", "c"), edges)
    for names in permutations(baseline.items):
        for rows in permutations(edges):
            candidate = order(names, rows)
            assert candidate.closure == baseline.closure
            assert candidate.incomparable_pairs == baseline.incomparable_pairs
            assert candidate.classes == baseline.classes


def test_lp_extrema_match_enumerated_vertices_of_small_order_polytope():
    poset = order(("a", "b", "c", "d"), (("a", "c"), ("b", "c")))
    admissible = domain(poset)
    vertices = [dict(zip(poset.items, coordinates)) for coordinates in np.ndindex(2, 2, 2, 2)]
    vertices = [point for point in vertices if admissible.contains(point)]
    rng = np.random.default_rng(37)
    for _ in range(20):
        coefficients = dict(zip(poset.items, rng.normal(size=4)))
        values = [sum(coefficients[name] * point[name] for name in poset.items) for point in vertices]
        result = extrema(admissible, coefficients)
        assert result["minimum"]["value"] == pytest.approx(min(values))
        assert result["maximum"]["value"] == pytest.approx(max(values))


@pytest.mark.parametrize("bounds", [[], [bound("a", 0, 1)],
                                    [bound("a", 0, None), bound("b", 0, 1)]])
def test_requesting_a_uniform_law_cannot_define_one_without_a_finite_explicit_box(bounds):
    admissible = po.MagnitudeDomain(
        order(("a", "b")), bounds, gap=0, gap_source=SOURCE,
    )
    batch = sample(admissible, n=5, budget=10)
    report = batch.report()
    assert report["status"] == "unidentified"
    assert "uniform" in report["requested_proposal"]
    assert report["proposal_status"] == "undefined_missing_bounds"
    assert report["proposal_distribution"] is report["distribution"] is None
    assert report["conditioning_status"] == "undefined_proposal"
    assert report["accepted"] == report["attempted"] == report["rejected"] == 0
    assert report["acceptance_rate"] is None
    assert report["proposal_source"]["reference"] == SOURCE.reference
    forward = po.evaluate_forward(batch, lambda point: 1, estimand="no proposal", unit="toy")
    assert forward["sampling"]["distribution"] is None
    assert forward["sign_stability"] == "inconclusive"
    json.dumps(report, allow_nan=False)


def test_order_propagated_finite_bounds_still_do_not_define_a_box_proposal():
    admissible = po.MagnitudeDomain(
        order(edges=(("a", "b"), ("b", "c"))),
        anchors=[po.Anchor("a", 0, SCOPE, SOURCE), po.Anchor("c", 1, SCOPE, SOURCE)],
        gap=0, gap_source=SOURCE,
    )
    assert extrema(admissible, {"b": 1})["status"] == "bounded"
    report = sample(admissible).report()
    assert report["proposal_status"] == "undefined_missing_bounds"
    assert report["proposal_distribution"] is report["distribution"] is None


def test_an_empty_explicit_box_defines_neither_base_nor_conditioned_law():
    admissible = po.MagnitudeDomain(
        order(("a",), ()), [bound("a", 0, 1)],
        anchors=[po.Anchor("a", 2, SCOPE, SOURCE)], gap=0, gap_source=SOURCE,
    )
    report = sample(admissible).report()
    assert report["status"] == "infeasible"
    assert report["proposal_status"] == "undefined_empty_box"
    assert report["proposal_distribution"] is report["distribution"] is None
    assert report["conditioning_status"] == "undefined_proposal"
    assert report["attempted"] == 0


@pytest.mark.parametrize("fixed", [False, True])
def test_a_valid_base_box_does_not_define_conditioning_on_an_empty_event(fixed):
    admissible = po.MagnitudeDomain(
        order(("a", "b")), [bound("a", 2, 2 if fixed else 3), bound("b", 0, 0 if fixed else 1)],
        gap=0, gap_source=SOURCE,
    )
    report = sample(admissible).report()
    assert report["status"] == "infeasible"
    assert report["proposal_status"] == "defined"
    assert "uniform" in report["proposal_distribution"]
    assert report["distribution"] is None
    assert report["conditioning_status"] == "empty"
    assert report["attempted"] == 0


@pytest.mark.parametrize("budget", [0, 29])
@pytest.mark.parametrize("case", ["anchored-upper", "anchored-lower", "gap", "path"])
def test_zero_probability_conditioning_is_not_defined_by_lp_feasibility_or_a_budget(case, budget):
    if case == "anchored-upper":
        admissible = po.MagnitudeDomain(
            order(("a", "b")), [bound("a", 0, 1)],
            anchors=[po.Anchor("b", 0, SCOPE, SOURCE)], gap=0, gap_source=SOURCE,
        )
    elif case == "anchored-lower":
        admissible = po.MagnitudeDomain(
            order(("a", "b")), [bound("b", 0, 1)],
            anchors=[po.Anchor("a", 1, SCOPE, SOURCE)], gap=0, gap_source=SOURCE,
        )
    elif case == "gap":
        admissible = domain(order(("a", "b")), gap=1)
    else:
        admissible = domain(order(edges=(("a", "b"), ("b", "c"))), gap=0.5)
    report = sample(admissible, n=3, seed=101, budget=budget).report()
    assert report["feasibility"]["status"] == "optimal"
    assert report["status"] == "budget_exhausted"
    assert report["attempted"] == report["rejected"] == budget
    assert report["accepted"] == 0
    assert report["proposal_status"] == "defined"
    assert report["proposal_distribution"] is not None
    assert report["conditioning_status"] == "zero_probability"
    assert report["distribution"] is None


def test_a_floating_point_boundary_hit_does_not_create_a_conditioned_law(monkeypatch):
    admissible = po.MagnitudeDomain(
        order(("a", "b")), [bound("a", 0, 1)],
        anchors=[po.Anchor("b", 0, SCOPE, SOURCE)], gap=0, gap_source=SOURCE,
    )
    assert admissible.contains({"a": 0, "b": 0})
    monkeypatch.setattr(po.np.random, "default_rng", lambda seed: SimpleNamespace(
        uniform=lambda lower, upper: lower.copy(),
    ))
    batch = sample(admissible, n=2, budget=4)
    assert batch.points == ()
    assert batch.attempted == 4
    assert batch.status == "budget_exhausted"
    assert batch.report()["distribution"] is None


@pytest.mark.parametrize("budget", [0, 10])
def test_no_acceptances_do_not_disprove_a_positive_probability_conditioning_event(budget):
    admissible = po.MagnitudeDomain(
        order(("a", "b")), [bound("a", 0, 1)],
        anchors=[po.Anchor("b", 2 ** -40, SCOPE, SOURCE)], gap=0, gap_source=SOURCE,
    )
    report = sample(admissible, n=20, seed=101, budget=budget).report()
    assert report["status"] == "budget_exhausted"
    assert report["accepted"] == 0
    assert report["attempted"] == budget
    assert report["proposal_status"] == "defined"
    assert report["conditioning_status"] == "positive_probability"
    assert report["distribution"] is not None


@pytest.mark.parametrize("gap", [0.0, 0.25])
def test_explicit_point_masses_are_not_mistaken_for_zero_probability_conditioning(gap):
    admissible = po.MagnitudeDomain(
        order(("a", "b")), anchors=[po.Anchor("a", 0, SCOPE, SOURCE),
                                    po.Anchor("b", gap, SCOPE, SOURCE)],
        gap=gap, gap_source=SOURCE,
    )
    report = sample(admissible, n=4, budget=4).report()
    assert report["status"] == "complete"
    assert report["accepted"] == 4
    assert report["proposal_free_dimensions"] == 0
    assert "point masses" in report["proposal_distribution"]
    assert report["conditioning_status"] == "positive_probability"
    assert report["distribution"] is not None


def test_declared_ties_define_the_base_measure_in_quotient_not_ambient_coordinates():
    admissible = domain(order(("a", "b"), (), ties=[po.DeclaredTie(("a", "b"), SCOPE, SOURCE)]))
    report = sample(admissible, n=4, budget=4).report()
    assert report["proposal_free_dimensions"] == 1
    assert report["conditioning_status"] == "positive_probability"
    assert report["distribution"] is not None
    assert all(point["a"] == point["b"] for point in report["points"])


@pytest.mark.parametrize("upper,expected", [(0.3, "empty"), (3 * 0.1, "positive_probability")])
def test_support_check_does_not_use_lp_tolerance_or_rounded_path_gaps(upper, expected):
    admissible = domain(
        order(("a", "b", "c", "d"), (("a", "b"), ("b", "c"), ("c", "d"))),
        upper=upper, gap=0.1,
    )
    report = sample(admissible, budget=0).report()
    assert report["feasibility"]["status"] == "optimal"
    assert report["conditioning_status"] == expected
    assert (report["distribution"] is not None) == (expected == "positive_probability")
    assert "binary-float" in report["measure_check"]


def test_conditioning_support_is_invariant_to_parameter_names_and_edge_order():
    for names in permutations(("a", "b", "c")):
        for left, right, expected in ((0.5, 0.5, "zero_probability"),
                                      (0.25, 0.75, "positive_probability")):
            poset = order(names, ((names[1], names[2]), (names[0], names[1])))
            admissible = po.MagnitudeDomain(
                poset, [bound(names[1], 0, 1)],
                anchors=[po.Anchor(names[0], left, SCOPE, SOURCE),
                         po.Anchor(names[2], right, SCOPE, SOURCE)],
                gap=0, gap_source=SOURCE,
            )
            report = sample(admissible, budget=0).report()
            assert report["conditioning_status"] == expected


def test_proposal_existence_is_independent_of_a_numerical_sampling_failure(monkeypatch):
    def broken_uniform(lower, upper):
        raise FloatingPointError("deliberate proposal execution failure")

    monkeypatch.setattr(po.np.random, "default_rng", lambda seed: SimpleNamespace(uniform=broken_uniform))
    report = sample(domain()).report()
    assert report["status"] == "numerical_failure"
    assert report["attempted"] == report["accepted"] == 0
    assert report["proposal_distribution"] is not None
    assert report["conditioning_status"] == "positive_probability"
    assert report["distribution"] is not None


def test_seeded_sampling_is_explicit_uniform_box_conditioning_and_reproducible():
    admissible = domain(order(("a", "b"), (("a", "b"),)))
    batch = sample(admissible, n=5000, budget=20000)
    assert batch.status == "complete"
    assert batch.points == sample(admissible, n=5000, budget=20000).points
    assert batch.points != sample(admissible, n=5000, seed=4, budget=20000).points
    report = batch.report()
    assert report["accepted"] == 5000
    assert report["attempted"] <= report["budget"]
    assert 0.46 < report["acceptance_rate"] < 0.54
    assert "uniform" in report["distribution"]
    assert "quotient" in report["distribution"]
    assert not report["biological_posterior"]
    assert report["proposal_source"]["reference"] == SOURCE.reference
    points = [dict(zip(admissible.order.items, row)) for row in batch.points]
    assert all(admissible.contains(point) for point in points)
    values = np.asarray(batch.points)
    assert values[:, 0].mean() == pytest.approx(1 / 3, abs=0.025)
    assert values[:, 1].mean() == pytest.approx(2 / 3, abs=0.025)


def test_tie_quotient_and_explicit_fixed_coordinates_are_not_zero_measure_rejection():
    poset = order(ties=[po.DeclaredTie(("a", "c"), SCOPE, SOURCE)])
    admissible = po.MagnitudeDomain(
        poset, [bound("a", 0, 1), bound("c", 0.2, 0.8)],
        anchors=[po.Anchor("b", 1, SCOPE, SOURCE)], gap=0.1, gap_source=SOURCE,
    )
    batch = sample(admissible)
    assert batch.status == "complete"
    for a, b, c in batch.points:
        assert a == c
        assert b == 1
        assert 0.2 <= a <= 0.8


def test_sampling_never_uses_tolerance_to_enlarge_the_conditioned_domain():
    poset = order(("a", "b"), (("a", "b"),))
    admissible = domain(poset, upper=1e-9, gap=1e-9, tolerance=1e-8)
    batch = sample(admissible, n=5, budget=100)
    assert batch.status == "budget_exhausted"
    assert batch.points == ()
    assert batch.attempted == 100


def test_budget_exhaustion_has_no_biased_fallback_or_discarded_denominator():
    batch = sample(domain(), n=10, budget=1)
    assert batch.status == "budget_exhausted"
    assert batch.attempted == 1
    assert len(batch.points) <= 1
    assert batch.report()["requested"] == 10
    empty = sample(domain(), n=10, budget=0)
    assert empty.status == "budget_exhausted" and not empty.points
    assert empty.report()["acceptance_rate"] is None


@pytest.mark.parametrize("kwargs", [{"n": 0}, {"n": True}, {"seed": -1},
                                    {"budget": -1}, {"budget": 1.5}])
def test_invalid_sampling_controls_fail_early(kwargs):
    with pytest.raises(ValueError):
        sample(domain(), **kwargs)


def test_forward_callback_reports_only_sampled_extrema_and_preserves_failures():
    batch = sample(domain(), n=100)

    def forward(values):
        if values["c"] > 0.8:
            raise RuntimeError("deliberate synthetic integration failure")
        return (values["b"] - values["a"]) ** 2 + 0.01

    report = po.evaluate_forward(batch, forward, estimand="neutral nonlinear response", unit="toy")
    assert not report["globally_certified"]
    assert report["sampled_minimum"] >= 0.01
    assert report["sampled_maximum"] <= 1.01
    assert report["failed"] > 0
    assert report["successful"] + report["failed"] == 100
    assert len(report["records"]) == report["accepted"] == 100
    assert report["sign_stability"] == "inconclusive"
    assert report["successful_sample_sign"] == "positive"
    failed = [row for row in report["records"] if row["status"] == "failed"]
    assert all(row["point"]["c"] > 0.8 for row in failed)
    assert all("deliberate" in row["error"] for row in failed)
    assert report["sampling"]["proposal_source"]["kind"] == "synthetic"
    json.dumps(report, allow_nan=False)


def test_forward_nonfinite_outputs_are_failures_and_all_failed_extrema_are_missing():
    batch = sample(domain(), n=10)
    report = po.evaluate_forward(batch, lambda p: math.nan, estimand="invalid output", unit="toy")
    assert report["successful"] == 0 and report["failed"] == 10
    assert report["sampled_minimum"] is report["sampled_maximum"] is None
    assert report["successful_sample_sign"] == "unidentified"
    assert report["sign_stability"] == "inconclusive"


def test_forward_callback_cannot_mutate_the_samples():
    batch = sample(domain(), n=2)
    before = batch.points

    def mutating(point):
        point["a"] = 100
        return 1

    report = po.evaluate_forward(batch, mutating, estimand="mutation attempt", unit="toy")
    assert batch.points == before
    assert report["failed"] == 2


def test_forward_sign_stability_is_sample_only_even_with_no_failures():
    batch = sample(domain(), n=100)
    positive = po.evaluate_forward(batch, lambda p: 1 + p["a"] ** 2,
                                   estimand="positive toy response", unit="toy")
    mixed = po.evaluate_forward(batch, lambda p: math.sin(2 * math.pi * p["c"]),
                                estimand="oscillating toy response", unit="toy")
    assert positive["sign_stability"] == "positive_on_samples"
    assert mixed["sign_stability"] == "crosses_zero_on_samples"
    assert not positive["globally_certified"] and not mixed["globally_certified"]
    incomplete = po.evaluate_forward(sample(domain(), n=10, budget=1), lambda p: 1,
                                     estimand="incomplete experiment", unit="toy")
    assert incomplete["sign_stability"] == "inconclusive"


def test_synthetic_cli_retains_budget_failure_in_its_results(tmp_path):
    root = Path(__file__).resolve().parents[1]
    script = root / "scripts" / "run_order_robustness.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output-dir", str(tmp_path), "--budget", "0"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 1
    report = json.loads((tmp_path / "partial_order_synthetic_demo.json").read_text())
    assert report["sampling"]["status"] == "budget_exhausted"
    assert report["sampling"]["accepted"] == report["sampling"]["attempted"] == 0
    assert report["sampling"]["conditioning_status"] == "positive_probability"
    assert report["sampling"]["distribution"] is not None
    assert report["order_only_sampling"]["distribution"] is None
    for row in report["nonlinear"].values():
        assert row["sampled_minimum"] is None
        assert row["sign_stability"] == "inconclusive"
        assert not row["globally_certified"]


def test_synthetic_cli_writes_strict_json_and_refuses_to_overwrite(tmp_path):
    root = Path(__file__).resolve().parents[1]
    script = root / "scripts" / "run_order_robustness.py"
    command = [sys.executable, str(script), "--output-dir", str(tmp_path),
               "--samples", "30", "--seed", "5", "--budget", "1000"]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    path = tmp_path / "partial_order_synthetic_demo.json"
    report = json.loads(path.read_text())
    assert report["schema_version"] == 2
    assert report["order_only_sampling"]["distribution"] is None
    assert report["order_only_sampling"]["proposal_distribution"] is None
    assert report["sampling"]["conditioning_status"] == "positive_probability"
    assert report["synthetic_only"] is True
    assert report["biology_data_used"] is False
    assert report["law_ready"] is False
    assert report["affine"]["robust"]["sign"] == "positive"
    assert report["affine"]["nonrobust"]["sign"] == "crosses_zero"
    assert report["order_only"]["status"] == "unbounded"
    assert all(not row["globally_certified"] for row in report["nonlinear"].values())
    with (tmp_path / "partial_order_synthetic_summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if row["method"] == "sampled nonlinear callback":
            assert row["requested_proposal"] == report["sampling"]["requested_proposal"]
            assert row["proposal_distribution"] == report["sampling"]["proposal_distribution"]
            assert row["proposal"] == report["sampling"]["distribution"]
            assert row["conditioning_status"] == "positive_probability"
            assert row["proposal_source"] == report["sampling"]["proposal_source"]["reference"]
    before = path.read_bytes()
    rerun = subprocess.run(command, capture_output=True, text=True, check=False)
    assert rerun.returncode != 0
    assert path.read_bytes() == before
