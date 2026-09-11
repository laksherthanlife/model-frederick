from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from itertools import combinations
from math import fsum, inf, isfinite
from numbers import Integral, Real
from types import MappingProxyType

import numpy as np
from scipy.optimize import linprog


__all__ = [
    "Source", "Calibration", "Scope", "OrderClaim", "DeclaredTie", "Conversion",
    "CycleError", "PartialOrder", "compile_order", "Bound", "Anchor", "MagnitudeDomain",
    "affine_extrema", "SampleBatch", "sample_uniform_box", "ForwardEvaluation", "evaluate_forward",
]

_MAGNITUDE_KINDS = {"magnitude", "effect", "importance"}
_RELATION_KINDS = _MAGNITUDE_KINDS | {"temporal_precedence", "causal_dependency"}
_CERTIFICATE = "numerical SciPy/HiGHS LP certificate, not a rational exact proof"


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _number(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return value


def _instance(value, cls, name):
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be {cls.__name__}")


def _scoped(value):
    _instance(value.scope, Scope, "scope")
    _instance(value.source, Source, "source provenance")


@dataclass(frozen=True)
class Source:
    kind: str
    reference: str
    method: str

    def __post_init__(self):
        if self.kind not in {"measured", "assumed", "simulated", "synthetic"}:
            raise ValueError("source kind must be measured, assumed, simulated, or synthetic")
        _text(self.reference, "source reference")
        _text(self.method, "assay or method")


@dataclass(frozen=True)
class Calibration:
    family: str
    gain: str

    def __post_init__(self):
        _text(self.family, "calibration family")
        if self.gain not in {"known", "shared_unknown_positive"}:
            raise ValueError("gain must be known or shared_unknown_positive within this family")


@dataclass(frozen=True)
class Scope:
    relation_kind: str
    context: str
    quantity: str
    unit: str
    calibration: Calibration
    estimand: str | None = None

    def __post_init__(self):
        if self.relation_kind not in _RELATION_KINDS:
            raise ValueError(f"unsupported relation kind: {self.relation_kind!r}")
        for name in ("context", "quantity", "unit"):
            _text(getattr(self, name), name)
        _instance(self.calibration, Calibration, "calibration")
        if self.relation_kind in {"effect", "importance"} or self.estimand is not None:
            _text(self.estimand, "named estimand")


@dataclass(frozen=True)
class OrderClaim:
    lower: str
    upper: str
    scope: Scope
    source: Source

    def __post_init__(self):
        _text(self.lower, "lower item")
        _text(self.upper, "upper item")
        _scoped(self)


@dataclass(frozen=True)
class DeclaredTie:
    members: tuple[str, ...]
    scope: Scope
    source: Source

    def __post_init__(self):
        if isinstance(self.members, str):
            raise ValueError("tie members must be a sequence of distinct item names")
        members = tuple(self.members)
        for member in members:
            _text(member, "tie member")
        if len(members) < 2 or len(set(members)) != len(members):
            raise ValueError("a declared tie needs at least two distinct members")
        object.__setattr__(self, "members", members)
        _scoped(self)


@dataclass(frozen=True)
class Conversion:
    origin: Scope
    target: Scope
    factor: float
    offset: float
    source: Source

    def __post_init__(self):
        _instance(self.origin, Scope, "conversion origin")
        _instance(self.target, Scope, "conversion target")
        _instance(self.source, Source, "conversion source")
        if self.origin.relation_kind != self.target.relation_kind:
            raise ValueError("conversion cannot change relation kind, including causal-as-magnitude")
        if self.origin.quantity != self.target.quantity:
            raise ValueError("conversion cannot change quantity")
        if self.origin.estimand != self.target.estimand:
            raise ValueError("conversion cannot change estimand")
        if self.origin.relation_kind == "causal_dependency":
            raise ValueError("causal dependencies do not have affine magnitude conversions")
        factor = _number(self.factor, "conversion factor")
        if factor <= 0:
            raise ValueError("conversion factor must be positive to preserve order")
        offset = _number(self.offset, "conversion offset")
        if self.origin == self.target and (factor != 1.0 or offset != 0.0):
            raise ValueError("a nonidentity conversion needs distinct origin and target scopes")
        object.__setattr__(self, "factor", factor)
        object.__setattr__(self, "offset", offset)


def _conversions(rows, target):
    rows = tuple(rows)
    origins = set()
    for row in rows:
        _instance(row, Conversion, "conversion")
        if row.target != target:
            raise ValueError("each conversion must explicitly target this poset's scope")
        if row.origin in origins:
            raise ValueError("ambiguous duplicate conversion origin")
        origins.add(row.origin)
    return rows


def _transform(origin, target, conversions):
    if origin == target:
        return 1.0, 0.0
    for conversion in conversions:
        if conversion.origin == origin and conversion.target == target:
            return conversion.factor, conversion.offset
    raise ValueError("incompatible scope: explicit direct conversion evidence is required")


class CycleError(ValueError):
    def __init__(self, cycle, claims, relation_kind):
        self.cycle = tuple(cycle)
        self.claims = tuple(claims)
        path = " -> ".join("=".join(group) for group in self.cycle)
        sources = ", ".join(sorted({claim.source.reference for claim in claims}))
        suffix = "; causal feedback is not a poset" if relation_kind == "causal_dependency" else ""
        super().__init__(f"strict cycle {path}; sources: {sources}{suffix}")


@dataclass(frozen=True)
class PartialOrder:
    items: tuple[str, ...]
    scope: Scope
    claims: tuple[OrderClaim, ...]
    ties: tuple[DeclaredTie, ...]
    conversions: tuple[Conversion, ...]
    classes: tuple[tuple[str, ...], ...]
    closure: tuple[tuple[str, str], ...]
    incomparable_pairs: tuple[tuple[str, str], ...]
    _edges: tuple[tuple[int, int], ...] = field(repr=False)

    def compare(self, left: str, right: str) -> str:
        if left not in self.items or right not in self.items:
            raise ValueError("unknown item in comparison")
        if left == right:
            return "same"
        if any(left in group and right in group for group in self.classes):
            return "tied"
        if (left, right) in self.closure:
            return "less"
        if (right, left) in self.closure:
            return "greater"
        return "incomparable"

    def report(self) -> dict:
        return {
            "scope": asdict(self.scope), "items": list(self.items),
            "claims": [asdict(row) for row in self.claims],
            "declared_ties": [asdict(row) for row in self.ties],
            "conversions": [asdict(row) for row in self.conversions],
            "classes": [list(group) for group in self.classes],
            "strict_closure": [list(pair) for pair in self.closure],
            "incomparable_pairs": [list(pair) for pair in self.incomparable_pairs],
            "source_kinds": sorted({row.source.kind for row in (*self.claims, *self.ties)}),
            "incomparable_meaning": "unknown, not a tie or an arbitrarily completed ranking",
            "claim_direction": "lower strictly precedes upper in the declared relation kind",
            "scale_note": "a shared unknown positive gain preserves order but supplies no range",
        }


def compile_order(
    items: Iterable[str], scope: Scope, claims: Iterable[OrderClaim] = (), *,
    ties: Iterable[DeclaredTie] = (), conversions: Iterable[Conversion] = (),
) -> PartialOrder:
    _instance(scope, Scope, "scope")
    if isinstance(items, str):
        raise ValueError("items must be a sequence of distinct names")
    items = tuple(items)
    for item in items:
        _text(item, "item")
    if not items or len(set(items)) != len(items):
        raise ValueError("items must be nonempty and distinct")
    items = tuple(sorted(items))
    claims, ties = tuple(claims), tuple(ties)
    conversions = _conversions(conversions, scope)
    parent = {item: item for item in items}

    def root(item):
        while item != parent[item]:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for rows, cls in ((claims, OrderClaim), (ties, DeclaredTie)):
        for row in rows:
            _instance(row, cls, "order input")
            _transform(row.scope, scope, conversions)
            names = (row.lower, row.upper) if isinstance(row, OrderClaim) else row.members
            if any(name not in parent for name in names):
                raise ValueError("unknown item in order input")
    for tie in ties:
        for member in tie.members[1:]:
            parent[root(member)] = root(tie.members[0])
    grouped = {}
    for item in items:
        grouped.setdefault(root(item), []).append(item)
    classes = tuple(sorted(tuple(group) for group in grouped.values()))
    index = {name: i for i, group in enumerate(classes) for name in group}
    edges = tuple(sorted({(index[row.lower], index[row.upper]) for row in claims}))
    adjacency = {i: [] for i in range(len(classes))}
    for left, right in edges:
        adjacency[left].append(right)
    done, postorder = set(), []
    for start in adjacency:
        if start in done:
            continue
        path, active = [start], {start: 0}
        stack = [(start, iter(adjacency[start]))]
        while stack:
            node, children = stack[-1]
            child = next(children, None)
            if child is None:
                done.add(node)
                postorder.append(node)
                active.pop(node)
                path.pop()
                stack.pop()
            elif child in active:
                cycle = path[active[child]:] + [child]
                cycle_edges = set(zip(cycle, cycle[1:]))
                evidence = [row for row in claims
                            if (index[row.lower], index[row.upper]) in cycle_edges]
                raise CycleError([classes[i] for i in cycle], evidence, scope.relation_kind)
            elif child not in done:
                active[child] = len(path)
                path.append(child)
                stack.append((child, iter(adjacency[child])))
    reachable = {i: set() for i in adjacency}
    for node in postorder:
        for child in adjacency[node]:
            reachable[node].add(child)
            reachable[node].update(reachable[child])
    closure = tuple((a, b) for a in items for b in items if index[b] in reachable[index[a]])
    incomparable = tuple((a, b) for a, b in combinations(items, 2)
                         if index[a] != index[b] and index[b] not in reachable[index[a]]
                         and index[a] not in reachable[index[b]])
    return PartialOrder(items, scope, claims, ties, conversions, classes, closure, incomparable, edges)


@dataclass(frozen=True)
class Bound:
    item: str
    lower: float | None
    upper: float | None
    scope: Scope
    source: Source

    def __post_init__(self):
        _text(self.item, "bounded item")
        _scoped(self)
        if self.lower is None and self.upper is None:
            raise ValueError("a bound must supply at least one endpoint; omit absent bounds")
        for name in ("lower", "upper"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _number(getattr(self, name), name))
        if self.lower is not None and self.upper is not None and self.lower > self.upper:
            raise ValueError("lower bound exceeds upper bound")


@dataclass(frozen=True)
class Anchor:
    item: str
    value: float
    scope: Scope
    source: Source

    def __post_init__(self):
        _text(self.item, "anchored item")
        _scoped(self)
        object.__setattr__(self, "value", _number(self.value, "anchor"))


@dataclass(frozen=True)
class MagnitudeDomain:
    order: PartialOrder
    bounds: tuple[Bound, ...] = ()
    anchors: tuple[Anchor, ...] = ()
    conversions: tuple[Conversion, ...] = ()
    gap: float = field(kw_only=True)
    gap_source: Source = field(kw_only=True)
    tolerance: float = field(default=1e-8, kw_only=True)
    _limits: tuple[tuple[float | None, float | None], ...] = field(init=False, repr=False)

    def __post_init__(self):
        _instance(self.order, PartialOrder, "compiled order")
        if self.order.scope.relation_kind not in _MAGNITUDE_KINDS:
            raise ValueError("only magnitude, effect, or named importance orders define magnitude domains")
        _instance(self.gap_source, Source, "gap source")
        gap = _number(self.gap, "nonnegative gap")
        if gap < 0:
            raise ValueError("embedding gap must be nonnegative")
        tolerance = _number(self.tolerance, "LP tolerance")
        if not 1e-10 <= tolerance <= 1e-4:
            raise ValueError("LP tolerance must be between 1e-10 and 1e-4 for HiGHS validation")
        object.__setattr__(self, "gap", gap)
        object.__setattr__(self, "tolerance", tolerance)
        object.__setattr__(self, "bounds", tuple(self.bounds))
        object.__setattr__(self, "anchors", tuple(self.anchors))
        object.__setattr__(self, "conversions", tuple(self.conversions))
        conversions = _conversions((*self.order.conversions, *self.conversions), self.order.scope)
        index = self._index()
        limits = [[None, None] for _ in self.order.classes]
        for rows, cls in ((self.bounds, Bound), (self.anchors, Anchor)):
            for row in rows:
                _instance(row, cls, "range input")
                if row.item not in index:
                    raise ValueError(f"unknown bounded or anchored item: {row.item!r}")
                factor, offset = _transform(row.scope, self.order.scope, conversions)
                endpoints = (row.lower, row.upper) if isinstance(row, Bound) else (row.value, row.value)
                for side, value in enumerate(endpoints):
                    if value is not None:
                        converted = _number(factor * value + offset, "converted bound or anchor")
                        current = limits[index[row.item]][side]
                        combine = max if side == 0 else min
                        limits[index[row.item]][side] = converted if current is None else combine(
                            current, converted,
                        )
        object.__setattr__(self, "_limits", tuple(tuple(pair) for pair in limits))

    def _index(self):
        return {name: i for i, group in enumerate(self.order.classes) for name in group}

    def _point(self, coordinates):
        return {name: float(coordinates[i]) for name, i in self._index().items()}

    def _matrix(self):
        matrix = np.zeros((len(self.order._edges), len(self.order.classes)))
        for row, (left, right) in enumerate(self.order._edges):
            matrix[row, left], matrix[row, right] = 1.0, -1.0
        return matrix, np.full(len(matrix), -self.gap)

    def violation(self, point: Mapping[str, float]) -> float:
        if set(point) != set(self.order.items):
            raise ValueError("point must contain every item exactly once and no unknown items")
        point = {name: _number(value, "point coordinate") for name, value in point.items()}
        coordinates = [point[group[0]] for group in self.order.classes]
        errors = [0.0]
        for group, value, (lower, upper) in zip(self.order.classes, coordinates, self._limits):
            errors.extend(abs(point[name] - value) for name in group)
            if lower is not None:
                errors.append(lower - value)
            if upper is not None:
                errors.append(value - upper)
        errors.extend(self.gap - (coordinates[right] - coordinates[left])
                      for left, right in self.order._edges)
        return max(errors)

    def contains(self, point: Mapping[str, float], *, tolerance: float = 0.0) -> bool:
        tolerance = _number(tolerance, "membership tolerance")
        if tolerance < 0:
            raise ValueError("membership tolerance must be nonnegative")
        return self.violation(point) <= tolerance

    def feasibility(self) -> dict:
        result = _solve(self, np.zeros(len(self.order.classes)))
        if result["status"] == "unbounded":
            result = _empty_endpoint("numerical_failure", "zero-objective feasibility cannot be unbounded")
        result.update(tolerance=self.tolerance, certificate=_CERTIFICATE)
        return result

    def report(self) -> dict:
        return {
            "order": self.order.report(),
            "bounds": [asdict(row) for row in self.bounds],
            "anchors": [asdict(row) for row in self.anchors],
            "additional_conversions": [asdict(row) for row in self.conversions],
            "effective_box": [
                {"members": list(group), "lower": lower, "upper": upper}
                for group, (lower, upper) in zip(self.order.classes, self._limits)
            ],
            "gap": {"value": self.gap, "unit": self.order.scope.unit,
                    "source": asdict(self.gap_source),
                    "meaning": "minimum separation on each declared edge; path gaps accumulate"},
            "embedding": "weak_order_preserving" if self.gap == 0 else "positive_gap_order_preserving",
            "positive_gap_exceeds_tolerance": self.gap > self.tolerance,
            "order_reflecting": False,
            "range_information": "explicit bounds/anchors; query-dependent identification"
            if self.bounds or self.anchors else "none; no absolute scale or finite range supplied",
            "lp_tolerance": self.tolerance,
            "probability_distribution": "none unless an explicit sampling proposal is requested",
        }


def _maxabs(values):
    return float(np.max(np.abs(values), initial=0.0))


def _empty_endpoint(status, message):
    return {"status": status, "value": None, "point": None, "validated": False,
            "residuals": None, "message": message}


@np.errstate(over="raise", invalid="raise", divide="raise")
def _solve(domain, coefficients):
    for lower, upper in domain._limits:
        if lower is not None and upper is not None and lower > upper:
            return _empty_endpoint("infeasible", "conflicting explicit bounds/anchors in a quotient class")
    matrix, rhs = domain._matrix()
    try:
        result = linprog(
            coefficients, A_ub=matrix if len(matrix) else None,
            b_ub=rhs if len(matrix) else None, bounds=domain._limits, method="highs",
            options={"primal_feasibility_tolerance": domain.tolerance,
                     "dual_feasibility_tolerance": domain.tolerance},
        )
    except (ValueError, RuntimeError, OverflowError, FloatingPointError) as exc:
        return _empty_endpoint("numerical_failure", f"LP solver error: {exc}")
    if result.status != 0:
        status = {2: "infeasible", 3: "unbounded"}.get(result.status, "numerical_failure")
        return _empty_endpoint(status, str(result.message))
    try:
        coordinates = np.asarray(result.x, dtype=float)
        dual = np.asarray(result.ineqlin.marginals, dtype=float)
        lower_dual = np.asarray(result.lower.marginals, dtype=float)
        upper_dual = np.asarray(result.upper.marginals, dtype=float)
        vectors = (coordinates, dual, lower_dual, upper_dual)
        expected = (len(coefficients), len(matrix), len(coefficients), len(coefficients))
        if any(vector.shape != (size,) or not np.isfinite(vector).all()
               for vector, size in zip(vectors, expected)):
            raise ValueError("invalid solver primal/dual vectors")
        point = domain._point(coordinates)
        primal = _number(float(coefficients @ coordinates), "LP objective")
        reported = _number(result.fun, "solver objective")
        lower = np.array([0.0 if pair[0] is None else pair[0] for pair in domain._limits])
        upper = np.array([0.0 if pair[1] is None else pair[1] for pair in domain._limits])
        has_lower = np.array([pair[0] is not None for pair in domain._limits])
        has_upper = np.array([pair[1] is not None for pair in domain._limits])
        dual_objective = _number(float(rhs @ dual + lower @ lower_dual + upper @ upper_dual),
                                 "dual objective")
        scale = max(1.0, abs(primal), abs(dual_objective))
        stationarity = matrix.T @ dual + lower_dual + upper_dual
        dual_scale = max(1.0, _maxabs(coefficients), _maxabs(stationarity))
        residuals = {
            "primal_absolute": domain.violation(point),
            "dual_sign_absolute": max(0.0, float(np.max(dual, initial=0.0)),
                                      float(np.max(-lower_dual, initial=0.0)),
                                      float(np.max(upper_dual, initial=0.0))),
            "absent_bound_dual_absolute": max(_maxabs(lower_dual[~has_lower]),
                                              _maxabs(upper_dual[~has_upper])),
            "stationarity_relative": _maxabs(coefficients - stationarity) / dual_scale,
            "duality_gap_relative": abs(primal - dual_objective) / scale,
            "objective_relative": abs(primal - reported) / max(1.0, abs(primal), abs(reported)),
            "complementarity_relative": max(
                _maxabs(dual * (rhs - matrix @ coordinates)),
                _maxabs(lower_dual[has_lower] * (coordinates - lower)[has_lower]),
                _maxabs(upper_dual[has_upper] * (upper - coordinates)[has_upper]),
            ) / scale,
        }
        if not all(isfinite(value) for value in residuals.values()):
            raise ValueError("nonfinite certificate residual")
        valid = all(value <= domain.tolerance for value in residuals.values())
        endpoint = _empty_endpoint("numerical_failure", "solver optimum failed residual validation")
        endpoint["residuals"] = residuals
        if valid:
            endpoint.update(status="optimal", value=primal, point=point,
                            validated=True, message=str(result.message))
        return endpoint
    except (AttributeError, TypeError, ValueError, OverflowError, FloatingPointError) as exc:
        return _empty_endpoint("numerical_failure", f"invalid solver certificate: {exc}")


def _sign(lower, upper, tolerance):
    if lower > tolerance:
        return "positive"
    if upper < -tolerance:
        return "negative"
    if lower >= -tolerance and upper <= tolerance:
        return "within_zero_tolerance"
    if lower >= -tolerance:
        return "nonnegative_within_tolerance"
    if upper <= tolerance:
        return "nonpositive_within_tolerance"
    return "crosses_zero"


def affine_extrema(
    domain: MagnitudeDomain, coefficients: Mapping[str, float], *, estimand: str,
    unit: str, offset: float = 0.0,
) -> dict:
    _instance(domain, MagnitudeDomain, "domain")
    _text(estimand, "query estimand")
    _text(unit, "query unit")
    offset = _number(offset, "query offset")
    if any(name not in domain.order.items for name in coefficients):
        raise ValueError("unknown item in affine coefficients")
    coefficients = {name: _number(value, "affine coefficient") for name, value in coefficients.items()}
    try:
        vector = np.array([_number(fsum(coefficients.get(name, 0.0) for name in group),
                                   "quotient coefficient") for group in domain.order.classes])
    except OverflowError as exc:
        raise ValueError("quotient coefficient must be finite") from exc
    feasibility = domain.feasibility()
    if feasibility["status"] != "optimal":
        minimum = _empty_endpoint(feasibility["status"], "no validated feasible domain")
        maximum = dict(minimum)
        status = feasibility["status"]
    else:
        minimum, maximum = _solve(domain, vector), _solve(domain, -vector)
        for endpoint, direction in ((minimum, 1.0), (maximum, -1.0)):
            if endpoint["status"] == "optimal":
                value = direction * endpoint["value"] + offset
                if isfinite(value):
                    endpoint["value"] = value
                else:
                    endpoint.update(_empty_endpoint("numerical_failure", "nonfinite query result"))
        states = {minimum["status"], maximum["status"]}
        if "numerical_failure" in states or "infeasible" in states:
            status = "numerical_failure"
        else:
            status = "unbounded" if "unbounded" in states else "bounded"
    sign = "unidentified"
    if status in {"bounded", "unbounded"}:
        lower = -inf if minimum["status"] == "unbounded" else minimum["value"]
        upper = inf if maximum["status"] == "unbounded" else maximum["value"]
        if lower <= upper + domain.tolerance:
            sign = _sign(lower, upper, domain.tolerance)
        else:
            status = "numerical_failure"
    return {
        "status": status, "estimand": estimand, "unit": unit,
        "coefficients": coefficients, "offset": offset,
        "minimum": minimum, "maximum": maximum, "sign": sign,
        "feasibility": feasibility, "tolerance": domain.tolerance,
        "certificate": _CERTIFICATE,
        "certificate_scope": "affine objective on the stated closed numeric domain only",
        "residual_convention": "primal/dual-sign residuals absolute; objective/KKT residuals relative",
        "sign_convention": "computed LP endpoints at the stated tolerance, not an exact strict-sign proof",
        "domain": domain.report(),
    }


def _uniform_box_measure(domain):
    requested = (
        "explicit uniform quotient-coordinate box conditioned on numeric orders; "
        "explicit fixed coordinates are point masses, not ambient tie conditioning"
    )
    result = {
        "requested_proposal": requested,
        "proposal_distribution": None, "proposal_status": "undefined_missing_bounds",
        "proposal_free_dimensions": None, "conditioning_status": "undefined_proposal",
        "distribution": None,
        "measure_check": "exact strict/non-strict bound propagation on the declared binary-float "
                         "bounds and gap, relative to free quotient coordinates and explicit point masses; "
                         "not an LP certificate or validation of biological assumptions",
    }
    limits = domain._limits
    if any(lower is not None and upper is not None and lower > upper for lower, upper in limits):
        result["proposal_status"] = "undefined_empty_box"
        return result
    if any(lower is None or upper is None for lower, upper in limits):
        return result
    free = [lower < upper for lower, upper in limits]
    result.update(
        proposal_status="defined", proposal_free_dimensions=sum(free),
        proposal_distribution="independent uniform distributions on finite explicit quotient-coordinate "
                              "intervals; explicit fixed coordinates are point masses",
    )
    lower = [(Fraction(pair[0]), is_free) for pair, is_free in zip(limits, free)]
    upper = [Fraction(pair[1]) for pair in limits]
    gap = Fraction(domain.gap)
    children, indegree = [[] for _ in limits], [0 for _ in limits]
    for left, right in domain.order._edges:
        children[left].append(right)
        indegree[right] += 1
    ready = [i for i, degree in enumerate(indegree) if degree == 0]
    conditioning = "positive_probability"
    for left in ready:
        value, strict = lower[left]
        if value > upper[left]:
            conditioning = "empty"
        elif conditioning != "empty" and value == upper[left] and (strict or free[left]):
            conditioning = "zero_probability"
        for right in children[left]:
            lower[right] = max(lower[right], (value + gap, strict or free[left] or free[right]))
            indegree[right] -= 1
            if indegree[right] == 0:
                ready.append(right)
    if len(ready) != len(limits):
        conditioning = "unidentified"
    result["conditioning_status"] = conditioning
    if conditioning == "positive_probability":
        result["distribution"] = requested
    return result


@dataclass(frozen=True)
class SampleBatch:
    domain: MagnitudeDomain
    points: tuple[tuple[float, ...], ...]
    requested: int
    attempted: int
    seed: int
    budget: int
    source: Source
    status: str
    message: str
    feasibility: dict

    def report(self, *, include_points: bool = True) -> dict:
        report = {
            "status": self.status, "message": self.message, "requested": self.requested,
            "accepted": len(self.points), "attempted": self.attempted, "budget": self.budget,
            "rejected": self.attempted - len(self.points), "seed": self.seed,
            "acceptance_rate": len(self.points) / self.attempted if self.attempted else None,
            **_uniform_box_measure(self.domain),
            "proposal_source": asdict(self.source), "biological_posterior": False,
            "membership_tolerance": 0.0,
            "feasibility": self.feasibility, "domain": self.domain.report(),
        }
        if include_points:
            report["points"] = [dict(zip(self.domain.order.items, point)) for point in self.points]
        return report


def sample_uniform_box(
    domain: MagnitudeDomain, *, n: int, seed: int, budget: int, source: Source,
) -> SampleBatch:
    _instance(domain, MagnitudeDomain, "domain")
    _instance(source, Source, "uniform proposal source")
    for name, value, minimum in (("n", n, 1), ("seed", seed, 0), ("budget", budget, 0)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    n, seed, budget = int(n), int(seed), int(budget)
    measure = _uniform_box_measure(domain)
    feasibility = domain.feasibility()
    points, attempted = [], 0
    status, message = feasibility["status"], "no validated feasible domain"
    if status == "optimal":
        if measure["proposal_status"] == "undefined_missing_bounds":
            status, message = "unidentified", "no finite explicit box for every quotient coordinate; " \
                "order-propagated bounds are not silently substituted for a proposal"
        elif measure["proposal_status"] == "undefined_empty_box":
            status, message = "infeasible", "an empty explicit box does not define a uniform proposal"
        else:
            rng = np.random.default_rng(seed)
            lower, upper = np.array(domain._limits).T
            try:
                while attempted < budget and len(points) < n:
                    coordinates = rng.uniform(lower, upper)
                    attempted += 1
                    point = domain._point(coordinates)
                    if measure["conditioning_status"] == "positive_probability" and domain.contains(point):
                        points.append(tuple(point[item] for item in domain.order.items))
                status = "complete" if len(points) == n else "budget_exhausted"
                message = "uniform rejection only; no biased fallback"
                if measure["conditioning_status"] != "positive_probability":
                    message += f"; conditioning is {measure['conditioning_status']}; " \
                        "no conditional law or accepted samples, even for floating-point boundary hits"
            except (ValueError, OverflowError, FloatingPointError) as exc:
                status, message = "numerical_failure", f"uniform proposal failed: {exc}"
    return SampleBatch(domain, tuple(points), n, attempted, seed, budget, source,
                       status, message, feasibility)


@dataclass(frozen=True)
class ForwardEvaluation:
    status: str
    value: float | None
    details: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if self.status not in {"optimal", "infeasible", "unbounded", "numerical_failure", "unsupported"}:
            raise ValueError("unknown forward solver status")
        if self.status == "optimal":
            object.__setattr__(self, "value", _number(self.value, "forward scalar output"))
        elif self.value is not None:
            raise ValueError("an unsolved forward evaluation cannot carry a numeric value")
        if not isinstance(self.details, Mapping):
            raise TypeError("forward details must be a mapping")
        object.__setattr__(self, "details", MappingProxyType(dict(self.details)))


def evaluate_forward(
    batch: SampleBatch, callback: Callable[[Mapping[str, float]], float | ForwardEvaluation], *,
    estimand: str, unit: str, sign_tolerance: float = 1e-8,
) -> dict:
    _instance(batch, SampleBatch, "sample batch")
    _text(estimand, "forward estimand")
    _text(unit, "forward unit")
    tolerance = _number(sign_tolerance, "sign tolerance")
    if tolerance < 0:
        raise ValueError("sign tolerance must be nonnegative")
    if not callable(callback):
        raise TypeError("forward callback must be callable")
    records, values = [], []
    for i, coordinates in enumerate(batch.points):
        point = dict(zip(batch.domain.order.items, coordinates))
        record = {"sample_index": i, "point": point}
        try:
            if len(coordinates) != len(batch.domain.order.items) or not batch.domain.contains(point):
                raise ValueError("sample is outside the stated domain")
            output = callback(MappingProxyType(point))
            if isinstance(output, ForwardEvaluation):
                record["details"] = dict(output.details)
                if output.status != "optimal":
                    record.update(status=output.status, value=None)
                    records.append(record)
                    continue
                output = output.value
            value = _number(output, "forward scalar output")
            values.append(value)
            record.update(status="ok", value=value)
        except Exception as exc:
            record.update(status="failed", value=None, error_type=type(exc).__name__, error=str(exc))
        records.append(record)
    lower, upper = (min(values), max(values)) if values else (None, None)
    sign = _sign(lower, upper, tolerance) if values else "unidentified"
    failed = len(records) - len(values)
    stability = f"{sign}_on_samples" if values and failed == 0 and batch.status == "complete" else "inconclusive"
    return {
        "estimand": estimand, "unit": unit, "globally_certified": False,
        "scope": "sampled forward evaluations only; never a global nonlinear worst-case certificate",
        "sampled_minimum": lower, "sampled_maximum": upper,
        "successful_sample_sign": sign, "sign_stability": stability, "sign_tolerance": tolerance,
        "accepted": len(batch.points), "successful": len(values), "failed": failed,
        "records": records, "sampling": batch.report(include_points=False),
    }
