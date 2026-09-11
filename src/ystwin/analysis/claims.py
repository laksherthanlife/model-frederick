"""Read-only evaluation of explicitly declared scientific claims and value markers.

Current generated evidence can bind to durable receipts without duplicating their hashes::

    "provenance": {"source": "artifact_registry", "run_id": "product-prediction"}
    "configuration": {"id": "nested-strain-validation", "binding": {
        "source": "artifact_registry", "checks": [
            {"field": "/parameters/mode", "value": "nested"},
            {"source": "main", "field": "validation_mode", "value": "nested"}
        ]}}

The run and every evidence path must belong to ``data/artifact_registry.json``. Adopt a
successful isolated reproduction through ``ystwin.artifacts.adopt_reproduction`` (or its
CLI), with exact-hash review for scientific changes. This reader uses the artifact
validator to check the adopted child-execution proof, accepted bytes, code, inputs, model,
parameters, runtime and semantic contract. It never reproduces, adopts or restamps.
Receipt-backed configuration checks default to JSON pointers into the verified receipt's
``identity``; an explicit evidence source selects output fields instead. Missing proof,
stale proof or mismatched configuration blocks support even if a table cell agrees.

Existing frozen-record pins remain a separate provenance mode. Historical claims and
explicitly retired selectors cannot authorize current support; active selectors must
still have a current surface occurrence regardless of the historical inventory counts.
"""
from __future__ import annotations

import ast
from collections import Counter
from copy import deepcopy
import csv
from decimal import Decimal, InvalidOperation
import fnmatch
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import shlex
import statistics
import tokenize
from typing import Any

from ystwin import artifacts


__all__ = ["ClaimRegistryError", "audit_claims", "load_registry", "selector_identity"]

_ROLES = {"current", "historical", "refused", "investigation_pending"}
_SURFACES = {"markdown", "python_docstrings", "python_field_comments"}
_KINDS = {"numeric", "qualitative", "refusal", "equivalence", "nis", "pending"}
_AGGREGATIONS = {"identity", "count", "nunique", "sum", "mean", "median", "min", "max", "weighted_mean"}
_MARKER_START = re.compile(r"<!--\s*audit:value\b")
_MARKER = re.compile(r"<!--\s*audit:value\b(?P<body>.*?)-->", re.DOTALL)
_LITERAL = re.compile(r"^[+\-−]?(?:(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?|\.\d+)(?:[eE][+\-]?\d+)?")
_REQUIRED = {
    "id", "family", "statement", "role", "method", "conditions", "evidence", "model",
    "configuration", "provenance", "uncertainty", "coverage", "independent_unit", "evaluation",
}


class ClaimRegistryError(ValueError):
    pass


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ClaimRegistryError(message)


def _relative(value: str) -> str:
    _require(isinstance(value, str) and bool(value) and value == value.strip(),
             "a nonempty repository-relative path is required")
    path = PurePosixPath(value)
    _require(not path.is_absolute() and ".." not in path.parts and value != "."
             and path.as_posix() == value and "\\" not in value and ":" not in value,
             f"unsafe repository-relative path: {value!r}")
    return value


def _file(root: Path, relative: str) -> Path:
    path = root / _relative(relative)
    _require(path.resolve().is_relative_to(root), f"evidence escapes repository: {relative}")
    _require(path.is_file(), f"missing evidence or surface: {relative}")
    return path


def _pairs(pairs: list[tuple[str, Any]]) -> dict:
    result = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise ClaimRegistryError(f"nonfinite JSON value: {value}")


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs,
                      parse_constant=_invalid_constant)


def _pointer(value: Any, pointer: str) -> Any:
    _require(isinstance(pointer, str) and (not pointer or pointer.startswith("/")),
             f"invalid JSON pointer: {pointer!r}")
    if not pointer:
        return value
    try:
        for part in pointer.split("/")[1:]:
            key = part.replace("~1", "/").replace("~0", "~")
            if isinstance(value, list):
                _require(key.isdigit() and str(int(key)) == key,
                         f"invalid JSON array index in {pointer}")
                value = value[int(key)]
            else:
                value = value[key]
        return value
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise ClaimRegistryError(f"missing or mistyped field: {pointer}") from exc


def _field(row: dict, field: str) -> Any:
    if field.startswith("/"):
        return _pointer(row, field)
    _require(field in row, f"missing column: {field}")
    return row[field]


def _missing(value: Any) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def _number(value: Any, label: str = "value") -> float:
    _require(not isinstance(value, bool) and not _missing(value), f"missing numeric {label}")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ClaimRegistryError(f"nonnumeric {label}: {value!r}") from exc
    _require(math.isfinite(result), f"nonfinite {label}: {value!r}")
    return result


def _boolean(value: Any) -> bool:
    if type(value) is bool:
        return value
    if value in ("true", "True"):
        return True
    if value in ("false", "False"):
        return False
    raise ClaimRegistryError(f"a declared boolean is required, got {value!r}")


def _key(value: Any) -> tuple[str, str]:
    if _missing(value):
        return "missing", ""
    if type(value) is bool or value in ("true", "True", "false", "False"):
        return "bool", str(_boolean(value))
    text = str(value).strip()
    try:
        number = Decimal(text)
        if number.is_finite():
            return "number", str(number.normalize()) if number else "0"
    except InvalidOperation:
        pass
    return "text", text


def _same(actual: Any, expected: Any) -> bool:
    if isinstance(actual, (dict, list)) or isinstance(expected, (dict, list)):
        return type(actual) is type(expected) and actual == expected
    return _key(actual) == _key(expected)


def selector_identity(path: str, column: str, where: dict, scale: float = 1.0) -> str:
    identity = [path, column, sorted((key, _key(value)) for key, value in where.items()), _number(scale, "selector scale")]
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]


def _marker_unit(specification: str | dict) -> tuple[str, float]:
    if isinstance(specification, str):
        unit, scale = specification, 1.0
    else:
        _require(isinstance(specification, dict), "marker unit declaration must be text or an object")
        unit = specification.get("unit")
        scale = _number(specification.get("scale"), "declared marker scale")
    _require(isinstance(unit, str) and bool(unit) and scale != 0, "marker unit and nonzero scale are required")
    return unit, scale


def _validate_source(source: dict) -> None:
    _require(isinstance(source, dict), "evidence source must be an object")
    for field in ("id", "path", "format", "row_keys", "unit"):
        _require(field in source, f"evidence source missing {field}")
    _relative(source["path"])
    _require(source["format"] in {"csv", "tsv", "json"}, "unsupported evidence format")
    keys = source["row_keys"]
    _require(isinstance(keys, list) and len(keys) == len(set(keys)), "duplicate or invalid row keys")
    _require(all(isinstance(key, str) and key for key in keys), "empty row key")
    _require(bool(keys) or source["format"] == "json", "tabular evidence needs authoritative row keys")
    _require(isinstance(source["unit"], str) and bool(source["unit"]), "evidence unit is required")
    _require(set(source.get("nullable_keys", [])) <= set(keys), "nullable key is not a row key")
    if source["format"] == "json":
        _require("pointer" in source, "JSON evidence needs an explicit pointer")
        _require(isinstance(source["pointer"], str) and
                 (not source["pointer"] or source["pointer"].startswith("/")), "invalid JSON pointer")


def _validate_claim(claim: dict, families: dict) -> None:
    _require(_REQUIRED <= set(claim), f"claim {claim.get('id')!r} missing {sorted(_REQUIRED - set(claim))}")
    _require(isinstance(claim["id"], str) and bool(claim["id"]), "stable claim id is required")
    _require(claim["role"] in _ROLES, f"unknown claim role: {claim['role']}")
    _require(claim["family"] in families, f"undeclared claim family: {claim['family']}")
    _require(isinstance(claim["statement"], str) and bool(claim["statement"].strip()), "statement is required")
    _require(isinstance(claim["method"], dict) and isinstance(claim["method"].get("description"), str)
             and claim["method"]["description"].strip(), "method is required")
    implementations = claim["method"].get("implementation", [])
    _require(isinstance(implementations, list), "method implementation must list repository-relative paths")
    for implementation in implementations:
        _relative(implementation)
    _require(isinstance(claim["conditions"], dict) and bool(claim["conditions"]), "conditions are required")
    _require(isinstance(claim["model"], str) and bool(claim["model"]), "model identity is required")
    for key in ("configuration", "provenance", "uncertainty", "coverage", "independent_unit"):
        _require(isinstance(claim[key], dict) and bool(claim[key]), f"{key} must be declared")
    _require(claim["configuration"].get("id") and "binding" in claim["configuration"],
             "configuration needs an id and an evidence binding")
    provenance = claim["provenance"]
    _require(provenance.get("source") in {None, "artifact_registry"}, "unknown provenance source")
    receipt_backed = provenance.get("source") == "artifact_registry"
    if receipt_backed:
        _require(isinstance(provenance.get("run_id"), str) and provenance["run_id"].strip(),
                 "artifact provenance needs an exact run_id")
        _require(not {"status", "artifacts"} & provenance.keys(),
                 "receipt provenance is derived, not a manual status or hash assertion")
    else:
        _require(provenance.get("status") in {"verified", "pending", "failed"},
                 "provenance must be verified, pending or failed")
        if provenance["status"] == "verified":
            _require(isinstance(provenance.get("artifacts", []), list), "provenance artifacts must be a list")
            for artifact in provenance.get("artifacts", []):
                _require(isinstance(artifact, dict) and "path" in artifact and "sha256" in artifact,
                         "provenance artifact needs a path and SHA-256")
                _relative(artifact["path"])
        else:
            _require(isinstance(provenance.get("reason"), str) and provenance["reason"].strip(),
                     "pending or failed provenance must name its blocker")
    binding = claim["configuration"]["binding"]
    _require(isinstance(binding, dict), "configuration binding must be an object")
    _require(binding.get("source") in {None, "artifact_registry"}, "unknown configuration binding source")
    if binding.get("source") == "artifact_registry":
        _require(receipt_backed and "status" not in binding,
                 "receipt configuration needs receipt provenance, not a manual status")
    else:
        _require(binding.get("status") in {"verified", "pending", "failed"},
                 "configuration binding must declare its verification state")
    checks = binding.get("checks", [])
    _require(isinstance(checks, list), "configuration checks must be a list")
    for check in checks:
        _require(isinstance(check, dict) and isinstance(check.get("field"), str) and check["field"]
                 and "value" in check, "configuration check needs a field and expected value")
        if binding.get("source") == "artifact_registry" and check.get("source", "receipt") == "receipt":
            _require(check["field"].startswith("/"), "receipt configuration fields must be JSON pointers into identity")
    columns = claim["independent_unit"].get("columns", [])
    _require(isinstance(columns, list) and all(isinstance(column, str) and column for column in columns),
             "independent-unit columns must be explicit field names")
    _require(len(columns) == len(set(columns)), "duplicate independent-unit column")
    _require(isinstance(claim["evidence"], list), "evidence must be a list")
    for source in claim["evidence"]:
        _validate_source(source)
    ids = [source["id"] for source in claim["evidence"]]
    _require(len(ids) == len(set(ids)), "duplicate evidence id")
    evaluation = claim["evaluation"]
    _require(isinstance(evaluation, dict) and evaluation.get("kind") in _KINDS, "unknown evaluation kind")
    _require("expected" in claim or "acceptance" in claim, "expected value or acceptance criteria are required")
    kind = evaluation["kind"]
    if kind != "pending":
        _require(bool(ids) and evaluation.get("source", "main") in ids, "evaluation needs declared evidence")
    if kind == "numeric":
        _require(evaluation.get("aggregation") in _AGGREGATIONS, "numeric aggregation must be explicit")
        if "expected" in claim:
            expected = claim["expected"]
            _require(isinstance(expected, dict) and "value" in expected, "numeric expected value is required")
            _number(expected["value"], "expected value")
            for key in ("atol", "rtol"):
                _require(_number(expected.get(key, 0), key) >= 0, f"negative {key}")
    if kind in {"qualitative", "refusal"}:
        predicates = claim.get("acceptance", {}).get("predicates")
        _require(isinstance(predicates, list) and bool(predicates), "qualitative evidence needs predicates")
        for predicate in predicates:
            _require(isinstance(predicate, dict) and isinstance(predicate.get("field"), str),
                     "predicate needs an authoritative field")
            _require(predicate.get("operator", "equals") in {"equals", "in", "missing", "nonempty", "lt", "gt"},
                     "unknown predicate operator")
    _require(claim["role"] != "refused" or kind == "refusal", "refused claims require a refusal evaluator")
    if kind == "pending":
        _require(claim.get("acceptance", {}).get("requires"), "pending evidence must name what is required")


def load_registry(root: str | Path, registry_path: str | Path | None = None) -> dict:
    root = Path(root).resolve()
    requested = Path(registry_path) if registry_path is not None else Path("data/current_claims.json")
    if requested.is_absolute():
        _require(requested.resolve().is_relative_to(root), "registry must be inside the declared repository")
        relative = requested.resolve().relative_to(root).as_posix()
    else:
        relative = requested.as_posix()
    registry = _json(_file(root, relative))
    _require(isinstance(registry, dict) and type(registry.get("schema_version")) is int
             and registry["schema_version"] == 1, "unsupported claim schema")
    families = registry.get("families")
    surfaces = registry.get("surfaces")
    _require(isinstance(families, dict) and bool(families), "scientific family scope is required")
    _require(isinstance(surfaces, list) and bool(surfaces), "surface scope must be explicitly declared")
    surface_ids = set()
    for surface in surfaces:
        _require(isinstance(surface, dict) and surface.get("id") not in surface_ids, "duplicate surface id")
        _require(surface.get("id") and surface.get("reason"), "surface id and scope reason are required")
        surface_ids.add(surface["id"])
        _relative(surface["glob"])
        for excluded in surface.get("exclude", []):
            _relative(excluded)
        _require(surface.get("kind") in _SURFACES, "unsupported claim surface")
        _require(surface.get("role") in _ROLES | {"forbidden"}, "surface role must be explicit")
        _require(surface.get("unmarked_policy") == "registered_only", "unmarked text must not be inferred from numbers")
    for name, family in families.items():
        _require(isinstance(family, dict) and family.get("boundary"), f"family {name} has no claim boundary")
        _require(isinstance(family.get("locations"), list) and bool(family["locations"]),
                 f"family {name} needs declared surface locations")
        for location in family["locations"]:
            _relative(location["path"])
            _require(location.get("scope") in surface_ids and location.get("anchor"),
                     f"family {name} has an undeclared surface or anchor")
    profiles = registry.get("profiles", {})
    _require(isinstance(profiles, dict), "profiles must be an object")
    claims, ids = [], set()
    for item in registry.get("claims", []):
        _require(isinstance(item, dict), "claim must be an object")
        profile = item.get("profile")
        _require(profile is None or profile in profiles, f"unknown claim profile: {profile}")
        claim = {**deepcopy(profiles.get(profile, {})), **deepcopy(item)}
        _validate_claim(claim, families)
        _require(claim["id"] not in ids, f"duplicate claim id: {claim['id']}")
        ids.add(claim["id"])
        claims.append(claim)
    _require(bool(claims), "claim inventory must not be empty")
    _require(set(families) <= {claim["family"] for claim in claims}, "declared family has no inventory records")
    registry["claims"] = claims
    marker_inventory = registry.get("marked_values", {})
    _require(isinstance(marker_inventory, dict), "marked-value inventory must be an object")
    active_selectors = set()
    for path, table in marker_inventory.get("tables", {}).items():
        _relative(path)
        _require(table.get("profile") in profiles, f"unbound marker table: {path}")
        _require(isinstance(table.get("selectors"), list) and bool(table["selectors"]),
                 f"marker table {path} needs authoritative selectors")
        seen = set()
        for selector in table["selectors"]:
            _require(isinstance(selector.get("where"), dict) and bool(selector["where"]), "empty marker selector")
            _require(isinstance(selector.get("columns"), dict) and bool(selector["columns"]), "marker units are required by column")
            for column, specification in selector["columns"].items():
                _, scale = _marker_unit(specification)
                key = selector_identity(path, column, selector["where"], scale)
                _require(key not in seen, f"duplicate marker selector: {path} {column}")
                seen.add(key)
                active_selectors.add(key)
    retired = marker_inventory.get("retired_selectors", [])
    _require(isinstance(retired, list), "retired selectors must be an explicit list")
    retired_ids = set()
    current_claims = {claim["id"]: claim for claim in claims if claim["role"] in {"current", "refused"}}
    for entry in retired:
        _require(isinstance(entry, dict), "retired selector must be an object")
        _relative(entry["path"])
        _require(isinstance(entry.get("reason"), str) and entry["reason"].strip(),
                 "selector retirement needs an explicit reason")
        family = entry.get("replacement_family")
        _require(family in families, "selector retirement needs a declared replacement family")
        replacements = entry.get("replacement_claims")
        _require(isinstance(replacements, list) and replacements and
                 all(identity in current_claims and current_claims[identity]["family"] == family
                     for identity in replacements), "selector retirement needs current replacement coverage")
        _require(isinstance(entry.get("where"), dict) and entry["where"], "retired selector needs its row identity")
        _require(isinstance(entry.get("columns"), dict) and entry["columns"], "retired selector needs its columns and units")
        for column, specification in entry["columns"].items():
            _, scale = _marker_unit(specification)
            key = selector_identity(entry["path"], column, entry["where"], scale)
            _require(key not in active_selectors, "a selector cannot be both active and retired")
            _require(key not in retired_ids, "duplicate retired selector")
            retired_ids.add(key)
    return registry


class _Evidence:
    def __init__(self, root: Path):
        self.root = root
        self.tables: dict[tuple, list[dict]] = {}
        self.hashes: dict[str, str] = {}
        self.artifact_registry: dict | None = None
        self.run_identities: dict[str, dict] = {}
        self.receipts: dict[tuple[str, str], dict] = {}

    def receipt(self, claim: dict) -> dict:
        """Validate adopted per-artifact proofs with the same validator as the artifact gate."""
        if self.artifact_registry is None:
            self.artifact_registry = artifacts.load_registry(self.root)
        registry = self.artifact_registry
        run_id = claim["provenance"]["run_id"]
        _require(run_id in registry["runs"], f"unregistered artifact run: {run_id}")
        run = registry["runs"][run_id]
        _require(run["lifecycle"] == "current", f"artifact run is not current: {run_id}")
        references = []
        for source in claim["evidence"]:
            path = source["path"]
            item = registry["artifacts"].get(path)
            _require(item is not None and item["run"] == run_id and item["lifecycle"] == "current",
                     f"evidence is not a current member of {run_id}: {path}")
            _require(isinstance(item.get("comparison"), dict), f"semantic artifact contract absent: {path}")
            reference = (run.get("verification") or {}).get("artifacts", {}).get(path)
            _require(isinstance(reference, dict) and reference.get("receipt_id"),
                     f"durable reproduction receipt absent: {path}; reproduce and explicitly adopt {run_id}, do not restamp")
            key = run_id, path
            if key not in self.receipts:
                if run_id not in self.run_identities:
                    self.run_identities[run_id] = artifacts.run_identity(self.root, run)
                issues = artifacts._receipt_issues(self.root, registry, run, item, self.run_identities[run_id])
                _require(not issues, f"artifact receipt invalid for {path}: {'; '.join(issues)}")
                receipt = registry["receipts"][reference["receipt_id"]]
                _require(receipt.get("baseline_policy") == registry["document"].get("baseline_ref"),
                         f"artifact receipt baseline policy changed: {path}")
                self.receipts[key] = receipt
            references.append({"path": path, **reference})
        _require(bool(references), "receipt provenance needs evidence artifacts")
        identity = self.receipts[run_id, references[0]["path"]]["identity"]
        code = {entry["path"] for entry in identity["code"]}
        for implementation in claim["method"].get("implementation", []):
            _require(implementation in code, f"unbound producer provenance: {implementation}")
        return {"source": "artifact_registry", "registry": artifacts.REGISTRY_PATH, "run_id": run_id,
                "status": "verified", "artifacts": references, "identity": identity}

    def digest(self, relative: str) -> str:
        if relative not in self.hashes:
            self.hashes[relative] = hashlib.sha256(_file(self.root, relative).read_bytes()).hexdigest()
        return self.hashes[relative]

    def read(self, source: dict) -> list[dict]:
        identity = source["path"], source["format"], source.get("pointer", "")
        if identity not in self.tables:
            path = _file(self.root, source["path"])
            if source["format"] == "json":
                value = _pointer(_json(path), source["pointer"])
                rows = value if isinstance(value, list) else [value]
                _require(all(isinstance(row, dict) for row in rows), "JSON table must contain objects")
            else:
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t" if source["format"] == "tsv" else ",")
                    fields = reader.fieldnames
                    _require(fields and all(fields) and len(fields) == len(set(fields)),
                             f"empty or duplicate columns in {source['path']}")
                    rows = list(reader)
                _require(all(None not in row and None not in row.values() for row in rows),
                         f"malformed row in {source['path']}")
            _require(bool(rows), f"empty evidence table: {source['path']}")
            self.tables[identity] = rows
        rows = self.tables[identity]
        keys = source["row_keys"]
        _require(keys or len(rows) == 1, "JSON arrays need unique authoritative row keys")
        nullable = set(source.get("nullable_keys", []))
        seen = set()
        for row in rows:
            values = [_field(row, key) for key in keys]
            _require(all(not _missing(value) or key in nullable for key, value in zip(keys, values)),
                     f"missing row key in {source['path']}")
            key = tuple(_key(value) for value in values)
            _require(key not in seen, f"duplicate row key in {source['path']}: {dict(zip(keys, values))}")
            seen.add(key)
        return rows

    def select(self, source: dict, where: dict) -> list[dict]:
        rows = self.read(source)
        _require(isinstance(where, dict), "row selector must be an object")
        selected = [row for row in rows if all(_same(_field(row, key), value) for key, value in where.items())]
        _require(bool(selected), f"selector matched no rows of {source['path']}: {where}")
        return selected


def _scope_matches(surface: dict, relative: str) -> bool:
    def matches(pattern: str) -> bool:
        return fnmatch.fnmatchcase(relative, pattern) or (
            "**/" in pattern and fnmatch.fnmatchcase(relative, pattern.replace("**/", "")))
    return matches(surface["glob"]) and not any(matches(pattern) for pattern in surface.get("exclude", []))


def _segments(root: Path, surface: dict):
    paths = sorted(root.glob(surface["glob"]))
    matched = 0
    for path in paths:
        relative = path.relative_to(root).as_posix()
        if not path.is_file() or not _scope_matches(surface, relative):
            continue
        matched += 1
        source = _file(root, relative).read_text(encoding="utf-8")
        if surface["kind"] == "markdown":
            yield relative, "document", source, 1
            continue
        tree = ast.parse(source, filename=relative)
        if surface["kind"] == "python_field_comments":
            for token in tokenize.generate_tokens(io.StringIO(source).readline):
                if token.type == tokenize.COMMENT and token.string.startswith("#:"):
                    yield relative, "field", token.string, token.start[0]
            continue

        def walk(node, parents):
            name = ".".join(parents) or "<module>"
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                body = node.body
                if (body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)
                        and isinstance(body[0].value.value, str)):
                    first = body[0].value
                    yield relative, name, ast.get_source_segment(source, first), first.lineno
            for child in ast.iter_child_nodes(node):
                names = parents + [child.name] if isinstance(child, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) else parents
                yield from walk(child, names)

        yield from walk(tree, [])
    _require(matched or not surface.get("required", False), f"declared surface is absent: {surface['glob']}")


def _markers(text: str):
    marks = list(_MARKER.finditer(text))
    _require(len(_MARKER_START.findall(text)) == len(marks), "malformed or incomplete audit:value marker")
    for index, mark in enumerate(marks):
        fields = {}
        for token in shlex.split(mark.group("body")):
            key, separator, value = token.partition("=")
            _require(separator and value and key in {"table", "column", "row", "scale"}, "invalid marker field")
            _require(key not in fields, f"duplicate marker field: {key}")
            fields[key] = value
        _require({"table", "column", "row"} <= set(fields), "marker needs table, column and row")
        _relative(fields["table"])
        where = {}
        for clause in fields["row"].split(";"):
            if not clause.strip():
                continue
            key, separator, value = clause.partition("=")
            key, value = key.strip(), value.strip()
            _require(separator and key and value, "invalid marker row selector")
            _require(key not in where, f"duplicate marker row key: {key}")
            where[key] = value
        _require(bool(where), "empty marker row selector")
        line_end = text.find("\n", mark.end())
        end = min(line_end if line_end >= 0 else len(text),
                  marks[index + 1].start() if index + 1 < len(marks) else len(text))
        literal = _LITERAL.match(text[mark.end():end].lstrip(" \t*`_"))
        _require(literal is not None, "marker must be followed by its numeric literal on the same line")
        written = literal.group().replace(",", "").replace("−", "-")
        number = Decimal(written)
        yield {
            "path": fields["table"], "column": fields["column"], "where": where,
            "scale": _number(fields.get("scale", 1), "marker scale"),
            "expected": {"value": _number(written), "atol": float(Decimal("0.5") * Decimal(10) ** number.as_tuple().exponent)},
            "literal": literal.group(), "line": text.count("\n", 0, mark.start()) + 1,
        }


def _numeric_acceptance(actual: float, claim: dict) -> bool:
    if "expected" in claim:
        expected = claim["expected"]
        return math.isclose(actual, _number(expected["value"]),
                            abs_tol=_number(expected.get("atol", 0)) * (1 + 1e-12),
                            rel_tol=_number(expected.get("rtol", 0)))
    acceptance = claim["acceptance"]
    _require("minimum" in acceptance or "maximum" in acceptance, "numeric acceptance needs a bound")
    return (("minimum" not in acceptance or actual >= _number(acceptance["minimum"])) and
            ("maximum" not in acceptance or actual <= _number(acceptance["maximum"])))


def _row_key(source: dict, row: dict) -> dict:
    return {key: _field(row, key) for key in source["row_keys"]} or {"json_pointer": source["pointer"]}


def _coverage(rows: list[dict], claim: dict) -> dict:
    evaluation = claim["evaluation"]
    coverage = claim["coverage"]
    expected = evaluation.get("expected_rows", coverage.get("expected_rows"))
    if expected is not None:
        _require(type(expected) is int and expected > 0, "expected row count must be positive")
        _require(len(rows) == expected, f"coverage mismatch: expected {expected} rows, found {len(rows)}")
    expected_keys = coverage.get("expected_keys")
    if expected_keys is not None:
        _require(isinstance(expected_keys, list) and bool(expected_keys) and
                 all(isinstance(key, dict) and key for key in expected_keys), "invalid expected population keys")
        columns = sorted(expected_keys[0])
        _require(all(sorted(key) == columns for key in expected_keys), "inconsistent expected population keys")
        expected_set = {tuple(_key(key[column]) for column in columns) for key in expected_keys}
        _require(len(expected_set) == len(expected_keys), "duplicate expected population key")
        observed_set = {tuple(_key(_field(row, column)) for column in columns) for row in rows}
        _require(observed_set == expected_set, "coverage mismatch: selected population differs from declared row identities")
    unit = claim["independent_unit"]
    columns = unit.get("columns", [])
    count = None
    if columns:
        units = set()
        for row in rows:
            values = [_field(row, column) for column in columns]
            _require(not any(_missing(value) for value in values), "missing independent unit")
            units.add(tuple(_key(value) for value in values))
        count = len(units)
        if "expected_count" in unit:
            _require(count == unit["expected_count"], f"independent-unit coverage mismatch: expected {unit['expected_count']}, found {count}")
    return {"selected_rows": len(rows), "independent_unit": unit.get("name"), "independent_units": count}


def _numeric(rows: list[dict], source: dict, claim: dict):
    evaluation = claim["evaluation"]
    aggregation = evaluation["aggregation"]
    values, results = [], []
    for row in rows:
        record = {"source": source["path"], "row_key": _row_key(source, row),
                  "unit": evaluation.get("unit", source["unit"]),
                  "input_unit": evaluation.get("input_unit", source["unit"]),
                  "context": {field: _field(row, field) for field in source.get("context_columns", [])}}
        if aggregation == "count":
            value = 1.0
        elif aggregation == "nunique":
            value = _field(row, evaluation["column"])
            _require(not _missing(value), "missing distinct value")
        else:
            value = _number(_field(row, evaluation["column"]), evaluation["column"])
            record["observed"] = value
            transform = evaluation.get("transform")
            if transform:
                reference = _number(_field(row, transform["reference_column"]), "reference")
                _require(reference > 0, "relative comparisons need a positive reference")
                record["reference"] = reference
                ratio = value / reference
                if transform["kind"] == "absolute_relative_error":
                    value = abs(ratio - 1)
                elif transform["kind"] == "symmetric_fold_error":
                    _require(ratio > 0, "fold error needs a positive prediction")
                    value = max(ratio, 1 / ratio)
                elif transform["kind"] == "ratio":
                    value = ratio
                else:
                    raise ClaimRegistryError("unknown numeric transform")
            value *= _number(evaluation.get("scale", 1), "scale")
        record["value"] = value
        results.append(record)
        values.append(value)
    if aggregation == "identity":
        _require(len(values) == 1, f"single-cell selector matched {len(values)} rows")
        actual = values[0]
    elif aggregation == "nunique":
        actual = len({_key(value) for value in values})
    elif aggregation == "weighted_mean":
        weights = [_number(_field(row, evaluation["weight_column"]), "weight") for row in rows]
        _require(all(weight >= 0 for weight in weights) and sum(weights) > 0, "invalid aggregation weights")
        actual = sum(value * weight for value, weight in zip(values, weights)) / sum(weights)
    else:
        aggregate = {"count": len, "sum": sum, "mean": statistics.mean, "median": statistics.median,
                     "min": min, "max": max}[aggregation]
        actual = aggregate(values)
    actual = _number(actual, "aggregate")
    return actual, _numeric_acceptance(actual, claim), results


def _predicates(rows: list[dict], source: dict, claim: dict):
    results, passed = [], True
    for row in rows:
        for predicate in claim["acceptance"]["predicates"]:
            actual = _field(row, predicate["field"])
            operator = predicate.get("operator", "equals")
            expected = predicate.get("value")
            if operator == "missing":
                ok = _missing(actual)
            elif operator == "nonempty":
                ok = not _missing(actual)
            elif operator == "in":
                _require(isinstance(expected, list) and expected, "in predicate needs declared alternatives")
                ok = any(_same(actual, value) for value in expected)
            elif operator == "equals":
                _require(not _missing(actual) and "value" in predicate, "missing qualitative evidence")
                ok = _same(actual, expected)
            elif operator == "lt":
                ok = _number(actual) < _number(expected)
            else:
                ok = _number(actual) > _number(expected)
            results.append({"source": source["path"], "row_key": _row_key(source, row),
                            "field": predicate["field"], "value": actual, "expected": expected,
                            "predicate_passed": ok, "unit": source["unit"]})
            passed &= ok
    return {"predicates_satisfied": passed, "checked": len(results)}, passed, results


def _equivalence(rows: list[dict], source: dict, claim: dict):
    acceptance = claim["acceptance"]
    margins = acceptance.get("margins")
    _require(isinstance(margins, list) and len(margins) == 2 and acceptance.get("margin_source"),
             "equivalence requires two declared margins and their scientific justification")
    lower_margin, upper_margin = map(_number, margins)
    null = _number(acceptance.get("null", 1))
    _require(lower_margin < null < upper_margin, "equivalence margins must bracket the null")
    confidence = _number(acceptance.get("confidence"), "interval confidence")
    _require(0 < confidence < 1, "invalid interval confidence")
    evaluation = claim["evaluation"]
    results = []
    for row in rows:
        if evaluation.get("estimable_column"):
            _require(_boolean(_field(row, evaluation["estimable_column"])), "interval was refused as unestimable")
        lower = _number(_field(row, evaluation["lower_column"]), "interval lower limit")
        upper = _number(_field(row, evaluation["upper_column"]), "interval upper limit")
        _require(lower <= upper, "reversed uncertainty interval")
        if evaluation.get("value_column"):
            estimate = _number(_field(row, evaluation["value_column"]), "point estimate")
            _require(lower <= estimate <= upper, "point estimate is outside its uncertainty interval")
        equivalent = lower_margin < lower and upper < upper_margin
        results.append({"source": source["path"], "row_key": _row_key(source, row), "unit": source["unit"],
                        "low": lower, "high": upper, "margins": margins, "confidence": confidence,
                        "equivalent": equivalent})
    passed = all(row["equivalent"] for row in results)
    return {"equivalent": passed, "margins": margins, "intervals": len(results)}, passed, results


def _nis(rows: list[dict], source: dict, claim: dict):
    from scipy.stats import chi2

    evaluation = claim["evaluation"]
    alpha = _number(evaluation.get("alpha"), "NIS alpha")
    ess_floor = _number(evaluation.get("ess_floor"), "ESS floor")
    white_threshold = _number(evaluation.get("white_fraction_threshold"), "whiteness fraction threshold")
    _require(0 < alpha < 1 and ess_floor > 0 and 0 < white_threshold <= 1, "invalid NIS criteria")
    results, tested = [], []
    for row in rows:
        ess = _number(_field(row, "min_ess"), "minimum ESS")
        _require(ess > 0, "ESS must be positive")
        eligible = ess >= ess_floor
        _require(_boolean(_field(row, "tested")) == eligible, "tested flag disagrees with declared ESS floor")
        record = {"source": source["path"], "row_key": _row_key(source, row), "unit": source["unit"],
                  "tested": eligible, "min_ess": ess}
        if eligible:
            steps = _number(_field(row, "k_steps"), "NIS steps")
            _require(steps >= 3 and steps.is_integer(), "NIS needs at least three integer steps per well")
            mean_nis = _number(_field(row, "mean_nis"), "mean NIS")
            correlation = _number(_field(row, "lag1_autocorr"), "lag-one autocorrelation")
            _require(mean_nis >= 0 and abs(correlation) <= 1, "invalid NIS or autocorrelation")
            low, high = (_number(chi2.ppf(alpha / 2, steps) / steps, "NIS lower limit"),
                         _number(chi2.ppf(1 - alpha / 2, steps) / steps, "NIS upper limit"))
            bound = 2 / math.sqrt(steps)
            record.update({"k_steps": int(steps), "mean_nis": mean_nis, "lag1_autocorr": correlation,
                           "band_low": low, "band_high": high, "white_bound": bound,
                           "in_band": low <= mean_nis <= high, "above_band": mean_nis > high,
                           "below_band": mean_nis < low, "white": abs(correlation) <= bound})
            tested.append(record)
        results.append(record)
    acceptance = claim["acceptance"]
    _require(set(acceptance) >= {"nis_verdict", "whiteness_verdict"}, "NIS and whiteness must be assessed together")
    if not tested:
        return {"tested_wells": 0, "inconclusive_wells": len(rows), "median_mean_nis": None,
                "median_lag1_autocorr": None, "nis_verdict": "inconclusive",
                "whiteness_verdict": "inconclusive"}, False, results

    def fraction(key):
        return sum(row[key] for row in tested) / len(tested)

    white_fraction = fraction("white")
    actual = {
        "tested_wells": len(tested), "inconclusive_wells": len(rows) - len(tested),
        "median_mean_nis": statistics.median(row["mean_nis"] for row in tested),
        "median_lag1_autocorr": statistics.median(row["lag1_autocorr"] for row in tested),
        "fraction_in_band": fraction("in_band"), "fraction_white": white_fraction,
        "nis_verdict": "over_confident" if fraction("above_band") > 0.5 else
                       "over_dispersed" if fraction("below_band") > 0.5 else "mixed",
        "whiteness_verdict": "not_white" if white_fraction < white_threshold else "mostly_white",
        "interpretation": "Per-well nominal chi-square diagnostics; serial dependence prevents an independent-draw calibration claim.",
    }
    passed = all(actual[key] == acceptance[key] for key in ("nis_verdict", "whiteness_verdict"))
    return actual, passed, results


def _trust(claim: dict, evidence: _Evidence) -> tuple[list[str], dict | None]:
    errors, proof = [], None
    provenance = claim["provenance"]
    if provenance.get("source") == "artifact_registry":
        try:
            proof = evidence.receipt(claim)
        except (OSError, ValueError, KeyError, TypeError, SyntaxError) as exc:
            errors.append(f"provenance unverified: {exc}")
    elif provenance["status"] != "verified":
        errors.append(f"provenance {provenance['status']}: {provenance.get('reason', 'verified regeneration is required')}")
    else:
        pinned = provenance.get("artifacts", [])
        if not pinned:
            errors.append("verified provenance requires pre-existing artifact identities, not a status assertion")
        paths = set()
        for artifact in pinned:
            try:
                expected = artifact["sha256"]
                _require(isinstance(expected, str) and len(expected) == 64 and
                         set(expected) <= set("0123456789abcdef"), "invalid provenance SHA-256")
                _require(artifact["path"] not in paths, "duplicate provenance artifact")
                paths.add(artifact["path"])
                _require(evidence.digest(artifact["path"]) == expected,
                         f"provenance mismatch: {artifact['path']}; regenerate and review, do not restamp")
            except (ClaimRegistryError, KeyError, OSError, ValueError) as exc:
                errors.append(str(exc))
        for source in claim["evidence"]:
            if source["path"] not in paths:
                errors.append(f"unbound evidence provenance: {source['path']}")
        if provenance.get("scope") == "reproduced_result":
            for implementation in claim["method"].get("implementation", []):
                if implementation not in paths:
                    errors.append(f"unbound producer provenance: {implementation}")
        elif provenance.get("scope") != "frozen_record":
            errors.append("provenance scope must distinguish a reproduced result from a frozen record")
    binding = claim["configuration"]["binding"]
    receipt_backed = binding.get("source") == "artifact_registry"
    if not receipt_backed and binding.get("status") != "verified":
        errors.append(f"configuration unverified: {binding.get('reason', 'missing configuration binding')}")
    else:
        checks = binding.get("checks", [])
        if not checks:
            errors.append("configuration needs authoritative field checks")
        if receipt_backed and not any(check.get("source", "receipt") == "receipt" for check in checks):
            errors.append("configuration needs checks against the verified receipt identity")
        for check in checks:
            try:
                source_id = check.get("source", "receipt" if receipt_backed else "main")
                if receipt_backed and source_id == "receipt":
                    _require(proof is not None, "configuration unverified: no valid artifact receipt")
                    rows = [proof["identity"]]
                else:
                    sources = {source["id"]: source for source in claim["evidence"]}
                    rows = evidence.select(sources[source_id], check.get("where", {}))
                _require(all(_same(_field(row, check["field"]), check["value"]) for row in rows),
                         f"configuration mismatch: {check['field']} != {check['value']!r}")
            except (ClaimRegistryError, KeyError, OSError, ValueError, TypeError) as exc:
                errors.append(str(exc))
    return list(dict.fromkeys(errors)), proof


def _record(claim: dict) -> dict:
    return {
        "id": claim["id"], "check": f"structured claim: {claim['id']}", "record_type": "claim",
        **{key: deepcopy(claim[key]) for key in ("family", "statement", "role", "method", "conditions", "evidence",
                                                 "model", "configuration", "provenance", "uncertainty", "coverage", "independent_unit")},
        "expected": deepcopy(claim.get("expected", claim.get("acceptance"))),
        "evaluation": deepcopy(claim["evaluation"]), "actual": None, "observed": None,
        "status": "FAIL", "verdict": "BLOCKED", "supported": False, "blocking": True,
        "comparison_status": "NOT_EVALUATED", "result_rows": [], "detail": "",
    }


def _audit_one(claim: dict, evidence: _Evidence) -> dict:
    record = _record(claim)
    if claim["role"] == "historical":
        record.update(status="SKIP", verdict="HISTORICAL", blocking=False,
                      detail=claim.get("withdrawal_reason", "Retained historical context; not asserted as a current scientific result."))
        return record
    if claim["evaluation"]["kind"] == "pending":
        record["detail"] = str(claim["acceptance"]["requires"])
        return record
    problems, proof = _trust(claim, evidence)
    if proof is not None:
        record["provenance_observed"] = proof
    if claim["role"] == "investigation_pending":
        problems.append("investigation pending: acceptance and regeneration are not yet authorized")
    try:
        evaluation = claim["evaluation"]
        source = next(source for source in claim["evidence"] if source["id"] == evaluation.get("source", "main"))
        rows = evidence.select(source, evaluation.get("where", {}))
        record["coverage_observed"] = _coverage(rows, claim)
        evaluators = {"numeric": _numeric, "qualitative": _predicates, "refusal": _predicates,
                      "equivalence": _equivalence, "nis": _nis}
        actual, passed, results = evaluators[evaluation["kind"]](rows, source, claim)
        if evaluation["kind"] == "nis" and actual["nis_verdict"] == "inconclusive":
            problems.append("NIS inconclusive: no well meets the declared ESS floor")
        record["observed"] = actual
        record["comparison_status"] = "PASS" if passed else "FAIL"
        record["result_rows"] = [dict(row, evidence_status="UNVERIFIED" if problems else "VERIFIED") for row in results]
        if problems:
            record["detail"] = "; ".join(problems)
        else:
            record.update(actual=actual, status="PASS" if passed else "FAIL", supported=passed,
                          blocking=not passed, verdict=("REFUSAL_CONFIRMED" if claim["role"] == "refused" else "SUPPORTED")
                          if passed else "NOT_SUPPORTED",
                          detail="Acceptance criteria satisfied within the declared scope." if passed else
                                 "Evidence does not satisfy the declared acceptance criteria.")
    except (ClaimRegistryError, OSError, ValueError, KeyError, TypeError, StopIteration, OverflowError) as exc:
        record["detail"] = "; ".join([*problems, str(exc)])
    return record


def _failure(identity: str, detail: str, record_type: str = "scope") -> dict:
    return {"id": identity, "check": identity, "record_type": record_type, "status": "FAIL",
            "verdict": "BLOCKED", "supported": False, "blocking": True, "actual": None,
            "observed": None, "expected": "valid declared scope and evidence", "detail": detail,
            "result_rows": [], "comparison_status": "NOT_EVALUATED"}


def _audit_markers(root: Path, registry: dict, evidence: _Evidence) -> list[dict]:
    inventory = registry.get("marked_values", {})
    tables = inventory.get("tables", {})
    declarations = {}
    for path, table in tables.items():
        for selector in table["selectors"]:
            for column, specification in selector["columns"].items():
                unit, scale = _marker_unit(specification)
                key = selector_identity(path, column, selector["where"], scale)
                declarations[key] = table, unit
    retired = {}
    for entry in inventory.get("retired_selectors", []):
        for column, specification in entry["columns"].items():
            unit, scale = _marker_unit(specification)
            key = selector_identity(entry["path"], column, entry["where"], scale)
            retired[key] = {"path": entry["path"], "column": column, "where": entry["where"],
                            "scale": scale, "unit": unit, "reason": entry["reason"],
                            "replacement_family": entry["replacement_family"],
                            "replacement_claims": entry["replacement_claims"]}
    rows, seen_selectors, seen_surfaces = [], set(), set()
    occurrences = Counter()
    for surface in registry["surfaces"]:
        try:
            for relative, anchor, text, first_line in _segments(root, surface):
                identity = relative, surface["kind"], anchor, first_line
                _require(identity not in seen_surfaces, f"overlapping surface declarations: {relative} {anchor}")
                seen_surfaces.add(identity)
                if not _MARKER_START.search(text):
                    continue
                for marker in _markers(text):
                    key = selector_identity(marker["path"], marker["column"], marker["where"], marker["scale"])
                    if surface["role"] not in {"historical", "forbidden"}:
                        seen_selectors.add(key)
                    occurrence = relative, anchor, key
                    occurrences[occurrence] += 1
                    digest = hashlib.sha256(json.dumps(occurrence).encode()).hexdigest()[:20]
                    claim_id = f"marked.{digest}.{occurrences[occurrence]}"
                    location = {"path": relative, "anchor": anchor, "line": first_line + marker["line"] - 1,
                                "scope": surface["id"]}
                    if surface["role"] == "forbidden":
                        result = _failure(claim_id, f"marker on forbidden surface: {relative}:{location['line']}", "marked_value")
                    elif surface["role"] == "historical" and key not in declarations:
                        result = _failure(claim_id, "Historical marker lies outside the current scientific inventory.", "marked_value")
                        result.update(role="historical", status="SKIP", verdict="HISTORICAL", blocking=False,
                                      expected=marker["expected"])
                    elif key in retired:
                        result = _failure(claim_id, f"retired selector reintroduced as current: {retired[key]['reason']}", "marked_value")
                    elif key not in declarations:
                        result = _failure(claim_id, f"unregistered authoritative selector: {marker['path']} {marker['column']} {marker['where']}", "marked_value")
                    else:
                        table, unit = declarations[key]
                        claim = deepcopy(registry["profiles"][table["profile"]])
                        claim.update({"id": claim_id, "role": surface["role"], "statement":
                                      f"Marked {marker['column']} = {marker['literal']} under its declared table selector.",
                                      "expected": marker["expected"], "evaluation": {
                                          "kind": "numeric", "source": "main", "where": marker["where"],
                                          "column": marker["column"], "scale": marker["scale"],
                                          "aggregation": "identity", "expected_rows": 1}})
                        claim["independent_unit"] = {"name": "one cited table row; not an independent experiment"}
                        claim["evidence"][0]["unit"] = unit
                        _require(claim["evidence"][0]["path"] == marker["path"], "marker profile points at another table")
                        _validate_claim(claim, registry["families"])
                        result = _audit_one(claim, evidence)
                        result["record_type"] = "marked_value"
                    result.update(selector_id=key, location=location,
                                  check=f"marked value: {relative}:{location['line']} ({marker['path']}[{marker['column']}])",
                                  selector={"path": marker["path"], "column": marker["column"],
                                            "where": marker["where"], "scale": marker["scale"]})
                    rows.append(result)
        except (ClaimRegistryError, OSError, SyntaxError, ValueError, KeyError, TypeError, tokenize.TokenError) as exc:
            rows.append(_failure(f"scope.{surface['id']}", str(exc)))
    for key in sorted(set(declarations) - seen_selectors):
        rows.append(_failure(f"marker_inventory.{key}", "registered selector has no remaining surface marker; review the declaration rather than silently shrinking coverage"))
    for key, entry in sorted(retired.items()):
        rows.append({"id": f"retired_selector.{key}", "check": f"retired selector: {key}",
                     "record_type": "marker_retirement", "role": "historical", "status": "SKIP",
                     "verdict": "RETIRED", "supported": False, "blocking": False, "actual": None,
                     "expected": "no current occurrence", "selector_id": key, "selector": entry,
                     "detail": f"{entry['reason']} Replacement coverage: {', '.join(entry['replacement_claims'])}.",
                     "result_rows": []})
    complete = not any(row["blocking"] and (row["record_type"] == "scope" or "family" not in row) for row in rows)
    summary = {
        "id": "registry.marker_inventory", "check": "declared audit:value inventory", "record_type": "scope",
        "status": "PASS" if complete else "FAIL", "verdict": "INVENTORIED" if complete else "BLOCKED",
        "supported": False, "blocking": not complete,
        "actual": {"occurrences": sum(occurrences.values()), "distinct_selectors": len(seen_selectors)},
        "expected": {"baseline_occurrences": inventory.get("baseline_occurrences"),
                     "baseline_distinct_selectors": inventory.get("baseline_distinct_selectors"),
                     "declared_selectors": len(declarations), "retired_selectors": len(retired)},
        "detail": "Active selectors require current occurrences. Baseline counts are historical; neither counts nor cell agreement establish scientific support.",
        "result_rows": [],
    }
    rows.append(summary)
    return rows


def audit_claims(root: str | Path, registry_path: str | Path | None = None) -> list[dict]:
    root = Path(root).resolve()
    try:
        registry = load_registry(root, registry_path)
    except (ClaimRegistryError, OSError, ValueError, KeyError, TypeError) as exc:
        return [_failure("registry.schema", str(exc), "registry")]
    evidence = _Evidence(root)
    rows = []
    scopes = {surface["id"]: surface for surface in registry["surfaces"]}
    for family, declaration in registry["families"].items():
        for index, location in enumerate(declaration["locations"]):
            try:
                _file(root, location["path"])
                _require(_scope_matches(scopes[location["scope"]], location["path"]),
                         "family location is outside its declared surface")
            except (ClaimRegistryError, OSError) as exc:
                rows.append(_failure(f"scope.family.{family}.{index}", str(exc)))
    rows.extend(_audit_one(claim, evidence) for claim in registry["claims"])
    rows.extend(_audit_markers(root, registry, evidence))
    ids = [row["id"] for row in rows]
    if len(ids) != len(set(ids)):
        rows.append(_failure("registry.result_ids", "duplicate audit record identities"))
    return rows
