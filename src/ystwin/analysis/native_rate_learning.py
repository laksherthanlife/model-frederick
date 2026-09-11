from __future__ import annotations

import copy
import hashlib
import itertools
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass, fields
from enum import Enum
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as f_distribution
from scipy.stats import t as t_distribution


class ContractError(ValueError):
    pass


class PhysicalType(str, Enum):
    ENZYME_CONCENTRATION = "enzyme_concentration"
    ENZYME_AMOUNT_PER_DRY_MASS = "enzyme_amount_per_dry_mass"
    RELATIVE_ENZYME = "relative_enzyme"
    METABOLITE_CONCENTRATION = "metabolite_concentration"
    RELATIVE_METABOLITE = "relative_metabolite"
    VOLUMETRIC_FLUX = "volumetric_flux"
    SPECIFIC_FLUX = "specific_flux"
    RELATIVE_FLUX = "relative_flux"


_CONCENTRATION = {
    "mol/L": 1.0, "mmol/L": 1e-3, "umol/L": 1e-6, "nmol/L": 1e-9,
    "M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9,
}
_AMOUNT = {f"{prefix}mol/gDW": scale for prefix, scale in (("", 1.0), ("m", 1e-3), ("u", 1e-6), ("n", 1e-9), ("p", 1e-12))}
_FLUX = {
    basis: {
        f"{prefix}mol/({basis}*{time})": scale / seconds
        for prefix, scale in (("", 1.0), ("m", 1e-3), ("u", 1e-6), ("n", 1e-9))
        for time, seconds in (("s", 1.0), ("min", 60.0), ("h", 3600.0))
    }
    for basis in ("L", "gDW")
}
_UNIT_TABLE = {
    PhysicalType.ENZYME_CONCENTRATION: ("mol/L", _CONCENTRATION),
    PhysicalType.ENZYME_AMOUNT_PER_DRY_MASS: ("mol/gDW", _AMOUNT),
    PhysicalType.METABOLITE_CONCENTRATION: ("mol/L", _CONCENTRATION),
    PhysicalType.VOLUMETRIC_FLUX: ("mol/(L*s)", _FLUX["L"]),
    PhysicalType.SPECIFIC_FLUX: ("mol/(gDW*s)", _FLUX["gDW"]),
    PhysicalType.RELATIVE_ENZYME: ("relative", {"relative": 1.0}),
    PhysicalType.RELATIVE_METABOLITE: ("relative", {"relative": 1.0}),
    PhysicalType.RELATIVE_FLUX: ("relative", {"relative": 1.0}),
}
_ENZYMES = {PhysicalType.ENZYME_CONCENTRATION, PhysicalType.ENZYME_AMOUNT_PER_DRY_MASS, PhysicalType.RELATIVE_ENZYME}
_METABOLITES = {PhysicalType.METABOLITE_CONCENTRATION, PhysicalType.RELATIVE_METABOLITE}
_RATES = {PhysicalType.VOLUMETRIC_FLUX, PhysicalType.SPECIFIC_FLUX, PhysicalType.RELATIVE_FLUX}
_RELATIVE = {PhysicalType.RELATIVE_ENZYME, PhysicalType.RELATIVE_METABOLITE, PhysicalType.RELATIVE_FLUX}
_FAMILIES = {"enzyme_only", "mass_action", "saturation", "activation", "inhibition", "log_linear"}
_ROLES = {"train", "development"}
_FORMULAS = {
    "enzyme_only": "v_oriented = coefficient * E",
    "mass_action": "v_oriented = coefficient * E * S",
    "saturation": "v_oriented = coefficient * E * S / (K_substrate + S)",
    "activation": "v_oriented = coefficient * E * S / (K_substrate + S) * R / (K_regulator + R)",
    "inhibition": "v_oriented = coefficient * E * S / (K_substrate + S) * K_regulator / (K_regulator + R)",
    "log_linear": "v_oriented = coefficient * E * exp(sum(beta_j * (log1p(X_j / scale_j) - center_j)))",
}


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError(f"{name} must be an explicit mapping")
    return value


def _list(value: Any, name: str) -> list:
    if not isinstance(value, (list, tuple)):
        raise ContractError(f"{name} must be an explicit list")
    return list(value)


def _text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"{name} must be a nonempty string")
    return value


def _number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ContractError(f"{name} must be a finite numeric value")
    return float(value)


def content_sha256(value: Any) -> str:
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    except (ValueError, TypeError) as error:
        raise ContractError("digest input must be finite JSON data") from error
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class QuantitySpec:
    physical_type: PhysicalType
    unit: str
    entity_id: str
    reference_id: str | None = None

    @classmethod
    def parse(cls, value: Mapping[str, Any]) -> QuantitySpec:
        value = _mapping(value, "quantity specification")
        if set(value) - {"physical_type", "unit", "entity_id", "reference_id"}:
            raise ContractError("unknown quantity specification fields; arbitrary conversions are not inferred")
        try:
            kind = PhysicalType(value.get("physical_type"))
        except ValueError as error:
            raise ContractError("quantity physical_type is not supported") from error
        unit = _text(value.get("unit"), "quantity unit")
        if unit not in _UNIT_TABLE[kind][1]:
            raise ContractError(f"unit {unit!r} is not explicitly defined for {kind.value}")
        entity_id = _text(value.get("entity_id"), "quantity entity_id")
        reference = value.get("reference_id")
        if kind in _RELATIVE:
            reference = _text(reference, "relative quantity reference_id")
        elif reference is not None:
            raise ContractError("absolute quantity cannot use an unidentified reference scale")
        return cls(kind, unit, entity_id, reference)

    @property
    def factor(self) -> float:
        return _UNIT_TABLE[self.physical_type][1][self.unit]

    @property
    def canonical_unit(self) -> str:
        return _UNIT_TABLE[self.physical_type][0]

    @property
    def relative(self) -> bool:
        return self.physical_type in _RELATIVE

    def to_dict(self) -> dict:
        result = {"physical_type": self.physical_type.value, "unit": self.unit, "entity_id": self.entity_id}
        if self.reference_id is not None:
            result["reference_id"] = self.reference_id
        return result


@dataclass(frozen=True)
class FitSettings:
    max_nfev: int = 200
    profile_points: int = 9
    profile_confidence: float = 0.95
    log_shape_bounds: tuple[float, float] = (-8.0, 8.0)
    coefficient_bounds: tuple[float, float] = (-4.0, 4.0)
    start_values: tuple[float, ...] = (-2.0, 0.0, 2.0)
    condition_limit: float = 1e6
    constant_target_rtol: float = 1e-12
    profile_floor: float = 1e-10
    max_profile_log_span: float = 8.0

    @classmethod
    def parse(cls, value: Mapping[str, Any]) -> FitSettings:
        value = dict(_mapping(value, "settings"))
        if set(value) - {item.name for item in fields(cls)}:
            raise ContractError("unknown optimizer/profile settings")
        for key in ("log_shape_bounds", "coefficient_bounds", "start_values"):
            if key in value:
                value[key] = tuple(_list(value[key], key))
        result = cls(**value)
        for key, lower, upper in (("max_nfev", 20, 1000), ("profile_points", 7, 17)):
            number = getattr(result, key)
            if isinstance(number, bool) or not isinstance(number, int) or not lower <= number <= upper:
                raise ContractError(f"{key} exceeds the bounded computation budget")
        if not 0.8 <= _number(result.profile_confidence, "profile_confidence") <= 0.99:
            raise ContractError("profile_confidence must be between 0.8 and 0.99")
        for name in ("log_shape_bounds", "coefficient_bounds"):
            bounds = getattr(result, name)
            if len(bounds) != 2 or not -12 <= _number(bounds[0], name) < _number(bounds[1], name) <= 12:
                raise ContractError("parameter bounds exceed the bounded search budget")
        if not 1 <= len(result.start_values) <= 3:
            raise ContractError("start_values exceeds the bounded multistart budget")
        for start in result.start_values:
            _number(start, "start_values")
            if any(not bounds[0] < start < bounds[1] for bounds in (result.log_shape_bounds, result.coefficient_bounds)):
                raise ContractError("every optimizer start must be interior to both declared bounds")
        if not 10 <= _number(result.condition_limit, "condition_limit") <= 1e8:
            raise ContractError("condition_limit must be finite and at most 1e8")
        for name in ("constant_target_rtol", "profile_floor"):
            if not 0 < _number(getattr(result, name), name) <= 1e-4:
                raise ContractError(f"{name} must lie in (0, 1e-4]")
        if not 0 < _number(result.max_profile_log_span, "max_profile_log_span") <= 12:
            raise ContractError("max_profile_log_span exceeds the bounded search budget")
        return result


@dataclass(frozen=True, order=True)
class MeasurementIdentity:
    source_id: str
    measurement_id: str

    @classmethod
    def parse(cls, value: Any, default_source: str | None = None) -> MeasurementIdentity:
        if isinstance(value, str) and default_source is not None:
            return cls(default_source, _text(value, "ancestry measurement_id"))
        value = _mapping(value, "measurement ancestry identity")
        if set(value) != {"source_id", "measurement_id"}:
            raise ContractError("measurement ancestry references require explicit source_id and measurement_id")
        return cls(_text(value["source_id"], "ancestry source_id"), _text(value["measurement_id"], "ancestry measurement_id"))

    def token(self) -> tuple[str, str, str]:
        return "measurement", self.source_id, self.measurement_id


_ANCESTRY_IDENTITY_FIELDS = {
    "measurement": ("source_id", "measurement_id"),
    "source_location": ("source_artifact_sha256", "source_locator"),
    "physical_measurement": ("physical_measurement_id",),
    "derivation": ("source_id", "derivation_id"),
}
_PROVENANCE_FIELDS = {
    "source_id", "measurement_id", "evidence_type", "source_artifact_sha256", "source_locator",
    "physical_measurement_id", "derivation_id", "derivation_method", "source_observation_ids", "ancestry",
}


def _serialize_identities(tokens: set | frozenset) -> list[dict]:
    return [{"identity_type": token[0], **dict(zip(_ANCESTRY_IDENTITY_FIELDS[token[0]], token[1:], strict=True))} for token in sorted(tokens)]


def _parse_identities(values: Sequence[Mapping]) -> set[tuple[str, ...]]:
    result = set()
    for value in _list(values, "target measurement ancestry identities"):
        value = _mapping(value, "target measurement ancestry identity")
        kind = value.get("identity_type")
        if not isinstance(kind, str) or kind not in _ANCESTRY_IDENTITY_FIELDS:
            raise ContractError("unknown target measurement ancestry identity type")
        names = _ANCESTRY_IDENTITY_FIELDS[kind]
        if set(value) != {"identity_type", *names}:
            raise ContractError("incomplete target measurement ancestry identity")
        result.add((kind, *(_text(value[name], "measurement ancestry identity") for name in names)))
    return result


def _measurement_ancestry(
    records: Sequence[Mapping], quantity_names: Sequence[str] | None, target: str | None,
    evidence_kind: str, *, forbidden_target_identities: Sequence[Mapping] = (), enforce_groups: bool = True,
) -> dict:
    nodes, roots, pending = {}, [], []
    for record in records:
        provenance = _mapping(record.get("provenance"), "measurement ancestry provenance")
        names = list(provenance) if quantity_names is None else list(quantity_names)
        if set(provenance) != set(names):
            raise ContractError("every projected quantity needs explicit measurement ancestry provenance")
        record_id = _text(record.get("record_id"), "ancestry record_id")
        group_id = _text(record.get("group_id"), "ancestry group_id")
        for name in names:
            origin = _mapping(provenance[name], "measurement ancestry provenance")
            identity = MeasurementIdentity.parse({key: origin.get(key) for key in ("source_id", "measurement_id")})
            roots.append((record_id, group_id, name, identity))
            pending.append((origin, 0))
    if len(pending) > 50000:
        raise ContractError("measurement ancestry exceeds the 50000-node declaration budget")
    declarations, edge_count = 0, 0
    while pending:
        origin, depth = pending.pop()
        declarations += 1
        if declarations > 50000 or depth > 64:
            raise ContractError("measurement ancestry exceeds the bounded declaration/depth budget")
        origin = _mapping(origin, "measurement ancestry node")
        if set(origin) - _PROVENANCE_FIELDS:
            raise ContractError("unknown measurement ancestry provenance fields cannot establish independence")
        identity = MeasurementIdentity.parse({key: origin.get(key) for key in ("source_id", "measurement_id")})
        evidence_type = _text(origin.get("evidence_type"), "ancestry evidence_type")
        if evidence_kind == "synthetic_engineering":
            if evidence_type != "synthetic_engineering":
                raise ContractError("engineering ancestry nodes cannot be labeled as biological measurements")
        elif evidence_kind == "native_measurement":
            if evidence_type not in {"direct_measurement", "measurement_derived"}:
                raise ContractError("native ancestry nodes must identify measured or measurement-derived observations")
        else:
            raise ContractError("measurement ancestry requires an explicit evidence kind")
        aliases = {identity.token()}
        if "physical_measurement_id" in origin:
            aliases.add(("physical_measurement", _text(origin["physical_measurement_id"], "physical measurement identity")))
        located = "source_artifact_sha256" in origin or "source_locator" in origin
        if located or evidence_kind == "native_measurement":
            artifact = _text(origin.get("source_artifact_sha256"), "ancestry source artifact digest")
            locator = _text(origin.get("source_locator"), "ancestry source locator")
            try:
                if len(artifact) != 64 or len(bytes.fromhex(artifact)) != 32:
                    raise ValueError("not SHA256")
            except ValueError as error:
                raise ContractError("ancestry source artifact identity requires a SHA256 digest") from error
            aliases.add(("source_location", artifact.lower(), locator))
        derivation = None
        if "derivation_id" in origin:
            derivation = _text(origin["derivation_id"], "ancestry derivation instance identity")
            aliases.add(("derivation", identity.source_id, derivation))
        if "derivation_method" in origin:
            _text(origin["derivation_method"], "ancestry derivation method")
        parents = None
        if "source_observation_ids" in origin:
            parents = tuple(sorted(MeasurementIdentity.parse(value, identity.source_id) for value in _list(origin["source_observation_ids"], "measurement ancestry source_observation_ids")))
            if len(parents) != len(set(parents)):
                raise ContractError("duplicate parent references in measurement ancestry")
        legacy = set(origin) == {"source_id", "measurement_id", "evidence_type"} and evidence_type == "synthetic_engineering"
        inventory = _list(origin.get("ancestry", []), "measurement ancestry inventory")
        if len(inventory) + declarations + len(pending) > 50000:
            raise ContractError("measurement ancestry exceeds the 50000-node declaration budget")
        pending.extend((node, depth + 1) for node in inventory)
        if identity in nodes:
            prior = nodes[identity]
            if prior["evidence_type"] != evidence_type or prior["derivation_id"] != derivation or (prior["parents"] is not None and parents is not None and prior["parents"] != parents):
                raise ContractError("conflicting definitions of measurement ancestry")
            prior["aliases"].update(aliases)
            if prior["parents"] is None and parents is not None:
                prior["parents"] = parents
                edge_count += len(parents)
            prior["legacy"] = prior["legacy"] and legacy
        else:
            nodes[identity] = {"aliases": aliases, "parents": parents, "evidence_type": evidence_type, "derivation_id": derivation, "legacy": legacy}
            edge_count += len(parents or ())
        if edge_count > 100000:
            raise ContractError("measurement ancestry exceeds the 100000-edge budget")
    representatives = {identity: identity for identity in nodes}

    def representative(identity):
        trail = []
        while identity in representatives and representatives[identity] != identity:
            trail.append(identity)
            identity = representatives[identity]
        for previous in trail:
            representatives[previous] = identity
        return identity

    alias_owners = {}
    for identity, node in nodes.items():
        for token in node["aliases"]:
            if token in alias_owners:
                first, second = representative(identity), representative(alias_owners[token])
                representatives[max(first, second)] = min(first, second)
            else:
                alias_owners[token] = identity
    target_identities = _parse_identities(forbidden_target_identities)
    target_components = {representative(identity) for _, _, name, identity in roots if name == target}
    target_components.update(representative(identity) for identity, node in nodes.items() if node["aliases"] & target_identities)
    for record_id, _, name, identity in roots:
        if name != target and representative(identity) in target_components:
            raise ContractError(f"target/feature measurement ancestry overlap for feature {name!r} in record {record_id!r}; physical aliases identify the same observation")
    joined = {}
    for identity, node in nodes.items():
        canonical = representative(identity)
        parents = None if node["parents"] is None else tuple(sorted({representative(parent) for parent in node["parents"]}))
        if parents and canonical in parents:
            raise ContractError("measurement ancestry contains a cycle through physical aliases")
        if canonical not in joined:
            joined[canonical] = {**node, "parents": parents, "aliases": set(node["aliases"])}
            continue
        prior = joined[canonical]
        if prior["evidence_type"] != node["evidence_type"] or prior["derivation_id"] != node["derivation_id"] or (prior["parents"] is not None and parents is not None and prior["parents"] != parents):
            raise ContractError("conflicting definitions of physical measurement ancestry")
        prior["aliases"].update(node["aliases"])
        if len(prior["aliases"]) > 4096:
            raise ContractError("measurement ancestry exceeds the physical-alias component budget")
        if prior["parents"] is None:
            prior["parents"] = parents
        prior["legacy"] = prior["legacy"] and node["legacy"]
    nodes = joined
    roots = [(record_id, group_id, name, representative(identity)) for record_id, group_id, name, identity in roots]
    edge_count = sum(len(node["parents"] or ()) for node in nodes.values())
    parent_ids = {parent for node in nodes.values() for parent in node["parents"] or ()}
    visiting, resolved, unknown, path_lengths = set(), {}, set(), {}
    closure_entries = 0

    def closure(identity, depth=0):
        nonlocal closure_entries
        if identity in visiting:
            raise ContractError("measurement ancestry contains a cycle")
        if depth > 64:
            raise ContractError("measurement ancestry exceeds the bounded graph-depth budget")
        if identity in resolved:
            return resolved[identity]
        visiting.add(identity)
        node = nodes.get(identity)
        tokens = {identity.token()} if node is None else set(node["aliases"])
        if node is None or node["parents"] is None:
            unknown.add(identity)
        path_length = 0
        if node is not None:
            for parent in node["parents"] or ():
                tokens.update(closure(parent, depth + 1))
                path_length = max(path_length, path_lengths[parent] + 1)
                if len(tokens) > 4096 or path_length > 64:
                    raise ContractError("measurement ancestry exceeds the bounded closure/depth budget")
        path_lengths[identity] = path_length
        visiting.remove(identity)
        closure_entries += len(tokens)
        if closure_entries > 1000000:
            raise ContractError("measurement ancestry exceeds the total closure-entry budget")
        resolved[identity] = frozenset(tokens)
        return resolved[identity]

    for identity in nodes:
        closure(identity)
    for _, _, name, identity in roots:
        if name == target:
            target_identities.update(resolved[identity])
    for record_id, _, name, identity in roots:
        if name != target and resolved[identity] & target_identities:
            raise ContractError(f"target/feature measurement ancestry overlap for feature {name!r} in record {record_id!r}; renaming or regrouping does not create independence")
    if enforce_groups:
        identity_groups = {}
        for _, group_id, _, identity in roots:
            for token in resolved[identity]:
                if token in identity_groups and identity_groups[token] != group_id:
                    raise ContractError("one source or ancestral measurement occurs in multiple groups")
                identity_groups[token] = group_id
    projected_roots = {identity for _, _, _, identity in roots}
    permitted_legacy = {identity for identity in unknown if identity in projected_roots and identity in nodes and nodes[identity]["legacy"] and identity not in parent_ids}
    if unknown - permitted_legacy or (unknown and evidence_kind != "synthetic_engineering"):
        raise ContractError("unresolved measurement ancestry; every derived observation requires a closed explicit source graph")
    for identity in resolved:
        node = nodes.get(identity)
        if node is None:
            continue
        parents, derivation = node["parents"], node["derivation_id"]
        if parents and (derivation is None or node["evidence_type"] == "direct_measurement"):
            raise ContractError("derived measurement ancestry requires a derivation instance and a derived evidence class")
        if derivation is not None and not parents:
            raise ContractError("derived measurement ancestry requires explicit nonempty source_observation_ids")
        if node["evidence_type"] == "measurement_derived" and (not parents or derivation is None):
            raise ContractError("measurement-derived ancestry is incomplete")
    graph = [{"identity": asdict(identity), "aliases": _serialize_identities(node["aliases"]), "parents": None if node["parents"] is None else [asdict(parent) for parent in node["parents"]], "evidence_type": node["evidence_type"], "derivation_id": node["derivation_id"]} for identity, node in sorted(nodes.items())]
    return {
        "schema_version": "native_rate_measurement_ancestry.v1",
        "resolution": "synthetic_identity_only_unresolved" if unknown else "declared_graph_closed",
        "authority": "declared_provenance_not_external_source_authorization",
        "graph_sha256": content_sha256(graph),
        "node_count": len(nodes),
        "edge_count": edge_count,
        "target_identities": _serialize_identities(target_identities),
        "unresolved_synthetic_roots": [asdict(identity) for identity in sorted(unknown)],
    }


class _FitFailure(Exception):
    def __init__(self, code: str, message: str, diagnostics: dict | None = None):
        super().__init__(message)
        self.code = code
        self.diagnostics = diagnostics or {}


@dataclass
class _Prepared:
    measurements: Mapping[str, Any]
    protocol: Mapping[str, Any]
    contract: Mapping[str, Any]
    settings: FitSettings
    quantities: dict[str, QuantitySpec]
    values: dict[str, np.ndarray]
    records: list[Mapping[str, Any]]
    groups: np.ndarray
    target: str
    reaction: Mapping[str, Any]
    evidence_kind: str
    folds: list[dict]
    fit_groups: list[str]
    measurement_ancestry: dict


def _approval_contract(measurements: Mapping, protocol: Mapping, catalogue: Mapping | None) -> Mapping:
    kind = measurements.get("evidence_kind")
    if kind == "synthetic_engineering":
        contract = protocol
        if contract.get("schema_version") != "native_rate_learning_protocol.v1":
            raise ContractError("synthetic fixtures require an explicit native_rate_learning_protocol.v1 contract")
        if contract.get("approval", {}).get("status") != "synthetic_engineering_only":
            raise ContractError("synthetic fixtures require synthetic_engineering_only approval scope")
    elif kind == "native_measurement":
        contract = protocol.get("learning_contract")
        if protocol.get("schema_version") != 1 or not isinstance(contract, Mapping):
            raise ContractError("native training approval requires the independent protocol's explicit learning_contract")
        if contract.get("schema_version") != "native_rate_learning_protocol.v1":
            raise ContractError("native approval has no supported learning contract")
        approval = _mapping(contract.get("approval"), "native training approval")
        if approval.get("status") != "approved_native_training":
            raise ContractError("native training approval is absent")
        if contract.get("protocol_id") != protocol.get("protocol_id"):
            raise ContractError("native training approval belongs to another protocol")
        _text(approval.get("approved_by"), "independent training approval issuer")
        _text(approval.get("external_registration_receipt_sha256"), "external registration receipt digest")
        catalogue = _mapping(catalogue, "approved biochemical catalogue")
        allowlist = _list(catalogue.get("access_policy", {}).get("learner_measurement_allowlist"), "catalogue learner allowlist")
        if not allowlist or approval.get("catalogue_allowlist_entry") not in allowlist:
            raise ContractError("native training approval is not in the custodian learner allowlist")
        if approval.get("catalogue_sha256") != content_sha256(catalogue):
            raise ContractError("native training approval catalogue digest mismatch")
        if approval.get("measurements_sha256") != content_sha256(measurements):
            raise ContractError("native training approval measurement digest mismatch")
    else:
        raise ContractError("evidence_kind must explicitly distinguish native_measurement from synthetic_engineering")
    if contract.get("evidence_kind") != kind:
        raise ContractError("protocol and measurement evidence kinds differ")
    _text(contract.get("protocol_id"), "protocol_id")
    return contract


def _prepare(measurements: Mapping, protocol: Mapping, catalogue: Mapping | None) -> _Prepared:
    measurements = _mapping(measurements, "measurements")
    protocol = _mapping(protocol, "protocol")
    if measurements.get("schema_version") != "native_rate_measurements.v1":
        raise ContractError("measurements require an explicit native_rate_measurements.v1 projection; no columns are inferred")
    if set(measurements) != {"schema_version", "evidence_kind", "reaction", "target", "quantities", "records"}:
        raise ContractError("measurement projection fields must follow the explicit typed schema")
    records = _list(measurements.get("records"), "measurement records")
    if not 1 <= len(records) <= 20000:
        raise ContractError("measurement record count exceeds the bounded computation budget")
    for record in records:
        _mapping(record, "measurement record")
        if not isinstance(record.get("role"), str) or record.get("role") not in _ROLES:
            raise ContractError("only explicitly permitted train/development roles may enter fitting")
        if set(record) - {"record_id", "group_id", "role", "grouping", "values", "provenance", "standard_uncertainties"}:
            raise ContractError("measurement record projection fields must follow the explicit typed schema")
    contract = _approval_contract(measurements, protocol, catalogue)
    permitted = _mapping(contract.get("permitted_groups"), "permitted_groups")
    if not 3 <= len(permitted) <= 64 or any(role not in _ROLES for role in permitted.values()):
        raise ContractError("permitted_groups must contain 3 to 64 train/development groups")
    if "train" not in permitted.values():
        raise ContractError("training groups are required")
    factors = _list(contract.get("independence_factors"), "independence_factors")
    if not factors or len(set(factors)) != len(factors):
        raise ContractError("explicit unique independence factors are required")
    factors = [_text(item, "independence factor") for item in factors]
    record_ids, unit_groups = set(), {}
    for record in records:
        record_id = _text(record.get("record_id"), "record_id")
        if record_id in record_ids:
            raise ContractError("duplicate measurement record_id")
        record_ids.add(record_id)
        group = _text(record.get("group_id"), "group_id")
        if group not in permitted:
            raise ContractError(f"group {group!r} is not permitted")
        if record["role"] != permitted[group]:
            raise ContractError("record role disagrees with the permitted group role")
        grouping = _mapping(record.get("grouping"), "record independence grouping")
        for factor in factors:
            identity = (factor, _text(grouping.get(factor), "independence identity"))
            if identity in unit_groups and unit_groups[identity] != group:
                raise ContractError("one independence component appears in multiple groups")
            unit_groups[identity] = group
    groups = np.asarray([record["group_id"] for record in records], dtype=object)
    if set(groups) != set(permitted):
        raise ContractError("projection does not contain exactly the permitted groups")
    selection = _mapping(contract.get("selection"), "selection")
    mode = selection.get("mode")
    if mode not in {"group_cross_validation", "development_groups"}:
        raise ContractError("selection mode must be explicitly grouped")
    if selection.get("metric") != "equal_group_scaled_mae":
        raise ContractError("the supported selection metric is equal_group_scaled_mae")
    if not 0 <= _number(selection.get("tie_tolerance"), "tie_tolerance") <= 1e-8:
        raise ContractError("tie_tolerance must be predeclared and at most 1e-8")
    if selection.get("refit_policy") not in {"train_only", "all_permitted_groups"}:
        raise ContractError("refit_policy must explicitly identify the permitted fitting roles")
    if mode == "development_groups" and selection["refit_policy"] != "train_only":
        raise ContractError("development selection requires train_only normalization and refitting")
    if mode == "development_groups" and "development" not in permitted.values():
        raise ContractError("development selection requires development groups")
    if measurements["evidence_kind"] == "native_measurement":
        if mode != "development_groups" or selection["tie_tolerance"] != 1e-12:
            raise ContractError("native protocol requires development-group MAE and its preregistered tie tolerance")
        if list(permitted.values()).count("train") < 3 or list(permitted.values()).count("development") < 2:
            raise ContractError("native protocol requires at least three training and two development groups")
    folds = _list(contract.get("folds"), "folds")
    if not 2 <= len(folds) <= 64:
        raise ContractError("fold count exceeds the bounded grouped-CV budget")
    seen_validation, fold_ids = Counter(), set()
    training_groups = {group for group, role in permitted.items() if role == "train"}
    development_groups = set(permitted) - training_groups
    development_folds = 0
    for fold in folds:
        _mapping(fold, "fold")
        fold_id = _text(fold.get("fold_id"), "fold_id")
        if fold_id in fold_ids:
            raise ContractError("fold_id must be unique")
        fold_ids.add(fold_id)
        fit_list = _list(fold.get("fit_groups"), "fold fit_groups")
        validation_list = _list(fold.get("validation_groups"), "fold validation_groups")
        fit, validation = set(fit_list), set(validation_list)
        if len(fit) != len(fit_list) or len(validation) != len(validation_list) or not fit.isdisjoint(validation):
            raise ContractError("fold fit and validation groups must be unique and disjoint")
        if len(fit) < 2 or not validation or not fit.union(validation) <= set(permitted):
            raise ContractError("fold contains missing or unpermitted groups")
        if mode == "group_cross_validation" and fit.union(validation) != set(permitted):
            raise ContractError("each CV fold must partition all permitted whole groups")
        if mode == "development_groups":
            if not fit <= training_groups:
                raise ContractError("development groups cannot enter fold fitting or normalization")
            if validation & development_groups:
                if validation != development_groups or fit != training_groups:
                    raise ContractError("development fold must hold out all whole development groups")
                development_folds += 1
            elif fit.union(validation) != training_groups:
                raise ContractError("training CV folds must partition training groups only")
        seen_validation.update(validation)
    if set(seen_validation) != set(permitted) or any(count != 1 for count in seen_validation.values()):
        raise ContractError("each permitted whole group must be validated exactly once")
    if mode == "development_groups" and development_folds != 1:
        raise ContractError("exactly one complete development fold is required")
    target = _text(measurements.get("target"), "target")
    quantity_values = _mapping(measurements.get("quantities"), "quantities")
    if not 2 <= len(quantity_values) <= 64:
        raise ContractError("quantity count exceeds the bounded feature budget")
    quantities = {_text(name, "quantity name"): QuantitySpec.parse(spec) for name, spec in quantity_values.items()}
    if target not in quantities or quantities[target].physical_type not in _RATES:
        raise ContractError("target must be an explicitly typed reaction-rate quantity")
    if any(spec.physical_type in _RATES for name, spec in quantities.items() if name != target):
        raise ContractError("additional outcomes cannot enter the feature projection")
    reaction = _mapping(measurements.get("reaction"), "reaction")
    _text(reaction.get("reaction_id"), "reaction_id")
    _text(reaction.get("source_id"), "reaction stoichiometry source_id")
    if quantities[target].entity_id != reaction["reaction_id"]:
        raise ContractError("target entity does not identify the declared reaction")
    direction = reaction.get("direction")
    if isinstance(direction, bool) or direction not in (-1, 1):
        raise ContractError("reaction direction must be explicitly +1 or -1")
    for species, coefficient in _mapping(reaction.get("stoichiometry"), "provided stoichiometry").items():
        _text(species, "stoichiometric species")
        if _number(coefficient, "stoichiometric coefficient") == 0:
            raise ContractError("provided stoichiometric coefficients cannot be zero")
    ancestry = _measurement_ancestry(records, list(quantities), target, measurements["evidence_kind"])
    numeric = {name: [] for name in quantities}
    native = measurements["evidence_kind"] == "native_measurement"
    permitted_evidence = contract.get("quantity_evidence_types", {})
    for record in records:
        values = _mapping(record.get("values"), "record values")
        provenance = _mapping(record.get("provenance"), "record provenance")
        uncertainties = _mapping(record.get("standard_uncertainties", {}), "source uncertainty mapping")
        if not set(uncertainties) <= set(quantities) or (native and set(uncertainties) != set(quantities)):
            raise ContractError("source uncertainty entries must identify projected quantities; native measurements need every uncertainty")
        if set(values) != set(quantities) or set(provenance) != set(quantities):
            raise ContractError("every explicitly projected quantity needs a value and measurement provenance")
        for name, spec in quantities.items():
            origin = _mapping(provenance[name], "quantity provenance")
            evidence_type = origin["evidence_type"]
            if native and evidence_type not in permitted_evidence.get(name, []):
                raise ContractError("native quantity provenance is not an approved measurement class")
            number = _number(values[name], f"{name} measurement") * spec.factor
            if name == target:
                number *= direction
                if number < 0:
                    raise ContractError("measured rate violates the declared one-direction law; signed/reversible laws are not inferred")
            elif number < 0:
                raise ContractError("enzyme and metabolite quantities must be nonnegative")
            if not math.isfinite(number):
                raise ContractError("canonical unit conversion produced a nonfinite quantity")
            numeric[name].append(number)
            if name in uncertainties:
                uncertainty = _number(uncertainties[name], f"{name} source standard uncertainty") * spec.factor
                if uncertainty < 0 or not math.isfinite(uncertainty):
                    raise ContractError("source uncertainty must convert to a finite nonnegative canonical quantity")
    settings = FitSettings.parse(contract.get("settings", {}))
    candidates = _list(contract.get("candidates"), "candidates")
    if not 1 <= len(candidates) <= 24:
        raise ContractError("candidate count exceeds the bounded computation budget")
    candidate_ids = [_text(_mapping(candidate, "candidate").get("candidate_id"), "candidate_id") for candidate in candidates]
    if len(set(candidate_ids)) != len(candidate_ids):
        raise ContractError("candidate_id must be unique")
    candidate_families = {item["family"] for item in candidates if isinstance(item.get("family"), str)}
    if native and ("enzyme_only" not in candidate_families or len(candidate_families & _FAMILIES) < 2):
        raise ContractError("native selection needs distinct functional forms and the enzyme-only dependency ablation")
    call_bound = 0
    for candidate in candidates:
        count = min(len(_shape_names(candidate)), 2)
        call_bound += (len(folds) + 1) * (len(settings.start_values) ** count + count * (settings.profile_points + 1))
    if call_bound > 5000:
        raise ContractError("declared candidates/folds exceed the 5000-solver-call computation budget")
    content_sha256(measurements)
    content_sha256(protocol)
    fit_groups = sorted(training_groups if selection["refit_policy"] == "train_only" else permitted)
    return _Prepared(measurements, protocol, contract, settings, quantities, {name: np.asarray(value) for name, value in numeric.items()}, records, groups, target, reaction, measurements["evidence_kind"], folds, fit_groups, ancestry)


def _shape_names(candidate: Mapping) -> list[str]:
    family = candidate.get("family")
    if not isinstance(family, str):
        return []
    if family in {"saturation", "activation", "inhibition"}:
        return ["K_substrate"] + (["K_regulator"] if family != "saturation" else [])
    if family == "log_linear":
        covariates = candidate.get("covariates", [])
        if not isinstance(covariates, (list, tuple)):
            return []
        return [f"beta:{name}" for name in covariates]
    return []


def _features(candidate: Mapping) -> list[str]:
    result = [candidate.get("enzyme")]
    if candidate.get("family") in {"mass_action", "saturation", "activation", "inhibition"}:
        result.append(candidate.get("substrate"))
    if candidate.get("family") in {"activation", "inhibition"}:
        result.append(candidate.get("regulator"))
    if candidate.get("family") == "log_linear":
        result.extend(candidate.get("covariates", []))
    return result


def _check_candidate(candidate: Mapping, data: _Prepared) -> None:
    family = candidate.get("family")
    if not isinstance(family, str) or family not in _FAMILIES:
        raise _FitFailure("unknown_family", "candidate family is not implemented")
    allowed = {"candidate_id", "family", "enzyme"}
    if family in {"mass_action", "saturation", "activation", "inhibition"}:
        allowed.add("substrate")
    if family in {"activation", "inhibition"}:
        allowed.add("regulator")
    if family == "log_linear":
        allowed.add("covariates")
        if not isinstance(candidate.get("covariates"), (list, tuple)) or not 1 <= len(candidate["covariates"]) <= 2:
            raise _FitFailure("feature_budget", "log_linear requires one or two predeclared covariates")
    if set(candidate) != allowed:
        raise _FitFailure("candidate_schema", "candidate has missing or unsupported settings")
    names = _features(candidate)
    if data.target in names:
        raise _FitFailure("target_as_feature", "the outcome cannot be used as an input feature")
    if any(not isinstance(name, str) or name not in data.quantities for name in names):
        raise _FitFailure("missing_feature", "candidate refers to an unprojected feature")
    if len(set(names)) != len(names):
        raise _FitFailure("duplicate_feature", "enzyme, substrate, regulator and covariates must be distinct")
    if data.quantities[names[0]].physical_type not in _ENZYMES:
        raise _FitFailure("enzyme_units", "enzyme must be a measured enzyme quantity, not gene copy or a free gain")
    if any(data.quantities[name].physical_type not in _METABOLITES for name in names[1:]):
        raise _FitFailure("metabolite_units", "substrates and candidate regulators must be measured metabolite quantities")
    if "substrate" in candidate and data.reaction["stoichiometry"]:
        entity = data.quantities[candidate["substrate"]].entity_id
        coefficient = data.reaction["stoichiometry"].get(entity)
        if coefficient is None or coefficient * data.reaction["direction"] >= 0:
            raise _FitFailure("substrate_stoichiometry", "selected substrate is not consumed in the declared stoichiometric direction")


def _weights(groups: np.ndarray) -> np.ndarray:
    counts = Counter(groups)
    return np.asarray([1.0 / (len(counts) * counts[group]) for group in groups])


def _constant(values: np.ndarray, tolerance: float) -> bool:
    maximum = float(np.max(np.abs(values)))
    return maximum == 0 or float(np.ptp(values)) <= tolerance * maximum


def _iqr(values: np.ndarray, groups: np.ndarray, tolerance: float) -> float:
    means = np.asarray([np.mean(values[groups == group]) for group in sorted(set(groups))])
    lower, upper = np.quantile(means, [0.25, 0.75], method="linear")
    result = float(upper - lower)
    if not math.isfinite(result) or result <= tolerance * float(np.max(np.abs(means))):
        raise _FitFailure("zero_training_iqr", "training-group mean IQR is not strictly positive; no replacement scoring scale is fitted")
    return result


def _indices(data: _Prepared, groups: Sequence[str]) -> np.ndarray:
    return np.flatnonzero(np.isin(data.groups, groups))


def _normalization(data: _Prepared, candidate: Mapping, indices: np.ndarray) -> dict:
    weights = _weights(data.groups[indices])
    scales = {name: float(weights @ data.values[name][indices]) for name in _features(candidate)}
    if any(value <= 0 or not math.isfinite(value) for value in scales.values()):
        raise _FitFailure("unidentifiable", "a required measured feature has no positive training support")
    rate_scale = float(weights @ data.values[data.target][indices])
    if rate_scale <= 0 or not math.isfinite(rate_scale):
        raise _FitFailure("unidentifiable", "no positive measured training rate is available")
    centers = {
        name: float(weights @ np.log1p(data.values[name][indices] / scales[name]))
        for name in candidate.get("covariates", [])
    }
    return {
        "features": scales,
        "rate": rate_scale,
        "log_feature_centers": centers,
        "training_group_mean_iqr": _iqr(data.values[data.target][indices], data.groups[indices], data.settings.constant_target_rtol),
        "fit_groups": sorted(set(data.groups[indices])),
        "scope": "fit_groups_only_canonical_units",
    }


def _base(candidate: Mapping, x: Mapping[str, np.ndarray], theta: np.ndarray, centers: Mapping) -> np.ndarray:
    family = candidate["family"]
    result = x[candidate["enzyme"]].copy()
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        if family == "mass_action":
            result *= x[candidate["substrate"]]
        elif family in {"saturation", "activation", "inhibition"}:
            substrate = x[candidate["substrate"]]
            result *= substrate / (math.exp(float(theta[0])) + substrate)
            if family in {"activation", "inhibition"}:
                regulator, half = x[candidate["regulator"]], math.exp(float(theta[1]))
                result *= regulator / (half + regulator) if family == "activation" else half / (half + regulator)
        elif family == "log_linear":
            exponent = np.zeros_like(result)
            for beta, name in zip(theta, candidate["covariates"], strict=True):
                exponent += beta * (np.log1p(x[name]) - centers[name])
            result *= np.exp(exponent)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite candidate response")
    return result


def _profile(
    candidate: Mapping, theta: np.ndarray, residual, bounds: tuple, objective: float,
    degrees: int, normalization: Mapping, settings: FitSettings,
) -> dict:
    profiles = {}
    threshold = objective + max(settings.profile_floor, objective * float(f_distribution.ppf(settings.profile_confidence, 1, degrees)) / degrees)
    names = _shape_names(candidate)
    for index, name in enumerate(names):
        grid = sorted(set(np.linspace(bounds[0], bounds[1], settings.profile_points).tolist() + [float(theta[index])]))
        points = []
        others = [j for j in range(len(theta)) if j != index]
        for fixed in grid:
            trial = theta.copy()
            trial[index] = fixed
            status, nfev, loss = "fixed_shape", 0, None
            try:
                if others:
                    def reduced(parameters):
                        varying = trial.copy()
                        varying[others] = parameters
                        return residual(varying)
                    optimized = least_squares(reduced, trial[others], bounds=bounds, max_nfev=settings.max_nfev, ftol=1e-12, xtol=1e-12, gtol=1e-12)
                    nfev = int(optimized.nfev)
                    if optimized.success:
                        loss = float(optimized.fun @ optimized.fun)
                        status = "converged"
                    else:
                        status = f"failed:{optimized.status}"
                else:
                    errors = residual(trial)
                    loss = float(errors @ errors)
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
                status = f"failed:{type(error).__name__}"
            feature = candidate.get("substrate") if name == "K_substrate" else candidate.get("regulator")
            value = fixed if candidate["family"] == "log_linear" else math.exp(fixed) * normalization["features"][feature]
            points.append({"value": value, "loss": loss, "optimizer_status": status, "nfev": nfev})
        accepted = [i for i, point in enumerate(points) if point["loss"] is not None and point["loss"] <= threshold]
        if any(point["loss"] is None for point in points) or not accepted:
            status, interval = "incomplete", None
        elif 0 in accepted or len(points) - 1 in accepted:
            status, interval = "reaches_optimization_bound", None
        else:
            interval = [points[min(accepted) - 1]["value"], points[max(accepted) + 1]["value"]]
            width = interval[1] - interval[0] if candidate["family"] == "log_linear" else math.log(interval[1] / interval[0])
            status = "bounded_on_grid" if width <= settings.max_profile_log_span else "weakly_identified"
        profiles[name] = {
            "status": status,
            "interval": interval,
            "threshold": threshold,
            "points": points,
            "interpretation": "outer_profile_grid_bracket_not_calibrated_coverage",
            "group_degrees_of_freedom": degrees,
        }
    return profiles


def _parameter_templates(candidate: Mapping) -> dict:
    return {name: {"value": None, "status": "no_identified_point_estimate"} for name in ["coefficient", *_shape_names(candidate)]}


def _fit_one(data: _Prepared, candidate: Mapping, indices: np.ndarray) -> dict:
    settings = data.settings
    target = data.values[data.target][indices]
    if _constant(target, settings.constant_target_rtol):
        raise _FitFailure("constant_training_target", "training target is constant; no kinetic identification is asserted")
    normalization = _normalization(data, candidate, indices)
    weights = _weights(data.groups[indices])
    root_weights = np.sqrt(weights)
    y = target / normalization["rate"]
    x = {name: data.values[name][indices] / scale for name, scale in normalization["features"].items()}
    names = _shape_names(candidate)
    bounds = settings.coefficient_bounds if candidate["family"] == "log_linear" else settings.log_shape_bounds
    centers = normalization["log_feature_centers"]

    def projected(theta):
        base = _base(candidate, x, np.asarray(theta), centers)
        denominator = float(weights @ (base * base))
        if denominator <= 0:
            raise FloatingPointError("no enzyme/substrate support for a positive rate coefficient")
        amplitude = max(0.0, float(weights @ (base * y)) / denominator)
        return amplitude, base

    def residual(theta):
        amplitude, base = projected(theta)
        return root_weights * (amplitude * base - y)

    starts = list(itertools.product(settings.start_values, repeat=len(names))) if names else [()]
    attempts, solutions = [], []
    for start in starts:
        try:
            if names:
                fit = least_squares(residual, np.asarray(start), bounds=bounds, max_nfev=settings.max_nfev, ftol=1e-12, xtol=1e-12, gtol=1e-12)
                theta, success, evaluations, status = fit.x, bool(fit.success), int(fit.nfev), int(fit.status)
                errors = fit.fun
                message = str(fit.message)
            else:
                theta = np.asarray([], dtype=float)
                errors = residual(theta)
                success, evaluations, status, message = True, 0, 1, "closed_form_nonnegative_coefficient"
            objective = float(errors @ errors)
            success = success and math.isfinite(objective)
            attempts.append({"initial_shape": list(start), "success": success, "nfev": evaluations, "status": status, "message": message, "objective": objective if math.isfinite(objective) else None})
            if success:
                solutions.append((objective, theta.copy(), len(attempts) - 1))
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
            attempts.append({"initial_shape": list(start), "success": False, "nfev": None, "status": "numerical_failure", "message": str(error), "objective": None})
    diagnostics = {
        "normalization": normalization,
        "optimization": {
            "method": "bounded_variable_projection_least_squares",
            "objective": "equal_group_mean_squared_canonical_rate_residual_scaled_for_optimization",
            "coefficient_solution": "nonnegative_closed_form_conditional_on_shape",
            "shape_coordinate": "linear_beta" if candidate["family"] == "log_linear" else "natural_log_half_saturation_over_fit_feature_scale",
            "bounds": list(bounds),
            "max_nfev_per_start": settings.max_nfev,
            "ftol": 1e-12, "xtol": 1e-12, "gtol": 1e-12,
            "starts": attempts,
        },
    }
    if not solutions:
        raise _FitFailure("optimization_failed", "every bounded optimization attempt failed", diagnostics)
    objective, theta, selected_start = min(solutions, key=lambda item: (item[0], item[2]))
    diagnostics["optimization"]["selected_start"] = selected_start
    amplitude, base = projected(theta)
    if amplitude <= 0:
        raise _FitFailure("unidentifiable", "positive rate coefficient is not supported", diagnostics)
    comparable = [item for item in solutions if item[0] <= objective + max(settings.profile_floor * 0.01, objective * 1e-6)]
    if names and any(float(np.max(np.abs(item[1] - theta))) > 0.1 for item in comparable):
        diagnostics["equivalent_separated_optima"] = len(comparable)
        raise _FitFailure("unidentifiable", "separated parameter solutions explain the data equally well; no arbitrary point is retained", diagnostics)
    parameter_count = 1 + len(names)
    group_count = len(set(data.groups[indices]))
    degrees = group_count - parameter_count
    if degrees <= 0:
        raise _FitFailure("unidentifiable", "too few independent training groups for parameter uncertainty", diagnostics)
    full_theta = np.r_[math.log(amplitude), theta]
    jacobian = np.empty((len(indices), parameter_count))
    step = 1e-5
    for j in range(parameter_count):
        plus, minus = full_theta.copy(), full_theta.copy()
        plus[j] += step
        minus[j] -= step
        upper = math.exp(float(plus[0])) * _base(candidate, x, plus[1:], centers)
        lower = math.exp(float(minus[0])) * _base(candidate, x, minus[1:], centers)
        jacobian[:, j] = root_weights * (upper - lower) / (2 * step)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else math.inf
    diagnostics["sensitivity"] = {
        "singular_values": singular.tolist(),
        "condition_number": condition if math.isfinite(condition) else None,
        "condition_limit": settings.condition_limit,
        "independent_group_count": group_count,
        "parameter_count": parameter_count,
        "degrees_of_freedom": degrees,
    }
    if not math.isfinite(condition) or condition > settings.condition_limit:
        raise _FitFailure("unidentifiable", "measured data do not separate the candidate parameters", diagnostics)
    if names and np.any(np.minimum(theta - bounds[0], bounds[1] - theta) < 1e-5):
        raise _FitFailure("unidentifiable", "optimization-bound point estimates are not retained as identified constants", diagnostics)
    profiles = _profile(candidate, theta, residual, bounds, objective, degrees, normalization, settings)
    diagnostics["profiles"] = profiles
    if any(profile["status"] != "bounded_on_grid" for profile in profiles.values()):
        raise _FitFailure("unidentifiable", "bounded residual profiles do not identify every shape parameter", diagnostics)
    inverse_information = np.linalg.inv(jacobian.T @ jacobian)
    errors = root_weights * (amplitude * base - y)
    meat = np.zeros((parameter_count, parameter_count))
    for group in set(data.groups[indices]):
        mask = data.groups[indices] == group
        score = jacobian[mask].T @ errors[mask]
        meat += np.outer(score, score)
    covariance = inverse_information @ meat @ inverse_information * group_count / degrees
    covariance = (covariance + covariance.T) / 2
    if not np.all(np.isfinite(covariance)):
        raise _FitFailure("uncertainty_failed", "cluster covariance is nonfinite", diagnostics)
    quantile = float(t_distribution.ppf((1 + settings.profile_confidence) / 2, degrees))
    coefficient = amplitude * normalization["rate"] / normalization["features"][candidate["enzyme"]]
    denominators = [candidate["enzyme"]]
    if candidate["family"] == "mass_action":
        coefficient /= normalization["features"][candidate["substrate"]]
        denominators.append(candidate["substrate"])
    coefficient_unit = f"({data.quantities[data.target].canonical_unit})" + "".join(f"/({data.quantities[name].canonical_unit})" for name in denominators)
    relative_coefficient = any(data.quantities[name].relative for name in [data.target, *denominators])
    coefficient_error = quantile * math.sqrt(max(0.0, float(covariance[0, 0])))
    if coefficient_error > 600:
        raise _FitFailure("unidentifiable", "coefficient uncertainty is not finite on a useful numerical scale", diagnostics)
    parameters = {
        "coefficient": {
            "value": coefficient,
            "unit": coefficient_unit,
            "interpretation": "relative_effective_coefficient" if relative_coefficient else "coefficient_per_measured_enzyme_not_intrinsic_kcat",
            "interval": [coefficient * math.exp(-coefficient_error), coefficient * math.exp(coefficient_error)],
            "uncertainty_method": "conditional_group_cluster_sandwich_log_parameter",
            "reference_ids": {name: data.quantities[name].reference_id for name in [data.target, *denominators] if data.quantities[name].relative},
        }
    }
    for j, name in enumerate(names):
        if candidate["family"] == "log_linear":
            value, unit, interpretation, reference = float(theta[j]), "dimensionless", "predictive_coefficient_not_kinetic_constant", None
        else:
            feature = candidate["substrate"] if name == "K_substrate" else candidate["regulator"]
            spec = data.quantities[feature]
            value = math.exp(float(theta[j])) * normalization["features"][feature]
            unit, reference = spec.canonical_unit, spec.reference_id
            interpretation = "relative_effective_half_saturation" if spec.relative else "conditional_measured_concentration_half_saturation"
        parameters[name] = {"value": value, "unit": unit, "interpretation": interpretation, "reference_id": reference, "interval": profiles[name]["interval"], "profile": profiles[name]}
    source_records = [data.records[int(i)] for i in indices]
    model = {
        "schema_version": "native_rate_model.v1",
        "candidate_id": candidate["candidate_id"],
        "specification": copy.deepcopy(dict(candidate)),
        "formula": _FORMULAS[candidate["family"]],
        "target": data.target,
        "reaction": copy.deepcopy(dict(data.reaction)),
        "quantities": {name: data.quantities[name].to_dict() for name in [data.target, *_features(candidate)]},
        "parameters": parameters,
        "normalization": normalization,
        "feature_choices": _features(candidate),
        "training_domain": {name: {"minimum": float(np.min(data.values[name][indices])), "maximum": float(np.max(data.values[name][indices])), "unit": data.quantities[name].canonical_unit} for name in _features(candidate)},
        "uncertainty": {
            "method": "group_cluster_sandwich_and_bounded_residual_profiles",
            "parameter_coordinates": ["log_coefficient", *[(name if name.startswith("beta:") else f"log_{name}") for name in names]],
            "covariance": covariance.tolist(),
            "quantile": quantile,
            "independent_group_degrees_of_freedom": degrees,
            "source_measurement_uncertainty": "reported_separately_not_propagated_by_this_estimator",
            "coverage_claim": "none",
        },
        "identifiability": "conditional_on_candidate_and_declared_search_bounds",
        "given_invariants": {"stoichiometry": "provided_not_learned", "direction": "provided_not_learned", "nonnegative_oriented_rate": "imposed_not_learned", "enzyme_proportionality": "candidate_assumption_not_learned"},
        "training_record_ids": [record["record_id"] for record in source_records],
        "training_records_sha256": content_sha256(source_records),
        "measurement_ancestry": copy.deepcopy(data.measurement_ancestry),
        "protocol_sha256": content_sha256(data.protocol),
        "evidence_kind": data.evidence_kind,
    }
    model["model_sha256"] = content_sha256(model)
    return {"model": model, "normalization": normalization, "parameters": parameters, "diagnostics": diagnostics}


def _feature_records(data: _Prepared, indices: np.ndarray) -> list[dict]:
    return [
        {"record_id": data.records[int(i)]["record_id"], "group_id": data.records[int(i)]["group_id"], "values": {name: value for name, value in data.records[int(i)]["values"].items() if name != data.target}, "provenance": {name: value for name, value in data.records[int(i)]["provenance"].items() if name != data.target}}
        for i in indices
    ]


def predict_rate_law(model: Mapping[str, Any], records: Sequence[Mapping[str, Any]], quantities: Mapping[str, Any]) -> list[dict]:
    model = _mapping(model, "frozen model")
    if model.get("schema_version") != "native_rate_model.v1":
        raise ContractError("prediction requires a complete native_rate_model.v1")
    digest = model.get("model_sha256")
    if digest != content_sha256({key: value for key, value in model.items() if key != "model_sha256"}):
        raise ContractError("frozen model digest mismatch")
    frozen_ancestry = _mapping(model.get("measurement_ancestry"), "frozen model measurement ancestry")
    if frozen_ancestry.get("schema_version") != "native_rate_measurement_ancestry.v1" or not frozen_ancestry.get("target_identities"):
        raise ContractError("frozen model lacks a complete target measurement ancestry registry")
    if model["evidence_kind"] == "native_measurement" and frozen_ancestry.get("resolution") != "declared_graph_closed":
        raise ContractError("native prediction cannot upgrade unresolved measurement ancestry")
    records = _list(records, "feature-only prediction records")
    if len(records) > 20000:
        raise ContractError("prediction record count exceeds the bounded computation budget")
    given = {name: QuantitySpec.parse(spec) for name, spec in _mapping(quantities, "prediction quantities").items()}
    if model["target"] in given or any(spec.physical_type in _RATES for spec in given.values()):
        raise ContractError("prediction inputs must be feature-only; outcomes are forbidden")
    candidate = model["specification"]
    required = model["feature_choices"]
    for name in required:
        if name not in given:
            raise ContractError(f"missing prediction feature {name!r}")
        expected = QuantitySpec.parse(model["quantities"][name])
        if given[name].physical_type != expected.physical_type or given[name].entity_id != expected.entity_id:
            raise ContractError("prediction feature physical type/entity differs from the fitted measurement")
        if given[name].reference_id != expected.reference_id:
            raise ContractError("prediction feature reference scale differs from the fitted measurement")
    target_spec = QuantitySpec.parse(model["quantities"][model["target"]])
    seen = set()
    for record in records:
        record = _mapping(record, "prediction record")
        record_id = _text(record.get("record_id"), "prediction record_id")
        _text(record.get("group_id"), "prediction group_id")
        if record_id in seen:
            raise ContractError("duplicate prediction record_id")
        seen.add(record_id)
        values = _mapping(record.get("values"), "feature-only prediction values")
        if model["target"] in values:
            raise ContractError("prediction values must be feature-only")
        if not set(values) <= set(given) or not set(required) <= set(values):
            raise ContractError("prediction contains unknown or missing features")
        origin = _mapping(record.get("provenance"), "prediction feature provenance")
        if set(origin) != set(values):
            raise ContractError("prediction measurement ancestry must cover every provided feature")
    ancestry = _measurement_ancestry(records, None, None, model["evidence_kind"], forbidden_target_identities=frozen_ancestry["target_identities"], enforce_groups=False)
    predictions = []
    covariance = np.asarray(model["uncertainty"]["covariance"])
    for record in records:
        record_id, group_id = record["record_id"], record["group_id"]
        values, origin = record["values"], record["provenance"]
        canonical = {}
        for name in required:
            _mapping(origin.get(name), "prediction feature provenance")
            number = _number(values[name], f"prediction feature {name}") * given[name].factor
            if number < 0 or not math.isfinite(number):
                raise ContractError("prediction features must convert to finite nonnegative quantities")
            canonical[name] = number
        parameters = model["parameters"]
        family = candidate["family"]
        try:
            rate = parameters["coefficient"]["value"] * canonical[candidate["enzyme"]]
            gradient = [1.0]
            if family == "mass_action":
                rate *= canonical[candidate["substrate"]]
            elif family in {"saturation", "activation", "inhibition"}:
                substrate, half = canonical[candidate["substrate"]], parameters["K_substrate"]["value"]
                rate *= substrate / (half + substrate)
                gradient.append(-half / (half + substrate))
                if family in {"activation", "inhibition"}:
                    regulator, half = canonical[candidate["regulator"]], parameters["K_regulator"]["value"]
                    rate *= regulator / (half + regulator) if family == "activation" else half / (half + regulator)
                    gradient.append(-half / (half + regulator) if family == "activation" else regulator / (half + regulator))
            elif family == "log_linear":
                exponent = 0.0
                for name in candidate["covariates"]:
                    transformed = math.log1p(canonical[name] / model["normalization"]["features"][name]) - model["normalization"]["log_feature_centers"][name]
                    exponent += parameters[f"beta:{name}"]["value"] * transformed
                    gradient.append(transformed)
                rate *= math.exp(exponent)
            standard = math.sqrt(max(0.0, float(np.asarray(gradient) @ covariance @ np.asarray(gradient))))
            half_width = model["uncertainty"]["quantile"] * standard
            interval = [rate * math.exp(-half_width), rate * math.exp(half_width)] if rate else [0.0, 0.0]
            direction = model["reaction"]["direction"]
            signed = direction * rate
            interval = sorted(direction * bound / target_spec.factor for bound in interval)
            reported = signed / target_spec.factor
            if not all(math.isfinite(value) for value in [signed, reported, *interval]):
                raise FloatingPointError("nonfinite extrapolated prediction")
        except (ValueError, OverflowError, FloatingPointError) as error:
            raise ContractError(f"prediction cannot be represented as finite physical quantities: {error}") from error
        predictions.append(
            {
                "record_id": record_id,
                "group_id": group_id,
                "predicted_rate": reported,
                "unit": target_spec.unit,
                "predicted_rate_canonical": signed,
                "canonical_unit": target_spec.canonical_unit,
                "model_sha256": digest,
                "measurement_ancestry_resolution": ancestry["resolution"],
                "measurement_ancestry_sha256": ancestry["graph_sha256"],
                "target_ancestry_registry_sha256": content_sha256(frozen_ancestry["target_identities"]),
                "input_features_sha256": content_sha256({"record_id": record_id, "group_id": group_id, "values": {name: values[name] for name in required}, "quantities": {name: given[name].to_dict() for name in required}, "provenance": {name: origin[name] for name in required}}),
                "outside_training_range": [name for name in required if not model["training_domain"][name]["minimum"] <= canonical[name] <= model["training_domain"][name]["maximum"]],
                "uncertainty": {"scope": "conditional_parameter_only_not_calibrated_coverage", "lower": interval[0], "upper": interval[1], "measurement_plus_prediction_interval": None},
                "stoichiometric_contributions": {species: coefficient * reported for species, coefficient in model["reaction"]["stoichiometry"].items()},
            }
        )
    content_sha256(predictions)
    return predictions


def _group_losses(data: _Prepared, indices: np.ndarray, predictions: list[dict], scale: float) -> list[dict]:
    by_id = {prediction["record_id"]: prediction["predicted_rate_canonical"] for prediction in predictions}
    result = []
    for group in sorted(set(data.groups[indices])):
        members = [int(i) for i in indices if data.groups[i] == group]
        errors = [abs(by_id[data.records[i]["record_id"]] - data.reaction["direction"] * data.values[data.target][i]) for i in members]
        result.append({"group_id": group, "role": data.contract["permitted_groups"][group], "observation_count": len(members), "mae_canonical": float(np.mean(errors)), "scaled_mae": float(np.mean(errors)) / scale, "training_iqr_canonical": scale})
    return result


def _r_squared(data: _Prepared, predictions: list[dict]) -> float | None:
    if len(predictions) != len(data.records):
        return None
    by_id = {prediction["record_id"]: prediction["predicted_rate_canonical"] for prediction in predictions}
    y = data.values[data.target] * data.reaction["direction"]
    if _constant(y, data.settings.constant_target_rtol):
        return None
    weights = _weights(data.groups)
    variance = float(weights @ ((y - weights @ y) ** 2))
    errors = np.asarray([by_id[record["record_id"]] for record in data.records]) - y
    return 1.0 - float(weights @ (errors * errors)) / variance


def _failure(error: _FitFailure, stage: str, fold_id: str | None = None) -> dict:
    return {"stage": stage, "fold_id": fold_id, "code": error.code, "message": str(error), "diagnostics": error.diagnostics}


def fit_rate_laws(measurements: Mapping[str, Any], protocol: Mapping[str, Any], *, catalogue: Mapping[str, Any] | None = None) -> dict:
    data = _prepare(measurements, protocol, catalogue)
    contract = data.contract
    result = {
        "schema_version": "native_rate_fit.v1",
        "status": "no_identified_candidate",
        "evidence_kind": data.evidence_kind,
        "claim_scope": "synthetic_engineering_only" if data.evidence_kind == "synthetic_engineering" else "permitted_measurement_training_only_not_independent_validation",
        "input_sha256": content_sha256(measurements),
        "protocol_sha256": content_sha256(protocol),
        "protocol_id": contract["protocol_id"],
        "measurement_ancestry": copy.deepcopy(data.measurement_ancestry),
        "measurement_schema": {"target": data.target, "quantities": {name: spec.to_dict() for name, spec in data.quantities.items()}, "canonical_units": {name: spec.canonical_unit for name, spec in data.quantities.items()}},
        "resolved_settings": asdict(data.settings),
        "selection": {**copy.deepcopy(dict(contract["selection"])), "selected_candidate_id": None, "tie_break": "fewer_parameters_then_model_id", "fit_groups": data.fit_groups, "folds": copy.deepcopy(data.folds)},
        "candidates": [],
        "model": None,
        "predictions": {"out_of_fold": [], "refit": []},
        "baseline": {"model": "training_group_weighted_constant", "folds": [], "r_squared": None},
        "measurement_uncertainty": {
            "scope": "reported_source_uncertainty_separate_from_scoring_scale",
            "records": [{"record_id": record["record_id"], "canonical_standard_uncertainties": {name: value * data.quantities[name].factor for name, value in record.get("standard_uncertainties", {}).items() if name in data.quantities}} for record in data.records],
        },
        "limitations": [
            "The learner does not authenticate custody, registration, lineage or scientific authorization; the independent verifier owns those decisions.",
            "Conditional parameter/profile intervals are not measurement-plus-prediction intervals and carry no calibrated coverage claim.",
            "No errors-in-variables correction, active-enzyme calibration, intrinsic kcat, gene-copy conversion or product-allocation parameter is inferred.",
            "Search bounds, nonnegative orientation, enzyme proportionality and any supplied stoichiometry are assumptions, not learned invariants.",
            "Model selection scores use permitted groups only and are not independent final-test performance.",
        ],
    }
    if _constant(data.values[data.target], data.settings.constant_target_rtol):
        result["status"] = "refused_constant_target"
        result["candidates"] = [{"candidate_id": candidate["candidate_id"], "specification": copy.deepcopy(candidate), "status": "not_run_constant_target", "model": None, "parameters": _parameter_templates(candidate), "folds": [], "failures": [{"stage": "preflight", "code": "constant_target", "message": "constant targets do not establish kinetic parameter recovery"}], "cv_score": None, "selection_score": None} for candidate in contract["candidates"]]
        content_sha256(result)
        return result
    prediction_quantities = {name: spec.to_dict() for name, spec in data.quantities.items() if name != data.target}
    baseline_predictions = []
    for fold in data.folds:
        train, validation = _indices(data, fold["fit_groups"]), _indices(data, fold["validation_groups"])
        try:
            scale = _iqr(data.values[data.target][train], data.groups[train], data.settings.constant_target_rtol)
            mean = float(_weights(data.groups[train]) @ data.values[data.target][train]) * data.reaction["direction"]
            predictions = [{"record_id": data.records[int(i)]["record_id"], "group_id": data.records[int(i)]["group_id"], "predicted_rate_canonical": mean} for i in validation]
            baseline_predictions.extend(predictions)
            result["baseline"]["folds"].append({"fold_id": fold["fold_id"], "status": "fitted", "fit_groups": fold["fit_groups"], "group_losses": _group_losses(data, validation, predictions, scale), "predictions": predictions})
        except _FitFailure as error:
            result["baseline"]["folds"].append({"fold_id": fold["fold_id"], "status": "failed", "failure": _failure(error, "baseline", fold["fold_id"])})
    result["baseline"]["r_squared"] = _r_squared(data, baseline_predictions)
    for candidate in contract["candidates"]:
        report = {"candidate_id": candidate["candidate_id"], "specification": copy.deepcopy(candidate), "status": "failed", "model": None, "parameters": _parameter_templates(candidate), "folds": [], "failures": [], "cv_score": None, "selection_score": None, "parameter_count": 1 + len(_shape_names(candidate))}
        result["candidates"].append(report)
        try:
            _check_candidate(candidate, data)
        except _FitFailure as error:
            report["failures"].append(_failure(error, "candidate_validation"))
            continue
        for fold in data.folds:
            fit_indices, validation_indices = _indices(data, fold["fit_groups"]), _indices(data, fold["validation_groups"])
            fold_report = {"fold_id": fold["fold_id"], "fit_groups": fold["fit_groups"], "validation_groups": fold["validation_groups"], "status": "failed", "parameters": _parameter_templates(candidate), "normalization": None, "predictions": [], "group_losses": []}
            report["folds"].append(fold_report)
            try:
                fitted = _fit_one(data, candidate, fit_indices)
                predictions = predict_rate_law(fitted["model"], _feature_records(data, validation_indices), prediction_quantities)
                fold_report.update(fitted)
                fold_report.update({"status": "fitted", "predictions": predictions, "group_losses": _group_losses(data, validation_indices, predictions, fitted["normalization"]["training_group_mean_iqr"])})
            except (_FitFailure, ContractError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
                if not isinstance(error, _FitFailure):
                    error = _FitFailure("numerical_or_prediction_failure", str(error))
                report["failures"].append(_failure(error, "cross_validation", fold["fold_id"]))
                fold_report["normalization"] = error.diagnostics.get("normalization")
                fold_report["diagnostics"] = error.diagnostics
        try:
            fitted = _fit_one(data, candidate, _indices(data, data.fit_groups))
            report.update(fitted)
        except (_FitFailure, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
            if not isinstance(error, _FitFailure):
                error = _FitFailure("numerical_failure", str(error))
            report["failures"].append(_failure(error, "refit"))
            report["diagnostics"] = error.diagnostics
        if all(fold["status"] == "fitted" for fold in report["folds"]):
            losses = [loss for fold in report["folds"] for loss in fold["group_losses"]]
            report["cv_score"] = float(np.mean([loss["scaled_mae"] for loss in losses]))
            selected_losses = losses if contract["selection"]["mode"] == "group_cross_validation" else [loss for loss in losses if loss["role"] == "development"]
            report["selection_score"] = float(np.mean([loss["scaled_mae"] for loss in selected_losses]))
            report["out_of_fold_r_squared"] = _r_squared(data, [prediction for fold in report["folds"] for prediction in fold["predictions"]])
        if not report["failures"] and report["model"] is not None:
            report["status"] = "eligible"
    eligible = [report for report in result["candidates"] if report["status"] == "eligible"]
    if eligible:
        best = min(report["selection_score"] for report in eligible)
        tied = [report for report in eligible if report["selection_score"] <= best + contract["selection"]["tie_tolerance"]]
        selected = min(tied, key=lambda report: (report["parameter_count"], report["candidate_id"]))
        result["status"] = "fitted"
        result["model"] = selected["model"]
        result["selection"]["selected_candidate_id"] = selected["candidate_id"]
        result["selection"]["score"] = selected["selection_score"]
        result["predictions"]["out_of_fold"] = [prediction for fold in selected["folds"] for prediction in fold["predictions"]]
        result["predictions"]["refit"] = predict_rate_law(selected["model"], _feature_records(data, _indices(data, data.fit_groups)), prediction_quantities)
    content_sha256(result)
    return result
