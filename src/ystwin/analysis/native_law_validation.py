from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
import hashlib
import hmac
import json
import math
from pathlib import Path
from statistics import NormalDist, mean
from typing import Literal


PROTOCOL_ID = "native-law-v2-preregistered-01"
GRADE = "independently_validated_scoped_empirical_native_response"
RECEIPT_KINDS = ("registry", "dataset", "lineage", "freeze", "observations", "predictions", "support", "controls")
MODELS = ("candidate", "training_constant", "affine_inputs", "clock_only", "source_template",
          "simple_mechanism", "dependency_ablated", "input_reassigned")
PROTOCOL_PATH = "data/native_law_v2/protocol.json"
DOMAIN_PATH = "data/native_law_v2/claim_domain.json"


@dataclass(frozen=True)
class SealedArtifact:
    path: str
    sha256: str


@dataclass(frozen=True)
class Authority:
    issuer: str
    secret: bytes
    kinds: frozenset[str]
    domain: Literal["real", "fixture"]


@dataclass(frozen=True)
class ValidationTrust:
    authorities: tuple[Authority, ...]
    registry_heads: Mapping[str, str]


@dataclass(frozen=True)
class EvidenceBundle:
    registry: SealedArtifact | None = None
    dataset: SealedArtifact | None = None
    lineage: SealedArtifact | None = None
    freeze: SealedArtifact | None = None
    observations: SealedArtifact | None = None
    predictions: SealedArtifact | None = None
    support: SealedArtifact | None = None
    controls: SealedArtifact | None = None


@dataclass(frozen=True)
class Gate:
    name: str
    state: Literal["pass", "fail", "missing"]
    reasons: tuple[str, ...] = ()
    metrics: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class Assessment:
    candidate_id: str | None
    grade: str
    authorized: bool
    status: str
    scope: str
    gates: tuple[Gate, ...]
    evidence_sha256: Mapping[str, str]


class EvidenceError(ValueError):
    def __init__(self, reason: str, metrics: Mapping[str, object] | None = None):
        super().__init__(reason)
        self.metrics = {} if metrics is None else metrics


class MissingEvidence(EvidenceError):
    pass


def _need(condition: bool, reason: str) -> None:
    if not condition:
        raise EvidenceError(reason)


def _number(value: object) -> float:
    _need(type(value) in (int, float) and math.isfinite(value), "a finite numeric measurement is required")
    return float(value)


def _text(value: object) -> str:
    _need(isinstance(value, str) and bool(value.strip()) and value == value.strip(), "nonempty unpadded identity required")
    return value


def _sha(value: object) -> str:
    _need(isinstance(value, str) and len(value) == 64 and set(value) <= set("0123456789abcdef"), "invalid SHA256 identity")
    return value


def _strings(values: object, *, empty: bool = False) -> set[str]:
    _need(isinstance(values, list) and (empty or bool(values)), "explicit identity list required")
    result = {_text(value) for value in values}
    _need(len(result) == len(values), "duplicate identity in inventory")
    return result


def _index(rows: object, key: str) -> dict[str, dict]:
    _need(isinstance(rows, list) and bool(rows), "complete nonempty record inventory required")
    result = {}
    for row in rows:
        _need(isinstance(row, dict), "typed record required")
        identifier = _text(row[key])
        _need(identifier not in result, f"duplicate {key}")
        result[identifier] = row
    return result


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _pairs(items: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in items:
        _need(key not in result, "duplicate JSON field")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise EvidenceError(f"nonfinite JSON constant: {value}")


def _parse(raw: bytes) -> dict:
    value = json.loads(raw, object_pairs_hook=_pairs, parse_constant=_invalid_constant)
    _need(isinstance(value, dict), "JSON artifact must be an object")
    _canonical(value)
    return value


def _ref(value: object) -> SealedArtifact:
    _need(isinstance(value, dict) and set(value) == {"path", "sha256"}, "exact sealed artifact reference required")
    return SealedArtifact(_text(value["path"]), _sha(value["sha256"]))


class _Reader:
    def __init__(self, root: Path, trust: ValidationTrust):
        self.root = root.resolve(strict=True)
        self.authorities = {}
        _need(isinstance(trust, ValidationTrust), "trusted verifier configuration required")
        _need(bool(trust.authorities), "independent authority configuration is missing")
        secrets = set()
        for authority in trust.authorities:
            _need(isinstance(authority, Authority), "typed independent authority required")
            _text(authority.issuer)
            _need(authority.issuer not in self.authorities, "duplicate independent authority identity")
            _need(isinstance(authority.secret, bytes) and len(authority.secret) >= 32
                  and authority.secret not in secrets, "independent authority keys must be distinct and at least 32 bytes")
            _need(authority.domain in {"real", "fixture"} and bool(authority.kinds)
                  and set(authority.kinds) <= set(RECEIPT_KINDS), "invalid independent authority capabilities")
            secrets.add(authority.secret)
            self.authorities[authority.issuer] = authority
        self.issuers = {}
        self.cache = {}
        self.source_tables = {}

    def raw(self, path: str) -> bytes:
        relative = Path(_text(path))
        _need(not relative.is_absolute() and relative.parts and ".." not in relative.parts,
              "artifact path must stay within the explicit verification root")
        current = self.root
        for part in relative.parts:
            current = current / part
            _need(not current.is_symlink(), "symlink artifact path is not authoritative")
        _need(current.is_file(), f"missing artifact file: {path}")
        return current.read_bytes()

    def artifact(self, ref: SealedArtifact | dict) -> bytes:
        ref = _ref(ref) if isinstance(ref, dict) else ref
        _need(isinstance(ref, SealedArtifact), "typed sealed artifact required")
        _sha(ref.sha256)
        key = (ref.path, ref.sha256)
        if key not in self.cache:
            raw = self.raw(ref.path)
            _need(hashlib.sha256(raw).hexdigest() == ref.sha256, f"artifact digest mismatch: {ref.path}")
            self.cache[key] = raw
        return self.cache[key]

    def receipt(self, ref: SealedArtifact, kind: str) -> dict:
        document = _parse(self.artifact(ref))
        _need(set(document) == {"schema_version", "issuer", "kind", "receipt_id", "payload", "mac"},
              "authentication requires an exact receipt envelope, not pass flags")
        _need(type(document["schema_version"]) is int and document["schema_version"] == 1
              and document["kind"] == kind, "receipt schema/kind mismatch")
        authority = self.authorities.get(document["issuer"])
        _need(authority is not None and kind in authority.kinds, f"authentication authority cannot issue {kind}")
        _text(document["receipt_id"])
        mac = _sha(document["mac"])
        unsigned = {key: value for key, value in document.items() if key != "mac"}
        expected = hmac.new(authority.secret, _canonical(unsigned), hashlib.sha256).hexdigest()
        _need(hmac.compare_digest(mac, expected), f"receipt authentication failed: {kind}")
        _need(isinstance(document["payload"], dict), "receipt payload must contain typed evidence")
        self.issuers[kind] = authority.issuer
        return document["payload"]


def _close(actual: object, expected: object, name: str) -> None:
    _need(math.isclose(_number(actual), _number(expected), rel_tol=1e-8, abs_tol=1e-10), f"{name} numerical mismatch")


def _pointer(value: object, pointer: object) -> object:
    _need(isinstance(pointer, str) and (pointer == "" or pointer.startswith("/")), "source locator must be an RFC6901 JSON pointer")
    if pointer == "":
        return value
    for token in pointer[1:].split("/"):
        _need(all(part and part[0] in "01" for part in token.split("~")[1:]), "invalid JSON pointer escape")
        token = token.replace("~1", "/").replace("~0", "~")
        if isinstance(value, list):
            _need(token.isascii() and token.isdigit() and str(int(token)) == token and int(token) < len(value),
                  "source array locator is outside its exact indexed domain")
            value = value[int(token)]
        else:
            _need(isinstance(value, dict), "source locator does not resolve a structured observation")
            value = value[token]
    return value


def _source_fields(reader: _Reader, source: dict, identifier: str, decoder_key: str = "decoder") -> tuple[dict, object]:
    decoder = source.get(decoder_key)
    if not isinstance(decoder, dict) or decoder.get("kind") != "json_records":
        raise MissingEvidence(f"no implemented source adapter for {source['source_id']} {decoder_key}; hashes do not bind measurement values")
    _need(isinstance(decoder["fields"], dict) and bool(decoder["fields"]), "source adapter requires explicit field locators")
    ref = _ref(source["artifact"])
    key = (ref.sha256, _canonical(decoder))
    if key not in reader.source_tables:
        document = _parse(reader.artifact(ref))
        origin = _pointer(document, decoder["origin_pointer"])
        custodian = reader.authorities[reader.issuers["dataset"]]
        _need(origin == "experimental_measurement" or (origin == "synthetic_software_fixture" and custodian.domain == "fixture"),
              "source origin is not a real experimental measurement; reclassifying authority metadata cannot authorize teacher/fixture labels")
        table = _pointer(document, decoder["table_pointer"])
        _need(isinstance(table, list) and bool(table), "source adapter requires an actual nonempty record table")
        indexed = {}
        for row in table:
            row_id = _text(_pointer(row, decoder["key_pointer"]))
            _need(row_id not in indexed, "duplicate physical source-table identity")
            indexed[row_id] = row
        unit_key = "units_pointer" if decoder_key == "input_decoder" else "unit_pointer"
        units = _pointer(document, decoder[unit_key])
        reader.source_tables[key] = (indexed, units)
    indexed, units = reader.source_tables[key]
    _need(identifier in indexed, "registered observation identity is absent from its actual source table")
    values = {_text(name): _pointer(indexed[identifier], pointer) for name, pointer in decoder["fields"].items()}
    return values, units


def _source_equal(actual: object, expected: object, name: str) -> None:
    if expected is None:
        _need(actual is None, f"{name} source/projection missingness mismatch")
    elif type(expected) in (int, float):
        _need(_number(actual) == _number(expected), f"{name} source/projection value mismatch")
    else:
        _need(type(actual) is type(expected) and actual == expected, f"{name} source/projection identity or status mismatch")


def _quantile(values: list[float], probability: float) -> float:
    values = sorted(values)
    position = (len(values) - 1) * probability
    index = int(position)
    fraction = position - index
    return values[index] * (1 - fraction) + values[min(index + 1, len(values) - 1)] * fraction


def _wilson(wins: int, count: int) -> tuple[float, float]:
    z = NormalDist().inv_cdf(0.975)
    proportion = wins / count
    denominator = 1 + z * z / count
    center = (proportion + z * z / (2 * count)) / denominator
    radius = z * math.sqrt(proportion * (1 - proportion) / count + z * z / (4 * count * count)) / denominator
    return center - radius, center + radius


def _sign(wins: int, count: int) -> dict:
    return {"wins": wins, "groups": count, "p_value": sum(math.comb(count, k) for k in range(wins, count + 1)) / (1 << count),
            "winning_fraction_wilson_95": list(_wilson(wins, count))}


def _dependencies(group: dict) -> set[str]:
    return {group["group_id"], *_strings(group["biological_ids"]), *_strings(group["dependency_ids"])}


def _partition(groups: dict[str, dict], roles: set[str]) -> set[str]:
    return {gid for gid, group in groups.items() if group["role"] in roles}


def _group_order(candidate: str, group: str) -> str:
    return hashlib.sha256((PROTOCOL_ID + candidate + group).encode()).hexdigest()


def _data_metadata(data: dict, domain: dict, protocol: dict) -> dict:
    _need(data["source_doi"] == domain["source_doi"], "publication provenance differs from the named candidate")
    measurement = data["measurement"]
    _need(measurement["origin"] == "experimental_observation", "generated teacher targets are not real native observations")
    _need(measurement["observable"] == domain["observable"] and measurement["unit"] == domain["canonical_unit"],
          "observable identity/unit is outside the preregistered claim")
    _need(measurement["uncertainty_kind"] == "standard_uncertainty", "source measurement uncertainty is missing or not standard uncertainty")
    sources = _index(data["sources"], "source_id")
    for source in sources.values():
        _ref(source["artifact"])
        _ref(source["methods"])
        _text(source["accession"])
        _text(source["assay_id"])
        _need(source["kind"] in {"native_imaging", "nuclear_marker_imaging", "native_flux_metabolite_protein", "native_activity_assay"},
              "custody must identify a real experimental measurement method")
    _need(measurement["calibration_source_id"] in sources, "source observation calibration is missing")
    context = data["context"]
    for key in ("species", "strain", "protein", "medium", "platform"):
        _text(context[key])
    _need(context["species"] == "Saccharomyces cerevisiae" and _number(context["temperature_K"]) > 0,
          "native species/temperature context is invalid")
    expected_protein = "Sfp1" if domain["observable"] == "Sfp1_max5_over_median" else "Cdc19"
    _need(context["protein"] == expected_protein, "named native protein/isoform context does not match")
    groups = _index(data["groups"], "group_id")
    _need(set(data["input_units"]) == set(domain["input_quantities"]), "source-bound input units are incomplete")
    for unit in data["input_units"].values():
        _text(unit)
    seen = set()
    for group in groups.values():
        _need(group["role"] in {"train", "development", "final", "support", "sham"}, "unknown biological split role")
        dependencies = _dependencies(group)
        _need(not seen.intersection(dependencies), "duplicated biological dependency: merge complete components before splitting")
        seen.update(dependencies)
        _need(group["context"] == context, "independent support and target require a matched verified context")
        _number(group["acquisition_clock_s"])
        _need(set(group["inputs"]) == set(domain["input_quantities"]), "registered input quantities are incomplete")
        for value in group["inputs"].values():
            _number(value)
        _need(set(group["input_locators"]) == set(group["inputs"]), "jointly observed input provenance is incomplete")
        for name, locator in group["input_locators"].items():
            _need(locator["source_id"] in sources and locator["unit"] == data["input_units"][name], "source-bound input units/provenance do not match")
            _text(locator["locator"])
        schedule = group["input_schedule"]
        _need(isinstance(schedule, list) and len(schedule) >= 2, "source-bound input schedule is missing")
        times = [_number(event["time_s"]) for event in schedule]
        _need(times == sorted(set(times)), "input schedule must be strictly time ordered")
        for event in schedule:
            _number(event["value"])
            _text(event["quantity"])
            _text(event["unit"])
        if group["role"] == "sham":
            _need(len({(event["quantity"], event["value"], event["unit"]) for event in schedule}) == 1,
                  "sham experiment contains an input change")
            if expected_protein == "Sfp1":
                _need(group["inputs"]["initial_glucose_percent"] == group["inputs"]["final_glucose_percent"],
                      "sham label cannot replace an actual no-input-change experiment")
    selection = protocol["selection"]
    for role, key in (("train", "minimum_train_independent_groups"), ("development", "minimum_development_independent_groups"),
                      ("final", "minimum_final_independent_groups"), ("support", "minimum_support_independent_groups")):
        _need(len(_partition(groups, {role})) >= selection[key], f"insufficient independent {role} groups")
    _need(len(_partition(groups, {"sham"})) >= 2, "independent no-input-change sham groups are missing")
    candidates = sorted(_partition(groups, {"train", "development", "final"}), key=lambda g: _group_order(data["candidate_id"], g))
    nfinal, ndev = max(6, math.ceil(0.25 * len(candidates))), max(2, math.ceil(0.20 * len(candidates)))
    _need(set(candidates[:nfinal]) == _partition(groups, {"final"})
          and set(candidates[nfinal:nfinal + ndev]) == _partition(groups, {"development"})
          and set(candidates[nfinal + ndev:]) == _partition(groups, {"train"}),
          "whole-experiment partition differs from the prespecified metadata-only split")
    clocks = {}
    for gid in _partition(groups, {"final"}):
        group = groups[gid]
        clocks.setdefault(group["acquisition_clock_s"], set()).add(_canonical(group["input_schedule"]))
    crossed = [values for values in clocks.values() if len(values) >= 2]
    _need(len(crossed) >= 2 and any(a & b for i, a in enumerate(crossed) for b in crossed[i + 1:]),
          "real input variation is confounded with absolute acquisition clock; replay invariance cannot establish response")
    records = _index(data["records"], "record_id")
    physical = set()
    for record in records.values():
        _need(record["group_id"] in groups and groups[record["group_id"]]["role"] != "support", "target record belongs to an unknown or support group")
        _need(record["source_id"] in sources and record["observable"] == domain["observable"], "source-bound target observable is incorrect")
        _text(record["locator"])
        identifier = _text(record["physical_observation_id"])
        _need(identifier not in physical, "duplicate physical observation, including a source-file copy")
        physical.add(identifier)
        _need(record["status"] in {"quantified", "missing", "censored"}, "source missingness/censoring must be declared before release")
        _need(record["source_unit"] == measurement["adapter"]["source_unit"], "inconsistent source measurement units")
        _number(record["time_s"])
    support = _index(data["support_records"], "record_id")
    _need(not set(support) & set(records), "support reuses a target observation identity")
    target_sources = {record["source_id"] for record in records.values()}
    target_artifacts = {sources[key]["artifact"]["sha256"] for key in target_sources}
    target_assays = {sources[key]["assay_id"] for key in target_sources}
    support_groups = []
    for record in support.values():
        gid, sid = record["group_id"], record["source_id"]
        _need(gid in groups and groups[gid]["role"] == "support", "support is not in a reserved independent group")
        _need(sid in sources and sid not in target_sources and sources[sid]["artifact"]["sha256"] not in target_artifacts
              and sources[sid]["assay_id"] not in target_assays, "support must be independently acquired, not the same target signal")
        expected_kind = "nuclear_marker_imaging" if expected_protein == "Sfp1" else "native_activity_assay"
        _need(sources[sid]["kind"] == expected_kind, "independent support measurement kind is not the preregistered assay")
        _text(record["low_locator"])
        _text(record["high_locator"])
        identifiers = _strings(record["physical_observation_ids"])
        _need(len(identifiers) == 2 and not identifiers & physical, "support physical observations are reused")
        physical.update(identifiers)
        support_groups.append(gid)
    _need(len(set(support_groups)) == len(support_groups)
          and set(support_groups) == _partition(groups, {"support"}), "support requires one registered paired contrast per complete biological group")
    return {"groups": groups, "records": records, "support_records": support, "sources": sources}


def _registry(registry: dict, meta: dict, freeze: dict, refs: dict[str, SealedArtifact], protocol_sha: str, domain_sha: str) -> dict:
    _need(registry["protocol_id"] == PROTOCOL_ID and registry["protocol_sha256"] == protocol_sha
          and registry["claim_domain_sha256"] == domain_sha, "pre-fit registration does not bind the current protocol and claim domain")
    _need(registry["artifacts"] == {key: ref.sha256 for key, ref in refs.items() if key != "registry"},
          "registry does not bind every sealed evidence artifact")
    events = registry["events"]
    _need(isinstance(events, list) and [event["sequence"] for event in events] == list(range(1, len(events) + 1)),
          "registry event sequence must be complete and monotone")
    required = ("protocol_registered", "dataset_sealed", "split_registered", "model_library_registered", "development_opened", "fit_started",
                "selection_frozen", "predictions_sealed", "final_outcomes_opened", "support_outcomes_opened", "evaluation_sealed")
    named = {}
    for event in events:
        kind = event["event"]
        _need(kind in {*required, "development_attempt"}, "unknown/unresolved registry event")
        if kind != "development_attempt":
            _need(kind not in named, "repeated registration/freeze/final release consumes the holdout")
            named[kind] = event
    _need(set(named) == set(required), "registry lacks complete registration/freeze/outcome-access history")
    positions = [named[kind]["sequence"] for kind in required]
    _need(positions == sorted(positions), "protocol, split and selection freeze must precede all final outcome access")
    _need(named["protocol_registered"]["protocol_sha256"] == protocol_sha
          and named["protocol_registered"]["claim_domain_sha256"] == domain_sha, "protocol was not registered before fitting")
    for event, artifact in (("dataset_sealed", "dataset"), ("selection_frozen", "freeze"), ("predictions_sealed", "predictions")):
        _need(named[event]["sha256"] == refs[artifact].sha256, f"{event} bound a different artifact")
    _need(named["model_library_registered"]["sha256"] == freeze["config_artifacts"]["model_library"]["sha256"],
          "functional library and comparator mechanisms were not fixed before fitting/data access")
    groups = meta["groups"]
    _need(named["split_registered"]["partitions"] == {gid: g["role"] for gid, g in groups.items()}, "registered biological split changed")
    for event, roles in (("development_opened", {"train", "development"}), ("final_outcomes_opened", {"final", "sham"}),
                         ("support_outcomes_opened", {"support"})):
        _need(_strings(named[event]["group_ids"]) == _partition(groups, roles), "outcome release inventory is partial or changed")
    attempts = [event for event in events if event["event"] == "development_attempt"]
    _need(len(attempts) == len(freeze["attempts"])
          and {e["attempt_id"] for e in attempts} == {a["attempt_id"] for a in freeze["attempts"]},
          "all development/selection attempts must be retained in the authoritative log")
    _need(all(named["fit_started"]["sequence"] < e["sequence"] < named["selection_frozen"]["sequence"] for e in attempts),
          "selection attempts occurred outside the pre-final development window")
    _need(named["evaluation_sealed"]["artifacts"] == {key: refs[key].sha256 for key in ("observations", "support", "controls")},
          "final evaluation is incomplete or rebound to a later result")
    _need(registry["run_id"] == freeze["run_id"], "more than one selected final run")
    history = registry["access_history"]
    _need(isinstance(history, list) and bool(history), "current complete ancestral access ledger required")
    for event, roles in (("development_opened", {"train", "development"}), ("final_outcomes_opened", {"final", "sham"}),
                         ("support_outcomes_opened", {"support"})):
        wanted = set().union(*(_dependencies(g) for g in groups.values() if g["role"] in roles))
        matching = []
        for access in history:
            accessed = _strings(access["dependency_ids"])
            if not wanted & accessed:
                continue
            current = (access["protocol_id"] == PROTOCOL_ID and access["candidate_id"] == registry["candidate_id"]
                       and access["run_id"] == registry["run_id"] and access["purpose"] == event
                       and access["event_sequence"] == named[event]["sequence"])
            if event != "development_opened":
                _need(current, "heldout biological dependency was previously exposed/consumed, including under another protocol")
            if current:
                matching.append(accessed)
        _need(len(matching) == 1 and matching[0] == wanted, "complete single-release custody record is missing or repeated")
    return {"final_releases": 1, "support_releases": 1, "development_attempts": len(attempts), "candidate_alpha": 0.025}


def _lineage(lineage: dict, freeze: dict, meta: dict, reader: _Reader) -> dict:
    dependencies = {_ref(ref).sha256 for ref in freeze["dependencies"]}
    _need(_strings(lineage["roots"]) == dependencies, "ancestral audit roots do not cover all frozen dependencies")
    nodes = {}
    for node in lineage["nodes"]:
        ref = _ref(node["artifact"])
        _need(ref.sha256 not in nodes, "duplicate ancestral artifact node")
        reader.artifact(ref)
        _need(node["role"] in {"code", "model_dependency", "measurement_release", "source_prior", "external_calibration"},
              "unknown or unaudited ancestral dependency")
        _text(node["source_id"])
        _strings(node["locators"])
        _strings(node["parents"], empty=True)
        _strings(node["exposed_dependency_ids"], empty=True)
        nodes[ref.sha256] = node
    _need(dependencies <= set(nodes), "ancestral dependency graph is not closed")
    forbidden = set().union(*(_dependencies(g) for g in meta["groups"].values() if g["role"] in {"final", "support", "sham"}))
    visited, active = set(), set()

    def visit(identifier: str) -> None:
        _need(identifier in nodes, "ancestral dependency graph has an unknown leaf")
        _need(identifier not in active, "ancestral dependency graph is cyclic")
        if identifier in visited:
            return
        active.add(identifier)
        node = nodes[identifier]
        _need(not set(node["exposed_dependency_ids"]) & forbidden, "ancestral fit/selection/normalization exposure overlaps reserved experiments")
        _need(node["parents"] or node["role"] in {"code", "measurement_release", "source_prior", "external_calibration"},
              "unresolved ancestral model leaf")
        if node["role"] == "measurement_release":
            _need(bool(node["exposed_dependency_ids"]), "ancestral measurement exposure inventory is missing")
        for parent in node["parents"]:
            visit(parent)
        active.remove(identifier)
        visited.add(identifier)

    for root in dependencies:
        visit(root)
    _need(visited == set(nodes), "ancestral audit contains disconnected/unbound evidence")
    return {"audited_artifacts": len(nodes), "reserved_ancestral_overlap": 0}


def _freeze(freeze: dict, data: dict, meta: dict, domain: dict, reader: _Reader) -> dict:
    _need(freeze["context"] == data["context"], "frozen claim context differs from source context")
    _need(freeze["gauge_scope"] == "observable_only" and bool(_strings(freeze["unidentified_directions"])),
          "latent measurement-scale gauge must remain explicitly unidentified")
    train = _partition(meta["groups"], {"train"})
    development = _partition(meta["groups"], {"development"})
    allowed = train | development
    normalization = freeze["normalization"]
    _need(normalization["method"] == "training_group_iqr" and _strings(normalization["fit_groups"]) == train
          and _number(normalization["offset"]) == 0 and _number(normalization["scale"]) > 0,
          "normalization must use the complete training partition only, with a positive frozen IQR")
    _need(set(normalization["training_group_means"]) == train, "normalization training means are incomplete")
    adapter = freeze["observation_adapter"]
    _need(adapter == data["measurement"]["adapter"] and adapter["canonical_unit"] == domain["canonical_unit"],
          "observation adapter differs from authoritative source units/normalization")
    conversions = {"dimensionless_ratio": 1.0} if domain["canonical_unit"] == "dimensionless_ratio" else {
        "mmol/gDW/h": 1.0, "mol/gDW/h": 1000.0, "umol/gDW/h": 0.001, "mmol/gDW/s": 3600.0,
    }
    _need(adapter["source_unit"] in conversions and _number(adapter["offset"]) == 0
          and _number(adapter["factor"]) == conversions[adapter["source_unit"]], "unknown units or outcome-derived observation normalization")
    _need(adapter["kind"] == ("ratio" if domain["canonical_unit"] == "dimensionless_ratio" else "linear"),
          "observation measurement adapter does not match the biological observable")
    support = freeze["support_adapter"]
    _need(_strings(support["fit_groups"]) <= allowed and _number(support["factor"]) > 0, "support calibration used reserved observations or an invalid scale")
    _number(support["offset"])
    _need(all(r["unit"] == support["unit"] for r in meta["support_records"].values()), "support assay units/calibration are unresolved")
    error = freeze["error_model"]
    _need(error["kind"] == "gaussian_prediction_interval" and _number(error["nominal_coverage"]) == 0.9
          and _strings(error["fit_groups"]) <= allowed, "a pre-final statistical prediction interval/error budget is required")
    for field_name in ("prediction_sd", "measurement_sd", "support_prediction_sd"):
        _need(_number(error[field_name]) > 0, "frozen uncertainty components must be positive")
    _need(set(freeze["input_domain"]) == set(domain["input_quantities"]), "identifiable input domain is incomplete")
    for name, bounds in freeze["input_domain"].items():
        _need(_text(bounds["unit"]) == data["input_units"][name], "frozen input units are not the authoritative measurement units")
        lower, upper = _number(bounds["lower"]), _number(bounds["upper"])
        _need(lower <= upper, "invalid frozen input domain")
        _need(all(lower <= _number(g["inputs"][name]) <= upper for g in meta["groups"].values()),
              "registered experiment lies outside the frozen claim domain; do not drop or expand it")
    dependencies = {_ref(ref).sha256 for ref in freeze["dependencies"]}
    _need(len(dependencies) == len(freeze["dependencies"]), "duplicate frozen dependency")
    codes = {_ref(ref).sha256 for ref in freeze["code_artifacts"]}
    _need(bool(codes) and codes <= dependencies, "frozen inference/scoring code dependencies are missing")
    source_spec = freeze["source_template_specification"]
    _need(source_spec["parameter_policy"] == "published_fixed_except_observation_gain", "source template requires a fixed published parameter policy")
    _text(source_spec["source_reference"])
    _text(source_spec["source_locator"])
    source_ref = _ref(source_spec["artifact"])
    _need(source_ref.sha256 in dependencies, "source template parameters/structure lack ancestral audit coverage")
    reader.artifact(source_ref)
    library = _index(freeze["library"], "form_id")
    signatures = {_canonical(entry["functional_form"]) for entry in library.values()}
    _need(len(signatures) == len(library) and len(signatures) >= 2, "functional forms must differ in structure, not merely in their names")
    model_library = {"library": freeze["library"], "source_template_specification": source_spec,
                     "comparators": [{key: model[key] for key in ("model_id", "form_id", "feature_names", "parameter_count", "functional_form", "optimizer_budget")}
                                     for model in freeze["models"] if model["model_id"] != "candidate"]}
    expected_configs = {"normalization": normalization, "observation_adapter": adapter, "support_adapter": support,
                        "error_model": error, "model_library": model_library,
                        "selection_log": {"library": freeze["library"], "attempts": freeze["attempts"]}}
    _need(set(freeze["config_artifacts"]) == set(expected_configs), "frozen configuration artifacts are incomplete")
    for key, expected in expected_configs.items():
        ref = _ref(freeze["config_artifacts"][key])
        _need(ref.sha256 in dependencies and _parse(reader.artifact(ref)) == expected, f"{key} artifact differs from the frozen evidence")
    models = _index(freeze["models"], "model_id")
    _need(set(models) == set(MODELS), "all prespecified simpler and competing models must be frozen")
    for name, model in models.items():
        ref = _ref(model["artifact"])
        _need(ref.sha256 in dependencies, "model checkpoint missing from ancestral audit")
        saved = _parse(reader.artifact(ref))
        for key in ("model_id", "form_id", "feature_names", "parameter_count", "functional_form"):
            _need(saved[key] == model[key], "frozen model definition differs from checkpoint bytes")
        if name == "source_template":
            _need(saved["source_artifact_sha256"] == source_ref.sha256, "source template checkpoint does not bind the fixed source parameters/structure")
        _need(type(model["parameter_count"]) is int and model["parameter_count"] >= 0
              and len(saved["parameters"]) == model["parameter_count"], "model parameter count is not bound")
        for value in saved["parameters"]:
            _number(value)
        _need(saved["implementation_sha256"] in codes, "model implementation is outside frozen code ancestry")
        features = _strings(model["feature_names"], empty=True)
        if name == "clock_only":
            _need(features == {"acquisition_clock_s"}, "clock comparator must use absolute experimental clock only")
        else:
            _need(features <= set(domain["input_quantities"]), "candidate/comparator may not read absolute clock, identities or outcomes")
        if name == "training_constant":
            _need(not features, "training-constant baseline cannot use experiment inputs")
        dependency = "final_glucose_percent" if domain["canonical_unit"] == "dimensionless_ratio" else "fructose_1_6_bisphosphate"
        if name == "candidate":
            _need(dependency in features, "named functional dependency is absent from candidate")
        if name == "dependency_ablated":
            _need(dependency not in features, "dependency ablation retained the tested input")
        _need(_strings(model["fit_groups"]) <= allowed, "model fit used reserved observations")
        _need(type(model["optimizer_budget"]) is int and model["optimizer_budget"] > 0, "comparator fitting budget missing")
    _need(len({model["optimizer_budget"] for model in models.values()}) == 1, "comparator optimizer budgets are not matched")
    library = _index(freeze["library"], "form_id")
    _need(len(library) >= 2, "functional dependency cannot be inferred by fitting multipliers in only one supplied form")
    attempts = _index(freeze["attempts"], "attempt_id")
    _need(freeze["selected_attempt_id"] in attempts, "selected attempt is absent from the full development log")
    successes = set()
    for attempt in attempts.values():
        _need(attempt["form_id"] in library, "post-registration functional form was introduced")
        _need(type(attempt["seed"]) is int and type(attempt["evaluations"]) is int
              and 0 <= attempt["evaluations"] <= models["candidate"]["optimizer_budget"], "attempt seed/budget is unbound")
        _need(attempt["status"] in {"ok", "failed"}, "all unsuccessful attempts must be retained")
        if attempt["status"] == "failed":
            _text(attempt["reason"])
            continue
        checkpoint = _ref(attempt["checkpoint"])
        _need(checkpoint.sha256 in dependencies, "development candidate checkpoint not covered by ancestry")
        saved = _parse(reader.artifact(checkpoint))
        _need(saved["form_id"] == attempt["form_id"] and saved["parameter_count"] == library[attempt["form_id"]]["parameter_count"]
              and saved["feature_names"] == library[attempt["form_id"]]["feature_names"]
              and saved["functional_form"] == library[attempt["form_id"]]["functional_form"], "functional library is not bound to attempted checkpoint")
        successes.add(attempt["form_id"])
    _need(len(successes) >= 2, "at least two distinct successful real-data functional fits are required")
    selected = attempts[freeze["selected_attempt_id"]]
    _need(selected["status"] == "ok" and selected["checkpoint"] == models["candidate"]["artifact"], "selected candidate is not the frozen development attempt")
    return {"models": len(models), "functional_forms": len(library), "gauge_scope": "observable_only", "context": data["context"],
            "input_domain": freeze["input_domain"], "unidentified_directions": freeze["unidentified_directions"]}


def _measurements(observations: dict, data: dict, meta: dict, freeze: dict, reader: _Reader) -> dict:
    for source in meta["sources"].values():
        reader.artifact(source["artifact"])
        reader.artifact(source["methods"])
    rows = _index(observations["rows"], "record_id")
    _need(set(rows) == set(meta["records"]), "complete authoritative measurement inventory is required; do not drop failed/missing records")
    for gid, group in meta["groups"].items():
        for name, locator in group["input_locators"].items():
            original, units = _source_fields(reader, meta["sources"][locator["source_id"]], gid, "input_decoder")
            for key in ("group_id", "biological_ids", "dependency_ids", "context", "acquisition_clock_s", "input_schedule"):
                _source_equal(group[key], original[key], "experiment identity/design")
            _source_equal(group["inputs"][name], original[name], "input")
            _need(units[name] == data["input_units"][name] == locator["unit"], "input source/projection units mismatch")
    adapter = freeze["observation_adapter"]
    values = {}
    for rid, record in meta["records"].items():
        row = rows[rid]
        original, unit = _source_fields(reader, meta["sources"][record["source_id"]], rid)
        _need(unit == record["source_unit"], "measurement source/projection unit mismatch")
        required = {"raw_value", "raw_sd"}
        if adapter["kind"] == "ratio":
            required |= {"numerator", "denominator"}
        if meta["groups"][record["group_id"]]["role"] == "sham":
            required |= {"baseline_raw_value", "baseline_raw_sd"}
        if not required <= set(original):
            raise MissingEvidence("source adapter does not resolve every measurement/uncertainty component")
        for key in required | (set(original) & set(row)):
            _source_equal(row[key], original[key], "measurement")
        for key, expected in (("record_id", rid), ("group_id", record["group_id"]),
                              ("physical_observation_id", record["physical_observation_id"]), ("status", record["status"])):
            if key in original:
                _source_equal(original[key], expected, "measurement")
        if record["status"] != "quantified":
            _need(row["raw_value"] is None and row["raw_sd"] is None, "missing/censored measurements may not become zero or exact values")
            continue
        raw, sd = _number(row["raw_value"]), _number(row["raw_sd"])
        _need(sd > 0, "source standard measurement uncertainty must be positive")
        if adapter["kind"] == "ratio":
            numerator, denominator = _number(row["numerator"]), _number(row["denominator"])
            _need(numerator >= 0 and denominator > 0, "localization observation has an invalid raw intensity denominator")
            _close(raw, numerator / denominator, "raw localization ratio")
        value = raw * adapter["factor"] + adapter["offset"]
        _need(value >= data["physical"]["lower"], "measured native observable violates its source-backed physical bound")
        if data["physical"]["upper"] is not None:
            _need(value <= data["physical"]["upper"], "measured native observable exceeds its source-backed bound")
        values[rid] = {"value": value, "sd": sd * abs(adapter["factor"]), "raw": row}
    for gid, group in meta["groups"].items():
        if group["role"] != "support":
            _need(any(rid in values and row["group_id"] == gid for rid, row in meta["records"].items()),
                  "a registered biological group has no quantifiable target observations")
    return values


def _selection(freeze: dict, meta: dict, values: dict, reader: _Reader) -> dict:
    groups, records = meta["groups"], meta["records"]
    train = _partition(groups, {"train"})
    means = {gid: mean(values[rid]["value"] for rid, r in records.items() if r["group_id"] == gid and rid in values) for gid in train}
    scale = _quantile(list(means.values()), 0.75) - _quantile(list(means.values()), 0.25)
    _need(scale > 0, "training IQR is zero; no outcome-derived fallback scale is permitted")
    _close(freeze["normalization"]["scale"], scale, "normalization IQR")
    for gid, value in means.items():
        _close(freeze["normalization"]["training_group_means"][gid], value, "normalization training group mean")
    development = _partition(groups, {"development"})
    wanted = {rid for rid, r in records.items() if r["group_id"] in development and rid in values}
    scores = []
    for attempt in freeze["attempts"]:
        if attempt["status"] != "ok":
            continue
        predictions = attempt["development_predictions"]
        _need(set(predictions) == wanted, "complete whole-development-group predictions are required for every attempted form")
        losses = {gid: mean(abs(_number(predictions[rid]) - values[rid]["value"]) / scale for rid in wanted if records[rid]["group_id"] == gid)
                  for gid in development}
        saved = _parse(reader.artifact(_ref(attempt["checkpoint"])))
        scores.append((mean(losses.values()), saved["parameter_count"], saved["model_id"], attempt["attempt_id"]))
    minimum = min(row[0] for row in scores)
    selected = min((row for row in scores if row[0] <= minimum + 1e-12), key=lambda row: row[1:])
    _need(selected[3] == freeze["selected_attempt_id"], "selected form differs from the predeclared development-only selection rule")
    return {"recomputed_training_iqr": scale, "attempt_scores": [{"attempt_id": row[3], "mean_development_scaled_mae": row[0]} for row in scores]}


def _prediction_rows(predictions: dict, meta: dict, values: dict, freeze: dict) -> dict:
    wanted = {rid for rid, r in meta["records"].items() if r["group_id"] in _partition(meta["groups"], {"final", "sham"}) and rid in values}
    result = {name: {} for name in MODELS}
    for row in predictions["rows"]:
        name, rid = row["model_id"], row["record_id"]
        _need(name in result and rid in wanted and rid not in result[name], "duplicate, unknown or nonfinal prediction")
        for key in ("value", "lower", "upper", "prediction_sd", "measurement_sd", "last_input_time_s"):
            _number(row[key])
        error = freeze["error_model"]
        _close(row["prediction_sd"], error["prediction_sd"], "frozen prediction uncertainty")
        _close(row["measurement_sd"], error["measurement_sd"], "frozen measurement uncertainty")
        _need(row["measurement_sd"] + 1e-10 >= values[rid]["sd"], "prediction error budget understates independently measured uncertainty")
        width = NormalDist().inv_cdf(0.95) * math.hypot(row["prediction_sd"], row["measurement_sd"])
        _close(row["lower"], row["value"] - width, "statistical prediction interval lower")
        _close(row["upper"], row["value"] + width, "statistical prediction interval upper")
        result[name][rid] = row
    _need(all(set(rows) == wanted for rows in result.values()), "complete predictions for every frozen comparator and biological experiment are required")
    return result


def _quality(rows: dict, meta: dict, values: dict, freeze: dict, protocol: dict) -> dict:
    failures = []

    def require(condition: bool, reason: str) -> None:
        if not condition:
            failures.append(reason)

    final = _partition(meta["groups"], {"final"})
    records, scale = meta["records"], freeze["normalization"]["scale"]
    scoring = protocol["scoring"]
    losses = {name: {gid: mean(abs(row["value"] - values[rid]["value"]) / scale for rid, row in predictions.items()
                              if records[rid]["group_id"] == gid) for gid in final} for name, predictions in rows.items()}
    candidate = losses["candidate"]
    require(mean(candidate.values()) <= scoring["max_mean_group_scaled_mae"], "candidate final mean error exceeds the preregistered observable-specific budget")
    require(max(candidate.values()) <= scoring["max_any_group_scaled_mae"], "a final experiment exceeds the error budget; do not average it away")
    comparisons = {}
    for name, other in losses.items():
        if name == "candidate":
            continue
        wins = sum(candidate[gid] < other[gid] - max(scoring["minimum_absolute_scaled_improvement"], scoring["minimum_relative_improvement"] * other[gid]) for gid in final)
        test = _sign(wins, len(final))
        test["mean_group_scaled_mae"] = mean(other.values())
        comparisons[name] = test
        require(mean(candidate.values()) <= (1 - scoring["minimum_relative_improvement"]) * mean(other.values())
              and mean(other.values()) - mean(candidate.values()) >= scoring["minimum_absolute_scaled_improvement"]
              and test["p_value"] <= protocol["multiple_testing"]["candidate_alpha"],
              f"candidate does not independently beat {name} by the prespecified margin and exact group sign test")
    intervals = scoring["prediction_intervals"]
    coverage, widths = {}, {}
    for gid in final:
        selected = [(rid, row) for rid, row in rows["candidate"].items() if records[rid]["group_id"] == gid]
        coverage[gid] = mean(row["lower"] <= values[rid]["value"] <= row["upper"] for rid, row in selected)
        widths[gid] = mean((row["upper"] - row["lower"]) / scale for _, row in selected)
        require(all((row["upper"] - row["lower"]) / scale <= intervals["maximum_any_scaled_width"] for _, row in selected), "prediction interval width is uninformative")
    covered_groups = sum(value >= intervals["minimum_equal_group_coverage"] for value in coverage.values())
    coverage_interval = _wilson(covered_groups, len(final))
    require(mean(coverage.values()) >= intervals["minimum_equal_group_coverage"]
          and coverage_interval[0] >= intervals["minimum_wilson_lower_coverage"]
          and mean(widths.values()) <= intervals["maximum_mean_scaled_width"], "heldout prediction interval calibration/precision gate failed")
    sham_losses = {}
    for gid in _partition(meta["groups"], {"sham"}):
        selected = [(rid, row) for rid, row in rows["candidate"].items() if records[rid]["group_id"] == gid]
        sham_losses[gid] = mean(abs(row["value"] - values[rid]["value"]) / scale for rid, row in selected)
        require(sham_losses[gid] <= 0.5, "independent sham experiment error exceeds budget")
        for rid, prediction in selected:
            raw = values[rid]["raw"]
            baseline = _number(raw["baseline_raw_value"]) * freeze["observation_adapter"]["factor"]
            baseline_sd = _number(raw["baseline_raw_sd"]) * abs(freeze["observation_adapter"]["factor"])
            require(baseline_sd > 0 and abs(prediction["value"] - baseline) <= 1.96 * math.hypot(baseline_sd, values[rid]["sd"]),
                    "predicted sham response exceeds the independently measured no-change error budget")
    metrics = {"mean_group_scaled_mae": mean(candidate.values()), "group_losses": losses, "comparisons": comparisons,
               "equal_group_coverage": mean(coverage.values()), "coverage_group_wilson_95": list(coverage_interval),
               "mean_scaled_interval_width": mean(widths.values()), "sham_group_losses": sham_losses}
    if failures:
        raise EvidenceError("; ".join(failures), metrics)
    return metrics


def _support(support: dict, predictions: dict, meta: dict, freeze: dict, domain: dict, protocol: dict, reader: _Reader) -> dict:
    _need(support["kind"] == domain["support_kind"], "independent support is not the named activity/observation assay")
    measured = _index(support["rows"], "record_id")
    predicted = _index(predictions["support_predictions"], "record_id")
    _need(set(measured) == set(predicted) == set(meta["support_records"]), "complete independently reserved support outcomes and presealed predictions are required")
    results, wins, failures = [], 0, []
    factor = freeze["support_adapter"]["factor"]
    for rid, row in measured.items():
        record = meta["support_records"][rid]
        original, unit = _source_fields(reader, meta["sources"][record["source_id"]], rid)
        _need(unit == record["unit"], "support source/projection unit mismatch")
        for key in ("low", "high", "low_sd", "high_sd", "covariance"):
            _source_equal(row[key], original[key], "support")
        low, high = _number(row["low"]), _number(row["high"])
        low_sd, high_sd, covariance = _number(row["low_sd"]), _number(row["high_sd"]), _number(row["covariance"])
        _need(low_sd > 0 and high_sd > 0 and abs(covariance) <= low_sd * high_sd, "support paired-measurement uncertainty is invalid")
        variance = (low_sd * low_sd + high_sd * high_sd - 2 * covariance) * factor * factor
        _need(variance > 0, "support contrast requires positive independent measurement uncertainty")
        contrast, sd = (high - low) * factor, math.sqrt(variance)
        point, prediction_sd = _number(predicted[rid]["value"]), _number(predicted[rid]["prediction_sd"])
        _close(prediction_sd, freeze["error_model"]["support_prediction_sd"], "frozen support prediction uncertainty")
        budget = 1.96 * math.hypot(sd, prediction_sd)
        if point == 0 or abs(point - contrast) > budget:
            failures.append(f"independent support assay {rid} contradicts or cannot resolve the frozen response within the declared error budget")
        agrees = point != 0 and math.copysign(1.0, point) * contrast > 1.96 * sd
        wins += agrees
        results.append({"record_id": rid, "measured_contrast": contrast, "measurement_sd": sd,
                        "predicted_contrast": point, "prediction_sd": prediction_sd, "error_budget": budget,
                        "direction_resolved": agrees})
    test = _sign(wins, len(results))
    if len(results) < protocol["selection"]["minimum_support_independent_groups"] or test["p_value"] > protocol["multiple_testing"]["candidate_alpha"]:
        failures.append("independent support lacks resolved directional replication")
    metrics = {"kind": support["kind"], "contrasts": results, "direction_test": test, "scope": "observable_response_not_absolute_latent_activity"}
    if failures:
        raise EvidenceError("; ".join(failures), metrics)
    return metrics


def _same_predictions(replay: dict, expected: dict, name: str) -> None:
    _need(set(replay) == set(expected), f"{name} replay has an incomplete prediction inventory")
    for rid, value in replay.items():
        _close(value, expected[rid]["value"], name)


def _controls(controls: dict, rows: dict, meta: dict, values: dict, freeze: dict, data: dict, domain: dict, protocol: dict) -> dict:
    registered = protocol["controls"]
    expected = rows["candidate"]
    groups = meta["groups"]
    final = sorted(_partition(groups, {"final"}), key=lambda g: _group_order(data["candidate_id"], g))
    _need(controls["input_reassignment"] == {g: final[(i + 1) % len(final)] for i, g in enumerate(final)},
          "input negative control must use the prespecified whole-group derangement")
    for name in ("absolute_clock_shift", "normalization_poison", "measurement_scale_gauge", "unit_conversion"):
        _same_predictions(controls[name]["predictions"], expected, name)
    clock = controls["absolute_clock_shift"]
    _close(clock["seconds"], registered["absolute_clock_shift_seconds"], "absolute_clock_shift")
    clocks = _index(clock["clocks"], "group_id")
    _need(set(clocks) == _partition(groups, {"final", "sham"}), "absolute_clock_shift replay omitted biological groups")
    for gid, row in clocks.items():
        _close(row["before"], groups[gid]["acquisition_clock_s"], "absolute_clock_shift before")
        _close(row["after"], row["before"] + registered["absolute_clock_shift_seconds"], "absolute_clock_shift after")
        _need(row["inputs_before"] == row["inputs_after"] == groups[gid]["inputs"], "absolute_clock_shift changed the actual input response")
    template = controls["source_template_recovery"]
    model = next(model for model in freeze["models"] if model["model_id"] == "source_template")
    _need(template["checkpoint_sha256"] == model["artifact"]["sha256"] and template["label_origin"] == "engineering_source_template",
          "source-template recovery must remain separately labelled engineering evidence")
    _same_predictions(template["teacher"], rows["source_template"], "source_template_recovery teacher")
    _same_predictions(template["recovered"], rows["source_template"], "source_template_recovery recovered")
    poison = controls["normalization_poison"]
    _close(poison["factor"], registered["normalization_poison_factor"], "normalization_poison")
    _need(poison["artifacts_before"] == poison["artifacts_after"] == freeze["dependencies"],
          "normalization_poison changed a fitted coefficient, selected form or normalization artifact")
    _need(set(poison["poisoned_values"]) == set(expected), "normalization_poison lacks actual altered reserved inputs")
    for rid, value in poison["poisoned_values"].items():
        _close(value, values[rid]["raw"]["raw_value"] * registered["normalization_poison_factor"], "normalization_poison altered value")
    for name, factor_name in (("measurement_scale_gauge", "measurement_scale_factor"), ("unit_conversion", "unit_conversion_factor")):
        replay = controls[name]
        factor = registered[factor_name]
        _close(replay["factor"], factor, name)
        converted = _index(replay["rows"], "record_id")
        _need(set(converted) == set(expected), f"{name} measurement replay is incomplete")
        for rid, row in converted.items():
            original = values[rid]["raw"]
            adapter_factor = freeze["observation_adapter"]["factor"]
            if name == "measurement_scale_gauge" and freeze["observation_adapter"]["kind"] == "ratio":
                _close(row["numerator"], original["numerator"] * factor, name)
                _close(row["denominator"], original["denominator"] * factor, name)
                _need(row["denominator"] > 0, "measurement_scale_gauge denominator is invalid")
                _close(row["raw_value"], row["numerator"] / row["denominator"], name)
                _close(row["adapter_factor"], adapter_factor, name)
            else:
                _close(row["raw_value"], original["raw_value"] * factor, name)
                _close(row["adapter_factor"], adapter_factor / factor, name)
            _close(row["raw_value"] * row["adapter_factor"], values[rid]["value"], name)
            _close(row["raw_sd"] * row["adapter_factor"], values[rid]["sd"], f"{name} uncertainty")
    physical = data["physical"]
    _need(physical["source_id"] in meta["sources"] and _number(physical["lower"]) == domain["physical_lower_bound"]
          and physical["upper"] == domain["physical_upper_bound"], "physical constraints lack the correct source-backed observable scope")
    evaluations = controls["physical_consistency"]["rows"]
    indexed = {}
    for row in evaluations:
        key = (row["model_id"], row["record_id"])
        _need(key not in indexed, "duplicate physical consistency witness")
        indexed[key] = row
    _need(set(indexed) == {(model, rid) for model, predictions in rows.items() for rid in predictions},
          "physical consistency must cover all registered model/experiment predictions")
    maximum_residual = 0.0
    for (model, rid), witness in indexed.items():
        prediction = rows[model][rid]
        _close(witness["value"], prediction["value"], "physical prediction binding")
        _need(prediction["value"] >= physical["lower"] and (physical["upper"] is None or prediction["value"] <= physical["upper"]),
              "native prediction violates physical bounds")
        _close(witness["last_input_time_s"], prediction["last_input_time_s"], "causal input binding")
        _need(witness["last_input_time_s"] <= meta["records"][rid]["time_s"], "response prediction used a future input")
        if domain["canonical_unit"] == "dimensionless_ratio":
            _need(physical["kind"] == "localization_ratio", "localization cannot inherit phosphorylation or flux physical constraints")
        else:
            _need(physical["kind"] == "stoichiometric_flux", "flux requires actual stoichiometric and capacity witnesses")
            vector = witness["fluxes"]
            _need(set(vector) == set(physical["bounds"]), "physical flux witness lacks reaction coverage")
            for reaction, value in vector.items():
                lower, upper = map(_number, physical["bounds"][reaction])
                _need(lower <= _number(value) <= upper, "flux witness violates direction/capacity constraints")
            _close(vector[physical["target_reaction"]], prediction["value"], "native flux reaction binding")
            _need(bool(physical["balances"]), "source-linked carbon/charge balances are absent")
            for balance in physical["balances"]:
                _need(set(balance["coefficients"]) <= set(vector) and _number(balance["scale"]) > 0,
                      "invalid physical conservation constraint")
                residual = abs(sum(_number(coefficient) * _number(vector[reaction]) for reaction, coefficient in balance["coefficients"].items())
                               - _number(balance["rhs"])) / balance["scale"]
                maximum_residual = max(maximum_residual, residual)
                _need(residual <= 1e-8, "native flux violates carbon/charge/stoichiometric consistency")
    return {"replays": registered["required_replays"], "physical_witnesses": len(indexed), "maximum_normalized_balance_residual": maximum_residual,
            "gauge_scope": "measured_observable_only", "source_template_labels_are_biological_evidence": False}


def assess_native_law(bundle: EvidenceBundle, *, root: str | Path, trust: ValidationTrust) -> Assessment:
    gates = []
    candidate = None
    scope = "No biological claim: evidence has not been assessed."
    refs = {}
    reader = None

    def finish() -> Assessment:
        passed = bool(gates) and all(gate.state == "pass" for gate in gates)
        real = reader is not None and set(reader.issuers) == set(RECEIPT_KINDS)
        if real:
            real = all(reader.authorities[issuer].domain == "real" for issuer in reader.issuers.values())
        authorized = passed and real
        status = ("authorized_scoped_empirical_response" if authorized else "fixture_only_not_biology" if passed
                  else "refused" if any(g.state == "fail" for g in gates) else "missing_evidence")
        return Assessment(candidate, GRADE, authorized, status, scope, tuple(gates), {key: ref.sha256 for key, ref in refs.items()})

    def check(name: str, action: Callable[[], object]) -> object | None:
        try:
            result = action()
            gates.append(Gate(name, "pass", metrics=result if isinstance(result, dict) else {}))
            return result
        except (KeyError, FileNotFoundError, MissingEvidence) as exc:
            gates.append(Gate(name, "missing", (f"missing required evidence for {name}: {exc}",)))
        except EvidenceError as exc:
            gates.append(Gate(name, "fail", (f"{name}: {exc}",), exc.metrics))
        except (ValueError, TypeError, AttributeError, IndexError, ArithmeticError, OSError) as exc:
            gates.append(Gate(name, "fail", (f"{name}: {exc}",)))
        return None

    if not isinstance(bundle, EvidenceBundle):
        gates.append(Gate("evidence_bundle", "fail", ("typed sealed EvidenceBundle required; dictionaries and pass flags are not authorization",)))
        return finish()
    missing = [kind for kind in RECEIPT_KINDS if getattr(bundle, kind) is None]
    if missing:
        gates.extend(Gate(kind, "missing", (f"missing authoritative sealed {kind} artifact",)) for kind in missing)
        return finish()
    invalid = [kind for kind in RECEIPT_KINDS if not isinstance(getattr(bundle, kind), SealedArtifact)]
    if invalid:
        gates.append(Gate("evidence_references", "fail", (f"typed sealed artifact reference required for: {', '.join(invalid)}",)))
        return finish()
    refs = {kind: getattr(bundle, kind) for kind in RECEIPT_KINDS}
    reader = check("authority_configuration", lambda: _Reader(Path(root), trust))
    if reader is None:
        return finish()
    registry = check("registry_authentication", lambda: reader.receipt(refs["registry"], "registry"))
    if registry is None:
        return finish()

    def policy() -> dict:
        nonlocal candidate, scope
        candidate = _text(registry["candidate_id"])
        head = trust.registry_heads.get(f"{PROTOCOL_ID}/{candidate}")
        _need(head == refs["registry"].sha256, "receipt is not the current external registry head; caller-recomputed hashes cannot establish independence")
        protocol_raw, domain_raw = reader.raw(PROTOCOL_PATH), reader.raw(DOMAIN_PATH)
        protocol, domains = _parse(protocol_raw), _parse(domain_raw)
        _need(protocol["protocol_id"] == domains["protocol_id"] == PROTOCOL_ID and candidate in domains["candidates"],
              "unsupported claim/protocol; no universal, causal, intrinsic enzyme or product claim is implemented")
        _need(protocol["multiple_testing"]["candidate_alpha"] == 0.025 and protocol["multiple_testing"]["maximum_candidates"] == 2
              and protocol["multiple_testing"]["maximum_final_evaluations_per_candidate"] == 1,
              "preregistered final-test budget/alpha must not drift")
        domain = domains["candidates"][candidate]
        scope = domain["scope"]
        return {"protocol": protocol, "domain": domain, "protocol_sha256": hashlib.sha256(protocol_raw).hexdigest(),
                "claim_domain_sha256": hashlib.sha256(domain_raw).hexdigest()}

    bound = check("registered_policy", policy)
    if bound is None:
        return finish()
    documents = {"registry": registry}
    for kind in ("dataset", "lineage", "freeze"):
        document = check(f"{kind}_authentication", lambda kind=kind: reader.receipt(refs[kind], kind))
        if document is None:
            return finish()
        documents[kind] = document
    data, lineage, freeze = (documents[key] for key in ("dataset", "lineage", "freeze"))
    meta = check("data_provenance_and_independence", lambda: _data_metadata(data, bound["domain"], bound["protocol"]))
    if meta is None:
        return finish()
    check("registration_and_single_final_release", lambda: _registry(registry, meta, freeze, refs, bound["protocol_sha256"], bound["claim_domain_sha256"]))
    check("ancestral_fit_selection_lineage", lambda: _lineage(lineage, freeze, meta, reader))
    check("selection_observation_and_scope_freeze", lambda: _freeze(freeze, data, meta, bound["domain"], reader))
    if any(g.state != "pass" for g in gates):
        return finish()
    for kind in ("observations", "predictions", "support", "controls"):
        document = check(f"{kind}_authentication", lambda kind=kind: reader.receipt(refs[kind], kind))
        if document is None:
            return finish()
        documents[kind] = document

    def roles() -> dict:
        for kind, document in documents.items():
            _need(document["candidate_id"] == candidate, f"{kind} belongs to a different candidate")
        _need(len({reader.issuers[key] for key in ("registry", "dataset", "lineage", "freeze")}) == 4,
              "independent custody, ancestry, execution and registry authorities must be separate")
        _need(reader.issuers["observations"] == reader.issuers["dataset"] == reader.issuers["support"],
              "measurement/support receipts must be issued by the source custodian")
        _need(reader.issuers["predictions"] == reader.issuers["controls"] == reader.issuers["freeze"],
              "predictions and control replays must come from the isolated frozen executor")
        _need(all(documents[kind]["run_id"] == registry["run_id"] for kind in ("freeze", "predictions", "controls")),
              "multiple or substituted final model runs are not allowed")
        return {"issuers": dict(reader.issuers)}

    check("authority_role_and_run_binding", roles)
    if any(g.state != "pass" for g in gates):
        return finish()
    values = check("source_bound_measurements", lambda: _measurements(documents["observations"], data, meta, freeze, reader))
    if values is None:
        return finish()
    check("recomputed_development_selection_and_normalization", lambda: _selection(freeze, meta, values, reader))
    rows = check("complete_presealed_predictions", lambda: _prediction_rows(documents["predictions"], meta, values, freeze))
    if rows is not None:
        check("prediction_quality", lambda: _quality(rows, meta, values, freeze, bound["protocol"]))
        check("negative_controls_gauge_and_physics", lambda: _controls(documents["controls"], rows, meta, values, freeze, data, bound["domain"], bound["protocol"]))
    check("independent_matched_support", lambda: _support(documents["support"], documents["predictions"], meta, freeze, bound["domain"], bound["protocol"], reader))
    return finish()
