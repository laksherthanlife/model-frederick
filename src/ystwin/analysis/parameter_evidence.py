from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

from ystwin.mech.params import Param, ParameterUncertainty

PROGRAMMES = (
    "carbon_pka_snf1", "nitrogen_tor", "oxidative_yap1", "heat_hsf1_hsp70",
    "osmotic_hog", "upr_hac1", "ph",
)
EVIDENCE_STATUSES = ("measured", "fitted", "prior", "refused")
PARAMETER_SCALES = ("biological", "assay", "physical", "input")
DEFAULT_INVENTORY = "data/parameter_evidence.json"
JALIHAL_SBML_SHA256 = "a5416c7785df8001044602f2c47cfa084109396a0ee165b36e844d35b00fb1ed"
JALIHAL_EXPORT_ADAPTER = "jalihal2021_ode_export_v1"
JALIHAL_NORMALIZED_INPUTS = ("Carbon", "ATP", "Glutamine_ext", "NH4", "Proline")


class EvidenceGap(NotImplementedError):
    pass


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def content_sha256(value) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _text(value, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} requires nonempty text")
    return value


def _digest(value, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{label} requires a lowercase SHA-256")
    return value


def _relative_path(value: str) -> str:
    _text(value, "source path")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value or "\\" in value:
        raise ValueError("source paths must be canonical repository-relative paths")
    return value


def source_path(root: str | Path, relative: str) -> Path:
    relative = _relative_path(relative)
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("source path escapes the repository through a symlink")
    return path


def _unique_strings(values, label: str, *, nonempty: bool = True) -> tuple[str, ...]:
    if not isinstance(values, (list, tuple)) or (nonempty and not values):
        raise ValueError(f"{label} requires a {'nonempty ' if nonempty else ''}list")
    result = tuple(_text(value, label) for value in values)
    if len(set(result)) != len(result):
        raise ValueError(f"duplicate {label}")
    return result


def _no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


@dataclass(frozen=True)
class EvidenceParameter:
    parameter: Param
    status: str
    scale: str
    source_id: str
    locator: str
    units_status: str
    context: str
    identifiability: str
    unresolved_directions: tuple[str, ...]
    uncertainty: ParameterUncertainty
    interpretation: str

    def to_dict(self) -> dict:
        return {
            "id": self.parameter.name, "value": self.parameter.value,
            "units": self.parameter.units, "legacy_tag": self.parameter.tag,
            "status": self.status, "scale": self.scale, "source_id": self.source_id,
            "locator": self.locator, "units_status": self.units_status,
            "context": self.context, "identifiability": self.identifiability,
            "unresolved_directions": list(self.unresolved_directions),
            "uncertainty": self.uncertainty.to_dict(),
            "interpretation": self.interpretation,
            "missing": self.parameter.missing,
        }


def _parameter(record: dict, sources: dict) -> EvidenceParameter:
    required = {
        "id", "programme", "source_id", "locator", "value", "units", "units_status",
        "status", "scale", "context", "uncertainty", "identifiability",
        "unresolved_directions", "interpretation", "missing",
    }
    if not isinstance(record, dict) or set(record) != required:
        raise ValueError("parameter evidence fields do not match schema version 1")
    for key in ("id", "source_id", "locator", "units", "context", "interpretation"):
        _text(record[key], f"parameter {key}")
    if record["programme"] not in (*PROGRAMMES, "native_support", "observation"):
        raise ValueError("unknown parameter programme")
    if record["status"] not in EVIDENCE_STATUSES or record["scale"] not in PARAMETER_SCALES:
        raise ValueError("unknown evidence status or parameter scale")
    if record["source_id"] not in sources:
        raise ValueError("unknown parameter source")
    if record["units_status"] not in ("source_declared", "unspecified_in_source", "required_not_available"):
        raise ValueError("unknown units status")
    if record["identifiability"] not in ("source_context_only", "not_assessed", "unresolved", "composite_only"):
        raise ValueError("unknown identifiability status")
    directions = _unique_strings(record["unresolved_directions"], "unresolved directions", nonempty=False)
    if record["identifiability"] in ("unresolved", "composite_only") and not directions:
        raise ValueError("unresolved identifiability must retain its unresolved directions")
    uncertainty = ParameterUncertainty(**record["uncertainty"])
    if uncertainty.units != record["units"]:
        raise ValueError("parameter and uncertainty units differ")
    value = record["value"]
    if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)):
        raise ValueError("parameter value must be a finite real number or null")
    if uncertainty.interval is not None and value is not None:
        if not uncertainty.interval[0] <= value <= uncertainty.interval[1]:
            raise ValueError("parameter value lies outside its uncertainty interval")
    source = sources[record["source_id"]]
    provenance = f"{record['source_id']}: {record['locator']}; {record['interpretation']}"
    if record["status"] == "refused":
        if value is not None or uncertainty.kind != "unresolved":
            raise ValueError("refused evidence cannot carry a value or resolved uncertainty")
        parameter = Param.refused(
            record["id"], record["units"], provenance,
            reason=record["interpretation"], missing=_text(record["missing"], "missing evidence"),
        )
    elif record["status"] == "measured":
        if source["verification"]["level"] not in ("primary_table", "local_measurement"):
            raise ValueError("measured values require checked primary-table or local-measurement evidence")
        if record["units_status"] != "source_declared":
            raise ValueError("a measured value requires verified source-declared units")
        checked = source["verification"].get("checked_measurements", {}).get(record["id"])
        if checked != {"value": value, "units": record["units"], "locator": record["locator"]}:
            raise ValueError("measured parameter does not match a checked source measurement")
        parameter = Param.measured(record["id"], value, record["units"], provenance)
    else:
        parameter = Param.asserted(
            record["id"], value, record["units"],
            f"{record['status'].upper()}, not an independently measured kinetic constant. {provenance}",
            missing=record["missing"],
        )
    return EvidenceParameter(
        parameter, record["status"], record["scale"], record["source_id"], record["locator"],
        record["units_status"], record["context"], record["identifiability"], directions,
        uncertainty, record["interpretation"],
    )


@dataclass(frozen=True)
class EvidenceInventory:
    _json: str

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self._json.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict:
        return json.loads(self._json)

    @property
    def parameters(self) -> tuple[EvidenceParameter, ...]:
        data = self.to_dict()
        return tuple(_parameter(record, data["sources"]) for record in data["parameters"])

    def programme(self, name: str) -> dict:
        if name not in PROGRAMMES:
            raise ValueError(f"unknown signalling programme {name!r}")
        return self.to_dict()["programmes"][name]

    def require_reference(self, programme: str) -> tuple[str, ...]:
        record = self.programme(programme)
        models = self.to_dict()["models"]
        available = tuple(key for key in record["model_ids"] if models[key]["availability"] == "local_verified")
        if not available:
            raise EvidenceGap(f"{programme}: " + "; ".join(record["blockers"]))
        return available


def validate_inventory(data: dict) -> EvidenceInventory:
    if not isinstance(data, dict) or set(data) != {
        "schema_version", "inventory_version", "scope", "firewall", "programmes",
        "sources", "models", "parameters", "native_data", "shared_assumptions",
    }:
        raise ValueError("unexpected evidence inventory schema")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise ValueError("unsupported evidence inventory schema version")
    _text(data["inventory_version"], "inventory version")
    _text(data["scope"], "inventory scope")
    if set(data["programmes"]) != set(PROGRAMMES):
        raise ValueError("inventory must cover all seven signalling programmes")
    firewall = data["firewall"]
    if firewall != {
        "domain": "native_only", "product_outcomes_allowed": False,
        "reserved_response_values_allowed": False,
        "freeze_scope": "prospective_local_content_binding_not_external_preregistration",
    }:
        raise ValueError("native/product-blind evidence firewall is mandatory")
    sources, models = data["sources"], data["models"]
    if not isinstance(sources, dict) or not sources or not isinstance(models, dict):
        raise ValueError("sources and models must be keyed records")
    for key, source in sources.items():
        _text(key, "source id")
        if source["domain"] != "native" or source["use"] not in ("reference", "evaluation", "inventory_only"):
            raise ValueError("source is outside the native allowlist")
        _text(source["title"], "source title")
        _text(source["verification"]["level"], "verification level")
        _text(source["verification"]["finding"], "verification finding")
        if not isinstance(source["identity"], dict) or not source["identity"]:
            raise ValueError("sources require a checked identity")
        seen_paths = set()
        for asset in source["local_artifacts"]:
            _relative_path(asset["path"])
            _digest(asset["sha256"], "source artifact")
            if asset["path"] in seen_paths:
                raise ValueError("duplicate source artifact path")
            seen_paths.add(asset["path"])
            if "bytes" in asset and (type(asset["bytes"]) is not int or asset["bytes"] <= 0):
                raise ValueError("source artifact byte counts must be positive integers")
    for key, model in models.items():
        _text(key, "model id")
        if model["source_id"] not in sources:
            raise ValueError("model has an unknown source")
        if model["availability"] not in (
            "local_verified", "local_source_only", "upstream_verified_not_vendored", "equations_only", "identity_only",
        ):
            raise ValueError("unknown model availability")
        if "execution_adapter" in model and (
            key != "jalihal2021" or model["execution_adapter"] != JALIHAL_EXPORT_ADAPTER
            or model.get("sha256") != JALIHAL_SBML_SHA256
        ):
            raise ValueError("execution adapter is restricted to the audited immutable Jalihal export")
        if model.get("execution_adapter") == JALIHAL_EXPORT_ADAPTER:
            audit = model.get("unit_audit", {})
            time = audit.get("time", {})
            if (time.get("status") != "conflicting_declarations" or time.get("sbml_declared") != "second"
                    or time.get("upstream_timecourse") != "minute" or time.get("numeric_conversion") != "none"
                    or "physical_seconds_per_native_unit" not in time or time["physical_seconds_per_native_unit"] is not None
                    or audit.get("states", {}).get("absolute_units_verified") is not False
                    or audit.get("normalized_input_ids") != list(JALIHAL_NORMALIZED_INPUTS)
                    or "physical_input_calibration" not in audit or audit["physical_input_calibration"] is not None):
                raise ValueError("the native-unit audit must retain unresolved time, state, and nutrient calibration")
        if model["availability"] == "local_source_only":
            _relative_path(model["archive_path"])
            _text(model["runtime"], "source-only runtime")
            declared = {asset["path"] for asset in sources[model["source_id"]]["local_artifacts"]}
            if model["archive_path"] not in declared:
                raise ValueError("a source-only archive must be explicitly pinned by its source")
        if not set(_unique_strings(model["programmes"], "model programmes")) <= set(PROGRAMMES):
            raise ValueError("unknown model programme")
        if model["availability"] == "local_verified":
            _relative_path(model["path"])
            _digest(model["sha256"], "model")
            if not isinstance(model["model_id"], str):
                raise ValueError("internal SBML identity must be explicit text, including empty if the source omits it")
            if sources[model["source_id"]]["use"] != "reference":
                raise ValueError("local reference model lacks source permission")
    for source_id in sources:
        _source_assets(data, source_id)
    for name, programme in data["programmes"].items():
        _text(programme["scope"], "programme scope")
        _unique_strings(programme["blockers"], "programme blockers")
        for key in _unique_strings(programme["model_ids"], "programme model ids"):
            if key not in models or name not in models[key]["programmes"]:
                raise ValueError("programme/model identity mismatch")
    parameters = [_parameter(record, sources) for record in data["parameters"]]
    if len({p.parameter.name for p in parameters}) != len(parameters):
        raise ValueError("duplicate parameter evidence id")
    if not set(PROGRAMMES) <= {p["programme"] for p in data["parameters"]}:
        raise ValueError("each signalling programme needs parameter evidence or a refusal")
    for record in data["native_data"]:
        if record["source_id"] not in sources:
            raise ValueError("native data has an unknown source")
        _text(record["units"], "native data units")
        _text(record["gap"], "native data gap")
        declared_paths = {asset["path"] for asset in sources[record["source_id"]]["local_artifacts"]}
        paths = _unique_strings(record["artifact_paths"], "native asset paths", nonempty=False)
        if not set(paths) <= declared_paths:
            raise ValueError("native data paths must be explicitly pinned by their registered source")
    _unique_strings(data["shared_assumptions"], "shared assumptions")
    return EvidenceInventory(canonical_json(data))


def load_parameter_evidence(path: str | Path | None = None) -> EvidenceInventory:
    path = repository_root() / DEFAULT_INVENTORY if path is None else Path(path)
    data = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_no_duplicates)
    return validate_inventory(data)


def _source_assets(data: dict, source_id: str) -> list[dict]:
    assets = {asset["path"]: dict(asset) for asset in data["sources"][source_id]["local_artifacts"]}
    for key, model in data["models"].items():
        if model["source_id"] != source_id or model["availability"] != "local_verified":
            continue
        asset = assets.setdefault(model["path"], {"path": model["path"], "sha256": model["sha256"]})
        if asset["sha256"] != model["sha256"]:
            raise ValueError("source artifact and model SHA-256 declarations disagree")
        if "bytes" in asset and "bytes" in model and asset["bytes"] != model["bytes"]:
            raise ValueError("source artifact and model byte count declarations disagree")
        asset["model_key"] = key
    return list(assets.values())


def _asset_check(asset: dict, payload: bytes | None) -> dict:
    actual = None if payload is None else hashlib.sha256(payload).hexdigest()
    size = None if payload is None else len(payload)
    expected_blob = asset.get("acquisition", {}).get("git_blob")
    blob = None
    if expected_blob is not None and payload is not None:
        blob = hashlib.sha1(f"blob {size}\0".encode() + payload, usedforsecurity=False).hexdigest()
    return {
        "actual_sha256": actual, "actual_bytes": size, "actual_git_blob": blob,
        "verified": actual == asset["sha256"] and ("bytes" not in asset or size == asset["bytes"])
        and (expected_blob is None or blob == expected_blob),
    }


def load_source_asset(source_id: str, relative: str, inventory: EvidenceInventory | None = None, *, root=None) -> bytes:
    inventory = load_parameter_evidence() if inventory is None else inventory
    root = repository_root() if root is None else root
    data = inventory.to_dict()
    if source_id not in data["sources"]:
        raise ValueError(f"unknown source {source_id!r}")
    assets = {asset["path"]: asset for asset in _source_assets(data, source_id)}
    if relative not in assets:
        raise ValueError("asset is not explicitly registered to this native source")
    payload = source_path(root, relative).read_bytes()
    if not _asset_check(assets[relative], payload)["verified"]:
        raise ValueError(f"source SHA-256, byte count, or Git blob mismatch: {relative}")
    return payload


def audit_local_sources(inventory: EvidenceInventory, root: str | Path | None = None) -> list[dict]:
    root = repository_root() if root is None else Path(root)
    data = inventory.to_dict()
    results = []
    for source_id in data["sources"]:
        for asset in _source_assets(data, source_id):
            path = source_path(root, asset["path"])
            payload = path.read_bytes() if path.is_file() else None
            results.append({"source_id": source_id, **asset, **_asset_check(asset, payload)})
    return results


def _jalihal_export_model(path: Path, payload: bytes, record: dict):
    import libsbml

    from ystwin.mech.kinetic_sbml import KineticModel

    if hashlib.sha256(payload).hexdigest() != JALIHAL_SBML_SHA256:
        raise ValueError("the Jalihal export adapter cannot repair an unrecognized source")
    document = libsbml.readSBMLFromString(payload.decode("utf-8"))
    source = document.getModel()
    formulas = {reaction.getId(): libsbml.formulaToL3String(reaction.getKineticLaw().getMath())
                for reaction in source.getListOfReactions()}
    min_nodes = 0

    def lower_min(node):
        nonlocal min_nodes
        children = [lower_min(node.getChild(index)) for index in range(node.getNumChildren())]
        if node.getType() == libsbml.AST_FUNCTION_MIN:
            if len(children) != 2:
                raise ValueError("only the audited binary min expressions may be lowered")
            min_nodes += 1
            left, right = children
            condition = libsbml.ASTNode(libsbml.AST_RELATIONAL_LEQ)
            condition.addChild(left.deepCopy())
            condition.addChild(right.deepCopy())
            result = libsbml.ASTNode(libsbml.AST_FUNCTION_PIECEWISE)
            result.addChild(left)
            result.addChild(condition)
            result.addChild(right)
            return result
        result = node.deepCopy()
        for index, child in enumerate(children):
            result.replaceChild(index, child)
        return result

    # This is NOT an SBML Level 3 default. The exporter leaves it undefined.
    # The pinned variables.txt and translator establish one derivative per product.
    supplied_stoichiometry = []
    for species, reaction in zip(source.getListOfSpecies(), source.getListOfReactions(), strict=True):
        if reaction.getNumReactants() or reaction.getNumProducts() != 1:
            raise ValueError("Jalihal source is not the audited one-ODE-per-product export")
        product = reaction.getProduct(0)
        if (product.getSpecies() != species.getId() or product.isSetStoichiometry()
                or not product.getConstant()):
            raise ValueError("unexpected source product stoichiometry; no implicit repair")
        product.setStoichiometry(1.0)
        supplied_stoichiometry.append(reaction.getId())
        reaction.getKineticLaw().setMath(lower_min(reaction.getKineticLaw().getMath()))
    if len(supplied_stoichiometry) != 25 or min_nodes != 5:
        raise ValueError("Jalihal export compatibility inventory changed")
    execution_bytes = libsbml.writeSBMLToString(document).encode("utf-8")
    audit_json = canonical_json({
        "adapter": JALIHAL_EXPORT_ADAPTER, "source_sha256": JALIHAL_SBML_SHA256,
        "execution_sbml_sha256": hashlib.sha256(execution_bytes).hexdigest(),
        "source_bytes_modified": False, "unit_product_coefficients_supplied": supplied_stoichiometry,
        "binary_min_nodes_lowered": min_nodes, "min_identity": "min(a,b) = piecewise(a,a<=b,b)",
        "source_reaction_formulas": formulas, "unit_audit": record["unit_audit"],
        "interpretation": "source-normalized ODE export compatibility, not corrected physical kinetics or biological validation",
    })

    class SourceNativeKineticModel(KineticModel):
        @property
        def metadata(self):
            metadata = super().metadata
            metadata["native_reference"] = json.loads(audit_json)
            metadata["units"]["physical_calibration"] = {"time": None, "states": None, "nutrient_inputs": None}
            metadata["units"]["interpretation"] = (
                "Literal SBML second and mole/litre declarations conflict with upstream minute/normalized-state semantics; "
                "model-native arithmetic only, not physical time, concentration, or nutrient calibration."
            )
            return metadata

    return SourceNativeKineticModel(document, path, JALIHAL_SBML_SHA256)


def load_reference_model(model_key: str, inventory: EvidenceInventory | None = None, *, root=None):
    from ystwin.mech.kinetic_sbml import KineticModel

    inventory = load_parameter_evidence() if inventory is None else inventory
    root = repository_root() if root is None else root
    record = inventory.to_dict()["models"].get(model_key)
    if record is None:
        raise ValueError(f"unknown model key {model_key!r}")
    if record["availability"] != "local_verified":
        raise EvidenceGap(f"{model_key}: {record['gap']}")
    path = source_path(root, record["path"])
    payload = load_source_asset(record["source_id"], record["path"], inventory, root=root)
    if hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError(f"source SHA-256 mismatch for {model_key}")
    if record.get("execution_adapter") == JALIHAL_EXPORT_ADAPTER:
        model = _jalihal_export_model(path, payload, record)
    else:
        model = KineticModel.from_sbml(path)
    metadata = model.metadata
    if model.model_id != record["model_id"] or metadata["sha256"] != record["sha256"]:
        raise ValueError(f"internal SBML identity mismatch for {model_key}")
    for field in ("sbml_level", "sbml_version"):
        if metadata[field] != record[field]:
            raise ValueError(f"SBML version mismatch for {model_key}")
    if (len(model.species_ids), len(model.reaction_ids), len(metadata["parameter_definitions"])) != (
        record["species"], record["reactions"], record["parameters"],
    ):
        raise ValueError(f"SBML inventory count mismatch for {model_key}")
    return model


def model_parameter_inventory(model_key: str, inventory: EvidenceInventory | None = None, *, root=None) -> dict:
    model = load_reference_model(model_key, inventory, root=root)
    metadata = model.metadata
    native = metadata.get("native_reference", {})
    records = []
    for name, value in model.parameter_values.items():
        is_normalized_input = bool(native) and name in JALIHAL_NORMALIZED_INPUTS
        records.append({
            "id": name, "value": value, "units": metadata["units"]["parameters"][name],
            "status": "prior", "source_set_status": (
                "published_normalized_fit_with_fixed_inputs_and_totals" if native else "published_joint_fit_mixed_attribution"
            ),
            "scale": "normalized_input" if is_normalized_input else "model_native_unresolved",
            "locator": f"listOfParameters/parameter[@id='{name}']",
            "uncertainty": ParameterUncertainty(
                "not_reported", metadata["units"]["parameters"][name],
                "The SBML supplies a point, not a per-parameter uncertainty or measured-versus-fitted attribution.",
            ).to_dict(),
            "identifiability": "not_assessed",
            "unresolved_directions": ["joint kinetic parameters and observation scales; no rank certificate"],
        })
    return {
        "model_key": model_key, "model_id": model.model_id, "sha256": metadata["sha256"],
        "parameters": records,
        "assignments": metadata["assignment_formulas"],
        "reaction_formulas": native.get("source_reaction_formulas", metadata["reaction_formulas"]),
        "local_parameters": metadata["local_parameters"], "units": metadata["units"],
        "native_reference": native,
        "gap": "Inline equation constants and mixed-source parameter attribution are retained in the source equations, not certified as measured kinetics.",
    }


@dataclass(frozen=True)
class NativeReferenceTrajectory:
    times_native: np.ndarray
    states: np.ndarray
    variables: dict[str, np.ndarray]
    reaction_rates: np.ndarray
    metadata: dict


def simulate_native_reference(
    model_key: str, times_native, *, normalized_inputs: Mapping[str, float],
    inventory: EvidenceInventory | None = None, root=None, initial_state=None,
    rtol=1e-10, atol=1e-12, method="Radau",
) -> NativeReferenceTrajectory:
    """Replay normalized source inputs/time; the audited stiff solver has no fallback.

    Rich-input late-time oscillations are tolerance-sensitive. Numerical agreement
    on an audited horizon does not calibrate physical time, dose, or extrapolation.
    """
    if model_key != "jalihal2021":
        raise EvidenceGap("this source-native normalized-input interface is only audited for jalihal2021")
    if not isinstance(normalized_inputs, Mapping) or set(normalized_inputs) != set(JALIHAL_NORMALIZED_INPUTS):
        raise ValueError("declare all five normalized source inputs: " + ", ".join(JALIHAL_NORMALIZED_INPUTS))
    inputs = {}
    for name in JALIHAL_NORMALIZED_INPUTS:
        value = normalized_inputs[name]
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"{name} must be an explicit finite normalized source input in [0, 1], not g/L")
        inputs[name] = float(value)
    model = load_reference_model(model_key, inventory, root=root)
    trajectory = model.simulate(times_native, parameters=inputs, initial_state=initial_state,
                                rtol=rtol, atol=atol, method=method)
    return NativeReferenceTrajectory(
        trajectory.times_s.copy(), trajectory.states, trajectory.variables, trajectory.reaction_rates,
        {
            "model_key": model_key, "source_sha256": model.metadata["sha256"],
            "time_axis": "source_native_numeric", "physical_seconds_per_native_unit": None,
            "normalized_inputs": inputs, "physical_input_calibration": None,
            "state_scale": "source-normalized activities and pools, not absolute concentrations or kinase assays",
            "native_reference": model.metadata["native_reference"],
            "sbml_execution": trajectory.metadata,
            "scope": "conditional source-native simulation; not independent biological validation",
            "independent_real_validation_established": False,
        },
    )


@dataclass(frozen=True)
class FrozenEvidence:
    _json: str

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self._json.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict:
        return json.loads(self._json)


def freeze_evidence(
    inventory: EvidenceInventory, *, source_ids, product_outcomes_seen: bool,
    root: str | Path | None = None,
) -> FrozenEvidence:
    if product_outcomes_seen is not False:
        raise ValueError("freeze requires an explicit no-product-outcomes declaration for this run")
    selected = _unique_strings(source_ids, "freeze source ids")
    data = inventory.to_dict()
    if any(key not in data["sources"] or data["sources"][key]["use"] != "reference" for key in selected):
        raise ValueError("freeze accepts only explicitly allowlisted native reference sources")
    root = repository_root() if root is None else root
    assets = []
    for key in selected:
        for asset in _source_assets(data, key):
            payload = source_path(root, asset["path"]).read_bytes()
            check = _asset_check(asset, payload)
            if not check["verified"]:
                raise ValueError(f"source SHA-256, byte count, or Git blob mismatch before freeze: {asset['path']}")
            assets.append({"source_id": key, "path": asset["path"], "sha256": check["actual_sha256"]})
    if not assets:
        raise EvidenceGap("no local verified reference assets to freeze")
    return FrozenEvidence(canonical_json({
        "schema_version": 1, "inventory_sha256": inventory.sha256,
        "source_ids": sorted(selected), "assets": sorted(assets, key=lambda a: (a["source_id"], a["path"])),
        "product_outcomes_seen": False,
        "scope": "local content freeze before this run; not proof of ancestral blindness or external preregistration",
    }))


def verify_evidence_freeze(frozen: FrozenEvidence, inventory: EvidenceInventory, *, root=None) -> None:
    payload = frozen.to_dict()
    current = freeze_evidence(
        inventory, source_ids=payload["source_ids"], product_outcomes_seen=payload["product_outcomes_seen"], root=root,
    )
    if current.sha256 != frozen.sha256:
        raise ValueError("evidence freeze no longer matches the inventory or source bytes")


def native_asset_inventory(inventory: EvidenceInventory, *, root=None) -> list[dict]:
    root = repository_root() if root is None else root
    records = []
    for dataset in inventory.to_dict()["native_data"]:
        for relative in dataset["artifact_paths"]:
            path = source_path(root, relative)
            payload = path.read_bytes() if path.is_file() else None
            records.append({
                "source_id": dataset["source_id"], "path": relative,
                "sha256": None if payload is None else hashlib.sha256(payload).hexdigest(),
                "bytes": None if payload is None else len(payload),
                "scope": "byte identity only; no response values parsed or released",
            })
    return records


def evidence_report(inventory: EvidenceInventory | None = None, *, root=None) -> dict:
    inventory = load_parameter_evidence() if inventory is None else inventory
    data = inventory.to_dict()
    parameters = inventory.parameters
    return {
        "inventory_version": data["inventory_version"], "inventory_sha256": inventory.sha256,
        "scope": data["scope"], "programmes": data["programmes"], "models": data["models"],
        "parameter_status_counts": {status: sum(p.status == status for p in parameters) for status in EVIDENCE_STATUSES},
        "parameter_scale_counts": {scale: sum(p.scale == scale for p in parameters) for scale in PARAMETER_SCALES},
        "source_checks": audit_local_sources(inventory, root),
        "native_data": data["native_data"], "shared_assumptions": data["shared_assumptions"],
        "full_biological_coverage": False, "independent_real_validation_established": False,
    }
