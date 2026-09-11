from __future__ import annotations

import ast
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import libsbml
import numpy as np
import pandas as pd

from .hog_data import load_native_hog_data
from .hog_learning import fit_native_hog_parameters
from .native_physiology import (
    NativeExchangeModel,
    evaluate_native_exchange_model,
    fit_native_exchange_model,
    load_native_chemostat_data,
)
from ..mech.hog import HogModel, HogProtocol, hog_glycerol_balance


CHECKPOINT_PATH = "outputs/native_training_run_02/checkpoint.json"
CHECKPOINT_SHA256 = "384175e19b2451f542ab08e347ee9d9021e22772fbc205e0c404bbc28416ec94"
_NATIVE_INPUTS = (
    "data/hog2013/model_wt.xml", "data/hog2013/observations.xls",
    "data/hog2013/methods.pdf", "data/hog2013/sources.json",
    "data/physiology/chemostatData_VanHoek1998.tsv", "data/gem/ecYeastGEM_batch.xml.gz",
    "data/hog2013/native_training_protocol.json",
)
_AUDITED_CODE = (
    "src/ystwin/analysis/hog_learning.py", "src/ystwin/analysis/hog_data.py",
    "src/ystwin/analysis/native_physiology.py", "src/ystwin/analysis/training_freeze.py",
    "src/ystwin/analysis/blind_transfer.py", "src/ystwin/mech/hog.py",
    "src/ystwin/mech/kinetic_sbml.py", "src/ystwin/fba/native_obligations.py",
    "src/ystwin/fba/chemical_task.py", "src/ystwin/fba/dynamic_rates.py",
    "scripts/run_native_training.py", "scripts/evaluate_frozen_transfer.py",
    "src/ystwin/generator/in_silico.py", "src/ystwin/analysis/in_silico.py",
    "src/ystwin/in_silico.py",
)
_WESTERNS = frozenset({"Hog1PP_measured", "Gpd1_measured"})
_SOURCE_KEY = ("source_doi", "supplement_id", "source_sheet", "source_cell")
_PLAN = {
    "schema_version": 1,
    "checkpoint_sha256": CHECKPOINT_SHA256,
    "scope": "Retrospective native-data falsification and software/source audit, not independent biological validation",
    "cutoff_model_time_s": 7200.0,
    "cutoff_rule": "Fit times <= cutoff, score all later times; source priors and historical selection already saw these data",
    "nulls": ["within_trace_shuffle", "within_trace_reverse_time"],
    "seed": 19052013,
    "null_group_columns": ["experiment_id", "observable_id"],
    "null_fit_max_nfev": 30,
    "prefix_fit_max_nfev": 30,
    "search_starts": 1,
    "search_fold_bounds": [0.2, 5.0],
    "baselines": ["static_by_protocol", "affine_time_by_protocol"],
    "baseline_inputs": ["observable_id", "nacl_molar", "time_relative_s"],
    "common_fast_rate_factors": [0.2, 1.0, 5.0],
    "common_rate_slice": "Multiply kv16f_1 and kv16r_1 together; profile Western gains only; other kinetics fixed; not a likelihood confidence interval",
    "local_log_step": 0.001,
    "western_scale_factors": [0.5, 2.0],
    "scale_identity_atol": 1e-10,
    "unit_probe": "Convert all reported mol/L values and SD to mmol/L, then verify rejection without an explicit unit adapter",
    "exchange_null": "Shuffle quantified training values within each observable, retaining censored support locations",
    "performance_thresholds": None,
    "selection_from_diagnostic_results": False,
    "significance_or_confidence_intervals": None,
    "failure_policy": "Retain every attempted fit, simulation, unscored readout and unresolved claim; never choose a successful subset",
    "forbidden_inputs_opened": [],
}
_CLAIMS = {
    "conditional_parameter_estimation": (
        "checkpoint_integrity", "native_cells_bound", "native_fit_bound", "supplied_hog_structure",
    ),
    "empirical_interpolation": (
        "checkpoint_integrity", "exchange_knots_bound", "piecewise_interpolation_verified",
    ),
    "supplied_kinetic_law": ("checkpoint_integrity", "supplied_hog_structure"),
    "independent_native_validation": (
        "independent_native_validation_rows", "whole_intervention_groups_reserved",
        "split_before_all_selection", "historical_lineage_disjoint", "independent_experimental_units",
    ),
    "independent_mechanistic_evidence": (
        "independent_native_validation_rows", "whole_intervention_groups_reserved",
        "split_before_all_selection", "historical_lineage_disjoint", "independent_experimental_units",
        "competing_mechanisms_discriminated", "independent_observation_calibration",
    ),
    "teacher_recovery_is_biological_law_learning": ("teacher_targets_are_empirical",),
    "absolute_phosphorylation_fraction_learned": ("independent_observation_calibration",),
    "individual_fast_rates_identified": ("interventional_fast_rate_identification",),
    "native_gem_coupling_validated": ("independent_cross_context_flux_test",),
    "genotype_copy_to_active_enzyme_validated": ("measured_copy_expression_activity_mapping",),
    "realized_allocation_learned": ("measured_allocation_rule", "independent_allocation_validation"),
    "integrity_established": ("checkpoint_integrity",),
    "portable_scientific_content": (
        "portable_scientific_content_integrity", "native_cells_bound", "native_fit_bound",
        "supplied_hog_structure", "exchange_knots_bound", "piecewise_interpolation_verified",
    ),
}
_SCOPES = {
    "conditional_parameter_estimation": "Three grouped kinetic multipliers and two observation gains, conditional on supplied structure, source priors, initial states and observation processing",
    "empirical_interpolation": "Stored training knots and piecewise-linear interpolation at an imposed growth rate; not discovery of growth or exchange mechanisms",
    "supplied_kinetic_law": "Execution of published SBML equations, not learning their functional form",
    "integrity_established": "File identity only; neither independence, identifiability nor causality follows",
    "portable_scientific_content": "Preserved native scientific model content, source-bound values and unchanged relative code/data identities from a verified public export; not original checkpoint-byte integrity, original serialization timing, or independent biological validation",
}
_VERIFIED_GRADES = frozenset({
    "conditional_parameter_estimation", "empirical_interpolation", "supplied_kinetic_law",
    "integrity_established", "portable_scientific_content",
})
_GATE_CONTRACT = {
    "version": 2,
    "advisory_interface": "assess_claim never authorizes; predicate/reference shape is caller-supplied metadata only",
    "verified_interface": "require_claim re-collects authoritative source-bound evidence at root, compares every collector-owned JSON field, and evaluates only implemented scoped grades",
    "root": "Explicit function argument or the installed audit module's repository; never taken from a caller ledger or reference",
    "implemented_grades": sorted(_VERIFIED_GRADES),
    "unimplemented_grades": sorted(set(_CLAIMS) - _VERIFIED_GRADES),
    "stronger_grades": "Fail closed until dedicated typed measurement, split, intervention, calibration and score assessors are implemented; True flags cannot substitute",
    "report_trust": "Pure JSON is a collector snapshot, not a bearer authorization. Extra diagnostic/report fields are not evidence used by the release gate",
    "digest_scope": "No caller-recomputable ledger digest authorizes a claim. The externally retained checkpoint digest binds file identity, not biological truth",
    "assumptions": "Trusted running audit code and dependencies, and a stable filesystem during validation; not protection against a compromised interpreter or concurrent file replacement",
}


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if value is pd.NA or (isinstance(value, float) and not np.isfinite(value)):
        return None
    return value


def _write(path, value):
    with path.open("x") as handle:
        json.dump(_jsonable(value), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field: {key}")
        result[key] = value
    return result


def _json(payload):
    return json.loads(payload, object_pairs_hook=_pairs)


def _bytes(root, relative):
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError("audit inputs must be explicit repository-relative files")
    path = root / candidate
    if path.resolve() != path or not path.is_file():
        raise ValueError(f"audit input is not a regular nonsymlink file: {relative}")
    return path.read_bytes()


def _digest(payload):
    return hashlib.sha256(payload).hexdigest()


def _ref(root, path, location):
    return {"path": path, "sha256": _digest(_bytes(root, path)), "location": location}


def _definition(tree, name):
    node = tree
    for part in name.split("."):
        matches = [item for item in node.body
                   if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                   and item.name == part]
        if len(matches) != 1:
            raise ValueError(f"missing or ambiguous audited source definition: {name}")
        node = matches[0]
    return node


def _contains(node, expression):
    expected = ast.dump(ast.parse(expression, mode="eval").body)
    return any(ast.dump(item) == expected for item in ast.walk(node))


def _symbols(node):
    own = {node.getName()} if node.getType() in (libsbml.AST_NAME, libsbml.AST_NAME_TIME) else set()
    for index in range(node.getNumChildren()):
        own.update(_symbols(node.getChild(index)))
    return own


def _fact(value, refs, basis):
    return {"value": value, "evidence": refs, "verification": basis}


def assess_claim(claim, evidence):
    if claim not in _CLAIMS:
        raise ValueError(f"unknown learning claim: {claim}")
    predicates = []
    for name in _CLAIMS[claim]:
        fact = evidence.get(name, {})
        valid = (isinstance(fact, dict) and fact.get("verification") in {
            "source_structure", "source_bound_values", "frozen_payload_equality", "source_qualification",
            "portable_content_equality",
        } and isinstance(fact.get("evidence"), list) and bool(fact["evidence"])
            and all(isinstance(ref, dict) and set(ref) == {"path", "sha256", "location"}
                    and isinstance(ref["sha256"], str) and len(ref["sha256"]) == 64
                    and set(ref["sha256"]) <= set("0123456789abcdef")
                    and isinstance(ref["path"], str) and bool(ref["path"])
                    and isinstance(ref["location"], dict) and bool(ref["location"]) for ref in fact["evidence"]))
        value = fact.get("value") if valid else None
        state = "pass" if value is True else "fail" if value is False else "unresolved"
        predicates.append({"predicate": name, "state": state, "evidence": fact if valid else None})
    consistent = all(item["state"] == "pass" for item in predicates)
    if claim not in _VERIFIED_GRADES:
        predicates.append({"predicate": "implemented_authoritative_assessor", "state": "unresolved",
                           "evidence": None, "reason": "No implemented assessor for this evidence grade"})
    return {
        "claim": claim, "supported": False, "authorized": False,
        "assessment_kind": "advisory_metadata", "metadata_consistent": consistent,
        "status": "rejected" if any(item["state"] == "fail" for item in predicates)
        else "advisory_consistent_not_verified" if consistent and claim in _VERIFIED_GRADES else "not_established",
        "scope": _SCOPES.get(claim, "Requires evidence not supplied by a good training score, a source-kind flag, or file hashes"),
        "predicates": predicates,
    }


def _collector_assessment(claim, evidence):
    assessment = assess_claim(claim, evidence)
    supported = claim in _VERIFIED_GRADES and assessment["metadata_consistent"]
    return {**assessment, "assessment_kind": "authoritative_collector_snapshot", "supported": supported,
            "status": "supported_with_scope" if supported else assessment["status"]}


def _canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def require_claim(ledger, claim, *, root=None, artifact_source="original", manifest_path=None, expected_sha256=None):
    if claim not in _CLAIMS:
        raise ValueError(f"unknown learning claim: {claim}")
    if claim not in _VERIFIED_GRADES:
        raise ValueError(f"unsupported claim {claim}: no implemented authoritative evidence assessor")
    root = Path(__file__).resolve().parents[3] if root is None else Path(root).resolve(strict=True)
    current = collect_learning_evidence(root, artifact_source=artifact_source,
                                        manifest_path=manifest_path, expected_sha256=expected_sha256)
    if not isinstance(ledger, dict):
        raise ValueError("ledger must match the current authoritative collector")
    try:
        changed = [key for key, value in current.items()
                   if key not in ledger or _canonical_json(ledger[key]) != _canonical_json(value)]
    except (TypeError, ValueError) as exc:
        raise ValueError("ledger is not finite JSON from the authoritative collector") from exc
    if changed:
        raise ValueError(f"ledger differs from the current authoritative collector: {', '.join(changed)}")
    assessment = _collector_assessment(claim, current["evidence_predicates"])
    if not assessment["supported"]:
        failed = [item["predicate"] for item in assessment["predicates"] if item["state"] != "pass"]
        raise ValueError(f"unsupported claim {claim}: {', '.join(failed)}")
    return {**assessment, "assessment_kind": "verified_release_gate", "authorized": True,
            "verification_root": "." if artifact_source == "portable" else str(root),
            "checkpoint_sha256": CHECKPOINT_SHA256,
            **({"artifact_source": "portable", "integrity_scope": current["integrity_scope"]}
               if artifact_source == "portable" else {})}


def bind_native_observations(rows, authoritative):
    required = set(authoritative.columns) - {"split"}
    if (not isinstance(rows, pd.DataFrame) or rows.empty or not rows.columns.is_unique
            or required - set(rows.columns)):
        return {"verified": False, "rows": 0, "issues": [{"reason": "incomplete source-bound observation fields"}]}
    source = authoritative.set_index(list(_SOURCE_KEY), verify_integrity=True)
    issues = []
    if rows.duplicated(list(_SOURCE_KEY)).any():
        issues.append({"reason": "duplicate physical source observation"})
    for row in rows.to_dict("records"):
        key = tuple(row[field] for field in _SOURCE_KEY)
        if key not in source.index:
            issues.append({"source_key": key, "reason": "unknown physical source observation"})
            continue
        original = source.loc[key]
        for field in sorted(required - set(_SOURCE_KEY)):
            actual, expected = row[field], original[field]
            both_missing = pd.isna(actual) and pd.isna(expected)
            equal = not pd.isna(actual) and not pd.isna(expected) and actual == expected
            if not (both_missing or equal):
                issues.append({"source_key": key, "field": field,
                               "observed": _jsonable(actual), "authoritative": _jsonable(expected)})
    return {"verified": not issues, "rows": len(rows), "issues": issues,
            "scope": "Every imported source field except the explicitly replaceable analysis split; metadata flags alone do not bind values"}


def _select(source, protocol):
    selection = protocol["calibration_selection"]
    rows = source.observations
    chosen = rows.loc[
        rows.supplement_id.eq(selection["supplement_id"])
        & rows.genotype.eq(selection["genotype"])
        & rows.strain_id.eq(selection["strain_id"])
        & rows.nacl_molar.isin(selection["nacl_molar"])
        & rows.observable_id.isin(selection["observable_ids"])
    ].copy().reset_index(drop=True)
    if len(chosen) != selection["expected_rows"]:
        raise ValueError("audited native selection differs from the frozen protocol")
    chosen["source_split"] = chosen["split"]
    chosen["split"] = "fit"
    return chosen


def _integrity(root, frozen, checkpoint_digest, *, checkpoint_reference=None):
    entries = frozen["training_manifest"]
    if set(entry["path"] for entry in entries["inputs"]) != set(_NATIVE_INPUTS):
        raise ValueError("checkpoint native inputs differ from the audited explicit allowlist")
    reference = checkpoint_reference or {"path": CHECKPOINT_PATH, "sha256": CHECKPOINT_SHA256}
    results = [{"path": reference["path"], "expected_sha256": reference["sha256"],
                "actual_sha256": checkpoint_digest, "matches": checkpoint_digest == reference["sha256"]}]
    for entry in (*entries["inputs"], *entries["code"]):
        if entry["role"] not in {"fit", "prior", "provenance", "selection", "training_code"}:
            raise ValueError("unexpected input role in native checkpoint")
        actual = _digest(_bytes(root, entry["path"]))
        results.append({"path": entry["path"], "role": entry["role"],
                        "expected_sha256": entry["sha256"], "actual_sha256": actual,
                        "matches": actual == entry["sha256"]})
    return results


def _collect(root, *, artifact_source="original", manifest_path=None, expected_sha256=None):
    root = Path(root).resolve(strict=True)
    portable = None
    if artifact_source == "original":
        if manifest_path is not None or expected_sha256 is not None:
            raise ValueError("portable manifest arguments require artifact_source='portable'")
        reference = {"path": CHECKPOINT_PATH, "sha256": CHECKPOINT_SHA256}
        checkpoint_bytes = _bytes(root, CHECKPOINT_PATH)
        if _digest(checkpoint_bytes) != CHECKPOINT_SHA256:
            raise ValueError("this audit requires the retained frozen checkpoint SHA-256")
        frozen = _json(checkpoint_bytes)
    elif artifact_source == "portable":
        from .portable_replay import load_portable_native_checkpoint

        frozen, portable = load_portable_native_checkpoint(
            root, manifest_path=manifest_path, expected_sha256=expected_sha256)
        reference = portable["checkpoint_reference"]
        checkpoint_bytes = _bytes(root, reference["path"])
    else:
        raise ValueError("artifact_source must be 'original' or 'portable'")
    integrity = _integrity(root, frozen, _digest(checkpoint_bytes), checkpoint_reference=reference)
    if not all(row["matches"] for row in integrity):
        changed = [row["path"] for row in integrity if not row["matches"]]
        raise ValueError(f"frozen dependencies changed; do not replay this checkpoint: {changed}")
    config = frozen["model"]["config"]
    protocol = config["protocol"]
    source = load_native_hog_data(root / "data/hog2013")
    selected = _select(source, protocol)
    model = HogModel.from_source(root / "data/hog2013")
    native = load_native_chemostat_data(root / protocol["native_exchange_fit"]["input"])
    fit_path = "outputs/native_training_run_02/native_hog_fit.json" if portable is None else portable["native_hog_fit_reference"]["path"]
    fit_bytes = _bytes(root, fit_path)
    fit = _json(fit_bytes)
    if portable is not None:
        if _digest(fit_bytes) != portable["native_hog_fit_reference"]["sha256"]:
            raise ValueError("portable native fit export changed during collection")
        fit = fit["payload"]
    parameters = {key: frozen["model"]["parameters"][f"hog.{key}"] for key in config["hog_override_keys"]}
    gains = {name: frozen["model"]["parameters"][f"observation_gain.{name}"] for name in _WESTERNS}
    evidence = {}
    payload_pointer = "" if portable is None else "/payload"
    cp_ref = _ref(root, reference["path"], {"json_pointer": payload_pointer + "/training_manifest"})
    protocol_ref = _ref(root, "data/hog2013/native_training_protocol.json", {"lines": [13, 84]})
    hog_ref = _ref(root, "data/hog2013/observations.xls", {
        "cells": selected.loc[:, list(_SOURCE_KEY)].to_dict("records"),
    })
    fit_ref = _ref(root, fit_path, {"json_pointer": payload_pointer + "/diagnostics"})
    source_ref = _ref(root, "data/hog2013/sources.json", {"json_pointer": "/source_assumptions"})
    binding = bind_native_observations(selected, source.observations)
    evidence["checkpoint_integrity"] = _fact(True if portable is None else None, [cp_ref], "frozen_payload_equality")
    if portable is not None:
        manifest_ref = _ref(root, portable["manifest_reference"]["path"], {"scope": "public export integrity and declared transforms"})
        evidence["portable_scientific_content_integrity"] = _fact(
            portable["scientific_content_preserved"], [cp_ref, manifest_ref], "portable_content_equality")
    evidence["native_cells_bound"] = _fact(binding["verified"], [hog_ref, protocol_ref], "source_bound_values")
    fit_bound = (fit["parameters"] == parameters and fit["observation_gains"] == gains
                 and fit["diagnostics"]["parameter_groups"] == protocol["candidate_parameter_groups"]
                 and fit["diagnostics"]["training_rows"] == len(selected)
                 and set(fit["diagnostics"]["training_experiments"]) == set(selected.experiment_id)
                 and set(protocol["parameter_selection"]["selected_groups"]) == set(protocol["candidate_parameter_groups"])
                 and set(parameters) == {name for members in protocol["candidate_parameter_groups"].values() for name in members}
                 and all(frozen["model"]["parameters"].get(f"hog.{key}") == value
                         for key, value in model.parameter_values.items() if key not in parameters)
                 and all(np.allclose([parameters[key] / model.parameter_values[key] for key in members],
                                     fit["diagnostics"]["parameter_multipliers"][group], rtol=1e-12, atol=0)
                         for group, members in protocol["candidate_parameter_groups"].items()))
    evidence["native_fit_bound"] = _fact(fit_bound, [fit_ref, cp_ref, protocol_ref], "frozen_payload_equality")
    trees, imports, references = {}, [], {}
    for path in _AUDITED_CODE:
        tree = ast.parse(_bytes(root, path), filename=path)
        trees[path] = tree
        for item in ast.walk(tree):
            if isinstance(item, (ast.Import, ast.ImportFrom)):
                imports.append({"path": path, "line": item.lineno, "statement": ast.unparse(item)})

    def definition(path, name):
        node = _definition(trees[path], name)
        references[name] = _ref(root, path, {"symbol": name, "lines": [node.lineno, node.end_lineno]})
        return node

    kinetic = definition("src/ystwin/mech/kinetic_sbml.py", "KineticModel.__init__")
    fit_node = definition("src/ystwin/analysis/hog_learning.py", "fit_native_hog_parameters")
    structure_fixed = (_contains(kinetic, "law.getMath()")
                       and _contains(fit_node, "model.simulate(grid, protocol=protocol, parameters=parameters, rtol=1e-8, atol=1e-11)"))
    evidence["supplied_hog_structure"] = _fact(
        structure_fixed, [references["KineticModel.__init__"], references["fit_native_hog_parameters"]], "source_structure")
    interpolation = definition("src/ystwin/analysis/native_physiology.py", "NativeExchangeModel.predict")
    payload = config["native_exchange_model"]
    if portable is None:
        reconstructed = fit_native_exchange_model(native.split()["train"]).to_dict()
    else:
        from .portable_replay import restore_native_exchange_model

        reconstructed = restore_native_exchange_model(payload, native).to_dict()
    evidence["exchange_knots_bound"] = _fact(payload == reconstructed, [cp_ref, _ref(
        root, protocol["native_exchange_fit"]["input"], {"lines": [1, 11]})], "source_bound_values")
    evidence["piecewise_interpolation_verified"] = _fact(
        payload["method"] == "piecewise_linear" and _contains(
            interpolation, "(1 - fraction) * left.value + fraction * right.value"),
        [references["NativeExchangeModel.predict"], cp_ref], "source_structure")
    teacher = definition("src/ystwin/generator/in_silico.py", "simulate_teacher")
    student = definition("src/ystwin/analysis/in_silico.py", "fit_student")
    teacher_verified = (_contains(teacher, "teacher_controls(latent)") and _contains(student, "episode.latent")
                        and _contains(student, "episode.latent_derivative") and _contains(student, "episode.controls"))
    definition("src/ystwin/generator/in_silico.py", "TeacherParameters")
    evidence["teacher_targets_are_empirical"] = _fact(
        False if teacher_verified else None,
        [references[name] for name in ("simulate_teacher", "fit_student", "TeacherParameters")], "source_structure")
    definition("src/ystwin/analysis/hog_data.py", "_provenance")
    original_fit_reuse = selected.original_fit_status.eq("published_fitting_dataset").all()
    for name, value, refs in (
        ("independent_native_validation_rows", False if original_fit_reuse else None, [hog_ref, source_ref, protocol_ref]),
        ("whole_intervention_groups_reserved", False, [protocol_ref]),
        ("split_before_all_selection", None, [protocol_ref, source_ref]),
        ("historical_lineage_disjoint", False if original_fit_reuse else None, [source_ref, hog_ref]),
        ("independent_experimental_units", None, [references["_provenance"]]),
        ("competing_mechanisms_discriminated", None, [references["fit_native_hog_parameters"], protocol_ref]),
        ("independent_observation_calibration", False, [source_ref, references["_provenance"]]),
        ("interventional_fast_rate_identification", None, [fit_ref, protocol_ref]),
    ):
        evidence[name] = _fact(value, refs, "source_qualification")
    carbon = definition("src/ystwin/mech/hog.py", "HogGlycerolBalance.carbon_commitment")
    coupling = definition("src/ystwin/fba/native_obligations.py", "impose_native_obligation")
    coupling_verified = (_contains(carbon, "0.5 * (synthesis / glucose)")
                         and _contains(coupling, "obligation.carbon_fraction * 6. / 3.")
                         and _contains(coupling, "terminal - coefficient * uptake"))
    evidence["independent_cross_context_flux_test"] = _fact(
        None, [references["HogGlycerolBalance.carbon_commitment"], references["impose_native_obligation"], protocol_ref], "source_structure")
    enzyme = definition("src/ystwin/fba/chemical_task.py", "EnzymeDefinition")
    definition("src/ystwin/fba/chemical_task.py", "install_chemical_task")
    enzyme_fields = [node.target.id for node in enzyme.body if isinstance(node, ast.AnnAssign)]
    evidence["measured_copy_expression_activity_mapping"] = _fact(
        None, [references["EnzymeDefinition"], references["install_chemical_task"]], "source_structure")
    envelope = definition("src/ystwin/analysis/blind_transfer.py", "predict_balanced_envelope")
    directions = [ast.literal_eval(item.iter) for item in ast.walk(envelope)
                  if isinstance(item, ast.For) and isinstance(item.target, ast.Name)
                  and item.target.id == "direction" and isinstance(item.iter, ast.Tuple)]
    envelope_verified = directions == [("min", "max")] and config["prediction_kind"] == "feasible_envelope"
    for name in ("measured_allocation_rule", "independent_allocation_validation"):
        evidence[name] = _fact(False if envelope_verified and name == "measured_allocation_rule" else None,
                               [references["predict_balanced_envelope"], cp_ref], "source_structure")
    km = model.kinetic_model.metadata
    od = km["assignment_formulas"]["OD"]
    od_symbols = sorted(_symbols(libsbml.parseL3Formula(od)))
    estimated = set(parameters)
    source_parameters = model.parameter_values
    units = km["units"]
    components = [
        {"id": "teacher_student", "classification": "synthetic_teacher_recovery",
         "structure_verified": teacher_verified, "biological_law_credit": False,
         "estimated": ["encoder", "dynamics", "controller"],
         "supplied": ["latent coordinate meanings", "teacher dynamics", "sensor loadings", "synthetic control targets"],
         "evidence": [references[name] for name in ("TeacherParameters", "simulate_teacher", "fit_student")]},
        {"id": "native_hog", "classification": "conditional_parameter_estimation",
         "estimated_parameter_groups": protocol["candidate_parameter_groups"], "estimated_observation_gains": gains,
         "kinetic_degrees_of_freedom": len(protocol["candidate_parameter_groups"]),
         "absolute_parameter_overrides": parameters, "source_global_parameter_count": len(source_parameters),
         "remaining_fixed_parameter_count_including_protocol_inputs": len(set(source_parameters) - estimated),
         "remaining_fixed_parameters": {key: source_parameters[key] for key in sorted(set(source_parameters) - estimated)},
         "source_species_count": len(model.kinetic_model.species_ids),
         "source_reaction_count": len(model.kinetic_model.reaction_ids),
         "assignment_rules_count": len(km["assignment_formulas"]),
         "historical_fit_status_counts": selected.original_fit_status.value_counts().to_dict(),
         "native_rows": len(selected), "unshocked_rows_used_for_fit": int(selected.nacl_molar.eq(0).sum()),
         "source_binding": binding, "dependencies": ["source_hog_structure", "source_observation_and_growth_assumptions"],
         "evidence": [fit_ref, protocol_ref, hog_ref]},
        {"id": "source_hog_structure", "classification": "supplied_kinetic_law",
         "functional_forms_learned": False, "structure_verified": structure_fixed,
         "reaction_formulas": km["reaction_formulas"], "assignment_formulas": km["assignment_formulas"],
         "evidence": [references["KineticModel.__init__"], source_ref]},
        {"id": "source_observation_and_growth_assumptions", "classification": "supplied_nuisance_assumptions",
         "source_assumptions": model.source_metadata["source_assumptions"],
         "OD_formula": od, "OD_symbol_dependencies": od_symbols,
         "OD_is_time_only_assignment": set(od_symbols) == {"time"},
         "source_initial_state_policy": fit["diagnostics"]["initial_state_policy"],
         "normalization_policy": source.metadata["normalization_policy"],
         "declared_sbml_volume": units["definitions"].get("volume"),
         "substance_convention": units["substance"], "numeric_unit_conversion": units["numeric_conversion"],
         "absolute_protein_calibration_independent": False,
         "evidence": [source_ref, references["_provenance"], _ref(root, "data/hog2013/model_wt.xml", {"lines": [2476, 2620]})]},
        {"id": "native_exchange", "classification": "empirical_interpolation", "method": payload["method"],
         "training_growth_rates_per_h": payload["training_growth_rates_per_h"],
         "validation_growth_rates_per_h": protocol["native_exchange_fit"]["native_validation_growth_rates_per_h"],
         "growth_role": payload["input"]["role"], "study": native.metadata["source"],
         "split_scope": "Retrospective within-study condition check; no independent-laboratory or ancestral-GEM fit independence established",
         "dependencies": ["van_hoek_measurements"], "evidence": evidence["exchange_knots_bound"]["evidence"]},
        {"id": "native_gem_bridge", "classification": "untested_coupling", "structure_verified": coupling_verified,
         "source_ratio": "0.5*(v6+v6b)/v1 on the common reference-volume rate basis",
         "gem_constraint": "native_terminal >= 2*carbon_fraction*net_glucose_input",
         "identified_from_observations": False, "absolute_flux_or_atom_origin_identified": False,
         "scenario_range_is_confidence_interval": False,
         "dependencies": ["native_hog", "native_exchange", "host_prior"],
         "evidence": evidence["independent_cross_context_flux_test"]["evidence"]},
        {"id": "host_prior", "classification": "supplied_host_prior", "historical_fit_lineage": "incomplete",
         "upstream_fit_dataset_overlap": "unresolved, not assumed absent",
         "repository_model_selection": "Existing manifest records native van Hoek comparison during host-model selection; exact upstream fit-row lineage is not recorded in the checkpoint",
         "evidence": [cp_ref, _ref(root, "data/gem/MANIFEST.md", {"lines": [49, 73]})]},
        {"id": "chemical_capacity", "classification": "untested_coupling", "enzyme_input_fields": enzyme_fields,
         "genotype_copy_to_activity_mapping": "absent from this interface",
         "source_string_is_independent_measurement_proof": False,
         "dependencies": ["supplied_task_stoichiometry", "supplied_kcat_and_active_abundance", "host_prior"],
         "evidence": evidence["measured_copy_expression_activity_mapping"]["evidence"]},
        {"id": "output", "classification": "feasible_envelope", "structure_verified": envelope_verified,
         "quantity": config["quantity"], "unit": config["unit"], "realized_allocation_identified": False,
         "dependencies": ["native_gem_bridge", "chemical_capacity"],
         "evidence": [references["predict_balanced_envelope"], cp_ref]},
        {"id": "independent_mechanistic_evidence", "classification": "independent_mechanistic_evidence",
         "available": False, "evidence": [], "reason": "No new independent, interventional, calibration-bound mechanism comparison is supplied"},
    ]
    lineage = [{"path": entry["path"], "role": entry["role"], "source": entry["source"]}
               for entry in frozen["training_manifest"]["inputs"]]
    ledger = {
        "schema_version": 2, "kind": "biology_learning_audit", "checkpoint_sha256": CHECKPOINT_SHA256,
        "scope": _PLAN["scope"], "components": components, "evidence_predicates": evidence,
        "claims": [_collector_assessment(claim, evidence) for claim in _CLAIMS],
        "gate_trust_contract": deepcopy(_GATE_CONTRACT),
        "additional_qualifications": _additional_qualifications(root, model),
        "native_input_lineage": lineage, "source_imports": imports, "integrity_before": integrity,
        "native_observation_inventory": {
            "imported_rows": len(source.observations), "fitted_rows": len(selected),
            "unselected_rows": len(source.observations) - len(selected),
            "original_fit_status_counts": source.observations.original_fit_status.value_counts().to_dict(),
            "series": source.observations.groupby(
                ["genotype", "nacl_molar", "observable_id", "original_fit_status"], dropna=False
            ).size().reset_index(name="rows").to_dict("records"),
            "unselected_is_not_independent": True,
            "s004_scope": "Not opened; existing source qualification says preliminary data participated in model development despite not entering the published parameter fit",
        },
        "qualification": "Predicates are constructed from this retained checkpoint, verified native cells and inspected source structure. Caller-authored source-kind flags are not evidence. Passing a necessary gate never proves causality.",
        "excluded": ["product outcomes", "product calibration tables", "S4 outcome values", "new source/model selection"],
    }
    if portable is not None:
        ledger["artifact_source"] = "portable"
        ledger["evidence_mode"] = portable["mode"]
        ledger["integrity_scope"] = portable
        ledger["excluded"] = ["product-outcome-based fitting or selection", "product calibration tables",
                              "S4 outcome values", "new source/model selection"]
        ledger["source_access_scope"] = (
            "Whole-public-bundle integrity validation also parses previously recorded outcomes and scores. "
            "They are not inputs to native parameter estimation or to this scoped native-content grade."
        )
        ledger["qualification"] = (
            "Public export bytes, preserved scientific payload and relative pinned dependencies are checked. "
            "The private original checkpoint bytes are not checked and checkpoint_integrity remains unresolved. "
            "This is replay of previously recorded content, not independent biological evidence or proof of original serialization timing."
        )
    return {"root": root, "ledger": ledger, "model": model, "selected": selected, "source": source,
            "native": native, "frozen": frozen, "parameters": parameters, "gains": gains, "fit": fit}


def collect_learning_evidence(root, *, artifact_source="original", manifest_path=None, expected_sha256=None):
    return _jsonable(_collect(root, artifact_source=artifact_source,
                             manifest_path=manifest_path, expected_sha256=expected_sha256)["ledger"])


def _additional_qualifications(root, model):
    qualifications = {}
    gauge_path = "outputs/biology_observability_gauge.json"
    test_path = "tests/test_native_observability.py"
    candidate_path = "data/validation_candidates/granados2018.json"
    if (root / gauge_path).is_file():
        gauge = _json(_bytes(root, gauge_path))
        qualifications["gpd1_source_gauge"] = {
            "classification": "reported_source_structural_counterexample_not_empirical_validation",
            "source_model_identity_matches": gauge["source_model"] == {
                "path": "data/hog2013/model_wt.xml", "sha256": model.kinetic_model.metadata["sha256"],
            },
            "gauge": gauge["gauge"], "reported_probe": gauge["actual_probe"],
            "within_restricted_three_group_fit": False,
            "numerical_probe_rerun_by_collector": False, "biological_validation_credit": False,
            "qualification": "Scaling initial Gpd1 and kv18f_1 by s while dividing kv6_1 and the Gpd1 observation gain by s preserves the fitted readouts. The restricted fit fixes or ties these directions; that is prior dependence, not identification by new measurements.",
            "recorded_test_count": gauge["tests"]["passed"],
            "test_count_qualification": "The gauge artifact records its earlier test count; the separately referenced parent test file also contains the pool-counterflow test. No tests are reexecuted by this collector.",
            "evidence": [_ref(root, gauge_path, {"json_pointer": "/gauge"}),
                         _ref(root, gauge_path, {"json_pointer": "/actual_probe"})],
        }
    if (root / test_path).is_file():
        tree = ast.parse(_bytes(root, test_path), filename=test_path)
        name = "test_pool_balances_cannot_resolve_common_inward_and_outward_counterflow"
        node = _definition(tree, name)
        qualifications["pool_counterflow"] = {
            "classification": "referenced_parent_source_balance_test_not_empirical_validation",
            "null_directions": [["v13a", "v13b"], ["v13aBatch", "v13bBatch"]],
            "kinetic_parameter_equivalence_established": False, "biological_validation_credit": False,
            "qualification": "Equal extra inward/outward stoichiometric counterflow leaves bulk pool balances unchanged. Supplied kinetic forms may restrict this direction; this is not proof of equivalent parameterizations within the restricted kinetic fit.",
            "evidence": [_ref(root, test_path, {"symbol": name, "lines": [node.lineno, node.end_lineno]})],
        }
    if (root / candidate_path).is_file():
        candidate = _json(_bytes(root, candidate_path))
        qualifications["granados2018_candidate"] = {
            "classification": "metadata_only_future_validation_candidate",
            "catalogue_status": candidate["status"], "source": candidate["source"],
            "measurement_contract": candidate["measurement_contract"],
            "reported_design": ["biological replicates", "rich-to-rich controls", "co-measured pairs", "dose variation"],
            "measurement_arrays_opened_by_audit": False, "validation_completed": False,
            "biological_validation_credit": False,
            "leakage_and_interpretation_gates": candidate["leakage_and_interpretation_gates"],
            "catalogue_discrepancies": candidate["catalogue_discrepancies"],
            "qualification": "Catalogue existence, later publication date and a metadata-verification flag do not establish independence or causality. Replicate-1 duplication, exact normalization, localization-to-model observation mapping, full experimental splits and ancestral lineage must be resolved before opening measurement arrays for a registered test.",
            "evidence": [_ref(root, candidate_path, {"json_pointer": "/measurement_contract"}),
                         _ref(root, candidate_path, {"json_pointer": "/leakage_and_interpretation_gates"})],
        }
    return qualifications


def _raw_hog(model, rows, parameters):
    prediction = np.empty(len(rows))
    groups = {}
    for index, row in rows.reset_index(drop=True).iterrows():
        key = (float(row.nacl_molar), float(row.time_model_s - row.time_relative_s))
        groups.setdefault(key, []).append(index)
    trajectories = []
    for (dose, shock), indices in groups.items():
        selected = rows.iloc[indices]
        grid = np.unique(np.r_[0.0, selected.time_model_s.to_numpy(dtype=float)])
        trajectory = model.simulate(grid, protocol=HogProtocol(dose, shock), parameters=parameters,
                                    rtol=1e-8, atol=1e-11)
        for index in indices:
            row = rows.iloc[index]
            prediction[index] = trajectory.variables[row.observable_id][np.searchsorted(grid, row.time_model_s)]
        trajectories.append(trajectory)
    if not np.isfinite(prediction).all():
        raise ValueError("native prediction contains nonfinite observations")
    return prediction, trajectories


def _loss(rows, raw, *, gains=None, scales=None):
    values = rows.value.to_numpy(dtype=float)
    prediction = np.asarray(raw, dtype=float).copy()
    if prediction.shape != values.shape or not np.isfinite(prediction).all():
        raise ValueError("complete finite native predictions are required")
    groups = {name: np.flatnonzero(rows.observable_id.to_numpy() == name) for name in rows.observable_id.unique()}
    fitted_gains, normalizers, errors = {}, {}, {}
    residual = np.empty(len(rows))
    for name, indices in groups.items():
        response, target = prediction[indices], values[indices]
        if name in _WESTERNS:
            gain = float(response @ target / (response @ response)) if gains is None else gains[name]
            if not np.isfinite(gain) or gain <= 0:
                raise ValueError("Western profiling requires a positive finite gain")
            fitted_gains[name] = gain
            prediction[indices] *= gain
        scale = float(np.sqrt(np.mean(target ** 2))) if scales is None else scales[name]
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("native loss requires a positive finite training scale")
        normalizers[name] = scale
        difference = prediction[indices] - target
        residual[indices] = difference / scale / np.sqrt(len(indices) * len(groups))
        errors[name] = {"rows": len(indices), "rmse": float(np.sqrt(np.mean(difference ** 2))),
                        "nrmse": float(np.sqrt(np.mean(difference ** 2)) / scale)}
    return {"nrmse": float(np.linalg.norm(residual)), "gains": fitted_gains,
            "training_scales": normalizers, "per_observable": errors}, prediction, residual


def _baselines(train, test):
    outputs = {}
    for name in _PLAN["baselines"]:
        predicted = np.full(len(test), np.nan)
        coefficients = []
        for (observable, dose), group in train.groupby(["observable_id", "nacl_molar"], sort=True):
            t = group.time_relative_s.to_numpy(dtype=float) / 3600
            x = np.ones((len(t), 1)) if name == "static_by_protocol" else np.column_stack([np.ones(len(t)), t])
            weights, _, rank, _ = np.linalg.lstsq(x, group.value.to_numpy(dtype=float), rcond=None)
            mask = test.observable_id.eq(observable) & test.nacl_molar.eq(dose)
            query = test.loc[mask].time_relative_s.to_numpy(dtype=float) / 3600
            design = np.ones((len(query), 1)) if name == "static_by_protocol" else np.column_stack([np.ones(len(query)), query])
            predicted[np.flatnonzero(mask)] = design @ weights
            coefficients.append({"observable_id": observable, "nacl_molar": dose, "rank": int(rank),
                                 "coefficients": weights.tolist()})
        scales = {obs: float(np.sqrt(np.mean(group.value.to_numpy(dtype=float) ** 2)))
                  for obs, group in train.groupby("observable_id")}
        metrics, _, _ = _loss(test, predicted, gains=dict.fromkeys(_WESTERNS, 1.0), scales=scales)
        outputs[name] = {"metrics": metrics, "prediction": predicted.tolist(), "coefficients": coefficients,
                         "parameter_count": sum(item["rank"] for item in coefficients),
                         "scope": "Simple curve baseline, not biological mechanism; count excludes no hidden historically fitted priors"}
    return outputs


def _fit_probe(context, rows, budget):
    protocol = context["frozen"]["model"]["config"]["protocol"]
    fit = fit_native_hog_parameters(
        context["model"], rows, parameter_groups=protocol["candidate_parameter_groups"],
        fold_bounds=tuple(_PLAN["search_fold_bounds"]), max_nfev=budget,
    )
    return fit, {"status": "converged" if fit.diagnostics["optimization_converged"] else "budget_or_convergence_failure",
                 "parameters": fit.parameter_values, "gains": fit.observation_gains,
                 "nrmse": fit.training_nrmse, "diagnostics": fit.diagnostics}


def _failure(exc):
    return {"status": "failed", "exception": type(exc).__name__, "reason": str(exc)}


def _hog_diagnostics(context):
    rows, model = context["selected"], context["model"]
    raw, trajectories = _raw_hog(model, rows, context["parameters"])
    frozen_metrics, predicted, residual = _loss(rows, raw, gains=context["gains"])
    source_raw, _ = _raw_hog(model, rows, {})
    source_metrics, _, _ = _loss(rows, source_raw)
    report = {
        "scope": "Measured native fitting-data diagnostics, conditional on the historically fitted source model; not independent validation",
        "frozen_replay": frozen_metrics, "published_reference_with_reprofiled_training_gains": source_metrics,
        "recorded_fit_nrmse": context["fit"]["diagnostics"]["training_nrmse"],
        "replay_minus_recorded_nrmse": frozen_metrics["nrmse"] - context["fit"]["diagnostics"]["training_nrmse"],
        "native_rows": [{"source_sheet": row.source_sheet, "source_cell": row.source_cell,
                         "experiment_id": row.experiment_id, "observable_id": row.observable_id,
                         "nacl_molar": row.nacl_molar, "time_model_s": row.time_model_s,
                         "value": row.value, "frozen_prediction": predicted[index]}
                        for index, row in enumerate(rows.itertuples(index=False))],
        "in_sample_baselines": _baselines(rows, rows), "null_controls": [],
    }
    common_times = sorted(set.intersection(*(set(trajectory.times_s) for trajectory in trajectories)))
    od = [np.interp(common_times, trajectory.times_s, trajectory.variables["OD"]) for trajectory in trajectories]
    report["OD_driver"] = {"times_s": common_times, "dose_trajectories": [values.tolist() for values in od],
                           "maximum_dose_difference": float(np.max(np.abs(od[0] - od[1]))),
                           "interpretation": "Dose-invariant imposed OD assignment, not learned growth adaptation"}
    rng = np.random.default_rng(_PLAN["seed"])
    for operation in _PLAN["nulls"]:
        derived = rows.copy(deep=True)
        donors = np.arange(len(rows))
        for _, group in rows.groupby(_PLAN["null_group_columns"], sort=True):
            indices = group.sort_values("time_model_s").index.to_numpy()
            donors[indices] = rng.permutation(indices) if operation == "within_trace_shuffle" else indices[::-1]
        derived["value"] = rows.value.to_numpy()[donors]
        derived["sd"] = rows.sd.to_numpy()[donors]
        record = {"operation": operation, "evidence_class": "counterfactual_native_value_control_not_empirical_validation",
                  "source_binding": bind_native_observations(derived, context["source"].observations),
                  "mapping": [{"recipient": [rows.iloc[i].source_sheet, rows.iloc[i].source_cell],
                               "donor": [rows.iloc[j].source_sheet, rows.iloc[j].source_cell]}
                              for i, j in enumerate(donors)],
                  "baselines": _baselines(derived, derived)}
        try:
            _, result = _fit_probe(context, derived, _PLAN["null_fit_max_nfev"])
            record.update(result)
            record["nrmse_minus_native_frozen"] = result["nrmse"] - frozen_metrics["nrmse"]
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            record.update(_failure(exc))
        report["null_controls"].append(record)
    prefix = rows.loc[rows.time_model_s <= _PLAN["cutoff_model_time_s"]].reset_index(drop=True)
    suffix = rows.loc[rows.time_model_s > _PLAN["cutoff_model_time_s"]].reset_index(drop=True)
    retrospective = {
        "scope": "Retrospective blocked-time stress test, not a dose/strain/intervention holdout; original priors and selection saw both partitions",
        "training_rows": len(prefix), "test_rows": len(suffix), "test_refitted": False,
        "baselines": _baselines(prefix, suffix), "source_priors_test_exposed": True,
        "test_source_cells": suffix.loc[:, ["source_sheet", "source_cell", "observable_id", "value"]].to_dict("records"),
    }
    try:
        prefix_fit, result = _fit_probe(context, prefix, _PLAN["prefix_fit_max_nfev"])
        retrospective["fit"] = result
        test_raw, _ = _raw_hog(model, suffix, prefix_fit.parameter_values)
        test_metrics, test_prediction, _ = _loss(
            suffix, test_raw, gains=prefix_fit.observation_gains,
            scales=prefix_fit.diagnostics["normalization_scales"],
        )
        retrospective.update(test_metrics=test_metrics, test_prediction=test_prediction.tolist())
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        retrospective["failure"] = _failure(exc)
    report["retrospective_time_split"] = retrospective
    report.update(_profile_diagnostics(context, rows, raw, residual))
    return report


def _profile_diagnostics(context, rows, raw, residual):
    model, parameters = context["model"], context["parameters"]
    profiled, _, base_residual = _loss(rows, raw)
    groups = dict(context["frozen"]["model"]["config"]["protocol"]["candidate_parameter_groups"])
    groups.update(fast_deactivation=["kv16r_1"], active_uptake=["kv13b_1"])
    probes, columns = [], []
    for name, members in groups.items():
        perturbed = dict(parameters)
        for member in members:
            perturbed[member] = parameters.get(member, model.parameter_values[member]) * np.exp(_PLAN["local_log_step"])
        try:
            prediction, _ = _raw_hog(model, rows, perturbed)
            metrics, _, shifted = _loss(rows, prediction)
            columns.append((shifted - base_residual) / _PLAN["local_log_step"])
            probes.append({"group": name, "status": "ok", "metrics": metrics})
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            columns.append(None)
            probes.append({"group": name, **_failure(exc)})
    sensitivities = {"probes": probes, "scope": "Numerical local sensitivity with gains profiled, not global or causal identifiability"}
    if all(column is not None for column in columns):
        for label, count in (("three_fitted_groups", len(groups) - 2), ("including_two_fixed_rate_priors", len(groups))):
            matrix = np.column_stack(columns[:count])
            _, singular, vectors = np.linalg.svd(matrix, full_matrices=False)
            sensitivities[label] = {
                "coordinates": list(groups)[:count], "singular_values": singular.tolist(),
                "condition_number": float(singular[0] / singular[-1]) if singular[-1] > 0 else None,
                "weakest_direction": dict(zip(list(groups)[:count], vectors[-1].tolist())),
                "numerical_rank": int(np.count_nonzero(singular > max(matrix.shape) * np.finfo(float).eps * singular[0])),
            }
    else:
        sensitivities["status"] = "incomplete; no successful-only sensitivity matrix reported"
    slices, scale_cases = [], []
    for factor in _PLAN["common_fast_rate_factors"]:
        trial = dict(parameters)
        trial["kv16f_1"] = parameters["kv16f_1"] * factor
        trial["kv16r_1"] = model.parameter_values["kv16r_1"] * factor
        record = {"common_rate_factor": factor, "activation_deactivation_ratio": trial["kv16f_1"] / trial["kv16r_1"],
                  "other_kinetics_refitted": False, "biological_validation": False}
        try:
            prediction, trajectories = _raw_hog(model, rows, trial)
            metrics, _, local_residual = _loss(rows, prediction)
            record.update(status="ok", metrics=metrics, nrmse_minus_frozen=metrics["nrmse"] - profiled["nrmse"])
            record["native_carbon_ratio"] = []
            for trajectory in trajectories:
                try:
                    fractions = hog_glycerol_balance(model, trajectory).carbon_commitment()
                    values = fractions[trajectory.times_s >= 3600.0]
                    record["native_carbon_ratio"].append({
                        "nacl_molar": trajectory.metadata["hog_protocol"]["nacl_molar"],
                        "minimum_at_observation_grid": float(values.min()), "maximum_at_observation_grid": float(values.max()),
                        "scope": "Conditional source-law diagnostic, not absolute flux, isotope origin or uncertainty coverage",
                    })
                except (ValueError, RuntimeError) as exc:
                    record["native_carbon_ratio"].append(_failure(exc))
            for scale_factor in _PLAN["western_scale_factors"]:
                perturbed = rows.copy(deep=True)
                mask = perturbed.observable_id.isin(_WESTERNS)
                perturbed.loc[mask, "value"] *= scale_factor
                perturbed.loc[mask, "sd"] *= scale_factor
                scaled_metrics, _, scaled_residual = _loss(perturbed, prediction)
                delta = float(np.max(np.abs(scaled_residual - local_residual)))
                scale_cases.append({"common_rate_factor": factor, "western_scale_factor": scale_factor,
                                    "status": "identity_holds" if delta <= _PLAN["scale_identity_atol"] else "identity_failed",
                                    "maximum_residual_difference": delta,
                                    "gain_ratios": {key: scaled_metrics["gains"][key] / metrics["gains"][key] for key in _WESTERNS},
                                    "absolute_biological_truth_redefined": False})
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            record.update(_failure(exc))
            scale_cases.extend({"common_rate_factor": factor, "western_scale_factor": scale_factor,
                                "status": "not_evaluable", "reason": str(exc)}
                               for scale_factor in _PLAN["western_scale_factors"])
        slices.append(record)
    converted = rows.copy(deep=True)
    converted["value"] *= 1000.0
    converted["sd"] *= 1000.0
    converted["unit"] = "mmol/L"
    try:
        _fit_probe(context, converted, 1)
        unit_probe = {"status": "unexpected_acceptance"}
    except ValueError as exc:
        unit_probe = {"status": "rejected", "reason": str(exc)}
    return {
        "local_sensitivity": sensitivities, "common_fast_rate_slices": slices,
        "western_scale_identity": {
            "cases": scale_cases,
            "equation": "y'=c*y; gain'=c*gain; RMS'=c*RMS, so every gain-profiled standardized residual is unchanged",
            "interpretation": "Absolute Western amplitude cannot identify phosphorylation fractions independently of its assumed calibration; these are counterfactual calibration tests, not new observations",
        },
        "explicit_unit_conversion_probe": {
            **unit_probe, "scope": "The current native fitting interface accepts only mol/L and must not silently fit differently scaled units",
            "source_binding": bind_native_observations(converted, context["source"].observations),
        },
        "profiled_vs_frozen_gain_residual_difference": float(np.max(np.abs(base_residual - residual))),
    }


def _exchange_diagnostics(context):
    parts = context["native"].split()
    train, holdout = parts["train"], parts["holdout"]
    payload = context["frozen"]["model"]["config"]["native_exchange_model"]
    model = NativeExchangeModel.from_dict(payload)
    evaluation = evaluate_native_exchange_model(model, holdout)
    training = evaluate_native_exchange_model(model, train)
    predictions = []
    for row in evaluation["readouts"]:
        native_rows = train.loc[train.observable_id.eq(row["observable_id"]) & train.observation_status.eq("quantified")]
        x = native_rows.growth_rate_per_h.to_numpy(dtype=float)
        y = native_rows.value.to_numpy(dtype=float)
        baselines = {}
        for method in ("static_training_mean", "affine_growth_only"):
            value = None
            if len(y) and (method == "static_training_mean" or len(y) >= 2):
                if method == "static_training_mean":
                    value = float(np.mean(y))
                else:
                    weights = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), y, rcond=None)[0]
                    value = float(np.array([1., row["growth_rate_per_h"]]) @ weights)
            observed = row["observed_value"]
            baselines[method] = {
                "prediction": value, "status": "ok" if value is not None else "insufficient_quantified_training_support",
                "residual": value - observed if value is not None and observed is not None else None,
                "common_with_interpolator_scored": bool(row["scored"] and value is not None),
                "negative_prediction": value is not None and value < 0,
            }
        predictions.append({**row, "baselines": baselines})
    comparisons = []
    for observable in holdout.observable_id.unique():
        subset = [row for row in predictions if row["observable_id"] == observable]
        for method in ("static_training_mean", "affine_growth_only"):
            pairs = [row for row in subset if row["baselines"][method]["common_with_interpolator_scored"]]
            comparisons.append({
                "observable_id": observable, "baseline": method, "n_total": len(subset), "n_common_scored": len(pairs),
                "interpolation_rmse_common": float(np.sqrt(np.mean([row["residual"] ** 2 for row in pairs]))) if pairs else None,
                "baseline_rmse_common": float(np.sqrt(np.mean([row["baselines"][method]["residual"] ** 2 for row in pairs]))) if pairs else None,
            })
    rng = np.random.default_rng(_PLAN["seed"])
    shuffled = train.copy(deep=True)
    mapping = []
    for observable, group in train.groupby("observable_id", sort=True):
        indices = group.loc[group.observation_status.eq("quantified")].index.to_numpy()
        donors = rng.permutation(indices)
        shuffled.loc[indices, "reported_value"] = train.loc[donors, "reported_value"].to_numpy()
        shuffled.loc[indices, "value"] = train.loc[donors, "value"].to_numpy()
        mapping.extend({"observable_id": observable, "recipient_line": int(train.loc[i, "source_line"]),
                        "donor_line": int(train.loc[j, "source_line"])} for i, j in zip(indices, donors))
    null_model = fit_native_exchange_model(shuffled)
    null_training = evaluate_native_exchange_model(null_model, shuffled)
    return {
        "scope": "Retrospective within-study native condition assessment; imposed growth is never scored; upstream host-fit overlap remains unresolved",
        "original_interpolation_training": training,
        "native_condition_assessment": {key: value for key, value in evaluation.items() if key != "readouts"},
        "all_condition_readouts": predictions, "baseline_comparisons_on_common_support": comparisons,
        "shuffled_training_knots": {
            "evidence_class": "counterfactual_native_value_control_not_empirical_validation",
            "mapping": mapping, "training_assessment": null_training,
            "source_values_changed": int(np.count_nonzero(shuffled.reported_value.to_numpy() != train.reported_value.to_numpy())),
            "scope": "An interpolator reproduces even shuffled training knots exactly. Zero training error is not mechanism evidence. Censored support is not replaced by zero targets.",
        },
    }


def _findings(ledger, diagnostics):
    components = {item["id"]: item for item in ledger["components"]}
    descriptions = (
        ("historical_reuse", "high", "Native training is real data, but not independent evidence for the inherited biological laws",
         ["native_hog", "host_prior"],
         "Reserve entire intervention/dose/strain experiments before mechanism, parameter-group, observation-transform and host-model selection. Record ancestral parameter-fit dataset identifiers; acquire new experimental units when independence is unknown."),
        ("coupling_not_learned", "high", "The HOG synthesis/current-glucose ratio becomes a GEM lower-bound obligation by a supplied law, not a measured causal transfer",
         ["native_gem_bridge", "source_hog_structure"],
         "Measure condition-matched isotope-resolved synthesis, transport, glucose uptake, accessible volume and dry mass; test predictions under controlled glycerol-pathway perturbations before using a cross-strain/medium coupling as biological truth."),
        ("allocation_and_activity_gap", "high", "Chemical feasibility and genotype copy information do not determine active enzyme or realized allocation",
         ["chemical_capacity", "output"],
         "Separate genotype, expression, protein amount, active fraction and realized allocation at typed interfaces; require independent dose-response/activity/flux measurements or return only explicitly conditional envelopes."),
        ("observation_and_fixed_prior_dependence", "high", "The fitted numbers depend on supplied nuisance assumptions and constrained parameter combinations",
         ["source_observation_and_growth_assumptions", "native_hog"],
         "Make observation calibration, OD forcing and initial-state assumptions explicit inputs with independent uncertainty. Use high-time-resolution perturbation/relaxation data and competing-law tests; do not infer global identifiability from a reduced local sensitivity rank."),
        ("interpolation_shortcut", "medium", "Native exchange curves interpolate imposed growth; training recovery cannot identify an exchange mechanism",
         ["native_exchange"],
         "Keep an empirical interpolator as a named baseline, preserve censored coverage, and compare a frozen mechanism on independently held-out input conditions without imposing the scored outputs."),
        ("teacher_and_integrity_scope", "medium", "Synthetic teacher recovery and byte integrity are useful engineering tests, not biological-law evidence",
         ["teacher_student"],
         "Attach claim-gate receipts to the next-version fit/evaluation interface; never promote synthetic supervision or hashes into empirical mechanistic credit."),
    )
    findings = []
    for rank, (identifier, severity, finding, ids, action) in enumerate(descriptions, 1):
        findings.append({"rank": rank, "id": identifier, "severity": severity, "finding": finding,
                         "evidence": [ref for key in ids for ref in components[key]["evidence"]],
                         "next_version_minimal_change": action})
    findings[3]["diagnostics"] = {
        "common_rate_slices": diagnostics.get("hog", {}).get("common_fast_rate_slices"),
        "local_sensitivity": diagnostics.get("hog", {}).get("local_sensitivity"),
        "western_scale_identity": diagnostics.get("hog", {}).get("western_scale_identity"),
    }
    qualifications = ledger["additional_qualifications"]
    for index, name in ((3, "gpd1_source_gauge"), (1, "pool_counterflow"), (0, "granados2018_candidate")):
        if name in qualifications:
            findings[index]["additional_qualification"] = name
            findings[index]["evidence"].extend(qualifications[name]["evidence"])
    return findings


def _reused_native_diagnostics(root):
    plan_path = "outputs/biology_learning_audit_01/diagnostic_plan.json"
    source_path = "outputs/biology_learning_audit_01/native_diagnostics.json"
    prior_plan = _json(_bytes(root, plan_path))
    if any(prior_plan.get(key) != value for key, value in _PLAN.items()):
        raise ValueError("prior diagnostic plan differs from the retained checkpoint or diagnostic design")
    payload = _bytes(root, source_path)
    if set(_json(payload)) != {"hog", "native_exchange"}:
        raise ValueError("prior native diagnostic report has unexpected fields")
    return payload, {
        "mode": "reused_unchanged_retrospective_evidence", "numerical_diagnostics_rerun": False,
        "checkpoint_sha256": CHECKPOINT_SHA256, "biological_validation_credit": False,
        "source": {"path": source_path, "sha256": _digest(payload), "location": {"scope": "entire prior native diagnostic report"}},
        "original_plan": _ref(root, plan_path, {"scope": "original design and original audit-code digests"}),
        "qualification": "Prior bytes are copied unchanged, not refitted or newly validated. This digest binds report identity only; the release gate does not use diagnostic scores or this receipt to authorize stronger biological claims.",
    }


def run_biology_learning_audit(root, output_dir, *, reuse_native_diagnostics=False):
    root = Path(root).resolve(strict=True)
    output = Path(output_dir)
    if not output.is_absolute():
        output = root / output
    if (output.parent.resolve() != root / "outputs" or output.resolve() != output
            or not output.name.startswith("biology_learning_audit_")):
        raise ValueError("use a fresh direct outputs/biology_learning_audit_* directory")
    if not isinstance(reuse_native_diagnostics, bool):
        raise ValueError("reuse_native_diagnostics must be an explicit boolean")
    output.mkdir(exist_ok=False)
    plan = deepcopy(_PLAN)
    plan["audit_code"] = [_ref(root, path, {"scope": "audit implementation"}) for path in (
        "src/ystwin/analysis/biology_learning_audit.py", "scripts/audit_biology_learning.py",
    )]
    plan["execution_mode"] = "reuse_prior_retrospective_diagnostics" if reuse_native_diagnostics else "run_native_diagnostics"
    _write(output / "diagnostic_plan.json", plan)
    try:
        context = _collect(root)
        if reuse_native_diagnostics:
            payload, origin = _reused_native_diagnostics(root)
            diagnostics = _json(payload)
            with (output / "native_diagnostics.json").open("xb") as handle:
                handle.write(payload)
        else:
            diagnostics = {}
            for name, operation in (("hog", _hog_diagnostics), ("native_exchange", _exchange_diagnostics)):
                try:
                    diagnostics[name] = operation(context)
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                    diagnostics[name] = _failure(exc)
            _write(output / "native_diagnostics.json", diagnostics)
            origin = {"mode": "new_retrospective_diagnostic_run", "numerical_diagnostics_rerun": True,
                      "biological_validation_credit": False}
        ledger = context["ledger"]
        ledger["diagnostics_origin"] = origin
        ledger["diagnostics_sha256"] = _digest((output / "native_diagnostics.json").read_bytes())
        if reuse_native_diagnostics and ledger["diagnostics_sha256"] != origin["source"]["sha256"]:
            raise ValueError("reused native diagnostic bytes changed during copying")
        ledger["findings"] = _findings(ledger, diagnostics)
        ledger["diagnostic_plan_sha256"] = _digest((output / "diagnostic_plan.json").read_bytes())
        ledger["integrity_after"] = _integrity(root, context["frozen"], _digest(_bytes(root, CHECKPOINT_PATH)))
        ledger["all_frozen_files_unchanged"] = all(row["matches"] for row in ledger["integrity_after"])
        ledger["diagnostics_path"] = str((output / "native_diagnostics.json").relative_to(root))
        ledger["diagnostic_results_do_not_upgrade_claims"] = True
        _write(output / "ledger.json", ledger)
        return _jsonable(ledger)
    except Exception as exc:
        _write(output / "audit_failure.json", _failure(exc))
        raise
