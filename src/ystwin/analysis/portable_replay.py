from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
import hashlib
import importlib.metadata
import importlib.util
import math
from pathlib import Path, PurePosixPath
import platform

from .training_freeze import (
    InputSpec,
    _entries,
    _interval_rows,
    _json_bytes,
    _number,
    _parse_json,
    _reject_label_overlap,
    _sha256,
    _unique_rows,
)


CHECKPOINT_PATH = "outputs/native_training_run_02/checkpoint.json"
CHECKPOINT_SHA256 = "384175e19b2451f542ab08e347ee9d9021e22772fbc205e0c404bbc28416ec94"
MODEL_SHA256 = "97924083113565e1a7db10b07bb735c2819d4fddc69b6a43279d5377195aed79"
REPLAY_KIND = "portable_replay_of_previously_recorded_content"
DEFAULT_MANIFEST_PATH = "data/frozen_evidence/native_v1/manifest.json"
DEFAULT_MANIFEST_SHA256 = "eff317f689fe096f6b3df35e050f557b51cd7e51cdb192bb95775eef8ef529ec"
PROTOCOL_PATH = "data/hog2013/native_training_protocol.json"
NATIVE_INPUT_ROLES = {
    "data/gem/ecYeastGEM_batch.xml.gz": "prior",
    "data/hog2013/methods.pdf": "provenance",
    "data/hog2013/model_wt.xml": "prior",
    PROTOCOL_PATH: "selection",
    "data/hog2013/observations.xls": "fit",
    "data/hog2013/sources.json": "provenance",
    "data/physiology/chemostatData_VanHoek1998.tsv": "fit",
}
FROZEN_CODE_PATHS = frozenset({
    "pyproject.toml", "scripts/evaluate_frozen_transfer.py", "scripts/run_native_training.py",
    "src/ystwin/__init__.py", "src/ystwin/paths.py", "src/ystwin/analysis/__init__.py",
    "src/ystwin/analysis/blind_transfer.py", "src/ystwin/analysis/hog_data.py",
    "src/ystwin/analysis/hog_learning.py", "src/ystwin/analysis/native_physiology.py",
    "src/ystwin/analysis/training_freeze.py", "src/ystwin/fba/__init__.py",
    "src/ystwin/fba/chemical_task.py", "src/ystwin/fba/dynamic_rates.py",
    "src/ystwin/fba/native_obligations.py", "src/ystwin/fba/physiology.py",
    "src/ystwin/fba/solver.py", "src/ystwin/mech/__init__.py", "src/ystwin/mech/hog.py",
    "src/ystwin/mech/kinetic_sbml.py",
})
_VARIANTS = ("nad", "oxygen_peroxide", "oxygen_water")
_WESTERNS = frozenset({"Hog1PP_measured", "Gpd1_measured"})
_ENTRY_FIELDS = {"path", "resolved_path", "sha256", "size_bytes", "role", "source"}
_NUMERICAL_RTOL = 1e-10
_NUMERICAL_ATOL = 1e-14


def _load_bundle(root, manifest_path, expected_sha256):
    from .portable_evidence import load_portable_evidence

    root = Path(root).resolve(strict=True)
    if manifest_path is None:
        manifest_path = DEFAULT_MANIFEST_PATH
        expected_sha256 = DEFAULT_MANIFEST_SHA256 if expected_sha256 is None else expected_sha256
    elif expected_sha256 is None:
        raise ValueError("an explicit portable manifest requires its externally retained SHA-256")
    path = Path(manifest_path)
    if path.is_absolute():
        if not path.is_relative_to(root):
            raise ValueError("portable manifest escapes the caller-supplied root")
        path = path.relative_to(root)
    manifest_path = path.as_posix()
    _relative_file(root, manifest_path)
    bundle = load_portable_evidence(root, manifest_path, _sha256(expected_sha256))
    return root, manifest_path, bundle


def _public_reference(bundle, manifest_path, origin):
    matches = [record for record in bundle.manifest["records"] if record["origin_path"] == origin]
    if len(matches) != 1:
        raise ValueError("required scientific artifact is missing or ambiguous in the public export")
    record = matches[0]
    return {
        "path": (PurePosixPath(manifest_path).parent / record["public_path"]).as_posix(),
        "sha256": record["public_sha256"],
        "payload_json_pointer": "/payload" if record["format"] == "json_envelope" else None,
    }


def _portable_checkpoint(root, manifest_path, bundle):
    artifact = bundle.fetch(CHECKPOINT_PATH, original_sha256=CHECKPOINT_SHA256)
    frozen = artifact.payload
    dependencies = _verify_native_checkpoint(root, frozen)
    _same(frozen["training_manifest"], bundle.fetch("outputs/native_training_run_02/training_manifest.json").payload,
          "portable checkpoint and training manifest")
    contract = bundle.fetch("outputs/native_training_run_02/contract.json").payload
    if set(contract) != {"root", "inputs", "code"}:
        raise ValueError("invalid portable training contract")
    for section in ("inputs", "code"):
        _same(sorted(contract[section], key=lambda row: row["path"]),
              [{key: entry[key] for key in ("path", "sha256", "role", "source")}
               for entry in frozen["training_manifest"][section]], "portable training contract dependencies")
    checkpoint_receipt = bundle.fetch("outputs/native_training_run_02/checkpoint_receipt.json").payload
    _same(checkpoint_receipt["path"], CHECKPOINT_PATH, "retained checkpoint receipt path")
    _same(checkpoint_receipt["sha256"], CHECKPOINT_SHA256, "retained checkpoint receipt identity")
    scope = artifact.integrity_scope
    if scope["public_content_verified"] is not True or scope["original_bytes_verified_on_load"] is not False:
        raise ValueError("portable checkpoint must distinguish public content from private original bytes")
    receipt = {
        "artifact_source": "portable", "mode": REPLAY_KIND,
        "checkpoint_sha256": CHECKPOINT_SHA256, "model_sha256": MODEL_SHA256,
        "checkpoint_reference": _public_reference(bundle, manifest_path, CHECKPOINT_PATH),
        "native_hog_fit_reference": _public_reference(bundle, manifest_path, "outputs/native_training_run_02/native_hog_fit.json"),
        "manifest_reference": {"path": manifest_path, "sha256": bundle.manifest_sha256},
        "scientific_content_preserved": True, "scientific_content_basis": scope["scientific_payload"],
        "transforms": scope["transforms"], "public_integrity_scope": bundle.manifest["integrity_scope"],
        "original_checkpoint_bytes": {
            "status": "not_checked", "expected_sha256": CHECKPOINT_SHA256,
            "actual_sha256": None, "matches": None,
            "reason": "Public mode does not open the private original checkpoint",
        },
        "all_relative_frozen_dependencies_verified": True,
        "input_count": len(frozen["training_manifest"]["inputs"]),
        "code_count": len(frozen["training_manifest"]["code"]),
        "dependencies": dependencies, "runtime": _runtime_identity(frozen),
        "serialization_timestamp_verified": False, "independent_biological_validation": False,
    }
    return frozen, receipt


def load_portable_native_checkpoint(root, *, manifest_path=None, expected_sha256=None):
    root, manifest_path, bundle = _load_bundle(root, manifest_path, expected_sha256)
    return _portable_checkpoint(root, manifest_path, bundle)


def _same(actual, expected, name):
    if _json_bytes(actual) != _json_bytes(expected):
        raise ValueError(f"{name} differs from the recorded scientific content")


def _relative_file(root, relative):
    if (not isinstance(relative, str) or not relative or "\\" in relative
            or ":" in relative or relative != PurePosixPath(relative).as_posix()):
        raise ValueError("replay inputs require canonical repository-relative paths")
    candidate = PurePosixPath(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError("replay inputs require canonical repository-relative paths")
    path = root / relative
    if path.resolve() != path or not path.is_file():
        raise ValueError(f"replay input must be an existing nonsymlink file: {relative}")
    return path


def _relative_json(root, relative):
    return _parse_json(_relative_file(root, relative).read_bytes())


def _verify_entries(root, entries, roles):
    if not isinstance(entries, list) or not entries:
        raise ValueError("replay requires a nonempty explicit dependency list")
    specs = []
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != _ENTRY_FIELDS:
            raise ValueError("invalid replay dependency fields")
        _relative_file(root, entry["path"])
        if entry["resolved_path"] != entry["path"]:
            raise ValueError("replay dependency aliases are not allowed")
        specs.append(InputSpec(**{key: entry[key] for key in ("path", "sha256", "role", "source")}))
    actual = _entries(root, specs, frozenset(roles))
    _same(actual, entries, "relative input/code/data identities")
    return actual


def _runtime_identity(frozen):
    expected = frozen["model"]["config"]["environment"]
    if set(expected) != {"python", "packages"}:
        raise ValueError("invalid recorded runtime identity")
    actual = {
        "python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in expected["packages"]},
    }
    _same(actual, expected, "runtime versions")
    loaded = []
    for entry in frozen["training_manifest"]["code"]:
        path = PurePosixPath(entry["path"])
        if path.parts[0] != "src" or path.suffix != ".py":
            continue
        parts = list(path.with_suffix("").parts[1:])
        if parts[-1] == "__init__":
            parts.pop()
        name = ".".join(parts)
        module = importlib.import_module(name)
        filename = getattr(module, "__file__", None)
        if filename is None:
            raise ValueError(f"running frozen module lacks source identity: {name}")
        source = Path(filename)
        if source.suffix == ".pyc":
            source = Path(importlib.util.source_from_cache(str(source)))
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        if digest != entry["sha256"]:
            raise ValueError(f"running frozen module differs from pinned code: {name}")
        loaded.append({"module": name, "path": entry["path"], "sha256": digest, "matches": True})
    return {
        "recorded_versions": expected, "actual_versions": actual, "versions_match": True,
        "loaded_frozen_modules": loaded,
        "scope": "Exact relative repository code/data bytes and source bytes of loaded frozen modules; trusted interpreter and installed distributions, not a runtime-binary attestation",
    }


def _native_checkpoint_manifest(frozen):
    if (not isinstance(frozen, dict) or set(frozen) != {"schema_version", "kind", "model", "training_manifest"}
            or frozen["kind"] != "native_checkpoint" or type(frozen["schema_version"]) is not int
            or frozen["schema_version"] != 1):
        raise ValueError("portable replay requires the retained native checkpoint schema")
    if hashlib.sha256(_json_bytes(frozen["model"])).hexdigest() != MODEL_SHA256:
        raise ValueError("native scientific model digest differs from the retained sealed model")
    manifest = frozen["training_manifest"]
    if (set(manifest) != {"schema_version", "kind", "root", "dependency_mode", "inputs", "code", "sha256"}
            or manifest["kind"] != "training_input_manifest" or manifest["dependency_mode"] != "explicit"):
        raise ValueError("invalid native training dependency manifest")
    _same({entry["path"]: entry["role"] for entry in manifest["inputs"]}, NATIVE_INPUT_ROLES,
          "native input roles")
    _same(sorted(entry["path"] for entry in manifest["code"]), sorted(FROZEN_CODE_PATHS),
          "frozen executable dependency roster")
    return manifest


def _verify_native_checkpoint(root, frozen):
    manifest = _native_checkpoint_manifest(frozen)
    inputs = _verify_entries(root, manifest["inputs"], {"fit", "prior", "selection", "provenance"})
    code = _verify_entries(root, manifest["code"], {"training_code"})
    identities = [(root / entry["path"]).stat() for entry in (*inputs, *code)]
    if len({(item.st_dev, item.st_ino) for item in identities}) != len(identities):
        raise ValueError("native source/code roles must identify distinct physical files")
    return [{"path": entry["path"], "role": entry["role"], "expected_sha256": entry["sha256"],
             "actual_sha256": entry["sha256"], "matches": True,
             "scope": "unchanged_relative_dependency_bytes"} for entry in (*inputs, *code)]


def restore_native_exchange_model(payload, native):
    from . import native_physiology

    rows = native_physiology._validated_rows(native.split()["train"])
    series = []
    for observable in native_physiology._UNITS:
        knots = payload.get("knots", {}).get(observable)
        if not isinstance(knots, list):
            raise ValueError("native exchange payload omits an observable")
        expected = [{
            "growth_rate_per_h": float(row.growth_rate_per_h),
            "value": float(row.value) if row.observation_status == "quantified" else None,
            "reported_value": float(row.reported_value),
            "observation_status": row.observation_status, "source_line": int(row.source_line),
        } for row in rows.loc[rows.observable_id.eq(observable)].itertuples(index=False)]
        _same(knots, expected, "native exchange source-bound knots")
        series.append(tuple(native_physiology._Knot(**knot) for knot in knots))
    model = native_physiology.NativeExchangeModel(tuple(series))
    _same(model.to_dict(), payload, "native exchange parameters, units and provenance")
    return model


def _native_state(root, frozen, recorded):
    from . import native_physiology
    from .hog_data import load_native_hog_data
    from .hog_learning import _native_parameter_groups, _native_training
    from ..mech.hog import HogModel

    config = frozen["model"]["config"]
    protocol = _relative_json(root, PROTOCOL_PATH)
    _same(config["protocol"], protocol, "native protocol")
    for key in ("prediction_kind", "prediction_scope", "quantity", "unit"):
        _same(config[key], protocol["prediction_contract"][key], f"native readout {key}")
    source = load_native_hog_data(root / "data/hog2013")
    selection = protocol["calibration_selection"]
    rows = source.observations
    selected = rows.loc[
        rows.supplement_id.eq(selection["supplement_id"])
        & rows.genotype.eq(selection["genotype"])
        & rows.strain_id.eq(selection["strain_id"])
        & rows.nacl_molar.isin(selection["nacl_molar"])
        & rows.observable_id.isin(selection["observable_ids"])
    ].copy().reset_index(drop=True)
    if len(selected) != selection["expected_rows"]:
        raise ValueError("native replay selection differs from the frozen protocol")
    selected["source_split"] = selected["split"]
    selected["split"] = "fit"
    selected = _native_training(selected)
    model = HogModel.from_source(root / "data/hog2013")
    fit = recorded("outputs/native_training_run_02/native_hog_fit.json")
    groups = _native_parameter_groups(protocol["candidate_parameter_groups"], model.parameter_values)
    names = sorted(name for members in groups.values() for name in members)
    _same(names, config["hog_override_keys"], "HOG override keys")
    _same(sorted(groups), sorted(protocol["parameter_selection"]["selected_groups"]), "HOG parameter groups")
    parameters = {name: frozen["model"]["parameters"][f"hog.{name}"] for name in names}
    gains = {name: frozen["model"]["parameters"][f"observation_gain.{name}"] for name in sorted(_WESTERNS)}
    _same(parameters, fit["parameters"], "HOG fitted parameters")
    _same(gains, fit["observation_gains"], "HOG fixed observation gains")
    expected = {f"hog.{name}": float(value) for name, value in {**model.parameter_values, **parameters}.items()}
    expected.update({f"observation_gain.{name}": value for name, value in gains.items()})
    _same(frozen["model"]["parameters"], expected, "complete HOG parameter payload")
    if any(_number(value, "fixed fitted parameter or gain") <= 0 for value in (*parameters.values(), *gains.values())):
        raise ValueError("HOG overrides and fixed observation gains must be positive")
    diagnostic = fit["diagnostics"]
    _same(diagnostic["parameter_groups"], protocol["candidate_parameter_groups"], "recorded parameter groups")
    _same(diagnostic["observation_gains"], gains, "recorded observation gains")
    _same(diagnostic["training_rows"], len(selected), "recorded training row count")
    _same(sorted(diagnostic["training_experiments"]), sorted(selected.experiment_id.unique()), "training experiments")
    for group, members in groups.items():
        for name in members:
            if not math.isclose(parameters[name] / model.parameter_values[name],
                                diagnostic["parameter_multipliers"][group], rel_tol=1e-12, abs_tol=0):
                raise ValueError("HOG grouped multipliers differ from the frozen fit")
    native = native_physiology.load_native_chemostat_data(root / protocol["native_exchange_fit"]["input"])
    exchange = restore_native_exchange_model(config["native_exchange_model"], native)
    _same(config["native_exchange_model"], recorded("outputs/native_training_run_02/native_exchange_fit.json"),
          "recorded native exchange fit")
    obligations = recorded("outputs/native_training_run_02/native_obligations.json")
    _same(obligations, config["native_obligation_scenarios"], "recorded native obligations")
    doses = protocol["native_obligation_summary"]["dose_support_molar"]
    _same(sorted(obligations), sorted(format(float(dose), ".17g") for dose in doses), "native dose support")
    for key, scenario in obligations.items():
        _same(key, format(_number(scenario["nacl_molar"], "native dose"), ".17g"), "native scenario dose")
        lower, upper = (_number(scenario[field], "native obligation") for field in ("minimum", "maximum"))
        if not 0 <= lower <= upper <= 1:
            raise ValueError("native obligation endpoints are outside their physical range")
        _same(scenario["balance_metadata"]["parameter_overrides"], parameters, "native balance parameters")
        _same(scenario["balance_metadata"]["source_model_sha256"], model.kinetic_model.metadata["sha256"],
              "native balance source identity")
    return {
        "root": root, "frozen": frozen, "config": config, "protocol": protocol,
        "source": source, "selected": selected, "model": model, "fit": fit,
        "parameters": parameters, "gains": gains, "native": native, "exchange": exchange,
        "obligations": obligations,
    }


def _near(actual, expected, name):
    actual = _number(actual, name)
    expected = _number(expected, name)
    if not math.isclose(actual, expected, rel_tol=_NUMERICAL_RTOL, abs_tol=_NUMERICAL_ATOL):
        raise ValueError(f"recomputed {name} differs from its recorded value")


def _replay_hog(state):
    import numpy as np

    from .hog_learning import _protocol

    frame, model = state["selected"], state["model"]
    times = frame.time_model_s.to_numpy(dtype=float)
    identifiers = frame.observable_id.to_numpy()
    raw = np.empty(len(frame))
    protocol_rows = {}
    for _, experiment in frame.groupby("experiment_id", sort=False):
        protocol_rows.setdefault(_protocol(experiment), []).extend(experiment.index.tolist())
    for protocol, indices in protocol_rows.items():
        indices = np.asarray(indices, dtype=int)
        grid = np.unique(np.r_[0.0, times[indices]])
        trajectory = model.simulate(grid, protocol=protocol, parameters=state["parameters"], rtol=1e-8, atol=1e-11)
        for name in np.unique(identifiers[indices]):
            selected = indices[identifiers[indices] == name]
            series = np.asarray(trajectory.variables[name], dtype=float)
            if series.shape != grid.shape or not np.isfinite(series).all():
                raise ValueError("HOG replay returned incomplete or nonfinite observations")
            raw[selected] = series[np.searchsorted(grid, times[selected])]
    predicted = raw.copy()
    values = frame.value.to_numpy(dtype=float)
    diagnostics = state["fit"]["diagnostics"]
    blocks = {name: group.index.to_numpy() for name, group in frame.groupby("observable_id", sort=False)}
    weights = {name: 1.0 / len(blocks) for name in blocks}
    _same(weights, diagnostics["observable_weights"], "native observable weights")
    scales, residual_weights, per_observable = {}, np.empty(len(frame)), {}
    for name, indices in blocks.items():
        if name in _WESTERNS:
            predicted[indices] *= state["gains"][name]
        scale = float(np.sqrt(np.mean(values[indices] ** 2)))
        _near(scale, diagnostics["normalization_scales"][name], f"{name} training scale")
        scales[name] = scale
        residual_weights[indices] = np.sqrt(weights[name] / len(indices)) / scale
        rmse = float(np.sqrt(np.mean((predicted[indices] - values[indices]) ** 2)))
        error = diagnostics["per_observable_errors"][name]
        _near(rmse, error["rmse"], f"{name} RMSE")
        _near(rmse / scale, error["nrmse"], f"{name} NRMSE")
        _same(error["unit"], "mol/L", "recorded HOG unit")
        _same(error["rows"], len(indices), "recorded HOG readout count")
        per_observable[name] = {"rmse": rmse, "nrmse": rmse / scale, "rows": len(indices), "unit": "mol/L"}
    nrmse = float(np.linalg.norm((predicted - values) * residual_weights))
    _near(nrmse, diagnostics["training_nrmse"], "HOG training NRMSE")
    rows = [{
        "experiment_id": row.experiment_id, "observable_id": row.observable_id,
        "source_sheet": row.source_sheet, "source_cell": row.source_cell,
        "time_model_s": float(row.time_model_s), "unit": row.unit,
        "observed": float(row.value), "prediction": float(prediction),
    } for row, prediction in zip(frame.itertuples(index=False), predicted, strict=True)]
    return {
        "verified": True, "rows": rows, "training_nrmse": nrmse,
        "per_observable": per_observable, "normalization_scales": scales,
        "fixed_observation_gains": state["gains"], "parameters_fitted": False,
        "observation_gains_fitted": False, "protocols_reexecuted": len(protocol_rows),
        "scope": "Previously fitted source-model predictions and training metrics, not independent validation",
    }


def _replay_obligations(state):
    import numpy as np

    from ..mech.hog import HogProtocol, hog_glycerol_balance

    window = state["protocol"]["native_obligation_summary"]
    times = np.unique(np.r_[np.arange(0.0, window["window_end_s"] + 0.1, window["time_grid_step_s"]), 3605.0])
    result = {}
    for dose in window["dose_support_molar"]:
        key = format(float(dose), ".17g")
        recorded = state["obligations"][key]
        trajectory = state["model"].simulate(times, protocol=HogProtocol(float(dose)), parameters=state["parameters"])
        balance = hog_glycerol_balance(state["model"], trajectory)
        _same(balance.metadata, recorded["balance_metadata"], "native balance metadata, units and roles")
        values = balance.carbon_commitment()[times >= window["window_start_s"]]
        replayed = {"minimum": float(values.min()), "maximum": float(values.max()), "mean": float(values.mean())}
        for name, value in replayed.items():
            _near(value, recorded[name], f"native obligation {key} {name}")
        result[key] = {**replayed, "nacl_molar": float(dose), "verified": True}
    return {"scenarios": result, "scope": window["summary"], "parameters_fitted": False}


def recompute_interval_metrics(sealed, observed):
    if not isinstance(sealed, Mapping) or not isinstance(observed, Mapping):
        raise ValueError("interval replay requires explicit prediction and observation mappings")
    predictions = _interval_rows(sealed["predictions"], set(observed))
    _same(predictions, sealed["predictions"], "canonical sealed interval rows")
    rows, widths, violations = [], [], []
    status_counts = dict.fromkeys(("ok", "failed", "infeasible", "unsupported"), 0)
    for prediction in predictions:
        label = observed[prediction["experiment_id"]]
        if not isinstance(label, dict) or set(label) != {"value", "source_path"}:
            raise ValueError("each replayed observation needs its value and distinct source path")
        value = _number(label["value"], "recorded outcome")
        status_counts[prediction["status"]] += 1
        width = contained = violation = side = None
        if prediction["status"] == "ok":
            lower, upper = prediction["lower"], prediction["upper"]
            width = _number(upper - lower, "interval width")
            if value < lower:
                violation, side = _number(lower - value, "interval violation"), "below"
            elif value > upper:
                violation, side = _number(value - upper, "interval violation"), "above"
            else:
                violation = 0.0
            contained = side is None
            widths.append(width)
            violations.append(violation)
        rows.append({
            **prediction, "observed": value, "label_source_path": label["source_path"],
            "contained": contained, "width": width, "violation": violation, "violation_side": side,
        })
    total, scorable = len(rows), len(widths)
    contained = sum(row["contained"] is True for row in rows)
    return {
        "quantity": sealed["task_manifest"]["quantity"], "unit": sealed["task_manifest"]["unit"],
        "prediction_kind": sealed["prediction_kind"], "prediction_scope": sealed["prediction_scope"],
        "assessment": "feasible_envelope_consistency_not_point_forecast", "rows": rows,
        "complete": scorable == total, "status_counts": status_counts,
        "n_total": total, "n_successful": status_counts["ok"], "n_scorable": scorable,
        "successful_fraction_total": status_counts["ok"] / total, "scorable_fraction_total": scorable / total,
        "n_contained": contained, "n_violated": scorable - contained,
        "containment_on_scorable": contained / scorable if scorable else None,
        "containment_fraction_total": contained / total,
        "mean_width_on_scorable": math.fsum(value / scorable for value in widths) if scorable else None,
        "mean_violation_on_scorable": math.fsum(value / scorable for value in violations) if scorable else None,
        "max_violation_on_scorable": max(violations, default=None),
    }


def _transfer_score(root, frozen, variant, recorded, original_digest):
    from .blind_transfer import BalancedCondition, HostFluxRoles
    from ..fba.chemical_task import ChemicalTask
    from ..fba.native_obligations import NativeObligationRoles

    folder = f"outputs/heldout_transfer_{variant}"
    sealed_path, task_path = f"{folder}/predictions.json", f"{folder}/task_manifest.json"
    score_path = f"outputs/heldout_score_{variant}/score.json"
    sealed, task, score = (recorded(path) for path in (sealed_path, task_path, score_path))
    if (set(sealed) != {"schema_version", "kind", "checkpoint_sha256", "model_sha256", "task_manifest",
                       "predictions", "prediction_kind", "prediction_scope"}
            or sealed["kind"] != "sealed_interval_predictions" or sealed["schema_version"] != 1):
        raise ValueError("invalid recorded interval prediction schema")
    for artifact in (sealed, task, score):
        _same(artifact["checkpoint_sha256"], CHECKPOINT_SHA256, "original checkpoint role link")
    _same(sealed["model_sha256"], MODEL_SHA256, "sealed scientific model link")
    _same(sealed["task_manifest"], task, "sealed task manifest content")
    expected_sha = _sha256(original_digest(sealed_path))
    _same(score["predictions_sha256"], expected_sha, "score-to-prediction role link")
    config = frozen["model"]["config"]
    for field in ("prediction_kind", "prediction_scope"):
        _same(sealed[field], config[field], f"sealed {field}")
    for field in ("quantity", "unit"):
        _same(task[field], config[field], f"task {field}")
    if (set(task) != {"schema_version", "kind", "root", "checkpoint_sha256", "inputs",
                     "experiment_inputs", "quantity", "unit", "sha256"}
            or task["kind"] != "evaluation_task_manifest" or task["schema_version"] != 1):
        raise ValueError("invalid recorded task schema")
    entries = _verify_entries(root, task["inputs"], {"task_structure", "task_conditions", "evaluation_code"})
    code = [{**entry, "role": "evaluation_code", "source": {**entry["source"], "scope": "evaluation_code"}}
            for entry in frozen["training_manifest"]["code"]]
    _same([entry for entry in entries if entry["role"] == "evaluation_code"], code,
          "pre-frozen evaluation code identities")
    chemistry_path = f"data/holdout_transfer/hypotheses/chemistry_{variant}.json"
    conditions_path = f"data/holdout_transfer/hypotheses/conditions_{variant}.json"
    task_entries = {entry["path"]: entry for entry in entries if entry["role"] != "evaluation_code"}
    _same({path: entry["role"] for path, entry in task_entries.items()},
          {chemistry_path: "task_structure", conditions_path: "task_conditions"}, "distinct task source roles")
    chemistry, conditions = (_relative_json(root, path) for path in (chemistry_path, conditions_path))
    for path, artifact, fields in ((chemistry_path, chemistry, {"schema_version", "source", "task"}),
                                   (conditions_path, conditions, {"schema_version", "source", "experiments"})):
        if set(artifact) != fields or artifact["schema_version"] != 1:
            raise ValueError("invalid public task input schema")
        _same(artifact["source"], task_entries[path]["source"], "structured task source identity")
    ChemicalTask.from_dict(chemistry["task"])
    role_values = dict(config["host_roles"])
    role_values["native"] = NativeObligationRoles(**role_values["native"])
    HostFluxRoles(**role_values)
    experiments = _unique_rows(conditions["experiments"], {
        "experiment_id", "chemistry_path", "chemistry_sha256", "condition", "added_nacl_molar",
    }, "portable task conditions")
    expected_inputs = []
    for identifier, row in sorted(experiments.items()):
        _same(row["chemistry_path"], chemistry_path, "experiment chemistry path")
        _same(row["chemistry_sha256"], task_entries[chemistry_path]["sha256"], "experiment chemistry identity")
        condition = BalancedCondition(**row["condition"])
        _same(asdict(condition), row["condition"], "condition payload")
        _same(condition.experiment_id, identifier, "condition experiment identity")
        _same(condition.readout_basis, config["quantity"], "condition readout basis")
        if _number(condition.growth_rate_per_h, "imposed growth rate") <= 0:
            raise ValueError("condition growth must be positive")
        if _number(row["added_nacl_molar"], "added NaCl") < 0:
            raise ValueError("condition NaCl must be nonnegative")
        expected_inputs.append({"experiment_id": identifier, "input_paths": sorted([chemistry_path, conditions_path])})
    _same(task["experiment_inputs"], expected_inputs, "complete experiment-to-input role links")
    _interval_rows(sealed["predictions"], set(experiments))
    label_manifest = score["label_manifest"]
    if (set(label_manifest) != {"schema_version", "kind", "root", "inputs", "sha256"}
            or label_manifest["kind"] != "outcome_labels_manifest" or label_manifest["schema_version"] != 1):
        raise ValueError("invalid recorded outcome-label schema")
    label_entries = _verify_entries(root, label_manifest["inputs"], {"outcome_labels"})
    _reject_label_overlap(label_entries, [*frozen["training_manifest"]["inputs"], *code, *entries])
    label_specs = recorded("data/holdout_transfer/labels_manifest.json")
    _same(label_specs, {"labels": [{key: entry[key] for key in ("path", "sha256", "role", "source")}
                                   for entry in label_entries]}, "label source manifest")
    observed = {}
    for entry in label_entries:
        labels = _relative_json(root, entry["path"])
        if set(labels) != {"quantity", "unit", "observations"}:
            raise ValueError("invalid outcome data schema")
        for field in ("quantity", "unit"):
            _same(labels[field], config[field], f"label {field}")
        for identifier, row in _unique_rows(labels["observations"], {"experiment_id", "value"}, "outcome labels").items():
            if identifier in observed:
                raise ValueError("duplicate observed experiment across label files")
            observed[identifier] = {"value": _number(row["value"], "outcome value"), "source_path": entry["path"]}
    metrics = recompute_interval_metrics(sealed, observed)
    identity_fields = {"schema_version", "kind", "checkpoint_sha256", "predictions_sha256", "label_manifest", "sha256"}
    if (set(score) != identity_fields | set(metrics) or score["kind"] != "interval_prediction_score"
            or score["schema_version"] != 1):
        raise ValueError("invalid recorded interval score schema")
    _same({key: value for key, value in score.items() if key not in identity_fields}, metrics,
          "recomputed interval score rows and metrics")
    return {
        "variant": variant, "score_verified": True, "computed_metrics": metrics,
        "sealed_prediction_content_preserved": True, "prediction_kernel_reexecuted": False,
        "original_predictions_sha256": expected_sha, "recorded_score_digest": score["sha256"],
        "scope": "Arithmetic verification of all recorded interval predictions against the pinned outcome labels; no new LP solve, fit, hypothesis selection or validation evidence",
    }


def _replay_population_metrics(bundle):
    import numpy as np

    from .portable_evidence import POPULATION_FREEZE_PATH, POPULATION_FREEZE_SHA256

    freeze = bundle.fetch(POPULATION_FREEZE_PATH, original_sha256=POPULATION_FREEZE_SHA256).payload
    run = "outputs/native_population_development/native_training_01_phase_c_01_scoring"
    report = bundle.fetch(f"{run}/development_report.json").payload
    summary = bundle.fetch(f"{run}/scoring_summary.json").payload
    protocol = bundle.fetch("data/native_law_v2/development_protocol.json").payload
    _same(report["prediction_freeze_sha256"], POPULATION_FREEZE_SHA256, "population report-to-prediction link")
    _same(summary["prediction_freeze"]["sha256"], POPULATION_FREEZE_SHA256, "population summary-to-prediction link")
    _same(report["protocol_id"], protocol["protocol_id"], "population report protocol")
    _same(freeze["protocol_id"], protocol["protocol_id"], "population freeze protocol")
    for artifact in (freeze, report, summary):
        _same(artifact["strong_claim_authorized"], False, "population strong-claim boundary")
    for artifact in (report, summary):
        _same(artifact["final_test_defined"], False, "population final-test boundary")
    for field in ("H_train", "S_train"):
        _same(summary[field], freeze[field], "population frozen training scale")
    scale = _number(freeze["S_train"], "population frozen training scale")
    if scale <= 0:
        raise ValueError("this recorded population replay requires its positive frozen diagnostic scale")
    models = [row["model_id"] for row in protocol["model_recipe"]["canonical_families"]]
    groups = [row["group_id"] for row in protocol["partition"]["development"]]
    start, end = protocol["response_coordinates"]["expected_absolute_columns_1based"]
    frames = list(range(start, end + 1))
    _same(len(frames), protocol["response_coordinates"]["frame_count"], "population frame inventory")
    expected = {(model, group, frame) for model in models for group in groups for frame in frames}
    predictions = {}
    for row in freeze["predictions"]:
        key = row["model_id"], row["group_id"], row["frame_1based"]
        if key in predictions or key not in expected:
            raise ValueError("duplicate or unexpected sealed population prediction")
        _same(row["status"], "ok", "recorded population prediction status")
        _same(row["reason"], None, "recorded population prediction reason")
        for field in ("quantity", "unit"):
            _same(row[field], protocol["estimand"][field], f"population {field}")
        _number(row["mean"], "sealed population mean")
        _number(row["time_min"], "sealed population time")
        predictions[key] = row
    if set(predictions) != expected:
        raise ValueError("sealed population predictions omit registered frames or models")
    seen, observations, computed = set(), {}, {model: {} for model in models}
    for group in report["group_evaluations"]:
        model, identifier = group["model_id"], group["group_id"]
        if model not in computed or identifier not in groups or identifier in computed[model]:
            raise ValueError("duplicate or unexpected population evaluation group")
        _same([row["frame_1based"] for row in group["frames"]], frames, "population evaluated frame inventory")
        errors = []
        for row in group["frames"]:
            key = row["model_id"], row["group_id"], row["frame_1based"]
            if key in seen or key[:2] != (model, identifier) or key not in predictions:
                raise ValueError("population score row has an invalid prediction role link")
            seen.add(key)
            prediction = predictions[key]
            _same(row["status"], "quantified", "population observation status")
            _same(row["predicted_mean"], prediction["mean"], "population sealed mean")
            _same(row["time_min"], prediction["time_min"], "population fixed time coordinate")
            observed = _number(row["observed_mean"], "recorded aggregate observation")
            coordinate = identifier, row["frame_1based"]
            observation = {field: row[field] for field in ("observed_mean", "time_min", "prefix_cells", "valid_cells",
                                                          "source_missing_cells", "invalid_cells", "cell_sample_sd")}
            if coordinate in observations:
                _same(observation, observations[coordinate], "population aggregate shared across model families")
            observations[coordinate] = observation
            error = _number(prediction["mean"] - observed, "population residual")
            _near(row["squared_error"], error * error, "population squared error")
            _near(row["absolute_error"], abs(error), "population absolute error")
            errors.append(error)
        mse = float(np.mean(np.square(errors)))
        metrics = {"complete": True, "mean_mse": mse, "rmse": math.sqrt(mse),
                   "mae": float(np.mean(np.abs(errors)))}
        _same(group["complete"], True, "population complete-frame coverage")
        for field in ("mean_mse", "rmse", "mae"):
            _near(group[field], metrics[field], f"population group {field}")
        computed[model][identifier] = metrics
    if seen != expected or set(summary["models"]) != set(models):
        raise ValueError("population metric replay omitted a frozen model or frame")
    model_metrics = {}
    for model, per_group in computed.items():
        if set(per_group) != set(groups):
            raise ValueError("population metric replay omitted a registered biological group")
        mse = float(np.mean([per_group[group]["mean_mse"] for group in groups]))
        mae = float(np.mean([per_group[group]["mae"] for group in groups]))
        metrics = {"complete": True, "groups": per_group, "mean_mse": mse, "rmse": math.sqrt(mse),
                   "mae": mae, "scaled_rmse": math.sqrt(mse) / scale, "scaled_mae": mae / scale}
        previous = summary["models"][model]
        _same(sorted(previous), sorted(metrics), "population summary metric fields")
        _same(previous["complete"], True, "population model coverage")
        _same(previous["groups"], per_group, "population summary group metrics")
        for field in ("mean_mse", "rmse", "mae", "scaled_rmse", "scaled_mae"):
            _near(previous[field], metrics[field], f"population model {field}")
        model_metrics[model] = metrics
    _same(summary["prediction_count"], len(predictions), "population prediction count")
    _same(summary["evaluated_frame_rows"], len(seen), "population evaluated frame count")
    _same(summary["selected_development_model"], report["selected_development_model"], "recorded population selection label")
    return {
        "verified": True, "prediction_count": len(predictions), "evaluated_frame_rows": len(seen),
        "quantity": protocol["estimand"]["quantity"], "unit": protocol["estimand"]["unit"],
        "models": model_metrics, "recorded_selected_development_model": summary["selected_development_model"],
        "prediction_kernel_reexecuted": False, "selection_reexecuted": False,
        "raw_source_verification": "not_checked", "strong_claim_authorized": False, "final_test_defined": False,
        "scope": "Arithmetic verification of sealed population means against the previously reported aggregates, not private raw-source reanalysis, a new model selection or a fresh test",
    }


def _check_private_originals(bundle, root):
    from .portable_evidence import _pointer, _set_pointer

    root = Path(root).resolve()
    if root.exists() and not root.is_dir():
        raise ValueError("the caller-supplied private root must be a directory")
    checks = []
    for record in bundle.manifest["records"]:
        origin = record["origin_path"]
        path = root / origin
        result = {"path": origin, "expected_sha256": record["original_sha256"],
                  "actual_sha256": None, "matches": None, "status": "not_checked",
                  "transform_equivalence_verified": False}
        if not path.exists() and not path.is_symlink():
            checks.append({**result, "reason": "original file not present at the caller-supplied private root"})
            continue
        original = _relative_file(root, origin).read_bytes()
        actual = hashlib.sha256(original).hexdigest()
        if actual != record["original_sha256"]:
            raise ValueError(f"available original byte digest mismatch: {origin}")
        public = bundle.fetch(origin).payload
        if record["format"] == "json_envelope":
            projected = _parse_json(original)
            for transform in record["transforms"]:
                value = _pointer(projected, transform["json_pointer"])
                _same(hashlib.sha256(_json_bytes(value)).hexdigest(), transform["original_value_sha256"],
                      "original storage value digest")
                _set_pointer(projected, transform["json_pointer"], transform["replacement"])
            _same(projected, public, "original-to-public declared transform equivalence")
        elif original != public:
            raise ValueError("original-byte companion differs from public content")
        checks.append({**result, "status": "verified", "actual_sha256": actual, "matches": True,
                       "transform_equivalence_verified": True})
    return checks


def replay_portable_evidence(root, *, manifest_path=None, expected_sha256=None, private_original_root=None):
    from . import native_physiology

    root, manifest_path, bundle = _load_bundle(root, manifest_path, expected_sha256)
    frozen, integrity = _portable_checkpoint(root, manifest_path, bundle)

    def recorded(path):
        return bundle.fetch(path).payload

    def original_digest(path):
        return bundle.fetch(path).integrity_scope["original_sha256"]

    state = _native_state(root, frozen, recorded)
    hog = _replay_hog(state)
    obligations = _replay_obligations(state)
    exchange = native_physiology.evaluate_native_exchange_model(state["exchange"], state["native"].split()["holdout"])
    _same(exchange, recorded("outputs/native_training_run_02/native_exchange_validation.json"),
          "native exchange predictions and score")
    hypotheses = recorded("data/holdout_transfer/hypotheses/manifest.json")
    _same(hypotheses["registration"]["declared_variants_in_order"], list(_VARIANTS), "complete recorded hypothesis set")
    scores = [_transfer_score(root, frozen, variant, recorded, original_digest) for variant in _VARIANTS]
    population = _replay_population_metrics(bundle)
    private_checks = None
    if private_original_root is not None:
        private_checks = _check_private_originals(bundle, private_original_root)
        integrity["original_checkpoint_bytes"] = next(row for row in private_checks if row["path"] == CHECKPOINT_PATH)
    _load_bundle(root, manifest_path, bundle.manifest_sha256)
    _same(_verify_native_checkpoint(root, frozen), integrity["dependencies"], "post-replay relative dependency identity")
    _same(_runtime_identity(frozen), integrity["runtime"], "post-replay runtime identity")
    return {
        "schema_version": 1, "kind": REPLAY_KIND, "verified": True,
        "integrity": integrity, "private_original_checks": private_checks,
        "public_exports_and_relative_dependencies_unchanged": True,
        "native_hog": hog, "native_obligations": obligations,
        "native_exchange": {"verified": True, "computed_report": exchange, "parameters_fitted": False},
        "transfer_scores": scores, "population_metrics": population,
        "new_fitting": False, "new_selection": False, "independent_test": False,
        "biological_validation": False, "original_serialization_timestamp_verified": False,
        "limitations": [
            "The public manifest authenticates the exported content and declared transform projection, not unavailable private original bytes",
            "Native predictions and metrics reuse the original learned parameters, fixed gains and source kernels; historical optimizer and sensitivity diagnostics are not rerun",
            "Transfer interval score arithmetic is recomputed without rerunning the LP envelope solver or selecting among the three hypotheses",
            "No original serialization timestamp, historical operator blindness, new biological evidence or stronger biology grade is established",
            "Lineage-only external references and private population/raw-source replays remain outside this replay scope",
        ],
    }
