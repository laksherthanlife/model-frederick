from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import sys
import sysconfig

import numpy as np

from ystwin import paths
from ystwin.analysis.hog_data import load_native_hog_data
from ystwin.analysis.hog_learning import fit_native_hog_parameters
from ystwin.analysis.native_physiology import evaluate_native_exchange_model, fit_native_exchange_model, load_native_chemostat_data
from ystwin.analysis.training_freeze import InputSpec, TrainingContract, build_training_manifest, freeze_checkpoint
from ystwin.mech.hog import HogModel, HogProtocol, hog_glycerol_balance


TRAINING_CODE = (
    "pyproject.toml", "scripts/run_native_training.py", "src/ystwin/__init__.py", "src/ystwin/paths.py",
    "src/ystwin/analysis/__init__.py", "src/ystwin/analysis/hog_data.py",
    "src/ystwin/analysis/hog_learning.py", "src/ystwin/analysis/native_physiology.py",
    "src/ystwin/analysis/training_freeze.py", "src/ystwin/mech/__init__.py",
    "src/ystwin/mech/hog.py", "src/ystwin/mech/kinetic_sbml.py",
)
PREDICTION_CODE = (
    "scripts/evaluate_frozen_transfer.py", "src/ystwin/analysis/blind_transfer.py", "src/ystwin/fba/__init__.py",
    "src/ystwin/fba/chemical_task.py", "src/ystwin/fba/native_obligations.py",
    "src/ystwin/fba/dynamic_rates.py", "src/ystwin/fba/physiology.py", "src/ystwin/fba/solver.py",
)
PROTOCOL_PATH = "data/hog2013/native_training_protocol.json"


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source(identifier, *, code=False):
    return {"authority": "explicit native training protocol", "source_id": identifier,
            "reference": PROTOCOL_PATH,
            "scope": "native_training_code" if code else "native_physiology"}


def build_contract(root, protocol):
    role_map = {
        "published_native_reaction_structure": "prior",
        "native_training_and_native_control_measurements": "fit",
        "measurement_units_and_source_assumptions": "provenance",
        "source_identity_and_checksums": "provenance",
        "native_growth_exchange_training_and_condition_validation": "fit",
        "native_host_prior_frozen_for_generic_inference": "prior",
    }
    inputs = [InputSpec(entry["path"], _digest(root / entry["path"]), role_map[entry["role"]],
                        _source(entry["role"])) for entry in protocol["training_inputs"]]
    inputs.append(InputSpec(PROTOCOL_PATH, _digest(root / PROTOCOL_PATH), "selection", _source("native_selection")))
    code = [InputSpec(name, _digest(root / name), "training_code", _source(name, code=True))
            for name in (*TRAINING_CODE, *PREDICTION_CODE)]
    return TrainingContract(root, tuple(inputs), tuple(code))


def install_project_read_guard(root, allowed_paths, output):
    root, output = root.resolve(), output.resolve()
    allowed = {Path(path).resolve() for path in allowed_paths}
    environment_roots = {Path(sysconfig.get_path(name)).resolve()
                         for name in ("purelib", "platlib", "stdlib", "platstdlib")}
    audit = {"phase": "native_modeling", "reads": {}, "scope":
             "Python audit events for project-local file reads/writes; interpreter libraries are separate, not an OS sandbox or a claim about arbitrary native-code I/O"}

    def hook(event, arguments):
        if event != "open" or not arguments or isinstance(arguments[0], int):
            return
        logical = Path(os.fsdecode(arguments[0])).absolute()
        path = logical.resolve()
        if (not logical.is_relative_to(root) and not path.is_relative_to(root)
                or path.is_relative_to(output)
                or any(path.is_relative_to(folder) for folder in environment_roots)):
            return
        flags = arguments[2] if len(arguments) > 2 else 0
        writable = isinstance(flags, int) and (flags & os.O_ACCMODE) != os.O_RDONLY
        if writable:
            raise PermissionError("native training cannot modify project source inputs")
        canonical = path
        if path.suffix == ".pyc":
            try:
                canonical = Path(importlib.util.source_from_cache(str(path))).resolve()
            except ValueError:
                pass
        if canonical not in allowed:
            display = path.relative_to(root) if path.is_relative_to(root) else logical.relative_to(root)
            raise PermissionError(f"project input is outside the explicit native allowlist: {display}")
        key = f"{audit['phase']}:{canonical.relative_to(root)}"
        audit["reads"][key] = audit["reads"].get(key, 0) + 1

    sys.addaudithook(hook)
    return audit


def _json(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def _loaded_project_modules(root, allowed):
    loaded = set()
    for module in tuple(sys.modules.values()):
        name = getattr(module, "__name__", "")
        if not (name == "__main__" or name == "ystwin" or name.startswith(("ystwin.", "scripts."))):
            continue
        filename = getattr(module, "__file__", None)
        if filename is None:
            continue
        path = Path(filename).resolve()
        if path.is_relative_to(root):
            if path.suffix == ".pyc":
                path = Path(importlib.util.source_from_cache(str(path))).resolve()
            loaded.add(path)
    extra = loaded - allowed
    if extra:
        raise RuntimeError(f"undeclared project modules entered native training: {sorted(str(p.relative_to(root)) for p in extra)}")
    return sorted(str(path.relative_to(root)) for path in loaded)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Fit native physiology and freeze a product-agnostic model before any downstream task")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-nfev", type=int, default=30)
    args = parser.parse_args(argv)
    root = paths.REPO_ROOT.resolve()
    output = args.output_dir.resolve()
    if not output.is_relative_to(root):
        raise ValueError("the audited run output must be inside this repository")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("use a new or empty native-training destination")
    protocol = json.loads((root / PROTOCOL_PATH).read_text())
    contract = build_contract(root, protocol)
    manifest = build_training_manifest(contract)
    output.mkdir(parents=True, exist_ok=True)
    _json(output / "contract.json", {"root": str(root), "inputs": [asdict(item) for item in contract.inputs],
                                     "code": [asdict(item) for item in contract.code]})
    _json(output / "training_manifest.json", manifest)
    allowed = {(root / item.path).resolve() for item in (*contract.inputs, *contract.code)}
    read_audit = install_project_read_guard(root, allowed, output)
    source = load_native_hog_data(root / "data/hog2013")
    selection = protocol["calibration_selection"]
    rows = source.observations
    selected = rows.loc[
        rows.supplement_id.eq(selection["supplement_id"])
        & rows.genotype.eq(selection["genotype"])
        & rows.strain_id.eq(selection["strain_id"])
        & rows.nacl_molar.isin(selection["nacl_molar"])
        & rows.observable_id.isin(selection["observable_ids"])
    ].copy()
    if len(selected) != selection["expected_rows"]:
        raise ValueError("native observations do not match the frozen selection protocol")
    selected["source_split"] = selected["split"]
    selected["split"] = "fit"
    selected.to_csv(output / "native_training_observations.csv", index=False)
    model = HogModel.from_source(root / "data/hog2013")
    print("Fitting native measured protein/glycerol trajectories; no downstream task inputs are available.", flush=True)
    fit = fit_native_hog_parameters(
        model, selected, parameter_groups=protocol["candidate_parameter_groups"],
        fold_bounds=tuple(protocol["fitting_rules"]["search_fold_bounds"]), max_nfev=args.max_nfev)
    if not fit.diagnostics["optimization_converged"]:
        raise RuntimeError("native fitting did not converge; no checkpoint was frozen")
    _json(output / "native_hog_fit.json", {"parameters": fit.parameter_values,
                                          "observation_gains": fit.observation_gains,
                                          "diagnostics": fit.diagnostics})
    native_parts = load_native_chemostat_data(root / protocol["native_exchange_fit"]["input"]).split()
    exchange = fit_native_exchange_model(native_parts["train"])
    exchange_payload = exchange.to_dict()
    exchange_report = evaluate_native_exchange_model(exchange, native_parts["holdout"])
    _json(output / "native_exchange_fit.json", exchange_payload)
    _json(output / "native_exchange_validation.json", exchange_report)
    window = protocol["native_obligation_summary"]
    times = np.unique(np.r_[np.arange(0.0, window["window_end_s"] + 0.1, window["time_grid_step_s"]), 3605.0])
    obligations = {}
    for dose in window["dose_support_molar"]:
        trajectory = model.simulate(times, protocol=HogProtocol(float(dose)), parameters=fit.parameter_values)
        balance = hog_glycerol_balance(model, trajectory)
        fractions = balance.carbon_commitment()
        values = fractions[times >= window["window_start_s"]]
        obligations[format(float(dose), ".17g")] = {
            "nacl_molar": float(dose), "minimum": float(values.min()), "maximum": float(values.max()),
            "mean": float(values.mean()), "scope": window["summary"],
            "source": "native HOG reaction-balance fit; concentration scale is not a gDW flux calibration",
            "source_medium": "native W303 YPD; transfer to a glucose-only host scenario is an explicit assumption",
            "balance_metadata": balance.metadata,
        }
    _json(output / "native_obligations.json", obligations)
    parameters = {f"hog.{key}": float(value) for key, value in {**model.parameter_values, **fit.parameter_values}.items()}
    parameters.update({f"observation_gain.{key}": float(value) for key, value in fit.observation_gains.items()})
    config = {
        "protocol": protocol, "hog_override_keys": sorted(fit.parameter_values),
        "native_exchange_model": exchange_payload, "native_obligation_scenarios": obligations,
        "prediction_kind": protocol["prediction_contract"]["prediction_kind"],
        "prediction_scope": protocol["prediction_contract"]["prediction_scope"],
        "quantity": protocol["prediction_contract"]["quantity"], "unit": protocol["prediction_contract"]["unit"],
        "host_model_path": "data/gem/ecYeastGEM_batch.xml.gz",
        "host_roles": {"biomass_reaction_id": "r_2111",
                       "oxygen_exchange_ids": ["r_1992", "r_1992_REV"],
                       "intracellular_compartment_ids": ["ce", "c", "m", "n", "p", "er", "g", "lp", "v", "erm", "vm", "gm", "mm"],
                       "extracellular_compartment_ids": ["e"],
                       "native": {"native_metabolite_id": "s_0765[c]", "glucose_metabolite_id": "s_0565[e]",
                                  "glucose_exchange_ids": ["r_1714", "r_1714_REV"],
                                  "native_export_ids": ["r_1808"], "resource_metabolite_ids": ["prot_pool[c]"]}},
        "host_medium_policy": "The native source already has glucose-only carbon supply. Preserve all existing bounds and impose additional stoichiometry-weighted native glucose/oxygen caps; no medium opening or hard-bound relaxation is authorized.",
        "environment": {"python": sys.version.split()[0],
                        "packages": {name: importlib.metadata.version(name) for name in
                                     ("numpy", "scipy", "pandas", "cobra", "python-libsbml", "xlrd")}},
        "max_nfev": args.max_nfev,
    }
    loaded = _loaded_project_modules(root, allowed)
    read_audit["phase"] = "freeze_integrity_check"
    checkpoint = freeze_checkpoint(output / "checkpoint.json", contract=contract,
                                   training_manifest=manifest, parameters=parameters, config=config)
    _json(output / "checkpoint_receipt.json", {"path": str(checkpoint.path.relative_to(root)),
                                              "sha256": checkpoint.sha256,
                                              "retention_rule": "retain this expected digest separately before introducing any downstream task"})
    _json(output / "read_audit.json", {**read_audit, "loaded_project_modules": loaded})
    print(json.dumps({"checkpoint": str(checkpoint.path.relative_to(root)), "sha256": checkpoint.sha256,
                      "native_training_nrmse": fit.training_nrmse,
                      "native_condition_validation_scored": exchange_report["n_scored"],
                      "native_condition_validation_total": exchange_report["n_readouts"]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
