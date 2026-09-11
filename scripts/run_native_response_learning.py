from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import platform
import sys
from dataclasses import asdict, fields
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

np = importlib.import_module("numpy")
learning = importlib.import_module("ystwin.analysis.native_response_learning")


def interface_contract():
    recipe = learning.SelectionRecipe()
    return {
        "status": "synthetic_engineering_interface_native_projection_adapter_pending",
        "default_recipe_status": "engineering_example_not_registered_for_native_data",
        "native_custody_policy": "Only exact custodian-approved train/development projections; current Granados custody schema remains blocked pending verified adapter and independent admission.",
        "readout_fields": [field.name for field in fields(learning.Readout)],
        "history_fields": [field.name for field in fields(learning.PhysicalHistory)],
        "experiment_fields": [field.name for field in fields(learning.NativeExperiment)],
        "values_axes": "cell_by_time; null or explicit false mask means missing, not zero",
        "prediction_target": "conditional_population_mean; score each registered cell-time value within equal-weight whole experiments",
        "default_recipe": asdict(recipe),
        "candidate_library": learning.candidate_library(recipe),
        "equations": learning.EQUATIONS,
        "time_units": learning.TIME_UNITS,
        "dose_units": learning.DOSE_UNITS,
        "verifier_gaps": learning.VERIFIER_GAPS,
        "freeze_authority": "Canonical SHA256 and exclusive file creation provide integrity only, not independent prior registration or scientific authorization.",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Bounded native-response learner: metadata gate by default, no reserved measurements or raw acquisition.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/native_law_v2/granados_manifest.json")
    parser.add_argument("--protocol", type=Path, default=ROOT / "data/native_law_v2/protocol.json")
    parser.add_argument("--contract", action="store_true")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--freeze", type=Path)
    args = parser.parse_args(argv)
    if args.contract:
        if args.fit or args.freeze is not None:
            parser.error("--contract cannot be combined with --fit or --freeze")
        print(json.dumps(interface_contract(), indent=2, sort_keys=True, allow_nan=False))
        return 0
    if args.fit != (args.freeze is not None):
        parser.error("--fit and a new --freeze path must be supplied together")
    gate = learning.check_data_gate(args.manifest, args.protocol)
    gate["scientific_authorization"] = "none"
    if not gate["ready"] or not args.fit:
        print(json.dumps(gate, indent=2, sort_keys=True, allow_nan=False))
        return 0 if gate["ready"] else 2
    opened = False
    try:
        if args.freeze.exists() or args.freeze.is_symlink():
            raise FileExistsError("freeze already exists; overwriting is forbidden")
        if not args.freeze.parent.is_dir():
            raise learning.ContractError("freeze parent directory must already exist")
        opened = True
        bundle = learning.load_verified_train_development(args.manifest, args.protocol)
        result = learning.fit_native_response(bundle["train"], bundle["development"], recipe=bundle["recipe"])
        result["data_status"] = "projection_integrity_checked_external_attestations_not_independently_authenticated"
        result["provenance"] = {
            **bundle["provenance"],
            "learner_source_sha256": hashlib.sha256(Path(learning.__file__).read_bytes()).hexdigest(),
            "runner_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
        }
        envelope = learning.write_freeze(result, args.freeze)
        print(json.dumps({
            "status": "training_development_frozen_not_scientifically_authorized",
            "freeze": str(args.freeze),
            "result_sha256": envelope["result_sha256"],
            "model_sha256": envelope["model_sha256"],
            "predictions_sha256": envelope["predictions_sha256"],
            "selected_family": result["selected_family"],
            "candidate_attempt_count": len(result["candidates"]),
            "failures": result["failures"],
            "verifier_gaps": result["verifier_gaps"],
            "final_test_accessed": False,
        }, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except (learning.ContractError, OSError) as error:
        print(json.dumps({
            "status": "refused",
            "failure": str(error),
            "measurement_access_attempted": opened,
            "final_test_accessed": False,
            "attempts": getattr(error, "attempts", []),
        }, indent=2, sort_keys=True, allow_nan=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
