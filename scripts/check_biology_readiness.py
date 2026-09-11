from __future__ import annotations

import argparse
import json
from pathlib import Path

from ystwin.analysis.biology_learning_audit import (
    CHECKPOINT_SHA256,
    collect_learning_evidence,
    require_claim,
)
from ystwin.fba.native_reconciliation import (
    condition_constraints,
    ec_native_observables,
    native_conditions,
    require_native_consistency,
    solve_native,
)
from ystwin.fba.solver import load_model


_MODEL_PATH = "data/gem/ecYeastGEM_batch.xml.gz"


def _condition_probe(model, condition, observables, mode, timeout_s):
    constraints = condition_constraints(condition, observables, mode=mode, impose_growth=True)
    row = {
        "condition_id": condition["condition_id"],
        "growth_rate_per_h": condition["growth_rate_per_h"],
        "constraint_mode": mode,
        "scope": "native observation consistency at imposed growth, not an independent growth forecast",
        "uncertainty_scope": ("printed rounding, not a confidence interval" if mode == "printed_rounding"
                              else "reported points, not a measurement uncertainty model"),
        "censored_unquantified_readouts": [name for name, value in condition["observations"].items()
                                           if name in observables and value["value"] is None],
        "quantified_constraints": len(constraints),
        "status": "failed", "consistent": False, "reason": None, "physical_audit": None,
    }
    try:
        result = solve_native(model, constraints=constraints, objective=observables["growth"], timeout_s=timeout_s)
        row["status"] = result["status"]
        row["physical_audit"] = result["audit"]["summary"] if result["audit"] is not None else None
        require_native_consistency(result, growth_key="growth", target_growth_per_h=condition["growth_rate_per_h"])
        row["consistent"] = True
    except (ValueError, RuntimeError) as exc:
        row["reason"] = str(exc)
    return row


def run_readiness(root, output, *, mode="printed_rounding", timeout_s=15):
    root, output = Path(root).resolve(), Path(output).resolve()
    if mode not in ("quantified_point", "printed_rounding"):
        raise ValueError("readiness requires an explicit native constraint mode")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("readiness output must be new or empty")
    ledger = collect_learning_evidence(root)
    require_claim(ledger, "integrity_established", root=root)
    model, settings = load_model(root / _MODEL_PATH)
    observables = ec_native_observables(model)
    rows = [_condition_probe(model, condition, observables, mode, timeout_s) for condition in native_conditions()]
    claims = []
    for claim in ("conditional_parameter_estimation", "independent_mechanistic_evidence"):
        try:
            assessment = require_claim(ledger, claim, root=root)
            claims.append({"claim": claim, "authorized": True, "assessment": assessment})
        except ValueError as exc:
            claims.append({"claim": claim, "authorized": False, "reason": str(exc)})
    current = collect_learning_evidence(root)
    if current != ledger:
        raise ValueError("native evidence changed during readiness evaluation")
    require_claim(current, "integrity_established", root=root)
    passed = sum(row["consistent"] for row in rows)
    ready = passed == len(rows) and all(claim["authorized"] for claim in claims)
    report = {
        "schema_version": 1,
        "assessment": "current_frozen_model_biological_learning_readiness",
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "source_model": _MODEL_PATH, "solver": str(settings),
        "constraint_mode": mode, "native_conditions": rows,
        "n_native_conditions": len(rows), "n_native_consistent": passed,
        "claims": claims, "ready_for_biological_law_claim": ready,
        "all_frozen_files_unchanged": True, "calibration_adopted": False,
        "limitations": [
            "Native observations and priors have historical exposure; consistency is necessary for transfer but not independent validation",
            "Rounding ranges are representation bounds, not measurement confidence intervals; censored readouts remain unquantified",
            "No fit, model-selection change, solver-policy relaxation, product outcome or new measurement array enters this gate",
            "Strong biological claims remain blocked until dedicated authoritative assessors and genuinely independent evidence exist",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    with (output / "readiness.json").open("x") as handle:
        json.dump(report, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check native physics and verified evidence before claiming biological law learning")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--native-mode", choices=("quantified_point", "printed_rounding"), default="printed_rounding")
    parser.add_argument("--timeout-s", type=int, default=15)
    args = parser.parse_args(argv)
    report = run_readiness(Path(__file__).resolve().parents[1], args.output_dir,
                           mode=args.native_mode, timeout_s=args.timeout_s)
    print(json.dumps({key: report[key] for key in
                      ("assessment", "n_native_conditions", "n_native_consistent",
                       "ready_for_biological_law_claim", "all_frozen_files_unchanged", "calibration_adopted")}, indent=2))
    return 0 if report["ready_for_biological_law_claim"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
