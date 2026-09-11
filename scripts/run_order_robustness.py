from __future__ import annotations

import argparse
import csv
import io
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ystwin.analysis.partial_orders import (
    Bound,
    Calibration,
    MagnitudeDomain,
    OrderClaim,
    Scope,
    Source,
    affine_extrema,
    compile_order,
    evaluate_forward,
    sample_uniform_box,
)


REPORT_NAME = "partial_order_synthetic_demo.json"
SUMMARY_NAME = "partial_order_synthetic_summary.csv"


def build_demo(*, samples: int, seed: int, budget: int, tolerance: float) -> dict:
    scope = Scope(
        "magnitude", "synthetic fixed engineering context", "toy magnitude", "toy_unit",
        Calibration("synthetic shared readout", "shared_unknown_positive"),
    )
    order_source = Source(
        "synthetic", "synthetic-order-v1", "declared toy a < b; no measurement or model sensitivity",
    )
    bounds_source = Source(
        "synthetic", "synthetic-box-v1",
        "chosen [0,1] bounds in target toy coordinates, not derived from assays or from order",
    )
    gap_source = Source(
        "synthetic", "synthetic-weak-embedding-v1", "declared zero gap; weak a <= b relaxation",
    )
    proposal_source = Source(
        "synthetic", "synthetic-uniform-proposal-v1",
        "independent uniform quotient-coordinate box conditioned on a <= b; not a biological posterior",
    )
    poset = compile_order(("a", "b", "c"), scope, [OrderClaim("a", "b", scope, order_source)])
    qualitative_only = MagnitudeDomain(poset, gap=0.0, gap_source=gap_source, tolerance=tolerance)
    domain = MagnitudeDomain(
        poset, [Bound(item, 0.0, 1.0, scope, bounds_source) for item in poset.items],
        gap=0.0, gap_source=gap_source, tolerance=tolerance,
    )
    affine = {
        "robust": affine_extrema(
            domain, {"b": 1, "a": -1}, offset=0.125,
            estimand="synthetic b - a + 0.125", unit="toy_response",
        ),
        "nonrobust": affine_extrema(
            domain, {"c": 1}, offset=-0.5,
            estimand="synthetic c - 0.5", unit="toy_response",
        ),
        "weak_order_contrast": affine_extrema(
            domain, {"b": 1, "a": -1}, estimand="synthetic b - a", unit="toy_response",
        ),
    }
    batch = sample_uniform_box(domain, n=samples, seed=seed, budget=budget, source=proposal_source)

    def positive(point):
        return (point["b"] - point["a"]) ** 2 + 0.01

    def oscillating(point):
        return math.sin(2 * math.pi * point["c"]) - (point["b"] - point["a"]) / 4

    def with_failure(point):
        if point["c"] > 0.8:
            raise RuntimeError("deliberate synthetic callback failure for c > 0.8")
        return positive(point)

    nonlinear = {
        "robust_on_samples": evaluate_forward(
            batch, positive, estimand="synthetic (b-a)^2 + 0.01", unit="toy_response",
        ),
        "nonrobust_on_samples": evaluate_forward(
            batch, oscillating, estimand="synthetic sin(2*pi*c) - (b-a)/4", unit="toy_response",
        ),
        "with_failures": evaluate_forward(
            batch, with_failure,
            estimand="synthetic (b-a)^2 + 0.01, with deliberate failure when c > 0.8",
            unit="toy_response",
        ),
    }
    return {
        "schema_version": 2,
        "synthetic_only": True, "biology_data_used": False, "law_ready": False,
        "description": "SYNTHETIC engineering demonstration only; no biological facts or frozen model edits",
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "solver": "scipy.optimize.linprog, HiGHS"},
        "limitations": [
            "a < b is a synthetic qualitative claim; c is incomparable, not tied or arbitrarily ranked",
            "the declared zero-gap numeric embedding permits a = b; it is not strict separation",
            "shared unknown positive gain preserves the order but does not identify an absolute range",
            "[0,1] bounds are extra synthetic engineering assumptions, never inferred from the poset",
            "the poset alone supplies no probability distribution; requesting a uniform box does not define one",
            "base uniform proposals require finite nonempty explicit boxes; conditioning requires positive base probability",
            "distribution metadata describes a defined mathematical law, not whether a finite sampling budget succeeded",
            "measure checks use exact ratios of the declared binary-float inputs, not exact biological magnitudes",
            "LP results are residual-validated numerical certificates, not rational exact proofs",
            "LP feasibility uses absolute and relative tolerances; near-boundary signs remain numerical",
            "nonlinear extrema and signs apply only to successful sampled evaluations, never global worst cases",
            "all callback failures and their sample coordinates remain in the denominator and records",
            "no biological data, missing-assay bounds, model parameter hierarchies, or law-ready claims are used",
        ],
        "order_only": affine_extrema(
            qualitative_only, {"a": 1}, estimand="synthetic a without range assumptions", unit="toy_unit",
        ),
        "order_only_sampling": sample_uniform_box(
            qualitative_only, n=samples, seed=seed, budget=budget, source=proposal_source,
        ).report(),
        "affine": affine, "sampling": batch.report(), "nonlinear": nonlinear,
    }


def summary_csv(report: dict) -> str:
    rows = []
    for name, result in {"order_only": report["order_only"], **report["affine"]}.items():
        rows.append({
            "query": name, "method": "affine LP", "status": result["status"],
            "minimum": result["minimum"]["value"], "maximum": result["maximum"]["value"],
            "sign": result["sign"], "tolerance": result["tolerance"],
            "range_assumptions": "none" if name == "order_only" else "a,b,c in [0,1], synthetic-box-v1",
            "guarantee": "numerical affine LP certificate on the declared domain, not rational exact proof",
        })
    for name, result in report["nonlinear"].items():
        rows.append({
            "query": name, "method": "sampled nonlinear callback", "status": result["sign_stability"],
            "minimum": result["sampled_minimum"], "maximum": result["sampled_maximum"],
            "sign": result["successful_sample_sign"], "tolerance": result["sign_tolerance"],
            "accepted": result["accepted"], "successful": result["successful"], "failed": result["failed"],
            "requested": report["sampling"]["requested"], "attempted": report["sampling"]["attempted"],
            "seed": report["sampling"]["seed"], "budget": report["sampling"]["budget"],
            "range_assumptions": "a,b,c in [0,1], synthetic-box-v1",
            "requested_proposal": result["sampling"]["requested_proposal"],
            "proposal_distribution": result["sampling"]["proposal_distribution"],
            "proposal": result["sampling"]["distribution"],
            "conditioning_status": result["sampling"]["conditioning_status"],
            "proposal_source": result["sampling"]["proposal_source"]["reference"],
            "guarantee": "sampled successful evaluations only; no global nonlinear certificate",
        })
    text = io.StringIO(newline="")
    columns = [
        "synthetic_only", "query", "method", "status", "minimum", "maximum", "sign", "tolerance",
        "accepted", "successful", "failed", "requested", "attempted", "seed", "budget",
        "range_assumptions", "order_assumptions", "requested_proposal", "proposal_distribution",
        "proposal", "conditioning_status", "proposal_source", "guarantee", "details_file",
    ]
    writer = csv.DictWriter(text, fieldnames=columns)
    writer.writeheader()
    for row in rows:
        writer.writerow({
            "synthetic_only": True,
            "order_assumptions": "synthetic-order-v1: a<b, c incomparable; synthetic-weak-embedding-v1: gap=0",
            "details_file": REPORT_NAME,
            **row,
        })
    return text.getvalue()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="SYNTHETIC partial-order robustness demo; no biology data")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs")
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--budget", type=int, default=10000)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args(argv)
    if not args.output_dir.is_dir():
        parser.error("output directory must already exist")
    paths = (args.output_dir / REPORT_NAME, args.output_dir / SUMMARY_NAME)
    if any(path.exists() for path in paths):
        parser.error("refusing to overwrite existing partial_order results; choose a fresh output directory")
    try:
        report = build_demo(samples=args.samples, seed=args.seed, budget=args.budget, tolerance=args.tolerance)
        payloads = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", summary_csv(report))
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))
    try:
        for path, payload in zip(paths, payloads):
            with path.open("x", encoding="utf-8", newline="") as stream:
                stream.write(payload)
    except FileExistsError as exc:
        parser.error(f"refusing to overwrite existing output: {exc.filename}")
    print("SYNTHETIC engineering demonstration only: no biology data, no law-ready claim.")
    sampling = report["sampling"]
    print(f"Uniform rejection: {sampling['status']}, {sampling['accepted']}/{sampling['requested']} accepted, "
          f"{sampling['attempted']} proposals; seed={sampling['seed']}, budget={sampling['budget']}.")
    for path in paths:
        print(path)
    return 0 if sampling["status"] == "complete" and all(
        result["status"] == "bounded" for result in report["affine"].values()
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
