from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ystwin.analysis.mechanistic_reference import (
    FrozenReference,
    default_multisource_reference_design,
    default_reference_design,
    evaluate_native_reference,
    freeze_reference,
    generate_reference,
    reference_design_from_dict,
    run_synthetic_benchmark,
    score_parameter_recovery,
    score_synthetic_reference,
    training_mean_baseline,
)
from ystwin.analysis.parameter_evidence import (
    EvidenceGap,
    evidence_report,
    load_parameter_evidence,
    model_parameter_inventory,
)
from ystwin.mech.kinetic_sbml import KineticSimulationError


def _json(value) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=_unique_pairs)


def _write_new(path: Path, value) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(_json(value))


def _export(directory: Path, reference, score: dict) -> None:
    if not directory.parent.is_dir():
        raise ValueError("the explicit output directory's parent must already exist")
    directory.mkdir()
    learner = directory / "learner"
    evaluator = directory / "evaluator"
    learner.mkdir()
    evaluator.mkdir()
    _write_new(learner / "inputs.json", reference.learner_inputs())
    _write_new(evaluator / "freeze.json", reference.frozen.to_dict())
    _write_new(evaluator / "truth.json", reference.evaluator_truth())
    _write_new(evaluator / "benchmark.json", score)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Frozen HOG and explicitly declared source-native nutrient challenges; no physical nutrient/time conversion or full biological validation."
    )
    commands = parser.add_subparsers(dest="command", required=True)
    inventory_parser = commands.add_parser("inventory", help="audit source hashes and programme blockers without fitting")
    inventory_parser.add_argument("--model-parameters", help="enumerate a locally verified SBML parameter vector and equations")
    inventory_parser.add_argument("--require-complete", action="store_true", help="exit 2 while full biological coverage is blocked")
    synthetic = commands.add_parser("synthetic", help="freeze then generate independent reference challenges")
    design_options = synthetic.add_mutually_exclusive_group()
    design_options.add_argument("--design", type=Path, help="explicit HOG ReferenceDesign or declared MultiSourceReferenceDesign JSON; no learner coefficients or outcomes")
    design_options.add_argument("--multi-source", action="store_true", help="HOG plus verified Jalihal PKA/Snf1/TOR; all five normalized nutrient axes, native clocks and source counterfactuals")
    synthetic.add_argument("--seed", type=int, default=0, help="seed for the default design only")
    synthetic.add_argument("--declare-no-product-outcomes", action="store_true", required=True)
    synthetic.add_argument("--output-dir", type=Path, help="new directory; learner inputs and evaluator truth are separate; never overwrite")
    score_parser = commands.add_parser("score-synthetic", help="regenerate frozen truth and score external learner predictions")
    score_parser.add_argument("--freeze", required=True, type=Path)
    score_parser.add_argument("--predictions", required=True, type=Path)
    parameter_parser = commands.add_parser("score-parameters", help="source-coordinate parameter differences, not an identifiability certificate")
    parameter_parser.add_argument("--freeze", required=True, type=Path)
    parameter_parser.add_argument("--estimates", required=True, type=Path)
    native_parser = commands.add_parser("native-replay", help="real native HPLC fit replay only, explicitly not independent validation")
    native_parser.add_argument("--freeze", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        inventory = load_parameter_evidence()
        if args.command == "inventory":
            report = evidence_report(inventory)
            if args.model_parameters:
                report["model_parameter_inventory"] = model_parameter_inventory(args.model_parameters, inventory)
            print(_json(report), end="")
            blocked = (args.require_complete and not report["full_biological_coverage"]) or not all(
                check["verified"] for check in report["source_checks"]
            )
            return 2 if blocked else 0
        if args.command == "synthetic":
            if args.output_dir is not None and args.output_dir.exists():
                raise FileExistsError(f"refusing to overwrite or restamp {args.output_dir}")
            design = (reference_design_from_dict(_read(args.design)) if args.design else
                      (default_multisource_reference_design(seed=args.seed) if args.multi_source else default_reference_design(seed=args.seed)))
            frozen = freeze_reference(design, product_outcomes_seen=False, inventory=inventory)
            reference = generate_reference(frozen, inventory=inventory)
            result = run_synthetic_benchmark(reference, training_mean_baseline)
            result["learner"] = "training_mean_baseline_not_a_mechanistic_learner"
            if args.output_dir is not None:
                _export(args.output_dir, reference, result)
            summary = {key: value for key, value in result.items() if key != "rows"}
            summary["n_experiments"] = len(design.experiments)
            summary["source_ids"] = frozen.to_dict()["evidence"]["source_ids"]
            summary["source_contracts"] = frozen.to_dict().get("source_contracts", {})
            summary["n_training_observations"] = len(reference.learner_inputs()["training_observations"])
            summary["output_dir"] = None if args.output_dir is None else str(args.output_dir)
            summary["custody_limit"] = "projection and directory separation only; evaluator truth is not an OS-protected secret"
            print(_json(summary), end="")
            return 0
        frozen = FrozenReference.from_dict(_read(args.freeze))
        if args.command in ("score-synthetic", "score-parameters"):
            reference = generate_reference(frozen, inventory=inventory)
            result = (score_synthetic_reference(reference, _read(args.predictions))
                      if args.command == "score-synthetic"
                      else score_parameter_recovery(reference, _read(args.estimates)))
        else:
            result = evaluate_native_reference(frozen, evaluation_kind="published_fit_replay", inventory=inventory)
        print(_json(result), end="")
        return 0
    except (EvidenceGap, KineticSimulationError, ValueError, TypeError, KeyError, OSError) as exc:
        print(_json({"status": "blocked", "error": str(exc)}), end="", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
