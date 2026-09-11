from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
artifacts = importlib.import_module("ystwin.artifacts")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Stage registered reproductions, compare semantics, and explicitly adopt fresh receipts. Auditing and staging never replace retained outputs.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--run", help="exact registered run id to reproduce")
    action.add_argument("--adopt-receipt", type=Path, help="explicitly persist a verified reproduction receipt in the artifact registry")
    action.add_argument("--review-template", type=Path, help="produce an UNAPPROVED exact-hash scientific-change review template")
    action.add_argument("--inspect", nargs="?", const="all", help="read-only contract inspection, optionally for one run")
    action.add_argument("--check-contracts", action="store_true", help="read-only validation of retained current artifact schemas")
    action.add_argument("--plan", action="store_true", help="read-only dependency-ordered current recipes and input blockers")
    parser.add_argument("--output-dir", type=Path, help="new or empty temporary staging directory outside the repository")
    parser.add_argument("--private-inputs-from-env", type=Path,
                        help="explicit read-only opt-in for policy-approved private inputs from this env data file; never shell-sourced, never copied into Git")
    parser.add_argument("--output", type=Path, help="new temporary JSON file for an unapproved review template")
    parser.add_argument("--artifact", action="append", help="adopt only these exact full-path artifacts; never silently skip changed siblings")
    parser.add_argument("--accept-scientific-changes", action="store_true")
    parser.add_argument("--review", type=Path, help="reviewed exact-hash acceptance record, including any negative-case changes")
    parser.add_argument("--install", action="store_true", help="explicitly replace accepted artifacts after all receipt, review and staleness checks")
    parser.add_argument("--timeout", type=float, default=600.0)
    args = parser.parse_args(argv)
    if args.run and args.output_dir is None:
        parser.error("--run requires a temporary --output-dir")
    if not args.adopt_receipt and (args.install or args.accept_scientific_changes or args.review or args.artifact):
        parser.error("adoption flags require --adopt-receipt")
    if args.output is not None and args.review_template is None:
        parser.error("--output is only for an unapproved --review-template")
    if args.private_inputs_from_env is not None and not (args.run or args.adopt_receipt or args.review_template):
        parser.error("private input access is only for an explicit run, receipt validation/adoption, or review template")
    try:
        if args.inspect:
            result = artifacts.inspect_contracts(ROOT, None if args.inspect == "all" else args.inspect)
            code = 0
        elif args.plan:
            result = artifacts.plan_current_runs(ROOT)
            code = 0
        elif args.check_contracts:
            registry = artifacts.load_registry(ROOT)
            result = []
            for run in registry["runs"].values():
                if run["lifecycle"] != "current":
                    continue
                for item in run["artifacts"]:
                    path = artifacts.artifact_path(ROOT, item["path"])
                    check = artifacts.semantic_comparison(path, path, item.get("comparison"))
                    status = "PASS" if check["matches"] else "SKIP" if item.get("planned") and not path.exists() else "FAIL"
                    result.append({"run": run["id"], "path": item["path"], "status": status, "differences": check["differences"]})
            failures = [item for item in result if item["status"] == "FAIL"]
            result = {"checked": len(result), "passed": sum(item["status"] == "PASS" for item in result),
                      "planned_absent": sum(item["status"] == "SKIP" for item in result), "failures": failures}
            code = int(bool(failures))
        elif args.review_template:
            result = artifacts.review_template(ROOT, args.review_template, private_inputs_from_env=args.private_inputs_from_env)
            if args.output is not None:
                path = args.output.resolve()
                if path.is_relative_to(ROOT.resolve()):
                    raise artifacts.ArtifactError("review templates must be written to a new temporary file outside the repository")
                with path.open("x", encoding="utf-8") as stream:
                    json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
                    stream.write("\n")
            code = 0
        elif args.adopt_receipt:
            result = artifacts.adopt_reproduction(ROOT, args.adopt_receipt, artifact_ids=args.artifact,
                                                 accept_scientific_changes=args.accept_scientific_changes,
                                                 review_path=args.review, install=args.install,
                                                 private_inputs_from_env=args.private_inputs_from_env)
            code = 0
        else:
            receipt = artifacts.reproduce_run(ROOT, args.run, args.output_dir, timeout=args.timeout,
                                               private_inputs_from_env=args.private_inputs_from_env)
            result = {"run_id": receipt["run_id"], "success": receipt["success"],
                      "execution_complete": receipt["execution_complete"], "exit_code": receipt["exit_code"],
                      "receipt": str(args.output_dir / "reproduction.json"),
                      "artifacts": {name: {"matches": check["matches"], "dimensions": check["dimensions"],
                                           "differences": check["differences"], "negative_case_changes": len(check.get("negative_case_changes", []))}
                                    for name, check in receipt["comparisons"].items()},
                      "undeclared_reads": receipt["execution"].get("undeclared_read_paths", [])}
            code = 0 if receipt["success"] else 1
    except (artifacts.ArtifactError, OSError, KeyError, TypeError, SyntaxError, ValueError) as error:
        reason = ("private reproduction filesystem operation failed; source locations withheld"
                  if args.private_inputs_from_env is not None and isinstance(error, OSError) else str(error))
        print(json.dumps({"success": False, "status": "investigation_pending", "reason": reason}, allow_nan=False))
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
