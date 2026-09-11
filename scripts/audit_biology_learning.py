from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from ystwin.analysis.frozen_runtime import HistoricalReplayError, audit_historical_evidence


def main(argv=None):
    parser = argparse.ArgumentParser(description="Audit native learning claims without changing or retraining the frozen checkpoint")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--require-claim", action="append", default=[])
    parser.add_argument("--reuse-native-diagnostics", action="store_true")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--historical", action="store_true",
                      help="Audit public evidence with verified frozen Git source and runtime; never rerun fitting")
    mode.add_argument("--frozen-code-ref", help="Explicit full registered historical Git commit")
    parser.add_argument("--source-repository", type=Path)
    parser.add_argument("--python-executable", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--manifest-sha256")
    args = parser.parse_args(argv)
    historical = args.historical or args.frozen_code_ref is not None
    if not historical and any(value is not None for value in (
        args.source_repository, args.python_executable, args.manifest, args.manifest_sha256,
    )):
        parser.error("historical source/runtime and manifest options require --historical or --frozen-code-ref")
    root = args.root.resolve(strict=True)
    if historical:
        try:
            report = audit_historical_evidence(
                root, source_repository=args.source_repository, frozen_code_ref=args.frozen_code_ref,
                python_executable=args.python_executable, manifest_path=args.manifest,
                expected_sha256=args.manifest_sha256, required_claims=args.require_claim,
                reuse_native_diagnostics=args.reuse_native_diagnostics, output_dir=args.output_dir,
            )
        except HistoricalReplayError as exc:
            print(json.dumps(exc.report, sort_keys=True, indent=2, allow_nan=False), file=sys.stderr)
            if exc.returncode is not None and exc.returncode < 0:
                return 128 - exc.returncode
            return exc.returncode or 1
        print(json.dumps(report, sort_keys=True, indent=2, allow_nan=False))
        return 0 if report["required_claims_satisfied"] else 2
    from ystwin.analysis.biology_learning_audit import require_claim, run_biology_learning_audit

    ledger = run_biology_learning_audit(root, args.output_dir, reuse_native_diagnostics=args.reuse_native_diagnostics)
    rejected, authorized = [], []
    for claim in args.require_claim:
        try:
            require_claim(ledger, claim, root=root)
            authorized.append(claim)
        except ValueError as exc:
            rejected.append(str(exc))
    print(json.dumps({
        "output_dir": str(args.output_dir), "checkpoint_sha256": ledger["checkpoint_sha256"],
        "all_frozen_files_unchanged": ledger["all_frozen_files_unchanged"],
        "claims": {item["claim"]: item["status"] for item in ledger["claims"]},
        "rejected_required_claims": rejected, "root_verified_required_claims": authorized,
        "claim_summary_kind": "collector_snapshot_not_bearer_authorization",
        "diagnostics_mode": ledger["diagnostics_origin"]["mode"],
        "scope": ledger["scope"],
    }, indent=2, allow_nan=False))
    return 2 if rejected or not ledger["all_frozen_files_unchanged"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
