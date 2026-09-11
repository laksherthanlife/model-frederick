from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
development = importlib.import_module("ystwin.analysis.native_population_development")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Non-authorizing, source-gated population-mean development under the separate fixed protocol.")
    parser.add_argument("command", nargs="?", choices=("status", "contract", "freeze-implementation", "request", "refresh-request", "fit", "score"), default="status")
    parser.add_argument("--protocol", type=Path, default=ROOT / "data/native_law_v2/development_protocol.json")
    parser.add_argument("--allowlist", type=Path)
    parser.add_argument("--executor-admission", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--exporter-code", type=Path, action="append", default=[])
    parser.add_argument("--requester")
    parser.add_argument("--freeze", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "contract":
            result = {"scientific_scope": development.SCIENTIFIC_SCOPE, "admission_contract": development.admission_contract(), "executor_refresh_contract": development.executor_refresh_contract(), "environment": development.environment_record(), "canonical_parameters": development.PARAMETERS, "strong_claim_authorized": False}
        elif args.command == "status":
            result = development.check_custody_gate(ROOT, args.protocol, args.allowlist, executor_admission_path=args.executor_admission)
        elif args.command in {"freeze-implementation", "request", "refresh-request", "fit"}:
            if args.output is None:
                parser.error("this command requires a new root-relative --output directory under outputs/native_population_development")
            if args.command == "freeze-implementation":
                result = development.write_implementation_configuration(ROOT, args.protocol, args.output)
            elif args.command == "request":
                if not args.exporter_code or not args.requester:
                    parser.error("request requires actual --exporter-code paths and an actual --requester identity")
                result = development.write_release_request(ROOT, args.protocol, args.output, args.exporter_code, args.requester)
            elif args.command == "refresh-request":
                if not args.requester:
                    parser.error("refresh-request requires an actual --requester identity")
                result = development.write_executor_refresh_request(ROOT, args.protocol, args.allowlist, args.output, args.requester)
            else:
                result = development.execute_training(ROOT, args.protocol, args.allowlist, args.output, executor_admission_path=args.executor_admission)
        else:
            if args.freeze is None:
                parser.error("score requires the exact previously sealed --freeze artifact")
            path = args.freeze.relative_to(ROOT) if args.freeze.is_absolute() else args.freeze
            development.require(path.parts[:2] == ("outputs", "native_population_development"), "only this development run's frozen predictions may be supplied to scoring")
            reference = development.artifact_reference(ROOT, path, "development_predictions")
            result = development.execute_scoring(ROOT, args.protocol, args.allowlist, reference)
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        return 2 if result.get("ready") is False else 0
    except (development.DevelopmentError, OSError, KeyError, TypeError, ValueError) as error:
        print(json.dumps({"scientific_scope": development.SCIENTIFIC_SCOPE, "status": "blocked_no_biological_authorization", "reason_code": getattr(error, "code", "dependency_or_partition_violation"), "reason": str(error), "affected_group_ids": getattr(error, "groups", []), "strong_claim_authorized": False, "further_exports_permitted": False}, indent=2, sort_keys=True, allow_nan=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
