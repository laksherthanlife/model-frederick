from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
scoring = importlib.import_module("ystwin.analysis.native_population_scoring")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Phase-C-only scoring adapter with separate training, executor and response approvals; no fitting.")
    parser.add_argument("command", nargs="?", choices=("status", "score", "verify-published"), default="status")
    parser.add_argument("--manifest", type=Path, default=ROOT / scoring.MANIFEST_PATH)
    parser.add_argument("--manifest-sha256", default=scoring.MANIFEST_SHA256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "score":
            if args.output is None:
                parser.error("score requires a fresh --output directory under outputs/native_population_development")
            result = scoring.score_once(ROOT, args.manifest, args.output, args.manifest_sha256)
        elif args.command == "verify-published":
            if args.output is None:
                parser.error("verify-published requires a fresh --output directory for fix/version lineage, not a replacement report")
            result = scoring.verify_published_results(ROOT, args.output)
        else:
            bundle = scoring.inspect_phase_c(ROOT, args.manifest, args.manifest_sha256)
            result = {"adapter_version": scoring.ADAPTER_VERSION, "scientific_scope": scoring.SCIENTIFIC_SCOPE, "status": "phase_c_metadata_lineage_verified", "training_data_approval": bundle.training.data_approval, "executor_approval": bundle.training.executor_approval, "prediction_freeze": bundle.prediction.reference, "response_release_approval": bundle.release.manifest["custody_decision"], "response_projection": bundle.release.manifest["response_projection"], "response_values_opened": False, "strong_claim_authorized": False}
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except (scoring.ScoringError, scoring.frozen.DevelopmentError, OSError, KeyError, TypeError, ValueError) as error:
        print(json.dumps({"adapter_version": scoring.ADAPTER_VERSION, "scientific_scope": scoring.SCIENTIFIC_SCOPE, "status": "refused_phase_c_scoring", "reason": str(error), "strong_claim_authorized": False, "final_test_defined": False}, indent=2, sort_keys=True, allow_nan=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
