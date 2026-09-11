import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy

from ystwin.analysis import native_rate_learning
from ystwin.analysis.native_rate_learning import ContractError, content_sha256, fit_rate_laws


ROOT = Path(__file__).resolve().parents[1]


def read_bytes(path):
    path = Path(path)
    if not 0 < path.stat().st_size <= 64 * 1024 * 1024:
        raise ContractError("artifact size exceeds the explicit 64 MiB read budget")
    return path.read_bytes()


def parse_object(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ContractError("duplicate JSON keys are not permitted")
            result[key] = value
        return result

    value = json.loads(raw, object_pairs_hook=pairs)
    if not isinstance(value, dict):
        raise ContractError("artifact must be an explicit JSON object; no tabular columns are inferred")
    content_sha256(value)
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description="Fit explicit native-rate projections only after a frozen learning contract; never assess biological claims.")
    parser.add_argument("--protocol", type=Path, default=ROOT / "data/native_law_v2/protocol.json")
    parser.add_argument("--catalogue", type=Path, default=ROOT / "data/native_law_v2/biochemical_catalogue.json")
    parser.add_argument("--measurements", type=Path)
    parser.add_argument("--engineering-only", action="store_true", help="Accept explicitly synthetic engineering fixtures, never biological training evidence.")
    args = parser.parse_args(argv)
    access = "no_projection_opened"
    try:
        protocol_path = args.protocol.resolve()
        protocol_raw = read_bytes(protocol_path)
        protocol = parse_object(protocol_raw)
        catalogue = None
        catalogue_raw = None
        if args.engineering_only:
            if protocol.get("evidence_kind") != "synthetic_engineering" or protocol.get("approval", {}).get("status") != "synthetic_engineering_only":
                raise ContractError("engineering-only execution requires an explicitly synthetic engineering protocol")
            if args.measurements is None:
                raise ContractError("engineering-only execution requires an explicit synthetic measurement projection")
            measurement_path = args.measurements.resolve()
            expected_kind = "synthetic_engineering"
        else:
            contract = protocol.get("learning_contract")
            if protocol.get("schema_version") != 1 or not isinstance(contract, dict):
                raise ContractError("independent protocol has no approved learning_contract; no measurement projection was opened")
            approval = contract.get("approval")
            if not isinstance(approval, dict) or approval.get("status") != "approved_native_training":
                raise ContractError("independent learning_contract has no native training approval")
            projection = approval.get("projection_path")
            if not isinstance(projection, str) or not projection:
                raise ContractError("approval must explicitly name projection_path relative to the protocol directory")
            projection = Path(projection)
            if projection.is_absolute():
                raise ContractError("approved projection_path must be relative to the protocol directory")
            measurement_path = (protocol_path.parent / projection).resolve()
            if not measurement_path.is_relative_to(protocol_path.parent):
                raise ContractError("approved projection must remain inside the protocol directory")
            if args.measurements is not None and args.measurements.resolve() != measurement_path:
                raise ContractError("requested projection is not the explicit approved projection")
            expected_file_digest = approval.get("projection_file_sha256")
            if not isinstance(expected_file_digest, str) or len(expected_file_digest) != 64:
                raise ContractError("approval must bind the exact projection file digest before access")
            catalogue_raw = read_bytes(args.catalogue.resolve())
            catalogue = parse_object(catalogue_raw)
            if content_sha256(catalogue) != approval.get("catalogue_sha256"):
                raise ContractError("catalogue digest does not match the independent training approval")
            allowlist = catalogue.get("access_policy", {}).get("learner_measurement_allowlist")
            if not isinstance(allowlist, list) or not allowlist or approval.get("catalogue_allowlist_entry") not in allowlist:
                raise ContractError("custodian catalogue does not allow the approved learner projection")
            expected_kind = "native_measurement"
        measurement_raw = read_bytes(measurement_path)
        access = "synthetic_projection_opened" if args.engineering_only else "approved_projection_bytes_opened"
        measurement_digest = hashlib.sha256(measurement_raw).hexdigest()
        if not args.engineering_only and measurement_digest != expected_file_digest:
            raise ContractError("projection file digest mismatch; measurement values were not parsed")
        measurements = parse_object(measurement_raw)
        if measurements.get("evidence_kind") != expected_kind:
            raise ContractError("projection evidence kind is not authorized for this execution mode")
        result = fit_rate_laws(measurements, protocol, catalogue=catalogue)
        result["execution"] = {
            "scope": "learning_only_not_independent_verification",
            "measurement_access": access,
            "measurement_projection": str(measurement_path),
            "measurement_file_sha256": measurement_digest,
            "protocol_file_sha256": hashlib.sha256(protocol_raw).hexdigest(),
            "catalogue_file_sha256": hashlib.sha256(catalogue_raw).hexdigest() if catalogue_raw is not None else None,
            "learner_code_sha256": hashlib.sha256(read_bytes(native_rate_learning.__file__)).hexdigest(),
            "runner_code_sha256": hashlib.sha256(read_bytes(__file__)).hexdigest(),
            "python_version": sys.version.split()[0],
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "authority_receipts": "not_authenticated_by_this_learning_process",
        }
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))
        return 0 if result["status"] == "fitted" else 1
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(json.dumps({"schema_version": "native_rate_runner.v1", "status": "blocked", "reason": str(error), "error_type": type(error).__name__, "measurement_access": access}, sort_keys=True, allow_nan=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
