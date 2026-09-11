from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from ystwin.analysis.frozen_runtime import HistoricalReplayError, replay_historical_evidence
from ystwin.analysis.portable_replay import replay_portable_evidence


def main(argv=None):
    parser = argparse.ArgumentParser(description="Verify public evidence and replay recorded native predictions and score arithmetic without fitting")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--private-original-root", type=Path)
    parser.add_argument("--output", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--historical", action="store_true",
                      help="Replay with the registered frozen Git source and exact recorded runtime in an isolated child")
    mode.add_argument("--frozen-code-ref", help="Explicit full registered historical Git commit; never a branch or current-code fallback")
    parser.add_argument("--source-repository", type=Path,
                        help="Git worktree containing the frozen runtime registry and commit; defaults to --root")
    parser.add_argument("--python-executable", type=Path,
                        help="Trusted child interpreter with the recorded Python and package versions; defaults to this interpreter")
    args = parser.parse_args(argv)
    historical = args.historical or args.frozen_code_ref is not None
    if not historical and (args.source_repository is not None or args.python_executable is not None):
        parser.error("--source-repository and --python-executable require --historical or --frozen-code-ref")
    root = args.root.resolve(strict=True)
    output = args.output
    if output is not None:
        output = output if output.is_absolute() else root / output
        if output.resolve() != output or not output.parent.is_dir():
            raise ValueError("replay output requires an existing nonsymlink parent and a canonical path")
        if output.exists():
            raise FileExistsError("portable replay output is write-once; existing records cannot be overwritten")
    try:
        if historical:
            report = replay_historical_evidence(
                root, source_repository=args.source_repository, frozen_code_ref=args.frozen_code_ref,
                python_executable=args.python_executable, manifest_path=args.manifest,
                expected_sha256=args.manifest_sha256, private_original_root=args.private_original_root,
            )
        else:
            report = replay_portable_evidence(
                root, manifest_path=args.manifest, expected_sha256=args.manifest_sha256,
                private_original_root=args.private_original_root,
            )
    except HistoricalReplayError as exc:
        print(json.dumps(exc.report, sort_keys=True, indent=2, allow_nan=False), file=sys.stderr)
        if exc.returncode is not None and exc.returncode < 0:
            return 128 - exc.returncode
        return exc.returncode or 1
    data = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if output is not None:
        with output.open("x") as handle:
            handle.write(data)
    print(data, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
