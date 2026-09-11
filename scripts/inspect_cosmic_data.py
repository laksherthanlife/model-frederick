from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ystwin.analysis.cosmic_data import (
    CosmicDataError,
    PUBLISHER_WORKBOOK_BYTES,
    PUBLISHER_WORKBOOK_SHA256,
    WORKBOOK_FILENAMES,
    load_cosmic_workbook,
    summarize_cosmic_data,
)


_REPO_ROOT = Path(__file__).resolve().parents[1]


def _input_workbooks(args: argparse.Namespace) -> list[Path]:
    if args.workbook is not None:
        path = args.workbook.expanduser()
        if not path.is_file():
            raise CosmicDataError(f"Explicit workbook is not a file: {path}")
        return [path]
    directory = args.input_dir.expanduser()
    if not directory.is_dir():
        raise CosmicDataError(f"Explicit input directory does not exist: {directory}")
    candidates = [path for path in directory.iterdir() if path.suffix.lower() == ".xlsx"]
    unknown = [path.name for path in candidates if path.name not in WORKBOOK_FILENAMES]
    if unknown:
        raise CosmicDataError(
            f"Ambiguous workbook directory contains unregistered XLSX names {sorted(unknown)!r}; "
            "select a single file with --workbook (its bytes must still match the pinned source)"
        )
    paths = [directory / name for name in WORKBOOK_FILENAMES if (directory / name).is_file()]
    if not paths:
        raise CosmicDataError("No recognized workbook in the explicit directory; use --workbook for a renamed source")
    if len(paths) != len(candidates):
        raise CosmicDataError("Ambiguous workbook directory includes a non-file XLSX entry")
    return paths


def _inside_directory(path: Path, directory: Path) -> bool:
    return any(parent.exists() and parent.samefile(directory) for parent in path.parents)


def _report_destination(destination: Path, sources: Sequence[Path]) -> Path:
    destination = destination.expanduser()
    if not destination.is_absolute():
        raise CosmicDataError("--report requires an explicit absolute external JSON destination")
    path = destination.resolve()
    if path.suffix.lower() != ".json":
        raise CosmicDataError("--report must name a JSON diagnostic, never a source workbook")
    if _inside_directory(path, _REPO_ROOT):
        raise CosmicDataError("Report destination must be outside the repository; raw-data publication is not approved")
    for source in sources:
        if _inside_directory(path, source.resolve().parent) or _inside_directory(path, source.absolute().parent):
            raise CosmicDataError("Report destination must be outside the source directory")
    if destination.is_symlink() or path.exists():
        raise CosmicDataError(f"Report destination already exists; nothing will be overwritten: {path}")
    if not path.parent.is_dir():
        raise CosmicDataError("Report parent directory must already exist")
    return path


def _source_snapshot(path: Path) -> tuple[str, int]:
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    if digest != PUBLISHER_WORKBOOK_SHA256 or len(payload) != PUBLISHER_WORKBOOK_BYTES:
        raise CosmicDataError(f"Publisher workbook SHA-256/byte-count mismatch: {path.name}")
    return digest, len(payload)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect the pinned COSMIC CHO workbook read-only, preserving source roles, units and cells. "
            "This is not model reproduction or independent validation. No default source directory is used."
        ),
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input-dir", type=Path, help="Explicit directory containing the publisher workbook and/or its identical alias")
    inputs.add_argument("--workbook", type=Path, help="Explicit workbook file; acceptance requires the pinned publisher SHA-256")
    parser.add_argument("--report", type=Path, help="Optional new absolute JSON path outside the repository and source directory; never overwritten")
    args = parser.parse_args(argv)
    try:
        sources = _input_workbooks(args)
        destination = _report_destination(args.report, sources) if args.report is not None else None
        before = [_source_snapshot(path) for path in sources]
        data = load_cosmic_workbook(sources[0])
        report = summarize_cosmic_data(data)
        after = [_source_snapshot(path) for path in sources]
        if before != after:
            raise CosmicDataError("Source bytes changed during read-only inspection")
        report["input_files"] = [
            {"filename": path.name, "bytes": initial[1], "sha256_before": initial[0],
             "sha256_after": final[0], "unchanged": initial == final}
            for path, initial, final in zip(sources, before, after, strict=True)
        ]
        text = json.dumps(report, indent=2, ensure_ascii=True, allow_nan=False) + "\n"
        if destination is not None:
            destination = _report_destination(destination, sources)
            with destination.open("x", encoding="utf-8") as stream:
                stream.write(text)
        print(text, end="")
    except (CosmicDataError, OSError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
