from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ystwin.analysis.portable_evidence import export_portable_evidence


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export path-neutral historical evidence without fitting or altering originals")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--destination", type=Path, default=Path("data/frozen_evidence/native_v1"))
    args = parser.parse_args(argv)
    root = args.root.resolve(strict=True)
    destination = args.destination if args.destination.is_absolute() else root / args.destination
    manifest = export_portable_evidence(root, destination)
    digest = hashlib.sha256((destination / "manifest.json").read_bytes()).hexdigest()
    print(json.dumps({
        "manifest": "manifest.json", "sha256": digest,
        "records": len(manifest["records"]),
        "storage_transforms": sum(len(record["transforms"]) for record in manifest["records"]),
        "lineage_only_references": len(manifest["external_references"]),
        "new_scientific_claim_authorized": False,
    }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
