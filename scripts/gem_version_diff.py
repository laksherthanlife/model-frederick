"""Compare the resolved reference GEM with an explicit candidate, without adopting it.

    python3 scripts/gem_version_diff.py data/gem/yeast-GEM-9.1.1.xml.gz

Writes ``gem_version_diff.csv`` under ``YSTWIN_OUTPUTS`` (default: ``outputs/``), leaving
both models and all model resolvers unchanged. Set ``YSTWIN_OUTPUTS`` to an external
scratch directory to avoid overwriting the retained comparison. Other XML or XML.gz
candidates can still be compared by supplying their paths explicitly.

``--verify-candidate`` only verifies the retained v9.1.1 comparison asset, offline and
read-only. ``--fetch-candidate`` explicitly permits creation of that asset if missing,
using its immutable GitHub commit/path/blob and both byte identities in
``data/public_inputs.json``. It never overwrites an existing file, refreshes a hash, or
writes a comparison CSV. Ordinary comparisons never fetch; the registered candidate is
verified again before loading. The reference remains v9.0.2 unless explicitly overridden.

The table checks counts, the selected identifiers in ``PINNED``, and maximum objectives
under three glucose bounds with unlimited oxygen. For the retained models the sole
objective is biomass with coefficient one. Identifier survival is not proof of unchanged
stoichiometry, pathway behaviour, or predictive validity. An upgrade still requires
regenerating and reviewing the affected analyses; this comparison does not adopt one.

The retained v9.0.2 -> v9.1.1 table records 4131 -> 4105 reactions, 2806 -> 2748
metabolites, and 1161 -> 1143 genes. All listed identifiers survive, and growth deltas
round to -5.53% at glucose -1.0, -1.5 and -10.0. It contains no wet-lab validation.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import io
import json
import math
import pathlib
import subprocess
import sys
import warnings

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.fba.solver import load_model

warnings.filterwarnings("ignore")

#: Every id this repository looks up by name. A release that drops one of these breaks a
#: specific layer, so the check is by name and not by count.
PINNED = {
    "r_2111": "biomass objective",
    "r_4046": "NGAM, the maintenance layer",
    "r_1714": "glucose exchange",
    "r_1992": "oxygen exchange",
    "r_1761": "ethanol exchange",
    "s_0373": "acetyl-CoA, PHB precursor",
    "s_0189": "GGPP, beta-carotene precursor",
    "s_1543": "UDP-glucose, glycogen precursor",
    "s_1427": "sedoheptulose-7-phosphate, gadusol precursor",
    "s_1212": "NADPH",
    "s_1207": "NADP+",
    "s_0529": "coenzyme A",
    "s_0794": "H+",
}

GLUCOSE_BOUNDS = (-1.0, -1.5, -10.0)
CANDIDATE = pathlib.Path("data/gem/yeast-GEM-9.1.1.xml.gz")
COMPRESSION = {"format": "gzip", "compresslevel": 9, "mtime": 0, "filename": ""}


def _checked_path(path):
    path = pathlib.Path(path).absolute()
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError(f"candidate provenance paths must not contain symlinks: {path}")
    return path


def candidate_entry(root=None):
    root = paths.REPO_ROOT if root is None else pathlib.Path(root)
    document = json.loads(_checked_path(root / "data/public_inputs.json").read_text(encoding="utf-8"))
    entries = [entry for entry in document["files"] if entry["path"] == CANDIDATE.as_posix()]
    if len(entries) != 1:
        raise ValueError("exactly one retained v9.1.1 candidate entry is required")
    entry = entries[0]
    upstream = entry["upstream"]
    if (entry["source"] != "checkout" or entry["compression"] != COMPRESSION
            or upstream["repository"] != "SysBioChalmers/yeast-GEM"
            or upstream["tag"] != "v9.1.1" or upstream["path"] != "model/yeast-GEM.xml"):
        raise ValueError("candidate entry must describe the explicit, checkout-only v9.1.1 comparison asset")
    digests = [(entry["sha256"], 64), (entry["uncompressed_sha256"], 64)]
    digests.extend((upstream[key], 40) for key in ("commit", "tree", "model_tree", "git_blob"))
    for value, length in digests:
        if not isinstance(value, str) or len(value) != length or any(c not in "0123456789abcdef" for c in value):
            raise ValueError("candidate provenance requires pinned hexadecimal digests, not mutable refs")
    for key in ("size_bytes", "uncompressed_size_bytes"):
        if type(entry[key]) is not int or entry[key] <= 0:
            raise ValueError("candidate provenance requires positive pinned byte counts")
    url = f"https://raw.githubusercontent.com/{upstream['repository']}/{upstream['commit']}/{upstream['path']}"
    if upstream["url"] != url:
        raise ValueError("candidate URL does not match its immutable commit/path identity")
    return entry


def _verify_bytes(payload, expected_sha256, expected_size, label):
    if len(payload) != expected_size:
        raise ValueError(f"{label} byte count mismatch: expected {expected_size}, got {len(payload)}")
    actual = hashlib.sha256(payload).hexdigest()
    if actual != expected_sha256:
        raise ValueError(f"{label} SHA-256 mismatch: expected {expected_sha256}, got {actual}")


def _verify_xml(payload, entry):
    _verify_bytes(payload, entry["uncompressed_sha256"], entry["uncompressed_size_bytes"], "candidate SBML")
    blob = hashlib.sha1(f"blob {len(payload)}\0".encode() + payload, usedforsecurity=False).hexdigest()
    if blob != entry["upstream"]["git_blob"]:
        raise ValueError("candidate SBML Git blob identity mismatch")


def _verify_compressed(payload, entry):
    _verify_bytes(payload, entry["sha256"], entry["size_bytes"], "compressed candidate")
    _verify_xml(gzip.decompress(payload), entry)


def verify_candidate(path, entry):
    """Read-only verification, including the uncompressed upstream Git blob identity."""
    _verify_compressed(_checked_path(path).read_bytes(), entry)


def _compress_candidate(payload):
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, compresslevel=9, mtime=0) as handle:
        handle.write(payload)
    return buffer.getvalue()


def _github_json(endpoint):
    result = subprocess.run(["gh", "api", "--method", "GET", endpoint],
                            capture_output=True, check=True, timeout=120)
    return json.loads(result.stdout)


def _download_candidate(entry):
    upstream = entry["upstream"]
    repository = upstream["repository"]
    metadata = _github_json(f"repos/{repository}/contents/{upstream['path']}?ref={upstream['commit']}")
    if (metadata.get("type") != "file" or metadata.get("path") != upstream["path"]
            or metadata.get("sha") != upstream["git_blob"]
            or metadata.get("size") != entry["uncompressed_size_bytes"]):
        raise ValueError("GitHub commit/path/blob identity differs from the pinned candidate")
    blob = _github_json(f"repos/{repository}/git/blobs/{upstream['git_blob']}")
    if (blob.get("sha") != upstream["git_blob"] or blob.get("size") != entry["uncompressed_size_bytes"]
            or blob.get("encoding") != "base64"):
        raise ValueError("GitHub blob metadata differs from the pinned candidate")
    payload = base64.b64decode("".join(blob["content"].split()), validate=True)
    _verify_xml(payload, entry)
    return payload


def provision_candidate(*, root=None, fetch=False):
    """Default to offline verification; explicit acquisition can only create this asset."""
    root = paths.REPO_ROOT if root is None else pathlib.Path(root)
    entry = candidate_entry(root)
    path = _checked_path(root / CANDIDATE)
    if not path.exists() and fetch:
        compressed = _compress_candidate(_download_candidate(entry))
        _verify_compressed(compressed, entry)
        with _checked_path(path).open("xb") as handle:
            handle.write(compressed)
    verify_candidate(path, entry)
    return path, entry


def _growth(model, glucose_lower_bound):
    with model:
        model.reactions.get_by_id("r_1714").lower_bound = glucose_lower_bound
        model.reactions.get_by_id("r_1992").lower_bound = -1000.0
        solution = model.optimize()
        if solution.status != "optimal":
            return float("nan")
        value = solution.objective_value
        if value is None or not math.isfinite(value):
            return float("nan")
        return float(value)


def _percent_change(reference, candidate):
    if not math.isfinite(reference) or not math.isfinite(candidate):
        return "undefined (non-finite growth)"
    if reference == 0.0:
        return "undefined (zero reference)"
    return f"{100 * (candidate - reference) / reference:+.2f}%"


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("candidate", nargs="?", type=pathlib.Path,
                        help="explicit comparison model (SBML XML or XML.gz)")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--verify-candidate", action="store_true",
                       help="verify the retained v9.1.1 asset offline, without writing a table")
    modes.add_argument("--fetch-candidate", action="store_true",
                       help="explicitly create only the missing pinned v9.1.1 comparison asset")
    args = parser.parse_args(argv[1:])
    if args.verify_candidate or args.fetch_candidate:
        if args.candidate is not None:
            parser.error("candidate verification/acquisition cannot target another model path")
        try:
            path, entry = provision_candidate(fetch=args.fetch_candidate)
        except (OSError, ValueError, KeyError, TypeError, EOFError, subprocess.SubprocessError) as error:
            print(f"candidate verification/acquisition failed: {error}", file=sys.stderr)
            return 1
        print(f"verified: {path} ({entry['size_bytes']} bytes; sha256={entry['sha256']})")
        return 0
    if args.candidate is None:
        parser.error("supply an explicit candidate path or a candidate verification/acquisition mode")
    candidate = args.candidate
    if not candidate.is_file():
        print(f"candidate not found: {candidate}", file=sys.stderr)
        return 1
    if candidate.resolve() == (paths.REPO_ROOT / CANDIDATE).resolve():
        try:
            verify_candidate(candidate, candidate_entry())
        except (OSError, ValueError, KeyError, TypeError, EOFError) as error:
            print(f"candidate verification failed: {error}", file=sys.stderr)
            return 1
    vendored = paths.yeast_gem()
    if vendored is None:
        print("no vendored yeast-GEM; set YSTWIN_YEAST_GEM", file=sys.stderr)
        return 1

    old, old_settings = load_model(vendored)
    new, new_settings = load_model(candidate)
    print(f"Reference solver: {old_settings}; candidate solver: {new_settings}")

    rows = []
    for label, a, b in (("reactions", len(old.reactions), len(new.reactions)),
                        ("metabolites", len(old.metabolites), len(new.metabolites)),
                        ("genes", len(old.genes), len(new.genes))):
        rows.append({"check": f"count: {label}", "vendored": a, "candidate": b,
                     "delta": b - a, "breaking": False})

    old_ids = {r.id for r in old.reactions} | {m.id for m in old.metabolites}
    new_ids = {r.id for r in new.reactions} | {m.id for m in new.metabolites}
    missing = [i for i in PINNED if i not in new_ids]
    for identifier, what in PINNED.items():
        present = identifier in new_ids
        rows.append({"check": f"pinned id: {identifier} ({what})",
                     "vendored": identifier in old_ids, "candidate": present,
                     "delta": "", "breaking": not present})

    for bound in GLUCOSE_BOUNDS:
        a, b = _growth(old, bound), _growth(new, bound)
        rows.append({"check": f"growth /h at glucose {bound:g}, O2 unlimited",
                     "vendored": round(a, 6), "candidate": round(b, 6),
                     "delta": _percent_change(a, b), "breaking": False})

    frame = pd.DataFrame(rows)
    out = paths.outputs_dir() / "gem_version_diff.csv"
    frame.to_csv(out, index=False)
    print(frame.to_string(index=False))
    print()
    if missing:
        print(f"BREAKING: {len(missing)} pinned id(s) absent from the candidate: {missing}")
        print("Do not adopt without remapping every one of them.")
    else:
        print("No listed pinned id is missing. Identifier survival alone does not establish compatibility.")
        print("It is still a REGENERATION EVENT: rerun scripts/maintenance_scale.py, "
              "scripts/gem_environment_predictions.py and scripts/audit_claims.py, and "
              "expect every marked GEM number to move.")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
