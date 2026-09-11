"""Fetch the two public datasets the external validation needs, and nothing else.

Both are downloaded raw and cached under ``data/external/``; neither is tracked. The
split between them is the whole point of the design:

``gasch2000/``
    GEO GSE18, the Gasch et al. 2000 environmental stress compendium. Supplies the
    *measurements* -- one log2 ratio per gene per array, across nine stressors.

``sgd/regulation/``
    SGD's curated transcription-factor regulation records, one JSON per factor. Supplies
    the *labels* -- which genes belong to which module's regulon. Kept separate from the
    measurements deliberately: if the module definitions came from expression data they
    would be partly Gasch's own clustering, and the validation would be scoring the panel
    against the paper that parameterised it. SGD records carry an evidence type and a
    reference PMID per row, so the analysis can restrict to binding evidence and audit
    that Gasch 2000 contributes nothing. See docs/EXTERNAL_VALIDATION.md.

Usage: python scripts/fetch_external_validation_data.py [--skip-gasch] [--skip-sgd]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
GASCH = ROOT / "data" / "external" / "gasch2000"
SGD = ROOT / "data" / "external" / "sgd" / "regulation"

# GSE18 was printed on nine array batches and GEO keeps one series matrix per batch. All
# nine are needed: the stressors do not line up with the batches, so dropping one loses
# timepoints out of the middle of a time course rather than a whole condition.
PLATFORMS = ("51", "52", "53", "54", "55", "56", "57", "63", "64")
# GPL63 carries one array and has no annotation table published, so its probes cannot be
# mapped to genes. Listed here rather than silently omitted above.
NO_ANNOTATION = ("63",)

GEO_MATRIX = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE18/matrix/"
              "GSE18-GPL{gpl}_series_matrix.txt.gz")
GEO_ANNOT = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL{gpl}/annot/"
             "GPL{gpl}.annot.gz")
SGD_REGULATION = "https://www.yeastgenome.org/backend/locus/{gene}/regulation_details"

# Every transcription factor named by a MODULES entry in generator/stress_panel.py. ACE1
# is requested as CUP2, which is the name SGD and the GEO platform tables both use.
FACTORS = (
    "MSN2", "MSN4", "HAC1", "YAP1", "SKN7", "HSF1", "SKO1", "HOT1", "RPN4", "AFT1",
    "AFT2", "RLM1", "RFX1", "CRZ1", "ADR1", "CAT8", "GLN3", "GAT1", "ROX1", "UPC2",
    "HAP1", "CUP2", "MAC1", "ZAP1", "MET4", "MET31", "MET32", "PDR1", "PDR3", "RTG1",
    "RTG3", "RIM101",
)


def _get(url: str, destination: pathlib.Path, timeout: float = 300.0) -> bool:
    """Fetch one URL to one path, skipping what is already on disk."""
    if destination.exists() and destination.stat().st_size > 0:
        print(f"  have {destination.name}")
        return True
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            payload = response.read()
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        print(f"  FAILED {destination.name}: {error}", file=sys.stderr)
        return False
    destination.write_bytes(payload)
    print(f"  got  {destination.name} ({len(payload):,} B)")
    return True


def fetch_gasch() -> bool:
    print(f"GSE18 -> {GASCH}")
    ok = True
    for gpl in PLATFORMS:
        ok &= _get(GEO_MATRIX.format(gpl=gpl), GASCH / f"GSE18-GPL{gpl}_series_matrix.txt.gz")
    for gpl in PLATFORMS:
        if gpl in NO_ANNOTATION:
            continue
        ok &= _get(GEO_ANNOT.format(gpl=gpl), GASCH / f"GPL{gpl}.annot.gz")
    return ok


def fetch_sgd() -> bool:
    print(f"SGD regulation -> {SGD}")
    ok = True
    for gene in FACTORS:
        destination = SGD / f"{gene}.json"
        if destination.exists() and destination.stat().st_size > 0:
            print(f"  have {destination.name}")
            continue
        got = _get(SGD_REGULATION.format(gene=gene), destination)
        ok &= got
        if got:
            # Confirm it parsed as a record list before the analysis depends on it.
            try:
                records = json.loads(destination.read_text())
            except json.JSONDecodeError as error:
                print(f"  FAILED {gene}: not JSON ({error})", file=sys.stderr)
                destination.unlink()
                ok = False
                continue
            print(f"       {gene}: {len(records)} regulation records")
        time.sleep(0.5)  # courtesy to a public API with no documented rate limit
    return ok


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-gasch", action="store_true")
    parser.add_argument("--skip-sgd", action="store_true")
    args = parser.parse_args(argv)

    ok = True
    if not args.skip_gasch:
        ok &= fetch_gasch()
    if not args.skip_sgd:
        ok &= fetch_sgd()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
