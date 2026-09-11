"""Gate pipeline on the real plates: G1 optical quality, then D2 on what survives.

Order matters. D2 asks how much of a reporter is growth rate; that question is
meaningless in a well whose OD is not a quantitative biomass measurement.

The D2 medians printed at the end of each plate carry the spread of the wells they were
taken over, via ``run_d2.well_interval``. The G1 counts above them deliberately do not:
"41 of 87 wells passed" is a census of this plate, not an estimate from a sample of it, and
attaching a sampling interval to a count of the things counted would invent uncertainty
rather than report it.

Usage: python scripts/run_gates.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from ystwin import paths
from ystwin.diagnostics.dilution_confound import detect_blank_wells, plate_dilution_report
from ystwin.gates.g1_optical import OpticalQualityGate, assess_plate
from ystwin.plate import replay
from ystwin.plate.layout import recorded_well_roles
from ystwin.plate.synergy import read_synergy_kinetic
from ystwin.reporter import ReporterKinetics

# The D2 summary here is the same statistic over the same per-well frame that run_d2.py
# reports, so it uses the same estimator rather than a second copy of it that could drift.
from run_d2 import well_interval
from ystwin.readings import RawOD, RawRFU

OUT = paths.outputs_dir()

# 96-well path is ~0.3-0.5 cm, so linearity ends below the cuvette figure of 1.0.
GATE = OpticalQualityGate(od_linear_max=1.0, max_decline_fraction=0.15)


def run(path: pathlib.Path) -> dict | None:
    run_ = read_synergy_kinetic(path)
    # Resolve through raw_channel, not by matching the channel name: a plate read more
    # than once names its blocks "OD600[1]", "mCitrine[2]", which an exact match misses
    # silently rather than loudly.
    try:
        od_block = run_.raw_channel("OD600")
        fl_block = run_.raw_channel("mCitrine")
    except (KeyError, ValueError) as exc:
        # Never return None here. An absent plate and a plate that passed nothing look
        # identical in a report that says nothing about either.
        return {"file": path.name, "note": f"no usable channel pair ({exc})"}
    aligned = run_.aligned(od_block.channel)
    od, reporter = aligned[od_block.channel], fl_block.channel
    rfu = aligned[reporter]
    roles = recorded_well_roles(od, rfu)
    shared = [w for w in od.columns if w in set(rfu.columns)]
    unrecorded = (sorted((set(od.columns) | set(rfu.columns)) - set(roles))
                  if roles is not None else [])
    if roles is not None:
        expected_cultures = {well for well, role in roles.items() if role == "culture"}
        missing_cultures = sorted(expected_cultures - set(shared))
        missing_blanks = sorted(well for well, role in roles.items()
                                if role == "blank" and well not in shared)
        if missing_cultures or missing_blanks:
            return {"file": path.name, "n_cultures": len(expected_cultures),
                    "missing_culture_wells": missing_cultures,
                    "missing_blank_wells": missing_blanks,
                    "note": f"recorded wells lack a paired channel: cultures {missing_cultures}, blanks {missing_blanks}"}
    od, rfu = od[shared], rfu[shared]

    blanks = detect_blank_wells(RawOD(od))
    if not blanks:
        return {"file": path.name, "note": "no media blank; needs an explicit plate map"}
    od_blank = float(od[blanks].to_numpy().mean())
    rfu_bg = float(rfu[blanks].to_numpy().mean())

    cultures = [w for w in od.columns if (roles.get(w) == "culture" if roles is not None else w not in blanks)]
    if not cultures:
        return {"file": path.name, "n_cultures": 0, "blank_wells": blanks,
                "unrecorded_wells": unrecorded, "note": "no culture wells available for G1"}
    g1 = assess_plate(RawOD(od[cultures]), od_blank=od_blank, rfu=RawRFU(rfu[cultures]),
                      rfu_background=rfu_bg, gate=GATE)
    tag = path.stem.replace(" ", "_").replace("&", "and")
    g1.to_csv(OUT / f"g1_{tag}.csv", index=False)

    reasons = (
        g1[~g1.passed].failures.str.split(",").explode().value_counts()
        if (~g1.passed).any() else pd.Series(dtype=int)
    )
    survivors = g1[g1.passed].well.tolist()

    print(f"\n{'=' * 78}\n{path.name}   reporter={reporter}\n{'=' * 78}")
    print(f"  cultures {len(cultures)}   media blanks {blanks}")
    print(f"  G1 PASS: {len(survivors)}/{len(cultures)} ({len(survivors) / len(cultures):.0%})")
    if len(reasons):
        print("  failure reasons:")
        for reason, n in reasons.items():
            print(f"     {reason:<26} {n:>3} wells")
    print(f"  raw OD range across plate: {od[cultures].to_numpy().min():.3f} - "
          f"{od[cultures].to_numpy().max():.3f}")

    d2 = None
    if survivors:
        sub = aligned.loc[:, (slice(None), survivors + blanks)]
        # The resolved block name, not the literal -- same reason as above.
        d2 = plate_dilution_report(sub, od_block.channel, reporter,
                                   kinetics=ReporterKinetics(k_deg=0.0), blank_wells=blanks)
        pw = d2.per_well
        pw.to_csv(OUT / f"d2_g1passed_{tag}.csv", index=False)
        print(f"  -- D2 scored {len(pw)} of {len(survivors)} G1-passing cultures:")
        if not pw.empty:
            r2 = pw.dilution_r2.dropna()
            print(f"     verdicts {pw.verdict.value_counts().to_dict()}")
            if len(r2):
                print(f"     dilution R^2 median {well_interval(r2)}")
            print(f"     negative-activity fraction median "
                  f"{well_interval(pw.negative_activity_fraction, '{:.2f}')}")
            k_deg = pw.implied_min_k_deg.median()
            print(f"     implied min k_deg (1/h) median "
                  f"{well_interval(pw.implied_min_k_deg, '{:.4f}')}"
                  f"  (half-life {np.log(2) / max(k_deg, 1e-9):.1f} h)")
        for excluded in d2.excluded.itertuples():
            print(f"     D2 excluded {excluded.well}: {excluded.reason}")
        print("     bounds resample wells within this plate; one plate is one biological")
        print("     replicate, so they carry no between-plate variance and are not a")
        print("     statement about reproducibility.")
    else:
        print("  -- D2 not run: no culture passed G1.")
    return {"file": path.name, "n_cultures": len(cultures), "n_pass": len(survivors),
            "blank_wells": blanks, "unrecorded_wells": unrecorded,
            "n_d2_attempted": len(survivors),
            "n_d2_scored": len(d2.per_well) if d2 is not None else 0,
            "d2_exclusions": d2.excluded.to_dict("records") if d2 is not None else []}


def _select_exports(exports, manifest):
    required = {"export", "source_set", "source_sha256"}
    if not required <= set(manifest.columns):
        raise ValueError(f"gate source manifest lacks {sorted(required - set(manifest.columns))}")
    declared = {}
    for name, group in manifest.groupby("export", sort=False):
        if group.source_set.nunique(dropna=False) != 1 or group.source_sha256.nunique(dropna=False) != 1:
            raise ValueError(f"conflicting source classification or identity for {name!r}")
        scope, digest = group.source_set.iloc[0], group.source_sha256.iloc[0]
        if not isinstance(scope, str) or not scope or not isinstance(digest, str) or not digest:
            raise ValueError(f"missing source classification or identity for {name!r}")
        declared[name] = (scope, digest)
    selected, excluded, locations, identities = [], [], {}, {}
    for entry in exports:
        path = pathlib.Path(entry)
        if path.name in locations:
            if locations[path.name] != path:
                raise ValueError(f"multiple source locations for registered export {path.name!r}")
            continue
        locations[path.name] = path
        if path.name not in declared:
            raise ValueError(f"unregistered gate source {path.name!r}; declare its source set and identity")
        scope, digest = declared[path.name]
        if scope not in {"newprotocol", "july_2026_07"}:
            excluded.append({"file": path.name, "source_set": scope, "source_sha256": digest,
                             "status": "superseded" if scope == "newprotocol_superseded" else "outside_scope"})
            continue
        if digest in identities:
            raise ValueError(f"{path.name!r} and {identities[digest]!r} declare the same source; not independent replicates")
        identities[digest] = path.name
        selected.append(path)
    return selected, excluded


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--committed", action="store_true")
    args = parser.parse_args(argv)
    manifest = replay.load_manifest()
    # Both plate sets, because the finding this script produces is a *comparison*: the
    # July exports are what G1 was written against, the NewProtocol replicates are what
    # the changed inoculation density produced. Gate only one and the improvement claim
    # has no table behind it.
    sources = [
        (None if args.committed else paths.igem_results(), "the 2026-07 Synergy exports", "YSTWIN_IGEM_RESULTS"),
        (None if args.committed else paths.biosensor_plates(), "the NewProtocol replicates", "YSTWIN_PLATES"),
    ]
    # BOTH plate sets are now committed as text, so the replay covers both and is installed
    # once rather than as a special case for one of them.
    #
    # This loop used to key the swap on `variable == "YSTWIN_PLATES"`, because only the
    # NewProtocol replicates had been exported and the July pair still had to reach a real
    # workbook. That asymmetry is what left `paths.igem_results()` resolving to a sibling
    # directory by convention -- reproducible on one machine -- and it is gone:
    # `plate_readings.py --export` now writes both sets and verifies both value-for-value.
    #
    # `replay.install` dispatches per export on whether the manifest covers the name and
    # falls through to the real reader otherwise, so a machine holding one set of workbooks
    # and not the other gets the raw reader for the first and the committed text for the
    # second, in one pass.
    exports = []
    replayed = []
    for resolved, what, variable in sources:
        found = sorted(resolved.rglob("*.xlsx")) if resolved is not None else []
        if found:
            exports += found
        else:
            replayed.append(what)
    if replayed:
        have = {path.name for path in exports}
        exports += [path for path in replay.install(globals()) if path.name not in have]
        print(f"replaying the committed text in {replay.committed_dir()} for: "
              f"{', '.join(replayed)}")
    exports, records = _select_exports(exports, manifest)
    if not exports:
        # An empty gate report reads like a plate that passed nothing, not like no plate.
        raise SystemExit("no current declared gate sources are available")
    metadata = manifest.drop_duplicates("export").set_index("export")
    for excluded in records:
        print(f"\n{excluded['file']}: {excluded['status']} (source_set={excluded['source_set']})")
    for f in exports:
        try:
            info = run(f)
            if info is None:
                raise ValueError("gate run returned no result or refusal")
            info = {**info, "status": "refused" if "note" in info else "processed"}
            if "note" in info:
                print(f"\n{f.name}: {info['note']}")
        except Exception as exc:  # noqa: BLE001 - survey script
            info = {"file": f.name, "status": "refused", "note": f"{type(exc).__name__}: {exc}"}
            print(f"\n{f.name}: refused ({info['note']})")
        records.append({**info, "source_set": metadata.loc[f.name, "source_set"],
                        "source_sha256": metadata.loc[f.name, "source_sha256"]})
    expected = set(manifest.loc[manifest.source_set.isin(["newprotocol", "july_2026_07"]), "export"])
    for missing in sorted(expected - {path.name for path in exports}):
        records.append({"file": missing, "status": "refused", "note": "declared current gate source was not available",
                        "source_set": metadata.loc[missing, "source_set"],
                        "source_sha256": metadata.loc[missing, "source_sha256"]})
    refused = sum(record["status"] == "refused" for record in records)
    report = {"schema_version": 1, "complete": refused == 0,
              "committed_only": args.committed,
              "independent_unit": "biological plate; file representations are not independent replicates",
              "n_sources_processed": sum(record["status"] == "processed" for record in records),
              "n_sources_refused": refused, "exports": records}
    with (OUT / "gates_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return 2 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
