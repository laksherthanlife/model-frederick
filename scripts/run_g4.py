"""G4: does each biosensor track its independent qPCR anchor?

Anchors: Hac1 for the ER reporters, TRX2 for the oxidative ones, both against UBC.
Reporter response comes from the per-construct dose-response sheets -- the signal as
the team reads it.

Caveat carried into the output: the qPCR (11, 13 Aug) and the biosensor plates
(22 Jul, 3, 4 Aug) are different cultures, so this is paired by condition, not by
culture. A same-culture design would be stronger.

Usage: python scripts/run_g4.py
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.gates.g4_anchor import anchor_agreement
from ystwin.gates.artefact import refuse_partial_rebuild
from ystwin.plate import replay
from ystwin.plate.dose_response import endpoint_response, read_dose_response
from ystwin.qpcr import (
    STRESSOR_FOR_CONSTRUCT,
    delta_delta_cq,
    read_annotated_cq,
    read_positional_cq,
    rt_minus_margin,
)

OUT = paths.outputs_dir()
CONSTRUCTS = ["UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"]
_REPLICATES = (
    ("qpcr_24Jul", "positional", "20260724_ER_Ox_rep1_Cq.xlsx"),
    ("qpcr_11Aug", "annotated",
     "11 August ER and Oxidative Stress/iGEM 11 August_ER&Ox 2nd Biological Replicate -  Quantification Cq Results.xlsx"),
    ("qpcr_13Aug", "annotated",
     "13 August ER and Oxidative Stress/iGEM 13 August_ER&Ox 3rd Biological Replicate -  Quantification Cq Results.xlsx"),
)
"""Canonical replicate order. Read positionally, so it is also the plate-map order."""


def qpcr_files() -> list[tuple[int, str, pathlib.Path, str]]:
    """Index, tag, Cq export and reader for each replicate, from wherever it lives.

    Two locations, not one, and the split is not arbitrary. Replicate 1's archive entries
    were lowercased by the instrument software, so it exists only as the rebuilt copy
    ``scripts/repair_qpcr_export.py`` writes into the repository; the other two open as
    exported and stay wherever the lab keeps them.

    A replicate whose directory is absent is left out rather than pointed at a path that
    cannot exist -- but the ones that remain keep their index from ``_REPLICATES``, so the
    ``replicate`` column in the output means 24 July on every machine rather than
    "whichever export happened to be found first here".
    """
    directories = {"positional": paths.qpcr_dir(), "annotated": paths.qpcr_raw_dir()}
    return [(index, tag, directories[how] / relative, how)
            for index, (tag, how, relative) in enumerate(_REPLICATES)
            if directories[how] is not None]


def load_anchor() -> tuple[pd.DataFrame, pd.DataFrame]:
    folds, qcs = [], []
    for i, tag, path, how in qpcr_files():
        if not path.exists():
            print(f"  {tag}: {path.name} not present; skipped")
            continue
        tidy = read_positional_cq(path) if how == "positional" else read_annotated_cq(path)
        qc = rt_minus_margin(tidy)
        qc["replicate_name"] = tag
        qcs.append(qc)
        out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)
        out["replicate"] = i
        folds.append(out)
    if not folds:
        raise SystemExit(
            "no qPCR replicate could be read, so there is no anchor to gate against. "
            "Set YSTWIN_QPCR_RAW to the Bio-Rad export directory, and YSTWIN_QPCR to the "
            "repaired replicate-1 copy (see scripts/repair_qpcr_export.py)."
        )
    return pd.concat(folds, ignore_index=True), pd.concat(qcs, ignore_index=True)


def load_reporter() -> pd.DataFrame:
    plates = paths.biosensor_plates()
    if plates is None:
        # The workbooks are not committed -- their document properties name a private
        # individual -- so without YSTWIN_PLATES this used to exit and the anchor could
        # not be checked by anyone but the person who ran it. The numbers ARE committed,
        # as text under data/plates, and replaying them reaches the same reader.
        exports = replay.install(globals())
        plates = replay.committed_dir()
    else:
        exports = sorted(plates.glob("*.xlsx"))
    rows = []
    for path in exports:
        try:
            tidy = read_dose_response(path)
        except ValueError:
            continue
        summary = endpoint_response(tidy, fraction=0.25, relative_to_dose=0.0)
        summary["plate"] = path.stem
        rows.append(summary)
    if not rows:
        raise SystemExit(f"no export in {plates} carries a readable dose-response sheet, "
                         "so there is no reporter response to compare the anchor with")
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    # Before anything is written. Two of the three qPCR exports live outside the
    # repository, and on a machine without them the anchor drops from three biological
    # replicates to one -- which turns a measured effect into NaN and a replicate
    # requirement into 0. That is indistinguishable, in the table, from "there is no
    # anchor", and it is a fact about the machine.
    refuse_partial_rebuild(
        OUT / "g4_verdicts.csv",
        found=len(qpcr_files()), expected=len(_REPLICATES), what="qPCR replicates",
        remedy="Set YSTWIN_QPCR_RAW to the Bio-Rad export directory to reach the other "
               "two (see scripts/repair_qpcr_export.py for replicate 1).")

    anchor, qc = load_anchor()
    reporter = load_reporter()
    anchor.to_csv(OUT / "g4_anchor_fold_change.csv", index=False)
    reporter.to_csv(OUT / "g4_reporter_response.csv", index=False)
    qc.to_csv(OUT / "g4_rt_minus_qc.csv", index=False)

    print("=" * 78)
    print("QC gate: no-reverse-transcriptase controls")
    print("=" * 78)
    print(f"  {qc.passed.sum()}/{len(qc)} readings clear a 3-cycle margin\n")
    print(f"  {'construct':<13}{'target':<7}{'pass':>7}{'median margin':>15}")
    for (construct, target), g in qc.groupby(["construct", "target"]):
        flag = "   <-- signal is largely genomic DNA" if g.passed.mean() < 0.5 else ""
        print(f"  {construct:<13}{target:<7}{g.passed.mean():>6.0%}"
              f"{g.margin_cycles.median():>14.2f}{flag}")

    print("\n" + "=" * 78)
    print("Reporter dose response (fold vs 0 mM, mean of last 25% of timecourse)")
    print("  NOTE: UPRE1/UPRE2 doses are DTT; NativeYap1/AlteredYap1 doses are H2O2.")
    print("        The two dose axes are different agents and are not comparable.")
    print("=" * 78)
    pivot = reporter.pivot_table(index="construct", columns="dose_mM", values="response")
    print(pivot.to_string(float_format=lambda v: f"{v:.2f}"))

    print("\n" + "=" * 78)
    print("G4 verdicts")
    print("=" * 78)
    shared_doses = sorted(set(anchor.dose_mM) & set(reporter.dose_mM))
    print(f"  doses shared between anchor and reporter: {shared_doses}\n")

    results = []
    for construct in CONSTRUCTS:
        a = anchor[(anchor.construct == construct) & anchor.dose_mM.isin(shared_doses)]
        r = reporter[(reporter.construct == construct) & reporter.dose_mM.isin(shared_doses)]
        if a.empty or r.empty:
            continue
        response = r.groupby("dose_mM").response.mean()
        frame = pd.DataFrame({
            "replicate": a.replicate.to_numpy(),
            "dose_mM": a.dose_mM.to_numpy(),
            "anchor_fold_change": a.fold_change.to_numpy(),
            "reporter_response": a.dose_mM.map(response).to_numpy(),
            # Held constant: with no variance across doses the gate reports growth as
            # untested rather than as cleared.
            "growth_rate": 0.30,
        })
        target = a.target.iloc[0]
        margin = qc[(qc.construct == construct) & (qc.target == target)]
        result = anchor_agreement(frame, anchor_qc_pass_rate=float(margin.passed.mean()))
        results.append({
            "construct": construct, "stressor": STRESSOR_FOR_CONSTRUCT[construct],
            "anchor": target, "verdict": result.verdict,
            "anchor_effect": result.anchor_effect, "anchor_spread": result.anchor_spread,
            "reporter_anchor_r": result.reporter_anchor_r,
            "rt_minus_pass_rate": margin.passed.mean(),
            "replicates_needed": result.replicates_needed,
        })
        gate = "" if margin.passed.mean() >= 0.5 else "  [anchor fails RT- QC]"
        agent = STRESSOR_FOR_CONSTRUCT[construct]
        print(f"  {construct} vs {target}   stressor: {agent}{gate}")
        print(f"    {result.verdict}: {result.reason}")
        print(f"    reporter fold at {shared_doses}: "
              f"{[round(response.get(d, float('nan')), 2) for d in shared_doses]}")
        print(f"    anchor fold per replicate at {max(shared_doses)} mM: "
              f"{[round(v, 2) for v in a[a.dose_mM == max(shared_doses)].fold_change]}\n")

    pd.DataFrame(results).to_csv(OUT / "g4_verdicts.csv", index=False)
    print("\n  Growth rate per dose is now available (the plate map was recovered and\n"
          "  confirmed against the logbook), but the gate never reaches the growth test:\n"
          "  every anchor stops at the QC or power check first. See\n"
          "  scripts/run_sensor_characterisation.py for the growth-corrected responses.")


if __name__ == "__main__":
    main()
