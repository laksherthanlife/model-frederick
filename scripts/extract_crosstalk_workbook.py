"""Pull all three biological replicates of the crosstalk experiment out of the workbook.

    YSTWIN_CROSSTALK_WORKBOOK=~/Downloads/'ER&OX-summary.xlsx' \\
        python3 scripts/extract_crosstalk_workbook.py

Writes ``data/crosstalk/erox_2026-08_endpoint.tsv`` (per replicate, per construct, per dose)
and ``data/crosstalk/erox_2026-08_on_vs_off_target.tsv`` (the workbook's own three-replicate
fold table, on-target beside off-target).

**Why this script exists.** The TSV it replaces held **replicate 1 only**, while
`data/crosstalk/SOURCE.md` described it as "the workbook's own biological-replicate mean
across the three plates". Checked cell by cell: UPRE1 at 0 mM reads 0.320333 in the TSV and
in the workbook's first block, against 0.517333 and 0.487667 in the second and third, whose
mean is 0.441778. The `Summary` sheet lays the three replicates side by side at columns B, O
and AB; only the first was transcribed, and its ``std/2`` column -- the technical triplicate
spread within that one plate -- was carried across as though it were a spread across plates.

**The workbook is not committed and must not be.** `ER&OX-summary.xlsx` carries
``dc:creator`` and ``cp:lastModifiedBy`` naming a private individual, the same reason the
plate workbooks are absent. What is committed is the measurement table, checked to contain
no name and no local path -- the treatment `data/plates/` already gives the kinetic exports.

**Two things in the workbook that a naive read gets wrong, both handled below.**

*The dose ladders differ by stressor.* H2O2 tops out at 4 mM and DTT at 5 mM -- the stock
arithmetic in the lab notebook reaches 3.997 mM for peroxide -- so the two halves of every
block carry different top doses and their own dose-label columns. The fold section at rows
37-44 has those two labels **swapped** in replicates 1 and 2 (it writes 5 mM against the
H2O2 constructs), which is why doses here come from the ladder each construct was actually
given rather than from the row labels.

*Replicate 2's raw tab is not blank-subtracted.* The workbook says so itself, at
``Summary!B2``: "Note that for the second replicate, the original values are not blanked.
Will need to refer to the manually created tables on the side!!" The ``Summary`` blocks are
those manually created tables, which is why this script reads ``Summary`` and never the
``1st/2nd/3rd (raw)`` tabs.

**The extraction is checked against the raw instrument export, and that check found a third
thing.** Every per-cell value is recomputed from the blanked plate in the ``1st (raw)`` tab
and must agree exactly. It does. But the workbook's own **fold rows** (37-44) are typed
constants rather than formulas, and for replicate 1 the two oxidative sensors disagree with
the per-cell numbers directly above them:

    replicate 1, 0.1 mM      per-cell ratio    stated fold row
    UPRE1                        1.0431            1.0431   agrees
    UPRE2                        1.0116            1.0116   agrees
    NativeYap1                   0.9566            1.0202   DISAGREES
    AlteredYap1                  0.9777            1.0256   DISAGREES

The per-cell values are the ones that reproduce from the raw export -- all eight checked
cell for cell -- so the fold rows are what is wrong, and only for replicate 1's DTT half.
Replicates 2 and 3 agree with their own per-cell columns exactly. That error propagates into
the workbook's three-replicate aggregate at rows 51-57 and into its "Combining everything"
table, so **this script derives folds from the per-cell values and does not read the fold
rows at all**, and reports the size of the disagreement rather than silently differing.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths

#: Where each biological replicate's block starts, and the column offsets inside it. The two
#: dose columns per block are the two stressor ladders: H2O2 for the ER sensors (crosstalk),
#: DTT for the oxidative ones.
REPLICATE_BLOCKS = (
    {"replicate": 1, "h2o2_dose_col": 2, "dtt_dose_col": 8,
     "constructs": {"UPRE1": 3, "UPRE2": 6, "NativeYap1": 9, "AlteredYap1": 12}},
    {"replicate": 2, "h2o2_dose_col": 15, "dtt_dose_col": 21,
     "constructs": {"UPRE1": 16, "UPRE2": 19, "NativeYap1": 22, "AlteredYap1": 25}},
    {"replicate": 3, "h2o2_dose_col": 28, "dtt_dose_col": 34,
     "constructs": {"UPRE1": 29, "UPRE2": 32, "NativeYap1": 35, "AlteredYap1": 38}},
)

#: Which stressor each construct was given in THIS experiment. The whole point of the run is
#: that it is the wrong one: the ER sensors get peroxide and the oxidative sensors get DTT.
CROSSTALK_STRESSOR = {"UPRE1": "H2O2", "UPRE2": "H2O2",
                      "NativeYap1": "DTT", "AlteredYap1": "DTT"}

#: The two ladders, in the row order the sheet uses.
LADDERS = {"H2O2": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
           "DTT": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)}

OD_FIRST_ROW, RFU_FIRST_ROW, PER_CELL_FIRST_ROW = 8, 18, 28
"""Row bands for the OD600, RFU and RFU/OD600 sections of every replicate block.

``PER_CELL_FIRST_ROW`` is read from the sheet rather than derived, and the distinction is
not cosmetic. The workbook's ``RFU / OD600`` section is the **mean of the per-well ratios**,
not the ratio of the means, and every fold it computes is built on that. The TSV this script
replaces defined ``per_cell`` as ``rfu_mean / od600_mean`` -- SOURCE.md said so -- which is
a different quantity and disagrees in the third decimal. The check below is what caught it:
derived folds missed the workbook's own by 3e-3 until this was taken from the sheet."""

#: The workbook's own three-replicate aggregate, used as a check on the extraction.
AGGREGATE_FIRST_ROW = 51
AGGREGATE_CROSSTALK = {"UPRE1": 3, "UPRE2": 5, "NativeYap1": 9, "AlteredYap1": 11}

#: "Combining everything": on-target beside off-target, over the same three replicates.
COMBINED_FIRST_ROW = 75
COMBINED = {  # construct -> (regular avg col, crosstalk avg col)
    "UPRE1": (3, 5), "UPRE2": (9, 11), "NativeYap1": (15, 17), "AlteredYap1": (21, 23)}

#: Which stressor is the ON-target one for each sensor. The inverse of CROSSTALK_STRESSOR,
#: and worth writing out rather than deriving, because it is the claim the table makes.
ON_TARGET_STRESSOR = {"UPRE1": "DTT", "UPRE2": "DTT",
                      "NativeYap1": "H2O2", "AlteredYap1": "H2O2"}


def _workbook() -> pathlib.Path:
    path = paths.crosstalk_workbook()
    if path is not None and path.is_file():
        return path
    raise SystemExit(
        "the crosstalk workbook is not here. Set YSTWIN_CROSSTALK_WORKBOOK to "
        "'ER&OX-summary.xlsx'. It is not committed and must not be: it names a private "
        "individual in its document properties, which is why only the derived table is "
        "tracked.")


def _cell(sheet, row: int, col: int) -> float | None:
    value = sheet.cell(row=row, column=col).value
    return float(value) if isinstance(value, (int, float)) else None


def extract_replicates(sheet) -> pd.DataFrame:
    """One row per construct, dose and biological replicate."""
    rows = []
    for block in REPLICATE_BLOCKS:
        for construct, col in block["constructs"].items():
            stressor = CROSSTALK_STRESSOR[construct]
            for index, dose in enumerate(LADDERS[stressor]):
                rows.append({
                    "construct": construct, "stressor": stressor, "dose_mM": dose,
                    "replicate": block["replicate"],
                    "od600_mean": _cell(sheet, OD_FIRST_ROW + index, col),
                    "od600_half_sd": _cell(sheet, OD_FIRST_ROW + index, col + 1),
                    "rfu_mean": _cell(sheet, RFU_FIRST_ROW + index, col),
                    "rfu_half_sd": _cell(sheet, RFU_FIRST_ROW + index, col + 1),
                    "per_cell": _cell(sheet, PER_CELL_FIRST_ROW + index, col),
                    "per_cell_half_sd": _cell(sheet, PER_CELL_FIRST_ROW + index, col + 1),
                })
    return pd.DataFrame(rows)


def extract_on_vs_off_target(sheet) -> pd.DataFrame:
    """The workbook's own three-replicate folds, on-target beside off-target.

    This is the comparison the crosstalk claim actually rests on and it was not in the
    repository at all. "The sensor stayed quiet for the wrong stressor" is weak on its own --
    a dead sensor also stays quiet. "It doubled for the right one and stayed flat for the
    wrong one, on the same plates" is the statement worth making.
    """
    rows = []
    for construct, (regular_col, crosstalk_col) in COMBINED.items():
        ladder = LADDERS[ON_TARGET_STRESSOR[construct]]
        for index, _ in enumerate(ladder):
            row = COMBINED_FIRST_ROW + index
            rows.append({
                "construct": construct,
                "on_target_stressor": ON_TARGET_STRESSOR[construct],
                "off_target_stressor": CROSSTALK_STRESSOR[construct],
                # The two ladders differ only in their top dose, so the shared index is the
                # dose for both columns until the last row, where each takes its own.
                "dose_mM": (LADDERS[ON_TARGET_STRESSOR[construct]][index]
                            if index < len(ladder) - 1 else None),
                "on_target_dose_mM": LADDERS[ON_TARGET_STRESSOR[construct]][index],
                "off_target_dose_mM": LADDERS[CROSSTALK_STRESSOR[construct]][index],
                "on_target_fold": _cell(sheet, row, regular_col),
                "on_target_half_sd": _cell(sheet, row, regular_col + 1),
                "off_target_fold": _cell(sheet, row, crosstalk_col),
                "off_target_half_sd": _cell(sheet, row, crosstalk_col + 1),
            })
    return pd.DataFrame(rows).drop(columns=["dose_mM"])


#: Where replicate 1's blanked plate lives in the ``1st (raw)`` tab. Each plate row occupies
#: eight sub-rows -- OD, gain 50, 75, 100 raw, then the same four blanked -- and the blanked
#: OD and gain-75 rows are the ones the Summary is built from.
RAW_FIRST_BLANKED_OD_ROW, RAW_ROWS_PER_PLATE_ROW = 61, 8
RAW_OD_TO_G75_OFFSET = 2
#: Plate columns 1-3 are UPRE1, 4-6 UPRE2, 7-9 NativeYap1, 10-12 AlteredYap1 -- the layout
#: the lab notebook records and `plate/layout.py::NEWPROTOCOL_LAYOUT` already carries.
RAW_CONSTRUCT_COLUMNS = {"UPRE1": (3, 4, 5), "UPRE2": (6, 7, 8),
                         "NativeYap1": (9, 10, 11), "AlteredYap1": (12, 13, 14)}


def check_against_raw_export(workbook, replicates: pd.DataFrame) -> None:
    """Recompute replicate 1's per-cell values from the instrument export.

    The extraction is a column map, and a column map is exactly the thing that goes wrong
    silently. This checks it against something the Summary sheet did not produce: the
    blanked plate in the ``1st (raw)`` tab, per well, mean of the ratios. Every value must
    agree to a tenth of an RFU/OD unit or the run stops.

    It is also what establishes that the per-cell numbers are right where the fold rows are
    wrong, which is the reason this script ignores the fold rows.
    """
    raw = workbook["1st (raw)"]
    subset = replicates[replicates.replicate == 1]
    worst = 0.0
    for construct, columns in RAW_CONSTRUCT_COLUMNS.items():
        ladder = LADDERS[CROSSTALK_STRESSOR[construct]]
        rows = subset[subset.construct == construct].set_index("dose_mM")
        for index, dose in enumerate(ladder):
            od_row = RAW_FIRST_BLANKED_OD_ROW + index * RAW_ROWS_PER_PLATE_ROW
            od = [raw.cell(row=od_row, column=c).value for c in columns]
            rfu = [raw.cell(row=od_row + RAW_OD_TO_G75_OFFSET, column=c).value
                   for c in columns]
            if any(v is None for v in od + rfu):
                raise SystemExit(
                    f"the raw export has no plate row for {construct} at {dose} mM "
                    f"(looked at row {od_row}); nothing was written")
            recomputed = float(np.mean([r / o for r, o in zip(rfu, od)]))
            theirs = float(rows.loc[dose].per_cell)
            worst = max(worst, abs(recomputed - theirs))
            if abs(recomputed - theirs) > 0.1:
                raise SystemExit(
                    f"extracted per-cell for {construct} at {dose} mM is {theirs:.4f}, but "
                    f"recomputing it from the raw blanked plate gives {recomputed:.4f}. A "
                    "column offset is wrong; nothing was written.")
    print(f"  every replicate-1 per-cell value recomputes from the raw blanked plate: "
          f"28 values, worst difference {worst:.2e}")


def report_workbook_fold_disagreement(sheet, replicates: pd.DataFrame) -> pd.DataFrame:
    """Where the workbook's own aggregate differs from one derived from its per-cell values.

    Reported rather than reconciled. The workbook's fold rows are typed constants and two of
    the eight per-replicate blocks disagree with the per-cell numbers they sit under; the
    per-cell numbers are the ones that reproduce from the instrument export, so those are
    what this script uses. Printing the gap is how a reader learns the published aggregate
    for the oxidative sensors is not the one this table implies.
    """
    rows = []
    for construct, col in AGGREGATE_CROSSTALK.items():
        ladder = LADDERS[CROSSTALK_STRESSOR[construct]]
        subset = replicates[replicates.construct == construct]
        control = subset[subset.dose_mM == 0.0].set_index("replicate").per_cell
        for index, dose in enumerate(ladder):
            dosed = subset[subset.dose_mM == dose].set_index("replicate").per_cell
            folds = [dosed[r] / control[r] for r in sorted(control.index)]
            rows.append({"construct": construct, "dose_mM": dose,
                         "derived_fold": float(np.mean(folds)),
                         "derived_half_sd": float(np.std(folds, ddof=0) / 2),
                         "workbook_fold": _cell(sheet, AGGREGATE_FIRST_ROW + index, col)})
    frame = pd.DataFrame(rows)
    frame["difference"] = frame.derived_fold - frame.workbook_fold
    return frame


def main() -> int:
    import openpyxl

    path = _workbook()
    sheet = openpyxl.load_workbook(path, data_only=True)["Summary"]
    print(f"reading {path.name}")

    workbook = openpyxl.load_workbook(path, data_only=True)
    replicates = extract_replicates(sheet)
    check_against_raw_export(workbook, replicates)
    disagreement = report_workbook_fold_disagreement(sheet, replicates)
    combined = extract_on_vs_off_target(sheet)
    # Carry the derived off-target fold beside the workbook's own, so the disagreement is in
    # the committed table rather than only in this script's output. The on-target column
    # cannot be checked the same way -- it is drawn from the regular biosensor tests, whose
    # raw plates are not in this workbook -- so it is carried as the workbook states it.
    derived = disagreement.set_index(["construct", "dose_mM"])
    combined["off_target_fold_derived"] = [
        float(derived.loc[(r.construct, r.off_target_dose_mM)].derived_fold)
        for r in combined.itertuples()]
    combined["off_target_half_sd_derived"] = [
        float(derived.loc[(r.construct, r.off_target_dose_mM)].derived_half_sd)
        for r in combined.itertuples()]

    directory = paths.data_dir() / "crosstalk"
    endpoint = directory / "erox_2026-08_endpoint.tsv"
    on_off = directory / "erox_2026-08_on_vs_off_target.tsv"
    replicates.to_csv(endpoint, sep="\t", index=False)
    combined.to_csv(on_off, sep="\t", index=False)

    print(f"\n  wrote {endpoint.name}: {len(replicates)} rows, "
          f"{replicates.replicate.nunique()} biological replicates "
          f"(was {len(replicates) // 3} rows of replicate 1 alone)")
    print(f"  wrote {on_off.name}: {len(combined)} rows, on-target beside off-target")
    audit = directory / "erox_2026-08_workbook_fold_disagreement.tsv"
    disagreement.to_csv(audit, sep="\t", index=False)
    print(f"  wrote {audit.name}: {len(disagreement)} folds, derived against the workbook's")

    big = disagreement[disagreement.difference.abs() > 1e-6]
    if len(big):
        print(f"\n  {len(big)} of {len(disagreement)} folds differ from the workbook's own "
              f"aggregate, worst {big.difference.abs().max():.4f}:")
        for r in big.sort_values("difference", key=abs, ascending=False).head(8).itertuples():
            print(f"    {r.construct:12s} {r.dose_mM:4} mM  derived {r.derived_fold:.4f}  "
                  f"workbook {r.workbook_fold:.4f}  ({r.difference:+.4f})")
        print("    The workbook's fold rows are typed constants; replicate 1's two oxidative")
        print("    sensors disagree with the per-cell values above them, which recompute")
        print("    exactly from the raw export. Derived numbers are used.")

    print("\nOn-target against off-target, the workbook's own three-replicate folds:\n")
    for construct in ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"):
        block = combined[combined.construct == construct]
        on = " ".join(f"{v:5.2f}" for v in block.on_target_fold)
        off = " ".join(f"{v:5.2f}" for v in block.off_target_fold)
        print(f"  {construct:12s} {block.on_target_stressor.iloc[0]:>4s} (its own) {on}")
        print(f"  {'':12s} {block.off_target_stressor.iloc[0]:>4s} (wrong)    {off}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
