# Sensor crosstalk, endpoint, August 2026

Do the ER stress sensors respond to oxidative stress, and do the oxidative sensors respond to
ER stress? The reciprocal experiment: UPRE1 and UPRE2 challenged with H2O2, NativeYap1 and
AlteredYap1 challenged with DTT.

**The workbook is not committed and must not be.** `ER&OX-summary.xlsx` carries
`dc:creator` and `cp:lastModifiedBy` naming a private individual, the same reason the
NewProtocol workbooks are absent from this repository. `erox_2026-08_endpoint.tsv` is the
derived measurement table, checked to contain no name and no local path.

## What was run

| | |
| --- | --- |
| reader | Synergy H1, serial 269921 |
| plate | 96-well, double-orbital shake |
| reads | OD600 absorbance, then mCitrine fluorescence at three filter/gain settings |
| **reading type** | **endpoint — `Reading Type: Reader`, NOT kinetic** |
| **temperature** | **24.9 °C recorded, not 30** |
| layout | SPL1–28 in technical triplicate, plus a blank column |
| plates | three biological replicates: `20260821`, `20260823`, `20260826` |
| gain | 75 for the summarised fluorescence channel |

## What is in the TSVs

Three files, all derived by `scripts/extract_crosstalk_workbook.py` and none of them
hand-transcribed.

| file | rows | what it holds |
| --- | ---: | --- |
| `erox_2026-08_endpoint.tsv` | 84 | one row per construct, dose and **biological replicate** |
| `erox_2026-08_on_vs_off_target.tsv` | 28 | each sensor's own stressor beside the wrong one |
| `erox_2026-08_workbook_fold_disagreement.tsv` | 28 | derived folds against the workbook's own |

`od600_mean`, `rfu_mean` and `per_cell` are the workbook's `average` columns; the `_half_sd`
columns are its `std/2`, the **technical** triplicate spread within one plate. `per_cell` is
read from the workbook's `RFU / OD600` section rather than computed, and that is deliberate:
the sheet forms it as the **mean of the per-well ratios**, not the ratio of the means, and
every fold it computes rests on that. They differ in the third decimal.

### This file used to describe one replicate as three

The previous TSV held **replicate 1 alone** while this file called it "the workbook's own
biological-replicate mean across the three plates". The `Summary` sheet lays the three
replicates side by side at columns B, O and AB, and only the first was copied. UPRE1 at 0 mM
reads `0.320333 / 1774.0` in that block, against `0.517333 / 2691` and `0.487667 / 2373.67`
in the other two, whose mean is `0.441778`. Its `std/2` came across too, as though a
within-plate spread were a between-plate one.

It mattered beyond the sample size. On replicate 1 UPRE1 is at **90%** viability at 1 mM
H2O2; on the other two it is at **40%** and **31%**. `score_sensor_crosstalk.py` passed that
dose because the one plate it could see was the one where the cells were alive. The scorer
now forms every fold within its own plate and takes the **worst** plate's viability, not the
mean.

### And the workbook's own fold rows are wrong for two of the eight blocks

Found by checking the extraction against the instrument export rather than against the
workbook. All 28 of replicate 1's per-cell values recompute **exactly** from the blanked
plate in the `1st (raw)` tab. The workbook's fold rows (37-44) are typed constants, not
formulas, and for replicate 1 the two oxidative sensors disagree with the per-cell numbers
sitting directly above them:

| replicate 1, 0.1 mM | per-cell ratio | stated fold row | |
| --- | ---: | ---: | --- |
| UPRE1 | 1.0431 | 1.0431 | agrees |
| UPRE2 | 1.0116 | 1.0116 | agrees |
| NativeYap1 | 0.9566 | 1.0202 | **disagrees** |
| AlteredYap1 | 0.9777 | 1.0256 | **disagrees** |

Replicates 2 and 3 agree with their own per-cell columns exactly. The error propagates into
the workbook's three-replicate aggregate and its "Combining everything" table, moving 12 of
28 crosstalk folds by up to **0.124** — all of them NativeYap1 or AlteredYap1. So the
extraction derives folds from the per-cell values and never reads the fold rows, and
`erox_2026-08_workbook_fold_disagreement.tsv` records every difference so the choice is
auditable rather than asserted.

**This is not a criticism of the experiment.** The measurements are sound and reproduce from
the instrument file. It is a spreadsheet whose derived column was pasted rather than
computed, which is the same class of thing this repository keeps finding in its own prose.

### Replicate 2's raw tab is not blank-subtracted

The workbook says so itself, at `Summary!B2`, verbatim: *"Note that for the second replicate,
the original values are not blanked. Will need to refer to the manually created tables on the
side!!"* The `Summary` blocks **are** those manually created tables, which is why the
extractor reads `Summary` and never the `1st/2nd/3rd (raw)` tabs — except for the integrity
check above, which reads the *blanked* half of the `1st (raw)` grid on purpose.

### The dose ladders differ, and one set of labels is swapped

H2O2 tops out at 4 mM and DTT at 5 mM, so the two halves of each block carry different top
doses and their own dose-label columns. In the fold section the two labels are **swapped**
for replicates 1 and 2 — it writes 5 mM against the H2O2 constructs. Doses in the TSVs come
from the ladder each construct was actually given, not from the row labels.

## The limitation that decides what can be claimed

**This is an endpoint read, so the repository's dilution correction cannot be applied to it.**
`reporter.py::promoter_activity` inverts the reporter ODE using a growth rate taken from an
OD *trace*, and there is no trace here — one timepoint per well. The best available
correction is the per-cell ratio, which removes how many cells there are and does **not**
remove the fact that a stressed culture accumulates reporter for longer per division.

So a per-cell fold change here is a weaker quantity than the folds in
`outputs/late_window_sensitivity.csv`, and the two should not be put in the same table. What
it is strong enough to support is the negative: a ratio that does not move cannot be hiding a
response.

The 24.9 °C is worth carrying too. Every other plate in this project ran at 30, and the
generator's cardinal-temperature model puts growth at 24.9 °C around 80% of its 30 °C value.
Nothing here depends on that, because every comparison is within this plate set.
