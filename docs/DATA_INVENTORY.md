# What the biosensor plates actually measured

The README describes the NewProtocol dataset as "7 DTT doses x 3 technical replicates x 4
constructs, 4 biological replicates". That is true of the **optical density** channel. It
is not true of the **fluorescence** channel, which is the one every reported result rests
on, and the difference decides how much of the dose-response table is supportable.

Everything below was read off the exports themselves, not the logbook.

## Fluorescence coverage, per plate

**Current answer first: every plate carries all four constructs and its own media blank, and
the panel is n=4.** The two tables immediately below are the original mis-read, kept because
this document is where the corrections are recorded. What the files hold is in *Where row H's
blanks actually are* and in the coverage table in `docs/FINDINGS.md`.

| Export | OD wells | Fluorescence wells | Constructs with fluorescence | Blank in fluorescence? |
| --- | ---: | ---: | --- | --- |
| `20260722_..._ANALYSED` | 87 | **87** | all four | yes (H1-H3) |
| `20260728_..._Replicate2` | 87 | **21** | **UPRE1 only** | no |
| `20260803_..._Replicate3` | 87 | **87** | all four | yes (H1-H3) |
| `20260804_..._Replicate4` | 21 | **21** | **UPRE1 only** | no |

Replicates 2 and 4 read fluorescence on columns 1-3 rows A-G only. Under the recovered and
logbook-confirmed layout (`plate/layout.py::NEWPROTOCOL_LAYOUT`) columns 1-3 are UPRE1. So
those two plates carry a full 7-dose UPRE1 ladder in triplicate and **nothing at all** for
UPRE2, NativeYap1 or AlteredYap1.

The biological replication actually available is therefore:

| Construct | Plates carrying fluorescence | Plates usable for a fold change |
| --- | ---: | ---: |
| UPRE1 | 4 | **2** |
| UPRE2 | 2 | **2** |
| NativeYap1 | 2 | **2** |
| AlteredYap1 | 2 | **2** |

**This table is superseded twice over: every construct is n=4, and both of the reasons
given below for n=2 were wrong.** It read n=2 until 2026-08-29 and n=3 until 2026-08-30.

<!-- audit:retracted every sensor is effectively down to two -->
<!-- audit:retracted neither has a blank well -->

Those two markers are load-bearing rather than decorative. This retraction was written, and
`README.md` went on asserting both halves of what it retracts for weeks — its section 2 kept
the n=2 table and kept blaming missing blanks, three paragraphs below its own section 1
saying n=3, and above a shipped table carrying `n_plates=3` in all 144 rows. Nothing could
see it: `audit_claims.py` guarded one test count and 64 marked numbers, and an unmarked table
was invisible to both. `check_retracted_claims` now fails if either phrase appears in any
document except this one and `docs/superseded/`, so the next retraction holds without anyone
remembering to sweep for it.

20260728 and 20260804 were said to carry no blank well in the fluorescence
channel. Both carry all 87 wells, blanks included. They write the plate down the page
instead of across it -- one block per group of plate columns, only the first carrying the
Time column -- so the reader stopped after 21 wells, and the blanks live in a later
block. 20260728 is recovered and is now analysable; its recorded blanks (H4-H6) are
exactly the wells detection finds in it.

20260804 used to stay out because every xlsx copy of it is already blank-subtracted, its
H1-H3 reading +-0.002. There is no background to recover *from the exports*.

**The panel is n=4 as of 2026-08-30, and the background came out of the instrument file.**
`plate/gen5.py` reads the Gen5 `.xpt` binary directly, and in
`20260804_ER_Oxidative_Replicate4.xpt` wells H1-H3 read **0.087/0.090/0.088 OD and
347/349/354 RFU** at the first timepoint, against **-3/-1/+4 RFU** in the committed export.
The archive is now committed as text alongside the subtracted export it was exported to --
`20260804_ER_Oxidative_Replicate4__{OD600,mCitrine}.csv`, 96 wells against the export's 87,
in `data/plates/manifest.csv` with the archive's own sha256.

**The arithmetic between the two committed copies is checked, not asserted, and it is now
checkable from the checkout alone.** Subtracting the per-timepoint mean of H1:H3 from the
archive reproduces the subtracted export: in fluorescence over all 2,175 readings every
disagreement is either zero or exactly a third, which is what rounding the mean of three
integers back to an integer must produce and nothing else would; in absorbance the worst
disagreement is 0.0005, the arithmetic maximum for a three-decimal export of a four-decimal
reading. `tests/test_gen5_xpt.py` runs that against the binary and
`tests/test_xpt_export.py` runs it against the committed text with no `.xpt` present.

The parser is proven rather than asserted more broadly too: it reproduces the committed text
of every plate that has both a `.xpt` and a CSV -- fluorescence **exactly**, absorbance to
5e-4. Elapsed times match string for string, and the archive's time axis is decoded through
`H:MM:SS` rather than divided out of the stored milliseconds, because those are the same
real number and different floats at 16 of this plate's 50 timestamps.

Read off the exports directly:

| Export | Fluorescence wells | Blank in the fluorescence channel? |
| --- | ---: | --- |
| `20260722_..._ANALYSED` | 87 | yes (H1-H3) |
| `20260728_..._Replicate2` | 21 | **no** |
| `20260803_..._Replicate3` | 87 | yes (H1-H3) |
| `20260804_..._Replicate4` | 21 | **no** |

`scripts/run_sensor_characterisation.py` was right to drop replicates 2 and 4, and its
stated reason ("no blanks in the shared wells") is correct. But the reason it gives is the
proximate one; the real cause is that the reader was configured to read three columns.

### Where row H's blanks actually are, per plate

The logbook records blank positions per plate and its entries for 20260803 and 20260804
were both copies of 20260728's. The instrument files settle it, because they carry all 96
wells while the exports carry 87 or fewer, and because the reporter channel separates medium
from water where absorbance cannot. Every well of row H that is not a culture sits at the
density floor, 0.086-0.102 OD, on all four plates; only the ones holding medium fluoresce,
and they do so by roughly a factor of six over the instrument's dark count.

Means over the whole run, from the `.xpt`:

| Plate | H1-H3 RFU | H4-H6 RFU | H7-H12 RFU | media blanks |
| --- | ---: | ---: | ---: | --- |
| 20260722 | 406 | 50 | 50 | H1-H3 |
| 20260728 | 1493 — **cultures**, OD 0.117-0.642 | 843 | 51 | H4-H6 |
| 20260803 | 402 | 49 | 47 | H1-H3 |
| 20260804 | 341 | 56 | 57 | **H1-H3**, not the recorded H4-H6 |

So 20260728's logbook entry is right, and independently corroborated: its H1-H3 grew.
20260803's was corrected against its export in an earlier pass. 20260804's is corrected here
on the same evidence plus one piece the others did not have -- Gen5's own blank subtraction
on that plate is the mean of H1:H3, reproduced to the last bit above. Committing the archive
made the correction necessary rather than merely tidy: the archive carries H4-H6, so the
stale entry would have selected them and subtracted 56 RFU where 341 was right.

### What still cannot be added back: the UPRE1-only halves of replicates 2 and 4

Recovering the full 87-well read of 20260728 and the 96-well archive of 20260804 makes both
plates fully usable, so the bounded-background sweep that used to be needed for them is not.
What remains unmeasured on every plate is not the media blank but the *per-biomass*
autofluorescence, which the media blank does not touch -- see the next section.
`analysis/uncertainty.py` supports the bounded shape of statement; nothing needs it here.

## The consequence for the reported dose-response table

**Current state: every row rests on four biological replicates, and seven of the twenty
estimable folds clear "no change" at the canonical window.** `README.md` section 1 and
`docs/FINDINGS.md` carry the current numbers, pinned cell by cell to
`outputs/late_window_sensitivity.csv`.

The n=2 era table, and why the fold estimator moved twice before it settled -- pooled ratio,
then plate-paired, then the total-signal inversion in which the growth rate cancels
algebraically -- is in
[`superseded/dose-response-folds.md`](superseded/dose-response-folds.md). It is not repeated
here. Every row in it was INCONCLUSIVE, because `analysis/uncertainty.py::fold_change`
refused every interval below three biological replicates.

One row from that era is worth carrying forward as a caution. `AlteredYap1` at 1.0 mM H2O2
now has four readings: 1.55 plate-paired at n=2, **1.44 and then 1.49 at n=3** — the same
three plates before and after the estimator's smoothing window was measured rather than
inherited — and **1.00 at n=4** with an interval of [0.24, 2.50]. Both n=3 numbers are quoted
on purpose: a row whose point estimate moves 0.05 on a change to the estimator and 0.49 on a
change to the plates is a row whose point estimate was never the thing to read. It was never
a call at any n, and the direction it moved as replicates were added is the direction
[`superseded/heldout-interpolation.md`](superseded/heldout-interpolation.md) predicted for an
effect size first measured at the smallest n an estimator will accept. The oxidative response
that *did* survive the fourth plate is at a different dose -- 0.5 mM H2O2, fold 1.37,
interval [1.06, 1.70] -- and it is the first per-dose call either oxidative sensor has
produced.

INCONCLUSIVE is not the same as refuted, and the majority of rows are still INCONCLUSIVE at
n=4. The mechanism behind the point estimates -- that a stable reporter accumulates as growth
slows -- is sound and independently supported. What four plates bought is that for seven of
them the evidence that they are distinguishable from 1.0 now exists, and for one it survives
being corrected for all twenty.

## Autofluorescence has never been measured, on any plate

`observation.py` writes the measurement model as

    RFU = background + [autofluorescence + gain * R] * biomass * exp(-eps * carotenoid)

and its docstring insists that autofluorescence be measured "on an isogenic reporter-free
strain, not by subtracting a guess". Two such strains were plated:

- `20260807_BY4741_ER_Oxidative.xlsx` -- 33 wells, **OD600 channel only**
- `20260814_BY4741_0mM_OD600 Data.xlsx` -- 6 wells, **OD600 channel only**

**Corrected 2026-08-30, and this passage was wrong twice.** It said neither was read in the
mCitrine channel, and it said two strains were plated.

`20260807_BY4741_ER_Oxidative.xpt` -- the raw instrument file, which `plate/gen5.py` now
reads -- carries `mCitrine:480,530` for all 96 wells across 25 timepoints, under
`4h-10min_mCitrine_OD600.prt`: the same protocol file, gain, interval and optics as
replicates 2, 3 and 4. **The measurement was taken. The `.xlsx` export dropped the channel.**
And the 2026-08-14 file is not a second plating: its six wells, 25 elapsed times and
temperature series match the 2026-08-07 run exactly. One plate, exported twice.

**Measured 2026-08-30 by `scripts/measure_autofluorescence.py`: a = +18.2 RFU/OD, 95% CI
[-6.1, +42.4], which includes zero.** As a fraction of a reporter well that is 0.42% with an
upper bound of 0.97%, against a sensitivity sweep whose lowest verdict-changing fraction is
5% -- so every call in the panel survives it. The first pass below is superseded: it took a
median of `RFU/OD`, which measures `a + background/OD` and is mostly background at these
optical densities. The 2026-08-29 report that read the recovered channel and then
concluded the sweep had to stay is in
[`superseded/autofluorescence-first-pass.md`](superseded/autofluorescence-first-pass.md); it
prints 239 RFU/OD for the ratio this section quotes as 249, over a different well set.

A first pass over it put the reporter-free per-cell signal at a median of ~249 RFU/OD
against ~4354 and ~3649 for the reporter constructs, of order 6% of the median reporter
well. That is a size, not a measurement -- it splices two plates, and the BY4741 plate has
no medium-only well that fluoresces, so the subtraction can only make 249 an over-estimate.
`docs/research/XPT_INVENTORY.md` §2.2 carries the evidence and the layout; deriving `a` per
construct against its own 0 mM well is what would turn it into a measurement, and that is
now a data-reduction task rather than a plate to run.

What the analysis subtracts instead is the **media blank**, which is the `background`
term -- optics and medium, no cells. `autofluorescence` is a per-biomass term and survives
that subtraction entirely, so the per-cell signal
`(RFU - media_blank) / (OD - od_blank)` still contains it.

This matters because the ODE inversion makes the correction exact and therefore easy to
sweep: an unremoved constant `a` in the per-cell signal contributes exactly `a * mu` to the
recovered activity, since a constant has zero derivative. Sweeping `a` as a fraction of
each construct's own 0 mM signal:

| Construct | Dose | a=0% | a=5% | a=10% | a=20% | a=30% |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| NativeYap1 | 0.1 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;autofluorescence_fraction=0.0" -->0.92 | 0.92 | 0.92 | 0.92 | 0.91 |
| AlteredYap1 | 0.1 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;autofluorescence_fraction=0.0" -->1.01 | 1.01 | 1.02 | 1.02 | 1.02 |
| AlteredYap1 | 1.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=1.0;autofluorescence_fraction=0.0" -->1.00 | 1.01 | 1.03 | 1.06 | 1.10 |
| UPRE1 | 2.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;autofluorescence_fraction=0.0" -->1.06 | 1.08 | 1.10 | 1.16 | **1.22** |
| UPRE1 | 5.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=5.0;autofluorescence_fraction=0.0" -->0.93 | 0.95 | 0.97 | 1.01 | **1.06** |
| UPRE2 | 5.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=5.0;autofluorescence_fraction=0.0" -->1.06 | 1.09 | 1.11 | 1.16 | **1.22** |

Regenerated at four biological replicates. The a=0% column is pinned cell by cell to
`outputs/autofluorescence_sensitivity.csv`, so it cannot drift away from the file without
`scripts/audit_claims.py` failing; the rest of each row is read off the same file.

The fourth plate moved every row and it did not move them all the same way. NativeYap1 at
0.1 mM went 0.90 -> 0.92 and AlteredYap1 at 1.0 mM went 1.49 -> 1.00, which is the largest
single move any fold made in this pass and is why that row is no longer quoted as an
induction anywhere. The high-DTT rows all rose: UPRE1 at 2.0 mM now sweeps 1.06 -> 1.22
where it swept 1.02 -> 1.18, so it crosses into DIFFERENT at a smaller assumed
autofluorescence than before rather than at a larger one.

**The low-dose point estimates barely move.** Across a = 0-30% they shift by at most 0.02,
and their intervals contain one at every point of the sweep: NativeYap1 at 0.1 mM stays
inside [0.74, 1.08] and AlteredYap1 at 0.1 mM inside [0.94, 1.11]. Every one of those rows
is INCONCLUSIVE at every `a`.

**The 2-5 mM DTT point estimates move.** UPRE1 at 5 mM runs 0.93 -> 1.06 and UPRE2 at 5 mM
runs 1.06 -> 1.22, and three of those rows cross from INCONCLUSIVE to DIFFERENT inside the
sweep: UPRE2 at 2 mM at a = 20%, UPRE1 at 2 mM and UPRE2 at 5 mM at a = 30%. The correction
raises every fold, so it also enlarges the rows that already resolve: UPRE2 at 1.0 mM goes
1.56 -> 1.71.

### What the sweep cannot decide

This section previously read *"The 'no induction at all' calls at 0.1-0.2 mM are robust"*
and *"The 'entirely dilution' calls at 2-5 mM DTT are not"*, as though a fold that is
insensitive to `a` were a measurement of no response. Both readings are formally withdrawn.
`data/current_claims.json` holds `biosensor.no_induction.native_low`, `.native_second`,
`.altered_low` and `biosensor.entirely_dilution.upre1` at role `historical`, each with
`"margins": null` and `"margin_source": null`, and each demanding justified
practical-equivalence margins *before* the interval is read as anything. Its own words:
**"no SESOI is invented by this registry"**. None is invented here.

An interval containing one is UNRESOLVED. It is not evidence that the construct did not
respond, a stable point estimate under an assumed background does not make it so, and a
nonsignificant fold is not equivalence. Proving absence of response needs the whole interval
strictly inside two margins fixed in advance; this panel has no such margins, so none of the
low-dose rows is decided in either direction.

<!-- audit:retracted "no induction at all" calls at 0.1-0.2 mM are robust -->
<!-- audit:retracted "entirely dilution" calls at 2-5 mM DTT are not -->

These markers are load-bearing for the same reason the two above are. Nothing in the audit
could see either sentence: the registry's scope policy states that *"unmarked prose is not
discovered by decimal, citation-year, test-count or keyword heuristics"*, so a withdrawn
inference written as ordinary prose stood for weeks beside a registry that had already
retired it. `check_retracted_claims` now fails if either phrase appears in any document
except this one and `docs/superseded/`.

The size of `a` is still not known, and one plate would fix that. It would sharpen the
2-5 mM rows. It would not resolve the low-dose rows, whose obstruction is a missing margin
rather than a missing background.

## What to do, in cost order

1. **Read the reporter-free strain in the mCitrine channel.** One plate of BY4741 across
   the same dose ladder, same gain, same protocol. This is the cheapest item here and it is
   now the *only* cheap one: it is what fixes where in the sweep the 2-5 mM DTT rows sit,
   and no amount of replication settles that, because the term it would measure is common to
   every plate. It bounds `a`; it does not supply an equivalence margin, so it cannot by
   itself resolve a row whose interval contains one.
2. **Export from the instrument file, or do not select a transform.** 20260804 cost this
   project a biological replicate for three weeks because the operator's export had a
   blank-subtraction transform selected and the raw numbers were only in the `.xpt`.
   Committing the `.xpt` recovered it; not selecting the transform would have been free.
3. **Include blanks inside every channel that is read**, not only in OD, and put them in the
   same wells every time. Three of the four plates blanked at H1-H3 and one at H4-H6, and
   two logbook entries recorded the wrong three.
4. ~~**Two more full biological replicates**~~ -- **done, and the second one was already on
   disk.** The panel is n=4 for all four constructs. Seven folds clear "no change" against
   five at n=3, one survives correcting for all twenty for the first time, and every one of
   the seven holds across the whole window sweep. A *fifth* plate is the next thing that
   would move the distribution-free floor, which is 0.125 at four and needs six to reach
   0.05.
5. ~~**Recover UPRE1's replicates 2 and 4 under a background sweep**~~ -- not needed. Both
   plates have measured backgrounds now, so a bound over an unknown one would be strictly
   weaker than what is here.

## Provenance

Read from the exports on 2026-08-26 with `plate/synergy.py`, and from
`20260804_ER_Oxidative_Replicate4.xpt` on 2026-08-30 with `plate/gen5.py`. Plate maps from
`plate/layout.py::RECORDED_PLATES`, which agrees with the logbook except on where 20260803
and 20260804 put their blanks; both entries are corrected against the instrument and carry
the correction in their notes. Fold changes from `analysis/uncertainty.py::fold_change`
against `outputs/sensor_characterisation.csv`, regenerated at n=4 by
`scripts/plate_readings.py --write --window`.
