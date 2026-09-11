# The auto feedback loop plates, and the half of the finding that a fourth plate broke

**2026-08-29. Extended 2026-08-30 with the plate layout from the lab notebook. Rewritten
2026-08-30 after a fourth plate, `20260705_AFL&control_1st.xpt`, was scored for the first
time — it changed the answer, and this document now says so before it says anything else.**

`AFL` is the team's **auto negative feedback loop** — a galactose-inducible circuit carrying
its own repressor, so raising the inducer raises the repressor which holds the output down. An
earlier automated pass read the abbreviation as "autofluorescence"; it is not, and the
glycerol-stock box settles it, listing `AFL`, `AFL control` and `No-LacI` together.

## What changed, in one paragraph

This document used to say: **naive `RFU/OD` falls from 0% to 1.8% galactose while the
dilution-corrected activity rises, in every strain, in all three biological replicates.** The
second half survives a fourth plate and the first half does not. On `20260705` the naive
readout **rises** for the same strain that falls on all three August plates, and on any single
plate the naive endpoint ratio moves — a long way — with how many hours the plate was read
for. The claim was stated as a ratio of two endpoints, and that statistic is not stable. Stated
instead as a **rank correlation with the dose**, which is what the data actually supports, the
finding holds on all four plates and is stronger than it was:

> Corrected activity ranks **positively** with galactose in **28 of 28** strain × plate ×
> read-length blocks (Spearman ρ from **+0.86 to +1.00**). Naive `RFU/OD` ranks **zero or
> negatively** in **27 of 28** (ρ from **−1.00 to +0.14**). The two readouts disagree about the
> sign of the dose response everywhere, on both experiments, at both read lengths.

The endpoint ratio was the weaker statement all along; it just took a plate from a different
month to show it.

## Two experiments, not four replicates

`20260705_AFL&control_1st.xpt` is **not** a fourth replicate of the August plates. The notebook
keeps two programmes apart and so does `scripts/score_afl_circuit.py`:

| series | plates | strains | why separate |
| --- | --- | --- | --- |
| `AFL-vs-control` | 2026-07-05 | AFL, Control | "AFL vs control (1st biological replicate)". `No-Gal` and `No-LacI` **did not exist yet** — the notebook builds them in the second half of July. |
| `AFL-debugging` | 2026-08-02, 08-03, 08-04 | AFL, Control, No-Gal, No-LacI | Run a month later, "to see what the problem might be with our AFL". |

They are scored side by side and never pooled. `outputs/afl_circuit.csv` carries an
`experiment` column so nothing downstream can pool them by accident.

Two more things the file itself corrected against the notebook and against this document's
previous version:

- **The July run samples every 10 minutes, not 30.** Its protocol is named
  `24h-30min_mCitrine_OD600+75100.prt`, and its elapsed table steps 600,000 ms throughout. The
  file wins; the name is wrong.
- **The three August runs are not "121 timepoints over 20 hours" each.** That is only
  `20260802`. `20260803` and `20260804` were stopped at **109 points / 18.12 h**. This document
  asserted otherwise for a day and was wrong.

## The plates

### 2026-08-02 / 03 / 04, from the notebook's own map, which heads its column blocks

| | |
| --- | --- |
| columns 1-3 | **AFL** — the loop, closed |
| columns 4-6 | **Control** — the open-loop control |
| columns 7-9 | **No-Gal** — the galactose UAS removed |
| columns 10-12 | **No-LacI** — the repressor deleted; the loop, opened |
| rows A-G | galactose 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8 %, with glucose made up to 2% |
| row H1-H3 | YNB medium blanks; H4-H12 empty |

Confirmed against the readings: H1-H3 read **172 RFU** against H4-H12's **59**, which is the
difference between medium and air.

### 2026-07-05 — same ladder, and the strain assignment is *inferred*, not recorded

The notebook's 2026-07-05 map gives the galactose ladder and the blanks and **stops**. Unlike
the August maps it carries **no strain header row**, so which half of the plate held which
circuit is not written down anywhere. The rows and blanks check out against the data; the
columns do not, because there is nothing to check them against.

- Blanks: H1-H3 read **123 RFU** against H4-H12's **55**. Confirmed.
- Rows: mean time to reach blank-corrected OD 1.0 lengthens down the plate, **7.05 h in row A
  against 10.05 h in row G**, the same 3-hour spread the August plates show (7.78 h → 11.37 h).
  That is the inducer slowing the cultures, and it confirms the ladder's direction. Confirmed.
- Columns: **not recorded.** The table marks these rows `strain_source =
  inferred-from-signal`, against `notebook` for the August plates.

Three things fix the assignment anyway, and `tests/test_afl_circuit.py` asserts the second and
third:

1. **Only two circuits existed.** The entry prepares "5mL master stocks for each circuit (so 2
   falcon tubes)" and aims "to compare our autonegative feedback loop with the open-loop
   circuit".
2. **The plate splits 6/6 and nowhere else.** Late `RFU/OD` runs 106-164 in columns 1-6 and
   202-387 in columns 7-12, and the right half reads above the left at *every* dose. No split
   inside either half does that.
3. **Those two ranges are the August AFL and Control ranges.** Read at matched length the July
   halves span 105-154 and 185-319, against the August AFL's 82-146 and Control's 136-341,
   while No-Gal and No-LacI span **1,187-3,698** — more than ten times anything on this plate,
   which is what "the strain is absent" looks like. The dim half is the closed loop because
   closing the loop is what makes it dim.

If the assignment were swapped the only thing that changes is which block is called the loop.
Every number below is computed per block and stands either way.

Two anomalies are reported rather than cleaned:

- **H4 reads OD 0.237, flat, for 24 hours, with air-level fluorescence (54 RFU).** The map says
  it is empty. It never grows and it is not used, but an empty well should not read 0.24 OD —
  most likely something on the plate bottom.
- **G10 is a culture that dies, and it carried a headline until 2026-08-30.** Its
  blank-corrected OD peaks at **1.36 at 10.6 h and falls to 0.65 by 24 h** while its
  neighbours climb past 1.77. It is not dropped — dropping a well after seeing its result is
  exactly the move this repository refuses.

  **This bullet used to end "and it is nearly gone at the matched 18.12 h read". That was
  false, and it was the load-bearing sentence.** Over the matched late window (13.62–18.12 h)
  G10's mean blank-corrected OD is **1.107** against its five same-dose neighbours' 1.40–1.62,
  and its late `RFU/OD` is **598.6** against their 241–281: a 2.2× outlier at the matched
  length, not a fading one. The denominator is collapsing, which is the opposite of the
  reporter-accumulation mechanism this section attributed the effect to.

  It is the single largest term in the July Control's top-dose mean, and removing it moves
  that ladder's endpoint ratio from **1.087 to 0.897** — from rising to falling, onto the same
  side of 1.0 as all three August Controls (0.855, 0.708, 0.641). Every claim that the July
  plate contradicts the August ones on the naive readout was this one well.

### 2026-07-12 and 2026-07-13 — read, and deliberately not placed

**2026-08-30.** The `AFL-vs-control` series has three instrument files and scores one. The
notebook entries for `20260712_AFL&control_2nd_AND_oxidative_stress_1st.xpt` and
`20260713_AFL&Control3_ER Stress1.xpt` carry the map of the stress panel they were co-run
with and no AFL column assignment, so the same recovery that placed 2026-07-05 was run on
their readings. **It fails on both, in two different ways, and neither plate is added.** This
section records what the readings do say, so that the attempt is not repeated blind.

Both files parse. 07-12 is **131 reads to 21.790 h**, 07-13 is **125 to 20.790 h**, both on
the same 10-minute grid: the first 109 elapsed values of each are equal to the millisecond to
the four scored plates', so a matched-length comparison would have been available had the
columns been assignable. 84 wells grew on each — rows A–G across all twelve columns.

#### The method, and the two controls that show it works

Every one of the 2,047 bipartitions of the twelve columns is scored on blank-corrected
`RFU/OD`, and a partition *separates* only if the two sides' ranges do not overlap at any of
the seven rows. Run on **2026-07-05** it returns exactly one — `{1-6} | {7-12}`, the split
this document already rests on. Run on **2026-08-02**, whose map is in the notebook, it
returns `{1-6} | {7-12}` among two or three hits, and inside each half it returns the
notebook's own blocks: `{1,2,3} | {4,5,6}` and `{7,8,9} | {10,11,12}`. Everything below was
computed four ways — final timepoint or late-window mean, as-run or truncated to 109 reads —
and the counts are given as a range where the four disagree.

| plate | over all 12 columns | inside columns 1-6 | inside columns 7-12 |
| --- | --- | --- | --- |
| 2026-07-05 (scored) | **1** — `{1-6}\|{7-12}` | 0 | 0 |
| 2026-08-02 (notebook) | 2–3, incl. `{1-6}\|{7-12}` | 1–2, incl. `{1,2,3}\|{4,5,6}` | **1** — `{7,8,9}\|{10,11,12}` |
| **2026-07-12** | **0** | **0** | 0–1 |
| **2026-07-13** | 1 — `{1,2,3}\|{4-12}` | **1** — `{1,2,3}\|{4,5,6}` | 0 |

The 08-02 control's extra hits are `{1,2,3} | {rest}` over the whole plate and `{2,3} | {1,4,5,6}`
inside the left half, both only under the late-window mean; the notebook's blocks are present
under all four. Nothing in the recovery is tuned to a statistic, and nothing below depends on
which of the four is used.

**07-12 has no layout in its readings at all.** Not one of the 2,047 partitions separates
under any of the four statistics, and none of the 31 inside columns 1-6 does either. The one
hit anywhere on that plate is `{7,8,9,10,11} | {12}`, inside the right half, at the final
timepoint as-run only — one column being brightest, gone under the other three. There is
nothing here to adopt.

**07-13 does return a unique split, and it is still not enough.** `{1,2,3} | {4,5,6}` is the
only separating partition inside columns 1-6, under all four statistics. But it is a split
*inside a block that has to be assumed first*, and — the part that settles it — the same
readings say those columns do not hold the AFL and Control strains.

#### Where the stress panel lives is recorded, and columns 1-6 carry an inducer ladder

The co-run panel's home is not a guess. `data/plates/20260709_ER_stress_1st_(RAW)__OD600.csv`
— the committed export of `20260708_ER_stress_1st.xpt`, proven value-for-value in
`docs/research/XPT_INVENTORY.md` — has exactly **45 columns: A7–G12 plus H7–H9**. In the
`.xpt`, that plate's columns 1-6 never move (largest OD change across 20 h: **0.001**, at
41–73 RFU, which is air). So an ER stress panel is six columns wide, sits in columns 7-12, and
puts its medium blanks in H7–H9.

Both unplaced plates agree:

- **07-12 has two medium blanks.** H1–H3 read **95–99 RFU** and H7–H9 read **269–274**,
  against H4–H6 and H10–H12 at **53–57**, which is air. Two media, one under each half.
- **07-13's columns 7-12 reproduce the 07-08 kill ladder.** Mean max blank-corrected OD by
  row is **1.46, 1.61, 1.51, 1.46, 0.95, 0.49, 0.45** — the cultures in rows F and G never
  reach OD 1.0 at all.
- **Columns 1-6 on both plates are delayed, not poisoned.** Mean hours to blank-corrected
  OD 1.0 lengthens from row A to row G by **+3.83 h** (07-12) and **+3.78 h** (07-13),
  against the four scored plates' **+3.28 to +4.42** — and the yield is untouched: mean max
  blank-corrected OD stays at 1.45–1.64 (07-12) and 1.61–1.69 (07-13) across all seven rows.
  Delay without loss of final OD is the galactose signature; loss of final OD is the stress
  one.

#### And the block is still not AFL and Control

The 2026-07-05 assignment rested on three things, and the third was that the two halves' `RFU/OD`
ranges *are* the August AFL and Control ranges. On these two plates that criterion does not
merely fail to help — **it excludes the assignment**. Matched to the first 109 reads,
late-window blank-corrected `RFU/OD` at row A → row G, with the corrected span beside it:

| plate | block | naive `RFU/OD` A → G | corrected span |
| --- | --- | --- | ---: |
| 2026-08-02 (notebook) | AFL, cols 1-3 | 127 → 136 | +14.5 |
| 2026-08-02 (notebook) | Control, cols 4-6 | 326 → 279 | +33.1 |
| 2026-08-02 (notebook) | No-Gal, cols 7-9 | 2822 → 2315 | +436 |
| 2026-08-02 (notebook) | No-LacI, cols 10-12 | 3263 → 1190 | +263 |
| 2026-07-05 (scored) | AFL, cols 1-6 | 149 → 154 | +14.1 |
| 2026-07-05 (scored) | Control, cols 7-12 | 293 → 319 | +25.3 |
| **2026-07-12** | cols 1-3 | **6084 → 4161** | **+480** |
| **2026-07-12** | cols 4-6 | **6170 → 4141** | **+539** |
| **2026-07-13** | cols 1-3 | **5431 → 2293** | **+441** |
| **2026-07-13** | cols 4-6 | **8973 → 3797** | **+650** |

**Well by well the two populations do not overlap at all.** Across the four scored plates —
168 AFL and Control wells — late `RFU/OD` runs **70 to 599**. Across columns 1-6 of the two
unplaced plates it runs **1,999 to 9,038**. The dimmest well in the candidate block is
**3.3×** the brightest AFL or Control well ever measured in this series, and the block means
are 19–48× the AFL and Control means at the same dose.

It is not a gain difference and not an instrument difference. `mCitrine:480,530` is the low
gain in **all six** files — the median ratio of the second mCitrine channel to the first is
8.12, 8.14, 8.15, 8.16, 8.17 and 8.18 — and the plates' own floors agree too: medium blanks
95–172 RFU and air 51–60 RFU on all six. It is not induction either: **at the first read**,
before the ladder has done anything, per-cell signal is **3,192–10,160 RFU/OD** in columns 1-6
of the two unplaced plates against **60–1,123** across every AFL and Control well of the four
scored ones. Disjoint again, before the experiment starts.

So whatever is in columns 1-6 of these two plates went into the well brighter per cell than
any AFL or Control culture in this series and then responded far harder across the ladder:
**+441 to +650** corrected span, against AFL's +11.6 to +14.5 and Control's +24.6 to +33.1.
Those numbers sit in the **No-Gal / No-LacI** band — 1,029–3,890 late `RFU/OD` and +249 to
+436 span on the August plates — and above it, not in the AFL / Control one. Naming them AFL
and Control would put a block into a series whose whole point is that its levels are
comparable across plates, off by a factor of twenty to fifty.

**What this costs and what it buys.** The `AFL-vs-control` series stays at **n=1**, below
`uncertainty.MIN_PLATES_FOR_INTERVAL = 3`, and gets no interval — the same place it was
before. What it buys is that the two files are no longer "not established" by default: they
are established as *not placeable*, one because its readings carry no column structure at all
and one because the structure they carry says the wrong strain. Recovering them needs the
notebook, and the notebook is the only thing that would.

## Read the plates at the same length or do not compare them

The four runs were stopped at **18.1, 18.1, 20.1 and 24.1 hours**. `RFU/OD` at the end of a run
is not a property of the strain: it is how much stable reporter has piled up since the culture
stopped dividing, so six extra hours answer a different question.

All four files record the same 10-minute grid from 445 s, so the **first 109 reads of each are
the same timepoints to the millisecond**. Every plate is therefore scored twice — at its own
length (`truncation = as-run`) and truncated to those 109 reads (`first-109`). Truncating the
two 109-point plates changes nothing, which is the check that this is a truncation and not a
resampling.

What the extra hours do, in the six blocks that have a longer read:

| | naive endpoint ratio 0% → 1.8% | corrected span 0% → 1.8% |
| --- | --- | --- |
| at 18.1 h | 0.36 – 1.09 | +14 to +436 |
| read on to 20.1 / 24.1 h | **0.42 – 1.42, higher in all six** | +2 to +400, **lower in all six, positive in all six** |

Reading longer pushes the naive ratio **up, toward and past 1**, in every block. It never
changes the sign of the corrected response. That is the whole argument in one table: one of
these two numbers is a rate and the other is a level that depends on when you looked.

## The numbers

All four plates truncated to the same 18.12 h, 109 reads, 3.00 h activity window, late window
from 13.62 h. `outputs/afl_circuit.csv`, `truncation = first-109`.

| experiment | strain | n plates | naive `RFU/OD` 0% → 1.8% | corrected activity 0% → 1.8% |
| --- | --- | --- | --- | --- |
| AFL-vs-control | AFL | 1 | 149 → 154 (**1.03×**) | +0.3 → **+14.4** (**+14.1**) |
| AFL-vs-control | Control | 1 | 293 → 319 (**1.09×**) | −5.9 → **+19.4** (**+25.3**) |
| AFL-debugging | AFL | 3 | 124 → 140 (**1.13×**) | −0.5 → **+13.0** (**+13.5**) |
| AFL-debugging | Control | 3 | 299 → 221 (**0.74×**) | −7.6 → **+21.1** (**+28.7**) |
| AFL-debugging | No-Gal | 3 | 2991 → 2137 (**0.71×**) | −136 → **+258** (**+394**) |
| AFL-debugging | No-LacI | 3 | 3475 → 1331 (**0.38×**) | −171 → **+88** (**+258**) |

The **corrected** column is the one that reproduces across a month, two experiments and two
comparator sets. Read at the same length, the July plate recovers **+14.1** for AFL against the
August plates' **+13.5**, and **+25.3** for Control against **+28.7** — agreement to 5% and
12% between experiments that share no reagent prep, no operator note and no week.

The **naive** column is the one that does not. For AFL it has always risen — the original
table showed 1.16× — which this document previously papered over with the sentence "every
strain shows the same reversal". That sentence was false when it was written, and the AFL
blocks alone are enough to refute it: at matched length they rise on all four plates
(1.076, 1.039, 1.300 in August, 1.030 in July).

**The Control row is a weaker piece of evidence than it looks and the difference is one
well.** This section said the naive Control ratio "says 1.09× on one plate and 0.74× on the
other three", inviting the reading that the July experiment contradicts the August ones.
It does not. The 1.09× is carried entirely by **G10**, the dying culture described above:
excluding it the July Control ratio is **0.897**, falling like all three August Controls.
The July plate adds no naive counterexample that the August plates did not already contain,
and the honest statement is that the naive readout is unreliable in a way the AFL blocks
already showed, not that a second experiment disagreed with the first.

## What is claimed, and what the fourth plate did to it

**Claimed, and now at n=4 across two experiments:** the dilution-corrected activity rises with
galactose, monotonically or near-monotonically, in every strain on every plate, and naive
`RFU/OD` does not follow it. A stable reporter accumulates as growth slows; galactose slows
these cultures; so `RFU/OD` is partly a measurement of 1/µ and the promoter's own behaviour is
underneath it. On the biosensor plates the correction *removes* an apparent response, which
always invites "the correction is just eating signal". Here it *recovers* one, four times, on
plates collected for another purpose.

**Broken by the fourth plate, and it is the more interesting outcome:** the *endpoint-ratio*
form of the claim. "Naive `RFU/OD` falls across the ladder" is not a property of these strains.
It is a property of the strain **and the hour you stopped reading**, and on 2026-07-05 it comes
out the other way. Anyone quoting a fold from a stable reporter at a single late timepoint is
quoting a number that moves when the run length moves.

## What is NOT claimed

- **Not that the circuit works.** Three of the four files are named `debugging`, and whether
  these numbers are the intended behaviour is the team's call against what they expected. What
  can be said is the arithmetic: at matched read length AFL's corrected span is **+13.5**
  against No-LacI's **+258** and No-Gal's **+394**, so the closed loop moves roughly **19-29×**
  less than the opened one across the same inducer range. That comparison is between strains
  that differ ~25× in absolute brightness, so it is a statement about these numbers and not a
  clean measure of loop gain.
- **Not that the corrected rise is GAL induction.** `No-Gal` has its galactose UAS **removed**,
  and its corrected activity still rises further with galactose than any other strain's. So
  what both readouts see is at least partly the carbon source — swapping glucose for galactose
  — and neither of them can separate that from promoter regulation. This is a limit of the
  experiment, not of the estimator, and no strain on either plate controls for it.
- **The negative activities are an artefact, not repression — and there are far more of them
  than this bullet used to admit.** Activity is `(dF/dt)/X`; a culture that has stopped
  growing with a slowly falling signal gives a negative derivative. They are reported rather
  than clipped because clipping them would hide the one place the estimator is visibly
  working at the edge of its assumptions.

  This said "at 0% galactose" until 2026-08-30, which understated it by a factor of three.
  **65 of the 196 committed rows are negative**, across four of the seven doses — 23 at 0%,
  23 at 0.3%, 17 at 0.6% and 2 at 0.9% — and 23 of the 28 strain × plate blocks contain one
  above zero dose. A third of the table is a regime where `dF/dt < 0` at `k_deg = 0`, which
  is the estimator reporting a promoter that destroys its own reporter. That is not a small
  edge.

- **The dose response is confounded with growth rate, and the confound is now in the table.**
  `mu_late` — the late-window `d(ln OD)/dt` — is a column of `outputs/afl_circuit.csv`
  because a galactose ladder is also a growth ladder here: pooled,
  Spearman(galactose, `mu_late`) is **+0.891**. So "activity rises with dose" and "activity
  rises with whatever is still growing" make nearly the same prediction, and until this
  column existed a reader could not tell them apart.

  Conditioning on it does not destroy the result, and it does not leave it intact either:

  | | Spearman(galactose, corrected) |
  | --- | ---: |
  | pooled, raw | **+0.839** |
  | pooled, partial given `mu_late` | **+0.494** |
  | within a strain × plate block, partial | median **+1.000**, 15 of 16 positive |

  **Read it as: roughly half the pooled association is shared with growth rate, and what
  survives within a block is not explained by growth rate alone.** The within-block figures
  rest on seven points each and saturate, so they are weaker than they look. Separating the
  two properly needs a strain whose promoter does not respond to galactose, grown on the same
  ladder — which is the control the previous bullet already says nobody has built.
- **No interval, and n=1 for the July experiment.** The `AFL-vs-control` series has three
  `.xpt` files (07-05, 07-12, 07-13) and only the first is scored here. The other two were
  parsed and put through the same layout recovery that placed 07-05, and it fails on both —
  07-12 has no separating column partition anywhere in its readings, and 07-13's one clean
  split sits in a block whose every well is brighter than every AFL and Control well in the
  series, by at least 3.3× and on the means 19–48×. See *"2026-07-12 and
  2026-07-13 — read, and deliberately not placed"* above. One plate is below
  `uncertainty.MIN_PLATES_FOR_INTERVAL = 3`, so the July series gets no interval at all, and
  the August series at n=3 sits exactly on that floor where the exact sign test bottoms out at
  **p = 0.25** — no single dose can be significant, by design, exactly as `docs/FINDINGS.md`
  records for the biosensor folds. What carries weight is the monotone rank correlation across
  seven doses, repeated on four plates from two experiments.
- **Gain 100 is unused.** Both gains were recorded on every run; only 75 is read here, matching
  the sensor pipeline. Two gains on one culture is a linearity check nobody has used.

## Why this dataset matters to the model, before any biology

Every biosensor plate this project analyses is **four hours** long.
`analysis/estimator_accuracy.py` established that the activity inversion's error is governed by
the smoothing window's **duration in hours** — so on a four-hour run the rule collapses onto its
floor and sits at ~8% single-well error.

The 2026-07-05 plate is a closer match to that table than any plate in the project: **145 points
over 24 h** is, exactly, the geometry `outputs/estimator_accuracy.csv` labels *"the regime the
table was measured in"*, where the auto window is 4.0 h and the measured median single-well
error is **1.31%** (row `n_points=145, duration_h=24, window_h=4.0, is_auto_window=True`). The
August plates, at 3.00-3.33 h windows, fall between that table's 2.0 h row (2.30%) and its 4.0 h
row (1.31%).

**These are the only plates in the project measured in the regime the estimator is good at**,
and they are now four rather than three.

## Files

`scripts/score_afl_circuit.py` → `outputs/afl_circuit.csv` (196 rows: 4 plates × 2 read lengths
× strains × 7 doses). Tests in `tests/test_afl_circuit.py`. Raw `.xpt` files are the
instrument's own and are not redistributable; set `YSTWIN_GEN5_XPT` to read them.
