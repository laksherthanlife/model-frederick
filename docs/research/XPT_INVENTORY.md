# All 22 Gen5 `.xpt` files, parsed — and the two the repository was wrong about

**2026-08-30.** `docs/research/DATA_SURVEY_2026_08.md` located twenty-two `.xpt` instrument
files and recorded them as *"not known — mentioned twice in prose, never located, never
parsed"*. `plate/gen5.py` was then written and proven against four of them. This is the
other eighteen.

Everything below was produced by reading the files, not by reading their names:

```bash
export YSTWIN_GEN5_XPT="…/iGEM 2026/Experiments"
python3 scripts/read_gen5_xpt.py "$YSTWIN_GEN5_XPT"/*.xpt
```

**All 22 parse.** No file raised `Gen5FormatError`, no block failed its stride check, and
every one holds all 96 wells whether or not the operator's `.xlsx` export did.

---

## 1. The inventory

`occupied` is wells whose optical density exceeds 0.12 at the last read — the medium blanks
and the empty wells sit at 0.086–0.10, so this is roughly "how much of the plate had cells in
it". *Roughly*: a handful of unused wells carry a fixed optical offset (six on each ATP plate,
one as high as 0.48) that never changes over the run, and this column counts them. Growth is
the change, not the level, and §3 uses the change wherever it matters.

| File | First read | Protocol | Channels | Reads | Span | Occupied | Repository status |
| --- | --- | --- | --- | ---: | --- | ---: | --- |
| `20260612.xpt` | 2026-06-12 | `RFP_30min-interval_24h.prt` | RFP:530,580 ×2, OD600 | 49 | 24.5 h | 36 | **not known** |
| `20260705_AFL&control_1st.xpt` | 2026-07-05 | `24h-30min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 145 | 24.1 h | 85 | **not known** |
| `20260707_oxidative_stress_preliminary.xpt` | 2026-07-08 | `endpoint_mCitrine_shake+OD600+5075100.prt` | OD600, mCitrine:480,530 ×3 | 1 | endpoint | 93 | **not known** |
| `20260708_ER_stress_1st.xpt` | 2026-07-08 | `24h-30min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 123 | 20.5 h | 42 | known — proven pair for `20260709_ER_stress_1st_(RAW)` |
| `20260712_AFL&control_2nd_AND_oxidative_stress_1st.xpt` | 2026-07-12 | `24h-30min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 131 | 21.8 h | 84 | **not known** |
| `20260713_AFL&Control3_ER Stress1.xpt` | 2026-07-13 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 125 | 20.8 h | 84 | **not known** |
| `20260728_ER_Oxidative_Replicate2.xpt` | 2026-07-28 | `4h-10min_mCitrine_OD600.prt` | OD600, mCitrine:480,530 | 25 | 4.1 h | 87 | known — proven pair, `data/plates/` |
| `20260802_AFL-debugging_1st.xpt` | 2026-08-02 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 121 | 20.1 h | 84 | known — `scripts/score_afl_circuit.py` |
| `20260803_AFL-debugging_2nd.xpt` | 2026-08-03 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 109 | 18.1 h | 84 | known — `scripts/score_afl_circuit.py` |
| `20260803_ER_Oxidative_Replicate3.xpt` | 2026-08-03 | `4h-10min_mCitrine_OD600.prt` | OD600, mCitrine:480,530 | 25 | 4.1 h | 84 | known — proven pair, `data/plates/` |
| `20260804_AFL-debugging_3rd.xpt` | 2026-08-04 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 109 | 18.1 h | 87 | known — `scripts/score_afl_circuit.py` |
| `20260804_ER_Oxidative_Replicate4.xpt` | 2026-08-04 | `4h-10min_mCitrine_OD600.prt` | OD600, mCitrine:480,530 | 25 | 4.1 h | 84 | known — the un-subtracted source of the blank-subtracted export (`DATA_INVENTORY.md`) |
| `20260805_ATP-sensor.xpt` | 2026-08-05 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 109 | 18.1 h | 54 | half known — the `.xlsx` is vendored, the `.xpt` was never named. **Now proven** (§3) |
| `20260807_BY4741_ER_Oxidative.xpt` | 2026-08-07 | `4h-10min_mCitrine_OD600.prt` | OD600, mCitrine:480,530 | 25 | 4.1 h | 31 | known as an OD plate. **The repository is wrong about it** (§2.2) |
| `20260808_ATP-sensors(ACS-ICL).xpt` | 2026-08-08 | `24h-10min_mCitrine_OD600+75100.prt` | OD600, mCitrine:480,530 ×2 | 109 | 18.1 h | 78 | half known, as above. **Now proven** (§3) |
| `24h_30min-interval_mCitrine-2.xpt` | 2026-07-01 | — | OD600, mCitrine:480,530 ×2 | 111 | 18.5 h | 90 | known — proven pair for `20260701_ER_preliminary_(RAW)` |
| `24h_30min-interval_mCitrine.xpt` | 2026-07-01 | — | OD600, **mCitrine:510,530 ×2** | 13 | 4.3 h | 90 | **not known** |
| `24h_continuous_plasmidloss.xpt` | 2026-07-02 | — | OD600, RFP:580,610 ×2 | 141 | 23.5 h | 78 | **not known** |
| `Experiment1.xpt` | no stamp | — | RFP:530,580 ×2, OD600 | 49 | none | 0 | **not known** |
| `RFP-read_20-min-interval_24h.xpt` | no stamp | — | RFP:530,580 ×2 | 73 | none | — | **not known** |
| `UPRE_Yap1_StressTest_Replicate1.xpt` | 2026-07-22 | — | OD600, mCitrine:480,530 | 25 | 4.1 h | 84 | **not known. It is biosensor replicate 1** (§2.1) |
| `endpoint_mCitrine5075-shake.xpt` | 2026-07-02 | — | OD600, mCitrine:480,530 ×3 | 1 | endpoint | 73 | **not known** |

A `.prt` protocol name is present in fifteen files and absent from seven; where it is absent
the archive simply carries no `.prt` path, and the reader reports `None` rather than guessing
from the filename. Five protocol files are named across the corpus, and the two that matter
here are `4h-10min_mCitrine_OD600.prt` — every four-hour biosensor plate *and* the BY4741
control — and `24h-10min_mCitrine_OD600+75100.prt`, the AFL and ATP-sensor runs.

### How "known" was decided

Not by name. Every channel of every file was compared value-for-value against every kinetic
block in `data/plates/`. Six files match committed text; the four already in
`tests/test_gen5_xpt.py::PROVEN` and two that were not:

| `.xpt` | committed text | fluorescence | absorbance |
| --- | --- | --- | --- |
| `20260708_ER_stress_1st` | `20260709_ER_stress_1st_(RAW)` | 5535/5535 exact | worst 5e-4 |
| `20260728_ER_Oxidative_Replicate2` | same name | 2175/2175 exact | worst 5e-4 |
| `20260803_ER_Oxidative_Replicate3` | `…Replicate3__mCitrine_2` | 2175/2175 exact | worst 5e-4 |
| `24h_30min-interval_mCitrine-2` | `20260701_ER_preliminary_(RAW)` | 10230/10230 exact | worst 5e-4 |
| **`UPRE_Yap1_StressTest_Replicate1`** | **`20260722_…_NewProtocol_ANALYSED__mCitrine_1`** | **2175/2175 exact** | **worst 5e-4** |
| **`20260807_BY4741_ER_Oxidative`** | **`20260814_BY4741_0mM_OD600_Data__OD600`** | (no committed channel) | **worst 5e-4** |

5e-4 is the arithmetic maximum for absorbance: the `.xpt` stores four decimals and the
`.xlsx` export keeps three. That is the same signature `tests/test_gen5_xpt.py` already
asserts.

---

## 2. The two identifications that change something

### 2.1 `UPRE_Yap1_StressTest_Replicate1.xpt` is biosensor replicate 1

It is the instrument file for `20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx` — the
first of the four NewProtocol replicates, and the plate the whole sensor-characterisation
pipeline starts from. Nothing in the repository connected the two, because the filename
carries neither the date convention nor the word "NewProtocol".

The evidence is not the name: all **2175** mCitrine readings match the committed text
exactly, the absorbance matches to 5e-4, and the elapsed-time axis matches to the second
(`0:09:08`, `0:19:08`, …).

**What it is worth.** A fifth entry for `tests/test_gen5_xpt.py::PROVEN`, at no cost — the
reader's evidence base goes from four plates to five, and the fifth is the panel's first
replicate rather than a preliminary. The `.xpt` also holds nine wells the export dropped
(H4–H12); all nine read OD 0.086–0.102 and 43–57 RFU, so they are empty and the committed
87 wells lose nothing.

### 2.2 The reporter-free control **was** read in the mCitrine channel

`docs/DATA_INVENTORY.md` says, of the two BY4741 workbooks:

> Neither was read in the mCitrine channel. The reporter-free control exists and the
> quantity it was made to measure was not measured. […] The size of `a` is not known, which
> is the point. It is one plate to find out.

and puts *"Read the reporter-free strain in the mCitrine channel"* at the top of its
cost-ordered list of what to do next.

**`20260807_BY4741_ER_Oxidative.xpt` holds that channel.** All 96 wells, 25 timepoints,
`mCitrine:480,530`, under `4h-10min_mCitrine_OD600.prt` — *the same protocol file as
replicates 2, 3 and 4*, so the same gain, the same interval and the same optics as the
plates the correction would be applied to. The measurement was taken on 7 August. The export
dropped the channel.

Two further corrections to the same passage, both from the file:

- The committed `20260814_BY4741_0mM_OD600 Data.xlsx` is **not a second plate**. Its six
  wells, its 25 elapsed times and its temperature series all match this 2026-08-07 run
  exactly. There is one BY4741 plate exported twice, not two BY4741 platings.
- The plate map is **A1–F3**, a descending density ladder: thirty wells across rows A–E
  reaching OD 0.24–0.70, then F1–F3 at OD 0.087–0.133 which barely grew. Columns 7–12 and
  F4–F6 are empty, at 44–60 RFU. That is 33 wells, which is exactly the count
  `DATA_INVENTORY.md` records — so the export kept the whole occupied block and dropped only
  the channel.

**What it is worth, and what it is not.** The thirty wells that grew read 164–217 RFU at
the last timepoint. Taking each plate's own floor off both sides, the reporter-free per-cell
signal is a median of **249 RFU/OD** against **4354** (replicate 2) and **3649**
(replicate 3) for the reporter constructs — so of order **6% of the median reporter well**,
2–3% of the brightest, and around 20% of the dimmest.

That is a size, not a measurement, and it is written here rather than in
`autofluorescence_sensitivity.csv` for three reasons. It splices two plates read eleven days
apart. The BY4741 plate carries no medium-only well that fluoresces — its unoccupied wells
read like air (54 RFU) while the reporter plates' medium blanks read 402 and 832 RFU — so
the arithmetic above subtracts the smaller floor, which can only make 249 an over-estimate.
And `a` in the published sweep is defined per construct against *its own* 0 mM well, which
needs the plate layout rather than a plate median.

**But the direction is the finding.** `DATA_INVENTORY.md`'s sweep runs `a` from 0% to 30%
and three of the 2–5 mM DTT readings change verdict inside that range. A first pass at the
measurement it says has never been made puts `a` near the **bottom** of the range. If that
survives being done properly, those rows stay INCONCLUSIVE rather than crossing to
DIFFERENT. That is a statement about where in the sweep they sit, and nothing more: the
"entirely dilution" reading of them is withdrawn in `data/current_claims.json`
(`biosensor.entirely_dilution.upre1`, role `historical`, `"margins": null`), so measuring
`a` cannot reinstate it. Absence of response would need equivalence margins declared in
advance, which no measurement of a background term supplies. **Doing it properly is now a
parsing job on a file that is already here, not a plate.**

---

## 3. The ATP-sensor figures: what got a writer, and what could not

`docs/REPRODUCING.md` recorded three `outputs/atp_sensor_*.png` panels that no script
produced. `scripts/analyse_atp_sensor.py --figures` is now their writer. It draws one
three-panel figure per plate block and it does **not** draw all of what was archived.

### The premise this track started from was wrong, and that is worth stating

The hypothesis was that the raw `.xpt` files might unlock the figures. They do not, because
nothing was locked. Both ATP `.xpt` files were compared against the vendored exports:

| | `20260805_ATP-sensor.xpt` | `20260808_ATP-sensors(ACS-ICL).xpt` |
| --- | --- | --- |
| fluorescence vs the export | **5559/5559 exact** | **8175/8175 exact** |
| absorbance vs the export | worst 5e-4 | worst 5e-4 |
| wells the export dropped | 45 | 21 |
| …any of them cells? | **no** — none moves more than 0.05 OD in 18 h, all at 44–69 RFU | **no** — same |
| gain-100 overflow cells | 1096 of 5559 | 127 of 8175 — **and the export marks exactly the same cells** |

So the exports already carried every number the panels need. What was missing was the
drawing code, and only that. The `.xpt` did settle one thing by exclusion: **neither plate
carries a promoterless or a constitutive strain**, which `docs/ATP_SENSOR.md` states and
which could until now only be inferred from an export that had dropped half the plate. The
45 and 21 unexported wells are empty, not controls.

### Panel by panel

| Archived panel | Drawn now? | Why |
| --- | --- | --- |
| `out_of_sync` A — RFU/OD traces, peak per condition | **yes** | plate data |
| `out_of_sync` B — peak hour vs growth-arrest hour | **yes** | plate data |
| `reasoning_chain` 1 — bar chart at 10 h | folded into the clock line on the other two panels | plate data, but a bar chart at one clock time is the reading the whole document exists to argue against |
| `reasoning_chain` 2 — same traces without the band | **yes** (the band is drawn) | plate data |
| `reasoning_chain` 3 — spread at one clock time vs at matched OD | **yes** | plate data |
| `reasoning_chain` 4 / `confidence` bottom-left — implied biomass vs assumed OD saturation | **no** | an analytic curve in a free parameter `k` that these plates do not measure. The archived panel says so on its face ("UNSETTLED — needs a dilution series") and drawing it again would not change that |
| `reasoning_chain` 5 / `confidence` top-right — 1.06× / 1.78× / 2.72× model against data | **partly** — the measured 2.72× is reproduced exactly (§3.1); the two model bars are not | the bars are an output of `generator/context.py`, and its input has moved. `docs/ATP_SENSOR.md` states `_CARBON["galactose"]["carbon"]` "has been put back to 0.45"; the table today reads **0.05**. Redrawing them would mean picking a value nobody can pin |
| `reasoning_chain` 6 / `confidence` top-middle — ICL1 >200× on non-fermentable carbon | **no** | a citation, not a measurement on these plates — and a citation `docs/SOURCE_AUDIT.md` §3.1 records as **broken**: PMID 8710508 is *"The complete nucleotide sequence of bacteriophage HP1 DNA"*, and a PubMed search for every Schüller + Entian co-authored paper found none from 1996 |
| `confidence` top-left — histogram of the hour each well reaches OD 1.0 | **no**, but its number is printed | the report already prints "conditions are 6.0 h apart at OD 1.0"; a sixteen-bar histogram of sixteen numbers adds nothing |
| `confidence` bottom-right — "two additions settle everything" text box | **no** | prose. It lives in `docs/ATP_SENSOR.md`, where it can be revised |

### 3.1 The archived numbers reproduce to the digit

This was checked before anything was drawn, because a writer that produced *different*
numbers would not be a writer for those panels — it would be a second, disagreeing analysis.

| archived value | recomputed |
| --- | --- |
| `reasoning_chain` 3, "one clock time: 3.7×" | **3.66×** at 10 h |
| `reasoning_chain` 3, matched-OD curve | **1.13, 1.40, 1.81, 2.47, 3.39** at OD 0.5→1.3 |
| `reasoning_chain` 5, "measured at OD 1.3: 2.72×" | **16888 / 6201 = 2.72×** |
| `out_of_sync` A, the three peaks | **24026 @ h12.62, 17324 @ h11.46, 12379 @ h8.96** |
| `out_of_sync` B, growth-arrest hours | **11.29, 10.46, 8.62** |

`docs/ATP_SENSOR.md`'s own peak/no-peak table is reproduced exactly by one rule, now written
down as `PEAK_DOUBLINGS`: a maximum counts only if it arrives after the culture has doubled.
Without it, eleven of the ICL plate's sixteen conditions have their maximum at t = 0, in the
settling transient.

### 3.2 What the new figure says that the archived one did not

Two things, both from drawing every condition rather than the glucose column alone.

**The "every peak sits about an hour after growth stops" claim has an exception.** The
archived panel plotted three points, all glucose, all above the identity line. There are
**nine** conditions across the three plate blocks whose maximum survives the peak rule, and
one of them — the ICL plate's mid-high sugar, gal 3:7 well — peaks at h9.12 against growth
arrest at h10.79, i.e. **1.7 h before** growth stops. The figure draws it in a different
colour and says so in the legend.

**The fold range at high density is over a subset.** At OD 1.3 only 11 of 16 conditions on
the ICL plate ever got there, and 8 of 12 on the first ACS/ICL block. The printed report has
always said so; the archived figure did not, and the new one carries the count under every
point.

---

## 4. The files the repository does not know, one line each

Ordered by what they could bear on. **None of these was analysed** — this is an inventory.

**Could add replicates to a claim that has three**

- `20260705_AFL&control_1st.xpt`, `20260712_AFL&control_2nd_AND_oxidative_stress_1st.xpt`,
  `20260713_AFL&Control3_ER Stress1.xpt` — three earlier runs of the auto feedback loop, 20–24 h
  at 96 wells, i.e. the same 20-hour geometry that makes
  [`AUTO_FEEDBACK_LOOP.md`](AUTO_FEEDBACK_LOOP.md)'s estimator accurate to ~1.3%. **They are
  not drop-in replicates.** Two of the three carry a second experiment on the same plate by
  their own names, and the readings agree. On `20260802` the four column-triplets separate
  cleanly — mean end-point mCitrine **344 / 627 / 4330 / 3097** across rows A–G, which is the
  notebook's AFL / Control / No-Gal / No-LacI map. On `20260705` they read **331 / 322 / 531 /
  528**: two blocks of six, not four of three. On `20260712` and `20260713` the whole plate
  sits ten to thirty times brighter (**7854 / 8034 / 10098 / 10781** and **6895 / 10494 /
  12590 / 14230**), which is the ER and oxidative-stress constructs their names announce.
  What does transfer is the blank convention — H1–H3 read 95–176 RFU against H4–H12's 51–63,
  medium against air, on all four plates — with one exception that is itself informative:
  `20260712` has H7–H9 at 269–274 RFU, so that plate's blank row is not laid out like the
  others'. Reading any of them needs the lab notebook page, not an assumption. *Could support* the
  AFL direction-reversal finding, which currently rests on three replicates from three
  consecutive days; *could contradict* it just as easily, which is the point of looking.

**Bears on the reporter-free / autofluorescence question**

- `20260707_oxidative_stress_preliminary.xpt` — a single endpoint read of a 93-well plate at
  three mCitrine gains (50/75/100) with shaking, 2026-07-08. The only file in the corpus with
  a three-gain fluorescence ladder on cells, so it is the one that could say how the gains
  relate; the gain-100 channel runs to 9.7e4, at the edge of the overflow sentinel.
- `endpoint_mCitrine5075-shake.xpt` — the same protocol shape on 2026-07-02, 73 occupied
  wells; the gain-50 channel reads 0–31 RFU, i.e. nothing, and gain 75 reads 51–608.

**A second reporter the model has never seen**

- `24h_continuous_plasmidloss.xpt` — 23.5 h, 141 reads, OD600 with RFP at 580/610 in two
  gains, 78 occupied wells. A plasmid-retention assay: the quantity is the *fraction* of
  cells still carrying the marker, which is a different observation model from the per-cell
  activity `observation.py` writes. *Could contradict* nothing directly; it is the only
  measurement in the project of a term — plasmid loss — that the reporter ODE assumes away.
- `20260612.xpt` — 24.5 h, 49 reads, RFP at 530/580 in two gains plus OD600, 36 occupied
  wells. The earliest dated file in the corpus and the first RFP time course.

**Instrument failures, worth recording so nobody re-derives them**

- `24h_30min-interval_mCitrine.xpt` — 4.3 h, 13 reads, 90 occupied wells, and **every one of
  its 2496 fluorescence readings is exactly −99999**, the reader's out-of-range sentinel. It
  is also the only file in the corpus reading `mCitrine:510,530` rather than `480,530`. A
  misconfigured excitation filter, caught by the instrument and recorded as a total loss; the
  OD600 channel is intact. Its sibling `24h_30min-interval_mCitrine-2.xpt`, four hours later
  the same day, is the successful re-run and is the committed `20260701_ER_preliminary` plate.
- `Experiment1.xpt` and `RFP-read_20-min-interval_24h.xpt` — well blocks present at full size
  (49 × 96 and 73 × 96), **every value exactly 0.0**, no elapsed table, no temperature series,
  no read-start timestamp. Experiment files saved from a protocol that was never run.
  `Experiment1.xpt` is 843,776 bytes against `20260612.xpt`'s 844,288 with the same three
  channels and the same 49 reads, so it is almost certainly that run's shell under Gen5's
  default name.

`−99999` is worth naming once: it is what the `.xlsx` export writes as an empty cell. On
`20260805_ATP-sensor.xpt` the sentinel appears in 1096 of 5559 gain-100 readings and the
vendored export has exactly 1096 blanks in the same channel, so the two agree cell for cell
and the `.xpt` recovers nothing extra there.

---

## 5. What was checked and found to be nothing

Said explicitly, because a survey that only reports hits is not a survey.

- **No promoterless or constitutive strain** hides in the unexported wells of either ATP
  plate. All 45 and all 21 sit at 44–69 RFU and none moves more than 0.05 OD over 18 hours;
  six on each plate carry a fixed optical offset that never changes, which is a mark on the
  plastic rather than a culture.
- **No extra data stream.** Every file holds `Contents`, `SUBSETS/1/HEADER` and
  `SUBSETS/1/DATA`, and only `24h_30min-interval_mCitrine.xpt` has a second subset — 150
  bytes, no well block.
- **No sixth committed plate.** `20260804_ER_Oxidative_Replicate4.xpt` matches nothing in
  `data/plates/` for the reason `DATA_INVENTORY.md` already gives: the committed text for
  that plate is blank-subtracted and the `.xpt` is not. That is the file whose re-export
  takes the biosensor panel from n=3 to n=4, and it remains so.
- **No `.prt` protocol file is in the corpus.** Five are named inside the archives; the files
  themselves live on the instrument PC and would carry the gain and filter definitions.

---

## 6. What this leaves undone

1. **Add `UPRE_Yap1_StressTest_Replicate1.xpt` to `tests/test_gen5_xpt.py::PROVEN`** against
   `20260722_…__OD600_1` and `__mCitrine_1`. The comparison is in §1 and it passes; the test
   file is outside this pass's paths.
2. **Derive `a` from `20260807_BY4741_ER_Oxidative.xpt` properly**, with that plate's layout
   and a stated blank convention, and correct the passage in `docs/DATA_INVENTORY.md` that
   says the channel was never read.
3. **Reconcile `docs/ATP_SENSOR.md` with `generator/context.py`** on
   `_CARBON["galactose"]["carbon"]` — the document says 0.45, the code says 0.05, and one of
   the archived figure panels was drawn from whichever it was at the time.
4. **Get the notebook pages for the three July AFL plates** before treating them as anything.
