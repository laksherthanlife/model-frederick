# Six measurements for a wet lab, and the ones this team can run

> ## THE RANKING BELOW ASSUMES EQUIPMENT THIS TEAM DOES NOT HAVE
>
> **Do not build a work plan on §"The ranking".** The six-item ordering in this file was
> written for a laboratory with an **HPLC**, an **LC-MS/MS** and **strains carrying a
> heterologous pathway** — carotenoid, PHB or a secreted protein. This team has a **plate
> reader**, **qPCR** and **four mCitrine stress biosensors** (UPRE1, UPRE2, NativeYap1,
> AlteredYap1) and none of those three things. **All six are out of reach as experiments.**
> Three fragments are not, and none is bench work: M4's **step 1** is a citation check, M8's
> **step 1** is a letter to another laboratory, and the "0.46 to 0.47" discrepancy is a
> documentation repair. All three can be done today, at a desk.
>
> Nothing is deleted. A laboratory with the equipment may run any of it, and this repository
> does not delete recorded work. The plan for *this* team is
> [§What this team can measure](#what-this-team-can-measure), which is a different list
> answering different questions.

Six measurements block this project. They were scattered in prose across
[`DISTANCE_TO_THE_VISION.md`](DISTANCE_TO_THE_VISION.md) §6,
[`WHAT_IS_LEFT.md`](WHAT_IS_LEFT.md), [`HARD_TESTS.md`](HARD_TESTS.md) and
[`EXTERNAL_PRODUCT_VALIDATION.md`](EXTERNAL_PRODUCT_VALIDATION.md). This file is the single
place they live, written for someone with a bioreactor, an HPLC and an LC-MS/MS rather than
for someone reading the source.

*(The heading said "Five" until 2026-09-09 while the table listed six; M9 was added
2026-09-01 and the title was never updated. It says six because six are specified.)*

> **A CHEMOSTAT IS NOT REQUIRED FOR ANY OF THEM, AND NEITHER IS A FERMENTER.** Every
> measurement below that says "chemostat run" means "a run at a *known* growth rate".
> [`PROTOCOLS.md` §P3](PROTOCOLS.md) replaced the chemostat with an exponential fed-batch on
> 2026-09-02 — one vessel, no effluent line, no weir. An audit on 2026-09-03 then found that
> §P3 itself needs **61.2 mmol O₂/L/h**, and therefore a sparged, stirred, DO-controlled
> vessel, which is not obviously easier to come by than a chemostat.
> [**`PROTOCOLS.md` §P4**](PROTOCOLS.md) supersedes both. The balance the model closes,
> `d[X]/dt = v_in − v_out(X) − μ[X]`, is *per gram of biomass* and **the vessel appears
> nowhere in it**, so a shake flask converges on exactly the steady state `pathway/solve.py`
> returns; the cost of having no chemostat is a **generation count**, not a different answer.
> A plate reader taking 145 reads over 24 h measures μ to **SE 0.00021 /h (0.21%)** against
> the fed-batch protocol's **0.00074 (0.73%)** — **3.6× more precisely**. `μ` must be known
> and roughly constant; it does not have to be held.
>
> **So "no chemostat" is not what blocks M5, M6 or M7 at this bench. The HPLC and the strain
> are.** The cultivation half of M5 and M7 is already reachable with a flask and the plate
> reader. M6 is the exception even there: a flask gives whatever μ the medium gives, and all
> six of §P4's in-range carbon × temperature combinations land at or above 0.101 /h, so
> **D = 0.05 has no flask route on the record**. The analytical half and the genotype are not
> reachable for any of them.

**The numbering M4–M8 is deliberate**, and matches items 4–8 of `DISTANCE_TO_THE_VISION.md`
§6 so that cross-references from the other documents still land. **M9 is outside that
mapping** — it was added 2026-09-01 for the saturating export outlet and does not correspond
to that document's item 9, which is the concentration-based flux law. The *run order* is
different from the numbering and is the table below.

Every specification states what would **confirm** and what would **refute**, because a
measurement that cannot come out the other way settles nothing. Where a number here was
computed for this document rather than read out of another one, the last section says how, and every such
computation was run rather than estimated.

---

## What each of the six needs, and whether this team can run it

Read this before the ranking. Cultivation is largely settled by §P4 above — a flask and the
plate reader. What is not settled is the **analytical method** and the **strain**.

| item | analytical method it needs | strain it needs | this team? |
| --- | --- | --- | :---: |
| **M8** phytoene | HPLC with a **285 nm** channel, phytoene authentic standard | access to Elizondo's stored chromatograms, or a β-carotene strain to re-grow | **NO** |
| **M5** β-carotene ceiling | HPLC, β-carotene authentic standard, gravimetric DCW, residual glucose + ethanol | a **high-producing CEN.PK carotenoid strain**, ideally extra *crtYB* | **NO** |
| **M4** acetyl-CoA | **LC-MS/MS** with authentic standards for four acyl-CoAs, one quench protocol | **SCKK006** (PHB) and its parent | **NO** — except step 1 |
| **M6** third growth rate | HPLC as M5, plus a route to hold μ = 0.05 /h | the three Elizondo β-car strains | **NO** |
| **M7** batch calibration state | HPLC (≈36 injections), gravimetric DCW, residual glucose + ethanol, RT-qPCR for CrtE | the Elizondo β-car strains | **NO** |
| **M9** export capacity | HPLC, or quantitative MS, or a calibrated immunoassay; paired pellet + supernatant | a **secreted-product** strain (TAL, farnesene, resveratrol, a recombinant protein) at two expression levels | **NO** |

**The premise was verified, not assumed.** [`DATA_INVENTORY.md`](DATA_INVENTORY.md) records
every plate this project holds as a Synergy H1 **OD600 + mCitrine** read and nothing else;
`qpcr.py::STRESSOR_FOR_CONSTRUCT` declares the four constructs and the two agents they were
challenged with (**DTT** for the UPRE pair, **H₂O₂** for the Yap1 pair);
`qpcr.py::TARGET_FOR_CONSTRUCT` declares the two endogenous transcripts (Hac1, TRX2). There
is no HPLC, no LC-MS/MS and no chromatogram anywhere in the data inventory, and no
carotenoid, PHB or secreted-product strain anywhere in it either.
[`PROTOCOLS.md` §P3](PROTOCOLS.md) says so in its own words, declining to specify the
carotenoid inner-filter control because it is *"not applicable to the current biosensor
strains, **which carry no carotenoid pathway**"*.

### M7 is not the no-chemostat option, and reads as though it were

M7 is titled "A batch calibration state", its effort line reads *"The cheapest cultivation
here: three flasks, one day, no chemostat"*, and it is the only item whose run steps never
mention a vessel. A reader planning around equipment will land on it. **It is still out of
reach, and by a wider margin than the items above it**, because its cost was never the
vessel:

- its step 3 measures **β-carotene, lycopene and phytoene** at four time points in three
  flasks — the effort line itself puts that at *"roughly 36 HPLC injections"*;
- it measures **gravimetric dry weight** on at least two points;
- it measures **residual glucose and ethanol**, so the fermentative state is recorded rather
  than assumed;
- it measures **relative CrtE mRNA**, which needs a strain carrying *crtE*;
- and every one of those is measured **on a carotenoid-producing strain**, which does not
  exist in this collection.

Of that list this team has the OD600 trace, μ from the slope of ln(OD) — measured better
than a fermenter would give it — and RT-qPCR. It has no assay for any of the four analytes
and no strain to run them on. After §P4 removed the chemostat for *every* item, M7's one
advantage over M5 and M6 disappeared: it is now simply M5's question asked in a mode the
calibration never covered.

### The one piece of the six this team can do today

**M4 step 1 is a desk task, and it is the stated first action of the item this table ranks
highest on value.**
The two acetyl-CoA numbers the whole M4 argument rests on — 0.0199 and 0.85 μmol/gDW, a
42.7-fold gap — are named in `scripts/thiolase_threshold.py` only as "Kolbeinsen & Bruheim
2021, supplementary Table S1" and "Kozak & van Rossum 2016", and **neither carries a PMID or
a DOI anywhere in this repository**. Resolving them needs a literature search and no
equipment at all. Its own specification says what to record: identifier, strain, medium,
dilution rate, quench and extraction protocol, and whether the value is whole-cell or
compartment-resolved. **If a source is paywalled, or its abstract does not state the quench
or the compartment, record it UNRESOLVED with that reason.** Do not infer a protocol from an
abstract — the audit note already sitting in
`data/carotenoid_batch/published_batch_titres.tsv` against Yamano 1994 (PMID 7765036),
*"paywalled; abstract names beta-carotene AND lycopene but describes no assay. Earlier
beta_carotene call unfounded"*, is the mistake this instruction exists to prevent.

Two smaller desk fragments sit in the same category and are recorded here so they are not
mistaken for bench work: **M8's step 1** asks for a re-integration of chromatograms held by
another laboratory, so the reachable act is a letter to the Elizondo/Saa group asking
whether `SysBioengLab/BcarGRASP`'s underlying runs stored the full spectrum — an answer of
"no" closes M8 for everyone and is worth recording; and the **"0.46 to 0.47" discrepancy**
in the last section of this file is a documentation repair needing no instrument.

---

## The ranking

**This ordering is for a laboratory with an HPLC, an LC-MS/MS and a pathway strain. It is
not this team's plan.** See the table above for why, and
[§What this team can measure](#what-this-team-can-measure) for what is.

Two orderings, and they disagree, so both are given.

| run order | item | what it settles | effort | value | value per unit effort |
| :---: | --- | --- | --- | :---: | :---: |
| **1** | **M8 phytoene** | whether "pathway flux" means pathway flux | one standard, one re-integration of chromatograms that already exist | medium | **highest** |
| **2** | **M5 β-carotene ceiling** | the model's sharpest falsifiable claim, asserted and never tested | one run at a known growth rate, if a high producer already exists | high | high |
| **3** | **M4 cytosolic acetyl-CoA** | originally the thiolase mechanism; that ordering is now refuted, while the cross-protocol pool splice is unresolved | two culture trains, one quench protocol, LC-MS/MS on four acyl-CoAs | **highest (historical ranking)** | medium |
| **4** | **M6 third dilution rate, at D = 0.05** | the growth exponent, and two refusals | three steady states at 100 h each | medium | medium |
| **5** | **M7 batch calibration state** | the most common question a user asks, which the model refuses | flask time course, ≈36 HPLC injections | medium | medium |
| **6** | **M9 export capacity of a secreted product** | the two constants `solve_pathway` refuses on, and whether a saturating outlet is the right shape | two steady states, paired pellet and supernatant assay | medium | medium |

**The ranking is historical, not a newly validated power or value calculation.** M4 was
ranked first on value and third to run because M8 uses existing samples and M5 uses one
vessel before an acyl-CoA method is ready. The later M5 literature result refuted the old
ceiling, and the energy-unit correction refuted M4's proposed threshold ordering. Those
failures remain below; neither the old ranking nor the withdrawn flux-null bar establishes
that M4 alone can turn the environment axis into a validated model. Same-laboratory pool
measurements still address the cross-protocol splice, which is a narrower question.

**M8 is first to run and is not the most valuable.** It is first because the sample already
exists: Chen 2016's result means the current declaration is known-wrong rather than
unverified, and one re-integration converts a documented error into a measurement.

---

## What this team can measure

A plate reader, qPCR, and UPRE1 / UPRE2 / NativeYap1 / AlteredYap1. Every item below runs on
exactly that, on strains already in the collection, with one exception that is labelled.

**Two things this list is deliberately not.** It is **not a promise to improve the
dose-response panel.** That panel is at n=4 across all four constructs, seven of twenty
estimable folds clear "no change" unadjusted, and none survives a family-wise correction —
and the reason is arithmetic rather than effort: a cluster bootstrap on four plates has a
smallest atom of `(1/4)⁴`, and the exact sign-test floor is `2·(1/2)⁴ = 0.125`, so **five
plates is the first count at which a family-wise question can be asked at all, six the
first at which a single fold could reach 0.05, and ten the first at which one could survive
correcting for twenty** (`FINDINGS.md`). More plates move that floor and nothing else — and
a fifth plate buys the *question*, not an answer, so do not plan around one.

The low-dose rows are worse than underpowered: their obstruction is a **missing
equivalence margin**, which no plate supplies — `data/current_claims.json` holds them at
role `historical` with `"margins": null`, and *"no SESOI is invented by this registry"*.
None is invented here either.

And it is **not a list of new claims.** Every item is a measurement the repository already
identifies as missing, sitting in `PROTOCOLS.md` or in `docs/research/`, that was never
ranked into a plan. Ranked by value per unit effort:

| run order | item | what it settles | equipment | strain | value per unit effort |
| :---: | --- | --- | --- | --- | :---: |
| **1** | **B1 OD600 linear range + `gdcw_per_od`** | how much of every plate this project has ever read is usable | plate reader, filter, drying oven, balance | any | **highest** |
| **2** | **B2 the ER anchor, re-run on `KAR2`** | whether the UPRE reporters report the UPR — the qPCR arm currently answers nothing | qPCR + plate reader, six primer pairs already specified | UPRE1, UPRE2, BY4741 | high |
| **3** | **B3 the oxidative anchor at 1.0 mM** | the same question for the Yap1 pair, at a dose nobody sampled | rides on B2's plate; existing TRX2 pair | NativeYap1, AlteredYap1 | high |
| **4** | **B4 `k_deg` for the actual constructs** | whether the panel measures induction or `1/μ` | plate reader + cycloheximide | all four | high |
| **5** | **B5 a growth-rate-only ladder** | the control that separates a response from a dilution artefact | plate reader | all four | medium |
| **6** | **B6 crosstalk, re-read kinetically at 30 °C** | puts the specificity claim on the same readout as the dose-response claim | plate reader | all four | medium |
| **7** | **B7 a STRE-general reporter** | breaks a designed-in degeneracy that removes a whole dimension from the generator | plate reader | **one construct that does not exist yet** | medium |

### Free, whatever else is run

None of these is an experiment; each is a setting on the next plate, and each has already
cost this project something.

1. **Read more often, and run longer — in that order.** The activity estimator's smoothing
   window is chosen in **hours**, not in points, and a four-hour run caps it at
   `duration / 3`. Read straight off `outputs/estimator_accuracy.csv`, at the window the rule
   actually picks (`is_auto_window=True`):

   | geometry | auto window | median single-well error |
   | --- | ---: | ---: |
   | the committed plates — 25 points / 4.00 h | 1.333 h | **3.25%** |
   | 145 points / 4.00 h | 1.333 h | 1.55% |
   | 145 points / 24.0 h | 4.0 h | **1.31%** |

   **The larger lever is the read interval, not the run length**, and this document had it
   the wrong way round: at fixed duration, 25 → 145 points buys **2.1×**; at fixed density,
   4 h → 24 h buys a further **1.18×**. The catch is that 145 points in 4 h is a read every
   1.7 minutes across both channels, which the schedule may not permit — so the honest advice
   is to shorten the interval as far as the protocol allows *and* to run past four hours,
   for **2.5× in total**. It does **not** rescue any fold: the folds are limited by plate
   count and by the missing margins above.

   **The "~8% single-well error on a four-hour run" figure in
   [`docs/research/AUTO_FEEDBACK_LOOP.md`](research/AUTO_FEEDBACK_LOOP.md) is the
   *superseded* one-sixth rule and must not be quoted against the current one.**
   `reporter.py::TARGET_ACTIVITY_WINDOW_H` says so in its own words — "one sixth is 0.667 h,
   which collapses onto the floor and lands at 8% where 1.33 h gives 3.2%" — and the shipped
   rule returns 1.333 h on that geometry, verified by calling
   `default_activity_window_h(np.linspace(0, 4, 25))`. This document quoted the 8% until
   2026-09-09, which made the available gain look like 6× where it is 2.5×.
2. **Put blanks inside every channel that is read**, in the same wells every time. Three of
   four plates blanked at H1–H3 and one at H4–H6, and two logbook entries recorded the wrong
   three.
3. **Do not select an export transform.** Plate `20260804` cost this project a biological
   replicate for three weeks because its `.xlsx` was exported blank-subtracted; only the
   `.xpt` still carried the background. Committing the `.xpt` recovered it. Not selecting the
   transform would have been free.
4. **Keep reporter-free BY4741 wells on every plate** (`PROTOCOLS.md` §P3). This is what
   makes `ReporterOptics.autofluorescence` a measurement rather than a guess, and it has
   already paid: `a` is **+18.2 RFU/OD, 95% CI [−6.1, +42.4]**, which includes zero and is
   **0.42% of a reporter well (upper bound 0.97%)** against a sweep whose lowest
   verdict-changing fraction is 5% — so every call in the panel survives it. The remaining
   refinement, deriving `a` per construct against its own 0 mM well, is a **data-reduction
   task, not a plate to run**. Do not schedule it as bench work.

### B1 — OD600 linear range, and `gdcw_per_od`, on this reader and this plate

*(`PROTOCOLS.md` §P1, never ranked; run order 1)*

**What is blocked.** `OpticalQualityGate.od_linear_max` carries a **guess**. It decides how
much of every plate is usable, and at one inoculation density it **rejected 81% of one
experiment**. In a 96-well plate the optical path is ~0.3–0.5 cm, so linearity ends well
below the cuvette figure of 1.0 — but *where* is a property of this reader, this plate type
and this fill volume, and a literature value does not transfer.

**Run.** §P1 as written: a 12-step two-fold dilution of a culture at OD ≈ 2.5–3 in the same
medium, loaded at the **exact fill volume, plate type and lid used in real runs**, three
technical replicates, ≥3 media-only blanks. `calib.od.fit_od_calibration` consumes it and
`calib.optical_gate(tolerance=0.05)` drops straight into G1. **Filter and dry a known volume
of the stock in the same session** for `gdcw_per_od`.

**Do not also run this on the fluorescence channel — that half is already answered, on data
in hand.** Seventeen of the 22 Gen5 `.xpt` files read mCitrine twice per step at gains 75 and
100, and `scripts/score_gain_linearity.py` used the pair as a dilution series in the light
over 96,459 paired readings: the multiple is **8.15** (per-plate 8.079–8.217), the median of
measured/predicted stays within 0.999–1.001 from 2,000 RFU to 99,998, and there is no
roll-off — the channel fails by writing an overflow sentinel, not by compressing. (This is
why `AUTO_FEEDBACK_LOOP.md`'s "a linearity check nobody has used" no longer holds.) It says
nothing about `od_linear_max`, which is an absorbance property and the only reason B1 exists.

**What confirms, what refutes.** The ceiling is exactly `tolerance / k` — verified by running
`ODCalibration.linear_range_max`, whose two-step form algebraically reduces to it. A measured `k`
materially larger than the one behind the current guess means plates already analysed
contain wells above the linear range, and the affected readings must be re-gated rather than
re-interpreted. A measured `k` smaller means the gate has been discarding usable wells.

**What it cannot settle.** Nothing about the biology. It is an instrument property, and it
belongs in a methods section — which is also why it is first: it is the cheapest item here
and every other item on this list is read through it.

### B2 — The ER anchor, re-run on `KAR2` with a real reference set

*(`docs/research/UPR_ANCHOR.md`, never ranked; run order 2)*

**What is blocked, and it is the whole ER half of the qPCR dataset.** The existing Hac1
assay is **91.4% genomic DNA** — median RT− margin **+0.13 cycles** on UPRE1 — and the only
published gDNA correction is validated to 60% gDNA, a margin of about 0.74 cycles, so the
data sits far outside the correctable regime; **7 of 36 Hac1 readings have a *negative*
margin**, which falsifies the additive model the correction rests on. Underneath that, total
`HAC1` is the wrong transcript: Ire1 converts the unspliced message *directly* into the
spliced one, so the sum barely moves (Kawahara 1997, Northern blot: *"the sum of the two
mRNA species remained almost constant"*; Leber 2004 confirms it under DTT). Two independent
fatal faults. `G4_ANCHOR.md`'s verdict is *"unsalvageable and must be re-run"*, and this file
had never said so.

**What to measure.** `KAR2` ΔΔCq against the geometric mean of **`TAF10`, `ALG9` and
`TFC1`**, with **`HAC1i` exon-junction** qPCR on the same cDNA as a gDNA-proof witness that
Ire1 fired. **Nothing needs designing, and one thing still needs measuring.**
`UPR_ANCHOR.md` §7.2 tabulates forward and reverse sequences and amplicon sizes for all six
pairs — `KAR2`, `HAC1i`, `HAC1` total, `TAF10`, `ALG9`, `TFC1` — each verified to occur
exactly once across the 6039 S288C ORFs, and the two junction primers to occur zero times,
which is what makes `HAC1i` gDNA-proof. **Published efficiencies exist only for the three
reference genes** (`TAF10` 96%, `ALG9` 93%, `TFC1` 91%, Teste 2009). `KAR2`'s and `HAC1i`'s
must be measured here, which `UPR_ANCHOR.md` §7 requires of every assay anyway — a 5-point,
4-fold dilution of pooled cDNA, redesigning any pair landing outside E ≈ 0.9–1.05.

**The design changes are worth more than the reagents.** Five doses to **2 mM**, sampled at
**45 min**, is worth about **20× the replicates**. Both halves matter and both were wrong
before: the old anchor was sampled at 0.2 and 0.5 mM, and in Pincus 2010's titration DTT at
1.5 mM and below produces **less than 10% of the maximal UPR response**; and the UPR at low
DTT is a **pulse**, deactivating within 4 h at 2.2 mM and 2 h at 1.5 mM, while mCitrine is
stable and the plate reader records the time-integral. **A flat anchor beside a risen
reporter is the expected result of a correctly working system sampled hours late.** Use
sacrifice wells — one per dose per construct, harvested at 45 min while the rest of the
plate keeps reading — so the pairing stays at the level of the individual well.

**What confirms, what refutes.**

- **The reporters report the UPR:** `KAR2` rising monotonically across five doses with
  `HAC1i` splicing up on the same cDNA.
- **Translational shutdown, not a dead reporter:** `HAC1i` up and `KAR2` flat. Geronimo 2025
  showed cycloheximide 10 min after DTT abolishes `KAR2`, `ERO1`, `SIL1` and `JEM1`
  induction despite successful splicing. Running `HAC1i` alongside `KAR2` detects this
  directly, and it is why both are run.
- **The consequential negative:** `KAR2` flat *and* `HAC1i` flat at doses where the
  fluorescence reporter rose. That would say the reporter's rise is not the UPR, which is
  the largest claim in the biosensor panel and has never been tested against an endogenous
  transcript that works.
- **Settles nothing:** any of the above at a dose below 1.5 mM, or with RNA taken after
  45–90 min, for the two reasons above.

**What it cannot settle.** Whether the model's stress state constrains anything downstream.
It settles what the reporter is a reporter *of*, and no more.

### B3 — The oxidative anchor, at a dose the sensors actually respond to

*(`docs/research/G4_ANCHOR.md` §3; run order 3, on B2's plate)*

**What is blocked, and this is a much better problem than B2's.** The TRX2 assay is
**clean**: an 8.60-cycle RT− margin, **0.3% gDNA**, passing QC on 100% of readings. It is
merely underpowered — and **it was measured at the wrong doses**. Anchor RNA was taken at 0,
0.2 and 0.5 mM H₂O₂ while the dose-response analysis puts the responsive window at roughly
**0.5–1.0 mM**, so two of three anchor points sit at or below the bottom of the range.

**Why this is not "run more plates".** `gates/g4_anchor.py` derives the requirement as
`n_reps × (target_snr / snr)²`, and **the quadratic is the whole story**: at SNR 0.39
NativeYap1 needs roughly **80** replicates and at SNR 1.19 AlteredYap1 needs about **9**.
Doubling the effect size is worth four times as many replicates. Sampling the anchor at
**1.0 mM** may recover more SNR than any number of extra plates, and it costs three wells.

**What confirms, what refutes.** TRX2 ΔΔCq at 1.0 mM rising clear of the between-replicate
scatter confirms that the Yap1 reporters track their regulon at the dose where the
fluorescence panel produced its one oxidative call — **AlteredYap1 at 0.5 mM H₂O₂, fold
1.37, [1.06, 1.70]**, the first per-dose call either oxidative sensor has produced. A flat
TRX2 at 1.0 mM, with the assay's 0.3% gDNA and a clean melt, is a real negative and belongs
in `NULL_RESULTS.md`.

**What it cannot settle.** NativeYap1 at its current effect size. At SNR 0.39 no dose this
team can safely run is going to reach 80 replicates' worth of power, and saying so is the
honest output.

### B4 — `k_deg` for the actual constructs

*(`PROTOCOLS.md` §P2, never ranked; run order 4)*

**What is blocked.** `k_deg` *"separates a reporter that reports from one that reports
`1/μ`"*, and it is currently **assumed zero**. That assumption is the textbook position for
a stable YFP and it holds on the newer plates; on the older, badly-conditioned ones it did
not. Every fold in the dose-response panel is conditioned on it, and the panel's entire
mechanism — *a stable reporter accumulates as growth slows* — is a statement about this
parameter.

**Run.** §P2 as written: mid-exponential culture split three ways — cycloheximide at
100 µg/mL, vehicle control, reporter-free parent — read every 10 min for 6 h.
`calib.kdeg.fit_kdeg_from_chase` subtracts residual growth and photobleaching and **reports a
negative `k_deg` rather than clamping it**, which means translation was not actually blocked.

**The part that is easy to get wrong.** Degradation runs with the clock; photobleaching runs
with the number of reads. **They separate only across read intervals**, so this needs two
identical plates — one read every 10 min, one every 30 — and `calib.kdeg.partition_loss`.
One read interval gives only the total, **however many plates are run**.

**A free lower bound exists now.** `minimum_consistent_kdeg(t, RawOD(od), RawRFU(rfu))`
returns the smallest loss rate under which inferred promoter activity never goes negative.
Use it as the prior until the chase is done; it is not a measurement.

**What confirms, what refutes.** `k_deg` indistinguishable from zero confirms the stable-
reporter reading and leaves the panel's interpretation standing. A materially positive
`k_deg` refutes it, and the correct response is to re-derive the panel with the measured
rate rather than to re-interpret the existing folds.

### B5 — A growth-rate-only ladder

*(new here; run order 5)*

**The question nobody has asked.** Every point estimate in the panel rests on the mechanism
*"a stable reporter accumulates as growth slows"* — which means a treatment that slows growth
**without** challenging the ER or the Yap1 regulon should move the reporter too. **No such
control has ever been run.** The one place this project has looked, it found exactly the
ambiguity: on the auto-feedback-loop plates, galactose is both the inducer and a growth
inhibitor, pooled Spearman(galactose, `mu_late`) is **+0.891**, and conditioning the dose
association on growth rate takes it from **+0.839 to +0.494** — roughly half the pooled
signal shared with growth rate. Within a single strain on a single plate it survives almost
intact (median partial ρ +1.000, 15 of 16 blocks positive), so growth rate does not account
for the ranking; but nothing separates the two, and `outputs/afl_circuit.csv` carries
`mu_late` so this is checkable rather than assumed.

**Run.** The four constructs on a ladder that spans the same μ range the DTT and H₂O₂
ladders produce, driven by **carbon source and temperature rather than by a proteotoxic
agent**. §P4's own table gives the dials and the range: of twelve carbon × temperature
combinations, six land inside [0.101, 0.2543] /h — glucose at 20 and 15 °C, galactose at 25
and 20 °C, ethanol at 30 and 25 °C. Read the standard OD600 + mCitrine kinetic; μ comes from
the ln(OD) slope on the same wells.

**What confirms, what refutes.** Reporter activity flat across a μ ladder that matches the
stressor ladders' μ range says the panel's folds are a response and not a dilution artefact
— the single most valuable thing this equipment can say about its own data. Reporter
activity tracking μ on the growth-only ladder says the confound is real and quantifies it,
and the panel's folds must then be reported conditional on μ.

**What it cannot settle, and this caveat is load-bearing.** **No growth-slowing treatment is
guaranteed inert to the stress promoters.** Msn2/4 activity is coupled to reduced growth rate
by definition of the ESR (Gasch 2000; Brauer 2008 show ESR expression tracking growth rate
across nutrient limitations), which is precisely why `UPR_ANCHOR.md` **rejects** CTT1 as a
qPCR anchor. UPRE and Yap1 elements are not STRE, and this project's own crosstalk panel
measures the off-target response rather than assuming it away — but a null here **bounds**
the confound, it does not eliminate it. Carbon source is also not a clean μ dial: it is the
one environment channel this project has evidence for
(`docs/research/` and `outputs/env_to_product.csv`), so a positive result on a carbon ladder
is ambiguous between μ and carbon in a way a temperature ladder is not. Run temperature as
the primary dial for that reason, and carbon as the second.

### B6 — The crosstalk panel, re-read kinetically at 30 °C

*(`data/crosstalk/SOURCE.md`; run order 6)*

**What is blocked.** Sensor specificity — do the ER sensors respond to H₂O₂, do the
oxidative sensors respond to DTT — currently rests on plates read as
**`Reading Type: Reader`, endpoint, not kinetic**, at a recorded **24.9 °C rather than 30**.
No kinetic trace means **no growth rate**, which means **the dilution correction cannot be
applied to the crosstalk data at all**. So the project's specificity claim rests on naive
`RFU/OD` — the one readout it has independently shown can invert a result: on the AFL plates
naive `RFU/OD` fell 3390 → 1354 (**0.40×**) across the same ladder on which the corrected
activity **rose**, −151 → +85.

**Run.** The same 28 conditions, the same three biological replicates, the same gain — as a
**kinetic** read at **30 °C** on the standard schedule. Nothing else changes, and the
existing endpoint plates stay on the record as what they are.

**What confirms, what refutes.** Specificity surviving the correction confirms the
on-target/off-target split at the same standard as the dose-response panel. Specificity
reversing under the correction refutes it — and on the evidence above that is not a remote
possibility.

**What it cannot settle.** It re-measures a claim rather than making a new one, which is why
it ranks below B5 despite costing less.

### B7 — A STRE-general reporter — the one item needing a construct

*(`generator/stress_panel.py`; run order 7. **Requires cloning one new reporter**, not new
equipment.)*

**What is blocked, and it is designed in rather than unlucky.** Three asserted constants in
`stress_panel.py` jointly cancel the agent out of the ESR arm. Verified by reading the code:
`_LETHAL_MULTIPLE = 3.0` sets `lethal_dose = 3.0 × ec50`; `_GENERAL_POTENCY = 2.5` sets
`ESR EC50 = 2.5 × ec50`; and `healthy_ladder` tops out at `ceiling = 0.55` of the lethal
dose. So

```
ladder_top / ESR_EC50 = 0.55 * 3.0 / 2.5 = 0.66      <- the agent's own EC50 cancels
```

and the ESR saturation at the top of every ladder is **identical for 23 of the 25
stressors**. Only **DTT and H₂O₂** differ, because their lethal doses are *measured* (1.55
and 1.0 mM, `stress_panel.py:366` and `:407`) and override the multiple. Four of the
twenty-four declared `REPORTERS` load on ESR — counted by importing the module, not by
reading the file: `STRE-general` has `module="ESR"` outright, and `HSE-heat`,
`STRE-osmotic` and `CSRE-carbon` carry it in `also_reads` at 0.35, 0.4 and 0.25. So across
the panel's own dose design a whole shared arm carries **no agent-specific information at
all**.
`SOURCE_AUDIT.md` already grades `_GENERAL_POTENCY` **8/40** and files it under *"asserted
but dressed up"*.

**What to measure.** Clone **HSP12- or CTT1-mCitrine** — the `STRE-general` reporter the
panel already declares at `stress_panel.py:506` — and run it on the **same DTT and H₂O₂
ladders** that produced `outputs/sensor_characterisation.csv`. Fit both with
`generator.panel_calibration.fit_dose_response` and set
`_GENERAL_POTENCY = ec50_STRE / ec50_agent` from the measurement.

**Free, today, with no construct:** assert `_GENERAL_POTENCY < ceiling × _LETHAL_MULTIPLE`.
As shipped that is `2.5 < 1.65`, which is false — so the guard fails immediately, and it
fails for the right reason: the ESR arm is unreachable by construction at the ladder the
generator itself recommends.

**What confirms, what refutes.** Two measured EC50 ratios that differ between DTT and H₂O₂
refute the single asserted `_GENERAL_POTENCY` and restore a dimension to the panel. Two that
agree, and agree with 2.5, confirm the constant — and would be the first evidence it has
ever had.

**What it cannot settle.** The other 23 stressors. Their EC50s are literature values and
their lethal doses are `3.0 × ec50` aliases rather than measurements, so a measured ratio on
two agents does not license the multiple on the rest; it licenses **one** ratio, and the
others stay asserted and should stay labelled as such.

---

## M9 — Export capacity, paired pellet and supernatant at two expression levels

*(added 2026-09-01, with the saturating export outlet; run order 6)*

> **NOT REACHABLE BY THIS TEAM.** Needs HPLC, quantitative MS or a calibrated immunoassay,
> paired pellet and supernatant, on a **secreted-product strain** at two expression levels.
> No such strain is in the collection. [§What each of the six needs](#what-each-of-the-six-needs-and-whether-this-team-can-run-it).

### What is blocked

`solve_pathway` carries a saturating export outlet on a terminal node and **refuses to use
it**, because both of its constants are unmeasured:

    secretion_vmax_mmol_per_gdcw_h    export capacity, mmol/gDCW/h
    secretion_km_mmol_per_gdcw        content at half capacity, mmol/gDCW

Until they exist, every secreted product — TAL, farnesene, resveratrol, any recombinant
protein — loads and then refuses, naming this file. That is better than the load-time
refusal it replaced, and it is not an answer.

### What to measure

Paired, **on the same culture, at two or more expression levels**:

1. **Intracellular content** of the product, absolute units per gram of dry cell weight.
   Washed pellet. HPLC for a metabolite, quantitative MS or a calibrated immunoassay for a
   protein.
2. **Supernatant accumulation rate**, mmol/gDCW/h. In chemostat this is supernatant
   concentration × D / biomass concentration, which is why a chemostat is the natural
   vessel: at steady state the rate is a division rather than a slope.
3. **mu**, recorded. The capacity is an absolute flux and does not carry its own growth
   scaling, so it is a claim at the dilution rate it was measured at.

Two expression levels is the **minimum, not a comfort margin**. One steady state fixes only
the ratio `vmax/(km + X)` and leaves both constants free along a ridge; two invert exactly:

    km   = X1*X2*(v2 - v1) / (v1*X2 - v2*X1)
    vmax = v1*v2*(X2 - X1) / (v1*X2 - v2*X1)

Two gene dosages, two promoters, or two dilution rates all supply the second state.

### What confirms, what refutes

The outlet's shape is a claim, and this measurement can refute it. Plot secreted flux
against intracellular content across the levels:

- **Saturating** — flux flattens while content climbs. Confirms the shape, and the two
  constants follow.
- **Straight through the origin** — flux proportional to content. Refutes it: the outlet
  should be first order, like the degradation term, and `secretion_km_mmol_per_gdcw` should
  be deleted rather than fitted.

The shape was chosen on Kastberg 2025, which measured intracellular and secreted product on
the same *K. phaffii* chemostat samples and found a 17× intracellular rise against a
secretome that did not move. That is a refutation of the first-order law rather than a fit
of the saturating one, so the saturating shape is **the surviving hypothesis, not a
confirmed one**.

### One point of this already exists

**Pfeffer 2011** (PMID 21703020) is exactly this experiment at one state: ³⁴S pulse labelling
on a *K. phaffii* chemostat at D = 0.10 /h, giving intracellular 91.8 µg/gYDM against qSec
50.7 µg/gYDM/h. Its three outlets close on their own arithmetic to 0.06%. **One more
dilution rate on that system completes the calibration**, which makes this the cheapest item
in the list to start and the only one with a published half.

### What this cannot settle

**It cannot separate export from degradation.** Pfeffer's own split is 58% degraded / 35%
secreted / 7% diluted, with degradation the faster half-life, so a capacity fitted to net
supernatant appearance absorbs intracellular proteolysis into its own `vmax`. Separating
them needs a labelled pulse-chase and a node this walk does not have.

## M8 — Phytoene, on HPLC runs that already exist

*(item 8 of `DISTANCE_TO_THE_VISION.md` §6; run order 1)*

> **NOT REACHABLE BY THIS TEAM.** Needs an HPLC with a **285 nm** channel and a phytoene
> authentic standard, run on either Elizondo's stored chromatograms or a β-carotene strain.
> Its **step 1 is a letter, not a culture**, and that part is reachable — see
> [§What each of the six needs](#what-each-of-the-six-needs-and-whether-this-team-can-run-it).

### What is blocked

`data/pathways/beta_carotene.toml` declares the phytoene node

```toml
rate_law = "passthrough"      # see notes: an assumption, and a visible one
measurable = false            # no phytoene number exists for the Elizondo strains
```

`passthrough` means the node holds no pool, so nothing dilutes out of it and the flux passes
through unchanged. **That is known to be wrong, not merely untested.** Chen et al. 2016
(PMID 27329233, PMC4915043, doi 10.1186/s12934-016-0509-4) quantify phytoene at **3.99%**
and neurosporene at **4.87%** of total carotenoid in an engineered *S. cerevisiae*, and
report CrtI conversion as rate-limiting. A node carrying 4% of the pathway's carbon holds a
pool.

Two consequences, both currently live:

1. **What the flux law predicts is the post-CrtI flux.** `BETA_CAROTENE_FLUX.alpha =
   1.008474e-3` mmol/gDCW/h per unit relative CrtE was fitted against
   `q_lycopene + q_betacarotene`, which excludes anything still sitting as phytoene. Any
   claim about *total pathway flux* is low by an unmeasured amount.
2. **`calibrated_kinetics()` raises rather than returning numbers.** Its message names the
   fix: phytoene and GGPP were never measured, so the synthase and the desaturase have a
   pinned rate and an unobserved substrate, and any `(vmax, km)` pair reproduces the data.
   `UNIDENTIFIABLE_STEPS = ("psy", "crti")`.

### What to measure

**Phytoene content, in the six Elizondo steady states or their equivalent, from the same
extract that gives lycopene and β-carotene.** Nothing new is grown.

### Run

1. Re-inject, or re-integrate the stored chromatograms from, the existing carotenoid
   extractions. **Phytoene is resolved by the same injection at 285 nm** — it is
   uncoloured, so a method reading only 450–478 nm never sees it. This is the whole
   experiment.
2. Obtain a phytoene authentic standard and build a calibration curve on the same column,
   gradient and detector. Without a standard the peak gives a ratio and not a content, and a
   ratio cannot enter `q / μ`.
3. Report phytoene as **nmol/gDCW**, matching the units of
   `data/carotenoid/elizondo2025_steady_states.tsv`, alongside lycopene and β-carotene from
   the same injection so the three sum to a total.
4. Do this in **every** state, not one. Six phytoene contents against six different pathway
   fluxes is what identifies two parameters; one content is one equation in two unknowns.
5. Report neurosporene if the method resolves it. Chen 2016 puts it at 4.87%, larger than
   phytoene, and it is a second undeclared node on the same branch.

### What confirms, what refutes

- **Refutes `passthrough`:** any phytoene content above the detection limit. The declaration
  says the pool is zero; a measurable pool falsifies it directly. On Chen's proportions, a
  phytoene content near 4% of the carotenoid total is expected.
- **Confirms `passthrough`:** phytoene below the limit of detection in all six states, with
  the limit reported in nmol/gDCW so it can be compared against the lycopene contents, which
  span 331–7533 nmol/gDCW. A detection limit above ~300 nmol/gDCW would confirm nothing —
  it would only fail to look.

### What the model does differently afterwards

- The phytoene node becomes `rate_law = "saturating"`, `measurable = true`, and the fitted
  `(vmax, km)` for CrtI join `pathway/calibrations.py` beside the cyclase pair.
- `calibrated_kinetics()` stops raising for `crti`. It keeps raising for `psy`, because GGPP
  is still unmeasured — **one measurement closes one of the two unidentified steps, not
  both**, and the refusal message should be narrowed rather than deleted.
- `α` is refitted against `q_phytoene + q_lycopene + q_betacarotene`, so "relative CrtE
  expression predicts pathway flux" becomes a statement about the whole pathway.
- The 14.2% leave-one-strain-out error on β-carotene is scored against measured β-carotene
  either way, so it is not expected to move. Predicting that it will not move, and then
  finding it does, would say the branch structure is wrong.

### Effort

One authentic standard, one detector channel, and analyst time. No cultivation. If the raw
chromatograms were stored with the full spectrum, step 1 costs an afternoon.

### What this cannot settle

Nothing about GGPP, and therefore nothing about the synthase. GGPP is a soluble prenyl
diphosphate at intracellular concentrations that
[`DISTANCE_TO_THE_VISION.md`](DISTANCE_TO_THE_VISION.md) §6 records as not existing in
usable form in any published work on this pathway.

---

## M5 — β-carotene content in a high-producing strain

*(item 5; run order 2)*

> **NOT REACHABLE BY THIS TEAM.** Needs HPLC with a β-carotene authentic standard,
> gravimetric dry weight, residual glucose and ethanol, on a **high-producing CEN.PK
> carotenoid strain**. Step 6 below explains why a plate reader cannot substitute: an
> absorbance assay reads *above* the ceiling in two of three strains while the ceiling is
> being obeyed. [§What each of the six needs](#what-each-of-the-six-needs-and-whether-this-team-can-run-it).

> ### SETTLED BY LITERATURE, 2026-08-28 — and the answer is that the ceiling is WRONG
>
> **This measurement no longer needs to be run to answer the question it was written for.**
> Arhar et al. 2024 (PMID 39215465) report **79 mg/gDCW of β-carotene** — HPLC on a C18
> column at 448 nm against an authentic standard, with lycopene affirmatively below the
> detection limit, on **gravimetric** dry weight. That is ~63× the 1.2483 mg/gDCW ceiling.
> Fathi 2021 (PMID 33332529) independently reports 46.5 mg/gDCW by HPLC on gravimetric DCW.
> `docs/HARD_TESTS.md` carries the quotes and the consequences.
>
> **Two rows that had looked like refutations are withdrawn or reclassified instead**, and
> both corrections matter for how this table is read in future:
> - **Bubphasawan 2024 (PMID 38710418), 71.80 and 60.50 mg/gDCW — STRUCK.** The paper was
>   **withdrawn by Elsevier**; Crossref carries an explicit retraction and the title is
>   prefixed "WITHDRAWN:". This repository was citing a retracted paper as evidence. The
>   authors' own replacement (PMID 40944773) reports 2.29 and 2.90 mg/gDCW — ~31× and ~21×
>   lower — by genuine HPLC.
> - **Olson 2016 (PMID 27423881) and López 2019 (PMID 31380362) — reclassified `A453-sum`.**
>   Both quantify by single-wavelength absorbance of a crude extract, so both report the sum
>   of the coloured carotenoids, and both derive dry weight from OD600.
>
> **What is still worth running, and it is a different experiment.** The ceiling's *level* is
> settled. Its **μ-dependence is not**: no paper in the refuting set reports a growth rate at
> all, and every one is a flask. The steps below therefore still stand as the way to get a
> chemostat number at a known μ — and running D = 0.18 would still convert
> `outputs/registered_prediction_D018.csv` from Tier 0 to Tier 3. But it would now be testing
> **whether the repaired model predicts content**, not whether the old bound holds. It does
> not.

### What is blocked

The model asserts a hard ceiling and has never tested it. For a terminal node fed by a
saturating step, `content = flux_out / μ` and `flux_out ≤ vmax_per_growth · μ`, so

```
content ≤ vmax_per_growth = 2.3252e-3 mmol/gDCW × 536.87 g/mol = 1.2483 mg/gDCW
```

and **μ cancels exactly**. The model therefore says no genotype at any growth rate exceeds
1.2483 mg/gDCW. A hundredfold increase in cassette dosage buys 0.06%:

| relative entry expression | predicted content, mg/gDCW |
| ---: | ---: |
| 1 | 1.1065 |
| 10 | 1.2369 |
| 100 | 1.2476 |
| ceiling | **1.2483** |

It is not a fitting artefact: leave-one-out over the six calibration states moves the
capacity between 2192 and 2695 nmol/gDCW (1.177 to 1.447 mg/gDCW), and the fit-scatter 95%
interval `capacity_ci95` is 1.7563e-3 to 3.0784e-3 mmol/gDCW, i.e. **0.9429 to 1.6527
mg/gDCW**. One calibration state sits at 96.7% of the point value. It is set by data, which
is what makes it the sharpest falsifiable claim in the package.

**Why the published literature has not already settled it**, which is the part to be careful
about. [`EXTERNAL_CAROTENOID_BOUND.md`](EXTERNAL_CAROTENOID_BOUND.md) §5 tabulates 28
published contents, and many sit far above 1.2483 mg/gDCW — Verwaal 2007 at 5.90 (PMID
17496128), Olson 2016 at 25.52 (PMID 27423881), López 2019 flask at 21.00 (PMID 31380362),
Arhar 2024 at 79.00 (PMID 39215465). None of them refutes *this* ceiling, for three separate
reasons, and all three must be closed by the new measurement:

1. **Analyte.** Four of those rows are *total carotenoid*, not β-carotene. Total carotenoid
   is an upper bound on the β-carotene fraction, so a total-carotenoid number above the
   ceiling is not a β-carotene number above the ceiling. **This repository disagrees with
   itself about López 2019 specifically**: `HARD_TESTS.md` §2 calls the 21 mg/gDCW "total
   carotenoid", while `EXTERNAL_CAROTENOID_BOUND.md` §5 records all four López 2019 rows
   with analyte `beta-car` and its §8 limitations list names only the *fed-batch* López rows
   as total. **RESOLVED 2026-08-28 by that one full-text read: `HARD_TESTS.md` was right and
   `EXTERNAL_CAROTENOID_BOUND.md` was wrong.** López 2019 quantified "the absorbance at 453
   nm of the hexane extracts" against a β-carotene standard, with **no chromatography**, so
   all four of its rows are a summed absorbance of the coloured carotenoids and none is a
   β-carotene number. They are relabelled `A453-sum` there. No culture needed to be
   inoculated, and the largest apparent refutation of the ceiling is withdrawn.
2. **Cultivation mode.** Every refuting row is flask or fed-batch. In those, μ at the moment
   of harvest is unknown, and 72–120 h harvests are late-exponential or post-diauxic. The
   model refuses to predict outside μ ∈ [0.101, 0.2543] at all, so it makes no claim about
   those cultures — the ceiling's algebra is μ-free but the calibration behind it is not.
3. **Genotype.** The capacity was fitted **flat** across three strains whose crtYB copy
   number is 1:2:3. crtYB is bifunctional — one protein running both the synthase step that
   sets the flux and the cyclase step that splits it — so a strain with more crtYB should
   have more cyclase capacity, and the model has no way to say so: `Genotype` carries only
   the entry-enzyme expression. This paragraph said until 2026-09-04 that "fitting a capacity
   per strain was tried and scored **−0.58** held out"; that was the score of a candidate that
   ALSO dropped the growth term, so it measured growth-rate-independence and not dosage.
   Fitted WITH the growth term, three capacities on six points give conditional state-LOO
   log RMSE
   <!-- audit:value table=outputs/carotenoid_fit.csv column="RMSE(log) LOO" row="model=capacity per strain x mu" -->0.205 — lower error than the constant baselines, a hair worse than the pooled
   fit. This conditions on measured lycopene, retains sibling-strain outcomes and compares
   families on the same cohort; it is not nested forward new-strain validation. The three
   full-data fitted capacities
   land <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=capacity spread across strains" -->1.1494x apart and
   non-monotone against a 3× copy span, correlating with measured CrtYB transcript at
   <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=pearson log mrna vs log capacity" -->+0.0593, and the two extra
   parameters do not survive a nested F test — F(2, 2) =
   <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F statistic" -->0.5666, p =
   <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F p" -->0.6383. Non-rejection with six states and two denominator degrees of freedom does
   not establish identical capacities, absence of a dosage effect or adequate power. The
   dosage series remains a proposal for independent measurement, not a conclusion supplied
   by these full-data conditional diagnostics.

### What to measure

**β-carotene content, specifically and not as total carotenoid, in a chemostat at a
dilution rate inside [0.101, 0.2543] /h under Elizondo's own cultivation conditions, in a
strain engineered well past 1.2483 mg/gDCW — ideally one carrying extra crtYB.**

### Run

1. **Strain.** A high producer in the CEN.PK background, so the host matches Elizondo's
   β-car2/3/4 (CEN.PK2-1c). If a strain with elevated crtYB copy number exists, use it: it
   tests the ceiling *and* the flat-capacity assumption in one vessel. Record the copy number
   or the relative crtYB mRNA — without it, a violation cannot be attributed.
2. **Medium — and this is where it is easiest to go wrong.** Match Elizondo's
   cultivation, which was glucose-**excess**: all six calibration states ferment at every
   rate, including μ = 0.101 where a glucose-limited culture makes no ethanol at all
   (`HARD_TESTS.md` §6 measures the model's stamped physiology at 4.1–6.4× off on uptake for
   exactly this reason). A glucose-**limited** Verduyn run is a different culture, and it
   breaks the comparison with the registered prediction in step 3 even though it is the more
   standard chemostat. Whichever is run, record the feed glucose concentration **and the
   residual glucose and ethanol at steady state**, so the regime is on the record rather
   than assumed.
3. **Dilution rate.** **D = 0.18 /h.** Inside the fitted window, and the rate at which
   `outputs/registered_prediction_D018.csv` already holds a written, unmeasured prediction —
   0.5354, 0.9095 and 1.0223 mg/gDCW for β-car2, β-car3 and β-car4. Running D = 0.18
   therefore settles M5 *and* converts the only registered prediction in the repository from
   Tier 0 to Tier 3, which nothing here has ever been.
4. **Steady state.** Sample after five volume changes, 5 / 0.18 = **27.8 h** at constant OD
   and constant off-gas. Confirm steadiness with two samples one residence time apart, not
   with one.
5. **Replicates.** Three independent vessels, not three samples from one vessel. One vessel
   sampled three times measures the pipette.
6. **Analysis. CHROMATOGRAPHY IS NOT OPTIONAL, and this is the step most likely to be
   economised away.** Extract and quantify by **HPLC** against a **β-carotene authentic
   standard**, reporting β-carotene, lycopene and (with M8) phytoene separately from the same
   injection.

   **A plate reader or spectrophotometer cannot do this job, even with a β-carotene
   standard.** Single-wavelength absorbance of a crude extract sums every *coloured*
   carotenoid — β-carotene **and lycopene**, plus γ-carotene and neurosporene. (Phytoene is
   colourless, λmax ≈ 286 nm, and is the one intermediate that does *not* interfere.) A pure
   standard only fixes the extinction coefficient; it does not confer specificity. Lycopene
   is the **direct precursor** here, and the model predicts it accumulating past the product:

   | strain | β-carotene | lycopene | an A453 assay reads | vs ceiling 1.2483 |
   | --- | ---: | ---: | ---: | :---: |
   | b-car2 | 0.5354 | 0.2408 | 0.7762 | below (0.62×) |
   | b-car3 | 0.9095 | 0.8609 | 1.7703 | **ABOVE (1.42×)** |
   | b-car4 | 1.0223 | 1.4503 | 2.4726 | **ABOVE (1.98×)** |

   So an absorbance run returns a number above the ceiling in two of three strains **while
   the ceiling is being obeyed**, and it reads as a refutation. That is exactly what happened
   to López 2019. Absorbance stays fine for *relative* work — tracking a time course, ranking
   strains, confirming a culture produces at all — but it cannot settle M5, and a mg/gDCW
   from it must never be compared against the ceiling. Report
   biomass as **gravimetric dry weight**, filtered and dried, not converted from OD.
   `EXTERNAL_CAROTENOID_BOUND.md` §4 declines to apply an OD-to-gDCW factor for exactly this
   reason, even though two independent papers imply 0.40–0.41 g/L per OD600.
7. Record the relative CrtE mRNA on the same sample, on the same normalisation as
   `data/carotenoid/elizondo2025_relative_mrna.tsv`. Without it the prediction has no input.

### What confirms, what refutes

- **Refutes the ceiling:** β-carotene, measured against an authentic standard, above
  **1.2483 mg/gDCW** at a dilution rate inside [0.101, 0.2543] with dry weight measured
  gravimetrically. Above **1.6527 mg/gDCW** it is outside even the upper fit-scatter bound
  and the refutation carries no argument about the fit's uncertainty.
- **Confirms it:** content at or below 1.2483 mg/gDCW in a strain whose entry expression is
  well above the fitted range [0.245, 1.000] — and the informative version of that is a
  strain with *more crtYB*, since extra crtE is exactly what the table above shows buys
  nothing.
- **Settles nothing either way:** a flask number, a total-carotenoid number, an OD-derived
  dry weight, or a content between 1.2483 and 1.6527 in a single vessel.

### What the model does differently afterwards

- **If refuted**, `NodeKinetics.vmax_per_growth` stops being a constant. The correct repair
  is not a larger constant but a capacity that reads the cyclase's own dosage, which means
  `Genotype` gains a second expression field and `beta_carotene.toml`'s note about crtYB
  bifunctionality stops being "flagged, not modelled". Note that this is a **structural**
  change: no re-fit of the present model can express it.
- **If confirmed**, the "within 10% of the hard ceiling" note that predictions already carry
  becomes a validated statement rather than a warning, and the model gains its first Tier 3
  claim.
- Either way the registered D = 0.18 prediction resolves.

### Effort

One chemostat train, ~28 h to steady state after the batch phase, three vessels. HPLC as in
M8. The dominant cost is having a high-producing CEN.PK strain to hand; if one does not
exist, this is a strain-construction project and drops below M4 in the run order.

### What this cannot settle

Whether the ceiling holds in flask or fed-batch, which is where every published high content
was made. That is M7.

---

## M4 — Cytosolic acetyl-CoA, one laboratory, one quench protocol

*(item 4; run order 3. Historically ranked highest value; the threshold lead was later refuted.)*

> **NOT REACHABLE BY THIS TEAM AS AN EXPERIMENT** — needs LC-MS/MS with authentic standards
> for four acyl-CoAs, on **SCKK006** and its parent. **Its step 1 is the exception**: a
> citation check at a desk, and the one piece of the six that can be done today.
> [§What each of the six needs](#what-each-of-the-six-needs-and-whether-this-team-can-run-it).

### What is blocked

Kocharin & Nielsen 2013 (PMID 23514405, PMC3610212) provides eleven condition summaries
with genotype fixed at SCKK006: four glucose, three ethanol and four mixed-feed states.
The feed contrast is descriptive; the released table has no cultivation/biological-replicate
identifiers or justified randomization scheme. It is not eleven independently held-out
biological experiments.

**The former universal null bar is withdrawn.** The old skills −0.100, −0.153, −0.140 and
−0.162 used a full-data mean comparator, while the alleged candidate null redrew random
covariates on every response shuffle. Its −0.046 / +0.077 percentiles did not establish
unidentifiability, a significance threshold or the number of vessels this experiment needs.
[`EXTERNAL_PRODUCT_VALIDATION.md` §3](EXTERNAL_PRODUCT_VALIDATION.md) retains the full
historical table and the corrected, not-yet-adopted candidate description: fold-training
mean comparison, a zero-skill constant, three negative nonconstant skills and all four
`inference_status=pending_exchangeability`. The random-design demonstration remains
separate; its new percentiles must not replace the old bar.

This measurement therefore needs traceable independent cultures and a declared analysis,
not permission obtained from a parameter-count threshold. The original proposed mechanism
was a **driving-force threshold**. Its failure below remains a separate result; withdrawing
the invalid null argument does not revive that lead. yeast-GEM's cytosolic thiolase
`r_0103` — the chemical analogue of PhaA, both condensing two acetyl-CoA into acetoacetyl-CoA
plus CoA — has, computed from the vendored ModelSEED formation energies through
`bridge/thermodynamic.py`,

```
dGr'0 = +38.07 kJ/mol          UPHILL
threshold = sqrt([acetoacetyl-CoA][CoA] / Keq) = 19,042 uM acetyl-CoA
                                at [CoA] = 100 uM, [acetoacetyl-CoA] = 1 uM
```

**and BOTH carbon feeds are below it, so this lead is refuted.** It read +15.61 kJ/mol and
221.3 µM until 2026-08-30, when `bridge/thermodynamic.py` was found to be reading a kcal/mol
table under `pytfa`'s kJ/mol default. eQuilibrator, unaffected by that bug and matching the
literature at +24.96 kJ/mol, puts the threshold at 1.4 mM and also leaves ethanol below. The
count of mechanisms that fail to explain the 4.24x is therefore six, not five. That is why
five earlier mechanisms
failed — FBA on the product, GEM precursor ceilings, GEM precursor throughput, E-Flux, and
an mRNA-reading law are all arguments about how much the network *can* carry, and a rate cap
is blind to a driving force.

**And it rests on two numbers from two laboratories.**

| feed | measured | source, as recorded in `scripts/thiolase_threshold.py` |
| --- | ---: | --- |
| glucose-limited | 0.0199 μmol/gDW | Kolbeinsen & Bruheim 2021, supplementary Table S1 |
| ethanol-limited | 0.85 μmol/gDW | Kozak & van Rossum 2016 |

A **42.7-fold** difference, across two quench protocols. **Neither source carries a PMID or
a DOI anywhere in this repository**, so the first action below is a citation check and not a
culture. Nothing in this specification depends on which paper is right; it depends on the
two numbers never having been made the same way.

### What to measure

**Acetyl-CoA, and the two other participants in the thiolase reaction, in glucose-limited
and ethanol-limited chemostats at matched dilution rate, in one laboratory, with one quench
protocol, in the same analytical batch.**

### Run

1. **Before anything: resolve the two citations.** Record PMID or DOI, the strain, the
   medium, the dilution rate, the quench and extraction protocol, and whether the reported
   value is a whole-cell or a compartment-resolved figure, for both. If either turns out not
   to be a chemostat at D = 0.05 in Verduyn medium in CEN.PK113-7D — which is what
   `thiolase_threshold.py` asserts in a comment — the splice is worse than documented and
   this experiment gets more urgent, not less.
2. **Strain.** SCKK006, Kocharin & Nielsen 2013's PHB strain, so the acetyl-CoA measurement
   is in the organism whose flux is being explained. Run its parent alongside: the threshold
   is a claim about the *host* pool, and the PHB pathway itself drains it.
3. **Medium and feeds.** Verduyn mineral medium, carbon matched at **0.666 Cmol/L** exactly
   as Kocharin ran it: glucose 20 g/L; ethanol 15.32 g/L. Adding the 1:2 mix (6.35 g/L
   glucose + 10.21 g/L ethanol) is worth the third train if vessels allow — it is the feed
   with the *highest* PHB flux, 0.01065 mmol/gDCW/h, and if it also has the highest
   acetyl-CoA the threshold gains a third point rather than a second.
4. **Dilution rate.** **D = 0.05 /h**, matched across feeds, because that is the rate the
   two published acetyl-CoA numbers are recorded at in `thiolase_threshold.py` — and the
   whole purpose of this experiment is to make the new pair comparable to the old one. It is
   also where the feed contrast in PHB flux is largest: glucose 0.00251 against ethanol
   0.00961, a factor of **3.83**, falling to 2.41 at D = 0.10. The cost is the residence
   time: 5 / 0.05 = **100 h** per vessel against 50 h at D = 0.10.
   **If vessel-hours are the binding constraint**, run D = 0.10 instead and say so — the
   comparison against the published pair then carries a growth-rate difference as well as a
   protocol difference, which is a weaker but still decisive test of the 42.7-fold gap.
   **If vessels allow, run both rates**: two dilution rates crossed with two feeds separates
   a feed effect on the acetyl-CoA pool from a growth-rate effect, which one rate cannot.
   D = 0.20 is excluded on ethanol either way, because the culture **washes out** there.
5. **Replicates.** **Three independent vessels per feed**, sampled once each at steady
   state. Three samples from one vessel measure the sampling, not the culture. The design
   must distinguish "these two feeds differ" from "these two runs differ", and with two
   feeds and three vessels each that is a two-sample comparison on n = 3. Preserve vessel/
   cultivation identifiers, treatment assignment and any randomization blocks with the
   measurements; feed labels and condition-summary rows alone do not establish those units.
6. **Quench.** One protocol, applied identically to both feeds, in the same week, by the
   same operator, with the two feeds interleaved in the analytical batch rather than run as
   two blocks. The specific protocol matters less than that it is one protocol — the entire
   defect being repaired is that the existing pair is not. Whatever is chosen, pin it in the
   methods with a citation, and report the acetyl-CoA recovery of a spiked internal standard
   so the next laboratory can compare.
7. **Measure four acyl-CoAs in the same extract, not one.** Acetyl-CoA, **free CoA**,
   **acetoacetyl-CoA**, and (for the PHB strain) 3-hydroxybutyryl-CoA. This is the single
   largest improvement available and it is nearly free once the method exists — see the box
   below.
8. **Report per gDCW, and report the biomass gravimetrically.** μmol/gDW is the unit both
   existing numbers are in and it is the unit in which the comparison can be made without a
   volume convention at all (step 4 of "what confirms").
9. **Record the PHB content on the same samples.** Then the threshold and the flux it is
   supposed to explain come from one culture, and the argument stops being a splice in a
   second place.

> **Why measuring CoA and acetoacetyl-CoA matters as much as acetyl-CoA.**
> The threshold is `sqrt([acetoacetyl-CoA][CoA] / Keq)`, so the other two participants set
> it. **Under the historical, erroneous energy**, the grid [CoA] 30/100/300 μM crossed with
> [acetoacetyl-CoA] 0.3/1/3 μM gave **66.4 to 663.8 μM**, a ten-fold band, and the
> glucose-below/ethanol-above ordering held in **8 of those 9 cells** rather than 9. It failed
> at [CoA] = 300 μM, [acetoacetyl-CoA] = 3 μM, where the 663.8 μM threshold exceeded the
> ethanol state at 425 μM. This is retained sensitivity history, not the current threshold.
> Measuring the cofactors narrows concentration uncertainty; it cannot undo the energy-unit
> correction or guarantee the proposed ordering.

### What confirms, what refutes

> **Outcome 2 already happened, on 2026-08-30, without the experiment.** This list was
> written to say what a measurement would have to show. Criterion 2 — *"`A_glc` and `A_eth`
> both above the threshold, or both below it"* — is now met by the existing numbers, because
> the threshold moved rather than the pools: correcting a kcal/kJ error in
> `bridge/thermodynamic.py` took it from 221 µM to 19 mM, and 10 µM and 425 µM are both
> below that. eQuilibrator, independent of that error, agrees. So the design below is kept
> as the record of what was asked, and its own words are what closes it: *"the mechanism is
> dead, and the environment axis stays refuted with one fewer lead — which is a real result
> and should be written up as one."* It is written up in `docs/FINDINGS.md`.

Let `A_glc` and `A_eth` be the measured acetyl-CoA in μmol/gDW, and let the threshold be
computed from the measured CoA and acetoacetyl-CoA in the same extract.

1. **Supports the specified threshold ordering, not an identified rate mechanism:**
   `A_glc` below the threshold and `A_eth` above it, with the gap larger than the combined
   replicate scatter, in one laboratory with one protocol. The original stronger criterion
   also asked that the mix land above and the acetyl-CoA ordering match PHB flux (glucose
   0.00251 < ethanol 0.00961 < mix 0.01065 mmol/gDCW/h at D = 0.05). Even that agreement
   would not establish that thermodynamic feasibility determines the realized rate.
2. **Refutes the proposed straddle:** `A_glc` and `A_eth` both above the threshold, or
   both below it. A difference unresolved against replicate scatter also does not establish
   a straddle. The original proposal additionally called any of these outcomes proof that
   the 42.7-fold splice was a protocol artefact and the whole environment axis refuted.
   Those conclusions do not follow. The threshold mechanism's failed explanation remains;
   the cause of the cross-laboratory difference needs its own evidence.
3. **Refutes it in the most likely way:** `A_eth / A_glc` far smaller than 42.7 but still
   above 1. Then the direction survives and the magnitude does not, and the question becomes
   whether the surviving ratio still straddles the threshold. This is the outcome the
   design must be able to resolve, which is why the cofactors are measured.
4. **The volume-free form of the test, which avoids one convention entirely.** Converting
   μmol/gDW to molar needs a cytosolic volume, and 1.0, 2.0 and 2.7 mL/gDW are all defensible
   — `thiolase_threshold.py` prints all three because the choice is a convention rather than
   a measurement. Stated in μmol/gDW, the threshold is **19.0, 38.1 and 51.4 μmol/gDW** at
   those three volumes. Neither feed reaches any of them: the larger pool, 0.85 μmol/gDW, is
   22× below the lowest. Under the pre-2026-08-30 energy the threshold was 0.221, 0.442 and
   0.597 μmol/gDW and the existing pair straddled all three, which is what made this
   measurement worth designing. It is no longer the test it was — a measurement here would
   have to move the pool by more than an order of magnitude to change any verdict.

### What the model does differently afterwards

- **The thermodynamic layer gets the one input it lacks.** `bridge/thermodynamic.py`
  already computes ΔGr'° for any reaction from vendored tables; what it cannot do is say
  whether a step *runs*, because that needs measured concentrations and there are none. With
  them, the solver refuses a step whose ΔG at the implied concentrations is positive instead
  of returning a number — a refusal, and one that says what would lift it: a higher measured
  precursor concentration. Whether such a gate is already in the chain is a code question;
  this measurement is what makes it mean anything, and without it a gate can only assert.
- Measured concentrations can test a specified thermodynamic condition; they do not by
  themselves identify a rate law or validate product predictions. The original proposal
  said the threshold could bypass a three-parameter +0.077 null bar by clearing a
  one-parameter bar. **Both bars were invalid for this purpose and are withdrawn.** A
  mechanistic candidate still needs its own declared target, independent-unit validation
  and any justified null, including selection if choices are data-driven. The failed
  threshold lead above is not rescued by having fewer fitted parameters.
- The concentration-based flux law of `DISTANCE_TO_THE_VISION.md` §6 item 9,
  `v = kcat·E·[S]ⁿ/(Kⁿ + [S]ⁿ)`, becomes fittable on the PHB axis. It stays blocked on the
  carotenoid axis, where intracellular FPP and GGPP concentrations do not exist in usable
  form.

### Effort

The largest item here. Two or three chemostat trains, three vessels each, 100 h to steady
state per vessel at D = 0.05 plus the batch phase; an LC-MS/MS acyl-CoA method with authentic
standards for four compounds; and one operator holding the quench constant. If the acyl-CoA
method does not already exist in the laboratory, method development dominates the cost — and
it should be validated on a single culture before any of the comparison vessels are
inoculated, since the entire result is the difference between two numbers made the same way.

### What this cannot settle

**Compartmentation.** A whole-cell extract cannot separate cytosolic acetyl-CoA from
mitochondrial, nuclear and peroxisomal pools, and the thiolase argument is about the
cytosolic one. Whether either existing number is compartment-resolved is **not recorded in
this repository**. If the extract is whole-cell, say so, and treat the result as an upper
bound on the cytosolic pool — which is the conservative direction for the glucose state and
the *un*conservative direction for the ethanol state.

**The energy's own uncertainty.** This section was written when the energy was +15.61 kJ/mol
and is kept because its sensitivity argument is what the refutation turned on. The threshold
scales as `exp(dG0 / 2RT)`, so `d ln(threshold)/d(dG0) = 0.1984` per kJ/mol at 30 °C — which
is why a 22.5 kJ/mol unit error moved the threshold by a factor of 86, from 221.3 μM to
19,042 μM, and took the lead with it. At the old energy an error of ±5 kJ/mol moved the
threshold from 221.3 μM to the range **82.1 to 596.6 μM**, and ±8 kJ/mol to **45.3 to
1081.8 μM**. That band is wider than the cofactor grid's and it
cannot be narrowed by any culture — it is narrowed by a measured equilibrium constant for
this reaction, or by accepting the ordering rather than the absolute value as the claim.

---

## M6 — A third dilution rate, at D = 0.05

*(item 6; run order 4)*

> **NOT REACHABLE BY THIS TEAM.** Needs HPLC as M5, the three Elizondo β-car strains, and a
> route to μ = 0.05 /h — which even §P4's flask does not supply, since all six of its
> in-range carbon × temperature combinations land at or above 0.101 /h. [§What each of the six needs](#what-each-of-the-six-needs-and-whether-this-team-can-run-it).

### What is blocked

The cyclase law is `vmax = capacity · μ^β`. Fitted freely on Elizondo's six states, β = 1.041
with residual RMSE 0.0777 in log rate, and the **95% profile-likelihood interval is
[0.53, 1.78]** — 1.25 wide, 3.6× the leave-one-out resampling spread of [0.85, 1.20], which
is not a confidence interval and must not be quoted as one. So μ^0.6 and μ^1.7 are both
inside; 0 and 3 are both outside.

The consequence a user meets is a refusal, quoted exactly:

```
node 'lycopene': growth rate 0.05 /h is outside the range its kinetics were fitted
over, [0.101, 0.2543]. The rate law is `vmax = capacity * mu` fitted over that window
and its growth exponent has a 95% interval of [0.53, 1.78] even inside it --
extrapolating a law that loose is not a small error. Measure at this rate, or widen
the range deliberately and record that you did
```

**The model refuses both D = 0.05 and D = 0.35, and that is exactly why they are
informative:** a measurement is worth most where the model declines to answer.

### What a point inside the anchors buys, and what one outside buys — measured

`DISTANCE_TO_THE_VISION.md` §6 and `WHAT_IS_LEFT.md` both state that an interior point moves
the 95% width "from 0.46 to 0.47" and that a point outside "halves it". **No script or test
in this repository computes 0.46 or 0.47**, and those documents' own profile interval,
[0.53, 1.78], is 1.25 wide rather than 0.46. So the claim was re-derived for this document
by simulation, using the same profile-likelihood machinery as
`tests/test_growth_exponent_identifiability.py`. Method in the last section.

Median 95% profile width on β, over Monte Carlo draws, noise sd 0.130 in log rate:

| design | states added | median width | vs 1.225 |
| --- | :---: | ---: | ---: |
| the existing six states | — | **1.225** | — |
| + D = 0.18, three strains *(between)* | 3 | 0.680 | 0.56× |
| + D = 0.35, one strain *(outside, high)* | 1 | 0.720 | 0.59× |
| + D = 0.35, three strains | 3 | 0.60–0.64 | ~0.51× |
| + D = 0.05, one strain *(outside, low)* | 1 | 0.640 | 0.52× |
| + **D = 0.05, three strains** | 3 | **0.40–0.44** | **~0.34×** |

Three things this says that the prose did not:

1. **The direction of the published advice is right and its magnitude is wrong.** A point
   outside does not merely halve the interval; three states at D = 0.05 cut it to about a
   third. And an interior point is *not* worthless — it narrows by 1.8×, because three more
   observations raise the residual degrees of freedom from 3 to 6 whatever their placement.
   "Worth nothing" overstates it.
2. **D = 0.05 and D = 0.35 are not interchangeable, and the repository's "0.05 or 0.35"
   should read "0.05".** D = 0.35 buys barely more than an interior point, despite a lever
   arm of |log(D/0.1603)| = 0.781 against 0.116 for D = 0.18. The reason is not statistical
   but mechanistic, and it is measurable directly — the sensitivity of the observable to β,
   `|∂ log q_βcar / ∂β|`, averaged over the three strains:

   | D | 0.05 | 0.101 | 0.18 | 0.254 | 0.35 |
   | --- | ---: | ---: | ---: | ---: | ---: |
   | sensitivity to β | **2.904** | 1.977 | 1.192 | 0.785 | **0.482** |
   | cyclase share of pathway flux | 0.12–0.36 | 0.23–0.59 | 0.39–0.71 | 0.51–0.75 | 0.60–0.77 |
   | [lycopene] / Km | 6.1–29.6 | 2.0–12.7 | 0.8–5.7 | 0.5–3.3 | 0.3–1.9 |

   At high μ the lycopene content (`flux/μ`) falls below Km, the cyclase runs unsaturated
   and converts nearly everything that arrives, so its rate reports the *flux* and stops
   reporting `vmax` — and therefore stops reporting β. At low μ the pool saturates the
   cyclase and the rate is `vmax` almost exactly. **A state at D = 0.05 carries six times
   the information about β that a state at D = 0.35 does.**
3. **Run all three strains, not one.** One state at D = 0.05 gives 0.640; three give ~0.42.
   The three strains sit at fluxes spanning 3.5×, and it is the flux axis rather than the
   growth axis that carries the identifiability.

### Run

1. Three strains — β-car2, β-car3, β-car4 or their equivalents — in glucose chemostats at
   **D = 0.05 /h**, matching Elizondo's cultivation in every other respect.
2. Steady state at five volume changes: 5 / 0.05 = **100 h**, so ~4.2 days per vessel after
   the batch phase. This is the price of the low rate and it is unavoidable.
3. Sample lycopene, β-carotene and (with M8) phytoene from the same injection, against
   authentic standards, in nmol/gDCW with gravimetric dry weight.
4. Record relative CrtE mRNA on the same normalisation as the existing states.
5. Replicate at the vessel level, not the sample level.

### What confirms, what refutes

- **Confirms:** the three new states fall within the LOSO band of the extrapolated fit — a
  typical fold error of `exp(0.2002) = 1.22×` on flux, `exp(0.1858) = 1.20×` on the cyclase
  rate — and the refitted β interval narrows around 1 without moving the point estimate
  outside [0.85, 1.20].
- **Refutes the power law:** the refit lands with β outside [0.53, 1.78], or the new states
  sit outside the extrapolation band, or the profile becomes bimodal. **This is the more
  valuable outcome and the simulation above cannot represent it**, because it generates data
  from the fitted law. A misfit at D = 0.05 says the μ-scaling is the wrong functional form
  outside the window, which is precisely what nobody knows.
- **A second, free test rides along.** At D = 0.05 the predicted content is highest, since
  `content = flux/μ` in the passthrough limit. If any state exceeds 1.2483 mg/gDCW, M5 is
  settled at the same time.

### What the model does differently afterwards

- `ELIZONDO2025.growth_rate_range` widens to `(0.05, 0.2543)` and the refusal at D = 0.05
  stops firing. D = 0.35 stays refused, correctly, until someone measures there.
- The interval quoted inside the refusal message itself narrows, so the message stops
  arguing from a width of 1.25.
- Whether `vmax_per_growth` is a capacity or a rate constant becomes testable: a *content*
  ceiling and a growth-proportional *vmax* are arithmetically the same statement over
  [0.101, 0.254], and they separate as the window widens.

### Effort

Three vessels × 100 h to steady state, plus the batch phase. No new analytical method beyond
M5's. The long residence time is the whole cost.

### What this cannot settle

Anything about batch. D = 0.05 extends the window downward; μ ≈ 0.4 /h is in the other
direction and is M7.

---

## M7 — A batch calibration state

*(item 7; run order 5)*

> **NOT REACHABLE BY THIS TEAM, AND IT READS AS THOUGH IT WERE.** "No chemostat" was never
> its cost: it needs ≈36 HPLC injections for β-carotene, lycopene and phytoene, gravimetric
> dry weight, residual glucose and ethanol, and relative CrtE mRNA — all on a
> **carotenoid-producing strain**. See
> [§M7 is not the no-chemostat option](#m7-is-not-the-no-chemostat-option-and-reads-as-though-it-were).

### What is blocked

Batch is the most common question a user has and the model refuses it. Asked with no
dilution rate, the context supplies μ = 0.4 /h and the same guard fires:

```
node 'lycopene': growth rate 0.4 /h is outside the range its kinetics were fitted
over, [0.101, 0.2543]. ...
```

Of eight questions `DISTANCE_TO_THE_VISION.md` §2 puts to the model, one gets a real answer,
three are refused, and shake flask is one of the three.

### Why this is not simply "run a flask"

The whole solver rests on a steady-state identity. For a diluted pool,

```
d[P]/dt = v − μ[P] = 0    ⇒    q = μ · [P]
```

and every content in `data/carotenoid/` is `q / μ` on that basis. **In batch that identity
does not hold**, because `d[P]/dt ≠ 0` while the culture accelerates. A single end-point
flask measurement therefore cannot calibrate this model, however carefully it is done: it
gives a content with no rate attached and no growth rate at the moment it was attached.

There are two honest ways out, and the specification takes the first:

1. **Balanced exponential growth**, where the *content* has stopped changing even though the
   culture has not. Then `q = μ[P]` is recovered with μ measured from the OD slope. This has
   to be demonstrated, not assumed — by showing content is flat across at least three
   sampling points inside exponential phase.
2. Fit the full ODE to a content time course. More information, much more work, and a
   different solver.

### Run

1. **Strain and medium.** The same strains as the chemostat states, in the same medium, so
   the only variable is the cultivation mode. Glucose in excess, which is what a flask is.
2. **Culture.** A shake flask with a headspace ratio and shaking rate recorded explicitly —
   `EXTERNAL_PRODUCT_VALIDATION.md` §2 records a case where the *only* thing separating a
   refuted constant from a corroborated one was reading the flask geometry out of the
   methods. 150 mL in a 500 mL flask at 140 rpm is an aerated culture, whatever the paper
   calls it.
3. **Sampling.** At least **four points inside exponential phase**, before glucose is
   exhausted and well before the diauxic shift. At each point measure: OD600, gravimetric
   dry weight on at least two of them, β-carotene, lycopene, phytoene, residual glucose and
   ethanol, and relative CrtE mRNA on at least two.
4. **μ.** From the slope of ln(OD) across the sampled window, with its standard error
   reported. μ is now an *estimate* rather than a pump setting, and that is the main thing
   batch loses relative to a chemostat.
5. **Expect fermentation and record it.** van Hoek 1998 (PMID 9797269, Table 1) has ethanol
   below detection up to D = 0.25 /h, appearing at D = 0.28 (0.11 mmol/gDW/h) and reaching
   13.9 mmol/gDW/h by D = 0.40. A batch culture at μ ≈ 0.4 is fully respiro-fermentative.
   `CRITICAL_GROWTH_RATE_PER_H = 0.28` is the model's name for that threshold.
6. **Replicates.** Three independent flasks. Batch scatter is larger than chemostat scatter,
   which is the second thing batch loses.

### What confirms, what refutes

- **Prerequisite, and it can fail:** β-carotene content flat across the exponential
  sampling points — within the assay's own replicate scatter — establishes that a
  pseudo-steady state exists at all. **If content is still rising at the last exponential
  point, the experiment has refuted the premise rather than produced a calibration state**,
  and the honest output is that batch needs the ODE route.
- **Confirms the extrapolation:** measured content at μ ≈ 0.4 within the LOSO band
  (`1.20×` on the cyclase rate) of the law extrapolated from [0.101, 0.2543].
- **Refutes it:** measured content outside that band. The interesting direction is *high*:
  `content ≤ vmax_per_growth` puts a hard 1.2483 mg/gDCW on batch too, and published flask
  contents in other strains reach 5.90 (Verwaal 2007), 21.00 (López 2019) and 79.00 (Arhar
  2024) mg/gDCW. **A flask content above 1.2483 mg/gDCW in these strains refutes the ceiling
  and settles M5 in the cheapest possible vessel** — with the caveat that it does so in a
  mode the calibration never covered, which is a weaker refutation than M5's chemostat.

### What the model does differently afterwards

- `growth_rate_range` widens upward, and `Environment(growth_rate_setpoint_per_h=None)` returns a number
  instead of raising — which converts three of the eight questions in
  `DISTANCE_TO_THE_VISION.md` §2 from refusals into answers.
- The `fermentative` flag becomes meaningful for the culture being predicted. Today it is
  read from van Hoek's **glucose-limited** table while the calibration states were
  glucose-**excess** and fermenting at every rate, so it reports `False` for cultures that
  were visibly fermenting. A batch state is glucose-excess by construction and is the right
  place to fix it.
- If the flat-content prerequisite fails, the correct change is the opposite of a widened
  window: `solve_pathway` should refuse *non-chemostat* environments explicitly, saying that
  the steady-state identity does not hold, instead of refusing them incidentally through a
  growth-rate guard.

### Effort

The cheapest cultivation here: three flasks, one day, no chemostat. The cost is analytical —
four time points × three flasks × three carotenoids × a dry weight, roughly 36 HPLC
injections plus RT-qPCR.

### What this cannot settle

The environment axis. A batch culture varies μ, carbon status and oxygen status all at once,
so it cannot separate them — which is why the chemostat exists.

---

## Where the numbers in this document came from

Everything quoted above is either read out of this repository or computed for this document.
Nothing is estimated.

**Read from the repository:**

| number | source |
| --- | --- |
| dGr'0 = +38.07 kJ/mol; threshold 19,042 μM; both feeds below | `python3 scripts/thiolase_threshold.py`, run |
| 0.0199 and 0.85 μmol/gDW | `scripts/thiolase_threshold.py::MEASURED_UMOL_PER_GDW` |
| ceiling 2.3252e-3 mmol/gDCW; `capacity_ci95`; window [0.100988, 0.254320] | `ystwin.kinetic.carotenoid.ELIZONDO2025`, printed |
| α = 1.008474e-3; LOSO rmse 0.2002; skill 0.6262; expression range [0.2452, 1.000] | `ystwin.pathway.calibrations.BETA_CAROTENE_FLUX` |
| refusal strings, and 1.0710 mg/gDCW at D = 0.18 | `predict_product` called at D = 0.05, 0.101, 0.18, 0.254, 0.35 and `None` |
| registered D = 0.18 predictions | `outputs/registered_prediction_D018.csv` |
| historical alleged null bars −0.046 (95th) and +0.077 (99th) at k = 3 — withdrawn, not replaced | old `scripts/flux_null_baseline.py`; full-data comparator and random-design ensemble, not a candidate significance test (M4) |
| PHB fluxes, feeds, and the washout at D = 0.20 on ethanol | `data/phb/kocharin2013_chemostat_states.tsv` |
| profile interval [0.53, 1.78]; resampling spread [0.85, 1.20] | `tests/test_growth_exponent_identifiability.py`; `docs/research/KINETIC_FIT.md` §6 |
| phytoene 3.99%, neurosporene 4.87% | `data/pathways/beta_carotene.toml` notes, citing PMID 27329233 |
| ethanol at D = 0.25 / 0.28 / 0.40 | `ystwin.generator.culture.CRITICAL_GROWTH_RATE_PER_H` docstring, citing PMID 9797269 |
| the four constructs and their two agents (DTT, H₂O₂), and the two anchor transcripts | `qpcr.py::STRESSOR_FOR_CONSTRUCT`, `::TARGET_FOR_CONSTRUCT` — the premise of §What this team can measure |
| every plate is OD600 + mCitrine on a Synergy H1, and no chromatogram exists | `DATA_INVENTORY.md`; `data/plates/manifest.csv` |
| 91.4% gDNA, +0.13-cycle margin, 7 of 36 negative; TRX2 at 8.60 cycles and 0.3% | `G4_ANCHOR.md`, `G4_STATISTICS.md` (B2, B3) |
| SNR 0.39 / 1.19 → 80 and 9 replicates; `unsalvageable and must be re-run` | `G4_ANCHOR.md`; the `replicates_needed` field of `gates/g4_anchor.py`, `ceil(n_reps · (target_snr/snr)²)` at its line 138 — **not** `analysis/power.py::replicates_needed`, an unrelated simulation (B2, B3) |
| primer sequences and amplicons for the six pairs; efficiencies 96 / 93 / 91% | `UPR_ANCHOR.md` §7.2, §5.3 (B2) |
| single-well error 3.25% / 1.55% / 1.31% at the auto window | `outputs/estimator_accuracy.csv`, rows `is_auto_window=True` (Free §1) |
| `a` = +18.2 RFU/OD, CI [−6.1, +42.4], 0.42% of a reporter well | `DATA_INVENTORY.md`, `scripts/measure_autofluorescence.py` (Free §4) |
| Spearman +0.891; partial +0.839 → +0.494; 15 of 16 blocks positive | `AUTO_FEEDBACK_LOOP.md`, `outputs/afl_circuit.csv` (B5) |
| endpoint read, 24.9 °C, 28 conditions × 3 biological replicates | `data/crosstalk/SOURCE.md` (B6) |
| naive `RFU/OD` 3390 → 1354 against corrected −151 → +85 | `DISTANCE_TO_THE_VISION.md` (B6) |
| `_LETHAL_MULTIPLE = 3.0`, `_GENERAL_POTENCY = 2.5`, `ceiling = 0.55`; `_GENERAL_POTENCY` graded 8/40 | `stress_panel.py:307–308`, `:692`; `SOURCE_AUDIT.md` (B7) |
| 20 estimable folds, 7 clearing unadjusted, `(1/4)⁴`, `2·(1/2)⁴`, five / six / ten plates | `FINDINGS.md` |
| six of twelve carbon × temperature combinations inside [0.101, 0.2543] /h | `PROTOCOLS.md` §P4 (B5) |
| gain multiple 8.15, 96,459 paired readings, no roll-off to 99,998 RFU | `calib/gain_linearity.py`, `scripts/score_gain_linearity.py` (B1) |
| `od_linear_max` is a placeholder at 1.0; the 81% rejection | `tests/test_gain_linearity.py::test_od_linear_max_is_still_the_placeholder`; `PROTOCOLS.md` §P1 (B1) |

**Computed for this document**, all by running rather than by argument:

- **42.7-fold**, as `0.85 / 0.0199`.
- **1.2483 mg/gDCW** and its interval **0.9429–1.6527**, as `capacity × 536.87` on the point
  estimate and on `capacity_ci95`, with the molar mass from `beta_carotene.toml`.
- **The threshold in μmol/gDW** — 19.0, 38.1, 51.4 at 1.0, 2.0, 2.7 mL/gDW — as
  `19042e-6 M × volume`, which is `thiolase_threshold.py`'s conversion run backwards so the
  comparison can be made in the unit the measurement arrives in. (0.221, 0.442, 0.597 under
  the pre-2026-08-30 energy.)
- **`d ln(threshold)/d(dG0) = 0.1984` per kJ/mol**, from `threshold ∝ exp(dG0/2RT)` at
  RT = 2.5204 kJ/mol, and the resulting 82.1–596.6 μM at ±5 kJ/mol.
- **Steady-state times** as `5 / D` hours, five volume changes being a convention rather
  than a measured requirement.
- **The profile-width table in M6.** The `_predict`, `_fit` and `_profile` functions of
  `tests/test_growth_exponent_identifiability.py` were reproduced verbatim in a scratch
  script — the shipped solver pins the exponent at 1 by construction, so a study of the
  exponent cannot go through it. On the six real states this machinery returns β = 1.041,
  capacity 2.4581e-3, Km 5.3453e-4, RMSE(log) 0.0777 and a 95% interval of width 1.225,
  reproducing `KINETIC_FIT.md` §6 exactly, which is the check that it was reproduced right.
  Three simulated states were then appended at a candidate D — one per strain, at that
  strain's geometric-mean measured flux, since the shipped flux law has **no** growth term —
  with `q_βcar` generated from the fitted law times lognormal noise, and the profile recomputed
  at n = 9, k = 3. Reported as the median over 25–60 draws; the residual seed-to-seed scatter
  is about ±0.04 in width, which is why two of the rows carry a range. Noise sd 0.130 is the
  residual sd `KINETIC_FIT.md` §5 reports; the ordering is unchanged at sd 0.0777 (the
  three-parameter train residual, widths 0.36 / 0.56 / 0.56) and at a pessimistic 0.200
  (0.52 / 0.96 / 0.96).
- **The sensitivity table** `|∂ log q_βcar / ∂β|`, by central difference at h = 1e-4 on the
  fitted parameters, with the cyclase share and `[lycopene]/Km` from the same solve.
- **B1's `ceiling = tolerance / k`**, by calling `ODCalibration.linear_range_max(0.05)` at
  k = 0.2, 0.5, 1.0 and 2.0 and comparing against `0.05 / k`: equal to eight decimals in all
  four. The shipped code takes two steps and the second cancels the first.
- **The window the estimator rule actually picks**, by calling
  `default_activity_window_h(np.linspace(0, 4, 25))` → 1.333 h and
  `(np.linspace(0, 24, 145))` → 4.0 h, which is what selects the two
  `is_auto_window=True` rows quoted in *Free* §1.
- **Four of twenty-four `REPORTERS` load on ESR**, by importing `stress_panel` and counting
  `module == "ESR"` plus `"ESR" in also_reads` — one and three respectively. Reading the file
  by eye gives one.

**Three discrepancies found while writing this, all recorded rather than repaired** — this
document edits neither `src/`, `scripts/` and `tests/` nor the other documents:

1. **The "0.46 to 0.47" figures have no reproducer.** They appear in
   `DISTANCE_TO_THE_VISION.md` §6 and `WHAT_IS_LEFT.md` and nowhere else; no script or test
   computes them; and the interval those same documents quote is 1.25 wide. The M6 table
   above replaces them with a computation that can be re-run.
2. **López 2019's analyte is labelled two ways.** `HARD_TESTS.md` §2 calls the 21 mg/gDCW
   "total carotenoid"; `EXTERNAL_CAROTENOID_BOUND.md` §5 records all four López rows as
   `beta-car` while its §8 names the *fed-batch* rows as total. One full-text read of
   PMID 31380362 settles it, and it should happen before M5 is designed around it.
3. **The "~8% single-well error" in `AUTO_FEEDBACK_LOOP.md` is the superseded window rule.**
   That document (§"Why this dataset matters to the model") says the current rule "collapses
   onto its floor and sits at ~8%" on a four-hour run. It does not:
   `default_activity_window_h` on a 25-point, 4.00 h grid returns **1.333 h**, whose measured
   error in `outputs/estimator_accuracy.csv` is **3.25%**, and `reporter.py`'s own docstring
   states the 8% belongs to the one-sixth rule at 0.667 h. Both files are correct about the
   24 h figure of 1.31%. The consequence is only for the size of the claim — a longer run is
   worth **2.5×**, not 6× — and the recommendation stands either way.

A third gap, which is a genuine absence rather than a contradiction: **the two acetyl-CoA
sources carry no PMID or DOI anywhere in this repository.** They are named only in
`scripts/thiolase_threshold.py` as "Kolbeinsen & Bruheim 2021, supplementary Table S1" and
"Kozak & van Rossum 2016". Step 1 of M4 exists for that reason.
