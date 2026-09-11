# The latent stress state on real multi-stressor data

Every transfer number in this repository was fitted to `generator/stress_panel.py` and
scored on data the same file generated. [`NULL_RESULTS.md`](NULL_RESULTS.md) already
records that the headline one fails its rotated-subspace null in simulation — observed
+0.0499 against a null of +0.0362, p = 0.268. This file runs the same test on measurements
nobody in this project made, and reports what happens.

The wet-lab data here cannot run this test at all. It carries two stressors, DTT and H2O2,
and each is read only by its own dedicated construct pair — UPRE1/UPRE2 read DTT,
NativeYap1/AlteredYap1 read H2O2. Holding out a stressor therefore removes every channel
that reads its private modules, which is the zero-peer case the redundancy law says
recovers nothing. A public multi-stressor dataset is the only route.

---

## 1. Verdict

**The real-data transfer result fails its `rotated_subspace` null, the same way the
simulated one does, and by a wider margin.** Median score over 18 (held-out stressor, own
target) pairs is **+0.688** against a null median of **+0.695**, p = **0.990**. Skill over
the null is **−0.025**: arbitrary axes of the same rank carry slightly *more* transfer
information than the fitted ones. The circularity-restricted six-stressor panel gives
+0.605 against +0.608, p = 0.726. This is a second independent failure, on real
measurements, of the claim that the latent structure is doing the work.

**The redundancy law survives, and then loses its interpretation.** Recovery tracks peer
count on real data — Spearman **+0.734**, p = 0.001, and 0-peer pairs mean 0.280 against
6+-peer pairs mean 0.790, which is the same shape simulation reported. It clears a
permutation null comfortably. It does **not** clear a rotated-subspace null on the
restricted panel (+0.786 against +0.786, p = 0.871). So the law is a true statement about
which modules this panel covers twice, and not evidence about the fitted basis. Both
things at once.

**Two by-products.** The known-loadings readout `x = y (L⁺)′` beats `shuffled_loadings` on
real data — median Spearman 0.465 against a null of 0.091, p = 0.005 — which is the first
time "these channels identify these modules" has had a null anywhere in this repository,
and it passes. And a transcriptome cannot test 5 of the 24 modules at all: `redox`,
`peroxide`, `atp`, `ph` and `nadh` are metabolite pools with no regulon, and those are
exactly what the ratiometric sensors read.

---

## 2. The dataset, and why it beat the others

**Chosen: Gasch et al. 2000, GEO [GSE18](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18)**
("Yeast Stress Response"), PMID 11102521. 155 mappable arrays over nine array print runs,
2.7 MB of series matrices plus 7 MB of platform annotation, no credentials.

| Candidate | Distinct stressors | Dose or timecourse | Maps to `MODULES`? | Downloadable | Verdict |
| --- | ---: | --- | --- | --- | --- |
| **Gasch 2000 (GSE18)** | **8 usable of ~12** | timecourses, 5–19 arrays each, plus a ΔT ladder for heat | yes, via regulon membership | yes, GEO FTP | **chosen** |
| Causton 2001 | ~6 | endpoint | yes in principle | **no** — not in GEO, MBoC returns 403 to automated fetches | rejected |
| Brauer 2008 | 6 nutrient limitations | growth-rate grid | only `carbon`, `nitrogen`, `sulfur` — 3 of 24 | yes | rejected: too few modules, and limitation is not stress |
| Hackett 2016 / Xia 2022 | 5 / 1 limitations | chemostat grid | same 3 modules | partly behind publisher pages | rejected, same reason |
| O'Duibhir 2014 | 0 | deletion mutants | no | yes | rejected: not stressors |
| Newman 2006 | 2 media | endpoint protein abundance | no | yes | rejected: 2 conditions |
| Keren 2013 | 4 stress conditions of 10 | one scalar per promoter per condition | yes | yes | rejected: too few stressors, no dynamics |
| omniplate (Swain lab) | 0 | reporter timecourses | no | yes | rejected here — nutrient shifts, not stressors. It is still the right dataset for the *forecast*, per [`research/EXTERNAL_DATA.md`](research/EXTERNAL_DATA.md) §3 |

[`research/EXTERNAL_DATA.md`](research/EXTERNAL_DATA.md) §4 item 1 already records the
governing fact: **no public dataset combines yeast stress reporters, a dose ladder, and OD
plus fluorescence time series.** So the readout has to be a transcriptome, and among
transcriptomes GSE18 is the only single-lab, many-stressor yeast compendium with
per-condition labels.

**One dead end, recorded so nobody repeats it.** Harbison et al. 2004 (PMID 15343339) was
the first choice for the module labels. Its genome-wide location data is live at the Young
lab, but both the text distribution (`pval_by_gene_abbr.txt`, 6229 × 352) and the MATLAB
one (`Results_for_upcoming_papers_abbr.mat`) publish the p-value matrix with **no row or
column labels**; the labelled copy is a 40 MB legacy `.xls` needing `xlrd`, which is not a
pinned dependency. MacIsaac et al. 2006's own supplementary deposit is four `.doc` files
and its MIT host 404s. SGD carries the same binding calls with labels, so SGD is used.

---

## 3. Circularity

This is the part that decides whether the numbers above mean anything, because
`stress_panel.py` is parameterised from this literature and cites Gasch 2000 by name.
Grepping the file, it does so in exactly four places.

| Where | What it takes from Gasch 2000 | Held out? |
| --- | --- | --- |
| `STRESSORS["diamide"]` | the agent, its **EC50 of 1.5 mM** — literally Gasch's dose — and its 0.60 ESR arm | **yes**, `--independent` drops the stressor entirely |
| `STRESSORS["sorbitol"]` | the agent and its **EC50 of 1.0 M** — also Gasch's dose | **yes**, same |
| `MODULES["ESR"]` | jointly with Martinez-Pastor 1996 (PMID 8641288), which supplies the STRE element independently | **no** — reported separately instead, see below |
| `STRESSORS["menadione"]` | nothing. The comment records that *Gasch used 1 mM where the panel uses 100 µM*, i.e. a disagreement | n/a |

Unattributed but pervasive: the *pattern* that every stressor carries an ESR arm is the
Gasch finding. It cannot be held out, because it is in all 25 stressors' target lists. So
ESR is reported as its own row everywhere below, and the reader can discount it.

What is left genuinely independent after holding diamide and sorbitol out:

- the whole of `reporter_loadings` — Kuge & Jones 1994, Morgan 1997, Petrezselyova 2016,
  Platara 2006. No Gasch anywhere in it, and this is what §8's null tests.
- every module→transcription-factor assignment except ESR's — Cox & Walter, Boy-Marcotte
  1999, Owsianik 2002, Hahn 2006, Yamaguchi-Iwai 1996, Zhao & Eide 1997, and the rest.
- the `driven_by` cascade: `proteasome ← oxidative, heat`; `iron ← oxidative`;
  `xenobiotic ← proteasome`.
- the DTT and H2O2 parameterisation — Pincus 2010, Goulev 2017, Kritsiligkou 2021,
  Ayer 2013, MacGilvray 2020.
- six of the eight mapped stressors, and four of the six own-target modules that are not
  ESR.

**The labels are independent of the measurements by construction, and this is asserted
rather than claimed.** Module gene sets come from SGD's curated regulation records
restricted to **MacIsaac et al. 2006** (PMID 16522208) — conserved binding motifs over
ChIP data, no expression evidence of any kind. `load_regulons` raises if Gasch 2000
contributes a single record to any of the 32 transcription factors. It contributes **zero
records to all 32**. That matters: had the module definitions come from expression data
they would be partly Gasch's own clustering, and the whole exercise would be scoring the
panel against the paper that parameterised it.

**Verdict: partly circular, and the circular parts are separable.** The two Gasch-derived
stressors change the numbers by little and the conclusion by nothing — the transfer result
fails its null on both panels. This is not a fully independent validation of
`stress_panel.py`, because the ESR arm cannot be held out. It is an independent test of
the *estimator and the architecture*, which is what G-d is a claim about.

---

## 4. The mapping

### Conditions → stressors

GSE18's only condition labels are free-text sample titles, so every rule is a regex over
the title. Eight groups map onto panel stressors.

| Gasch condition | Panel stressor | Arrays used | Note |
| --- | --- | ---: | --- |
| 25→37 °C shock, plus the variable pre-shift ladder (17→37 … 37→37) | `heat` | 17 | the pre-shift ladder is a genuine ΔT dose series |
| growth to stationary phase in YPD, 2 h → 13 d | `glucose_starvation` | 17 | glucose exhaustion; the early rungs are near-zero dose |
| two DTT timecourses (2.5 mM, and a second series) | `DTT` | 10 | includes `dtt 000 min`, a zero-dose rung |
| 1 mM menadione | `menadione` | 9 | the panel's EC50 is 100 µM, so this is a high rung |
| nitrogen depletion, 30 min → 3 d | `rapamycin` | 8 | **matched on module, not on agent** — see below |
| constant 0.32 mM H2O2 | `H2O2` | 7 | |
| 1.5 mM diamide | `diamide` | 8 | Gasch-derived; dropped by `--independent` |
| 1 M sorbitol | `sorbitol` | 5 | Gasch-derived; dropped by `--independent` |

`rapamycin` is the one mapping that is not chemistry. Nitrogen depletion and rapamycin are
different agents; they are the same *module profile* — `nitrogen` 1.0, `ESR` 0.60, `carbon`
0.30 — and the panel has no nitrogen-starvation stressor. The positive control in §5 says
this mapping works: nitrogen depletion is the argmax for both the `GATA-nitrogen` channel
(+3.83, by a factor of four) and the `nitrogen` module (+0.56). Every other reading in this
file that depends on it is flagged.

### Conditions dropped

74 of 155 arrays, all with a reason recorded per array in
`outputs/external_conditions.csv`.

| Reason | Arrays |
| --- | ---: |
| mild or combined shift with no matching panel dose (25 °C, 29→33 °C, 33 °C ± sorbitol) | 15 |
| a channel gene has no probe on this array's print run | 15 |
| adapted steady state, not a shift; cold has no module (15–36 °C ladder, "N deg growth") | 13 |
| genotype perturbation, not a stressor (`msn2 msn4`, `yap1`, MSN2/MSN4 overexpression) | 13 |
| reference is a pool, not a paired unstressed control (the six carbon-source arrays) | 6 |
| no panel module: general amino-acid starvation is Gcn4, which `MODULES` does not carry | 5 |
| no panel stressor: hypo-osmotic shock has no module | 5 |
| flagged by the depositors ("1M sorbitol - 45 min (* centrifugation problem)") | 1 |
| adapted steady state, not a shock ("1M sorbitol vs. YPD - basal") | 1 |

The carbon-source arrays are the most interesting loss. Galactose, raffinose, ethanol,
sucrose and fructose against glucose would genuinely exercise the `carbon` module, and
re-referencing all six to the glucose array would recover them. They are dropped because
the deposited reference is a pool rather than a paired control, and because the panel has
no "alternative carbon source" stressor to map them to.

### Readouts → channels

The 19 transcriptional `REPORTERS` name their marker gene. That gene's log2 ratio *is* the
channel. Where the panel writes a slash — `HSP12/CTT1`, `PIR3/CRH1`, `PMC1/CMK2`,
`ADH2/ICL1`, `DAL5/MEP2`, `ANB1/DAN1` — it means either promoter would serve, so both
genes are read and averaged. `CUP1` is a tandem duplication and is read as `CUP1-1` plus
`CUP1-2`.

Two channels are dropped: **`Xbox-dna`** (RNR3 is absent from 55% of these arrays) and
**`PACE-proteasome`** (RPT1 from 16%), both because whole print runs omit the probe.
`train_stress_model` centres with a plain mean, so one missing entry poisons a whole
channel; dropping the two keeps 81 arrays where keeping them would keep 35. Losing the
proteasome channel is a real cost, and the `DTT → proteasome` pair scores 0.000 below.

### Readouts → modules

Module activity is the mean log2 ratio of the module's regulon, where the regulon is the
union of its transcription factors' MacIsaac 2006 targets present on the arrays, **with
every channel marker gene removed from every regulon** — not just from its own module's.
That is what stops a channel predicting a module because the module's truth contains that
channel.

| Module | Factors | Regulon genes | Measurable |
| --- | --- | ---: | --- |
| ESR | Msn2/Msn4 | 76 | yes |
| oxidative | Yap1/Skn7 | 119 | yes |
| hypoxia | Rox1/Upc2/Hap1 | 93 | yes |
| iron | Aft1/Aft2 | 72 | yes |
| proteasome | Rpn4 | 64 | yes |
| nitrogen | Gln3/Gat1 | 49 | yes |
| retrograde | Rtg1/Rtg3 | 31 | yes |
| sulfur | Met4/Met31/Met32 | 29 | yes |
| heat | Hsf1 | 28 | yes |
| cell_wall | Rlm1 | 22 | yes |
| carbon | Adr1/Cat8 | 15 | yes |
| osmotic | Sko1/Hot1 | 14 | yes |
| copper | Ace1 (Cup2)/Mac1 | 14 | yes |
| dna_damage | Rfx1/Crt1 | 7 | yes |
| UPR | Hac1 | 6 | yes |
| zinc | Zap1 | 6 | yes |
| xenobiotic | Pdr1/Pdr3 | 6 | yes |
| calcium | Crz1 | 1 | **no**, below the 5-gene floor |
| alkaline_ph | Rim101 | 2 | **no**, below the floor |
| redox, peroxide, atp, ph, nadh | — | 0 | **no**, no regulon exists |

The five pool modules are not a data gap. A metabolite pool has no transcription factor and
no regulon, so **no transcriptome can measure them, ever.** They are supplied as zero
columns, which `module_transfer`'s variance guard scores as 0.0, and they are exactly what
`roGFP2-Grx1`, `HyPer7`, `QUEEN-2m`, `pHluorin` and `Peredox` exist to read. Seven of 24
modules are outside this test's reach.

**Why MacIsaac rather than every binding record SGD holds.** The wider set runs to 500+
promiscuous ChIP hits per factor, and a regulon mean over hundreds of genes on a
median-normalised array *is* the array mean. Measured, not asserted: under the wide set the
`oxidative` module reads **+0.01 under H2O2 and −0.08 under menadione** and `dna_damage`
reads ~0.00 everywhere. Under MacIsaac they read +0.27 and +0.01. The tighter map is
strictly better as a positive control and no less independent.

---

## 5. Positive controls

A transfer number computed on a mapping that does not respond to its own stimulus is a
number about nothing, so this comes before the result. Full tables in
`outputs/external_channel_response.csv` and `outputs/external_module_response.csv`.

**Channels: 4 clean passes of 6 testable, 2 partial.** Only six of the 17 channels have a
matching stressor in this panel; the rest have no agent to respond to.

| Channel | Marker gene | Should peak under | Actually peaks under | |
| --- | --- | --- | --- | --- |
| `GATA-nitrogen` | DAL5/MEP2 | nitrogen depletion | **nitrogen depletion, +3.83** | pass, by a factor of four |
| `TRX2-oxidative` | TRX2 | H2O2 or menadione | **menadione, +2.60** (H2O2 +2.11) | pass |
| `STRE-osmotic` | GPD1 | sorbitol | **sorbitol, +2.06** | pass |
| `CSRE-carbon` | ADH2/ICL1 | glucose exhaustion | **glucose_starvation, +1.10** | pass |
| `HSE-heat` | HSP104 | heat | glucose_starvation +2.23, diamide +2.30, heat +1.84 | partial — diamide is a real HSP inducer |
| `UPRE-ER` | KAR2 | DTT | diamide +1.70, heat +0.80, DTT +0.66 | partial — see below |
| `STRE-general` | HSP12/CTT1 | everything | glucose_starvation +4.46, all eight ≥ +1.0 except DTT | the ESR, behaving as an ESR |

`Rbox-retrograde` (CIT2) has no matching agent — antimycin A is absent — and peaks at
+3.16 under glucose starvation, which is the known glucose derepression of CIT2. That is
not a scored control but it is not noise either.

**Modules: weaker, and one outright failure.**

| Module | Should peak under | Actually peaks under | |
| --- | --- | --- | --- |
| `nitrogen` | nitrogen depletion | **nitrogen depletion, +0.56** | pass |
| `carbon` | glucose exhaustion | **glucose_starvation, +0.43** | pass |
| `heat` | heat | diamide +1.15, heat +0.82 | partial |
| `UPR` | DTT | diamide +0.65, DTT +0.41 | partial |
| `ESR` | everything | diamide +0.82, all eight ≥ +0.12 | pass in shape |
| `oxidative` | H2O2 or menadione | **H2O2 +0.27**, but menadione +0.01 | weak: right argmax, wrong magnitude |
| `osmotic` | sorbitol | glucose_starvation +0.93, sorbitol +0.39 | **fail** |
| `proteasome` | DTT | diamide +0.54, DTT +0.05 | **fail** |

Three things to read out of this, all of which limit what follows.

**The regulon mean is a much weaker readout than the marker gene.** `TRX2` moves +2.60
under menadione while the whole Yap1/Skn7 regulon moves +0.01. Binding is not regulation,
and 119 bound promoters average away the ~10 genes that actually carry the oxidative
response.

**The panel's target lists are incomplete, and the data says where.** Diamide is the argmax
for `UPR`, `heat` and `proteasome`, and `STRESSORS["diamide"]` lists none of the three. Its
targets are `oxidative`, `redox` and `ESR`. A thiol oxidant that induces KAR2 1.7-fold and
HSP104 2.3-fold is doing more than the panel records. This is a proposed correction, not a
finding about the model — see §9.

**Zero-dose arrays are inside the blocks.** `dtt 000 min` and three `Heat Shock 000
minutes` arrays are the zero rungs of their ladders. Keeping them is right — they are the
bottom of the dose response — and it pulls the DTT and heat block means down, which is part
of why DTT loses its own channel to diamide.

---

## 6. Transfer, with nulls

Leave one stressor out, fit the latent basis and the readout on the rest, then score the
held-out stressor's **own target** modules — the modules `STRESSORS[s].targets` names,
restricted to the measurable 17 — with the fixed `module_transfer` metric and its
`_RECOVERABLE_SPAN` guard.

Latent width is **3**, chosen by `select_dimension` on the channels alone, which is blind
to the module labels. The width sweep is reported for transparency and deliberately not
used to choose: widths 1–6 give 0.562 / 0.676 / **0.688** / 0.674 / 0.689 / 0.665, so
choosing on it would move the headline by 0.001 and would be selecting on the test.
Channel variance is 0.509 / 0.157 / 0.089 / 0.055 / 0.051 in the first five components —
one dominant axis, which is the shape that makes the null below inevitable.

### Per held-out stressor, full eight-stressor panel

| held out | pairs | median | best | pairs > 0.25 |
| --- | ---: | ---: | ---: | ---: |
| heat | 2 | 0.806 | 0.930 | 2 |
| menadione | 2 | 0.796 | 0.872 | 2 |
| diamide | 2 | 0.788 | 0.937 | 2 |
| DTT | 3 | 0.717 | 0.740 | 2 |
| glucose_starvation | 2 | 0.485 | 0.746 | 1 |
| sorbitol | 2 | 0.346 | 0.692 | 1 |
| H2O2 | 2 | 0.341 | 0.590 | 1 |
| rapamycin | 3 | **0.000** | 0.810 | 1 |

**12 of 18 pairs clear the 0.25 recovery threshold; median +0.688.** Every pair, with its
peer counts and the span of its held-out truth, is in
`outputs/external_transfer_pairs.csv`:

| held out | module | score | peers (panel) | peers (measured) | truth SD |
| --- | --- | ---: | ---: | ---: | ---: |
| diamide | ESR | 0.937 | 7 | 5 | 0.227 |
| heat | ESR | 0.930 | 7 | 5 | 0.907 |
| menadione | ESR | 0.872 | 7 | 5 | 0.182 |
| rapamycin | ESR | 0.810 | 7 | 6 | 0.283 |
| glucose_starvation | ESR | 0.746 | 7 | 5 | 0.621 |
| DTT | ESR | 0.740 | 7 | 6 | 0.269 |
| menadione | oxidative | 0.719 | 2 | 1 | 0.085 |
| DTT | UPR | 0.717 | 0 | 3 | 0.495 |
| sorbitol | ESR | 0.692 | 7 | 5 | 0.231 |
| heat | heat | 0.683 | 0 | 5 | 1.598 |
| diamide | oxidative | 0.640 | 2 | 1 | 0.054 |
| H2O2 | ESR | 0.590 | 7 | 5 | 0.254 |
| glucose_starvation | carbon | 0.223 | 1 | 1 | 0.463 |
| H2O2 | oxidative | 0.092 | 2 | 0 | 0.102 |
| DTT | proteasome | 0.000 | 0 | 2 | 0.129 |
| rapamycin | nitrogen | 0.000 | 0 | 0 | 0.096 |
| rapamycin | carbon | 0.000 | 1 | 1 | 0.106 |
| sorbitol | osmotic | 0.000 | 0 | 4 | 0.122 |

### And is that better than no structure at all?

| | median own-target score | null median | p | verdict |
| --- | ---: | ---: | ---: | --- |
| **observed** | **+0.6876** | — | — | — |
| null: random axes, same rank (`rotated_subspace`) | +0.6876 | **+0.6953** | **0.990** | **DOES NOT BEAT** |
| null: no cross-channel structure (`matched_marginals`) | +0.6876 | +0.0003 | 0.005 | beats |

Skill over the rotated null: **−0.025**. Negative, so arbitrary axes carry marginally more
transfer information than the fitted ones.

Restricting to the six stressors with nothing traceable to Gasch — dropping diamide and
sorbitol, 68 arrays, width 2 — moves the number and not the verdict:

| | median | null median | p | verdict |
| --- | ---: | ---: | ---: | --- |
| observed, six-stressor panel | +0.6048 | — | — | — |
| `rotated_subspace` | +0.6048 | +0.6077 | 0.726 | **DOES NOT BEAT** |
| `matched_marginals` | +0.6048 | +0.0000 | 0.005 | beats |

Skill **−0.007**. 10 of 14 pairs clear the threshold.

**What this means.** Clearing `matched_marginals` says only that these channels are
correlated at all, which any low-rank projection of correlated data exploits. Failing
`rotated_subspace` says the *identity* of the axes contributes nothing. One component
carries 51% of the channel variance, and that component is the general stress response —
which every one of these eight stressors drives, and which every regulon mean partly
tracks. A random rotation of 17 correlated channels retains that axis, so it predicts the
regulon means just as well. **The scores are real and the latent structure is not what
produces them.**

This is the same conclusion the simulation-only null reached, arrived at independently, on
data the model never saw, with the sign of the gap unchanged and the p-value worse.

---

## 7. The redundancy law

The claim under test, from
[`superseded/module-transfer-inflation.md`](superseded/module-transfer-inflation.md): across
all 76 simulated (held-out stressor, own target) pairs, recovery tracks how many *other*
stressors drive the same module — 0 peers recovered 0 of 6, 6+ peers recovered 84%,
Spearman +0.591.

Real data, peers counted from the panel's target lists over the 8-stressor panel:

| peers | pairs | mean score | recovered |
| ---: | ---: | ---: | ---: |
| 0 | 5 | 0.280 | 40% |
| 1–2 | 5 | 0.335 | 40% |
| 3–5 | 0 | — | — |
| 6+ | 8 | 0.790 | **100%** |

Spearman **+0.734**, p = 0.001, n = 18. Against a permutation null on the peer counts:
observed +0.734, null median +0.027, **p = 0.005**.

The panel's own bookkeeping is not the best peer definition available here, because §5
showed the target lists are incomplete. Counting peers from the data instead — a stressor
drives a module when its regulon mean moves more than 0.25 in log2, computed only over the
*other* stressors, all of which are in training — fills in the middle of the range and
tightens the bottom:

| peers, measured | pairs | mean score | recovered |
| ---: | ---: | ---: | ---: |
| 0 | 2 | **0.046** | **0%** |
| 1–2 | 5 | 0.316 | 40% |
| 3–5 | 9 | 0.685 | 89% |
| 6+ | 2 | 0.775 | 100% |

Spearman **+0.680**, p = 0.002, against a permutation null of +0.005, p = 0.005.

**The zero-peer edge holds, and the two apparent exceptions confirm the mechanism rather
than breaking it.** `DTT → UPR` (0.717) and `heat → heat` (0.683) both score well with zero
peers *by the panel's count*. By the measured count they have 3 and 5 peers: diamide raises
KAR2 1.7-fold and the Hsf1 regulon +1.15, and the panel's `diamide` entry records neither.
The pairs that genuinely have no peer under either definition — `rapamycin → nitrogen`
(0.000) and `H2O2 → oxidative` (0.092) — recover nothing. Two of two, mean 0.046, in the
same direction as simulation's 0 of 6.

**And then the interpretation collapses.** A permutation null says the relation is not
chance. It does not say the relation is evidence about the *fitted* basis, because the
scores it is built from do not beat arbitrary axes. Running the law itself under rotated
readings:

| panel | observed Spearman | rotated null median | p | verdict |
| --- | ---: | ---: | ---: | --- |
| eight stressors | +0.734 | +0.664 | 0.030 | beats, narrowly |
| six stressors (independent) | +0.786 | +0.786 | 0.871 | **DOES NOT BEAT** |

On the panel with nothing traceable to Gasch, arbitrary axes reproduce the law exactly. So
the honest statement is: **redundant coverage predicts what is recoverable from these
channels, and that is a fact about the design of the panel rather than about the latent
state.** Which is, read carefully, what the retraction already said — "interpolation across
redundant coverage, not generalisation to unseen biology." Real data agrees, and adds that
you do not need the fitted basis to get it.

**One caveat that limits this whole section.** On the eight-stressor panel every pair in the
6+ bin is ESR, and ESR is the one module whose panel encoding traces to Gasch 2000. The
high-peer end of the relation is a single module measured eight times. The measured-peer
definition spreads the range across `UPR`, `heat`, `osmotic` and `proteasome` as well, which
is why both definitions are reported.

---

## 8. A null the repository has never run

[`NULL_RESULTS.md`](NULL_RESULTS.md) §"Still without a null" lists *"These channels identify
these modules"* as needing `shuffled_loadings`. It can be run here, on real readings, and
it is the one thing in this file that passes.

The known-loadings route — Path A in [`CLAIM_BOUNDARY.md`](CLAIM_BOUNDARY.md), no fitting at
all — estimates module activity as `x = y (L⁺)′` with `L = reporter_loadings(channels)`, 17
channels × 24 modules, rank 17. Scored as the median Spearman correlation across the 17
measurable modules between that estimate and the regulon truth:

| | median Spearman | null median | p | verdict |
| --- | ---: | ---: | ---: | --- |
| eight-stressor panel | **0.465** | 0.091 | 0.005 | **beats** |
| six-stressor panel | **0.503** | 0.079 | 0.005 | **beats** |

`shuffled_loadings` permutes every entry of `L`, keeping the multiset of loading strengths
and destroying which reporter reads which module. The literature matrix beats it by a
factor of five. **The crosstalk topology in `reporter_loadings` carries real information
about which channel reads which regulon**, measured against data none of it came from. It
is still Tier 0 — transcribed, not measured — but it is transcribed correctly.

Note the asymmetry with §6, which is the most useful thing in this file. The *given* loading
matrix beats its null; the *fitted* latent basis does not beat its own. That is the
`IDENTIFIABILITY.md` §1–2 argument turning up as a measurement: Path A works because `L` is
known, and Path B — which every transfer number in the README goes through — is trying to
estimate what Path A is handed.

---

## 9. Limitations

**The module truth is a proxy, and a weak one.** A regulon mean over ChIP-derived binding
targets is not module activity. Binding is not regulation; `TRX2` moves +2.60 under
menadione while its whole regulon moves +0.01. Three modules fail their positive control
outright (`osmotic`, `proteasome`, and `oxidative` on magnitude). Everything in §6 and §7
inherits that.

**Channels and truth are the same measurement.** Both are log2 ratios from the same array,
so array-level artefacts — dye bias, normalisation, the print run — inflate the
channel-to-truth relation for every module alike. This is precisely what `rotated_subspace`
controls for, which is why that null is the load-bearing one and why failing it is the
result rather than a technicality.

**Eight stressors is small, and the peer range is thin.** 18 own-target pairs against
simulation's 76. The 3–5 peer bin is empty under the panel definition. Nothing here has the
power to separate a Spearman of +0.59 from +0.73.

**Seven modules of 24 are untestable.** Five pools have no regulon at all; `calcium` (1
gene) and `alkaline_ph` (2) fall below the 5-gene floor. Nine of the 17 measurable modules have no
driving stressor in this panel -- `cell_wall`, `copper`, `dna_damage`, `hypoxia`, `iron`,
`retrograde`, `sulfur`, `xenobiotic`, `zinc` -- so they contribute to the fit and never to
a score.

**Two channels are missing.** No `Xbox-dna`, no `PACE-proteasome`, both because whole print
runs omit the probe. `DTT → proteasome` scoring 0.000 is partly the absence of the channel
that reads it.

**One mapping is by module profile, not by agent.** Nitrogen depletion stands in for
`rapamycin`. Its positive controls are the strongest in §5, but it is a substitution.

**Zero-dose arrays are inside the treatment blocks**, which lowers the DTT and heat block
means and is why DTT loses its own channel to diamide.

**No dose is fitted anywhere.** `doses` is supplied as zeros. Gasch's timecourses are a
time axis, not a concentration ladder, and the panel's dose-response machinery is not
exercised by any of this.

**This is not Tier 3.** Nothing was registered before the measurement. Under
[`CLAIM_BOUNDARY.md`](CLAIM_BOUNDARY.md) the results here are the second half of Tier 2 —
agreement, or in this case disagreement, with a measurement the model never saw — and they
belong in the same row as the Vos 2016 retentostat check and the van Hoek 1998 failure.

### Two proposed changes, neither made here

Both are for the owner, not for this file, and neither touches a result above.

1. **`STRESSORS["diamide"]` is missing three arms.** It lists `oxidative`, `redox` and
   `ESR`. On these arrays diamide is the argmax for `UPR` (+0.65), `heat` (+1.15) and
   `proteasome` (+0.54), and raises KAR2 1.7-fold and HSP104 2.3-fold. A thiol oxidant
   causing ER and cytosolic misfolding is mechanistically unsurprising. Adding `UPR`,
   `heat` and `proteasome` arms would change peer counts across the panel, so it should be
   done deliberately and with its own sources — Gasch 2000 alone is the wrong warrant given
   §3.
2. **`train_stress_model` cannot take a channel with a missing entry.** It centres with
   `readings.mean(axis=0)`, so one NaN poisons a whole channel, while `fit_latent`
   underneath it masks missingness properly and documents that it does. Switching that one
   call to `np.nanmean` would have kept `Xbox-dna` and `PACE-proteasome` and taken this
   analysis from 81 arrays to 96 on 19 channels. It is a one-word change inside `src/`,
   which this work is not permitted to make.

---

## Reproduce

```bash
# 23 MB, no credentials: GSE18 from GEO, plus one JSON per transcription factor from SGD.
python3 scripts/fetch_external_validation_data.py

# The full eight-stressor panel. ~6 min at 200 null draws.
python3 scripts/run_external_validation.py --draws 200

# The circularity-restricted panel: diamide and sorbitol dropped entirely.
python3 scripts/run_external_validation.py --independent --draws 200
```

Writes, with `_independent` appended for the restricted run:

| File | Contents |
| --- | --- |
| `outputs/external_conditions.csv` | every one of the 155 arrays, its title, its assigned stressor, and why it was excluded if it was |
| `outputs/external_module_map.csv` | each module, its factors, its regulon size, whether it is measurable |
| `outputs/external_channel_response.csv` | mean log2 per channel per stressor, and the argmax — the §5 positive control |
| `outputs/external_module_response.csv` | the same for module regulon means |
| `outputs/external_transfer_pairs.csv` | all 18 (held-out stressor, own target) pairs with both peer counts and the truth span |
| `outputs/external_redundancy.csv` | the peer-count bins under both definitions |
| `outputs/external_nulls.csv` | every null in this file: statistic, null median, p, draws, verdict |

Both runs are deterministic: `compare_to_null(seed=0)`, `select_dimension` and
`module_transfer` are all seeded, and the only stochastic input is the surrogate draws.

## Provenance

Data fetched and analysis run 2026-08-26. All four identifiers verified against PubMed
esummary on the same day, titles quoted:

- **Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown
  PO.** PMID 11102521 — *"Genomic expression programs in the response of yeast cells to
  environmental changes."* Mol Biol Cell 11:4241-57, 2000. GEO GSE18, public domain.
- **MacIsaac KD, Wang T, Gordon DB, Gifford DK, Stormo GD, Fraenkel E.** PMID 16522208 —
  *"An improved map of conserved regulatory sites for Saccharomyces cerevisiae."* BMC
  Bioinformatics 7:113, 2006. The module labels, via SGD's curated regulation records.
- **Harbison CT, Gordon DB, Lee TI, et al.** PMID 15343339 — *"Transcriptional regulatory
  code of a eukaryotic genome."* Nature 431:99-104, 2004. Attempted as the label source and
  rejected: the deposit publishes its matrix without labels. §2.
- **Martinez-Pastor MT, Marchler G, Schüller C, Marchler-Bauer A, Ruis H, Estruch F.** PMID
  8641288 — the STRE element behind `MODULES["ESR"]`, and the half of that entry which is
  not Gasch. §3.

Regulation records were fetched from `www.yeastgenome.org/backend/locus/{gene}/
regulation_details` on 2026-08-26, one request per factor, and are cached under
`data/external/sgd/regulation/`. SGD's records carry an evidence type and a reference PMID
per row; `load_regulons` asserts that PMID 11102521 supplies none of them.
