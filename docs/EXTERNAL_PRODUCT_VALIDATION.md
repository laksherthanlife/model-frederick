# Product comparisons and their validation scope

The deliverable this file exists for, stated once: *based on this paper, under these
conditions, they measured X; the model predicts Y.*

These comparisons do not all have the same holdout scope. Elizondo is the calibration
cohort, with historical fixed-candidate strain holdouts in §1; Torello Pianale is an
external condition comparison in §2; Kocharin supplies source-cohort descriptions and
retrospective state holdouts in §3; **§4 is neither — it is a pair of within-paper matched
contrasts, which is a weaker design than a holdout and a stronger one than a between-paper
comparison, and it is named that way throughout.** The former claim that nothing here entered
a fit or model choice was wrong. Corrected candidates await governed artifact regeneration; a
calculation inspected outside `outputs/` is not an adopted result or new biological data.

---

## 1. Elizondo 2025 — three beta-carotene strains, two dilution rates

PMID 40891387, doi 10.1021/acssynbio.5c00256. Six chemostat steady states, three strains.

**Historical fixed-candidate strain holdout.** The paper is not withheld — these strains
are the calibration set. A **whole strain** is withheld while refitting the expression-to-flux
scalar and branch kinetics on the other two. This excludes the held-out product values
from those fits, but CrtE and the model family were chosen after examining this cohort.
The table below is therefore a selection-unadjusted comparison, not evidence that every
model choice was independent of its test outcomes. The corrected nested workflow selects
within outer-training strains; do not substitute the better historical fixed-CrtE numbers
for that selected-pipeline evaluation.

**What goes in:** one relative CrtE mRNA value, and the dilution rate. That is the entire
input.

| state | strain | CrtE | q_p measured | q_p predicted | fold | lycopene measured | lycopene predicted | fold |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D01 | b-car2 | 0.271 | 1.513e-4 | 1.612e-4 | 1.07× | 9.345e-4 | 1.126e-3 | 1.21× |
| 2D025 | b-car2 | 0.245 | 1.842e-4 | 1.915e-4 | 1.04× | 3.314e-4 | 2.251e-4 | 0.68× |
| 3D01 | b-car3 | 0.615 | 2.271e-4 | 2.010e-4 | 0.89× | 4.674e-3 | 4.426e-3 | 0.95× |
| 3D025 | b-car3 | 0.562 | 2.948e-4 | 3.509e-4 | 1.19× | 4.974e-4 | 9.514e-4 | 1.91× |
| 4D01 | b-car4 | 1.000 | 1.853e-4 | 2.464e-4 | 1.33× | 7.533e-3 | 7.056e-3 | 0.94× |
| 4D025 | b-car4 | 0.644 | 4.536e-4 | 3.766e-4 | 0.83× | 1.545e-3 | 9.480e-4 | 0.61× |

**Product: median 14.2% error, worst 1.33×. Intermediate: 26.3%, worst 1.91×.**

Rates are mmol/gDCW/h; lycopene is a content in mmol/gDCW. This historical table is not the
current nested headline from `scripts/predict_product.py`; retain its fixed-candidate scope
when comparing it with that runner's outputs. Follow the governed recipe in
[`REPRODUCING.md`](REPRODUCING.md), not an in-place rerun to overwrite the old evidence.

**What this does and does not establish.** The historical fixed-candidate calculation
predicts from cassette expression rather than the held-out strain's product measurement.
That addresses direct target leakage for that calculation, not model-selection bias or
external biological validation. It does **not** establish transfer to another host, promoter
set or pathway: three strains from one paper are one source cohort, not a survey. The
corrected nested companion is an evaluation of selection within that same cohort, not
additional independent biology.

---

## 2. Torello Pianale 2022 — a growth rate, and a condition label

doi 10.3389/fmicb.2021.802169, Front Microbiol 12:802169. Open access. Nothing here was
fitted to it.

**Measured (their Table 4):** CEN.PK113-7D parental, **0.379 ± 0.000 /h**, in Delft
(Verduyn) mineral medium, 20 g/L glucose, 30 °C, 150 mL in a 500 mL flask at 140 rpm, with
N₂ flushed for ten seconds. The paper's own word for the condition is **microaerobic**.

| the model asked as | predicts | vs 0.379 |
| --- | ---: | ---: |
| aerobic, air | 0.400 | **+5.5%** |
| strictly anaerobic, supplemented | 0.310 | −18.2% |
| microaerobic, 5% O₂ | 0.309 | −18.6% |
| strictly anaerobic, **not** supplemented | **0.000** | — |

A shaken 500 mL flask holding 150 mL, flushed once and then left, is an aerated culture: the
headspace is 350 mL and 140 rpm keeps transferring oxygen. Asked as air, the model is 5.5%
high, in the direction it should be.

**The last row is the result, not the first.** The paper reports no Tween-80 and no
ergosterol. This model holds that an unsupplemented anaerobic glucose culture does not grow
**at all**, because yeast cannot make a sterol without molecular oxygen (Andreasen & Stier
1953, PMID 13034889). They measured 0.379 /h. Those are compatible only if the culture had
oxygen — so the model reads the condition correctly, and independently agrees with the
authors' own label. A second route says the same thing: 0.379 is above
`ANAEROBIC_MU_MAX_PER_H` (Verduyn 1990's measured 0.31), so it is unreachable anaerobically
under this model whatever the supplements.

**A near miss worth recording.** The first pass at this used a summary of the paper that
called the cultures "anaerobic". Read that way the model is 18% low and the constant looks
refuted. The methods say microaerobic. The difference between the two readings is reading
the methods, and it is the difference between a refuted constant and a corroborated one.

Pinned by `tests/test_external_torello_pianale.py`.

### What their sensor design says about ours

Their five sensors split along exactly the axis of this repository's headline finding:

- **Ratiometric, one fluorophore, two excitations** — QUEEN-2m (410/480, em 520) and
  sfpHluorin (390/470, em 512). Both channels report the *same molecule*, so anything
  scaling that molecule — growth dilution above all — cancels exactly in the ratio. These
  are **structurally immune** to the confound.
- **Intensiometric against a constitutive reference** — GlyRNA, OxPro, RibPro, normalised to
  mCherry from pTEFmut8. Dilution cancels only insofar as two different proteins share
  maturation and degradation. Partial protection.
- **Single channel, no reference** — this repository's UPRE1, UPRE2, NativeYap1 and
  AlteredYap1, all single mCitrine. Fully exposed.

That ordering explains something the paper does not comment on: it reports **no analysis of
growth-rate dependence at all**, and does not address signal dilution quantitatively. Its
ratiometric sensors did not need it to. Ours did, which is why the correction had to be done
in software and why it removed most of the apparent response.

**The actionable form:** the fix for this project's largest measurement problem is a
construct change, not an analysis change. A ratiometric sensor gets for free what
`reporter.py::promoter_activity` has to model.

---

---

## 3. Kocharin & Nielsen 2013 — PHB environment contrasts within one source cohort

PMID 23514405, PMC3610212, CC-BY. Eleven aerobic chemostat steady states, **one genotype**
(SCKK006), four dilution rates crossed with three carbon feeds matched at 0.666 Cmol/L,
content in mg/gDW. Every one of the eleven rows was checked against the Europe PMC full
text before use; that source check is not independent biological replication.

Holding genotype fixed while varying environment exposes a limitation of a μ-only
prediction: it cannot express the reported between-feed contrast at matched μ. Fitting
new environment terms on these summaries is retrospective development within this cohort.

**Reported flux spans 4.24× across feeds at identical dilution rate.**

| q_PHB (mmol/gDCW/h) | glucose | ethanol | glucose:ethanol 1:2 |
| ---: | ---: | ---: | ---: |
| D = 0.05 | 0.00251 | 0.00961 | 0.01065 |
| D = 0.10 | 0.00649 | 0.01567 | 0.01442 |
| D = 0.15 | 0.01035 | 0.00923 | 0.01307 |
| D = 0.20 | 0.01299 | *washes out* | 0.01059 |

A μ-only model returns the same flux across feeds at a fixed μ and cannot express that
observed contrast. This challenges that restricted model, not every possible environment
law. The separately fitted PHB content-response terms below are retrospective additions,
not independent validation of the environment channel.

**Descriptive growth-rate slopes differ by feed.** Glucose gives q ~ μ^1.20, ethanol μ^0.04
and the mix μ^0.03 — a spread of 1.18; two point estimates lie outside the [0.53, 1.78]
interval fitted to the β-carotene exponent. These few same-cohort estimates are not
independent biological tests of distinct exponents or proof that no common law can work.

**The GEM does not rescue it, and the measured flux is what has to be explained.** q_PHB
moves <!-- audit:value table=outputs/phb_environment_score.csv column=q_span_over_all_states row="state_id=glc_D005" -->6.2309×
across all eleven states, on one genotype.

*This paragraph read "a ceiling that moves **1.3×** across all eleven states while the
measured flux moves **5.2×**. Pearson against the measurement is **−0.016**" until 2026-09-04.
Its source, `scripts/score_phb_environment.py`, imports no cobra and computed none of the
three. Two of them came from a CEILING formulation that `pathway/gem_environment.py` has since
measured to be the wrong question and abandoned — the chain's predicted flux sits 12.8–195×
below the ceiling so it never binds — so they were deleted rather than recomputed under a
retracted method. The third was computable and wrong: 5.2× is the span with `etoh_D010` and
`mix_D010` removed, a partial state set labelled as the whole. It is now computed, written to
`outputs/phb_environment_score.csv`, and marked above.*

The verdict is the same wall as everywhere else: FBA returns **ceilings**, the biology runs
far below them, and a ceiling therefore carries no information about the flux. It is why
`fba/audit.py` audits rather than supplies, and the PHB data says the same thing about
precursor supply that FVA said about the product reaction. What replaced the ceiling question
is a YIELD question, scored on these same eleven states — see `pathway/gem_environment.py`.

### Flux candidates: descriptive scores, inference pending

Carbon source removes 28.7% of the variance in log(flux) **in sample**. That description
of this cohort does not validate a predictive law. The former stronger claim that nothing
in the data predicts flux, and that eleven states cannot identify an environment effect,
was not established by the calculation cited for it.

**Historical table, retained to identify the withdrawn comparison:**

| model | parameters | old LOO skill | old alleged 95th-percentile bar | old verdict (withdrawn) |
| --- | ---: | ---: | ---: | --- |
| constant only | 1 | −0.100 | −0.100 | — |
| log(μ) | 2 | −0.153 | −0.046 | inside the null |
| carbon source | 3 | −0.140 | −0.046 | inside the null |
| carbon source + log(μ) | 4 | −0.162 | −0.056 | inside the null |

There were two different errors. The comparator used the **full-data mean**, including
the response being held out, making even an intercept-only fit score negative when it
should match a training-only mean comparator. And `scripts/flux_null_baseline.py` redrew
Gaussian covariates on each response shuffle: it did **not** refit each real candidate
under a justified null. The old −0.046 95th-percentile and +0.077 (roughly 0.08)
99th-percentile bars and their claimed identifiability/sample-size consequence are
withdrawn, not replaced by new random-design percentiles.

The corrected script separates two artifacts. `flux_null_baseline.csv` retains all five
parameter-count rows as `random_design_demonstration_not_inference`, 2000 draws per size,
seed 0. `flux_candidate_scores.csv` retains the four fixed candidates. Both model and
comparator are fitted on the ten training states in each state-LOO fold. The score is
`1 - RMSE(model) / RMSE(fold-training-mean log flux)`, not R².

**Corrected external candidate inspected during integration; not yet an adopted root
artifact.** Its common comparator log-flux RMSE is 0.5313503:

| model | parameters | descriptive state-LOO RMSE skill | score status | inference status |
| --- | ---: | ---: | --- | --- |
| constant only | 1 | 0 (roundoff only) | ok | pending_exchangeability |
| log(μ) | 2 | −0.0485 | ok | pending_exchangeability |
| carbon source | 3 | −0.0364 | ok | pending_exchangeability |
| carbon source + log(μ) | 4 | −0.0561 | ok | pending_exchangeability |

All three nonconstant candidates remain worse than the training-mean comparator on this
cohort. All four scores are defined; their empty failure reasons mean no failed fits, not
missing candidates. All eleven condition summaries remain: four glucose, three ethanol
and four mixed-feed states across four dilution rates. The table supplies no cultivation,
biological-replicate or randomization identifiers. Neither row count nor feed labels
establishes independent or exchangeable units, so **no candidate p-value, null percentile,
significance verdict or selection-adjusted claim is issued**. Rank-deficient or incomplete
fits refuse the whole score; failed permutations retain their positions, and a single
failed demonstration draw withholds that size's percentile rather than filtering it out.
The fixed-design permutation helper requires explicit response permutations; it does not
certify their exchangeability. `scripts/environment_channel.py` instead relabels a carbon
descriptor across feeds for **log-content** candidates. Its six relabellings do not justify
response shuffles or significance claims for these log-flux candidates.

The source-cohort flux descriptions come from `scripts/score_phb_environment.py` and
`tests/test_phb_environment.py`; comparator and failure semantics are checked by
`tests/test_flux_null_baseline.py`. Reproduction and adoption are separate operations under
[`REPRODUCING.md`](REPRODUCING.md).

### The content comparison is a different target, not a rescue of the flux null

`scripts/env_to_product.py` retains all eleven PHB measurements and scores log-content
regressions leave-one-state-out: μ only versus μ plus ethanol-carbon fraction and their
interaction. Median fold error is 1.505× versus 1.320×, but the environment fit is more than
2× wrong on **two** states (`glc_D005`, `etoh_D005`) rather than the μ-only fit's one
(`glc_D005`); its worst error is about 2.55×. The shipped environment factor is reported
separately, not presented as a held-out refit. These are retrospective descriptions of the
same Kocharin cohort, not an independent paper or new biological replication.

The companion β-carotene comparison retains **six** fixed-gene requests, glucose/ethanol
at μ = 0.101, 0.15 and 0.2543 /h: **four answered, two refused**. Ethanol at the two higher
setpoints raises `SetpointUnreachable` under the context maximum of 0.140 /h. Those rows
retain their request keys and reasons, with missing returned growth/content and
`environment_sets_the_number`, and `environment_layer_state=not-run`. Answered rows have
finite μ-only content, `reported` layers and `environment_sets_the_number=False`. Only
0.101 /h has a complete answered carbon pair. Its equal model contents do not validate
biological invariance; at the other rates that comparison is unavailable. The old six-answer
artifact and its all-`reported`/`nunique` interpretation must not be used as current evidence.

### The companion result, which is worse

Kocharin et al. **2012** (PMID 23009357, PMC3519744) supplies the other half. All four
strains carry the same `phaA` on the same PGK1 promoter — **entry-enzyme expression is 1.0
across the panel by construction** — and PHB specific productivity still moves **16.5-fold**,
on host acetyl-CoA supply alone.

So for this pathway the flux law's only input does not vary and the flux does. The
β-carotene result — where the genotype *did* vary and expression tracked it at skill +0.63 —
does not transfer here unexamined, and `data/pathways/phb.toml` carries no calibration entry
for exactly that reason.

Two further independent papers point the same way and are worth reading before anyone fits
an mRNA-reading law again: Hazelwood 2009 (PMID 19734328) concludes that transcriptional
regulation "could not quantitatively account for" storage-carbohydrate changes and implicates
post-transcriptional control; Boender 2009 (PMID 19592533) has glycogen roughly doubling
while GSY1 and GSY2 transcripts do not move.

---

## 4. Verwaal 2007 and Ukibe 2009 — two within-paper matched contrasts, 2026-09-09

PMIDs 17496128 and 19801484, both open access (PMC1932764, PMC2786542), both read in full
on 2026-09-09. `scripts/matched_contrast.py`, pinned by `tests/test_matched_contrast.py`.

**The honest name first, because the name is the entire question.** A within-paper
before/after is a **matched contrast**, not a held-out prediction. Nothing was withheld and
nothing was forecast. The model issues **no number for either arm of either pair**:
`predict.py` takes a relative entry expression, `data/carotenoid_batch/published_batch_titres.tsv`
states an expression for no row, and this package holds no genotype-to-content law. Any
sentence here of the form "the model predicts X for this strain" would be false.

**What a matched contrast can still test.** `pathway/solve.py` asserts one thing that needs
no genotype input at all — for a terminal node fed by a saturating step, `content <=
vmax_per_growth` with `mu` cancelling exactly, for **every** genotype at **every** growth
rate. That is the one model output a table with no expression column can falsify, and it is
the whole reason these two pairs are worth running. The refutation already on the books,
Arhar 2024 at 79 mg/gDCW, is a **between-paper** comparison, and every between-paper
comparison carries the confounds `pathway/capacity.py` names: host, promoter, medium,
cultivation, assay and biomass basis all move together with the cassette. A within-paper
contrast removes all of them at once. It is the design `capacity.py`'s refusal says the
literature lacks — and the literature turns out to hold two of them.

### 4.1 The one-sided ceiling bracket, over all nine `dfba_ready` rows

Computed from the shipped calibration, not from a literal: `content_ceiling(load_pathway(
"beta_carotene"), BETA_CAROTENE_KINETICS).mg_per_gdcw` = **1.2483 mg/gDCW**.

| strain | mg/gDCW | × ceiling | position | measurand |
| --- | ---: | ---: | :---: | --- |
| Ukibe 2009 | 0.018 | 0.014 | inside | β-carotene |
| Li 2013 | 0.390 | 0.312 | inside | β-carotene |
| Bubphasawan 2025 | 2.150 | 1.722 | **outside** | β-carotene |
| Lange 2011 | 3.897 | 3.122 | **outside** | β-carotene |
| Verwaal 2007 | 5.900 | 4.726 | **outside** | β-carotene |
| Xie 2014 | 7.410 | 5.936 | **outside** | β-carotene (of `both_reported`) |
| Bu 2022 | 11.400 | 9.132 | **outside** | β-carotene |
| Fathi 2021 | 46.500 | 37.250 | **outside** | β-carotene |
| Arhar 2024 | 79.000 | 63.285 | **outside** | β-carotene |

**2 inside, 7 outside, worst 63.285×.** The prior read-only count is confirmed exactly. Every
one of the nine is a specific assay — eight `beta_carotene`, one `both_reported` scored on its
HPLC β-carotene number and not on its 453 nm total — so none of these is the López 2019 trap
that §5 of `EXTERNAL_CAROTENOID_BOUND.md` records.

**"Inside" is not evidence for the bound.** The bracket is one-sided by construction: a
content above an upper bound refutes it, and a content below is compatible with the bound and
with every smaller value. Both inside rows sit two orders of magnitude below it. Nothing here
supports the ceiling; two rows merely fail to attack it.

A detail worth writing down once so nobody conflates them: Ukibe's 0.018 and Li 2013's 0.390
are two different papers, and **0.390 mg/gDCW is also the value of Ukibe's own +BTS1 arm.**
Substituting that arm leaves the bracket unchanged at 2 inside / 7 outside, because 0.390 is
still inside.

### 4.2 Verwaal 2007 — one crtI copy, and the bound is crossed inside one paper

Table 4 of Appl Environ Microbiol 73:4342, both values from one table, one HPLC diode-array
method, one balance:

| arm | strain | β-carotene | × ceiling |
| --- | --- | ---: | ---: |
| control | CEN.PK113-6B + YIplac211 YB/I/E\* + YIplac204 *tHMG1* | 501 µg/gDCW | 0.401 |
| derivative | the same, **+ YIplac128 *crtI*** | 5,918 µg/gDCW | **4.741** |

**11.81×**, on one added *crtI* copy. Held constant across the pair: host CEN.PK113-6B, TDH3p
on every *crt* gene, YNB + 2% glucose, 72 h at 30 °C and 225 rpm, one HPLC method, one table.
"Single-copy integration of the constructs was confirmed by Southern blotting by standard
laboratory procedures", so the dosage step is measured rather than inferred from a construct.

**The contrast straddles the bound.** The control is inside at 0.401×; the derivative is
4.741× outside. There is no appeal available to a different assay class, a different biomass
basis, a different host, a different medium or storage engineering, because both numbers come
from the same strain background measured the same way on the same day.

**And here the model's answer is not "unsupported" — it is a number, and the number is
wrong.** `data/pathways/beta_carotene.toml` declares the phytoene node — the node CrtI drains
— as `passthrough`: flux in equals flux out, no capacity, no fitted parameter.
`BETA_CAROTENE_KINETICS` holds an entry for the **lycopene node alone**. So the shipped model
predicts a fold of **exactly 1.000 at any crtI dosage**, and Verwaal measured 11.81. Three
independent routes into the answer are all closed for this gene: CrtI is not `spec.entry_enzyme`
(that is CrtE), it is not `CAPACITY_SETTING_ENZYME["beta_carotene"]` (that is crtYB), and its
node carries no rate law that can vary. `scripts/matched_contrast.py::model_response` derives
that from the spec rather than looking it up, so a spec edit moves the verdict.

**The sting is in the citation.** `beta_carotene.toml`'s own header reads "Chemistry from
Verwaal et al. 2007 … which built this pathway in *S. cerevisiae*", and its notes already
concede that the phytoene passthrough is "A LUMPED APPROXIMATION, NOT EVIDENCE OF AN ABSENT
POOL", citing Chen 2016 for limiting CrtI conversion **in another strain**. The paper the file
cites for its chemistry contains the matched contrast that settles the same question in the
lineage the file is about. It did not need another strain.

### 4.3 Ukibe 2009 — 21.7×, and it tests nothing here

| arm | strain | β-carotene | × ceiling |
| --- | --- | ---: | ---: |
| control | INVSc1 + pTV-*crtI* + pUV-*crtYB* | ~18 µg/gDCW | 0.014 |
| derivative | the same, **+ pHV-*BTS1*** (GGPP synthase) | ~390 µg/gDCW | 0.312 |

**Both arms are inside the ceiling**, so this pair cannot refute a one-sided upper bound and
is recorded as `not_a_ceiling_test`. Counting it as a second refutation would be counting a
null as a result, and `tests/test_matched_contrast.py` asserts that exactly one of the two
contrasts refutes.

It is reported for two other reasons. First, it is the **reproduction check**: the paper's own
words are "a prominent **22-fold** increase", and the two printed contents give 21.67× — 1.5%
apart. Verwaal prints no fold at all for its step, so there `paper_stated_fold` is `None` and
nothing is reproduced against an invented target; its arms are instead cross-checked against
the vendored survey. Second, BTS1 is the *S. cerevisiae* GGPP synthase, so it acts **upstream
of the chain's first node**, on the precursor `s_0189`. That is the one route the model does
expose — the entry flux — which makes Ukibe the pair the model could in principle express and
cannot, for want of an expression number the table does not carry. It is a null for the
ceiling and an unmet input requirement for the flux law.

### 4.4 What neither contrast can establish

Assembled from the record rather than written beside it, so a better-reported pair stops
carrying the caveat:

- **Neither is a held-out prediction.** Nothing was withheld; the model forecasts neither arm.
- **Neither carries a dispersion.** Verwaal's Table 4 is "averages of two independent
  cultures, except for YB/I/E+tHMG1+I, which are triplicates", with no standard deviations
  printed; Ukibe's Figure 2b shows no error bars and states no replicate count. **No interval
  attaches to either fold**, and none is manufactured.
- **Neither pair is isogenic.** In both papers the derivative arm gains a vector and its
  selection marker — YIplac128 (LEU2) over a *leu2,3-112* host in Verwaal, pHV-BTS1 (HIS3) in
  Ukibe — so each pair also differs in one complemented auxotrophy and in the supplementation
  its medium requires. That is a real confound and it is a field on the record, not a footnote.
- **One paper is one laboratory.** Removing between-paper confounds is not replication.
- **Ukibe's values are the paper's own approximations** ("approximately 18", "approximately
  390 µg per g"), so its 21.67× is coarse by construction.

**What survives all of that** is narrow and it is the point: a *single* Southern-confirmed
added copy moves β-carotene content **11.81×** with host, promoter, medium, cultivation and
assay held constant, and carries the strain from inside the model's flat genotype-independent
ceiling to 4.74× outside it. A ceiling that no genotype can exceed is not compatible with that,
and the incompatibility is not attributable to any of the confounds that a between-paper
comparison would leave open.

### 4.5 Provenance, and the level these numbers get

**The control arms are not in any vendored file.** They were transcribed here, on 2026-09-09,
from the open-access full text, and they are **this repository's artifact** — every `Arm`
carries `provenance = "transcribed_here_from_open_access_full_text"`. That is deliberately
*not* the level a checksummed upstream copy gets. A transcription must never be recorded as a
verified local file.

What makes them checkable rather than asserted: each contrast declares **which of its two arms
the vendored survey holds** — Verwaal's row is the *derivative*, Ukibe's row is the *control* —
and `survey_agreement` binds the contrast to the survey by PMID and **raises** if the two
independent transcriptions disagree by more than the survey's 1% rounding. Verwaal agrees to
0.30% (5.9 stored against 5,918 µg/g printed), Ukibe to 0.00%. A silent edit to the survey
breaks this script instead of quietly moving a result.

### 4.6 A third within-paper contrast, declined

Xie 2014's row records a four-rung internal copy ladder in its `copy_evidence`
(YXWP-14 0.03 → YXWP-40 0.959 → YXWP-47 1.947 → YXWP-53 8.005 mg/gDCW). It is **declined**,
with the reason computed from the row rather than remembered: the rungs were never transcribed
as rows, and that row's `measurand` is `both_reported`, so which of its two assay classes each
rung belongs to is unrecorded. A ladder scored across two assay classes is precisely the
mistake `EXTERNAL_CAROTENOID_BOUND.md` §5 records López 2019 falling into. Transcribing the
four rungs with their per-rung assay and genotype would make it usable, and would be the
third matched contrast in the literature.

### 4.7 Why the between-paper alternatives cannot replace these — verified, with one correction

**Cassette dosage is unidentified between papers.** Eight of the nine `dfba_ready` rows state
a numeric *crtE* dosage, over three levels (0, 1, 2). **79.5% of the log-content variance sits
WITHIN dosage levels.** The prior pass reported "R² = 0.205"; that figure is right but is the
one-way **factor** R² (1 − 0.795), not a linear regression. Regressed **linearly** in dosage
the R² is **0.173**, lower, because the highest content in the whole survey — Arhar 2024 at
79 mg/gDCW — belongs to a strain carrying **zero *crtE***, BTS1 having replaced it. Both
numbers are reported so the distinction cannot be lost again. No pair of the nine rows shares
host, promoter, carbon species and measurand.

**The published table cannot run a leave-one-published-strain-out validation, and this was
run rather than argued.** `pathway/flux.py::score_product_validation` was given the most
generous frame the survey can produce: one strain per paper, the content read as a rate, and a
**fabricated** μ of 0.1777 /h — the table states no growth rate for any row — so that the
failure cannot be blamed on a missing μ. Result: **9/9 outer folds `failed`**, all 58
candidates offered and none scored, `failure_reason = "no candidate has a complete finite
inner validation score"`. The 14 expression channels and the lycopene channel do not exist in
this table and no generosity conjures them.

*Correction to the prior read-only pass: it reported this as **7/7**. Over all nine
`dfba_ready` rows it is **9/9**. The reason and the candidate count are unchanged.*

### Reproduce

```
python3 scripts/matched_contrast.py                       # the report above
python3 scripts/matched_contrast.py --out <path>.csv      # bracket + contrasts as one table
python3 -m pytest tests/test_matched_contrast.py
```

The script writes no artifact unless `--out` is given, so nothing in `outputs/` is adopted by
running it; adoption is a separate operation under [`REPRODUCING.md`](REPRODUCING.md).

## Still to do
- ~~**Scalcinati 2012** (PMID 22938570): α-santalene, FPP-derived, two dilution rates, six
  strains. Chemically nearest external product.~~ **Closed 2026-09-01: out of scope.**
  The cultures carry a dodecane overlay and santalene is assayed in the organic phase, so
  the product is secreted and `PathwaySpec` refuses it at load. Farnesol is figure-only.
  `docs/SECOND_DATASET_HUNT.md` carries the three blockers and the quoted method.
- **Jouhten 2008** (PMID 18613954): five oxygen levels at fixed D. The only clean test of
  an environment axis this model claims and has never been scored on.
