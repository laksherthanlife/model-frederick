# Is there a second dataset that can calibrate the flux law?

**Searched 2026-08-28.** Eight product classes, independently, with adversarial verification
on the shortlist. 32 unique candidates, 5 shortlisted, 4 survived. Every PMID below was
resolved against PubMed by this repository, not taken on an agent's word.

## A second negative, and it closes the cheapest statistical fix

**Searched 2026-08-28.** Before hunting for a second *dataset*, the cheaper question is
whether the dataset already used has a fourth strain nobody pulled. It does not.

The whole of the flux law rests on **three strains**, and three strains cap what any
analysis of it can ever show: there are 3! = 6 assignments of strain label to expression
profile, so the smallest attainable permutation p is **1/6 = 0.167**. The observed
leave-one-strain-out skill of +0.626 is rank 1 of those 6 — the best possible result — and
therefore sits exactly at the floor. No re-analysis moves it. A **fourth** strain takes the
floor to 1/24 = 0.042, which is the first value that can cross 0.05, and it is by a wide
margin the cheapest thing that would make the headline claim significant rather than merely
best-of-six.

So: is there a fourth strain in the authors' own release? `github.com/SysBioengLab/BcarGRASP`
was cloned and read. **There is not.**

| file | strain conditions it carries |
| --- | --- |
| `biomass_and_metabolites_data/BiomassMeasurements.xlsx` | sheets `2D01 2D025 3D01 3D025 4D01 4D025` |
| `biomass_and_metabolites_data/MetabolitesMeasurements.xlsx` | columns `2D01 … 4D025` on every sheet |
| `data_pre-processing/fluxes/output/FluxesBiomassStatistics.xlsx` | rows `2D01 … 4D025` |
| `data_pre-processing/transcripts/output/images/heatmaps/` | `D010_2D010`, `D010_4D010`, `D025_2D025`, `D025_4D025` |

Six conditions, three strains, two dilution rates — the same six this repository already
vendors. `ERG9BiomassStatistics.xlsx` and `estimateERG9Flux.m` look like a fourth strain and
are not: they are an **estimated squalene-synthase flux** across the same six conditions, the
native branch competing for FPP. Worth having for a different reason — nothing here currently
models the competing drain — but it is not a strain.

**Consequence.** The fourth strain is a chemostat run, not a data pull, and it is not on any
of the five specifications in [`MEASUREMENTS_NEEDED.md`](MEASUREMENTS_NEEDED.md), all of
which target the environment axis, the ceiling or the growth exponent. It belongs beside
them.

## A candidate the 2026-08-28 search missed, found 2026-09-01

**It searched eight product classes and all of them were small molecules.** Recombinant
PROTEIN was not among them, and that is where the design this document asks for actually
exists.

**Kastberg et al. 2025**, *FEMS Yeast Res* foaf007, PMID `39971732`, PMC11881926, open
access. *Komagataella phaffii*, glucose-limited **chemostat at D = 0.1/h**, **five production
strains** — I1G, I1S, I6G, M1G, M6G — varying the heterologous gene (human insulin
precursor, Mambalgin-1), the promoter (P_GAP_ / P_SPI1_) and the **gene dosage (one cassette
or six)**. RT-qPCR of the heterologous transcript and secretome proteomics **on the same
chemostat samples**. Raw data deposited: BioProject `PRJNA1144800`, PRIDE `PXD055378` and
`PXD055382`.

Compare that against the ask at the foot of this document — five-plus producing strains, one
dilution rate, glucose-limited chemostat, expression sampled from the same steady state as
the product. It matches, in a different genus and for a protein rather than a small molecule.

**And what it reports falsifies the flux law for secreted protein.** The numbers below are
pulled from the supplementary tables themselves (`foaf007_supplemental_file.docx`, tables
S7-S9 for transcript, S13-S14 intracellular, S15-S16 secretome), not from the running text.
Every value is a label-free log2 fold change. The supplement anchors each strain against the
GS115 reference, so pairwise contrasts subtract exactly, and that is a usable internal check:
13.50 - 9.39 = 4.11 and 7.77 - 3.61 = 4.16 and 11.23 - 9.45 = 1.78 reproduce the three
figures the text quotes as 4.12, 4.15 and 1.8.

| contrast | transcript | intracellular | secreted | ic/tx | **sec/tx** |
| --- | ---: | ---: | ---: | ---: | ---: |
| I1S vs I1G (promoter) | +3.35 | +4.11 | -0.72 | 1.23 | **-0.21** |
| I6G vs I1G (dosage) | +2.90 | +2.75 | +1.07 | 0.95 | **0.37** |
| M6G vs M1G (dosage) | +2.25 | +4.16 | +1.78 | 1.85 | **0.79** |

`flux = alpha * expression` asserts an exponent of **1**. Transcript to INTRACELLULAR protein
sits near it — 0.95, 1.23, 1.85 — so translation tracks transcription about as the law
expects. Transcript to SECRETED protein does not: **-0.21, 0.37, 0.79**, one of them
negative. More message makes more protein, and the extra protein does not leave the cell.

**That is a statement about the OUTLET, and it is why `solve_pathway` gained a saturating
export term on 2026-09-01 rather than a bent flux law.** The distinguishing observable is the
intracellular pool, and this paper is one of the few that measured it in the same samples as
the secretome. A first-order export, `v_sec = k*X`, secretes in proportion to the pool by
construction, so no rate constant reproduces a 17x pool beside a flat supernatant. A capacity
does. See `NodeKinetics.secretion_km_mmol_per_gdcw` and `tests/test_secretion_outlet.py`.

**One number in this paper contradicts itself, and it changes an exponent.** Supplementary
table S9 gives the M6G vs M1G transcript contrast as **2.25** log2FC; the Results text gives
**4.74** for what reads as the same comparison. The table is used above, giving an exponent of
0.79; on the text's figure it would be 0.38. Sublinear either way, and not pinned.

**What this data still cannot do.** It cannot calibrate the outlet it motivated. Kastberg
reports **no absolute quantity of product anywhere** — no mg/L, no g/L, no mg/gDCW, no qp, no
ELISA, no product HPLC, in the body, the PDF or any of the sixteen supplementary tables; the
HPLC ran on glucose, ethanol, acetate and CO2 only. Every product value is a relative fold
change, and a relative change fits neither `vmax` nor `km`. The deposited data does not lift
this: **neither PRIDE deposition contains a processed quantitative matrix** (42 files each,
all raw/peak/mzid — verified against the v3 REST file listing and the FTP directory), and
**no GEO series exists** for the BioProject, so only raw reads are deposited on that side.

What Kastberg does supply is per-strain physiology at the same steady state, in absolute
units, which is the part a model can use (table S4, D = 0.1/h, two timepoints):

| strain | residual glucose (g/L) | biomass yield (gCDW/g glucose) | CO2 yield (gCO2/gCDW) |
| --- | ---: | ---: | ---: |
| GS115 (reference) | 0.140 / 0.135 | 0.58 / 0.60 | 0.72 / 0.67 |
| I1G | 0.143 / 0.144 | 0.62 / 0.54 | 0.70 / 0.78 |
| I6G | 0.142 / 0.140 | 0.60 / 0.56 | 0.68 / 0.72 |
| I1S | 0.152 / 0.149 | 0.54 / 0.56 | 0.74 / 0.71 |
| M1G | 0.155 / 0.161 | 0.51 / 0.51 | 0.78 / 0.79 |
| M6G | 0.138 / 0.139 | 0.52 / 0.55 | 0.77 / 0.73 |

**The calibration that is missing, and where one point of it exists.** Fitting a saturating
export needs paired intracellular content and supernatant accumulation rate on the same
culture at two or more expression levels. **Pfeffer 2011** (PMID `21703020`) has exactly that
shape and exactly one point of it: 34-S pulse labelling on a *K. phaffii* chemostat at the
same D = 0.1/h, giving intracellular 91.8 ug/gYDM against qSec 50.7, qDeg 83.3 and qDil 9.18
ug/gYDM/h. Those close on their own arithmetic — the three outlets sum to 143.18 against a
measured formation rate of 143.1, and mu times the pool is 9.18 exactly — so it is a real
three-outlet balance and not a set of separately reported figures. It is one state, so it
fixes the ratio `vmax/(km + X)` and neither constant. **One more dilution rate closes it.**

Its headline is also a warning about what the lumped constant covers: **58% of the product
made is degraded inside the cell, 35% is secreted, 7% is inherited by daughter cells**, and
degradation is the FASTER half-life (45.8 min against 75.3). The dominant sink for a
secretory protein is proteolysis, not export.

---

## The headline was a negative, and on 2026-09-02 it was falsified

**This section used to open: "No second chemostat expression series exists, in any product
class."** That is false, and the counterexample is in the same product class as the chain's
own calibration.

**GEO `GSE8451`** — "Transcriptional profiling of carotenoid producing S. cerevisiae cells",
published as **Verwaal et al. 2010, PMID `20632327`** (*Yeast*, "Heterologous carotenoid
production in Saccharomyces cerevisiae induces the pleiotropic drug resistance stress
response"). Verified against the GEO record directly, not inferred: **14 arrays**,
*S. cerevisiae*, **carbon-limited continuous culture**, three genotypes — CEN.PK113-7D as the
reference plus two carotenoid producers, Orange01 and Orange02 — under two feed regimes,
glucose-limited and glucose/ethanol-limited.

That is a chemostat expression series on carotenoid producers. The old sentence claimed no
such thing existed anywhere.

**It clears the Vos 2015 bar for ONE of its two producers and not the other, and the design
is nested.** This paragraph said something stronger until 2026-09-04, and the correction is
recorded rather than overwritten because of how the error was made: the claim was read off
GEO's SERIES SUMMARY, and the refutation is in the PER-SAMPLE characteristics of the same
record. A dataset assessment must quote the per-sample field it read, per blocking check.
That is the durable protection here — this section describes an external dataset with no
in-repo artefact, so no `audit:value` marker can reach it.

It read: *"the entry step is GGPP synthase, and the producers overexpress BTS1, which is the
native S. cerevisiae gene. A genome-wide yeast array carries a probe for a native gene by
construction."* Pulled from `GSE8451_family.soft.gz` directly, the two producers are not the
same construct:

| strain | genotype, from the per-sample characteristics | GGPP synthase |
| --- | --- | --- |
| Orange01 | `ura3-52::PTDH3-crtYB; PTDH3-crtI; PTDH3-BTS1` | native **BTS1** |
| Orange02 | `PTDH3-crtYB; PTDH3-crtI; PTDH3-crtE` | heterologous **crtE**, *X. dendrorhous* |

So blocking check 1 is answered **no for half the panel**. Orange02's entry enzyme is exactly
the kind of heterologous gene that is absent from GPL90 by construction — the same objection
that rejected Vos 2015, arriving by the same route.

**And the design is NESTED, not crossed, which is the sharper limit.** All four CEN.PK113-7D
and all three Orange01 arrays are **glucose-limited**; all four CEN.PK113-7D and all three
Orange02 arrays are **glucose/ethanol-limited**. The only strain measured at both feeds is the
NON-PRODUCING reference. Carbon feed is therefore perfectly confounded with producer genotype:
no carotenoid producer in this series is measured at two feeds, and the two producers cannot
be compared to each other at all. The "two feed regimes" in the summary are not a
carbon-source contrast in the sense `pathway/environment_flux.py` fits — the second is 7.5 g/L
glucose plus a **9.38 mM ethanol anti-oscillation supplement**, at D = 0.10 h⁻¹.

**What is still NOT verified, and must be before Orange01 is used:**

1. That BTS1's signal on GPL90 separates Orange01 from the reference. (Its PRESENCE on a
   genome-wide array is no longer in question; its dynamic range is.)
2. That per-strain carotenoid content is reported in absolute units at these steady states.

**What survives.** Not "five producers across two laboratories". **Orange01 alone, at one
feed, at one dilution rate** — at most ONE additional calibration point beside Elizondo's
three strains, and one whose entry-enzyme axis has two levels (reference and producer) rather
than a dosage series. Orange02 is out on the array objection; the feed axis is out on the
confounding. That is a partial retraction of a partial rescue, and it leaves the second-dataset
question closer to where it started than the previous version of this section suggested.

The structural pattern below is still what the search kept returning, and it still explains
why generalising has been hard. It is a strong tendency, not the absolute negative it was
written as:

- **Chemostat papers vary the ENVIRONMENT on one genotype.** That is the axis this repository
  already refuted (11 states, every model negative leave-one-out skill, inside the
  permutation null).
- **Genotype-panel papers are shake flask or microplate.** They report titre, not `q`, and
  no growth rate.

Elizondo 2025 (PMID 40891387) sits at the intersection of the two, and that intersection is
close to empty. The dataset the chain is calibrated on is **near-unique**, not merely the
first one we happened to find.

## Why that bites this specific architecture

The chain is `genotype -> flux -> steady-state pools -> product`, and the pools layer is
`X = v_in / mu`. **Growth rate is not optional** — it is the denominator. A titre-only
dataset can exercise `expression -> something` but not the law we actually fitted, because a
titre is time-integrated and confounded by biomass and duration, where `q` is neither.

So every candidate below tests a **weaker claim than the one β-carotene supports**, unless an
assumption about mu is added and declared.

## The four that survived

Ranked here for **fit to this chain**, which is not the same as the ranking on data quality.

| | Product / host | Identifier | n | Product units | mu | The catch |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Triacetic acid lactone, *K. marxianus* | PMID 32995271 | 18 (6 genotypes × 3 temps) | **specific, biomass-normalised** | per temperature (0.28/0.34/0.35 /h), not per strain | promoter strengths measured in separate EGFP strains — paired by promoter identity, not on the producing cell |
| 2 | Caffeic acid, *S. cerevisiae* | PMID 35300487 | 21 | titre (mg/L) | none | no numeric table; RT-qPCR is Fig 6A/6B and titre is Fig 5B, both bitmaps. Six δ-strains span only ~1.4-fold |
| 3 | (2S)-Naringenin, *E. coli* | PMID 35346204 | 34 | titre (mg/L) | none | **n = 1 per titre, no replicates**. Expression axis needs a cross-paper join to Bervoets 2018 (PMID 29361130) |
| 4 | Lycopene, *E. coli* | PMID 16123130 | ~27 | titre | qualitative only | figure digitisation only, and it is another carotenoid, so it is weakest for "beyond β-carotene" |

**TAL is first for us** because it is the only one whose product is biomass-normalised — a
quantity in the same class as `q` — and the only one with any mu at all. Its two-assumption
cost (batch mu ≈ mu_max per temperature; promoter strength transfers from an EGFP reporter
strain) is smaller than the cost of having no denominator.

**Caffeic acid is first on data quality**: same host, real RT-qPCR (2^-ΔΔCT, ALG9 reference,
triplicate), a two-enzyme expression vector (coC3H + coCPR1) that would test whether α is
per-enzyme or needs a multi-gene term, and an orthogonal copy-number column. It is titre-only,
and near-collinear copy number means the two coefficients will not separate.

## The two that miss by exactly one measurement

Both are glucose-limited chemostats in *S. cerevisiae* — the right shape, one column short:

- **Scalcinati 2012, PMC3527295** — α-santalene, six strains at steady state.
  **Closed 2026-09-01, and not by the column this entry named.** The source was read
  directly and it fails the scope rule before it fails the expression one, three ways:

  1. **The product is secreted.** The cultures are two-phase with a dodecane overlay —
     *"The product is continuously captured in the organic phase due to its high
     hydrophobicity"* — and santalene is assayed by GC-MS **on the organic layer**, not on
     cells. `PathwaySpec` refuses `Fate.SECRETED` at load, correctly: the balance closes on
     `mu*[X]`, and a product partitioning into dodecane leaves at a rate the solvent sets.
  2. **`q` is not printed.** Table 1's columns are `Strain, D, Yxs, rs, rCO2, rO2, retoh,
     racet, RQ, Cs, Cbalance, Totsant` — the last in mg 24h⁻¹ l⁻¹, a volumetric rate. The
     specific productivity appears only in Figures 3 and 4. It is *derivable* from `Totsant`
     and `Yxs`, which would be arithmetic rather than digitisation — but the point is moot
     given (1).
  3. **Farnesol, the intermediate that made this the interesting candidate, has no printed
     number at all** — rates and yields in Figures 3–5 only. So the two-node branch test
     this dataset seemed to offer does not exist without digitising a figure, which is what
     `tests/test_tal_pathway.py` records this repository refusing.

  Zero expression measurement remains true and is now the third reason, not the first.
- **Vos 2015, PMID 26369953** — 15 resveratrol chemostats, `q` per steady state, coumaric /
  cinnamic / phloretic acid pools. One genotype, and the microarray carries **no probes for
  the heterologous PAL/4CL/VST1**.

## The wet-lab ask, now exact

> RT-qPCR (or targeted proteomics) of the **pathway entry enzyme** on **6 producing strains**
> held at **one dilution rate** in **glucose-limited chemostat**, sampled from the same steady
> state as the product measurement.

**Six, not "≥5", and the number is now computed rather than guessed.**
`scripts/design_flux_experiment.py` simulates the panel at the scatter the real CrtE fit
produced (leave-one-strain-out RMSE 0.2002 in log space, 1.22x typical error) and reports the
power of a strain-label permutation test at each size:

| strains | permutation floor 1/n! | power at alpha = 0.05 |
| ---: | ---: | ---: |
| 3 | 0.167 | **unreachable** -- the floor is above alpha |
| 4 | 0.042 | 0.44 |
| 5 | 0.008 | 0.68 |
| **6** | **0.001** | **0.93** |
| 8 | <0.001 | 0.98 |

Five was the number this paragraph asked for and it reaches about two thirds, so a third of
such experiments would come back inconclusive having cost a chemostat run. Six is the first
size clearing 80%.

**What that simulation does not do is validate the law**, and the script says so in its own
docstring: its generator and its model are the same equation, so a good fit measures the
fitter. It is Tier 1 by this repository's grading, it sizes the experiment, and it cannot
replace it.

Nothing else. Every candidate above already has the product side; what none has is a measured
entry-enzyme level on a strain panel at a set mu. The cheapest version: take an existing
strain panel into chemostat at a single D and qPCR the entry enzyme — that converts a titre
series into a `q` series and restores the denominator.

## Metadata corrections found during verification

Record these; do not propagate the originals.

- **Qi 2022 (PMID 35300487)** — host is **BY4742**-derived, not BY4741. Paired n is **21**,
  not 24. The 7.9 mg/L for D9 is **YPD**, not SC-Ura.
- **Van Brempt 2022 (PMID 35346204)** — usable n is **34**, not 35 (one strain failed
  sequencing). The promoter TIF values are *not* absent from the paper: Additional file 1,
  Fig. S2A/S2C.
- **Alper 2005 (PMID 16123130)** — n is **~27**, not 18, and there is a published correction
  in PNAS 2006 that any use of this dataset must apply.

**No candidate was rejected for a fabricated identifier.** All five shortlisted PMIDs resolved
to the claimed paper with matching title, journal, year and first author.
