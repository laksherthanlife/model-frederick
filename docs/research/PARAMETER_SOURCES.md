# Parameter sources — what the literature can actually supply for `stress_panel.py`

**Compiled 2026-08-26.** Follows `docs/SOURCE_AUDIT.md` §3.8 and §6, which established that
`generator/stress_panel.py` holds **136 numeric parameters** (76 stressor target weights, 25
EC50s, 25 lethal doses, 6 `also_reads` crosstalk weights, 4 cascade weights), of which two are
measured and none carries a per-number citation. This document does not repeat the audit. It
asks the next question: **which of those numbers does the published literature contain, and at
what confidence.**

Every PMID below was resolved through NCBI E-utilities `esummary.fcgi` and the returned title,
first author, year and journal are quoted in §7. Nothing is cited that was not resolved.
Where the literature does not have the number, this document says so.

---

## 1. Verdict

**Counted strictly against the 136: eight carry a published numeric value, two more are
measured on this project's own plates, twenty-three are not parameters at all, and the
remaining 103 have no number in the literature in any form a reader could cite.**

The largest structural finding is arithmetic: **23 of the 25 `lethal_dose` values are exactly
`_LETHAL_MULTIPLE × ec50`** (verified by running the module). They encode one number, `3.0`.
The free parameter count is 113, not 136.

Of the 113 free numbers: **7 EC50s** have a real published anchor (H2O2, DTT, diamide, sorbitol,
NaCl, menadione, heat), **1 crosstalk weight** has a directly measured value
(`RIM101-alkaline.also_reads["calcium"]`, and the measurement says the encoded value is wrong by
5×), **2 lethal doses** are measured in-repo, **5 more crosstalk weights and all 4 cascade
weights** have solid primary sources for the *edge* and none for the *magnitude*, and the
**76 stressor target weights have nothing** — not one published number.

**No published quantitative regulon-overlap matrix exists that could be dropped in wholesale.**
Every candidate in the brief was chased. Gasch 2000 and Causton 2001 are condition × gene
expression matrices, not regulon × regulon; Harbison 2004 is a TF × promoter binding matrix
(useful, and the closest available object for `also_reads`); YEASTRACT+ returns a curated,
unweighted TF × gene matrix; Balaji 2006 is the one published *TF × TF co-regulation network*
for yeast, but it scores association significance, not signed loadings. The weights are
**derivable** from public data — the recipe is §3.2 — but somebody has to run it. That is a
half-day of data reduction, not a literature lookup.

One correction is available immediately and is not a magnitude fix but a topology fix:
`TRX2-oxidative.also_reads["ESR"]` attributes to Msn2/4 an effect the literature attributes to
**Skn7**, which the panel already places *inside* the oxidative module (§3.3).

Outside the 136, the culture and reporter constants fare much better: `_MU_MAX`,
`carrying_capacity`, `integrates_hours`, `reporter_dilution_exponent` and the design noise model
can all be grounded today, four of them from BioNumbers entries already vendored in this
repository (§2.5, §4).

---

## 2. The triage table

Category key: **1** = a real paper reports this quantity for *S. cerevisiae*; **2** = measurable,
and the data to measure it is public, but nobody has published the number in this form;
**3** = structural — it encodes a topology, a normalisation or a units choice, and asking for
its "measured value" is a category error.

### 2.1 The parameters with the largest blast radius

`docs/ARCHITECTURE.md` §2.4 ranks `stress_panel.py` first (fan-in 8 in-package / 40 total).
Within it, `reporter_loadings()` builds `L`, so the six `also_reads` weights are the whole
off-diagonal of the matrix every identifiability claim rests on.

| # | Parameter | Current | Cat | Measured value found | Source | Confidence |
|---|---|---|---|---|---|---|
| 1 | `REPORTERS["RIM101-alkaline"].also_reads["calcium"]` (ENA1) | **0.30** | **1** | **Crz1/calcineurin supplies ≈ 60 % of the early ENA1 alkaline-pH response**; two stress-responsive Crz1 sites in the promoter; the residual is Rim101(→Nrg1) + Snf1 | Petrezsélyová 2016, PMID 27362362; topology from Platara 2006, PMID 17023428 | **High** — the paper is titled "A Quantitative Study" and models the two inputs explicitly |
| 2 | `REPORTERS["TRX2-oxidative"].also_reads["ESR"]` | **0.30** | **3 → wrong module** | The Yap1-independent arm at the TRX2 promoter is **Skn7**, which binds TRX2 directly and co-operates with Yap1 — not Msn2/4 | Morgan 1997, PMID 9118942; Kuge & Jones 1994, PMID 8313910 | **High** on the attribution, **none** on the number |
| 3 | `REPORTERS["HSE-heat"].also_reads["ESR"]` (HSP104) | **0.35** | **2** | HSEs and STREs *cooperate* for maximal induction; in `msn2Δmsn4Δ` proper heat induction is obtained **exclusively through HSEs**; in `ras2Δ` derepression is **exclusively through STREs** | Grably 2002, PMID 11967066 (β-gal values are in the figures, not the abstract) | Direction **high**, magnitude **not published in extractable form** |
| 4 | `REPORTERS["UPRE-ER"].also_reads["heat"]` (KAR2) | **0.25** | **2** | KAR2's promoter carries a functional **HSE**; HSE deletion abolishes induction on cytosolic Hsp70 loss, and either the KAR2 HSE or the SSA1 HSE alone drives an HSE-CYC1-lacZ fusion | Kohno 1993, PMID 8423809; Oka 1997, PMID 9133628 | Topology **high**; the heat:DTT induction *ratio* is computable from GSE18 but unpublished |
| 5 | `REPORTERS["STRE-osmotic"].also_reads["ESR"]` (GPD1) | **0.40** | **2** | 186 genes ≥ 3-fold induced by NaCl/sorbitol; Hot1p carries a subset, Msn2/4 a different subset, and Msn2/4-dependent induction is *itself* reduced in `hog1Δ` | Rep 2000, PMID 10722658; Rep 1999, PMID 10409737 | Topology **high**; per-gene split is in the JBC tables |
| 6 | `REPORTERS["CSRE-carbon"].also_reads["ESR"]` (ADH2/ICL1) | **0.25** | **2** | Snf1 co-regulates Adr1 and Cat8 targets; glucose withdrawal is also a core ESR trigger | Young 2003, PMID **12676948** (note: the file still cites 12482873, a cyclin D3 paper — see §5) | Topology **high**, magnitude **unpublished** |
| 7 | `MODULES["proteasome"].driven_by` = `{oxidative: 0.45, heat: 0.40}` | 0.45 / 0.40 | **2** | YRE **and** HSE in the RPN4 promoter, both confirmed by site mutagenesis; Hsf1 co-ordinates proteasome expression | Hahn 2006, PMID 16556235; Owsianik 2002, PMID 11918814 | Both edges **high**; the near-equality (0.45 vs 0.40) is an assertion |
| 8 | `MODULES["iron"].driven_by["oxidative"]` = 0.10 | 0.10 | **2** | One conserved YRE 147 bp upstream of *AFT2*'s ATG, plus Yap1 ChIP-chip binding at the AFT2 promoter | Salin 2008, PMID 18627600 (verified in the audit) | Edge **high**, weight **unpublished** |
| 9 | `MODULES["xenobiotic"].driven_by["proteasome"]` = 0.20 | 0.20 | **3** | A positive transcriptional loop connects *RPN4* and *PDR1*. A loop **gain** is not a promoter-content weight and no paper reports one | Salin 2008, PMID 18627600 | Edge **high**, the number is **structural** |
| 10 | `NoiseModel.relative_cv` (`analysis/design.py:42` = 0.02) and `_PLATE_NOISE` (`sensor_selection.py:180` = 0.05) | 0.02 / 0.05 | **1 — in-repo** | The repo's own plates measure **0.146** (`panel_experiment.OBSERVED_ACTIVITY_CV`) and **0.14** (`MEASURED_ACTIVITY_CV`) | This repository, 28 matched conditions on two real plates | **High.** See §4 — this is the cheapest correction in the file |

### 2.2 The 25 EC50s

Running `module_ec50` shows the encoded EC50 is not one number per stressor but three: the
agent's `ec50` for its specific regulon, `× _POOL_POTENCY (0.35)` for a metabolite pool, and
`× _GENERAL_POTENCY (2.5)` for the ESR. So the 25 EC50s generate 76 module-level half-maximal
doses through two multipliers.

| Stressor | `ec50` | Cat | What the literature says | Source | Confidence |
|---|---|---|---|---|---|
| **H2O2** | **0.5 mM** | **1** | **Yap1-GFP nuclear relocation is partial at 0.1 mM and saturates above 0.2 mM**; 0.1 mM has no effect on growth rate; cell-cycle arrest emerges ~0.2 mM, growth arrest ~0.5 mM, **full growth arrest at 0.6 mM**. S288C in SCD, microfluidic sustained delivery | Goulev 2017, PMID 28418333 | **High**, with a readout caveat (Yap1 nuclear entry is signalling, transcription may lag) |
| H2O2 (conventional dose) | — | **1** | **0.30 mM** (Gasch 2000, DBY7286/YPD); **0.4 mM** (Causton 2001, S288c/YPD; Godon 1998) | PMIDs 11102521, 11179418, 9712873 | High — but a *dose a paper chose* is an upper bound on an EC50, not an EC50 |
| **DTT** | **1.0 mM** | **1** | **Half-maximal UPR at ≈ 2.2 mM** (≤ 1.5 mM gives < 10 % of max; 3.3 mM ≈ 75 %; 5 mM near-saturating). W303a in 2×SDC, UPRE-GFP + HAC1 splicing | Pincus 2010, PMID 20625545 | **High** for W303/SDC; see §5 for why it may not transfer to BY4741 |
| DTT (conventional dose) | — | **1** | **2.5 mM**, in **BY4741** in YPD — "enough stress to reduce cell growth but enables acclimation at a new growth rate" | MacGilvray 2020, PMID 32597660 (= the `PMC7646510` the file cites) | High |
| diamide | 1.5 mM | **1** | **1.5 mM** — exactly the Gasch 2000 dose | PMID 11102521 | High (dose, not EC50) |
| sorbitol | 1.0 M | **1** | **1 M** (Gasch 2000); 1.5 M (Causton 2001) | PMIDs 11102521, 11179418 | High (dose) |
| NaCl | 0.5 M | **1** | **1.0 M** (Causton 2001); Gasch used 0.4–1 M | PMIDs 11179418, 11102521 | Medium — 0.5 M sits inside the conventional band |
| menadione | 100 µM | **1** | Gasch 2000 used **1 mM menadione bisulfite** — 10× the encoded EC50 | PMID 11102521 | Medium; salt form and cell density both matter |
| heat | 6 °C above 30 (= 36 °C) | **1** | Every ESR dataset shifts from **25 → 37 °C** (Gasch, Causton) or **25 → 38 °C** (Boy-Marcotte). Gasch also shows 29 → 33 °C works | PMIDs 11102521, 11179418, 10411744 | Medium — see §5 on the parameterisation |
| MMS | 0.02 % | **2** | The Hahn 2006 dose remains **UNVERIFIED** (paywalled; the audit tried too) | PMID 16556235 | Low |
| tunicamycin | 1 µg/mL | **2** | Standard doses span 0.5–5 µg/mL and are strongly strain- and medium-dependent; no EC50 found | — | Low |
| rapamycin | 20 nM | **2** | 200 nM is the conventional dose; no EC50 for the Gln3 arm found | — | Low |
| The other 14 (BPS, MG132, hydroxyurea, congo red, caffeine, CaCl2, cobalt, copper sulfate, TPEN, methionine withdrawal, fluconazole, acetic acid, alkaline shift, antimycin A, glucose withdrawal) | various | **2** | Conventional doses exist in the literature for all of them; **published EC50s do not.** Each would need a dose ladder with a regulon-specific reporter | — | Low. These are also the stressors nobody in this project will dose — lowest priority by design |

### 2.3 The 25 lethal doses

| Parameter | Current | Cat | Finding |
|---|---|---|---|
| `STRESSORS["DTT"].lethal_dose` | 1.55 mM | **1 — in-repo** | Measured on this project's plates (growth halves at 1.45 and 1.65 mM). **But see §5**: MacGilvray 2020 grew BY4741 at 2.5 mM DTT with only a growth-rate reduction |
| `STRESSORS["H2O2"].lethal_dose` | 1.0 mM | **1 — in-repo** | Measured (halves at 0.79 and 1.24 mM). Goulev 2017 reports **full growth arrest at 0.6 mM** under sustained microfluidic delivery — the discrepancy is mechanistically explained by Tomalin 2016 (PMID 26944189): a bolus is consumed by the cell's thiol buffering, a perfusion is not |
| The other **23** | `3.0 × ec50` | **3** | **Not parameters.** Verified by execution: 23 of 25 `lethal_dose` values equal `_LETHAL_MULTIPLE × ec50` to machine precision. They encode `_LETHAL_MULTIPLE = 3.0` and nothing else. The audit's count of 136 should be read as **113 free numbers + 23 aliases** |

### 2.4 The 76 stressor target weights

All 76 are **category 2**. No paper publishes them; the data to compute them is public and
identified in §3. Two structural notes that reduce the problem:

- **~24 of the 76 are the `ESR` entry**, which appears in 24 of the 25 stressors with values
  0.20–1.00. These are a single quantity — *how strongly does this agent drive the common
  environmental stress response* — and Gasch 2000's 283 induced ESR genes give exactly that,
  per condition, as a mean log2 ratio, from one dataset (GEO **GSE18**).
- **Signs are structural, magnitudes are not.** `DTT → redox: −0.90` and `H2O2 → redox: +0.70`
  encode that a reductant and an oxidant displace the same pool in opposite directions. The
  sign is a fact; the 0.90 and 0.70 are not measurements.

### 2.5 Culture, reporter and design constants (`culture.py`, `literature.py`, `design.py`)

| Parameter | Current | Cat | Measured value | Source | Confidence |
|---|---|---|---|---|---|
| `_MU_MAX` (`stress_panel.py:267`) / `mu_max_glucose` (`literature.py:37`) | 0.40 h⁻¹ | **1** | **0.37 h⁻¹** on minimal medium at 30 °C (**BNID 106359**); doubling **~90 min in YPD → µ = 0.46 h⁻¹** and **~140 min in synthetic → µ = 0.30 h⁻¹** (**BNID 108255**); typical 90–120 min (**BNID 110545**) | `data/kaggle/BioNumbers_Nov2024.csv`, already vendored | **High.** 0.40 is the YPD figure. If the plates are synthetic, 0.30 is the right anchor |
| `_CONSTRUCTS[*]["mu_max"]` (`literature.py:92–103`) | 0.31–0.40 | **1** | Same BNIDs. Also relevant: **BY4741 and W303 grow poorly on synthetic complete because `leu2` cannot take up leucine** — a strain-specific growth penalty that is not a stress effect and will be read as one | Cohen & Engelberg 2007, PMID 17573937; Brachmann 1998, PMID 9483801 | **High**, and directly actionable |
| `carrying_capacity` (`literature.py:94–103`) | 3.0 / 2.6 g/L | **1** | **Max OD600 = 5.2 in YPD** and **6.87 in YPGal** (**BNID 106303 / 106302**); **1 OD600 = 0.62 gDCW/L** (**BNID 108275**) → **3.2 g/L**. Total biomass yield 54 g/mol glucose = 0.30 gDW/g (**BNID 110953**) → 6.0 g/L ceiling on 2 % glucose if fully respired | vendored BioNumbers | **Medium-high.** 3.0 g/L is well inside the range. Caveat: a 200 µL well is O2-limited, so full respiration of ethanol will not happen — the fermentative-only yield (~0.15 g/g → 3.0 g/L) is the better match, and that is where 3.0 lands |
| `growth_hill`, `promoter_hill` (`literature.py:92–103`) | 1.1–2.2 | **2** | **No published Hill coefficients for yeast stress growth-inhibition or promoter induction were found.** They are fitted shape parameters. `panel_calibration.fit_dose_response` already exists to fit them from a ladder | — | **None** |
| `_LETHAL_HILL` (`stress_panel.py:268`) | 2.5 | **2** | Same. Not published; fittable from the repo's own growth ladders, which `lethal_dose_from_growth` already does | — | **None** |
| `promoter_basal` = 1.0e-3, `promoter_peak` = 1.45–1.60e-3 | arbitrary units | **2** | The right published anchor is **Keren 2013** (~900 *S. cerevisiae* and ~1800 *E. coli* promoters, absolute activities by fluorescent reporter). For units, **BNID 106766**: median transcription rate 0.12 mRNA/min; **BNID 106767**: 90 % of genes between 2.33 and 29.7 mRNA/h | Keren 2013, PMID 24169404 | **Medium** — the mapping to the repo's units is a modelling choice |
| `_DEFAULT_BASAL` (`stress_panel.py:37`) | 0.9 | **1 — in-repo** | Fitted to this project's DTT ladder (1.38–1.47 fold), verified arithmetically by the audit. **The general statement behind it is published**: Keren 2013 finds **60–90 % of promoters change between conditions by a global scaling factor that depends only on the condition, not the promoter.** A high constitutive floor with a small specific component is the normal case, in both yeast and *E. coli* | Keren 2013, PMID 24169404 | **High** for the local fit; the Keren result is the generalisable justification and is currently uncited |
| `Reporter.integrates_hours` (`stress_panel.py:164`) | 3.0 h | **3 — derived** | This is **1/(µ + k_deg)**. At the BNID synthetic-medium µ = 0.30 h⁻¹ it is 3.3 h; at YPD µ = 0.46 h⁻¹ it is 2.2 h. It is not an independent parameter and should be computed from µ, not fixed — which also makes it host-transferable for free | BNID 108255; Hintsche & Klumpp 2013, PMID 24041253 (= the `PMC3847955` the file cites) | **High** |
| `reporter_dilution_exponent` = 1.0 (`literature.py:64`) | 1.0 | **1** | Correct: steady-state concentration goes as α/µ. The cited PMC3847955 resolves to **Hintsche & Klumpp 2013, *J Biol Eng* 7:22**, which is the right paper for this claim | PMID 24041253 | **High** |
| `h2o2_sublethal` = 0.5 mM (`literature.py:60`) | 0.5 mM | **1** | Goulev 2017: at 0.5 mM only **22 % of cells adapt normally**, 36 % show prolonged cell-cycle arrest and **42 % permanently arrest**. Calling 0.5 mM "sublethal" is defensible for a bolus but false for a sustained dose | PMID 28418333 | **High** — and it contradicts the label |
| `min_reporter_above_background` = 50.0 (`gates/g1_optical.py`) | 50 | **1, indirectly** | **BNID 115096**: lower limit for reliable GFP detection ≈ **1,400 molecules/cell** | vendored BioNumbers | Medium — converting molecules/cell to RFU needs the instrument's gain |

---

## 3. The crosstalk matrix — can the hand-assigned weights be replaced by published data?

**Short answer: not by citation, but yes by computation, and the computation is small.**

### 3.1 What was checked, and what each candidate actually is

| Candidate | What it is | Can it replace the weights? |
|---|---|---|
| **Gasch 2000** (PMID 11102521, GEO **GSE18**, SGD dataset page `yeastgenome.org/dataset/GSE18`) | ~6,150 genes × ~150 arrays of log2 ratios across ~12 stresses and time courses. Defines the ESR: ~900 genes, ~300 induced / ~600 repressed | **Not as published.** It is condition × *gene*, not condition × *regulon*. But it is the numerator of the recipe in §3.2 |
| **Causton 2001** (PMID 11179418) | Independent replication: heat, acid, alkali, H2O2 0.4 mM, NaCl 1 M, sorbitol 1.5 M. Defines the CER: **499 genes, 216 induced, 283 repressed** | Same. Its value is as an *independent* second dataset — a weight that agrees across Gasch and Causton is real; one that does not is protocol-specific |
| **Harbison 2004** (PMID 15343339) | Genome-wide binding location for **203 TFs**, in YPD **and 12 other environmental conditions**, as per-promoter binding p-values | **This is the closest published object to `also_reads`.** `also_reads` is *promoter content* — which TFs bind this reporter's promoter — and that is exactly what Harbison measures, condition-resolved. It gives presence and confidence, not a loading magnitude |
| **Lee 2002** (PMID 12399584) | The 106-TF predecessor, YPD only | Subsumed by Harbison |
| **YEASTRACT+** (PMIDs 31586406, 36350610) | Curated TF→target regulations (~175k associations, 183 TFs with documented binding sites) with a **"Regulation Matrix"** tool that emits a TF × gene matrix | **Gives the regulon gene sets** the recipe needs. The matrix is documented/binary, not weighted |
| **Kemmeren 2014** (PMID 24766815) | Expression profiles of **1,484 deletion strains** — for every TF deletion, a genome-wide response with M-values. Distributed as the Deleteome | **The best weighted TF→gene object available.** A TF's regulon becomes a real-valued vector, not a set |
| **Balaji 2006** (PMID 16762362) | The **published TF × TF co-regulatory network** for yeast: all significant associations among TFs regulating common targets | The only paper that publishes a *yeast TF-overlap matrix*. It scores association significance, not signed loading. **Use it as the null/structure check, not as `L`** |
| **Ihmels 2002** (PMID 12134151), **Segal 2003** (PMID 12740579) | Module decompositions of the yeast expression compendium, with condition-specific regulators | Modules are data-derived, not the seven named regulons. Useful as an independent partition to test the encoded one against |
| **Chasman 2014** (PMID 25411400), **Ho & Gasch 2015** (PMID 25957506) | The inferred yeast stress-activated *signalling* network, integrating interactions, fitness, mutant transcriptomes and phosphoproteomes — for **salt stress** | Signalling topology, one stress. Confirms crosstalk hubs; supplies no loadings |
| **yStreX** (PMID 25024351), **PROPHECY** (PMID 15608218) | Curated cross-study yeast stress expression; and quantitative microplate growth phenotypes for the deletion collection under environmental challenge | yStreX is a convenience layer over the same primary data. PROPHECY is the right source for §2.5's growth constants, not for crosstalk |

### 3.2 The recipe that would replace 76 guessed weights with measured ones

This is the highest-value single piece of work identified in this document.

1. **Define each module's regulon as a gene set.** YEASTRACT+ Regulation Matrix, restricted to
   *documented DNA binding + expression evidence*, for the seven TFs the panel names
   (Msn2/Msn4, Hac1, Yap1/Skn7, Hsf1, Hog1→Sko1/Hot1, Rpn4, Aft1/Aft2). Sanity-check the sizes
   against the primary regulon papers — Solís 2016 (PMID 27320198) puts the *core* Hsf1 regulon
   at **18 genes**, ten-fold fewer than the >160 previously assumed, which is a real and
   load-bearing correction to any Hsf1 weight.
2. **Score each stressor against each regulon** as the mean log2 ratio of that regulon's genes
   at the peak time point in GSE18 (Gasch), normalised so the agent's own regulon is 1.0. That
   is a measured `targets` vector, in exactly the model's units.
3. **Replicate on Causton 2001** for the six conditions it shares. Weights that agree across
   both are reportable; weights that do not are protocol artefacts and should be widened, not
   averaged.
4. **Replace the regulon *sets* with Kemmeren 2014 weighted vectors** for the second pass —
   a TF-deletion response is a direct measurement of that TF's contribution per gene, which
   removes the binary in/out decision entirely.
5. **For `also_reads` specifically, do not use expression at all — use binding.** `also_reads`
   is promoter content. Take Harbison 2004's binding p-values at the eight reporter promoters
   (`HSP12`/`CTT1`, `KAR2`, `TRX2`, `HSP104`, `GPD1`, `RPT1`, `ADH2`/`ICL1`, `ENA1`) and set
   the weight from the secondary TF's binding confidence relative to the primary's. This
   matches the semantics of the field, which expression data does not.
6. **Then run the null the repo already owns.** `analysis/nulls.shuffled_loadings` exists for
   exactly this question. Comparing the identifiability of the derived `L` against the encoded
   one converts a provenance argument into evidence, which is what `SOURCE_AUDIT.md` §7 item 8
   asks for.

### 3.3 The correction that must happen before any of that

**`REPORTERS["TRX2-oxidative"].also_reads = {"ESR": 0.30}` names the wrong module.** The
Yap1-independent arm at the TRX2 promoter is **Skn7**, which binds the TRX2 promoter directly
and co-operates with Yap1 on it (Morgan 1997, PMID 9118942; Kuge & Jones 1994, PMID 8313910).
Skn7 is *already inside* the `oxidative` module — `MODULES["oxidative"].transcription_factor`
is literally `"Yap1/Skn7"`. So the residual TRX2 induction in `yap1Δ` is **on-diagonal, not
off-diagonal**: it belongs in the 1.0 self-loading, not in a 0.30 crosstalk term to the ESR.

This is a topology error, not a magnitude error, and it changes `L`'s rank structure rather
than its scale. The source string in the file ("a Yap1-independent arm reads at this promoter")
is true; the module it is assigned to is not the one the literature names.

### 3.4 What the literature says about crosstalk that the model does *not* encode

**Keren 2013 (PMID 24169404) is the most important uncited result for this project.** Measuring
~900 yeast and ~1800 *E. coli* promoters with fluorescent reporters, it finds that **60–90 % of
promoters change expression between conditions by a global scaling factor that depends only on
the condition and not on the promoter's identity.** Specific regulation is the *deviation* from
that scale line.

This bears directly on `L`: a large fraction of what an uncorrected reporter panel would read as
crosstalk is not promoter content at all, it is the global scale factor. This repository already
does the right thing structurally — growth-dilution correction removes the dominant part of it —
but the result deserves to be cited where `_DEFAULT_BASAL` and the crosstalk weights are
justified, because it is the published evidence that the modelling choice is correct.

---

## 4. Ready-to-apply values

Ordered by (credibility restored) ÷ (effort). Each names the exact location.

| # | Location | Change | Basis |
|---|---|---|---|
| **1** | `src/ystwin/analysis/sensor_selection.py:180` — `_PLATE_NOISE = NoiseModel(relative_cv=0.05)` | **→ `relative_cv=0.146`** (or expose it and pass `OBSERVED_ACTIVITY_CV`) | The repo's own measurement, `generator/panel_experiment.py:38`. Verified by execution: module standard errors under `RECOMMENDED_BUILD` scale **exactly linearly** with this constant — ESR/UPR/oxidative/ATP go from 0.050/0.049/0.054/0.050 at cv=0.05 to **0.146/0.142/0.159/0.146** at the measured value. The *ranking* is invariant, so the build recommendation survives; every absolute precision claim and every `identifiable_by_precision(min_effect=…)` verdict is **3× optimistic** |
| **2** | `src/ystwin/analysis/design.py:42` — `relative_cv: float = 0.02` | Make it **required**, or default to 0.146 | 0.02 is **7.3×** tighter than measured, and it is the value a new caller silently gets. Same failure mode as `_MAINTENANCE_PER_ACTIVITY` in `SOURCE_AUDIT.md` §3.5 |
| **3** | `src/ystwin/generator/stress_panel.py:414–417` — `RIM101-alkaline`, `also_reads={"calcium": 0.30}` | **→ 1.5**, or renormalise to `{alkaline_ph: 0.4, calcium: 0.6}` | Petrezsélyová 2016, PMID 27362362: Crz1 supplies **≈ 60 %** of the early ENA1 alkaline response. As encoded, the model has the dominant arm as the minor one |
| **4** | `src/ystwin/generator/stress_panel.py:377–380` — `TRX2-oxidative` | Delete `also_reads={"ESR": 0.30}`; rewrite the source string to name **Skn7** and cite Morgan 1997, PMID 9118942 | §3.3. The current weight assigns a Skn7 effect to Msn2/4 |
| **5** | `src/ystwin/generator/stress_panel.py:289` — `H2O2` `ec50=0.5` | **→ 0.12–0.15 mM** for the `oxidative` arm; keep `peroxide`/`redox` where they are | Goulev 2017, PMID 28418333: Yap1 nuclear entry partial at 0.1 mM, saturating above 0.2 mM, S288C in SCD. **Consequence:** `healthy_ladder("H2O2", 5)` currently returns **0.034–0.55 mM**, so the top four rungs all sit at or above where the literature says the oxidative response is already flat. With the corrected EC50 the informative ladder is ~0.02–0.3 mM |
| **6** | `src/ystwin/generator/stress_panel.py:277` — `DTT` source string | Strike "the literature 1.0 mM from PMC7646510". **That paper contains no 1.0 mM value** — MacGilvray 2020 (PMID 32597660) used **2.5 mM in BY4741 in YPD**. Cite Pincus 2010 (PMID 20625545) for the only published DTT dose-response: **half-maximal UPR ≈ 2.2 mM** in W303a/2×SDC | §5 |
| **7** | `src/ystwin/generator/stress_panel.py:425–429` — `HyPer7` | Add explicit `demonstrated_in_yeast=True` with a citation: **Kritsiligkou 2023, PNAS 120:e2314043120, PMID 37991942** — a *S. cerevisiae* library with HyPer7 fused to every protein-coding ORF | Closes the audit's open item (§3.2, "UNVERIFIED for *S. cerevisiae*") with a definitive source |
| **8** | `src/ystwin/generator/literature.py:37` and `stress_panel.py:267` | `mu_max` **0.40 → 0.46 h⁻¹ if YPD, 0.30 h⁻¹ if synthetic**; state which | **BNID 108255** (~90 min YPD, ~140 min synthetic) and **BNID 106359** (0.37 h⁻¹ minimal, 30 °C), both in `data/kaggle/BioNumbers_Nov2024.csv`. Wire into `tests/test_constants_against_bionumbers.py`, which already opens that file |
| **9** | `src/ystwin/generator/literature.py:94–103` — `carrying_capacity` | Keep **3.0 g/L**, but attach the derivation: **BNID 106303** (max OD600 = 5.2 in YPD) × **BNID 108275** (0.62 gDCW/L per OD600) = 3.2 g/L | Moves it from bare literal to derived-and-checkable |
| **10** | `src/ystwin/generator/stress_panel.py:164` — `integrates_hours: float = 3.0` | Compute it as `1/(mu + k_deg)` instead of fixing it | It already *is* 1/µ at µ = 0.33 h⁻¹. Deriving it makes the constant transfer to any host for free, which is the stated goal |
| **11** | `src/ystwin/generator/literature.py:60` — `h2o2_sublethal` note | Change "sublethal" to what Goulev 2017 measured: at 0.5 mM, **22 % of cells adapt, 42 % arrest permanently** | PMID 28418333 |
| **12** | Wherever the 25 `lethal_dose` values are described | State that **23 of 25 are `3.0 × ec50`**, not independent values | Verified by execution. It removes 23 numbers from the "unsourced parameters" count honestly rather than by fiat |
| **13** | `src/ystwin/generator/literature.py` | Attach **Keren 2013, PMID 24169404** to `_DEFAULT_BASAL` and `promoter_basal`/`promoter_peak` | It is the published basis for a high basal floor and a small specific component, in both hosts this project might target |
| **14** | Documentation of the leucine effect | Note that **BY4741 grows poorly on synthetic complete because `leu2` blocks leucine uptake** (Cohen & Engelberg 2007, PMID 17573937) | A strain-specific growth penalty that the growth-dilution correction will otherwise absorb as if it were stress |

**Still outstanding from `SOURCE_AUDIT.md` §3.1** (verified again on 2026-08-26 against the
current file, and *not* yet fixed): the 14 broken module PMIDs are still present —
`stress_panel.py:210` (10361303 → Tovar 1999), `:213` (9663389 → Galant 1998), `:250`
(32160542 → **Laver 2020, *Cell Rep*, RNA-binding protein Rasputin/G3BP**; Pak 2020 is
**32130885**), `:253` (24815987 → **Vasefi 2014, dermoscopy**; Yaginuma is **25283467**),
`:259` (21982710 → **Krahmer 2011, lipid droplets**; Hung 2011 is **21982714**). The QUEEN-2m
*reporter* entry was repaired; the `MODULES` entries were not. One further slip: the QUEEN-2m
source cites "Takaine 2021, PMID 33654827" — that record resolves to **Takaine 2019**,
*Bio Protoc* 9:e3320. Right paper, wrong year.

---

## 5. What cannot be grounded, and what rests on it

### 5.1 The 76 target weights and the identifiability result

No published source gives them. §3.2 is a derivation, not a citation, and until it is run the
honest statement is the one `SOURCE_AUDIT.md` §7 item 8 proposes: *the topology is
literature-transcribed and PMID-cited; the weights are expert priors, and the identifiability
result is conditional on them.*

What actually rests on it: `reporter_loadings` → `fisher_information` → `module_standard_errors`
→ `RECOMMENDED_BUILD` → the README's "2 → 2 of 7, 4 → 4 of 7, 7 → 7 of 7" table and the
5→7-module transfer claim. Two mitigations are real and worth stating alongside the caveat:
(i) the weights enter `L` only through the six `also_reads` terms — the 76 `targets` weights
drive the *simulated data*, not `L` itself; and (ii) module standard errors were verified above
to scale linearly in the noise constant, so the *ordering* of which modules are identifiable is
robust to the one parameter we can check.

### 5.2 The DTT window: a genuine, unresolved contradiction

Three measurements disagree, and the disagreement is not a rounding error:

| Source | Strain | Medium | Result |
|---|---|---|---|
| This repository's plates | BY4741 constructs | (not stated in `DATA_INVENTORY.md`) | Growth **halves at 1.45–1.65 mM** |
| MacGilvray 2020, PMID 32597660 | **BY4741** | **YPD** | **2.5 mM** reduces growth but the culture acclimates to a new rate |
| Pincus 2010, PMID 20625545 | W303a | 2×SDC | UPR **< 10 % of max at ≤ 1.5 mM**; half-maximal at **2.2 mM**; usable to 5 mM |

The repo's plates put the *lethal* dose below where two published studies see cells still
growing, and its whole `healthy_ladder("DTT", 5)` (0.053–0.853 mM) lies inside the window where
Pincus reports almost no UPR — yet the plates measured 1.38–1.47 fold induction there.

Both cannot be right about the same organism. The candidate explanations, in order of how
cheaply they can be tested:

1. **Medium.** MacGilvray used YPD, which supplies cysteine, methionine and glutathione
   precursors that buffer a thiol reductant. A synthetic medium does not. **`DATA_INVENTORY.md`
   does not record the plate medium** — recording it is the single cheapest step here.
2. **Genotype.** BY4741 is `met15Δ0` (Brachmann 1998, PMID 9483801) and therefore methionine
   auxotrophic with a broken sulfate-assimilation branch. A strain that cannot make its own
   sulfur amino acids under a thiol-reductant challenge is a plausible sensitiser. *This is a
   hypothesis, not a literature finding — I found no paper testing met15Δ0 DTT sensitivity.*
3. **Readout.** The repo measures a promoter fusion over hours with growth-dilution correction;
   Pincus measured HAC1 splicing and UPRE-GFP. A 1.4-fold promoter-fusion change and "10 % of
   maximal HAC1 splicing" are not the same observable.

Until this is resolved, any statement that the DTT EC50 is "the literature 1.0 mM" is doubly
wrong: the cited paper does not contain 1.0 mM, and the only paper that does report a DTT
dose-response reports 2.2 mM.

### 5.3 Hill coefficients

`growth_hill` (1.1–1.6), `promoter_hill` (1.8–2.2) and `_LETHAL_HILL` (2.5) have **no published
values**. Searches for yeast stress dose-response Hill coefficients returned nothing usable.
They are fitted shape parameters and should be labelled as such. The repo already has the tool
to fit them (`panel_calibration.fit_dose_response`, `lethal_dose_from_growth`) and already
knows that one ladder fixes the shape but confounds the amplitude — which is the right
statement and is already in that module's docstring.

### 5.4 The heat parameterisation

`STRESSORS["heat"]` uses "degC above 30", but every ESR dataset shifts from **25 °C**, and Gasch
2000 showed 29 → 33 °C also elicits the response. The magnitude of a heat shock is a function of
the *shift*, not the absolute temperature, so a parameter denominated against a fixed 30 °C
baseline will not transfer to a lab that pre-grows at 25 °C, and cannot be compared with any of
the source data. This is a units problem, not a value problem, and it is invisible in the
current encoding.

### 5.5 What I could not find at all

| Wanted | Tried | Status |
|---|---|---|
| A published regulon × regulon overlap matrix for yeast stress | Gasch 2000, Causton 2001, Harbison 2004, Lee 2002, YEASTRACT+, Kemmeren 2014, Balaji 2006, Ihmels 2002, Segal 2003, Chasman 2014, iModulon/ICA searches | **Does not exist.** Balaji 2006 is the nearest published object; it scores TF-TF association significance, not loadings |
| Published EC50s for 14 of the 25 stressors | PubMed by agent + "EC50"/"IC50"/"dose-response" + *S. cerevisiae* | **Not found.** Conventional doses exist; half-maximal values do not |
| Hill coefficients for yeast stress dose-response | Same | **Not found** |
| Grably 2002 per-construct β-galactosidase values for HSP104 HSE vs STRE | Abstract via E-utilities; full text paywalled (Wiley) | **In the figures, not extracted** |
| Rep 2000 gene-level GPD1 split between Hot1 and Msn2/4 | Abstract via E-utilities; JBC full text returned HTTP 403 | **In the tables, not extracted** |
| Hahn 2006 relative YRE vs HSE contribution to RPN4 | Paywalled — the audit hit the same wall | **UNVERIFIED** |
| A quantitative *Komagataella phaffii* stress-regulon resource | PubMed + web | **Does not exist** at the level yeast and *E. coli* have. See §6 |

---

## 6. Porting to another host

The brief asks which parameters are host-specific and which transfer. The split is clean.

### 6.1 Transfers unchanged (structural)

- **`reporter_dilution_exponent = 1.0`** — steady-state concentration goes as α/µ. This is
  algebra, not biology (Hintsche & Klumpp 2013, PMID 24041253). It holds in any growing cell.
- **`integrates_hours`** — once it is computed as `1/(µ + k_deg)` rather than fixed at 3.0, it
  transfers for free. In *E. coli* at µ ≈ 1.0 h⁻¹ it becomes ~1 h, which changes the whole
  design: a promoter fusion in *E. coli* integrates over a third of the window it does in yeast.
- **The ratiometric-vs-transcriptional distinction** (`Kind`), the null-space check in
  `design.py`, and the growth-dilution correction itself. These are the deep parts.
- **The basal-floor structure.** Keren 2013 (PMID 24169404) measured both hosts with one method
  and found the same global-scaling result in each. A `_DEFAULT_BASAL` near 0.9 is not a yeast
  fact; it is a promoter fact.

### 6.2 Must be re-measured (host-specific)

Everything else: all 25 EC50s, all 25 lethal doses, all 76 target weights, all 6 `also_reads`
weights, all 4 cascade weights, `_MU_MAX`, `carrying_capacity`, and every Hill coefficient.
The *module list itself* does not survive: Hac1/Ire1 is conserved into *E. coli*'s… nothing.
There is no UPR in a prokaryote; the nearest analogue is σE/Cpx envelope stress, which is a
different mechanism with a different reporter set.

### 6.3 *Escherichia coli* — better resourced than yeast, for this specific model

This is worth stating plainly because it inverts the usual assumption. The one object yeast
does **not** have — a measured, quantitative module decomposition with both a loading matrix and
a condition-activity matrix — **exists for *E. coli***:

| Yeast parameter | *E. coli* equivalent source | What it gives |
|---|---|---|
| `MODULES` + `also_reads` (i.e. `L`) | **Sastry 2019, PMID 31797920** — ICA of 278 RNA-seq samples yielding **iModulons**, each a gene-weight vector; and **iModulonDB**, Rychel 2021, PMID 33045728 | ICA returns exactly the decomposition this repo builds by hand: a gene × module weight matrix **M** and a module × condition activity matrix **A**. `L` is measured, not assigned |
| Regulon membership / topology | **RegulonDB v12.0**, Salgado 2024, PMID 37971353 | The curated equivalent of YEASTRACT+, with far denser mechanistic annotation |
| `promoter_basal`, `promoter_peak`, reporter panel | **Zaslaver 2006, PMID 16862137** — 1,820 promoter-GFP fusions with measured absolute activities; **Zaslaver 2009, PMID 19851443** — the invariant distribution of those activities | Absolute promoter activity in real units, for the whole genome. Yeast has no equal-coverage equivalent |
| Cross-host anchor for basal/peak | **Keren 2013, PMID 24169404** | The *same* experiment in both hosts — ~900 yeast and ~1800 *E. coli* promoters. The only paper that lets a basal/peak parameter be ported with a measured scale factor rather than a guess |
| `_MU_MAX` | µ ≈ 0.6–1.0 h⁻¹ on glucose minimal, well above yeast's 0.3–0.46 | Changes `integrates_hours` by ~3×, which is a design change, not a re-parameterisation |
| Stress modules | OxyR (peroxide), SoxRS (superoxide), σ32/rpoH (heat), σE + Cpx (envelope), σS/rpoS (general stress), Fur (iron), stringent response (nutrient) | Structurally analogous to the panel's seven: a general-stress hub (σS ≈ Msn2/4), specific redox arms, an iron arm with its own logic, and a protein-quality arm |

**What would have to be re-measured to run this on *E. coli*:** every dose parameter, and the
reporter panel entirely (a yeast STRE-HSP12 fusion means nothing in *E. coli*). What would
*not*: the loading matrix, which for the first time could be estimated from data rather than
assigned — Sastry 2019's **M** matrix is `L`, measured, for 92 modules.

### 6.4 *Komagataella phaffii* (*Pichia pastoris*) — honestly, the pieces are not there

Ire1/Hac1 is conserved and characterised (HAC1 splicing, constitutive partial UPR activity even
without an ER stressor; `ire1Δ` and `hac1Δ` give only partially overlapping expression changes).
Promoter engineering toolboxes exist (Vogl 2022, PMID 35781205). **But there is no *P. pastoris*
equivalent of Gasch 2000, no YEASTRACT-grade regulatory database, and no measured stress-regulon
decomposition.** Porting there means measuring the loading matrix from scratch — which, given
§3.2, is arguably no worse off than yeast is today, just more honest about it.

### 6.5 The generalisable claim this supports

The transferable object is not the numbers; it is the **method**: define modules, measure a
loading matrix from perturbation data rather than assigning it, correct reporter readings for
growth dilution, and score sensor sets on Fisher information with a noise model taken from the
instrument rather than from hope. Two of those four are already right in this repository. The
loading matrix is assigned rather than measured (§3), and the noise model is 3–7× tighter than
the repo's own plates (§4, items 1–2). Fixing those two is what makes the claim general.

---

## 7. Annotated bibliography

Every entry below was resolved on 2026-08-26 through
`eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=<PMID>&retmode=json`, and the
title shown is the title that record returned.

### The environmental stress response and its datasets

- **PMID 11102521** — Gasch AP, 2000, *Mol Biol Cell* 11:4241-57 — "Genomic expression programs
  in the response of yeast cells to environmental changes." *Doses (verified from full text,
  PMC15070): H2O2 0.30 mM; menadione 1 mM bisulfite; DTT 2.5 mM; diamide 1.5 mM; sorbitol 1 M;
  heat 25→37 °C and 29→33 °C; nitrogen limitation 0.025 % ammonium sulfate. Strain DBY7286, YPD,
  30 °C. ESR ≈ 900 genes, ~300 induced / ~600 repressed. Data: GEO **GSE18**.*
- **PMID 11179418** — Causton HC, 2001, *Mol Biol Cell* 12:323-37 — "Remodeling of yeast genome
  expression in response to environmental changes." *Doses (verified, PMC30946): H2O2 0.4 mM;
  NaCl 1.0 M; sorbitol 1.5 M; heat 25→37 °C; acid pH 6→4 (succinic); alkali pH 6→7.9 (Tris).
  Strain ATCC-201388 (S288c), YPD. CER = **499 genes, 216 induced, 283 repressed**.*
- **PMID 17605132** — Gasch AP, 2007, *Yeast* 24:961-76 — "Comparative genomics of the
  environmental stress response in ascomycete fungi." *The cross-species framing for §6.*
- **PMID 42165638** — Gasch AP, 2026, *Yeast* — "Twenty-Five Years of the Environmental Stress
  Response and the Enduring Power of Yeast in Stress Biology." *Current review; the right
  citation for the ESR as a concept.*
- **PMID 25957506** — Ho YH, 2015, *Curr Genet* 61:503-11 — "Exploiting the yeast
  stress-activated signaling network to inform on stress biology and disease signaling."
- **PMID 30078561** — Ho YH, 2018, *Curr Biol* 28:2673-2680.e4 — "Decoupling Yeast Cell Division
  and Stress Defense Implicates mRNA Repression in Translational Reallocation during Stress."
- **PMID 22102822** — Berry DB, 2011, *PLoS Genet* 7:e1002353 — "Multiple means to the same end:
  the genetic basis of acquired stress resistance in yeast." *Relevant caution: pre-treatment
  with salt, DTT or heat produces H2O2 tolerance by **non-overlapping** genetic routes — i.e.
  cross-stress history changes the dose-response, which no single EC50 can carry.*
- **PMID 25024351** — Wanichthanarak K, 2014, *Database (Oxford)* 2014 — "yStreX: yeast stress
  expression database."

### Regulatory network resources (the crosstalk-matrix candidates)

- **PMID 15343339** — Harbison CT, 2004, *Nature* 431:99-104 — "Transcriptional regulatory code
  of a eukaryotic genome." *203 TFs, YPD + 12 other conditions. The best published object for
  `also_reads`, because it measures promoter content.*
- **PMID 12399584** — Lee TI, 2002, *Science* 298:799-804 — "Transcriptional regulatory networks
  in Saccharomyces cerevisiae."
- **PMID 24766815** — Kemmeren P, 2014, *Cell* 157:740-52 — "Large-scale genetic perturbations
  reveal regulatory networks and an abundance of gene-specific repressors." *1,484 deletion
  strains; the best weighted TF→gene object for yeast.*
- **PMID 16762362** — Balaji S, 2006, *J Mol Biol* 360:213-27 — "Comprehensive analysis of
  combinatorial regulation using the transcriptional regulatory network of yeast." *The one
  published yeast **TF × TF co-regulatory network**. Association significance, not loadings.*
- **PMID 31586406** — Monteiro PT, 2020, *Nucleic Acids Res* 48:D642-D649 — "YEASTRACT+: a
  portal for cross-species comparative genomics of transcription regulation in yeasts."
- **PMID 36350610** — Teixeira MC, 2023, *Nucleic Acids Res* 51:D785-D791 — "YEASTRACT+: a
  portal for the exploitation of global transcription regulation and metabolic model data in
  yeast biotechnology and pathogenesis." *Current citation; provides the Regulation Matrix tool.*
- **PMID 12134151** — Ihmels J, 2002, *Nat Genet* 31:370-7 — "Revealing modular organization in
  the yeast transcriptional network."
- **PMID 12740579** — Segal E, 2003, *Nat Genet* 34:166-76 — "Module networks: identifying
  regulatory modules and their condition-specific regulators from gene expression data."
- **PMID 25411400** — Chasman D, 2014, *Mol Syst Biol* 10:759 — "Pathway connectivity and
  signaling coordination in the yeast stress-activated signaling network." *Salt stress only;
  signalling, not loadings.*

### Per-edge primary sources for the crosstalk topology

- **PMID 11967066** — Grably MR, 2002, *Mol Microbiol* 44:21-35 — "HSF and Msn2/4p can
  exclusively or cooperatively activate the yeast HSP104 gene." *`HSE-heat.also_reads["ESR"]`.
  HSEs and STREs cooperate for maximal induction; in `msn2Δmsn4Δ` induction runs exclusively
  through HSEs; in `ras2Δ` derepression runs exclusively through STREs.*
- **PMID 8423809** — Kohno K, 1993, *Mol Cell Biol* 13:877-90 — "The promoter region of the yeast
  KAR2 (BiP) gene contains a regulatory domain that responds to the presence of unfolded proteins
  in the endoplasmic reticulum." *`UPRE-ER.also_reads["heat"]`.*
- **PMID 9133628** — Oka M, 1997, *J Biochem* 121:578-84 — "Saccharomyces cerevisiae KAR2 (BiP)
  gene expression is induced by loss of cytosolic HSP70/Ssa1p through a heat shock
  element-mediated pathway." *The functional demonstration that KAR2's HSE works.*
- **PMID 9118942** — Morgan BA, 1997, *EMBO J* 16:1035-44 — "The Skn7 response regulator controls
  gene expression in the oxidative stress response of the budding yeast Saccharomyces
  cerevisiae." *Skn7 binds the TRX2 promoter directly and co-operates with Yap1 on it. This is
  the correction in §3.3.*
- **PMID 8313910** — Kuge S, 1994, *EMBO J* 13:655-64 — "YAP1 dependent activation of TRX2 is
  essential for the response of Saccharomyces cerevisiae to oxidative stress by hydroperoxides."
- **PMID 10722658** — Rep M, 2000, *J Biol Chem* 275:8290-300 — "The transcriptional response of
  Saccharomyces cerevisiae to osmotic shock. Hot1p and Msn2p/Msn4p are required for the induction
  of subsets of high osmolarity glycerol pathway-dependent genes." *`STRE-osmotic.also_reads`.
  186 genes ≥3-fold induced; Msn2/4-dependent induction is itself reduced in `hog1Δ`, which is
  the mechanistic basis for the panel's Hog1→Msn2 routing.*
- **PMID 10409737** — Rep M, 1999, *Mol Cell Biol* 19:5474-85 — "Osmotic stress-induced gene
  expression in Saccharomyces cerevisiae requires Msn1p and the novel nuclear factor Hot1p."
- **PMID 27362362** — Petrezsélyová S, 2016, *PLoS One* 11:e0158424 — "Regulation of the
  Na+/K+-ATPase Ena1 Expression by Calcineurin/Crz1 under High pH Stress: A Quantitative Study."
  *The single best-quantified crosstalk edge in the panel: **Crz1 contributes ≈ 60 % of the early
  ENA1 alkaline response**, via two stress-responsive Crz1 sites; a second input of similar
  kinetics is attributed to Snf1.*
- **PMID 17023428** — Platara M, 2006, *J Biol Chem* 281:36632-42 — "The transcriptional response
  of the yeast Na(+)-ATPase ENA1 gene to alkaline stress involves three main signaling pathways."
  *Rim101(→Nrg1), Snf1(→Nrg1/Mig2) and calcineurin. In `snf1 rim101` + FK506, induction is
  abolished entirely — the three arms are exhaustive.*
- **PMID 27320198** — Solís EJ, 2016, *Mol Cell* 63:60-71 — "Defining the Essential Function of
  Yeast Hsf1 Reveals a Compact Transcriptional Program for Maintaining Eukaryotic Proteostasis."
  *The core Hsf1 regulon is **18 genes**, ~10× fewer than the >160 previously assumed; genome-wide
  heat-shock response in Hsf1-depleted vs wild-type cells correlates at R² = 0.81, i.e. most of
  the heat response is *not* Hsf1. HSP104 is among the 18. Anyone deriving an Hsf1 weight from a
  legacy regulon list will overestimate it by an order of magnitude.*
- **PMID 12676948** — Young ET, 2003, *J Biol Chem* 278:26146-58 — "Multiple pathways are
  co-regulated by the protein kinase Snf1 and the transcription factors Adr1 and Cat8."
  *The paper `MODULES["carbon"]` intends; the file still cites 12482873.*
- **PMID 16556235** — Hahn JS, 2006, *Mol Microbiol* 60:240-51 — "A stress regulatory network for
  co-ordinated activation of proteasome expression mediated by yeast heat shock transcription
  factor." *`MODULES["proteasome"].driven_by`. Verified in the audit; the quantitative split
  between YRE and HSE remains paywalled.*

### Dose-response

- **PMID 20625545** — Pincus D, 2010, *PLoS Biol* 8:e1000415 — "BiP binding to the ER-stress
  sensor Ire1 tunes the homeostatic behavior of the unfolded protein response." *The only
  published DTT dose-response for the yeast UPR. Titration 0.13–7.5 mM; ≤ 1.5 mM gives < 10 % of
  maximal activity, **2.2 mM ≈ 50 %**, 3.3 mM ≈ 75 %, 5 mM near-saturating. Strain YDP001
  (W303a), 2×SDC. In the `ire1` BiP-less mutant the curve shifts left and saturates by 1.5 mM.*
- **PMID 19026441** — Merksamer PI, 2008, *Cell* 135:933-47 — "Real-time redox measurements during
  endoplasmic reticulum stress reveal interlinked protein folding functions." *Companion source
  for DTT effects on ER redox.*
- **PMID 32597660** — MacGilvray ME, 2020, *J Proteome Res* 19:3405-3417 — "Phosphoproteome
  Response to Dithiothreitol Reveals Unique Versus Shared Features of Saccharomyces cerevisiae
  Stress Responses." ***This is `PMC7646510`.*** *Uses **2.5 mM DTT in BY4741 in YPD**; the dose
  reduces growth but permits acclimation to a new rate. It contains no 1.0 mM value, so the DTT
  `ec50` source string misattributes it.*
- **PMID 28418333** — Goulev Y, 2017, *eLife* 6:e23971 — "Nonlinear feedback drives homeostatic
  plasticity in H2O2 stress response." *The best H2O2 dose-response for yeast. S288C in SCD,
  microfluidics: 0.1 mM gives partial Yap1-GFP nuclear relocation with **no** growth-rate effect;
  nuclear localisation **saturates above 0.2 mM**; cell-cycle arrest emerges ~0.2 mM and growth
  arrest ~0.5 mM; **0.6 mM causes full growth arrest**. At 0.5 mM the population splits 22 %
  adapted / 36 % prolonged arrest / 42 % permanent arrest. Hormesis on replicative lifespan at
  10–25 µM.*
- **PMID 26944189** — Tomalin LE, 2016, *Free Radic Biol Med* 95:333-48 — "Increasing extracellular
  H2O2 produces a bi-phasic response in intracellular H2O2, with peroxiredoxin hyperoxidation only
  triggered once the cellular H2O2-buffering capacity is overwhelmed." *Explains why a bolus dose
  and a sustained dose give different lethal thresholds — the mechanistic reconciliation of the
  repo's 1.0 mM with Goulev's 0.6 mM.*
- **PMID 23740476** — Hoffmann GR, 2013, *Environ Mol Mutagen* 54:384-96 — "Adaptive response to
  hydrogen peroxide in yeast: induction, time course, and relationship to dose-response models."
  *Priming doses 0.000975–2 mM against a 1 mM challenge; adaptive response maximal at a priming
  dose of 0.125–0.25 mM in exponential culture.*
- **PMID 9712873** — Godon C, 1998, *J Biol Chem* 273:22480-9 — "The H2O2 stimulon in
  Saccharomyces cerevisiae." *The 0.4 mM convention.*
- **PMID 9885153** — Jamieson DJ, 1998, *Yeast* 14:1511-27 — "Oxidative stress responses of the
  yeast Saccharomyces cerevisiae." *Standing review.*

### Culture, growth and promoter activity

- **PMID 12489126** — Warringer J, 2003, *Yeast* 20:53-67 — "Automated screening in environmental
  arrays allows analysis of quantitative phenotypic profiles in Saccharomyces cerevisiae."
  ***The methodological match for this project's assay***: 350 µL microcultivation with frequent
  optical reads, 98 conditions from 33 growth inhibitors at three concentrations each, and
  explicit extraction of *rate of growth* and *stationary-phase OD increment* — i.e. µ_max and
  carrying capacity in a microplate. Strains W303, FY1679, CEN.PK.2. Also documents that
  micro-scale and 10 mL cultures differ in respiratory potential and stress-protein expression,
  which is directly relevant to a 96-well reporter assay.*
- **PMID 14676322** — Warringer J, 2003, *Proc Natl Acad Sci U S A* 100:15724-9 —
  "High-resolution yeast phenomics resolves different physiological features in the saline
  response."
- **PMID 15608218** — Fernandez-Ricaud L, 2005, *Nucleic Acids Res* 33:D369-73 — "PROPHECY--a
  database for high-resolution phenomics." *Quantitative microplate growth data for the deletion
  collection under environmental challenge — the queryable form of the above.*
- **PMID 21965291** — Zakrzewska A, 2011, *Mol Biol Cell* 22:4435-46 — "Genome-wide analysis of
  yeast stress survival and tolerance acquisition to analyze the central trade-off between growth
  rate and cellular robustness." *The growth-vs-defence trade-off this model's growth correction
  is built around.*
- **PMID 17573937** — Cohen R, 2007, *FEMS Microbiol Lett* 273:239-43 — "Commonly used
  Saccharomyces cerevisiae strains (e.g. BY4741, W303) are growth sensitive on synthetic complete
  medium due to poor leucine uptake." *A `leu2`-linked growth penalty on SC medium, rescued by
  BAP2, TAT1 or LEU2 overexpression. Not a stress effect, and a growth-dilution correction will
  absorb it as though it were.*
- **PMID 9483801** — Brachmann CB, 1998, *Yeast* 14:115-32 — "Designer deletion strains derived
  from Saccharomyces cerevisiae S288C: a useful set of strains and plasmids for PCR-mediated gene
  disruption and other applications." *BY4741's genotype of record.*
- **PMID 24169404** — Keren L, 2013, *Mol Syst Biol* 9:701 — "Promoters maintain their relative
  activity levels under different growth conditions." ***The most important uncited paper for this
  project.*** *~900 *S. cerevisiae* and ~1800 *E. coli* promoters measured with fluorescent
  reporters; **60–90 % of promoters change between conditions by a global scaling factor that
  depends only on the condition, not the promoter**, and specific regulation is the deviation from
  that scale line. This is the published justification for both the growth-corrected estimator and
  a high `_DEFAULT_BASAL`, and it is the only measurement that spans yeast and *E. coli* with one
  method.*
- **PMID 26355006** — Keren L, 2015, *Genome Res* 25:1893-902 — "Noise in gene expression is
  coupled to growth rate." *Bears on `DEFAULT_WELL_CV` and `biological_cv`: expression noise is
  not a constant, it is a function of µ — so a single CV across a dose ladder that also changes µ
  is structurally wrong, in a direction the repo can measure.*
- **PMID 24041253** — Hintsche M, 2013, *J Biol Eng* 7:22 — "Dilution and the theoretical
  description of growth-rate dependent gene expression." ***This is `PMC3847955`***; correct
  support for `reporter_dilution_exponent = 1.0`.

### Sensors

- **PMID 37991942** — Kritsiligkou P, 2023, *Proc Natl Acad Sci U S A* 120:e2314043120 —
  "Proteome-wide tagging with an H2O2 biosensor reveals highly localized and dynamic redox
  microenvironments." ***Settles the audit's open HyPer7 question***: a *S. cerevisiae* library
  with HyPer7 fused to the C-terminus of every protein-coding ORF, plus a SypHer7 redox-insensitive
  control library. HyPer7 is demonstrated in *S. cerevisiae*.
- **PMID 32130885** — Pak VV, 2020, *Cell Metab* 31:642-653.e6 — "Ultrasensitive Genetically
  Encoded Indicator for Hydrogen Peroxide Identifies Roles for the Oxidant in Cell Migration and
  Mitochondrial Function." *The real HyPer7 paper. `stress_panel.py:250` and `:428` cite 32160542,
  which is **Laver JD, 2020, *Cell Rep* 30:3353** — "The RNA-Binding Protein Rasputin/G3BP…".*
- **PMID 25283467** — Yaginuma H, 2014, *Sci Rep* 4:6522 — "Diversity in ATP concentrations in a
  single bacterial cell population revealed by quantitative single-cell imaging." *Applied
  correctly at the QUEEN-2m reporter; `MODULES["atp"]` (`:253`) still cites 24815987 =
  **Vasefi F, 2014, *Sci Rep* 4:4924**, a dermoscopy paper.*
- **PMID 30858198** — Takaine M, 2019, *J Cell Sci* 132 — "Reliable imaging of ATP in living
  budding and fission yeast." *The yeast demonstration of QUEEN.*
- **PMID 33654827** — Takaine M, **2019**, *Bio Protoc* 9:e3320 — "QUEEN-based Spatiotemporal ATP
  Imaging in Budding and Fission Yeast." *The file dates this 2021.*
- **PMID 21982714** — Hung YP, 2011, *Cell Metab* 14:545-54 — "Imaging cytosolic NADH-NAD(+) redox
  state with a genetically encoded fluorescent biosensor." *The real Peredox paper;
  `stress_panel.py:259` and `:451` cite 21982710 = **Krahmer N, 2011, *Cell Metab* 14:504**, on
  lipid-droplet phosphatidylcholine synthesis.*

### Other hosts

- **PMID 31797920** — Sastry AV, 2019, *Nat Commun* 10:5536 — "The Escherichia coli transcriptome
  mostly consists of independently regulated modules." *ICA of 278 RNA-seq samples into iModulons;
  the measured equivalent of this repo's hand-assigned `L`, plus a condition-activity matrix.*
- **PMID 33045728** — Rychel K, 2021, *Nucleic Acids Res* 49:D112-D120 — "iModulonDB: a
  knowledgebase of microbial transcriptional regulation derived from machine learning."
- **PMID 37971353** — Salgado H, 2024, *Nucleic Acids Res* 52:D255-D264 — "RegulonDB v12.0: a
  comprehensive resource of transcriptional regulation in E. coli K-12."
- **PMID 16862137** — Zaslaver A, 2006, *Nat Methods* 3:623-8 — "A comprehensive library of
  fluorescent transcriptional reporters for Escherichia coli."
- **PMID 19851443** — Zaslaver A, 2009, *PLoS Comput Biol* 5:e1000545 — "Invariant distribution of
  promoter activities in Escherichia coli."
- **PMID 35781205** — Vogl T, 2022, *Methods Mol Biol* 2513:153-177 — "Engineering of Promoters for
  Gene Expression in Pichia pastoris."

### BioNumbers entries usable today

All from `data/kaggle/BioNumbers_Nov2024.csv`, organism *Saccharomyces cerevisiae*, and therefore
wirable straight into `tests/test_constants_against_bionumbers.py`, which already opens that file.

| BNID | Quantity | Value | Bears on |
|---|---|---|---|
| 106359 | Growth rate on minimal medium at 30 °C | **0.37 h⁻¹** | `_MU_MAX`, `mu_max_glucose` |
| 108255 | Doubling time, "normal" laboratory haploid | **~90 min YPD (µ = 0.46 h⁻¹); ~140 min synthetic (µ = 0.30 h⁻¹)** | `_MU_MAX`, `_CONSTRUCTS[*].mu_max` |
| 110545 | Typical population doubling time | 90–120 min | same |
| 106303 | Maximum OD600 in YPD | **5.2** | `carrying_capacity` |
| 106302 | Maximum OD600 in YPGal | 6.87 | `carrying_capacity` |
| 108275 | OD600 to dry weight | **1 OD600 = 0.62 gDCW/L** | `GDCW_PER_OD` (0.3 vs 0.42 fork, audit §6 item 5) |
| 100986 | Cells/mL at OD600 = 1 | 3.0 × 10⁷ | same |
| 110953 | Total biomass yield (glucose + ethanol) | 54 g/mol = 0.30 gDW/g | `carrying_capacity` ceiling |
| 109863 | Cytosolic pH, exponential growth on glucose | **7.2** | `CYTOSOLIC_PH` — supports the *prose* ("resting is about 7.2") that the audit flagged as inconsistent with the tuple's 7.4. Both figures are in the literature; the range covers them |
| 111465 | Glutathione redox potential by Grx1-roGFP2 | −310 to −320 mV | The Kojer/Ayer 40 mV gap the panel already documents; BNID sits between them |
| 115096 | Lower limit for reliable GFP detection | **~1,400 molecules/cell** | `g1_optical.min_reporter_above_background` |
| 106766 / 106767 | Median transcription rate / 90 % range | 0.12 mRNA/min; 2.33–29.7 mRNA/h | Absolute units for `promoter_basal` / `promoter_peak` |
| 102974, 108677 | YFP maturation time | 39 min | `unmixing.PANEL` (audit §3.4) |
| 106883, 108678 | CFP maturation time | 49 min | same |

---

## 8. Method note

PMIDs were resolved through NCBI E-utilities `esummary.fcgi` and, where an abstract was needed,
`efetch.fcgi` with `rettype=abstract`. Full texts were read through PubMed Central where open
(PMC15070, PMC30946, PMC2897766, PMC4938784) and through the publisher where open (eLife 23971).
Where a full text was paywalled it is marked as such and the claim is not asserted. Several
plausible-looking PMIDs guessed from memory during this work resolved to unrelated papers and
were discarded before use — which is the same failure mode the source audit found in the file,
and the reason every identifier here was resolved before it was written down.
